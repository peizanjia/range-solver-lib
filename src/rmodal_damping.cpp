#include "rmodal_damping.h"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

#include <rbl_error.h>
#include <rml_eigen_value_solver_conf.h>

#include "reigenvaluesolver.h"

namespace
{
constexpr double EigenResidualTolerance = 1.0e-8;
constexpr double MassOrthogonalityTolerance = 1.0e-8;
constexpr double ClusterRelativeGap = 1.0e-6;
constexpr double EqualRatioTolerance = 1.0e-12;
constexpr double PositivePivotTolerance = 1.0e-12;

double norm(const RRVector &values)
{
    return RRVector::euclideanNorm(values);
}

RRVector product(const RSparseMatrix &matrix, const RRVector &values)
{
    RRVector result;
    RSparseMatrix::mlt(matrix,values,result);
    return result;
}

bool sameCluster(double left, double right)
{
    return std::fabs(left-right) <= ClusterRelativeGap*std::max(left,right);
}

void requirePositiveDefinite(const RSparseMatrix &matrix, const char *label)
{
    const uint n = matrix.getNRows();
    RRVector diagonal(n,0.0);
    for (uint i=0;i<n;i++)
    {
        diagonal[i] = matrix.findValue(i,i);
        if (!std::isfinite(diagonal[i]) || diagonal[i] <= 0.0)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Modal damping requires positive definite %s; a zero/negative diagonal may indicate a rigid-body or massless DOF.",label);
        }
    }

    // Work on diagonally equilibrated values so mixed translation/rotation
    // units do not turn a valid pivot into a false absolute-zero test.
    RRMatrix factor(n,n,0.0);
    for (uint i=0;i<n;i++)
    {
        for (uint j=0;j<=i;j++)
        {
            const double aij = matrix.findValue(i,j);
            const double aji = matrix.findValue(j,i);
            const double scale = std::sqrt(diagonal[i])*std::sqrt(diagonal[j]);
            if (!std::isfinite(aij) || !std::isfinite(aji) ||
                !std::isfinite(scale) || scale <= 0.0 ||
                std::fabs(aij-aji) > 1.0e-10*std::max({scale,std::fabs(aij),std::fabs(aji)}))
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Modal damping requires a finite symmetric %s.",label);
            }
            double pivot = aij/scale;
            for (uint k=0;k<j;k++) pivot -= factor[i][k]*factor[j][k];
            if (i == j)
            {
                if (!std::isfinite(pivot) || pivot <= PositivePivotTolerance)
                {
                    throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                                 "Modal damping requires positive definite %s; free rigid-body or near-mechanism modes are not supported.",label);
                }
                factor[i][j] = std::sqrt(pivot);
            }
            else
            {
                factor[i][j] = pivot/factor[j][j];
            }
        }
    }
}
}

RModalDampingResult RModalDamping::synthesize(const RSparseMatrix &mass,
                                              const RSparseMatrix &stiffness,
                                              const std::vector<double> &ratios,
                                              const RMatrixSolverConf &linearSolver)
{
    return synthesize(mass,stiffness,ratios,linearSolver,Options{});
}

RModalDampingResult RModalDamping::synthesize(const RSparseMatrix &mass,
                                              const RSparseMatrix &stiffness,
                                              const std::vector<double> &ratios,
                                              const RMatrixSolverConf &linearSolver,
                                              Options options)
{
    const uint n = mass.getNRows();
    const uint p = uint(ratios.size());
    if (n == 0 || stiffness.getNRows() != n || n > MaxEquations || p > MaxModes || p > n)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Modal damping requires 1..512 reduced equations and no more than 16 configured modes or available equations.");
    }
    if (options.eigenIterations == 0 || !std::isfinite(options.eigenConvergence) ||
        options.eigenConvergence <= 0.0)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Modal damping eigen solver settings must be finite and positive.");
    }
    for (double ratio : ratios)
    {
        if (!std::isfinite(ratio) || ratio < 0.0 || ratio >= 1.0)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Modal damping ratios must be finite and within [0,1).");
        }
    }
    RModalDampingResult result;
    result.damping.setNRows(n);
    if (p == 0) return result;

    requirePositiveDefinite(mass,"mass matrix");
    requirePositiveDefinite(stiffness,"stiffness matrix");

    // The current Lanczos path has no eigenvectors. Always use the block
    // subspace method, including one extra pair to detect a split cluster.
    const uint wanted = std::min(n,p+1);
    REigenValueSolverConf eigen(REigenValueSolverConf::SubspaceIteration);
    eigen.setNEigenValues(wanted);
    eigen.setNIterations(options.eigenIterations);
    eigen.setSolverCvgValue(options.eigenConvergence);
    eigen.setOutputFrequency(0);
    REigenValueSolver solver(eigen,linearSolver);
    RRVector eigenvalues;
    RRMatrix rawModes;
    solver.solve(mass,stiffness,eigenvalues,rawModes);
    if (eigenvalues.size() != wanted || rawModes.getNRows() != wanted ||
        rawModes.getNColumns() != n)
    {
        throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                     "Modal damping eigen solver returned fewer eigenpairs than needed, including the guard mode.");
    }
    for (uint i=0;i<wanted;i++)
    {
        if (!std::isfinite(eigenvalues[i]) || eigenvalues[i] <= 0.0 ||
            (i > 0 && eigenvalues[i] < eigenvalues[i-1]))
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Modal damping requires finite positive ordered eigenvalues; rigid-body modes are unsupported.");
        }
    }
    for (uint i=1;i<wanted;i++)
    {
        if (!sameCluster(eigenvalues[i-1],eigenvalues[i])) continue;
        if (i == p || std::fabs(ratios[i]-ratios[i-1]) > EqualRatioTolerance)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Modal damping cannot split a repeated/near-repeated eigenvalue cluster or assign different damping ratios within it.");
        }
    }

    std::vector<RRVector> modes;
    std::vector<RRVector> massModes;
    modes.reserve(wanted);
    massModes.reserve(wanted);
    for (uint i=0;i<wanted;i++)
    {
        RRVector mode(n,0.0);
        for (uint row=0;row<n;row++)
        {
            mode[row] = rawModes[i][row];
            if (!std::isfinite(mode[row]))
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Modal damping eigen solver returned a non-finite mode shape.");
            }
        }
        // Only equal-eigenvalue clusters may be mixed without changing their
        // physical meaning. Two M-inner-product passes stabilize that basis.
        for (uint pass=0;pass<2;pass++)
        {
            for (uint j=0;j<i;j++)
            {
                if (!sameCluster(eigenvalues[i],eigenvalues[j])) continue;
                const double projection = RRVector::dot(mode,massModes[j]);
                for (uint row=0;row<n;row++) mode[row] -= projection*modes[j][row];
            }
        }
        RRVector massMode = product(mass,mode);
        const double modalMass = RRVector::dot(mode,massMode);
        if (!std::isfinite(modalMass) || modalMass <= 0.0)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Modal damping eigenvector has non-positive modal mass.");
        }
        const double inverseNorm = 1.0/std::sqrt(modalMass);
        for (uint row=0;row<n;row++)
        {
            mode[row] *= inverseNorm;
            massMode[row] *= inverseNorm;
        }
        const RRVector stiffnessMode = product(stiffness,mode);
        RRVector residual = stiffnessMode;
        for (uint row=0;row<n;row++) residual[row] -= eigenvalues[i]*massMode[row];
        const double denominator = norm(stiffnessMode)+std::fabs(eigenvalues[i])*norm(massMode);
        if (!std::isfinite(denominator) || denominator <= 0.0 ||
            !std::isfinite(norm(residual)) || norm(residual)/denominator > EigenResidualTolerance)
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Modal damping eigenpair failed the generalized K*phi=lambda*M*phi residual check.");
        }
        for (uint j=0;j<i;j++)
        {
            const double overlap = RRVector::dot(mode,massModes[j]);
            if (!std::isfinite(overlap) || std::fabs(overlap) > MassOrthogonalityTolerance)
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Modal damping mode shapes failed mass orthogonality.");
            }
        }
        modes.push_back(std::move(mode));
        massModes.push_back(std::move(massMode));
    }

    result.eigenvalues.resize(p,0.0);
    result.modes.resize(p,n,0.0);
    RRMatrix denseDamping(n,n,0.0);
    for (uint i=0;i<p;i++)
    {
        result.eigenvalues[i] = eigenvalues[i];
        for (uint row=0;row<n;row++) result.modes[i][row] = modes[i][row];
        const double coefficient = 2.0*ratios[i]*std::sqrt(eigenvalues[i]);
        if (!std::isfinite(coefficient))
        {
            throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                         "Modal damping coefficient overflowed.");
        }
        for (uint row=0;row<n;row++)
        {
            for (uint column=0;column<=row;column++)
            {
                denseDamping[row][column] += coefficient*massModes[i][row]*massModes[i][column];
            }
        }
    }
    for (uint row=0;row<n;row++)
    {
        for (uint column=0;column<=row;column++)
        {
            const double value = denseDamping[row][column];
            if (!std::isfinite(value))
            {
                throw RError(RError::Type::InvalidInput,R_ERROR_REF,
                             "Modal damping matrix overflowed.");
            }
            if (value == 0.0) continue;
            result.damping.addValue(row,column,value);
            if (row != column) result.damping.addValue(column,row,value);
        }
    }
    return result;
}

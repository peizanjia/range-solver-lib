#include "rstructural_cable_equilibrium.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <rbl_rvector.h>
#include <rml_sparse_matrix.h>

namespace
{
constexpr std::size_t noFreeDof = std::numeric_limits<std::size_t>::max();
constexpr double armijoFraction = 1.0e-4;
constexpr double lineSearchContraction = 0.5;

struct NetworkState
{
    std::vector<double> internalForces;
    std::vector<double> freeResidual;
    std::vector<double> tensions;
    RSparseMatrix tangent;
    double residualNorm = 0.0;
    double freeLoadNorm = 0.0;
    double energy = 0.0;
    double freeDofPotential = 0.0;
    //! Sum of magnitudes contributing to the variable part of potential,
    //! used only to bound line-search roundoff acceptance.
    double freePotentialScale = 0.0;
};

double stableNorm(const std::vector<double> &values)
{
    double norm = 0.0;
    for (double value : values)
    {
        if (!std::isfinite(value))
            throw std::runtime_error("Cable equilibrium encountered a non-finite vector value");
        norm = std::hypot(norm,value);
    }
    if (!std::isfinite(norm))
        throw std::runtime_error("Cable equilibrium vector norm is non-finite");
    return norm;
}

double dot(const std::vector<double> &left, const std::vector<double> &right)
{
    if (left.size() != right.size())
        throw std::logic_error("Cable equilibrium vector dimensions do not match");
    double value = 0.0;
    for (std::size_t i=0;i<left.size();i++)
    {
        value += left[i]*right[i];
        if (!std::isfinite(value))
            throw std::runtime_error("Cable equilibrium dot product is non-finite");
    }
    return value;
}

void validateInputs(
    const std::vector<RStructuralCableVector3> &referenceNodes,
    const std::vector<RStructuralCableEquilibriumMember> &members,
    const std::vector<bool> &prescribedMask,
    const std::vector<double> &prescribedDisplacements,
    const std::vector<double> &deadNodalLoads,
    const RStructuralCableEquilibriumOptions &options)
{
    if (referenceNodes.empty() || referenceNodes.size() >
        std::numeric_limits<std::size_t>::max()/3)
        throw std::invalid_argument("Cable equilibrium requires a valid nonempty node array");
    const std::size_t dofCount = referenceNodes.size()*3;
    if (dofCount > std::numeric_limits<uint>::max() ||
        prescribedMask.size() != dofCount || prescribedDisplacements.size() != dofCount ||
        deadNodalLoads.size() != dofCount)
        throw std::invalid_argument("Cable equilibrium prescribed/load vectors must have size 3N");
    if (options.maxIterations == 0 || options.maxLineSearchIterations == 0 ||
        !std::isfinite(options.absForceTol) || options.absForceTol < 0.0 ||
        !std::isfinite(options.relForceTol) || options.relForceTol < 0.0 ||
        (options.absForceTol == 0.0 && options.relForceTol == 0.0))
        throw std::invalid_argument("Invalid cable equilibrium iteration or force tolerance options");

    for (const auto &node : referenceNodes)
    {
        for (double value : node)
        {
            if (!std::isfinite(value))
                throw std::invalid_argument("Cable equilibrium reference coordinates must be finite");
        }
    }
    for (std::size_t dof=0;dof<dofCount;dof++)
    {
        if (!std::isfinite(prescribedDisplacements[dof]) || !std::isfinite(deadNodalLoads[dof]))
            throw std::invalid_argument("Cable equilibrium displacement and load vectors must be finite");
    }
    for (const auto &member : members)
    {
        if (member.node1 >= referenceNodes.size() || member.node2 >= referenceNodes.size() ||
            member.node1 == member.node2 || !std::isfinite(member.referenceLength) ||
            member.referenceLength <= 0.0 || !std::isfinite(member.initialTension) ||
            member.initialTension < 0.0 || !std::isfinite(member.axialRigidity) ||
            member.axialRigidity <= 0.0)
            throw std::invalid_argument("Cable equilibrium member has invalid nodes or material data");
        const double axialStiffness = member.axialRigidity/member.referenceLength;
        if (!std::isfinite(axialStiffness) || axialStiffness <= 0.0)
            throw std::invalid_argument("Cable equilibrium member EA/L0 must be positive and finite");
    }
}

RStructuralCableVector3 currentPosition(
    const std::vector<RStructuralCableVector3> &referenceNodes,
    const std::vector<double> &displacements,
    std::size_t node)
{
    RStructuralCableVector3 position{};
    for (std::size_t axis=0;axis<3;axis++)
    {
        const std::size_t dof = node*3+axis;
        position[axis] = referenceNodes[node][axis]+displacements[dof];
        if (!std::isfinite(position[axis]))
            throw std::invalid_argument("Cable equilibrium current coordinates are non-finite");
    }
    return position;
}

NetworkState evaluateNetwork(
    const std::vector<RStructuralCableVector3> &referenceNodes,
    const std::vector<RStructuralCableEquilibriumMember> &members,
    const std::vector<double> &displacements,
    const std::vector<double> &loads,
    const std::vector<std::size_t> &globalToFree,
    std::size_t freeDofCount,
    bool assembleTangent)
{
    NetworkState state;
    state.internalForces.assign(displacements.size(),0.0);
    state.freeResidual.assign(freeDofCount,0.0);
    state.tensions.reserve(members.size());
    if (assembleTangent) state.tangent.setNRows(uint(freeDofCount));

    for (const auto &member : members)
    {
        const auto response = RStructuralCableKernel::evaluate(
            currentPosition(referenceNodes,displacements,member.node1),
            currentPosition(referenceNodes,displacements,member.node2),
            member.referenceLength,member.initialTension,member.axialRigidity);
        state.tensions.push_back(response.tension);
        state.energy += response.energy;
        state.freePotentialScale += response.energy;
        if (!std::isfinite(state.energy) || !std::isfinite(state.freePotentialScale))
            throw std::invalid_argument("Cable network strain energy is non-finite");

        const std::array<std::size_t,6> globalDofs = {{
            member.node1*3,member.node1*3+1,member.node1*3+2,
            member.node2*3,member.node2*3+1,member.node2*3+2
        }};
        for (std::size_t local=0;local<6;local++)
        {
            const std::size_t global = globalDofs[local];
            state.internalForces[global] += response.internalForce[local];
            if (!std::isfinite(state.internalForces[global]))
                throw std::invalid_argument("Cable network internal force is non-finite");
        }

        if (assembleTangent)
        {
            for (std::size_t row=0;row<6;row++)
            {
                const std::size_t freeRow = globalToFree[globalDofs[row]];
                if (freeRow == noFreeDof) continue;
                for (std::size_t column=0;column<6;column++)
                {
                    const std::size_t freeColumn = globalToFree[globalDofs[column]];
                    if (freeColumn == noFreeDof) continue;
                    const double value = response.tangent[row][column];
                    if (value == 0.0) continue;
                    const double sum = state.tangent.findValue(uint(freeRow),uint(freeColumn))+value;
                    if (!std::isfinite(sum))
                        throw std::invalid_argument("Cable network tangent assembly is non-finite");
                    state.tangent.addValue(uint(freeRow),uint(freeColumn),value);
                }
            }
        }
    }

    double externalWork = 0.0;
    double freeExternalWork = 0.0;
    for (std::size_t dof=0;dof<displacements.size();dof++)
    {
        const double loadWork = loads[dof]*displacements[dof];
        externalWork += loadWork;
        if (!std::isfinite(externalWork))
            throw std::invalid_argument("Cable network dead-load potential is non-finite");
        const std::size_t freeIndex = globalToFree[dof];
        if (freeIndex != noFreeDof)
        {
            state.freeResidual[freeIndex] = state.internalForces[dof]-loads[dof];
            state.freeLoadNorm = std::hypot(state.freeLoadNorm,loads[dof]);
            freeExternalWork += loadWork;
            state.freePotentialScale += std::fabs(loadWork);
            if (!std::isfinite(freeExternalWork) || !std::isfinite(state.freePotentialScale))
                throw std::invalid_argument("Cable network free-load potential is non-finite");
        }
    }
    state.freeDofPotential = state.energy-freeExternalWork;
    state.energy -= externalWork;
    if (!std::isfinite(state.energy) || !std::isfinite(state.freeDofPotential) ||
        !std::isfinite(state.freeLoadNorm))
        throw std::invalid_argument("Cable network total potential or load norm is non-finite");
    state.residualNorm = stableNorm(state.freeResidual);
    return state;
}

std::vector<double> sparseMultiply(const RSparseMatrix &matrix,
                                  const std::vector<double> &vector)
{
    RRVector x(uint(vector.size()),0.0);
    for (std::size_t i=0;i<vector.size();i++) x[uint(i)] = vector[i];
    RRVector y;
    RSparseMatrix::mlt(matrix,x,y);
    return std::vector<double>(y.begin(),y.end());
}

std::vector<double> solveJacobiCg(const RSparseMatrix &matrix,
                                  const std::vector<double> &rhs)
{
    const std::size_t n = rhs.size();
    std::vector<double> solution(n,0.0);
    if (n == 0) return solution;

    std::vector<double> diagonal(n,0.0);
    for (std::size_t i=0;i<n;i++)
    {
        diagonal[i] = matrix.findValue(uint(i),uint(i));
        if (!std::isfinite(diagonal[i]) || diagonal[i] <= 0.0)
            throw std::runtime_error("Cable tangent is singular or has a nonpositive Jacobi diagonal");
    }

    std::vector<double> residual = rhs;
    const double rhsNorm = stableNorm(rhs);
    if (rhsNorm == 0.0) return solution;
    // Keep the linear tolerance relative to this Newton rhs. An absolute
    // floor could make a small but physically required correction look solved.
    const double linearTolerance = 1.0e-10*rhsNorm;

    std::vector<double> preconditioned(n,0.0);
    std::vector<double> direction(n,0.0);
    for (std::size_t i=0;i<n;i++)
    {
        preconditioned[i] = residual[i]/diagonal[i];
        direction[i] = preconditioned[i];
        if (!std::isfinite(preconditioned[i]))
            throw std::runtime_error("Cable tangent Jacobi preconditioner is non-finite");
    }
    double rho = dot(residual,preconditioned);
    if (!std::isfinite(rho) || rho <= 0.0)
        throw std::runtime_error("Cable tangent has no positive preconditioned residual");

    const std::size_t maxCgIterations = n >= 2500 ? 10000 : std::max<std::size_t>(100,4*n);
    for (std::size_t iteration=0;iteration<maxCgIterations;iteration++)
    {
        const std::vector<double> product = sparseMultiply(matrix,direction);
        const double curvature = dot(direction,product);
        if (!std::isfinite(curvature) || curvature <= 0.0)
            throw std::runtime_error("Cable tangent is singular or has nonpositive CG curvature");
        const double alpha = rho/curvature;
        if (!std::isfinite(alpha))
            throw std::runtime_error("Cable tangent CG step is non-finite");
        for (std::size_t i=0;i<n;i++)
        {
            solution[i] += alpha*direction[i];
            residual[i] -= alpha*product[i];
            if (!std::isfinite(solution[i]) || !std::isfinite(residual[i]))
                throw std::runtime_error("Cable tangent CG iterate is non-finite");
        }
        if (stableNorm(residual) <= linearTolerance) break;

        for (std::size_t i=0;i<n;i++)
        {
            preconditioned[i] = residual[i]/diagonal[i];
            if (!std::isfinite(preconditioned[i]))
                throw std::runtime_error("Cable tangent CG preconditioned residual is non-finite");
        }
        const double nextRho = dot(residual,preconditioned);
        if (!std::isfinite(nextRho) || nextRho <= 0.0)
            throw std::runtime_error("Cable tangent CG lost positive residual norm");
        const double beta = nextRho/rho;
        if (!std::isfinite(beta))
            throw std::runtime_error("Cable tangent CG recurrence is non-finite");
        for (std::size_t i=0;i<n;i++) direction[i] = preconditioned[i]+beta*direction[i];
        rho = nextRho;
    }

    // Never trust only the recursive CG residual: recalculate b-A*x from the
    // assembled sparse tangent before accepting a Newton direction.
    const std::vector<double> product = sparseMultiply(matrix,solution);
    std::vector<double> freshResidual(n,0.0);
    for (std::size_t i=0;i<n;i++) freshResidual[i] = rhs[i]-product[i];
    const double freshNorm = stableNorm(freshResidual);
    if (freshNorm > 10.0*linearTolerance)
        throw std::runtime_error("Cable tangent CG failed fresh residual certification");
    return solution;
}
}

RStructuralCableEquilibriumResult RStructuralCableEquilibriumSolver::solve(
    const std::vector<RStructuralCableVector3> &referenceNodes,
    const std::vector<RStructuralCableEquilibriumMember> &members,
    const std::vector<bool> &prescribedMask,
    const std::vector<double> &prescribedDisplacements,
    const std::vector<double> &deadNodalLoads,
    const RStructuralCableEquilibriumOptions &options)
{
    validateInputs(referenceNodes,members,prescribedMask,prescribedDisplacements,
                   deadNodalLoads,options);
    const std::size_t dofCount = referenceNodes.size()*3;
    std::vector<std::size_t> globalToFree(dofCount,noFreeDof);
    std::vector<std::size_t> freeToGlobal;
    for (std::size_t dof=0;dof<dofCount;dof++)
    {
        if (!prescribedMask[dof])
        {
            globalToFree[dof] = freeToGlobal.size();
            freeToGlobal.push_back(dof);
        }
    }

    std::vector<double> displacements(dofCount,0.0);
    for (std::size_t dof=0;dof<dofCount;dof++)
    {
        if (prescribedMask[dof]) displacements[dof] = prescribedDisplacements[dof];
    }

    std::size_t iterations = 0;
    for (;;)
    {
        const NetworkState state = evaluateNetwork(
            referenceNodes,members,displacements,deadNodalLoads,globalToFree,
            freeToGlobal.size(),true);
        const double forceTolerance = options.absForceTol+
                                      options.relForceTol*state.freeLoadNorm;
        if (!std::isfinite(forceTolerance))
            throw std::invalid_argument("Cable equilibrium force tolerance is non-finite");
        if (state.residualNorm <= forceTolerance)
        {
            RStructuralCableEquilibriumResult result;
            result.displacements = displacements;
            result.internalForces = state.internalForces;
            result.elementTensions = state.tensions;
            result.residualNorm = state.residualNorm;
            result.energy = state.energy;
            result.iterations = iterations;
            result.reactions.resize(dofCount,0.0);
            for (std::size_t dof=0;dof<dofCount;dof++)
            {
                result.reactions[dof] = state.internalForces[dof]-deadNodalLoads[dof];
                if (!std::isfinite(result.reactions[dof]))
                    throw std::runtime_error("Cable equilibrium reaction is non-finite");
            }
            return result;
        }
        if (iterations >= options.maxIterations)
            throw std::runtime_error("Cable network failed to converge within the Newton iteration limit");

        std::vector<double> rhs(state.freeResidual.size(),0.0);
        for (std::size_t i=0;i<rhs.size();i++) rhs[i] = -state.freeResidual[i];
        const std::vector<double> increment = solveJacobiCg(state.tangent,rhs);
        const double descentSlope = dot(state.freeResidual,increment);
        if (!std::isfinite(descentSlope) || descentSlope >= 0.0)
            throw std::runtime_error("Cable tangent solve did not produce a descent direction");

        bool accepted = false;
        double stepLength = 1.0;
        for (std::size_t lineSearch=0;lineSearch<options.maxLineSearchIterations;lineSearch++)
        {
            std::vector<double> trial = displacements;
            for (std::size_t free=0;free<freeToGlobal.size();free++)
            {
                const std::size_t global = freeToGlobal[free];
                trial[global] += stepLength*increment[free];
                if (!std::isfinite(trial[global]))
                    throw std::runtime_error("Cable equilibrium trial displacement is non-finite");
            }
            try
            {
                const NetworkState trialState = evaluateNetwork(
                    referenceNodes,members,trial,deadNodalLoads,globalToFree,
                    freeToGlobal.size(),false);
                // Prescribed-load work is constant along admissible trials;
                // omit it from the Armijo arithmetic to avoid cancellation.
                const double armijoBound = state.freeDofPotential+
                                            armijoFraction*stepLength*descentSlope;
                const bool armijoAccepted = std::isfinite(armijoBound) &&
                    trialState.freeDofPotential <= armijoBound;

                // Near a converged equilibrium the predicted potential drop
                // can be much smaller than the rounding error in summing a
                // large prestress energy. Permit only a bounded 32-epsilon
                // energy shortfall, and only when the independently evaluated
                // physical residual drops by at least 10% or reaches the
                // unchanged outer force tolerance. This is a local
                // roundoff-resolution rule, not a relaxation of convergence.
                bool roundoffAccepted = false;
                if (!armijoAccepted && std::isfinite(armijoBound))
                {
                    const double energyScale = std::max(
                        state.freePotentialScale,trialState.freePotentialScale);
                    const double energyRoundoff = 32.0*std::numeric_limits<double>::epsilon()*energyScale;
                    const double energyShortfall = trialState.freeDofPotential-armijoBound;
                    const bool residualReduced = trialState.residualNorm <= forceTolerance ||
                        trialState.residualNorm <= 0.9*state.residualNorm;
                    roundoffAccepted = std::isfinite(energyRoundoff) &&
                        std::isfinite(energyShortfall) && energyShortfall >= 0.0 &&
                        energyShortfall <= energyRoundoff && residualReduced;
                }
                if (armijoAccepted || roundoffAccepted)
                {
                    displacements = std::move(trial);
                    accepted = true;
                    break;
                }
            }
            catch (const std::invalid_argument &)
            {
                // A finite backtrack may recover from a degenerate/overflowing
                // trial geometry; the current accepted state remains untouched.
            }
            catch (const std::runtime_error &)
            {
                // Trial residual norm overflow is also a rejected trial; the
                // accepted state and its convergence certificate are unchanged.
            }
            stepLength *= lineSearchContraction;
        }
        if (!accepted)
            throw std::runtime_error("Cable network Armijo line search exhausted without an acceptable step");
        iterations++;
    }
}

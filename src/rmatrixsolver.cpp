#include <cmath>
#include <mutex>
#include <unordered_map>
#include <vector>

#include <omp.h>

#include <rbl_utils.h>

#include "rmatrixsolver.h"

namespace
{

struct SparseMatrixCache
{
    uint nRows = 0;
    std::vector<uint> rowPtr;
    std::vector<uint> columnIndexes;
    std::vector<double> values;
    double norm = 0.0;
};

SparseMatrixCache buildSparseMatrixCache(const RSparseMatrix &A)
{
    SparseMatrixCache cache;
    cache.nRows = A.getNRows();
    cache.rowPtr.resize(cache.nRows+1,0);

    for (uint i=0;i<cache.nRows;i++)
    {
        cache.rowPtr[i+1] = cache.rowPtr[i] + A.getVector(i).size();
    }

    cache.columnIndexes.resize(cache.rowPtr.back());
    cache.values.resize(cache.rowPtr.back());

    double norm2 = 0.0;
#pragma omp parallel for default(shared) reduction(+:norm2)
    for (int64_t i=0;i<int64_t(cache.nRows);i++)
    {
        const RSparseVector<double> &row = A.getVector(uint(i));
        uint begin = cache.rowPtr[uint(i)];
        for (uint j=0;j<row.size();j++)
        {
            double value = row.getValue(j);
            cache.columnIndexes[begin+j] = row.getIndex(j);
            cache.values[begin+j] = value;
            norm2 += value*value;
        }
    }
    cache.norm = std::sqrt(norm2);

    return cache;
}

bool refreshSparseMatrixCache(const RSparseMatrix &A, SparseMatrixCache &cache)
{
    if (cache.nRows != A.getNRows() || cache.rowPtr.size() != A.getNRows()+1)
    {
        return false;
    }

    double norm2 = 0.0;
    int samePattern = 1;
#pragma omp parallel for default(shared) reduction(+:norm2) reduction(&:samePattern)
    for (int64_t i=0;i<int64_t(A.getNRows());i++)
    {
        const RSparseVector<double> &row = A.getVector(uint(i));
        uint begin = cache.rowPtr[uint(i)];
        uint end = cache.rowPtr[uint(i)+1];
        if (row.size() != end-begin)
        {
            samePattern = 0;
        }
        else
        {
            for (uint j=0;j<row.size();j++)
            {
                double value = row.getValue(j);
                if (cache.columnIndexes[begin+j] != row.getIndex(j))
                {
                    samePattern = 0;
                }
                cache.values[begin+j] = value;
                norm2 += value*value;
            }
        }
    }

    if (!samePattern)
    {
        return false;
    }

    cache.norm = std::sqrt(norm2);
    return true;
}

//! Return CSR cache for given matrix.
//! When allowRefresh is false an existing cache with a matching shape is
//! returned as-is — used by solveCG/solveGMRES which are only called from
//! solve() right after the cache has been refreshed.
const SparseMatrixCache &findSparseMatrixCache(const RSparseMatrix &A, bool allowRefresh = true)
{
    static std::mutex cachesMutex;
    static std::unordered_map<const RSparseMatrix*,SparseMatrixCache> caches;
    static const size_t maxCachedMatrices = 16;

    std::lock_guard<std::mutex> lock(cachesMutex);

    if (caches.find(&A) == caches.end() && caches.size() >= maxCachedMatrices)
    {
        // Bound memory usage - drop stale entries (matrices are keyed by
        // address, so destroyed matrices would otherwise accumulate forever).
        caches.clear();
    }

    SparseMatrixCache &cache = caches[&A];
    if (!allowRefresh && cache.nRows == A.getNRows() && cache.rowPtr.size() == A.getNRows()+1)
    {
        return cache;
    }
    if (!refreshSparseMatrixCache(A,cache))
    {
        cache = buildSparseMatrixCache(A);
    }
    return cache;
}

inline double multiplyRow(const SparseMatrixCache &matrix, uint row, const RRVector &x)
{
    double value = 0.0;
    for (uint j=matrix.rowPtr[row];j<matrix.rowPtr[row+1];j++)
    {
        value += matrix.values[j] * x[matrix.columnIndexes[j]];
    }
    return value;
}

}

void RMatrixSolver::_init(const RMatrixSolver *pMatrixSolver)
{
    this->iterationInfo.setConvergenceValue(this->matrixSolverConf.getSolverCvgValue());
    this->iterationInfo.setNIterations(this->matrixSolverConf.getNOuterIterations());
    this->iterationInfo.setOutputFrequency(this->matrixSolverConf.getOutputFrequency());
    this->iterationInfo.setOutputFileName(this->matrixSolverConf.getOutputFileName());

    if (pMatrixSolver)
    {
        this->matrixSolverConf = pMatrixSolver->matrixSolverConf;
        this->iterationInfo = pMatrixSolver->iterationInfo;
    }
}

RMatrixSolver::RMatrixSolver(const RMatrixSolverConf &matrixSolverConf)
    : matrixSolverConf(matrixSolverConf)
{
    this->_init();
}

RMatrixSolver::RMatrixSolver(const RMatrixSolver &matrixSolver)
{
    this->_init(&matrixSolver);
}

RMatrixSolver::~RMatrixSolver()
{
}

RMatrixSolver & RMatrixSolver::operator =(const RMatrixSolver &matrixSolver)
{
    this->_init(&matrixSolver);
    return (*this);
}

void RMatrixSolver::solve(const RSparseMatrix &A, const RRVector &b, RRVector &x, RMatrixPreconditionerType matrixPreconditionerType, unsigned int blockSize)
{
    RMatrixPreconditioner P(A,matrixPreconditionerType,blockSize);
    RRVector y(b);
    const SparseMatrixCache &matrix = findSparseMatrixCache(A);

    double An = matrix.norm;
    double bn = RRVector::euclideanNorm(b);
    double equationScale = 1.0;

    if (bn != 0.0 && An != 0.0)
    {
        equationScale = 1.0e9 / std::abs(bn/An);
    }

    RLogger::info("Unknowns = %u\n",b.size());
    RLogger::info("||A|| = %13e\n",An);
    RLogger::info("||b|| = %13e\n",bn);

    this->iterationInfo.printHeader(RMatrixSolverConf::getName(this->matrixSolverConf.getType()));

    this->iterationInfo.setEquationScale(equationScale);

    x.resize(b.getNRows(),0.0);

    y *= equationScale;
    x *= equationScale;

    switch (this->matrixSolverConf.getType())
    {
        case RMatrixSolverConf::CG:
            this->solveCG(A,y,x,P);
            break;
        case RMatrixSolverConf::GMRES:
            this->solveGMRES(A,y,x,P);
            break;
        default:
            throw RError(RError::Type::Application,R_ERROR_REF,"Unknown matrix solver type \'%d\'",this->matrixSolverConf.getType());
    }

    x *= 1.0/equationScale;

    this->iterationInfo.printFooter();
}

void RMatrixSolver::disableConvergenceLogFile()
{
    this->iterationInfo.setOutputFileName(QString());
}

void RMatrixSolver::solveCG(const RSparseMatrix &A, const RRVector &b, RRVector &x, RMatrixPreconditioner &P)
{
    unsigned int m = A.getNRows();
    const SparseMatrixCache &matrix = findSparseMatrixCache(A,false);

    RRVector r(m);
    RRVector z(m);
    RRVector p(m);
    RRVector q(m);

    r.fill(0.0);
    z.fill(0.0);
    p.fill(0.0);
    q.fill(0.0);

    double ro0 = 0.0;
    double ro1 = 0.0;
    double beta = 0.0;
    double dot = 0.0;
    double alpha = 0.0;
    double rn = 0.0;
    double bn = 0.0;

#pragma omp parallel default(shared)
    {
        // Compute initial residual.
#pragma omp for reduction(+:bn)
        for (int64_t i=0;i<int64_t(m);i++)
        {
            r[uint(i)] = b[uint(i)] - multiplyRow(matrix,uint(i),x);
            bn += b[uint(i)]*b[uint(i)];
        }
#pragma omp barrier
#pragma omp master
        {
            bn = std::sqrt(bn);
        }

        // Iterate and look for the solution
        for (unsigned int it=0;it<this->matrixSolverConf.getNOuterIterations();it++)
        {
#pragma omp barrier
#pragma omp master
            {
                this->iterationInfo.setIteration(it);
            }

#pragma omp barrier
            P.compute(r,z);

#pragma omp master
            ro0 = 0.0;
#pragma omp barrier
#pragma omp for reduction(+:ro0)
            for (int64_t i=0;i<int64_t(m);i++)
            {
                ro0 += r[uint(i)] * z[uint(i)];
            }

#pragma omp master
            {
                if (it > 0)
                {
                    if (std::abs(ro1) < RConstants::eps)
                    {
                        ro1 = (ro1 < 0.0) ? -RConstants::eps : RConstants::eps;
                    }
                    beta = ro0 / ro1;
                }
            }
#pragma omp barrier
#pragma omp for
            for (int64_t i=0;i<int64_t(m);i++)
            {
                if (it == 0)
                {
                    p[i] = z[i];
                }
                else
                {
                    p[i] = p[i]*beta + z[i];
                }
            }

            // q = A*p
#pragma omp master
            dot = 0.0;
#pragma omp barrier
#pragma omp for reduction(+:dot)
            for (int64_t i=0;i<int64_t(m);i++)
            {
                q[uint(i)] = multiplyRow(matrix,uint(i),p);
                dot += p[uint(i)]*q[uint(i)];
            }

#pragma omp master
            {
                if (dot > 0.0)
                {
                    dot = std::max(dot,RConstants::eps);
                }
                else if (dot < 0.0)
                {
                    dot = std::min(dot,-RConstants::eps);
                }
                else
                {
                    dot = RConstants::eps;
                }

                alpha = ro0 / dot;
                ro1 = ro0;
                rn = 0.0;
            }

#pragma omp barrier
#pragma omp for reduction(+:rn)
            for (int64_t i=0;i<int64_t(m);i++)
            {
                x[uint(i)] += alpha * p[uint(i)];
                r[uint(i)] -= alpha * q[uint(i)];
                rn += r[uint(i)]*r[uint(i)];
            }

#pragma omp master
            {
                rn = std::sqrt(rn);

                // A backward-error denominator (||A||*||x||+||b||) can
                // accept a large force residual in badly scaled systems.
                // Track the conventional relative residual instead.
                const double residualScale = std::max(bn,RConstants::eps);
                this->iterationInfo.setError(rn/residualScale);
                this->iterationInfo.printIteration();
            }

#pragma omp barrier
            if (this->iterationInfo.hasConverged() || this->iterationInfo.hasDiverged())
            {
                break;
            }
        }
    }
}

void RMatrixSolver::solveGMRES(const RSparseMatrix &A, const RRVector &b, RRVector &x, RMatrixPreconditioner &P)
{
    uint mA = A.getNRows();
    const SparseMatrixCache &matrix = findSparseMatrixCache(A,false);
    uint nouter = this->matrixSolverConf.getNOuterIterations();
    uint ninner = this->matrixSolverConf.getNInnerIterations();

    RRVector ro(mA,0.0);
    RRVector w(mA,0.0);
    RRVector p(ninner+1,0.0);
    RRVector y(ninner,0.0);
    RRVector c(ninner,0.0);
    RRVector s(ninner,0.0);
    RRMatrix v(ninner+1,mA,0.0);
    RRMatrix z(ninner+1,mA,0.0);
    RRMatrix h(ninner+1,ninner,0.0);

    double bn = RRVector::euclideanNorm(b);

    double beta = 0.0;
    double wNorm = 0.0;
    RRVector hValues(ninner+1,0.0);

#pragma omp parallel default(shared)
    {
        uint iti = 0;

        // Outer iteration
        for (uint ito=0;ito<nouter;ito++)
        {
#pragma omp master
            {
                this->iterationInfo.setIteration(ito);

                // Initialize matrix h, Krylov
                h.fill(0.0);

                // Compute initial residual ro and norm beta
                beta = 0.0;
            }
#pragma omp barrier
#pragma omp for reduction(+:beta)
            for (int64_t i=0;i<int64_t(mA);i++)
            {
                ro[uint(i)] = b[uint(i)] - multiplyRow(matrix,uint(i),x);
                beta += ro[uint(i)]*ro[uint(i)];
            }

#pragma omp master
            {
                beta = std::sqrt(beta);

                // Check convergence
                const double residualScale = std::max(bn,RConstants::eps);
                this->iterationInfo.setError(beta/residualScale);
                this->iterationInfo.printIteration();
            }
#pragma omp barrier
            if (this->iterationInfo.hasConverged() || this->iterationInfo.hasDiverged())
            {
                break;
            }

            // Define first krylov vector v[0]
#pragma omp for
            for (int64_t i=0;i<int64_t(mA);i++)
            {
                v[0][i] = (beta == 0.0) ? 0.0 : ro[i] / beta;
            }
#pragma omp master
            {
                // Initial rhs define p
                p.fill(0.0);
                p[0] = beta;
            }
            // Inner iteration
            for (iti=0;iti<ninner;iti++)
            {
#pragma omp barrier
                // Preconditioning
                P.compute(v[iti],z[iti]);
//                RLogger::warning("%u v %13g\n",iti,RRVector::norm(v[iti])-1.0);
//                RLogger::warning("%u z %13g\n",iti,RRVector::norm(z[iti]));
#pragma omp barrier
                // A multiplied by the last krylov vector at present
#pragma omp for
                for (int64_t i=0;i<int64_t(mA);i++)
                {
                    w[uint(i)] = multiplyRow(matrix,uint(i),z[iti]);
                }
                // Orthogonalize the candidate vector against the current basis.
#pragma omp single
                hValues.fill(0.0);
#pragma omp barrier
                std::vector<double> hLocal(iti+1,0.0);
#pragma omp for nowait
                for (int64_t j=0;j<int64_t(mA);j++)
                {
                    double wj = w[uint(j)];
                    for (uint i=0;i<=iti;i++)
                    {
                        hLocal[i] += wj * v[i][uint(j)];
                    }
                }
#pragma omp critical
                {
                    for (uint i=0;i<=iti;i++)
                    {
                        hValues[i] += hLocal[i];
                    }
                }
#pragma omp barrier
#pragma omp single
                {
                    for (uint i=0;i<=iti;i++)
                    {
                        h[i][iti] = hValues[i];
                    }
                }
#pragma omp barrier
#pragma omp for
                for (int64_t j=0;j<int64_t(mA);j++)
                {
                    double value = w[uint(j)];
                    for (uint i=0;i<=iti;i++)
                    {
                        value -= hValues[i] * v[i][uint(j)];
                    }
                    w[uint(j)] = value;
                }

#pragma omp single
                wNorm = 0.0;
#pragma omp barrier
#pragma omp for reduction(+:wNorm)
                for (int64_t j=0;j<int64_t(mA);j++)
                {
                    wNorm += w[uint(j)]*w[uint(j)];
                }
#pragma omp single
                {
                    // Finding norm of the krylov vector
                    h[iti+1][iti] = std::sqrt(wNorm);
//                    RLogger::warning("%u w %13g\n",iti,RRVector::norm(w));
                }
#pragma omp barrier
                // New krylov vector formed
#pragma omp for
                for (int64_t i=0;i<int64_t(mA);i++)
                {
                    v[iti+1][i] = (h[iti+1][iti] == 0.0) ? 0.0 : w[i] / h[iti+1][iti];
                }
#pragma omp barrier
#pragma omp master
                {
                    for (uint i=1;i<=iti;i++)
                    {
                        double hsave = h[i-1][iti];
                        h[i-1][iti] =  c[i-1]*hsave + s[i-1]*h[i][iti];
                        h[i]  [iti] = -s[i-1]*hsave + c[i-1]*h[i][iti];
                    }

                    double g = std::sqrt(  h[iti][iti]   * h[iti][iti]
                                         + h[iti+1][iti] * h[iti+1][iti]);
//                    RLogger::warning("%u g %13g\n",iti,g);
                    if (g == 0.0)
                    {
                        c[iti] = s[iti] = 0.0;
                    }
                    else
                    {
                        c[iti] = h[iti][iti]/g;
                        s[iti] = h[iti+1][iti]/g;
                    }
                    h[iti]  [iti] = g;
                    h[iti+1][iti] = 0.0;
                    p[iti+1]      = -s[iti]*p[iti];
                    p[iti]        =  c[iti]*p[iti];
                }
#pragma omp barrier
                if (fabs(p[iti+1]) < beta*1.0e-4)
                {
                    iti++;
                    break;
                }
            } // inner iteration

#pragma omp barrier
            // Backward substitution
#pragma omp master
            {
                y[iti-1] = (h[iti-1][iti-1] == 0.0) ? 0.0 : p[iti-1] / h[iti-1][iti-1];
                for (int i=iti-2;i>=0;i--)
                {
                    for (int j=iti-1;j>i;j--)
                    {
                        p[i] -= y[j] * h[i][j];
                    }
                    y[i] = (h[i][i] == 0.0) ? 0.0 : p[i] / h[i][i];
//                    RLogger::warning("%u %u y %13g\n",iti,i,y[i]);
                }
            }

#pragma omp barrier
            // Update prod
#pragma omp for
            for (int64_t i=0;i<int64_t(mA);i++)
            {
                for (uint j=0;j<=iti-1;j++)
                {
                    x[i] += z[j][i] * y[j];
                }
            }
        } // outer iteration
    }
}

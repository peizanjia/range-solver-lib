#include <cmath>
#include <limits>

#include "reigenvaluesolver.h"

#define R_EIS_PYTHAG(a,b) (std::sqrt ((a)*(a) + (b)*(b)))
#define R_ASSERT(_condition) { if (!(_condition)) { RLogger::unindent(); R_ERROR_ASSERT(_condition); } }

void REigenValueSolver::_init(const REigenValueSolver *pEigenValueSolver)
{
    if (pEigenValueSolver)
    {
        this->eigenValueSolverConf = pEigenValueSolver->eigenValueSolverConf;
        this->matrixSolverConf = pEigenValueSolver->matrixSolverConf;
    }
}

REigenValueSolver::REigenValueSolver(const REigenValueSolverConf &eigenValueSolverConf, const RMatrixSolverConf &matrixSolverConf)
    : eigenValueSolverConf(eigenValueSolverConf)
    , matrixSolverConf(matrixSolverConf)
{
    this->_init();
}

REigenValueSolver::REigenValueSolver(const REigenValueSolver &eigenValueSolver)
{
    this->_init(&eigenValueSolver);
}

REigenValueSolver::~REigenValueSolver()
{

}

REigenValueSolver &REigenValueSolver::operator =(const REigenValueSolver &eigenValueSolver)
{
    this->_init(&eigenValueSolver);
    return (*this);
}

void REigenValueSolver::solve(const RSparseMatrix &M, const RSparseMatrix &K, RRVector &d, RRMatrix &ev)
{
    switch (this->eigenValueSolverConf.method)
    {
        case REigenValueSolverConf::Lanczos:
        {
            // Find multiple eigen values.
            try
            {
                this->solveLanczos(M,K,d);
            }
            catch (const RError &error)
            {
                throw RError(RError::Type::Application,R_ERROR_REF,"Lanczos method failed. %s",error.getMessage().toUtf8().constData());
            }
            break;
        }
        case REigenValueSolverConf::SubspaceIteration:
        {
            // Find several of the lowest eigen values.
            try
            {
                this->solveSubspaceIteration(M,K,d,ev);
            }
            catch (const RError &error)
            {
                throw RError(RError::Type::Application,R_ERROR_REF,"Subspace iteration failed. %s",error.getMessage().toUtf8().constData());
            }
            break;
        }
        case REigenValueSolverConf::InversePowerIteration:
        {
            // Find the single lowest eigen value.
            try
            {
                this->solveInversePowerIteration(M,K,d,ev);
            }
            catch (const RError &error)
            {
                throw RError(RError::Type::Application,R_ERROR_REF,"Inverse power iteration failed. %s",error.getMessage().toUtf8().constData());
            }
            break;
        }
        default:
        {
            throw RError(RError::Type::Application,R_ERROR_REF,"Unknown eigen-value solver method.");
        }
    }

    // Every method returns lambda of K * phi = lambda * M * phi directly. Sort
    // the values ascending so that d[0] is always the lowest one, and carry the
    // eigen vectors along.
    try
    {
        if (d.getNRows() > 1)
        {
            std::vector<uint> indexes;
            RUtil::qSort(d,indexes);

            if (ev.getNRows() == d.getNRows())
            {
                // indexes[i] holds the original row of the i-th sorted value.
                RRMatrix evOld(ev);
                for (uint i=0;i<ev.getNRows();i++)
                {
                    for (uint j=0;j<ev.getNColumns();j++)
                    {
                        ev[i][j] = evOld[indexes[i]][j];
                    }
                }
            }
        }
    }
    catch (const RError &error)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,"Failed to reorder eigen values and vectors. %s",error.getMessage().toUtf8().constData());
    }
}

void REigenValueSolver::solveLanczos(const RSparseMatrix &M, const RSparseMatrix &K, RRVector &d)
{
    RRVector e;
    RRVector v(M.getNRows(),0.0);
    RRVector vo(M.getNRows(),0.0);
    RRVector w(M.getNRows(),0.0);

    for (uint i=0;i<M.getNRows();i++)
    {
        v[i] = double(std::rand()) / double(RAND_MAX);
    }
    v.normalize();

    RMatrixSolver solver(this->matrixSolverConf);

    uint ne = std::min(this->eigenValueSolverConf.getNEigenValues(),K.getNRows());

    d.resize(ne,0.0);
    e.resize(ne,0.0);

    // Lanczos iteration
    for (uint i=0;i<ne;i++)
    {
        RLogger::info("Lanczos iteration %u of %u\n",i+1,ne);
        RLogger::indent();

        // w = A*v - e(i)*vo
        // w = (M/K)*v - e(i)*vo
        // K*w = M*v - K*e(i)*vo

        RRVector b;

        RSparseMatrix::mlt(M,v,b);

        if (i > 0)
        {
            RRVector bTmp;
            RSparseMatrix::mlt(K,vo,bTmp);
            bTmp *= e[i];
            for (uint j=0;j<b.size();j++)
            {
                b[j] -= bTmp[j];
            }
        }

        try
        {
            solver.solve(K,b,w,R_MATRIX_PRECONDITIONER_JACOBI);
        }
        catch (const RError &error)
        {
            RLogger::unindent();
            throw RError(RError::Type::Application,R_ERROR_REF,"Failed to solve matrix system. %s", error.getMessage().toUtf8().constData());
        }

        // d(i) = w dot v
        d[i] = RRVector::dot(w,v);

        if (i+1 < ne)
        {
            // w = w - d(i)*v
            for (uint j=0;j<w.size();j++)
            {
                w[j] -= d[i] * v[j];
            }

            // e(i+1) = || w ||
            e[i+1] = RRVector::euclideanNorm(w);
            R_ASSERT(e[i+1] != 0.0);

            // vo = v
            vo = v;

            // v = w / e(i+1)
            for (uint j=0;j<v.size();j++)
            {
                v[j] = w[j] / e[i+1];
            }
        }
        RLogger::unindent();
    }

    try
    {
        REigenValueSolver::qlDecomposition(d,e);
    }
    catch (const RError &error)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,"QL decomposition failed. %s",error.getMessage().toUtf8().constData());
    }

    // The iteration runs on K^-1 * M, so its values are the reciprocals of the
    // eigen values of K * phi = lambda * M * phi. A value which stayed at zero
    // was not resolved and stands for an infinite eigen value.
    for (uint i=0;i<d.getNRows();i++)
    {
        double value = std::fabs(d[i]);
        d[i] = (value < RConstants::eps)
             ? std::numeric_limits<double>::infinity()
             : 1.0/value;
    }
}


void REigenValueSolver::solveSubspaceIteration(const RSparseMatrix &M, const RSparseMatrix &K, RRVector &d, RRMatrix &ev)
{
    uint n = K.getNRows();
    uint nEigen = std::min(this->eigenValueSolverConf.getNEigenValues(),n);

    if (nEigen == 0)
    {
        d.resize(0);
        ev.resize(0,0);
        return;
    }

    // A few extra vectors make the wanted eigen pairs converge markedly faster,
    // but every vector of the subspace costs a linear solve in every iteration.
    // q = min(2p,p+8) is the usual compromise between the two.
    uint m = std::min(n,std::min(2*nEigen,nEigen+8));

    RRMatrix X(n,m,0.0);
    RRMatrix Y(n,m,0.0);

    // Start from a set of unit vectors spread over the system plus a random
    // perturbation, so that the block is never accidentally deficient and the
    // result does not depend on the random sequence alone.
    for (uint j=0;j<m;j++)
    {
        for (uint i=0;i<n;i++)
        {
            X[i][j] = 1.0e-3 * (double(std::rand()) / double(RAND_MAX) - 0.5);
        }
        X[(j*n)/m][j] += 1.0;
    }

    RMatrixSolver solver(this->matrixSolverConf);

    RRVector theta(m,0.0);
    RRVector thetaOld(m,0.0);

    uint nIterations = std::max(this->eigenValueSolverConf.getNIterations(),uint(1));
    double convergenceValue = this->eigenValueSolverConf.getSolverCvgValue();

    RRVector x(n,0.0);
    RRVector y(n,0.0);
    RRVector b(n,0.0);

    // State the cost before starting - it is set by the number of requested
    // eigen values and is the first thing to reduce when a run takes too long.
    RLogger::info("Extracting %u eigen values from a subspace of %u vectors - up to %u iterations of %u linear solves each.\n",
                  nEigen,m,nIterations,m);

    for (uint it=0;it<nIterations;it++)
    {
        RLogger::info("Subspace iteration %u of %u\n",it+1,nIterations);
        RLogger::indent();

        // Y = K^-1 * M * X, one column at a time. Each column is a full linear
        // solve, so the progress is reported - a subspace iteration of a large
        // model takes a long time and is otherwise silent.
        for (uint j=0;j<m;j++)
        {
            RLogger::info("Subspace vector %u of %u\n",j+1,m);
            for (uint i=0;i<n;i++)
            {
                x[i] = X[i][j];
                y[i] = 0.0;
            }
            RSparseMatrix::mlt(M,x,b);
            try
            {
                solver.solve(K,b,y,R_MATRIX_PRECONDITIONER_JACOBI);
            }
            catch (const RError &error)
            {
                RLogger::unindent();
                throw RError(RError::Type::Application,R_ERROR_REF,"Failed to solve matrix system. %s",error.getMessage().toUtf8().constData());
            }
            for (uint i=0;i<n;i++)
            {
                Y[i][j] = y[i];
            }
        }

        uint mKept = REigenValueSolver::orthonormalizeColumns(Y);
        if (mKept == 0)
        {
            RLogger::unindent();
            throw RError(RError::Type::Application,R_ERROR_REF,"Subspace collapsed - no independent direction is left.");
        }

        // Rayleigh-Ritz projection onto the subspace spanned by Y.
        RRMatrix Kr(mKept,mKept,0.0);
        RRMatrix Mr(mKept,mKept,0.0);

        for (uint j=0;j<mKept;j++)
        {
            for (uint i=0;i<n;i++)
            {
                x[i] = Y[i][j];
            }

            RRVector kx,mx;
            RSparseMatrix::mlt(K,x,kx);
            RSparseMatrix::mlt(M,x,mx);

            for (uint l=0;l<mKept;l++)
            {
                double kSum = 0.0;
                double mSum = 0.0;
                for (uint i=0;i<n;i++)
                {
                    kSum += Y[i][l] * kx[i];
                    mSum += Y[i][l] * mx[i];
                }
                Kr[l][j] = kSum;
                Mr[l][j] = mSum;
            }
        }

        // Reduce K_r*z = theta*M_r*z to a standard symmetric problem using a
        // Cholesky decomposition of M_r.
        RRMatrix L;
        if (!REigenValueSolver::choleskyDecomposition(Mr,L))
        {
            RLogger::unindent();
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Projected mass matrix is not positive definite - check that every computable entity has a density.");
        }

        RRMatrix G(mKept,mKept,0.0);
        RRVector column(mKept,0.0);
        RRVector solution(mKept,0.0);

        // G = L^-1 * K_r
        for (uint j=0;j<mKept;j++)
        {
            for (uint i=0;i<mKept;i++)
            {
                column[i] = Kr[i][j];
            }
            REigenValueSolver::forwardSubstitution(L,column,solution);
            for (uint i=0;i<mKept;i++)
            {
                G[i][j] = solution[i];
            }
        }

        // A = L^-1 * G^T, symmetric.
        RRMatrix A(mKept,mKept,0.0);
        for (uint j=0;j<mKept;j++)
        {
            for (uint i=0;i<mKept;i++)
            {
                column[i] = G[j][i];
            }
            REigenValueSolver::forwardSubstitution(L,column,solution);
            for (uint i=0;i<mKept;i++)
            {
                A[i][j] = solution[i];
            }
        }

        RRVector thetaRitz;
        RRMatrix Z;
        REigenValueSolver::jacobiEigen(A,thetaRitz,Z);

        // Order the Ritz pairs by ascending eigen value.
        std::vector<uint> order(mKept);
        for (uint i=0;i<mKept;i++)
        {
            order[i] = i;
        }
        std::sort(order.begin(),order.end(),[&thetaRitz](uint a, uint b) { return thetaRitz[a] < thetaRitz[b]; });

        // Back transform the Ritz vectors and rebuild the subspace.
        X.resize(n,mKept,0.0);
        for (uint j=0;j<mKept;j++)
        {
            for (uint i=0;i<mKept;i++)
            {
                column[i] = Z[i][order[j]];
            }
            REigenValueSolver::backwardSubstitution(L,column,solution);

            for (uint i=0;i<n;i++)
            {
                double sum = 0.0;
                for (uint l=0;l<mKept;l++)
                {
                    sum += Y[i][l] * solution[l];
                }
                X[i][j] = sum;
            }
        }

        theta.resize(mKept,0.0);
        for (uint j=0;j<mKept;j++)
        {
            theta[j] = thetaRitz[order[j]];
        }
        m = mKept;
        Y.resize(n,mKept,0.0);

        // Ritz values may stabilize before their vectors do. Require both the
        // requested eigenvalue change and a generalized residual for every
        // wanted pair; modal damping depends on the vectors themselves.
        bool converged = (it > 0);
        uint nWanted = std::min(nEigen,mKept);
        double maxChange = 0.0;
        double maxResidual = 0.0;
        const double residualTolerance = std::max(1.0e-8,convergenceValue);
        for (uint j=0;j<nWanted;j++)
        {
            double scale = std::max(std::fabs(theta[j]),RConstants::eps);
            double change = std::fabs(theta[j]-thetaOld[j]) / scale;
            maxChange = std::max(maxChange,change);
            if (change > convergenceValue)
            {
                converged = false;
            }
            for (uint i=0;i<n;i++) x[i] = X[i][j];
            RRVector kx,mx;
            RSparseMatrix::mlt(K,x,kx);
            RSparseMatrix::mlt(M,x,mx);
            RRVector residual(kx);
            for (uint i=0;i<n;i++) residual[i] -= theta[j]*mx[i];
            const double denominator = RRVector::euclideanNorm(kx) +
                std::fabs(theta[j])*RRVector::euclideanNorm(mx);
            const double relativeResidual = denominator > 0.0
                ? RRVector::euclideanNorm(residual)/denominator
                : std::numeric_limits<double>::infinity();
            if (!std::isfinite(relativeResidual) || relativeResidual > residualTolerance)
            {
                converged = false;
            }
            maxResidual = std::max(maxResidual,relativeResidual);
        }
        if (it > 0)
        {
            RLogger::info("Convergence rate = %g (required %g); maximum generalized residual = %g (required %g)\n",
                          maxChange,convergenceValue,maxResidual,residualTolerance);
        }

        thetaOld = theta;

        RLogger::unindent();

        if (converged)
        {
            break;
        }
        if (it + 1 == nIterations)
        {
            RLogger::warning("Subspace iteration did not converge in %u iterations - eigenvalue change %g (required %g), "
                             "generalized residual %g (required %g). Extracting fewer modes or raising the iteration limit may be needed.\n",
                             nIterations,maxChange,convergenceValue,maxResidual,residualTolerance);
        }
    }

    uint nFound = std::min(nEigen,uint(theta.size()));

    if (nFound < nEigen)
    {
        RLogger::warning("Only %u of the %u requested eigen values could be extracted - the subspace lost independent directions.\n",nFound,nEigen);
    }

    d.resize(nFound,0.0);
    ev.resize(nFound,n,0.0);

    for (uint j=0;j<nFound;j++)
    {
        d[j] = theta[j];
        for (uint i=0;i<n;i++)
        {
            ev[j][i] = X[i][j];
        }
    }
}


void REigenValueSolver::solveInversePowerIteration(const RSparseMatrix &M, const RSparseMatrix &K, RRVector &d, RRMatrix &ev)
{
    uint n = K.getNRows();

    RRVector x(n,0.0);
    for (uint i=0;i<n;i++)
    {
        x[i] = double(std::rand()) / double(RAND_MAX) - 0.5;
    }
    if (x.normalize() == 0.0)
    {
        throw RError(RError::Type::Application,R_ERROR_REF,"Failed to generate a start vector.");
    }

    RMatrixSolver solver(this->matrixSolverConf);

    uint nIterations = std::max(this->eigenValueSolverConf.getNIterations(),uint(1));
    double convergenceValue = this->eigenValueSolverConf.getSolverCvgValue();

    double lambda = 0.0;
    RRVector b(n,0.0);
    RRVector y(n,0.0);

    for (uint it=0;it<nIterations;it++)
    {
        RLogger::info("Inverse power iteration %u of %u\n",it+1,nIterations);
        RLogger::indent();

        // Solving K*y = M*x drives x towards the eigen vector of the lowest
        // eigen value, which is the dominant one of K^-1 * M.
        RSparseMatrix::mlt(M,x,b);
        y.fill(0.0);
        try
        {
            solver.solve(K,b,y,R_MATRIX_PRECONDITIONER_JACOBI);
        }
        catch (const RError &error)
        {
            RLogger::unindent();
            throw RError(RError::Type::Application,R_ERROR_REF,"Failed to solve matrix system. %s",error.getMessage().toUtf8().constData());
        }

        if (y.normalize() == 0.0)
        {
            RLogger::unindent();
            throw RError(RError::Type::Application,R_ERROR_REF,"Inverse power iteration collapsed to a zero vector.");
        }
        x = y;

        RRVector kx,mx;
        RSparseMatrix::mlt(K,x,kx);
        RSparseMatrix::mlt(M,x,mx);

        double numerator = RRVector::dot(x,kx);
        double denominator = RRVector::dot(x,mx);

        if (std::fabs(denominator) < RConstants::eps)
        {
            RLogger::unindent();
            throw RError(RError::Type::Application,R_ERROR_REF,
                         "Projected mass is zero - check that every computable entity has a density.");
        }

        double lambdaNew = numerator / denominator;
        double change = std::fabs(lambdaNew-lambda) / std::max(std::fabs(lambdaNew),RConstants::eps);

        RLogger::info("Eigen value = %g, convergence rate = %g\n",lambdaNew,change);
        RLogger::unindent();

        lambda = lambdaNew;

        if (it > 0 && change < convergenceValue)
        {
            break;
        }
    }

    d.resize(1,0.0);
    d[0] = lambda;

    ev.resize(1,n,0.0);
    for (uint i=0;i<n;i++)
    {
        ev[0][i] = x[i];
    }
}


uint REigenValueSolver::orthonormalizeColumns(RRMatrix &X)
{
    uint n = X.getNRows();
    uint m = X.getNColumns();
    uint kept = 0;

    // A column is dependent when the orthogonalization leaves almost nothing of
    // it. What is left has to be measured against the length the column started
    // with: the columns arrive as K^-1 * M * x, so their length carries the
    // magnitude of the matrices and says nothing about linear independence.
    const double dependencyTolerance = 1.0e-8;

    for (uint j=0;j<m;j++)
    {
        double norm = 0.0;
        for (uint i=0;i<n;i++)
        {
            norm += X[i][j] * X[i][j];
        }
        norm = std::sqrt(norm);

        if (norm == 0.0)
        {
            // Nothing to normalize - the column is exactly zero.
            continue;
        }

        for (uint i=0;i<n;i++)
        {
            X[i][j] /= norm;
        }

        // Modified Gram-Schmidt against the columns already kept, run twice.
        // A single pass loses orthogonality when the columns are nearly
        // parallel, which is what subspace iteration drives them towards.
        for (uint pass=0;pass<2;pass++)
        {
            for (uint l=0;l<kept;l++)
            {
                double dot = 0.0;
                for (uint i=0;i<n;i++)
                {
                    dot += X[i][l] * X[i][j];
                }
                for (uint i=0;i<n;i++)
                {
                    X[i][j] -= dot * X[i][l];
                }
            }
        }

        // The column is a unit vector, so what is left of it after the
        // projection is the sine of the angle to the space already spanned.
        norm = 0.0;
        for (uint i=0;i<n;i++)
        {
            norm += X[i][j] * X[i][j];
        }
        norm = std::sqrt(norm);

        if (norm < dependencyTolerance)
        {
            // Dependent direction - drop it.
            continue;
        }

        for (uint i=0;i<n;i++)
        {
            X[i][kept] = X[i][j] / norm;
        }
        kept++;
    }

    return kept;
}


bool REigenValueSolver::choleskyDecomposition(const RRMatrix &A, RRMatrix &L)
{
    uint n = A.getNRows();

    L.resize(n,n,0.0);
    L.fill(0.0);

    for (uint i=0;i<n;i++)
    {
        for (uint j=0;j<=i;j++)
        {
            double sum = A[i][j];
            for (uint k=0;k<j;k++)
            {
                sum -= L[i][k] * L[j][k];
            }

            if (i == j)
            {
                if (sum <= 0.0)
                {
                    return false;
                }
                L[i][i] = std::sqrt(sum);
            }
            else
            {
                L[i][j] = sum / L[j][j];
            }
        }
    }

    return true;
}


void REigenValueSolver::forwardSubstitution(const RRMatrix &L, const RRVector &b, RRVector &x)
{
    uint n = L.getNRows();

    x.resize(n,0.0);

    for (uint i=0;i<n;i++)
    {
        double sum = b[i];
        for (uint j=0;j<i;j++)
        {
            sum -= L[i][j] * x[j];
        }
        x[i] = sum / L[i][i];
    }
}


void REigenValueSolver::backwardSubstitution(const RRMatrix &L, const RRVector &b, RRVector &x)
{
    uint n = L.getNRows();

    x.resize(n,0.0);

    for (uint i=n;i>0;i--)
    {
        uint row = i-1;
        double sum = b[row];
        for (uint j=row+1;j<n;j++)
        {
            sum -= L[j][row] * x[j];
        }
        x[row] = sum / L[row][row];
    }
}


void REigenValueSolver::jacobiEigen(RRMatrix &A, RRVector &d, RRMatrix &V)
{
    uint n = A.getNRows();

    V.resize(n,n,0.0);
    V.setIdentity(n);
    d.resize(n,0.0);

    double scale = 0.0;
    for (uint i=0;i<n;i++)
    {
        for (uint j=0;j<n;j++)
        {
            scale += A[i][j] * A[i][j];
        }
    }
    double threshold = std::max(scale,1.0) * 1.0e-30;

    for (uint sweep=0;sweep<100;sweep++)
    {
        double off = 0.0;
        for (uint p=0;p<n;p++)
        {
            for (uint q=p+1;q<n;q++)
            {
                off += A[p][q] * A[p][q];
            }
        }
        if (off <= threshold)
        {
            break;
        }

        for (uint p=0;p<n;p++)
        {
            for (uint q=p+1;q<n;q++)
            {
                if (std::fabs(A[p][q]) <= 0.0)
                {
                    continue;
                }

                double theta = (A[q][q] - A[p][p]) / (2.0 * A[p][q]);
                double t = ((theta >= 0.0) ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta*theta + 1.0));
                double c = 1.0 / std::sqrt(t*t + 1.0);
                double s = t * c;

                for (uint k=0;k<n;k++)
                {
                    double akp = A[k][p];
                    double akq = A[k][q];
                    A[k][p] = c*akp - s*akq;
                    A[k][q] = s*akp + c*akq;
                }
                for (uint k=0;k<n;k++)
                {
                    double apk = A[p][k];
                    double aqk = A[q][k];
                    A[p][k] = c*apk - s*aqk;
                    A[q][k] = s*apk + c*aqk;
                }
                for (uint k=0;k<n;k++)
                {
                    double vkp = V[k][p];
                    double vkq = V[k][q];
                    V[k][p] = c*vkp - s*vkq;
                    V[k][q] = s*vkp + c*vkq;
                }
            }
        }
    }

    for (uint i=0;i<n;i++)
    {
        d[i] = A[i][i];
    }
}


void REigenValueSolver::qlDecomposition(RRVector &d, RRVector &e)
{
    uint n = d.getNRows();

    // Convenient to renumber the elements of e
    for (uint i=1;i<n;i++)
    {
        e[i-1] = e[i];
    }
    e[n-1] = 0.0;

    for (uint l=0;l<n;l++)
    {
        RLogger::info("QL decomposition %u of %u\n",l+1,n);
        RLogger::indent();

        uint iter = 0;
        uint m = 0;
        while (true)
        {
            // Look for a single small subdiagonal element to split the matrix
            for (m=l;m<n-1;m++)
            {
                double dd = std::fabs(d[m]) + std::fabs(d[m+1]);
                if (std::fabs(e[m]+dd) == dd)
                {
                    break;
                }
            }
            if (m == l)
            {
                break;
            }
            if (iter ++ >= 30)
            {
                RLogger::warning("Too many iterations %u.\n",iter);
            }
            // Form shift
            double g = (d[l+1] - d[l]) / (2.0*e[l]);
            double r = R_EIS_PYTHAG(g,1.0);
            // This is dm - ks
            g = d[m] - d[l] + e[l] / (g+R_SAME_SIGN(r,g));

            double s = 1.0, c = 1.0, p = 0.0;
            // A plane rotation as in the original QL, followed by Givens rotations to restore tridiagonal form
            int gi;
            for (gi=int(m)-1;gi>=int(l);gi--)
            {
                double f = s * e[gi];
                double b = c * e[gi];
                e[gi+1] = (r = R_EIS_PYTHAG(f,g));

                // Recover from underflow
                if (r == 0.0)
                {
                    d[gi+1] -= p;
                    e[m] = 0.0;
                    break;
                }

                s = f/r;
                c = g/r;
                g = d[gi+1] - p;
                r = (d[gi] - g) * s + 2.0 * c * b;
                d[gi+1] = g + (p = s * r);
                g = c * r - b;
            }
            if (r == 0.0 && gi >= int(l))
            {
                continue;
            }
            d[l] -= p;
            e[l] = g;
            e[m] = 0.0;
        }
        RLogger::unindent();
    }
}


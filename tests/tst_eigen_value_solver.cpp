#include <QtTest>

#include <cmath>
#include <cstdlib>
#include <limits>

#include "reigenvaluesolver.h"
#include "rml_matrix_solver_conf.h"
#include "rml_sparse_matrix.h"

namespace
{

RMatrixSolverConf makeCgConf()
{
    RMatrixSolverConf conf(RMatrixSolverConf::CG);
    conf.setNInnerIterations(200);
    conf.setNOuterIterations(2000);
    conf.setSolverCvgValue(1.0e-14);
    conf.setOutputFrequency(0);
    return conf;
}

REigenValueSolverConf makeEigenConf(REigenValueSolverConf::Method method, uint nEigenValues)
{
    REigenValueSolverConf conf(method);
    conf.setNEigenValues(nEigenValues);
    conf.setNIterations(500);
    conf.setSolverCvgValue(1.0e-14);
    conf.setOutputFrequency(0);
    return conf;
}

//! K = [[2,-1],[-1,1]] * 1e6, M = identity.
//!
//! The generalised eigenvalues lambda are around 3.8e5 and 2.6e6, so lambda and
//! its reciprocal are more than ten orders of magnitude apart. That makes it
//! unambiguous which of the two a returned value is on, without relying on the
//! accuracy of the iteration itself.
void buildSystem(RSparseMatrix &m, RSparseMatrix &k)
{
    const double scale = 1.0e6;

    m.clear();
    m.setNRows(2);
    m.addValue(0,0,1.0);
    m.addValue(1,1,1.0);

    k.clear();
    k.setNRows(2);
    k.addValue(0,0,2.0*scale);
    k.addValue(0,1,-1.0*scale);
    k.addValue(1,0,-1.0*scale);
    k.addValue(1,1,1.0*scale);
}

void buildUnderiteratedSystem(RSparseMatrix &m, RSparseMatrix &k)
{
    m.clear();
    m.setNRows(3);
    m.addValue(0,0,1.0);
    m.addValue(1,1,1.5);
    m.addValue(2,2,0.75);

    k.clear();
    k.setNRows(3);
    k.addValue(0,0,4.0);
    k.addValue(0,1,-1.0);
    k.addValue(1,0,-1.0);
    k.addValue(0,2,0.25);
    k.addValue(2,0,0.25);
    k.addValue(1,1,3.0);
    k.addValue(1,2,-0.5);
    k.addValue(2,1,-0.5);
    k.addValue(2,2,2.0);
}

bool hasAccurateGeneralizedEigenpairs(const RSparseMatrix &m,
                                      const RSparseMatrix &k,
                                      const RRVector &eigenvalues,
                                      const RRMatrix &eigenvectors)
{
    if (eigenvalues.getNRows() != eigenvectors.getNRows())
    {
        return false;
    }

    for (uint mode=0;mode<eigenvalues.getNRows();mode++)
    {
        const double lambda = eigenvalues[mode];
        if (!std::isfinite(lambda))
        {
            return false;
        }

        RRVector phi(eigenvectors.getNColumns(),0.0);
        for (uint i=0;i<eigenvectors.getNColumns();i++)
        {
            phi[i] = eigenvectors[mode][i];
        }

        RRVector kPhi,mPhi;
        RSparseMatrix::mlt(k,phi,kPhi);
        RSparseMatrix::mlt(m,phi,mPhi);

        double kSquaredNorm = 0.0;
        double mSquaredNorm = 0.0;
        double residualSquaredNorm = 0.0;
        for (uint i=0;i<phi.getNRows();i++)
        {
            const double residual = kPhi[i] - lambda*mPhi[i];
            if (!std::isfinite(kPhi[i]) || !std::isfinite(mPhi[i]) || !std::isfinite(residual))
            {
                return false;
            }
            kSquaredNorm += kPhi[i]*kPhi[i];
            mSquaredNorm += mPhi[i]*mPhi[i];
            residualSquaredNorm += residual*residual;
        }

        const double denominator = std::sqrt(kSquaredNorm) +
                                   std::fabs(lambda)*std::sqrt(mSquaredNorm);
        const double relativeResidual = denominator > 0.0
                                      ? std::sqrt(residualSquaredNorm)/denominator
                                      : std::numeric_limits<double>::infinity();
        if (!std::isfinite(relativeResidual) || relativeResidual > 1.0e-8)
        {
            return false;
        }
    }

    return true;
}

}

class TestEigenValueSolver : public QObject
{
    Q_OBJECT

private slots:

    void multipleValuesAreGeneralisedEigenValues();
    void singleValueUsesTheSameScale();
    void inversePowerIterationFindsTheLowestValue();
    void inversePowerIterationAcceptsExactZeroTolerance();
    void valuesAreSortedAscending();
    void subspaceIterationFailsWhenUnderiterated();
    void inversePowerIterationFailsWhenUnderiterated();
};

//! Several extracted values have to come back as lambda of
//! K * phi = lambda * M * phi, not as its reciprocal.
void TestEigenValueSolver::multipleValuesAreGeneralisedEigenValues()
{
    RSparseMatrix m,k;
    buildSystem(m,k);

    RRVector d;
    RRMatrix ev;

    REigenValueSolver solver(makeEigenConf(REigenValueSolverConf::SubspaceIteration,2),makeCgConf());
    solver.solve(m,k,d,ev);

    QCOMPARE(d.getNRows(),2u);

    double lambdaLow = (3.0 - std::sqrt(5.0)) / 2.0 * 1.0e6;
    double lambdaHigh = (3.0 + std::sqrt(5.0)) / 2.0 * 1.0e6;

    QVERIFY2(std::fabs(d[0]-lambdaLow) < 1.0e-6*lambdaLow,
             qPrintable(QString("Lowest eigen value = %1, expected %2").arg(d[0]).arg(lambdaLow)));
    QVERIFY2(std::fabs(d[1]-lambdaHigh) < 1.0e-6*lambdaHigh,
             qPrintable(QString("Highest eigen value = %1, expected %2").arg(d[1]).arg(lambdaHigh)));
    QVERIFY2(hasAccurateGeneralizedEigenpairs(m,k,d,ev),"Generalized residual exceeded 1e-8.");
}

//! Extracting a single value used to skip the inversion, handing back a value
//! on the reciprocal scale. One and several values have to agree in scale.
void TestEigenValueSolver::singleValueUsesTheSameScale()
{
    RSparseMatrix m,k;
    buildSystem(m,k);

    RRVector d;
    RRMatrix ev;

    REigenValueSolver solver(makeEigenConf(REigenValueSolverConf::SubspaceIteration,1),makeCgConf());
    solver.solve(m,k,d,ev);

    QCOMPARE(d.getNRows(),1u);

    double lambdaLow = (3.0 - std::sqrt(5.0)) / 2.0 * 1.0e6;
    QVERIFY2(std::fabs(d[0]-lambdaLow) < 1.0e-6*lambdaLow,
             qPrintable(QString("Single eigen value = %1, expected %2").arg(d[0]).arg(lambdaLow)));
    QVERIFY2(hasAccurateGeneralizedEigenpairs(m,k,d,ev),"Generalized residual exceeded 1e-8.");
}

//! The inverse power iteration has to converge to the same lowest eigen value.
void TestEigenValueSolver::inversePowerIterationFindsTheLowestValue()
{
    RSparseMatrix m,k;
    buildSystem(m,k);

    RRVector d;
    RRMatrix ev;

    REigenValueSolver solver(makeEigenConf(REigenValueSolverConf::InversePowerIteration,1),makeCgConf());
    solver.solve(m,k,d,ev);

    QCOMPARE(d.getNRows(),1u);

    double lambdaLow = (3.0 - std::sqrt(5.0)) / 2.0 * 1.0e6;
    QVERIFY2(std::fabs(d[0]-lambdaLow) < 1.0e-6*lambdaLow,
             qPrintable(QString("Eigen value = %1, expected %2").arg(d[0]).arg(lambdaLow)));

    QCOMPARE(ev.getNRows(),1u);
    QCOMPARE(ev.getNColumns(),2u);
    QVERIFY2(hasAccurateGeneralizedEigenpairs(m,k,d,ev),"Generalized residual exceeded 1e-8.");
}

//! The lowest eigen value has to come first, so that mode 0 is the fundamental.
void TestEigenValueSolver::valuesAreSortedAscending()
{
    RSparseMatrix m,k;
    buildSystem(m,k);

    RRVector d;
    RRMatrix ev;

    REigenValueSolver solver(makeEigenConf(REigenValueSolverConf::SubspaceIteration,2),makeCgConf());
    solver.solve(m,k,d,ev);

    QCOMPARE(d.getNRows(),2u);
    QVERIFY2(d[0] <= d[1],
             qPrintable(QString("Eigen values %1, %2 are not sorted ascending").arg(d[0]).arg(d[1])));
    QVERIFY2(hasAccurateGeneralizedEigenpairs(m,k,d,ev),"Generalized residual exceeded 1e-8.");
}

//! A scalar eigenpair is exact, so a zero convergence tolerance can be met.
void TestEigenValueSolver::inversePowerIterationAcceptsExactZeroTolerance()
{
    RSparseMatrix m,k;
    m.setNRows(1);
    m.addValue(0,0,1.0);
    k.setNRows(1);
    k.addValue(0,0,2.0);

    REigenValueSolverConf eigenConf = makeEigenConf(REigenValueSolverConf::InversePowerIteration,1);
    eigenConf.setNIterations(2);
    eigenConf.setSolverCvgValue(0.0);
    REigenValueSolver solver(eigenConf,makeCgConf());

    RRVector d;
    RRMatrix ev;
    solver.solve(m,k,d,ev);

    QCOMPARE(d.getNRows(),1u);
    QCOMPARE(d[0],2.0);
    QVERIFY2(hasAccurateGeneralizedEigenpairs(m,k,d,ev),"Exact scalar eigenpair has a nonzero generalized residual.");
}

//! Exhausted subspace iteration must not publish an approximate mode.
void TestEigenValueSolver::subspaceIterationFailsWhenUnderiterated()
{
    RSparseMatrix m,k;
    buildUnderiteratedSystem(m,k);

    RRVector d;
    d.resize(1,0.0);
    d[0] = 2468.0;
    RRMatrix ev;
    ev.resize(1,3,0.0);
    ev[0][0] = 11.0;
    ev[0][1] = 12.0;
    ev[0][2] = 13.0;

    REigenValueSolverConf eigenConf = makeEigenConf(REigenValueSolverConf::SubspaceIteration,1);
    eigenConf.setNIterations(2);
    std::srand(41);
    REigenValueSolver solver(eigenConf,makeCgConf());

    bool didNotConverge = false;
    try
    {
        solver.solve(m,k,d,ev);
    }
    catch (const RError &error)
    {
        didNotConverge = error.getMessage().contains("did not converge");
    }

    QVERIFY2(didNotConverge,"Underiterated subspace solve did not fail with a convergence error.");
    QCOMPARE(d.getNRows(),1u);
    QCOMPARE(d[0],2468.0);
    QCOMPARE(ev.getNRows(),1u);
    QCOMPARE(ev.getNColumns(),3u);
    QCOMPARE(ev[0][0],11.0);
    QCOMPARE(ev[0][1],12.0);
    QCOMPARE(ev[0][2],13.0);
}

//! Exhausted inverse power iteration must not publish an approximate mode.
void TestEigenValueSolver::inversePowerIterationFailsWhenUnderiterated()
{
    RSparseMatrix m,k;
    buildUnderiteratedSystem(m,k);

    RRVector d;
    d.resize(1,0.0);
    d[0] = 2468.0;
    RRMatrix ev;
    ev.resize(1,3,0.0);
    ev[0][0] = 11.0;
    ev[0][1] = 12.0;
    ev[0][2] = 13.0;

    REigenValueSolverConf eigenConf = makeEigenConf(REigenValueSolverConf::InversePowerIteration,1);
    eigenConf.setNIterations(2);
    std::srand(41);
    REigenValueSolver solver(eigenConf,makeCgConf());

    bool didNotConverge = false;
    try
    {
        solver.solve(m,k,d,ev);
    }
    catch (const RError &error)
    {
        didNotConverge = error.getMessage().contains("did not converge");
    }

    QVERIFY2(didNotConverge,"Underiterated inverse power solve did not fail with a convergence error.");
    QCOMPARE(d.getNRows(),1u);
    QCOMPARE(d[0],2468.0);
    QCOMPARE(ev.getNRows(),1u);
    QCOMPARE(ev.getNColumns(),3u);
    QCOMPARE(ev[0][0],11.0);
    QCOMPARE(ev[0][1],12.0);
    QCOMPARE(ev[0][2],13.0);
}

QTEST_APPLESS_MAIN(TestEigenValueSolver)

#include "tst_eigen_value_solver.moc"

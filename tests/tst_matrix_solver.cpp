#include <QtTest>

#include <array>
#include <atomic>
#include <new>
#include <thread>

#include <omp.h>

#include "rmatrixsolver.h"
#include "rml_matrix_solver_conf.h"
#include "rml_sparse_matrix.h"
#include "rbl_rvector.h"

namespace
{
constexpr double tol = 1.0e-6;

RMatrixSolverConf makeCgConf()
{
    RMatrixSolverConf conf(RMatrixSolverConf::CG);
    conf.setNInnerIterations(1000);
    conf.setNOuterIterations(1000);
    conf.setSolverCvgValue(1.0e-13);
    conf.setOutputFrequency(0);
    return conf;
}
}

class TestMatrixSolver : public QObject
{
    Q_OBJECT

private slots:

    void solveSymmetricPositiveDefinite();
    void solveDiagonalSystem();
    void refreshesMutatedAndReusedMatrixAddress();
    void concurrentIndependentSolvers();
};

void TestMatrixSolver::solveSymmetricPositiveDefinite()
{
    // A = [[4,1],[1,3]], b = [1,2]  ->  x = [1/11, 7/11].
    RSparseMatrix a;
    a.setNRows(2);
    a.addValue(0, 0, 4.0);
    a.addValue(0, 1, 1.0);
    a.addValue(1, 0, 1.0);
    a.addValue(1, 1, 3.0);

    RRVector b(2, 0.0);
    b[0] = 1.0;
    b[1] = 2.0;

    RRVector x;
    RMatrixSolver solver(makeCgConf());
    solver.solve(a, b, x);

    QCOMPARE(x.getNRows(), 2u);
    QVERIFY(qAbs(x[0] - 1.0 / 11.0) < tol);
    QVERIFY(qAbs(x[1] - 7.0 / 11.0) < tol);
}

void TestMatrixSolver::solveDiagonalSystem()
{
    // Diagonal A trivially gives x[i] = b[i] / A[i][i].
    RSparseMatrix a;
    a.setNRows(3);
    a.addValue(0, 0, 2.0);
    a.addValue(1, 1, 4.0);
    a.addValue(2, 2, 5.0);

    RRVector b(3, 0.0);
    b[0] = 6.0;
    b[1] = 8.0;
    b[2] = 10.0;

    RRVector x;
    RMatrixSolver solver(makeCgConf());
    solver.solve(a, b, x);

    QVERIFY(qAbs(x[0] - 3.0) < tol);
    QVERIFY(qAbs(x[1] - 2.0) < tol);
    QVERIFY(qAbs(x[2] - 2.0) < tol);
}

void TestMatrixSolver::refreshesMutatedAndReusedMatrixAddress()
{
    // The cache must not treat a pointer as an immutable matrix identity.
    alignas(RSparseMatrix) unsigned char storage[sizeof(RSparseMatrix)];
    RSparseMatrix *a = new (storage) RSparseMatrix;
    a->setNRows(2);
    a->addValue(0,0,2.0);
    a->addValue(1,1,4.0);
    RMatrixSolver solver(makeCgConf());
    RRVector b(2,0.0), x;
    b[0] = 2.0; b[1] = 8.0;
    solver.solve(*a,b,x);
    QVERIFY(qAbs(x[0]-1.0) < tol && qAbs(x[1]-2.0) < tol);

    // Same pattern, different coefficients: refresh cached values.
    a->addValueAtPosition(0,0,3.0);
    b[0] = 5.0;
    solver.solve(*a,b,x);
    QVERIFY(qAbs(x[0]-1.0) < tol && qAbs(x[1]-2.0) < tol);

    // Same object, changed sparsity pattern: rebuild CSR indexes.
    a->addValue(0,1,1.0);
    a->addValue(1,0,1.0);
    b[0] = 7.0; b[1] = 9.0; // [[5,1],[1,4]] * [1,2]
    solver.solve(*a,b,x);
    QVERIFY(qAbs(x[0]-1.0) < tol && qAbs(x[1]-2.0) < tol);
    a->~RSparseMatrix();

    // New object at precisely the old address; neither shape nor pointer
    // alone may authorize reuse of the previous values/pattern.
    a = new (storage) RSparseMatrix;
    a->setNRows(2);
    a->addValue(0,0,11.0);
    a->addValue(1,1,13.0);
    b[0] = 11.0; b[1] = 26.0;
    solver.solve(*a,b,x);
    QVERIFY(qAbs(x[0]-1.0) < tol && qAbs(x[1]-2.0) < tol);
    a->~RSparseMatrix();
}

void TestMatrixSolver::concurrentIndependentSolvers()
{
    std::atomic<int> ready{0};
    std::atomic<bool> good{true};
    auto worker = [&](int workerId)
    {
        // Do not create a nested OpenMP team for each tiny test matrix.
        omp_set_num_threads(1);
        std::array<RSparseMatrix,12> matrices;
        for (size_t i=0;i<matrices.size();i++)
        {
            const double d0 = 2.0 + double(workerId) + double(i);
            const double d1 = 3.0 + double(workerId) + double(i);
            matrices[i].setNRows(2);
            matrices[i].addValue(0,0,d0);
            matrices[i].addValue(1,1,d1);
        }
        RMatrixSolver solver(makeCgConf());
        ready.fetch_add(1,std::memory_order_release);
        while (ready.load(std::memory_order_acquire) != 2) std::this_thread::yield();
        for (int pass=0;pass<3;pass++)
        {
            for (size_t i=0;i<matrices.size();i++)
            {
                const double d0 = 2.0 + double(workerId) + double(i);
                const double d1 = 3.0 + double(workerId) + double(i);
                RRVector b(2,0.0), x;
                b[0] = d0;
                b[1] = 2.0*d1;
                solver.solve(matrices[i],b,x);
                if (x.getNRows() != 2 || qAbs(x[0]-1.0) >= tol ||
                    qAbs(x[1]-2.0) >= tol)
                    good.store(false,std::memory_order_relaxed);
            }
        }
    };
    std::thread first(worker,0), second(worker,1);
    first.join();
    second.join();
    QVERIFY(good.load(std::memory_order_relaxed));
}

QTEST_APPLESS_MAIN(TestMatrixSolver)

#include "tst_matrix_solver.moc"

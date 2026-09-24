#include <QtTest>

#include <algorithm>
#include <cmath>
#include <initializer_list>
#include <limits>

#include <omp.h>

#include <rbl_error.h>

#include "rmodal_damping.h"
#include "rstructural_dynamics_integrator.h"

namespace
{
RSparseMatrix matrix(std::initializer_list<std::initializer_list<double>> rows)
{
    RSparseMatrix result;
    result.setNRows(uint(rows.size()));
    uint row = 0;
    for (const auto &values : rows)
    {
        uint column = 0;
        for (double value : values)
        {
            if (value != 0.0) result.addValue(row,column,value);
            column++;
        }
        row++;
    }
    return result;
}

RMatrixSolverConf cg()
{
    RMatrixSolverConf result(RMatrixSolverConf::CG);
    result.setNInnerIterations(200);
    result.setNOuterIterations(2000);
    result.setSolverCvgValue(1.0e-14);
    result.setOutputFrequency(0);
    return result;
}

bool close(double actual, double expected, double tolerance = 1.0e-9)
{
    return std::fabs(actual-expected) <= tolerance*std::max(1.0,std::fabs(expected));
}
}

class TestModalDamping : public QObject
{
    Q_OBJECT

private slots:
    void initTestCase();
    void twoModeAnalyticMatrix();
    void partialModesLeaveHigherModeUndamped();
    void nonidentityMassAndTransientDecay();
    void repeatedModesNeedUniformDampingAndWholeCluster();
    void repeatedExtractionHasStableDampingOperator();
    void unconvergedEigenvectorsAreRejected();
    void rejectsInvalidSystemsAndSettings();
};

void TestModalDamping::initTestCase()
{
    omp_set_num_threads(1);
}

void TestModalDamping::twoModeAnalyticMatrix()
{
    const RSparseMatrix M = matrix({{1.0,0.0},{0.0,1.0}});
    const RSparseMatrix K = matrix({{2.0,-1.0},{-1.0,2.0}});
    const RModalDampingResult result = RModalDamping::synthesize(M,K,{0.05,0.10},cg());
    QCOMPARE(result.eigenvalues.size(),2u);
    QVERIFY(close(result.eigenvalues[0],1.0));
    QVERIFY(close(result.eigenvalues[1],3.0));
    const double diagonal = (0.1+0.2*std::sqrt(3.0))/2.0;
    const double offDiagonal = (0.1-0.2*std::sqrt(3.0))/2.0;
    QVERIFY(close(result.damping.findValue(0,0),diagonal));
    QVERIFY(close(result.damping.findValue(1,1),diagonal));
    QVERIFY(close(result.damping.findValue(0,1),offDiagonal));
    QVERIFY(close(result.damping.findValue(1,0),offDiagonal));
}

void TestModalDamping::partialModesLeaveHigherModeUndamped()
{
    const RSparseMatrix M = matrix({{1.0,0.0},{0.0,1.0}});
    const RSparseMatrix K = matrix({{2.0,-1.0},{-1.0,2.0}});
    const RModalDampingResult result = RModalDamping::synthesize(M,K,{0.05},cg());
    for (uint row=0;row<2;row++)
    {
        for (uint column=0;column<2;column++)
        {
            QVERIFY(close(result.damping.findValue(row,column),0.05));
        }
    }
    RRVector highMode(2,0.0);
    highMode[0] = 1.0;
    highMode[1] = -1.0;
    RRVector dampingForce;
    RSparseMatrix::mlt(result.damping,highMode,dampingForce);
    QVERIFY(close(dampingForce[0],0.0));
    QVERIFY(close(dampingForce[1],0.0));
}

void TestModalDamping::nonidentityMassAndTransientDecay()
{
    const RSparseMatrix M = matrix({{2.0}});
    const RSparseMatrix K = matrix({{18.0}});
    const RModalDampingResult result = RModalDamping::synthesize(M,K,{0.05},cg());
    QVERIFY(close(result.eigenvalues[0],9.0));
    QVERIFY(close(result.modes[0][0]*result.modes[0][0]*2.0,1.0));
    QVERIFY(close(result.damping.findValue(0,0),0.6));

    RStructuralDynamicsIntegrator::Options options;
    RStructuralDynamicsIntegrator transient(M,result.damping,K,cg(),options,
        [](double, const RStructuralDynamicsState &, RRVector &force) { force[0] = 0.0; });
    transient.initialize(RRVector(1,1.0),RRVector(1,0.0));
    for (int i=0;i<1000;i++)
    {
        transient.trialStep(0.001);
        transient.acceptStep();
    }
    const double ratio = 0.05;
    const double omega = 3.0;
    const double dampedOmega = omega*std::sqrt(1.0-ratio*ratio);
    const double exact = std::exp(-ratio*omega)*
        (std::cos(dampedOmega)+ratio*omega/dampedOmega*std::sin(dampedOmega));
    QVERIFY(close(transient.acceptedState().q[0],exact,2.0e-6));
    QVERIFY(transient.acceptedState().dampingDissipation > 0.0);
}

void TestModalDamping::repeatedModesNeedUniformDampingAndWholeCluster()
{
    const RSparseMatrix M = matrix({{1.0,0.0},{0.0,1.0}});
    const RSparseMatrix K = matrix({{4.0,0.0},{0.0,4.0}});
    const RModalDampingResult equal = RModalDamping::synthesize(M,K,{0.05,0.05},cg());
    QVERIFY(close(equal.damping.findValue(0,0),0.2));
    QVERIFY(close(equal.damping.findValue(1,1),0.2));
    QVERIFY(close(equal.damping.findValue(0,1),0.0));
    QVERIFY_EXCEPTION_THROWN(RModalDamping::synthesize(M,K,{0.05},cg()),RError);
    QVERIFY_EXCEPTION_THROWN(RModalDamping::synthesize(M,K,{0.05,0.10},cg()),RError);
}

void TestModalDamping::repeatedExtractionHasStableDampingOperator()
{
    RSparseMatrix M;
    RSparseMatrix K;
    M.setNRows(9);
    K.setNRows(9);
    for (uint row=0;row<9;row++)
    {
        M.addValue(row,row,1.0);
        K.addValue(row,row,2.0);
        if (row > 0) K.addValue(row,row-1,-1.0);
        if (row+1 < 9) K.addValue(row,row+1,-1.0);
    }
    RModalDampingResult reference;
    try
    {
        reference = RModalDamping::synthesize(M,K,{0.03,0.07},cg());
    }
    catch (const RError &error)
    {
        QFAIL(qPrintable(QString("Reference extraction rejected: %1").arg(error.getMessage())));
    }
    for (int run=0;run<3;run++)
    {
        RModalDampingResult repeated;
        try
        {
            repeated = RModalDamping::synthesize(M,K,{0.03,0.07},cg());
        }
        catch (const RError &error)
        {
            QFAIL(qPrintable(QString("Repeated extraction rejected: %1").arg(error.getMessage())));
        }
        for (uint row=0;row<9;row++)
        {
            for (uint column=0;column<9;column++)
            {
                QVERIFY(close(repeated.damping.findValue(row,column),
                              reference.damping.findValue(row,column),1.0e-7));
            }
        }
    }
}

void TestModalDamping::unconvergedEigenvectorsAreRejected()
{
    RSparseMatrix M;
    RSparseMatrix K;
    M.setNRows(9);
    K.setNRows(9);
    for (uint row=0;row<9;row++)
    {
        M.addValue(row,row,1.0);
        K.addValue(row,row,2.0);
        if (row > 0) K.addValue(row,row-1,-1.0);
        if (row+1 < 9) K.addValue(row,row+1,-1.0);
    }
    RModalDamping::Options options;
    options.eigenIterations = 1;
    QVERIFY_EXCEPTION_THROWN(RModalDamping::synthesize(M,K,{0.03,0.07},cg(),options),RError);
}

void TestModalDamping::rejectsInvalidSystemsAndSettings()
{
    const RSparseMatrix M = matrix({{1.0,0.0},{0.0,1.0}});
    const RSparseMatrix K = matrix({{2.0,-1.0},{-1.0,2.0}});
    QVERIFY_EXCEPTION_THROWN(RModalDamping::synthesize(M,K,{0.05,0.10,0.05},cg()),RError);
    QVERIFY_EXCEPTION_THROWN(RModalDamping::synthesize(M,K,{-0.01},cg()),RError);
    QVERIFY_EXCEPTION_THROWN(RModalDamping::synthesize(M,K,{std::numeric_limits<double>::infinity()},cg()),RError);
    QVERIFY_EXCEPTION_THROWN(RModalDamping::synthesize(M,matrix({{0.0,0.0},{0.0,1.0}}),{0.05},cg()),RError);
    RModalDamping::Options options;
    options.eigenIterations = 0;
    QVERIFY_EXCEPTION_THROWN(RModalDamping::synthesize(M,K,{0.05},cg(),options),RError);

    RSparseMatrix oversized;
    oversized.setNRows(RModalDamping::MaxEquations+1);
    QVERIFY_EXCEPTION_THROWN(RModalDamping::synthesize(oversized,oversized,{0.05},cg()),RError);
}

QTEST_APPLESS_MAIN(TestModalDamping)

#include "tst_modal_damping.moc"

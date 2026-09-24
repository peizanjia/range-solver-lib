#include <QtTest>

#include <cmath>
#include <stdexcept>

#include <omp.h>

#include "rstructural_dynamics_integrator.h"

namespace
{
RSparseMatrix scalarMatrix(double value)
{
    RSparseMatrix matrix;
    matrix.setNRows(1);
    matrix.addValue(0,0,value);
    return matrix;
}

RMatrixSolverConf solverConf()
{
    RMatrixSolverConf conf(RMatrixSolverConf::CG);
    conf.setNInnerIterations(20);
    conf.setNOuterIterations(20);
    conf.setSolverCvgValue(1.0e-13);
    conf.setOutputFrequency(0);
    return conf;
}

using Integrator = RStructuralDynamicsIntegrator;

Integrator oscillator(double damping,
                      Integrator::Method method,
                      double forceMagnitude = 0.0,
                      double rhoInfinity = 1.0)
{
    Integrator::Options options;
    options.method = method;
    options.rhoInfinity = rhoInfinity;
    return Integrator(scalarMatrix(1.0),scalarMatrix(damping),scalarMatrix(4.0),
                      solverConf(),options,
                      [forceMagnitude](double, const RStructuralDynamicsState &, RRVector &force)
                      { force[0] = forceMagnitude; });
}

void advance(Integrator &integrator, double dt, int steps)
{
    for (int i=0; i<steps; ++i)
    {
        integrator.trialStep(dt);
        integrator.acceptStep();
    }
}
}

class TestStructuralDynamicsIntegrator : public QObject
{
    Q_OBJECT

private slots:
    void initTestCase();
    void freeVibrationAndEnergy();
    void forceStepAndWork();
    void dampedVibrationAndDissipation();
    void generalizedAlpha();
    void centralDifference();
    void centralDifferenceStabilityGuard();
    void trialRejectAndCheckpoint();
    void invalidStepAndForce();
};

void TestStructuralDynamicsIntegrator::initTestCase()
{
    // A one-DOF system should not pay for a full OpenMP team on every solve.
    omp_set_num_threads(1);
}

void TestStructuralDynamicsIntegrator::freeVibrationAndEnergy()
{
    auto coarse = oscillator(0.0,Integrator::Method::NewmarkAverageAcceleration);
    auto fine = oscillator(0.0,Integrator::Method::NewmarkAverageAcceleration);
    coarse.initialize(RRVector(1,1.0),RRVector(1,0.0));
    fine.initialize(RRVector(1,1.0),RRVector(1,0.0));
    advance(coarse,0.04,25);
    advance(fine,0.02,50);
    const double exact = std::cos(2.0);
    const double coarseError = std::abs(coarse.acceptedState().q[0]-exact);
    const double fineError = std::abs(fine.acceptedState().q[0]-exact);
    QVERIFY(fineError < coarseError/3.8);
    QVERIFY(std::abs(fine.acceptedState().kineticEnergy+
                     fine.acceptedState().strainEnergy-2.0) < 1.0e-11);
    QVERIFY(std::abs(fine.acceptedState().energyBalanceResidual) < 1.0e-11);
    QVERIFY(std::abs(fine.acceptedState().v[0]+2.0*std::sin(2.0)) < 0.001);
}

void TestStructuralDynamicsIntegrator::forceStepAndWork()
{
    auto integrator = oscillator(0.0,Integrator::Method::NewmarkAverageAcceleration,1.0);
    integrator.initialize(RRVector(1,0.0),RRVector(1,0.0));
    QCOMPARE(integrator.acceptedState().a[0],1.0);
    advance(integrator,0.01,100);
    const auto &state = integrator.acceptedState();
    QVERIFY(std::abs(state.q[0]-(1.0-std::cos(2.0))/4.0) < 4.0e-5);
    QVERIFY(std::abs(state.externalWork-state.q[0]) < 1.0e-11);
    QVERIFY(std::abs(state.energyBalanceResidual) < 1.0e-10);
}

void TestStructuralDynamicsIntegrator::dampedVibrationAndDissipation()
{
    auto integrator = oscillator(0.4,Integrator::Method::NewmarkAverageAcceleration);
    integrator.initialize(RRVector(1,1.0),RRVector(1,0.0));
    advance(integrator,0.005,200);
    const double wd = std::sqrt(4.0-0.2*0.2);
    const double exact = std::exp(-0.2)*(std::cos(wd)+0.2/wd*std::sin(wd));
    const auto &state = integrator.acceptedState();
    QVERIFY(std::abs(state.q[0]-exact) < 2.0e-5);
    QVERIFY(state.dampingDissipation > 0.0);
    QVERIFY(state.kineticEnergy+state.strainEnergy < 2.0);
    QVERIFY(std::abs(state.energyBalanceResidual) < 2.0e-4);
}

void TestStructuralDynamicsIntegrator::generalizedAlpha()
{
    auto noDamping = oscillator(0.0,Integrator::Method::GeneralizedAlpha,0.0,1.0);
    auto numericalDamping = oscillator(0.0,Integrator::Method::GeneralizedAlpha,0.0,0.3);
    noDamping.initialize(RRVector(1,1.0),RRVector(1,0.0));
    numericalDamping.initialize(RRVector(1,1.0),RRVector(1,0.0));
    advance(noDamping,0.02,50);
    advance(numericalDamping,0.02,50);
    QVERIFY(std::abs(noDamping.acceptedState().q[0]-std::cos(2.0)) < 0.001);
    QVERIFY(std::abs(noDamping.acceptedState().energyBalanceResidual) < 1.0e-10);
    QVERIFY(std::abs(numericalDamping.acceptedState().q[0]-std::cos(2.0)) < 0.003);
    QVERIFY(std::isfinite(numericalDamping.acceptedState().energyBalanceResidual));
}

void TestStructuralDynamicsIntegrator::centralDifference()
{
    auto coarse = oscillator(0.0,Integrator::Method::CentralDifference);
    auto fine = oscillator(0.0,Integrator::Method::CentralDifference);
    coarse.initialize(RRVector(1,1.0),RRVector(1,0.0));
    fine.initialize(RRVector(1,1.0),RRVector(1,0.0));
    advance(coarse,0.04,25);
    advance(fine,0.02,50);
    const double coarseError = std::abs(coarse.acceptedState().q[0]-std::cos(2.0));
    const double fineError = std::abs(fine.acceptedState().q[0]-std::cos(2.0));
    QVERIFY(fineError < coarseError/3.8);
    QVERIFY(std::abs(fine.acceptedState().v[0]+2.0*std::sin(2.0)) < 0.001);
    QVERIFY(std::abs(fine.acceptedState().energyBalanceResidual) < 0.001);
}

void TestStructuralDynamicsIntegrator::centralDifferenceStabilityGuard()
{
    auto integrator = oscillator(0.0,Integrator::Method::CentralDifference);
    integrator.initialize(RRVector(1,1.0),RRVector(1,0.0));
    // For M=1 and K=4, omega=2 and the exact central-difference limit is dt=1.
    QVERIFY_EXCEPTION_THROWN(integrator.trialStep(0.95),std::invalid_argument);
    QVERIFY(!integrator.hasCandidate());
    integrator.trialStep(0.5);
    integrator.acceptStep();
    QCOMPARE(integrator.acceptedState().time,0.5);
}

void TestStructuralDynamicsIntegrator::trialRejectAndCheckpoint()
{
    auto integrator = oscillator(0.0,Integrator::Method::NewmarkAverageAcceleration);
    integrator.initialize(RRVector(1,1.0),RRVector(1,0.0));
    const auto firstCandidate = integrator.trialStep(0.02);
    QVERIFY(integrator.hasCandidate());
    QCOMPARE(integrator.acceptedState().time,0.0);
    integrator.rejectStep();
    QVERIFY(!integrator.hasCandidate());
    QVERIFY(!integrator.hasCandidate());
    integrator.trialStep(0.02);
    QCOMPARE(integrator.candidateState().q[0],firstCandidate.q[0]);
    integrator.acceptStep();
    const auto saved = integrator.checkpoint();
    advance(integrator,0.02,5);
    integrator.trialStep(0.02);
    integrator.restore(saved);
    QCOMPARE(integrator.acceptedState().time,0.02);
    QVERIFY(!integrator.hasCandidate());
    QCOMPARE(integrator.acceptedState().q[0],saved.state.q[0]);
}

void TestStructuralDynamicsIntegrator::invalidStepAndForce()
{
    auto integrator = oscillator(0.0,Integrator::Method::NewmarkAverageAcceleration);
    integrator.initialize(RRVector(1,1.0),RRVector(1,0.0));
    QVERIFY_EXCEPTION_THROWN(integrator.trialStep(0.0),std::invalid_argument);
    integrator.trialStep(0.01);
    integrator.acceptStep();
    QVERIFY_EXCEPTION_THROWN(integrator.trialStep(0.02),std::invalid_argument);
    auto badForce = Integrator(scalarMatrix(1.0),scalarMatrix(0.0),scalarMatrix(4.0),
                               solverConf(),{},
                               [](double, const RStructuralDynamicsState &, RRVector &force)
                               { force.resize(2); });
    QVERIFY_EXCEPTION_THROWN(badForce.initialize(RRVector(1,0.0),RRVector(1,0.0)),
                             std::invalid_argument);
}

QTEST_APPLESS_MAIN(TestStructuralDynamicsIntegrator)

#include "tst_structural_dynamics_integrator.moc"

#include <QtTest>

#include <cmath>
#include <utility>
#include <vector>

#include <omp.h>

#include "rstructural_dynamics_integrator.h"

namespace
{
using Integrator = RStructuralDynamicsIntegrator;

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

Integrator oscillator(double damping, Integrator::ForceCallback force,
                      Integrator::Options options = {})
{
    // The analytic reference system is q'' + damping*q' + 4*q = f(t).
    return Integrator(scalarMatrix(1.0),scalarMatrix(damping),scalarMatrix(4.0),
                      solverConf(),options,std::move(force));
}

void advance(Integrator &integrator, double dt, int nSteps)
{
    for (int i=0; i<nSteps; ++i)
    {
        integrator.trialStep(dt);
        integrator.acceptStep();
    }
}

double restStartSineResponse(double time)
{
    // Unit sinusoidal force: q''+4q=sin(t), q(0)=q'(0)=0.
    return (std::sin(time)-0.5*std::sin(2.0*time))/3.0;
}
}

class TestStructuralDynamicsLoads : public QObject
{
    Q_OBJECT

private slots:
    void initTestCase();
    void sinusoidalSteadyState();
    void sinusoidalTransientAndTimeStepConvergence();
    void finiteSmoothPulse();
    void generalizedAlphaSamplesIntermediateTime();
};

void TestStructuralDynamicsLoads::initTestCase()
{
    // Keep the scalar tests independent of a machine's default OpenMP team size.
    omp_set_num_threads(1);
}

void TestStructuralDynamicsLoads::sinusoidalSteadyState()
{
    // A harmonic particular solution, including its phase lag, is exact for
    // these initial conditions; no start-up transient obscures the comparison.
    constexpr double damping = 0.4;
    const double amplitude = 1.0/std::sqrt(3.0*3.0+damping*damping);
    const double phase = std::atan2(damping,3.0);
    auto integrator = oscillator(damping,
                                 [](double time, const RStructuralDynamicsState &, RRVector &force)
                                 { force[0] = std::sin(time); });
    integrator.initialize(RRVector(1,-amplitude*std::sin(phase)),
                          RRVector(1,amplitude*std::cos(phase)));
    advance(integrator,0.02,100);

    const auto &state = integrator.acceptedState();
    QVERIFY(std::abs(state.time-2.0) < 1.0e-12);
    // Average-acceleration Newmark is second order; at dt=0.02 its phase
    // error is O(dt^2), comfortably below this 1e-4 displacement bound.
    QVERIFY(std::abs(state.q[0]-amplitude*std::sin(2.0-phase)) < 1.0e-4);
    QVERIFY(std::abs(state.v[0]-amplitude*std::cos(2.0-phase)) < 1.0e-4);
}

void TestStructuralDynamicsLoads::sinusoidalTransientAndTimeStepConvergence()
{
    auto force = [](double time, const RStructuralDynamicsState &, RRVector &result)
                 { result[0] = std::sin(time); };
    auto coarse = oscillator(0.0,force);
    auto fine = oscillator(0.0,force);
    coarse.initialize(RRVector(1,0.0),RRVector(1,0.0));
    fine.initialize(RRVector(1,0.0),RRVector(1,0.0));
    advance(coarse,0.04,50);
    advance(fine,0.02,100);

    const double exact = restStartSineResponse(2.0);
    const double coarseError = std::abs(coarse.acceptedState().q[0]-exact);
    const double fineError = std::abs(fine.acceptedState().q[0]-exact);
    QVERIFY(coarseError > 1.0e-7); // A nontrivial ratio, not numerical noise.
    QVERIFY(fineError < coarseError/3.7); // Expected ratio: 4 for order two.
    QVERIFY(fineError < 1.0e-4);
}

void TestStructuralDynamicsLoads::finiteSmoothPulse()
{
    // f(t)=sin(pi*t) for 0<=t<=1, then zero. Its continuous endpoints avoid
    // a step-force ambiguity at the boundary of a time step.
    const double pi = std::acos(-1.0);
    auto integrator = oscillator(0.0,
                                 [pi](double time, const RStructuralDynamicsState &, RRVector &force)
                                 { force[0] = time >= 0.0 && time <= 1.0 ?
                                              std::sin(pi*time) : 0.0; });
    integrator.initialize(RRVector(1,0.0),RRVector(1,0.0));
    // Halving dt from 0.01 to 0.005 reduces the pulse-end displacement
    // error from 4.0e-5 to 1.0e-5, as expected for a second-order scheme.
    advance(integrator,0.005,200);
    const double qAtPulseEnd = (std::sin(pi)-0.5*pi*std::sin(2.0))/(4.0-pi*pi);
    const double vAtPulseEnd = pi*(std::cos(pi)-std::cos(2.0))/(4.0-pi*pi);
    QVERIFY(std::abs(integrator.acceptedState().q[0]-qAtPulseEnd) < 1.5e-5);
    const double workAtPulseEnd = integrator.acceptedState().externalWork;

    advance(integrator,0.005,200);
    const double exactAfterPulse = qAtPulseEnd*std::cos(2.0)+
                                   0.5*vAtPulseEnd*std::sin(2.0);
    const auto &state = integrator.acceptedState();
    QVERIFY(std::abs(state.q[0]-exactAfterPulse) < 1.0e-5);
    QVERIFY(std::abs(state.externalWork-workAtPulseEnd) < 1.0e-11);
    QVERIFY(std::abs(state.kineticEnergy+state.strainEnergy-workAtPulseEnd) < 1.0e-5);
}

void TestStructuralDynamicsLoads::generalizedAlphaSamplesIntermediateTime()
{
    std::vector<std::pair<double,double>> samples;
    Integrator::Options options;
    options.method = Integrator::Method::GeneralizedAlpha;
    options.rhoInfinity = 0.5; // alpha_f=rho_inf/(1+rho_inf)=1/3.
    auto integrator = oscillator(0.0,
                                 [&samples](double time, const RStructuralDynamicsState &accepted,
                                            RRVector &force)
                                 {
                                     samples.emplace_back(time,accepted.time);
                                     force[0] = std::sin(time);
                                 },options);
    integrator.initialize(RRVector(1,0.0),RRVector(1,0.0),1.25);
    integrator.trialStep(0.1);

    QCOMPARE(samples.size(),std::size_t(3));
    QVERIFY(std::abs(samples[0].first-1.25) < 1.0e-12);
    QVERIFY(std::abs(samples[1].first-(1.25+0.1*2.0/3.0)) < 1.0e-12);
    QVERIFY(std::abs(samples[2].first-1.35) < 1.0e-12);
    for (const auto &sample : samples)
        QVERIFY(std::abs(sample.second-1.25) < 1.0e-12);
    QCOMPARE(integrator.acceptedState().time,1.25);
    integrator.acceptStep();
    QVERIFY(std::abs(integrator.acceptedState().time-1.35) < 1.0e-12);
}

QTEST_APPLESS_MAIN(TestStructuralDynamicsLoads)

#include "tst_structural_dynamics_loads.moc"

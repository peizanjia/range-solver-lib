#include <QtTest>

#include <cmath>
#include <algorithm>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <utility>
#include <array>

#include <omp.h>

#include "rstructural_dynamics_integrator.h"
#include "rparallel_execution_context.h"
#include <rbl_error.h>

namespace
{
RSparseMatrix scalarMatrix(double value)
{
    RSparseMatrix matrix;
    matrix.setNRows(1);
    matrix.addValue(0,0,value);
    return matrix;
}

RSparseMatrix diagonalMatrix(double first, double second)
{
    RSparseMatrix matrix;
    matrix.setNRows(2);
    matrix.addValue(0,0,first);
    matrix.addValue(1,1,second);
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

Integrator twoModeOscillator(double rhoInfinity)
{
    Integrator::Options options;
    options.method = Integrator::Method::GeneralizedAlpha;
    options.rhoInfinity = rhoInfinity;
    return Integrator(diagonalMatrix(1.0,1.0),diagonalMatrix(0.0,0.0),
                      diagonalMatrix(4.0,4.0e6),solverConf(),options,
                      [](double, const RStructuralDynamicsState &, RRVector &force)
                      {
                          force[0] = 0.0;
                          force[1] = 0.0;
                      });
}

double modalEnergy(const RStructuralDynamicsState &state, int mode, double stiffness)
{
    return 0.5*(state.v[mode]*state.v[mode]+stiffness*state.q[mode]*state.q[mode]);
}

void advance(Integrator &integrator, double dt, int steps)
{
    for (int i=0; i<steps; ++i)
    {
        integrator.trialStep(dt);
        integrator.acceptStep();
    }
}
void checkpointBitsEqual(const RStructuralDynamicsCheckpoint &actual,const RStructuralDynamicsCheckpoint &expected)
{
    for(const auto &vectors:{std::make_pair(&actual.state.q,&expected.state.q),std::make_pair(&actual.state.v,&expected.state.v),
                            std::make_pair(&actual.state.a,&expected.state.a),std::make_pair(&actual.state.externalForce,&expected.state.externalForce)}){
        QCOMPARE(vectors.first->size(),vectors.second->size());
        QVERIFY(std::memcmp(vectors.first->data(),vectors.second->data(),vectors.first->size()*sizeof(double))==0);
    }
    QCOMPARE(actual.state.time,expected.state.time);QCOMPARE(actual.fixedStep,expected.fixedStep);
    QCOMPARE(actual.initialMechanicalEnergy,expected.initialMechanicalEnergy);
    QCOMPARE(actual.state.kineticEnergy,expected.state.kineticEnergy);QCOMPARE(actual.state.strainEnergy,expected.state.strainEnergy);
    QCOMPARE(actual.state.externalWork,expected.state.externalWork);QCOMPARE(actual.state.dampingDissipation,expected.state.dampingDissipation);
    QCOMPARE(actual.state.energyBalanceResidual,expected.state.energyBalanceResidual);
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
    void generalizedAlphaLowFrequencyConvergence_data();
    void generalizedAlphaLowFrequencyConvergence();
    void generalizedAlphaAttenuatesHighModeAndPreservesLowMode();
    void centralDifference();
    void centralDifferenceStabilityGuard();
    void trialRejectAndCheckpoint();
    void acceptedPowerLedgerIsCommittedAndRestored();
    void heldControlJumpRefreshesBoundaryAcceleration();
    void invalidStepAndForce();
    void canonicalGridTwentyThousandStepsAndRestore();
    void canonicalGridRejectsInvalidClockAndStep();
    void scaleSeparatedNewmarkKeepsSmallForceRows_data();
    void scaleSeparatedNewmarkKeepsSmallForceRows();
    void cancellingOriginalMassRows_data();
    void cancellingOriginalMassRows();
    void correctionEquilibrationPreservesOriginalRows();
    void tinyCouplingsReachExactlyUnloadedFortyRowChain();
};

void TestStructuralDynamicsIntegrator::initTestCase()
{
    // A one-DOF system should not pay for a full OpenMP team on every solve.
    omp_set_num_threads(1);
}

void TestStructuralDynamicsIntegrator::cancellingOriginalMassRows_data()
{
    QTest::addColumn<int>("solverType");
    QTest::newRow("GMRES")<<int(RMatrixSolverConf::GMRES);
    QTest::newRow("CG")<<int(RMatrixSolverConf::CG);
}

void TestStructuralDynamicsIntegrator::cancellingOriginalMassRows()
{
    QFETCH(int,solverType);
    // SPD coupled inertia: the last load is the remainder of opposing finite
    // contributions, while inertia spans twelve orders. Certify the ORIGINAL
    // mass rows independently of the production norm-stopped Krylov solve.
    const std::array<double,9> entries{{1.,1e-4,1e-4,1e-4,1e6,1e-4,1e-4,1e-4,1e-6}};
    const std::array<double,3> seed{{.3,-3e-11,-30.}};
    RSparseMatrix mass,zero;mass.setNRows(3);zero.setNRows(3);RRVector force(3,0.);
    for(uint i=0;i<3;++i)for(uint j=0;j<3;++j){mass.addValue(i,j,entries[3*i+j]);force[i]+=entries[3*i+j]*seed[j];}
    QVERIFY(force[2]!=0.);QVERIFY(std::abs(force[2])<1e-10*std::abs(force[0]));
    Integrator integrator(mass,zero,zero,RMatrixSolverConf(static_cast<RMatrixSolverType>(solverType)),{},
        [force](double,const RStructuralDynamicsState &,RRVector &out){out=force;});
    try{integrator.initialize(RRVector(3,0.),RRVector(3,0.));}
    catch(const RError &e){QFAIL(qPrintable(e.getMessage()));}
    catch(const std::exception &e){QFAIL(e.what());}
    const auto &acceleration=integrator.acceptedState().a;
    double norm=0;
    for(uint i=0;i<3;++i){
        // Neumaier sum, with separate exact product tails. Its summation
        // algorithm is independent of the production TwoSum residual helper.
        double sum=-force[i],tail=0,scale=std::abs(force[i]);
        const auto add=[&](double value){const double next=sum+value;
            tail+=std::abs(sum)>=std::abs(value)?(sum-next)+value:(value-next)+sum;sum=next;};
        for(uint j=0;j<3;++j){const double value=entries[3*i+j]*acceleration[j];
            add(value);add(std::fma(entries[3*i+j],acceleration[j],-value));scale+=std::abs(value);}
        const double residual=sum+tail;norm=std::hypot(norm,residual);
        QVERIFY2(std::abs(residual)<=1e-12*scale,qPrintable(QString("original mass row %1: residual=%2 scale=%3")
            .arg(i).arg(residual,0,'g',17).arg(scale,0,'g',17)));
    }
    QVERIFY(norm<=1e-8*RRVector::euclideanNorm(force)+1e-12);
    integrator.refreshAcceptedForce();
    QVERIFY(!integrator.hasCandidate());
}

void TestStructuralDynamicsIntegrator::correctionEquilibrationPreservesOriginalRows()
{
    // Six orders of coordinate scaling create twelve orders of matrix scale.
    // G is SPD tridiagonal, A=D G D. Keep more coordinates than the default
    // restarted GMRES basis, and check every ORIGINAL row with a distinct
    // compensated summation, never a norm of the equilibrated equation.
    constexpr uint n=16;RSparseMatrix mass,zero;mass.setNRows(n);zero.setNRows(n);
    std::array<double,n> scales{},seed{};RRVector force(n,0.);
    for(uint i=0;i<n;++i){scales[i]=i%2?1e3:1e-3;seed[i]=(i%3?-.2:.3)/scales[i];}
    for(uint i=0;i<n;++i)for(uint j=0;j<n;++j){
        const double g=i==j?4.:(i+1==j || j+1==i?-1.:0.);
        if(g){const double value=scales[i]*g*scales[j];mass.addValue(i,j,value);force[i]+=value*seed[j];}}
    for(auto kind:{RMatrixSolverConf::CG,RMatrixSolverConf::GMRES}){
        Integrator integrator(mass,zero,zero,RMatrixSolverConf(kind),{},
            [force](double,const RStructuralDynamicsState &,RRVector &out){out=force;});
        try{integrator.initialize(RRVector(n,0.),RRVector(n,0.));}
        catch(const RError &e){QFAIL(qPrintable(e.getMessage()));}
        catch(const std::exception &e){QFAIL(e.what());}
        const auto &a=integrator.acceptedState().a;double norm=0;
        for(uint i=0;i<n;++i){double sum=-force[i],tail=0,scale=std::abs(force[i]);
            const auto add=[&](double value){const double next=sum+value;
                tail+=std::abs(sum)>=std::abs(value)?(sum-next)+value:(value-next)+sum;sum=next;};
            for(uint j=0;j<n;++j){const double entry=mass.findValue(i,j),value=entry*a[j];
                add(value);add(std::fma(entry,a[j],-value));scale+=std::abs(value);}
            const double residual=sum+tail;norm=std::hypot(norm,residual);
            QVERIFY2(std::abs(residual)<=1e-12*scale,qPrintable(QString("solver %1 original row %2: residual=%3 scale=%4")
                .arg(int(kind)).arg(i).arg(residual,0,'g',17).arg(scale,0,'g',17)));
        }
        QVERIFY(norm<=1e-8*RRVector::euclideanNorm(force)+1e-12);
        const auto before=integrator.checkpoint();integrator.trialStep(.001);integrator.rejectStep();
        checkpointBitsEqual(integrator.checkpoint(),before);
    }
}

void TestStructuralDynamicsIntegrator::tinyCouplingsReachExactlyUnloadedFortyRowChain()
{
    // Two independently specified SPD spring chains, coupled by finite 1e-14
    // entries. Coordinate scales produce twelve orders of matrix scale. The
    // first chain has EXACTLY zero RHS; its induced acceleration is nevertheless
    // nonzero. A global norm can hide errors in all twenty unloaded rows.
    constexpr uint n=40;RSparseMatrix mass,zero;mass.setNRows(n);zero.setNRows(n);
    std::array<double,n> coordinateScale{};RRVector force(n,0.0);
    for(uint i=0;i<n;++i)coordinateScale[i]=i<20?1e3:1e-3;
    for(uint base:{0u,20u})for(uint i=0;i<20;++i){
        mass.addValue(base+i,base+i,(.001+(i==0 || i==19?1.:2.))*coordinateScale[base+i]*coordinateScale[base+i]);
        if(i<19){const double value=-coordinateScale[base+i]*coordinateScale[base+i+1];
            mass.addValue(base+i,base+i+1,value);mass.addValue(base+i+1,base+i,value);}}
    for(uint i=0;i<20;++i){const double coupling=1e-14*coordinateScale[i]*coordinateScale[i+20];
        QVERIFY(coupling!=0.);mass.addValue(i,i+20,coupling);mass.addValue(i+20,i,coupling);
        force[i+20]=coordinateScale[i+20]*(i%3==0?.3:-.2);QCOMPARE(force[i],0.);}
    for(auto kind:{RMatrixSolverConf::CG,RMatrixSolverConf::GMRES}){
        bool invalidForce=false;Integrator integrator(mass,zero,zero,RMatrixSolverConf(kind),{},
            [&](double,const RStructuralDynamicsState &,RRVector &out){out=force;
                if(invalidForce)out[20]=std::numeric_limits<double>::infinity();});
        try{integrator.initialize(RRVector(n,0.),RRVector(n,0.));}
        catch(const RError &e){QFAIL(qPrintable(e.getMessage()));}
        catch(const std::exception &e){QFAIL(e.what());}
        const auto &a=integrator.acceptedState().a;double norm=0,unloadedMotion=0;
        for(uint i=0;i<n;++i){double sum=-force[i],tail=0,scale=std::abs(force[i]);
            const auto add=[&](double value){const double next=sum+value;
                tail+=std::abs(sum)>=std::abs(value)?(sum-next)+value:(value-next)+sum;sum=next;};
            for(uint j=0;j<n;++j){const double entry=mass.findValue(i,j),value=entry*a[j];
                add(value);add(std::fma(entry,a[j],-value));scale+=std::abs(value);}
            const double residual=sum+tail;norm=std::hypot(norm,residual);
            QVERIFY2(scale>0 && std::abs(residual)<=1e-12*scale,qPrintable(QString("solver %1 original row %2: residual=%3 scale=%4 zeroRHS=%5")
                .arg(int(kind)).arg(i).arg(residual,0,'g',17).arg(scale,0,'g',17).arg(force[i]==0)));
            if(i<20)unloadedMotion=std::hypot(unloadedMotion,a[i]);
        }
        QVERIFY(unloadedMotion>0);QVERIFY(norm<=1e-8*RRVector::euclideanNorm(force)+1e-12);
        const auto accepted=integrator.checkpoint();invalidForce=true;bool rejected=false;
        try{integrator.trialStep(.001);}catch(const RError &){rejected=true;}catch(const std::exception &){rejected=true;}
        QVERIFY(rejected);QVERIFY(!integrator.hasCandidate());checkpointBitsEqual(integrator.checkpoint(),accepted);
        invalidForce=false;integrator.trialStep(.001);integrator.rejectStep();checkpointBitsEqual(integrator.checkpoint(),accepted);
    }
}

void TestStructuralDynamicsIntegrator::scaleSeparatedNewmarkKeepsSmallForceRows_data()
{
    QTest::addColumn<int>("solverType");
    QTest::newRow("GMRES-default-budget")<<int(RMatrixSolverConf::GMRES);
    QTest::newRow("CG-default-budget")<<int(RMatrixSolverConf::CG);
}

void TestStructuralDynamicsIntegrator::scaleSeparatedNewmarkKeepsSmallForceRows()
{
    QFETCH(int,solverType);
    // SPD coupled coordinates with six orders of inertia separation and a
    // small loaded row obtained from cancelling finite terms, not a zero-force
    // workaround. The independent 2x2 solve verifies the complete trajectory.
    using V2=std::array<double,2>;using M2=std::array<double,4>;
    const M2 M{{1.,2e-4,2e-4,1e6}};M2 C{},K{};
    for(size_t i=0;i<4;++i){C[i]=.04*M[i];K[i]=4*M[i];}
    const auto multiply=[](const M2 &matrix,const V2 &x){return V2{{matrix[0]*x[0]+matrix[1]*x[1],matrix[2]*x[0]+matrix[3]*x[1]}};};
    const auto solve2=[](const M2 &matrix,const V2 &b){const double determinant=matrix[0]*matrix[3]-matrix[1]*matrix[2];
        return V2{{(matrix[3]*b[0]-matrix[1]*b[1])/determinant,(matrix[0]*b[1]-matrix[2]*b[0])/determinant}};};
    const auto sparse2=[](const M2 &matrix){RSparseMatrix result;result.setNRows(2);
        for(uint i=0;i<2;++i)for(uint j=0;j<2;++j)result.addValue(i,j,matrix[2*i+j]);return result;};
    const auto force=multiply(M,V2{{.3,-6e-11+3e-15}});QVERIFY(force[1]!=0);QVERIFY(std::abs(force[1]/force[0])<2e-8);
    RMatrixSolverConf policy(static_cast<RMatrixSolverType>(solverType));bool overflow=false;
    const auto callback=[&](double,const RStructuralDynamicsState &,RRVector &f){
        for(uint i=0;i<2;++i)f[i]=overflow?std::numeric_limits<double>::max():force[i];};
    const auto make=[&]{return Integrator(sparse2(M),sparse2(C),sparse2(K),policy,{},callback);};
    auto integrator=make();integrator.initialize(RRVector(2,0.),RRVector(2,0.));
    V2 q{{0,0}},v{{0,0}},a=solve2(M,force);const double h=.01;RStructuralDynamicsCheckpoint saved;
    for(std::uint64_t index=1;index<=20;++index){
        if(index==2){const auto before=integrator.checkpoint();overflow=true;bool failed=false;
            try{integrator.trialStepOnGrid(0,index,h);}catch(const RError &){failed=true;}catch(const std::exception &){failed=true;}
            overflow=false;QVERIFY(failed);QVERIFY(!integrator.hasCandidate());checkpointBitsEqual(integrator.checkpoint(),before);
        }
        V2 qp{},vp{};for(size_t i=0;i<2;++i){qp[i]=q[i]+h*v[i]+.25*h*h*a[i];vp[i]=v[i]+.5*h*a[i];}
        const auto cv=multiply(C,vp),kq=multiply(K,qp);V2 rhs{};M2 effective{};
        for(size_t i=0;i<2;++i)rhs[i]=force[i]-cv[i]-kq[i];
        for(size_t i=0;i<4;++i)effective[i]=M[i]+.5*h*C[i]+.25*h*h*K[i];
        a=solve2(effective,rhs);for(size_t i=0;i<2;++i){q[i]=qp[i]+.25*h*h*a[i];v[i]=vp[i]+.5*h*a[i];}
        integrator.trialStepOnGrid(0,index,h);const auto &candidate=integrator.candidateState();
        for(size_t i=0;i<2;++i){
            const auto relative=[](double actual,double expected){return std::abs(actual-expected)<=2e-10*std::max(std::abs(expected),1e-25);};
            QVERIFY2(relative(candidate.q[i],q[i]),"Scale-separated q differs from independent Newmark reference");
            QVERIFY2(relative(candidate.v[i],v[i]),"Scale-separated v differs from independent Newmark reference");
            QVERIFY2(relative(candidate.a[i],a[i]),"Scale-separated a differs from independent Newmark reference");
        }
        const V2 actualQ{{candidate.q[0],candidate.q[1]}},actualV{{candidate.v[0],candidate.v[1]}},actualA{{candidate.a[0],candidate.a[1]}};
        const auto ma=multiply(M,actualA),cvActual=multiply(C,actualV),kqActual=multiply(K,actualQ);
        for(size_t i=0;i<2;++i){const double scale=std::abs(force[i])+std::abs(M[2*i]*actualA[0])+std::abs(M[2*i+1]*actualA[1])
            +std::abs(C[2*i]*actualV[0])+std::abs(C[2*i+1]*actualV[1])+std::abs(K[2*i]*actualQ[0])+std::abs(K[2*i+1]*actualQ[1]);
            QVERIFY(std::abs(ma[i]+cvActual[i]+kqActual[i]-force[i])<=2e-12*scale);}
        integrator.acceptStep();if(index==7)saved=integrator.checkpoint();
    }
    const auto expected=integrator.checkpoint();auto restored=make();restored.restore(saved);
    for(std::uint64_t index=8;index<=20;++index){restored.trialStepOnGrid(0,index,h);restored.acceptStep();}
    checkpointBitsEqual(restored.checkpoint(),expected);
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

void TestStructuralDynamicsIntegrator::generalizedAlphaLowFrequencyConvergence_data()
{
    QTest::addColumn<double>("rhoInfinity");
    QTest::newRow("rho-0.0") << 0.0;
    QTest::newRow("rho-0.3") << 0.3;
    QTest::newRow("rho-0.7") << 0.7;
    QTest::newRow("rho-1.0") << 1.0;
}

void TestStructuralDynamicsIntegrator::generalizedAlphaLowFrequencyConvergence()
{
    QFETCH(double,rhoInfinity);
    auto coarse = oscillator(0.0,Integrator::Method::GeneralizedAlpha,0.0,rhoInfinity);
    auto fine = oscillator(0.0,Integrator::Method::GeneralizedAlpha,0.0,rhoInfinity);
    coarse.initialize(RRVector(1,1.0),RRVector(1,0.0));
    fine.initialize(RRVector(1,1.0),RRVector(1,0.0));
    // Use the same fine step as the separate two-mode accuracy benchmark.
    // The rho=0 scheme has a larger O(dt^2) error constant; dt=0.02 does
    // not meet the 5e-4 absolute bound even though it converges quadratically.
    advance(coarse,0.02,50);
    advance(fine,0.01,100);

    const double exact = std::cos(2.0);
    const double coarseError = std::abs(coarse.acceptedState().q[0]-exact);
    const double fineError = std::abs(fine.acceptedState().q[0]-exact);
    QVERIFY(fineError > 0.0);
    const double refinementRatio = coarseError/fineError;
    QVERIFY2(refinementRatio > 3.7,
             qPrintable(QString("rhoInfinity=%1, error refinement ratio=%2")
                            .arg(rhoInfinity).arg(refinementRatio)));
    QVERIFY2(fineError < 5.0e-4,
             qPrintable(QString("rhoInfinity=%1, fine-step displacement error=%2")
                            .arg(rhoInfinity).arg(fineError)));
}

void TestStructuralDynamicsIntegrator::generalizedAlphaAttenuatesHighModeAndPreservesLowMode()
{
    auto numericalDamping = twoModeOscillator(0.3);
    auto noNumericalDamping = twoModeOscillator(1.0);
    RRVector initialDisplacement(2,0.0);
    initialDisplacement[0] = 1.0;
    initialDisplacement[1] = 0.001;
    RRVector initialVelocity(2,0.0);
    numericalDamping.initialize(initialDisplacement,initialVelocity);
    noNumericalDamping.initialize(initialDisplacement,initialVelocity);

    const double initialHighModeEnergy = modalEnergy(numericalDamping.acceptedState(),1,4.0e6);
    advance(numericalDamping,0.01,100);
    advance(noNumericalDamping,0.01,100);

    const auto &dampedState = numericalDamping.acceptedState();
    const auto &undampedState = noNumericalDamping.acceptedState();
    const double attenuatedHighModeEnergy = modalEnergy(dampedState,1,4.0e6);
    const double preservedHighModeEnergy = modalEnergy(undampedState,1,4.0e6);
    const double lowModeError = std::abs(dampedState.q[0]-std::cos(2.0));

    QVERIFY(lowModeError < 5.0e-4);
    QVERIFY(attenuatedHighModeEnergy/initialHighModeEnergy < 1.0e-6);
    QVERIFY(std::abs(preservedHighModeEnergy/initialHighModeEnergy-1.0) < 1.0e-8);
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

void TestStructuralDynamicsIntegrator::acceptedPowerLedgerIsCommittedAndRestored()
{
    auto integrator = oscillator(0.4,Integrator::Method::NewmarkAverageAcceleration,1.0);
    QVERIFY_EXCEPTION_THROWN(integrator.acceptedPowerSample(),std::logic_error);
    integrator.initialize(RRVector(1,0.0),RRVector(1,2.0));
    const auto initial = integrator.acceptedPowerSample();
    QCOMPARE(initial.time,0.0);
    QCOMPARE(initial.effectiveForcePower,2.0);
    QCOMPARE(initial.viscousDissipationPower,1.6);
    QCOMPARE(initial.cumulativeEffectiveWork,0.0);
    QCOMPARE(initial.cumulativeDampingDissipation,0.0);

    integrator.trialStep(0.01);
    QCOMPARE(integrator.acceptedPowerSample().effectiveForcePower,initial.effectiveForcePower);
    integrator.rejectStep();
    const auto afterReject = integrator.acceptedPowerSample();
    QCOMPARE(afterReject.time,initial.time);
    QCOMPARE(afterReject.effectiveForcePower,initial.effectiveForcePower);
    QCOMPARE(afterReject.viscousDissipationPower,initial.viscousDissipationPower);
    QCOMPARE(afterReject.cumulativeEffectiveWork,initial.cumulativeEffectiveWork);
    QCOMPARE(afterReject.cumulativeDampingDissipation,initial.cumulativeDampingDissipation);
    QCOMPARE(afterReject.energyBalanceResidual,initial.energyBalanceResidual);

    integrator.trialStep(0.01);
    integrator.acceptStep();
    const auto accepted = integrator.acceptedPowerSample();
    const auto &state = integrator.acceptedState();
    QVERIFY(std::abs(accepted.effectiveForcePower-state.v[0]) < 1.0e-12);
    QVERIFY(std::abs(accepted.viscousDissipationPower-0.4*state.v[0]*state.v[0]) < 1.0e-12);
    // The applied force is exactly one, so its physical work over this step
    // is the displacement increment, independent of the ledger implementation.
    QVERIFY(std::abs(accepted.cumulativeEffectiveWork-state.q[0]) < 1.0e-12);
    const double expectedDissipation = 0.005*(1.6+0.4*state.v[0]*state.v[0]);
    QVERIFY(std::abs(accepted.cumulativeDampingDissipation-expectedDissipation) < 1.0e-12);
    QCOMPARE(accepted.cumulativeEffectiveWork,state.externalWork);
    QCOMPARE(accepted.cumulativeDampingDissipation,state.dampingDissipation);
    QCOMPARE(accepted.energyBalanceResidual,state.energyBalanceResidual);
    const auto saved = integrator.checkpoint();
    integrator.trialStep(0.01);
    integrator.acceptStep();
    integrator.restore(saved);
    const auto restored = integrator.acceptedPowerSample();
    QCOMPARE(restored.time,accepted.time);
    QCOMPARE(restored.effectiveForcePower,accepted.effectiveForcePower);
    QCOMPARE(restored.viscousDissipationPower,accepted.viscousDissipationPower);
    QCOMPARE(restored.cumulativeEffectiveWork,accepted.cumulativeEffectiveWork);
    QCOMPARE(restored.cumulativeDampingDissipation,accepted.cumulativeDampingDissipation);
    QCOMPARE(restored.energyBalanceResidual,accepted.energyBalanceResidual);
}

void TestStructuralDynamicsIntegrator::heldControlJumpRefreshesBoundaryAcceleration()
{
    double command = 0.0;
    Integrator integrator(scalarMatrix(1.0),scalarMatrix(0.0),scalarMatrix(4.0),
                          solverConf(),{},
                          [&command](double, const RStructuralDynamicsState &, RRVector &force)
                          { force[0] = command; });
    integrator.initialize(RRVector(1,0.0),RRVector(1,0.0));
    QCOMPARE(integrator.acceptedState().a[0],0.0);
    command = 1.0;
    integrator.refreshAcceptedForce();
    QCOMPARE(integrator.acceptedState().a[0],1.0);
    QCOMPARE(integrator.acceptedState().externalForce[0],1.0);
    // Average-acceleration Newmark under constant force has the refreshed
    // a(t0+)=1 predictor term plus the solved a1 contribution.
    integrator.trialStep(0.1);
    QVERIFY(std::abs(integrator.candidateState().q[0]-0.005/1.01) < 1.0e-12);
    QVERIFY(std::abs(integrator.candidateState().externalWork-
                     integrator.candidateState().q[0]) < 1.0e-12);
    QVERIFY_EXCEPTION_THROWN(integrator.refreshAcceptedForce(),std::logic_error);
    integrator.rejectStep();
    const auto checkpoint = integrator.checkpoint();
    command = std::numeric_limits<double>::quiet_NaN();
    QVERIFY_EXCEPTION_THROWN(integrator.refreshAcceptedForce(),std::invalid_argument);
    QCOMPARE(integrator.acceptedState().a[0],checkpoint.state.a[0]);
    QCOMPARE(integrator.acceptedState().externalForce[0],checkpoint.state.externalForce[0]);
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

void TestStructuralDynamicsIntegrator::canonicalGridTwentyThousandStepsAndRestore()
{
    RParallelExecutionConfig policy;
    policy.setMode(RParallelExecutionMode::SerialDeterministic);
    const RParallelExecutionContext execution(policy);
    const double time0 = 0.37,nominalDt = 0.0001;
    const auto make = [&] {
        return Integrator(scalarMatrix(1.0),scalarMatrix(0.4),scalarMatrix(4.0),
            solverConf(),{},[time0](double t,const RStructuralDynamicsState &,RRVector &f)
            {f[0] = 0.3*(t-time0);});
    };
    auto uninterrupted = make();
    uninterrupted.initialize(RRVector(1,1.0),RRVector(1,0.2),time0);
    double q = 1.0,v = 0.2,a = -4.08,dissipation = 0.0;
    RStructuralDynamicsCheckpoint split;
    bool exceedsOldStepTolerance = false;
    for (std::uint64_t index=1;index<=20000;++index)
    {
        const double left = index==1?time0:time0+double(index-1)*nominalDt;
        const double right = time0+double(index)*nominalDt,h = right-left;
        exceedsOldStepTolerance = exceedsOldStepTolerance ||
            std::fabs(h-nominalDt)>1.0e-12*std::max(h,nominalDt);
        // Independent scalar Newmark recurrence with the actual endpoint gap.
        const double qPredictor = q+h*v+0.25*h*h*a;
        const double vPredictor = v+0.5*h*a;
        const double nextA = (0.3*(right-time0)-0.4*vPredictor-4.0*qPredictor)/
                             (1.0+0.2*h+h*h);
        const double nextQ = qPredictor+0.25*h*h*nextA;
        const double nextV = vPredictor+0.5*h*nextA;
        dissipation += 0.2*h*(v*v+nextV*nextV);
        q=nextQ; v=nextV; a=nextA;
        uninterrupted.trialStepOnGrid(time0,index,nominalDt);
        uninterrupted.acceptStep();
        const auto &state = uninterrupted.acceptedState();
        QCOMPARE(state.time,right);
        QCOMPARE(state.externalForce[0],0.3*(right-time0));
        if (index==13781) split=uninterrupted.checkpoint();
    }
    QVERIFY(exceedsOldStepTolerance);
    const auto expected = uninterrupted.checkpoint();
    QCOMPARE(expected.fixedStep,nominalDt);
    QVERIFY(std::fabs(expected.state.q[0]-q)<2.0e-11);
    QVERIFY(std::fabs(expected.state.v[0]-v)<2.0e-11);
    QVERIFY(std::fabs(expected.state.a[0]-a)<2.0e-10);
    QVERIFY(std::fabs(expected.state.dampingDissipation-dissipation)<2.0e-11);
    auto resumed=make(); resumed.restore(split);
    for (std::uint64_t index=13782;index<=20000;++index)
    {
        resumed.trialStepOnGrid(time0,index,nominalDt); resumed.acceptStep();
    }
    const auto actual=resumed.checkpoint();
    QCOMPARE(actual.fixedStep,expected.fixedStep);
    QCOMPARE(actual.initialMechanicalEnergy,expected.initialMechanicalEnergy);
    QCOMPARE(actual.state.time,expected.state.time);
    for (const auto &vectors : {std::make_pair(&actual.state.q,&expected.state.q),
                               std::make_pair(&actual.state.v,&expected.state.v),
                               std::make_pair(&actual.state.a,&expected.state.a),
                               std::make_pair(&actual.state.externalForce,&expected.state.externalForce)})
    {
        QCOMPARE(vectors.first->size(),vectors.second->size());
        QVERIFY(std::memcmp(vectors.first->data(),vectors.second->data(),
                            vectors.first->size()*sizeof(double))==0);
    }
    QCOMPARE(actual.state.kineticEnergy,expected.state.kineticEnergy);
    QCOMPARE(actual.state.strainEnergy,expected.state.strainEnergy);
    QCOMPARE(actual.state.externalWork,expected.state.externalWork);
    QCOMPARE(actual.state.dampingDissipation,expected.state.dampingDissipation);
    QCOMPARE(actual.state.energyBalanceResidual,expected.state.energyBalanceResidual);
}

void TestStructuralDynamicsIntegrator::canonicalGridRejectsInvalidClockAndStep()
{
    auto integrator=oscillator(0.0,Integrator::Method::NewmarkAverageAcceleration);
    const double time0=1.0e10,dt=0.01;
    integrator.initialize(RRVector(1,1.0),RRVector(1,0.0),time0);
    for (auto index : {std::uint64_t(0),std::uint64_t(2),std::uint64_t(1)<<54})
        QVERIFY_EXCEPTION_THROWN(integrator.trialStepOnGrid(time0,index,dt),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(integrator.trialStepOnGrid(time0,1,0.0),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(integrator.trialStepOnGrid(time0,1,std::numeric_limits<double>::infinity()),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(integrator.trialStepOnGrid(std::numeric_limits<double>::quiet_NaN(),1,dt),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(integrator.trialStepOnGrid(time0,1,1.0e-20),std::invalid_argument);
    integrator.trialStepOnGrid(time0,1,dt); integrator.acceptStep();
    const auto accepted=integrator.checkpoint(); QCOMPARE(accepted.fixedStep,dt);
    QVERIFY_EXCEPTION_THROWN(integrator.trialStepOnGrid(time0,1,dt),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(integrator.trialStepOnGrid(time0,3,dt),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(integrator.trialStepOnGrid(std::nextafter(time0,std::numeric_limits<double>::infinity()),2,dt),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(integrator.trialStepOnGrid(time0,2,std::nextafter(dt,std::numeric_limits<double>::infinity())),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(integrator.trialStep(dt*2),std::invalid_argument);
    auto wrong=accepted; wrong.fixedStep=0; integrator.restore(wrong);
    QVERIFY_EXCEPTION_THROWN(integrator.trialStepOnGrid(time0,2,dt),std::invalid_argument);
    wrong=accepted; wrong.state.time=std::nextafter(wrong.state.time,std::numeric_limits<double>::infinity()); integrator.restore(wrong);
    QVERIFY_EXCEPTION_THROWN(integrator.trialStepOnGrid(time0,2,dt),std::invalid_argument);
    integrator.restore(accepted); integrator.trialStepOnGrid(time0,2,dt); integrator.rejectStep();
    QCOMPARE(integrator.acceptedState().time,accepted.state.time);
    QCOMPARE(integrator.acceptedState().q,accepted.state.q);
    QVERIFY_EXCEPTION_THROWN(integrator.trialStepOnGrid(std::numeric_limits<double>::max(),2,std::numeric_limits<double>::max()),std::invalid_argument);
}

QTEST_APPLESS_MAIN(TestStructuralDynamicsIntegrator)

#include "tst_structural_dynamics_integrator.moc"

#include <QtTest>
#include <QDataStream>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <rbl_error.h>
#include "rstructural_truss_dynamics.h"

namespace {
using System = RStructuralTrussDynamicsSystem;
using Inputs = RStructuralTrussDynamicsInputs;
using State = RStructuralTrussDynamicsState;
bool near(double a,double b,double tol=2e-9) { return std::abs(a-b)<=tol*(1+std::abs(b)); }
System axial()
{
    System s;
    s.referenceNodes={{{0,0,0}},{{1,0,0}}};
    s.members={{1,0,1,RStructuralTrussPassiveLaw{1,40}}};
    s.nodalLumpedMass={3,2};
    s.prescribedMask={true,true,true,false,true,true};
    s.fixedPrescribedDisplacements.assign(6,0);
    return s;
}
Inputs held(RStructuralTrussInput input=RStructuralTrussNoInput{},double force=0)
{
    Inputs i; i.memberInputs={input}; i.externalForces.assign(6,0); i.externalForces[3]=force; return i;
}
RStructuralTrussDynamicsOptions tight()
{
    RStructuralTrussDynamicsOptions o;
    o.absForceTol=1e-11; o.relForceTol=1e-11; o.linearRelativeTolerance=1e-12;
    return o;
}
State advance(const RStructuralTrussDynamics &d,State s,const Inputs &i,double dt,int n)
{
    for(int j=0;j<n;++j) s=d.trialStep(s,i,dt,tight()).candidate;
    return s;
}
// Serialize every public field (including variant discriminator) to verify the
// accepted object exactly, without comparing uninitialized struct padding.
QByteArray snapshot(const State &s)
{
    QByteArray b; QDataStream out(&b,QIODevice::WriteOnly);
    const auto vector=[&](const std::vector<double>&v) { out<<quint64(v.size()); for(double x:v) out<<x; };
    out<<s.time; vector(s.displacements); vector(s.velocities); vector(s.accelerations);
    vector(s.appliedInputs.externalForces);
    out<<quint64(s.appliedInputs.memberInputs.size());
    for(const auto &i:s.appliedInputs.memberInputs) {
        out<<quint64(i.index());
        std::visit([&](const auto &v) {
            using T=std::decay_t<decltype(v)>;
            if constexpr(!std::is_same_v<T,RStructuralTrussNoInput>) out<<v.value;
        },i);
    }
    vector(s.reactions);
    out<<s.storedEnergy<<s.kineticEnergy<<s.cumulativeExternalWork<<s.cumulativeSupportWork
       <<s.cumulativeElectricalWork<<s.cumulativeLengthWork<<s.initialTotalEnergy<<s.systemIdentity;
    return b;
}
template<class F> bool throws(F f)
{
    try { f(); } catch(const std::exception &) { return true; } catch(const RError &) { return true; }
    return false;
}
RStructuralPrescribedMotion polynomial()
{
    // u=.1t+.2t^2+.05t^3, v=.1+.4t+.15t^2, a=.4+.3t.
    RStructuralPrescribedMotion m;
    m.id="polynomial"; m.target.nodeId=1; m.target.dof=RDofType::UX;
    m.knots={{0,0,.1,.4},{1,.35,.65,.7}};
    return m;
}
System nonlinearSystem()
{
    System s;
    s.referenceNodes={{{-1,0,0}},{{1,0,0}},{{0,.5,0}}};
    s.members={{1,0,2,RStructuralTrussPassiveLaw{1,50}},
               {2,1,2,RStructuralTrussPiezoLaw{1,20,2,2}},
               {3,0,2,RStructuralTrussCableLaw{0,30}}};
    s.nodalLumpedMass={1,1,2}; s.prescribedMask.assign(9,true); s.prescribedMask[7]=false;
    s.fixedPrescribedDisplacements.assign(9,0); return s;
}
Inputs nonlinearInputs()
{
    Inputs i; i.memberInputs={RStructuralTrussNoInput{},RStructuralTrussVoltage{.2},RStructuralTrussCableLength{3}};
    i.externalForces.assign(9,0); return i;
}
// Independently reduced scalar equation: two oblique rods share one free y,
// the cable is slack throughout. No production element/assembly oracle.
double restoring(double y)
{
    const double l=std::sqrt(1+y*y);
    return (70*(l-1)+.4)*y/l;
}
std::array<double,2> rk4(double dt,int n)
{
    double y=.55,v=.07;
    for(int j=0;j<n;++j) {
        const double k1y=v,k1v=-restoring(y)/2;
        const double k2y=v+dt*k1v/2,k2v=-restoring(y+dt*k1y/2)/2;
        const double k3y=v+dt*k2v/2,k3v=-restoring(y+dt*k2y/2)/2;
        const double k4y=v+dt*k3v,k4v=-restoring(y+dt*k3y)/2;
        y+=dt*(k1y+2*k2y+2*k3y+k4y)/6;
        v+=dt*(k1v+2*k2v+2*k3v+k4v)/6;
    }
    return {y,v};
}
}

class TestStructuralTrussDynamics : public QObject
{
    Q_OBJECT
private slots:
    void smallStepAtFiniteDisplacement();
    void axialOscillators_data();
    void axialOscillators();
    void analyticPrescribedMotionAndSupportPower();
    void commandJumpsAndFixedMass();
    void fixedCableSlackTaut();
    void mixedIndefiniteEffectiveTangent();
    void mixedNonlinearRefinement();
    void malformedSystemsAndInputs();
    void staleStateAndFailedTrialAreTransactional();
    void decimalPrescribedEndpointAndExtrapolation();
    void heldChargeJumpAndBoundaryKindRejection();
    void singularLinearCertificateIsTransactional();
    void freeCableTautToSlack();
    void freeMassWithMovingSupport();
    void acceptedStateCertificationReturnsValidatedAssembly();
};

void TestStructuralTrussDynamics::axialOscillators_data()
{
    QTest::addColumn<int>("kind");
    QTest::newRow("passive-constant-load")<<0;
    QTest::newRow("held-voltage-physical-U")<<1;
    QTest::newRow("held-open-charge")<<2;
}
void TestStructuralTrussDynamics::axialOscillators()
{
    QFETCH(int,kind);
    System s=axial(); Inputs i=held({},1.2);
    double k=40,offset=0;
    if(kind) {
        s.members[0].law=RStructuralTrussPiezoLaw{1,40,2,.5};
        if(kind==1) { i.memberInputs[0]=RStructuralTrussVoltage{.7}; offset=1.4; }
        else { i.memberInputs[0]=RStructuralTrussCharge{.15}; k=48; offset=.6; }
    }
    RStructuralTrussDynamics d(s);
    std::vector<double> q(6,0),v(6,0); q[3]=.035; v[3]=.08;
    const State initial=d.initialize(0,q,v,i);
    const double eq=(1.2-offset)/k,w=std::sqrt(k/2),t=.4;
    const double exact=eq+(.035-eq)*std::cos(w*t)+.08/w*std::sin(w*t);
    const auto coarse=advance(d,initial,i,.02,20),fine=advance(d,initial,i,.01,40);
    const double ec=std::abs(coarse.displacements[3]-exact),ef=std::abs(fine.displacements[3]-exact);
    // Newmark phase error is O(dt^2); ratio window admits higher-order terms.
    QVERIFY(ec>1e-8); QVERIFY(ef<.28*ec); QVERIFY(ef<3e-5);
    for(const auto &end:{coarse,fine}) {
        const double u=end.displacements[3];
        const double voltage=kind==1?.7:(kind==2?(.15+2*u)/.5:0);
        const double physicalU=20*u*u+(kind?.25*voltage*voltage:0);
        QVERIFY(near(end.storedEnergy,physicalU,2e-11));
        QVERIFY(near(end.cumulativeExternalWork,1.2*(u-.035),2e-11));
        const double electrical=kind==1?-1.4*(u-.035):0;
        QVERIFY(near(end.cumulativeElectricalWork,electrical,2e-11));
        const double balance=end.kineticEnergy+physicalU-initial.initialTotalEnergy
            -end.cumulativeExternalWork-end.cumulativeElectricalWork;
        QVERIFY(near(balance,0,2e-10));
        QVERIFY(near(end.accelerations[3],(1.2-offset-k*u)/2,2e-10));
    }
}

void TestStructuralTrussDynamics::analyticPrescribedMotionAndSupportPower()
{
    System s=axial(); s.prescribedMask[3]=true; s.prescribedMotions={polynomial()};
    RStructuralTrussDynamics d(s); Inputs i=held({},.3);
    std::vector<double> q(6,0),v(6,0); v[3]=.1;
    const auto initial=d.initialize(0,q,v,i);
    QVERIFY(near(initial.accelerations[3],.4));
    // Binary-exact steps reach the frozen curve endpoint without roundoff extrapolation.
    const auto c=advance(d,initial,i,.125,8),f=advance(d,initial,i,.0625,16);
    // Integral R*v = Delta(T+U)-F*Delta(u), including support inertia.
    const double exact=1*(.65*.65-.1*.1)+20*.35*.35-.3*.35;
    for(const auto &end:{c,f}) {
        QVERIFY(near(end.displacements[3],.35,2e-12));
        QVERIFY(near(end.velocities[3],.65,2e-12));
        QVERIFY(near(end.accelerations[3],.7,2e-12));
        QVERIFY(near(end.reactions[3],2*.7+40*.35-.3,2e-11));
        QVERIFY(near(end.reactions[0],-40*.35,2e-11));
        QVERIFY(near(end.cumulativeExternalWork,.3*.35,2e-11));
    }
    const double ec=std::abs(c.cumulativeSupportWork-exact),ef=std::abs(f.cumulativeSupportWork-exact);
    // Endpoint derivative difference of R*v is 27.2; leading trapezoid
    // error at dt=.0625 is 27.2*dt^2/12=.008854 (higher term is O(dt^4)).
    QVERIFY(ec>1e-5); QVERIFY(ef<.26*ec); QVERIFY(ef<.009);
}

void TestStructuralTrussDynamics::commandJumpsAndFixedMass()
{
    System s=axial(); s.members={{1,0,1,RStructuralTrussPiezoLaw{1,40,2,.5}},
        {2,0,1,RStructuralTrussCableLaw{0,60}}};
    Inputs old=held(RStructuralTrussVoltage{.2},1); old.memberInputs.push_back(RStructuralTrussCableLength{1.2});
    Inputs next=held(RStructuralTrussVoltage{.8},3); next.memberInputs.push_back(RStructuralTrussCableLength{1});
    RStructuralTrussDynamics d(s); std::vector<double> q(6,0),v(6,0); q[3]=.1; v[3]=.05;
    const auto a=d.initialize(0,q,v,old); const auto saved=snapshot(a);
    const auto r=d.refresh(a,next);
    QVERIFY(near(a.accelerations[3],(1-4-.4)/2));
    QVERIFY(near(r.candidate.accelerations[3],(3-4-1.6-6)/2));
    QVERIFY(near(r.jumpWork.electricalWork,.25*(.8*.8-.2*.2)));
    QVERIFY(near(r.jumpWork.lengthWork,.5*60*.1*.1));
    QCOMPARE(r.jumpWork.externalWork,0.0); QCOMPARE(r.jumpWork.supportWork,0.0);
    QVERIFY(near(r.candidate.storedEnergy-a.storedEnergy,r.jumpWork.electricalWork+r.jumpWork.lengthWork));
    QCOMPARE(r.candidate.displacements,a.displacements); QCOMPARE(r.candidate.velocities,a.velocities);
    QCOMPARE(r.candidate.kineticEnergy,a.kineticEnergy); QCOMPARE(snapshot(a),saved);
    const auto trial=d.trialStep(a,next,.005,tight());
    QVERIFY(near(trial.jumpWork.lengthWork,r.jumpWork.lengthWork));
    QVERIFY(near(trial.jumpWork.electricalWork,r.jumpWork.electricalWork));
}

void TestStructuralTrussDynamics::fixedCableSlackTaut()
{
    System s=axial(); s.prescribedMask[3]=true; s.members[0].law=RStructuralTrussCableLaw{0,80};
    RStructuralTrussDynamics d(s); std::vector<double> zero(6,0);
    auto a=d.initialize(0,zero,zero,held(RStructuralTrussCableLength{1.25}));
    auto r=d.refresh(a,held(RStructuralTrussCableLength{.8}));
    // k=100, N=20, U=N^2/(2k)=2.
    QVERIFY(near(r.jumpWork.lengthWork,2)); QVERIFY(near(r.candidate.reactions[3],20));
    auto end=d.trialStep(r.candidate,r.candidate.appliedInputs,.1).candidate;
    QVERIFY(near(end.storedEnergy,2)); QCOMPARE(end.cumulativeSupportWork,0.0);
    r=d.refresh(end,held(RStructuralTrussCableLength{1})); // Exact kink is valid.
    QVERIFY(near(r.jumpWork.lengthWork,-2)); QCOMPARE(r.candidate.storedEnergy,0.0);
    QCOMPARE(r.candidate.reactions[3],0.0);
}

void TestStructuralTrussDynamics::mixedIndefiniteEffectiveTangent()
{
    System s=axial(); s.nodalLumpedMass[1]=.1; s.prescribedMask[4]=false; s.prescribedMask[5]=false;
    s.members={{1,0,1,RStructuralTrussPassiveLaw{2,100}},
        {2,0,1,RStructuralTrussPiezoLaw{1,20,1,2}}, {3,0,1,RStructuralTrussCableLaw{0,30}}};
    Inputs i=held({},-49.9); i.externalForces[4]=.01;
    i.memberInputs={RStructuralTrussNoInput{},RStructuralTrussVoltage{0},RStructuralTrussCableLength{2}};
    RStructuralTrussDynamics d(s); std::vector<double> q(6,0),v(6,0);
    // At q=0, free effective tangent eigenvalues are 70+.4, -50+.4,
    // -50+.4 for dt=1. A positive-curvature-only solver cannot certify it.
    QVERIFY(70+.4>0); QVERIFY(-50+.4<0);
    const auto a=d.initialize(0,q,v,i); const auto r=d.trialStep(a,i,1,tight());
    const auto &e=r.candidate;
    double l=0; for(int j=0;j<3;++j) l+=std::pow((j==0?1:0)+e.displacements[3+j],2);
    l=std::sqrt(l);
    const double n=50*(l-2)+20*(l-1); // Cable remains slack (l<2).
    QVERIFY(l<2); QVERIFY(r.iterations>0); QVERIFY(r.residualNorm<1e-8);
    for(int j=0;j<3;++j) {
        const double internal=n*((j==0?1:0)+e.displacements[3+j])/l;
        QVERIFY(near(.1*e.accelerations[3+j]+internal-i.externalForces[3+j],0,1e-8));
        QVERIFY(near(e.reactions[3+j],0,1e-8));
    }
    QVERIFY(near(e.accelerations[3],4*e.displacements[3]-a.accelerations[3],1e-8));
}

void TestStructuralTrussDynamics::mixedNonlinearRefinement()
{
    const auto s=nonlinearSystem(); const auto i=nonlinearInputs(); RStructuralTrussDynamics d(s);
    std::vector<double> q(9,0),v(9,0); q[7]=.05; v[7]=.07;
    const auto initial=d.initialize(0,q,v,i);
    const auto reference=rk4(.0001,4000),check=rk4(.0002,2000);
    QVERIFY(near(reference[0],check[0],1e-12)); QVERIFY(near(reference[1],check[1],1e-12));
    const auto c=advance(d,initial,i,.02,20),f=advance(d,initial,i,.01,40);
    const double ec=std::abs(.5+c.displacements[7]-reference[0]),ef=std::abs(.5+f.displacements[7]-reference[0]);
    // O(dt^2) nonlinear integration; reference RK4 error independently bounded above.
    QVERIFY(ec>1e-7); QVERIFY(ef<.3*ec); QVERIFY(ef<1e-4);
    QVERIFY(near(2*f.accelerations[7]+restoring(.5+f.displacements[7]),0,1e-8));
    const double bc=std::abs(c.kineticEnergy+c.storedEnergy-initial.initialTotalEnergy-c.cumulativeElectricalWork);
    const double bf=std::abs(f.kineticEnergy+f.storedEnergy-initial.initialTotalEnergy-f.cumulativeElectricalWork);
    QVERIFY(bc>1e-8); QVERIFY(bf<.35*bc);
}

void TestStructuralTrussDynamics::malformedSystemsAndInputs()
{
    const auto base=axial(); std::vector<double> zero(6,0);
    for(double mass:{0.0,-1.0,std::numeric_limits<double>::quiet_NaN()}) {
        auto s=base; s.nodalLumpedMass[1]=mass; QVERIFY(throws([&]{ RStructuralTrussDynamics d(s); }));
    }
    { auto s=base; s.members.push_back(s.members[0]); QVERIFY(throws([&]{ RStructuralTrussDynamics d(s); })); }
    { auto s=base; s.fixedPrescribedDisplacements[3]=.1; QVERIFY(throws([&]{ RStructuralTrussDynamics d(s); })); }
    for(int kind=0;kind<7;++kind) {
        auto s=base; s.prescribedMask[3]=true; auto m=polynomial();
        if(kind==0) s.prescribedMask[3]=false;
        if(kind==1) m.target.dof=RDofType::RX;
        if(kind==2) m.target.nodeId=2;
        if(kind==3) m.knots[1].time=0;
        s.prescribedMotions={m};
        if(kind>=4) {
            auto duplicate=m;
            if(kind==5) duplicate.id="other-id-same-dof";
            if(kind==6) duplicate.target.dof=RDofType::UY;
            s.prescribedMotions.push_back(duplicate);
        }
        QVERIFY(throws([&]{ RStructuralTrussDynamics d(s); }));
    }
    RStructuralTrussDynamics d(base);
    { auto i=held(); i.externalForces.pop_back(); QVERIFY(throws([&]{d.initialize(0,zero,zero,i);})); }
    { auto i=held(RStructuralTrussVoltage{0}); QVERIFY(throws([&]{d.initialize(0,zero,zero,i);})); }
    { auto q=zero; q[3]=-1; QVERIFY(throws([&]{d.initialize(0,q,zero,held());})); }
    QVERIFY(throws([&]{d.initialize(std::numeric_limits<double>::quiet_NaN(),zero,zero,held());}));
    auto s=base; s.members[0].law=RStructuralTrussCableLaw{0,40}; RStructuralTrussDynamics cable(s);
    for(double l:{0.0,-1.0,std::numeric_limits<double>::quiet_NaN()})
        QVERIFY(throws([&]{cable.initialize(0,zero,zero,held(RStructuralTrussCableLength{l}));}));
    s=base; s.prescribedMask[3]=true; s.prescribedMotions={polynomial()}; RStructuralTrussDynamics moving(s);
    auto v=zero; v[3]=.1;
    QVERIFY(throws([&]{moving.initialize(-.1,zero,v,held());}));
    QVERIFY(throws([&]{moving.initialize(0,zero,zero,held());}));
    const auto a=moving.initialize(0,zero,v,held());
    QVERIFY(throws([&]{moving.trialStep(a,held(),1.1);}));
}

void TestStructuralTrussDynamics::staleStateAndFailedTrialAreTransactional()
{
    const auto s=nonlinearSystem(); const auto i=nonlinearInputs(); RStructuralTrussDynamics d(s);
    std::vector<double> q(9,0),v(9,0); q[7]=.05; v[7]=.07;
    const auto accepted=d.initialize(0,q,v,i); const auto bytes=snapshot(accepted);
    auto options=tight(); options.maxIterations=1;
    QVERIFY(throws([&]{d.trialStep(accepted,i,.4,options);}));
    QCOMPARE(snapshot(accepted),bytes);
    const auto retry=d.trialStep(accepted,i,.01,tight());
    QVERIFY(retry.residualNorm<1e-8); QCOMPARE(snapshot(accepted),bytes);
    auto stale=accepted; stale.accelerations[7]+=.1;
    QVERIFY(throws([&]{d.trialStep(stale,i,.01);}));
    stale=accepted; stale.appliedInputs.memberInputs[1]=RStructuralTrussVoltage{.8};
    QVERIFY(throws([&]{d.refresh(stale,i);}));
    auto other=s; other.nodalLumpedMass[2]=3; RStructuralTrussDynamics wrongMass(other);
    QVERIFY(throws([&]{wrongMass.trialStep(accepted,i,.01);}));
    other=s; std::get<RStructuralTrussPiezoLaw>(other.members[1].law).coupling=3;
    RStructuralTrussDynamics wrongLaw(other); QVERIFY(throws([&]{wrongLaw.refresh(accepted,i);}));
    for(double dt:{0.0,-.01,std::numeric_limits<double>::quiet_NaN()})
        QVERIFY(throws([&]{d.trialStep(accepted,i,dt);}));
    QCOMPARE(snapshot(accepted),bytes);
}

void TestStructuralTrussDynamics::decimalPrescribedEndpointAndExtrapolation()
{
    System s=axial(); s.prescribedMask[3]=true; s.prescribedMotions={polynomial()};
    RStructuralTrussDynamics d(s); const auto i=held({},.3);
    std::vector<double> q(6,0),v(6,0); v[3]=.1;
    // Decimal accumulation typically produces 1+several ulps after 100 steps.
    // Curve endpoint tolerance must accept that roundoff, but not extrapolation.
    State end=d.initialize(0,q,v,i);
    for(int j=0;j<100;++j) end=d.trialStep(end,i,.01,tight()).candidate;
    QVERIFY(near(end.time,1,2e-14));
    QVERIFY(near(end.displacements[3],.35,2e-12));
    QVERIFY(near(end.velocities[3],.65,2e-12));
    QVERIFY(near(end.accelerations[3],.7,2e-12));
    QVERIFY(near(end.reactions[3],15.1,2e-11));
    const auto bytes=snapshot(end);
    QVERIFY(throws([&]{ d.trialStep(end,i,.01,tight()); }));
    QCOMPARE(snapshot(end),bytes);
}

void TestStructuralTrussDynamics::heldChargeJumpAndBoundaryKindRejection()
{
    System s=axial(); s.members[0].law=RStructuralTrussPiezoLaw{1,40,2,.5};
    RStructuralTrussDynamics d(s);
    const auto old=held(RStructuralTrussCharge{.15},1),next=held(RStructuralTrussCharge{.35},1);
    std::vector<double> q(6,0),v(6,0); q[3]=.1; v[3]=.05;
    const auto a=d.initialize(0,q,v,old); const auto bytes=snapshot(a);
    // Fixed delta=.1: V=(Q+g*delta)/Ce gives .7 -> 1.1;
    // N=k*delta+g*V gives 5.4 -> 6.2, U=.2+.25*V^2.
    QVERIFY(near(a.storedEnergy,.3225)); QVERIFY(near(a.accelerations[3],-2.2));
    const auto r=d.refresh(a,next);
    QVERIFY(near(r.candidate.storedEnergy,.5025));
    QVERIFY(near(r.candidate.accelerations[3],-2.6));
    QVERIFY(near(r.candidate.reactions[0],-6.2));
    QVERIFY(near(r.jumpWork.electricalWork,.18));
    QCOMPARE(r.jumpWork.lengthWork,0.0); QCOMPARE(r.jumpWork.externalWork,0.0);
    QCOMPARE(r.candidate.kineticEnergy,a.kineticEnergy);
    QCOMPARE(r.candidate.displacements,a.displacements); QCOMPARE(r.candidate.velocities,a.velocities);
    const auto switched=held(RStructuralTrussVoltage{1.1},1);
    QVERIFY(throws([&]{d.refresh(a,switched);}));
    QVERIFY(throws([&]{d.trialStep(a,switched,.01,tight());}));
    // The reverse switch is equally unavailable despite matching V and Q.
    const auto voltageState=d.initialize(0,q,v,held(RStructuralTrussVoltage{.7},1));
    QVERIFY(throws([&]{d.refresh(voltageState,old);}));
    QCOMPARE(snapshot(a),bytes);
}

void TestStructuralTrussDynamics::singularLinearCertificateIsTransactional()
{
    System s=axial(); s.members[0].law=RStructuralTrussPassiveLaw{2,100};
    s.nodalLumpedMass[1]=3.125; s.prescribedMask[4]=false;
    RStructuralTrussDynamics d(s); auto i=held({},-50); i.externalForces[4]=1;
    std::vector<double> q(6,0),v(6,0); v[4]=-.04;
    const auto a=d.initialize(0,q,v,i); const auto bytes=snapshot(a);
    QVERIFY(near(a.accelerations[4],.32));
    // dt=.5: predictor_y=.5*(-.04)+(.5^2/4)*.32=0.
    // N/l=-50, 4M/dt^2=50: K_eff=diag(100,0), rhs=(0,1).
    // Inconsistent linear equations require rejection before any candidate commit.
    bool rejected=false;
    try { d.trialStep(a,i,.5,tight()); }
    catch(const std::runtime_error &e) {
        rejected=true;
        QVERIFY2(QString::fromUtf8(e.what()).contains("linear residual certification"),e.what());
    }
    QVERIFY(rejected); QCOMPARE(snapshot(a),bytes);
}

void TestStructuralTrussDynamics::freeCableTautToSlack()
{
    System s=axial(); s.members[0].law=RStructuralTrussCableLaw{0,80};
    RStructuralTrussDynamics d(s); const auto i=held(RStructuralTrussCableLength{1});
    std::vector<double> q(6,0),v(6,0); q[3]=.1; v[3]=-.3;
    const auto initial=d.initialize(0,q,v,i);
    const double w=std::sqrt(40.0),tc=std::atan2(.1*w,.3)/w;
    // Before tc: .1*cos(wt)-.3/w*sin(wt). At tc speed is -.7;
    // afterwards the tension vanishes and u=-.7*(t-tc), v=-.7.
    const double exact=-.7*(.4-tc);
    QVERIFY(tc>0 && tc<.4); QVERIFY(1+exact>0);
    QVERIFY(near(initial.initialTotalEnergy,.49));
    const auto c=advance(d,initial,i,.02,20);
    const auto f=advance(d,initial,i,.005,80);
    const auto ff=advance(d,initial,i,.0025,160);
    for(const auto &end:{c,f,ff}) {
        QVERIFY(end.displacements[3]<0); QVERIFY(1+end.displacements[3]>0);
        QCOMPARE(end.storedEnergy,0.0); QCOMPARE(end.accelerations[3],0.0);
        QCOMPARE(end.reactions[0],0.0); QCOMPARE(end.cumulativeLengthWork,0.0);
        QCOMPARE(end.cumulativeExternalWork,0.0); QCOMPARE(end.cumulativeSupportWork,0.0);
        QCOMPARE(end.cumulativeElectricalWork,0.0);
        QVERIFY(end.kineticEnergy>.49);
    }
    const double ec=std::abs(c.displacements[3]-exact);
    const double ef=std::abs(f.displacements[3]-exact);
    const double eff=std::abs(ff.displacements[3]-exact);
    // Kink crossing has grid-dependent constants. Eightfold dt refinement
    // must reduce analytic error decisively; no smooth global order is assumed.
    QVERIFY(ec>1e-6); QVERIFY(ef<.5*ec); QVERIFY(eff<.15*ec);
    QVERIFY(eff<3e-5); QVERIFY(near(ff.velocities[3],-.7,5e-5));
    const double bc=std::abs(c.kineticEnergy-.49),bf=std::abs(f.kineticEnergy-.49);
    const double bff=std::abs(ff.kineticEnergy-.49);
    // Only the crossing step gains energy: DeltaE=-(EA/2L0)*u_left*u_right.
    // This independent kink identity gives O(dt^2), with phase-dependent factor.
    QVERIFY(bc>1e-6); QVERIFY(bf<.5*bc); QVERIFY(bff<.15*bc); QVERIFY(bff<5e-5);
}

void TestStructuralTrussDynamics::freeMassWithMovingSupport()
{
    System s=axial(); RStructuralPrescribedMotion m;
    m.id="quadratic-anchor"; m.target.nodeId=0; m.target.dof=RDofType::UX;
    // g=A*t^2, A=.05: exact quintic Hermite representation of a quadratic.
    m.knots={{0,0,0,.1},{1,.05,.1,.1}}; s.prescribedMotions={m};
    RStructuralTrussDynamics d(s); const auto i=held(); std::vector<double> z(6,0);
    const auto initial=d.initialize(0,z,z,i);
    const double t=.5,w=std::sqrt(20.0),g=.05*t*t,vg=.1*t;
    // m*x''+k*(x-g)=0 with x(0)=x'(0)=0.
    const double x=g-.1/(w*w)*(1-std::cos(w*t));
    const double velocity=.1*t-.1/w*std::sin(w*t);
    const double acceleration=.1*(1-std::cos(w*t));
    const double exactWork=velocity*velocity+20*(x-g)*(x-g)+1.5*vg*vg;
    const auto c=advance(d,initial,i,.025,20),f=advance(d,initial,i,.0125,40);
    const double ec=std::abs(c.displacements[3]-x),ef=std::abs(f.displacements[3]-x);
    QVERIFY(ec>1e-8); QVERIFY(ef<.28*ec); QVERIFY(ef<1e-5);
    QVERIFY(near(f.velocities[3],velocity,2e-5));
    QVERIFY(near(f.accelerations[3],acceleration,2e-4));
    for(const auto &end:{c,f}) {
        QVERIFY(near(end.displacements[0],g,2e-12));
        QVERIFY(near(end.velocities[0],vg,2e-12));
        QVERIFY(near(end.accelerations[0],.1,2e-12));
        QVERIFY(near(end.reactions[0],.3-40*(end.displacements[3]-g),2e-10));
        QVERIFY(near(end.reactions[3],0,1e-8));
    }
    const double wc=std::abs(c.cumulativeSupportWork-exactWork);
    const double wf=std::abs(f.cumulativeSupportWork-exactWork);
    // Quadrature plus Newmark error are O(dt^2) for this analytic smooth path.
    QVERIFY(wc>1e-8); QVERIFY(wf<.3*wc); QVERIFY(wf<1e-5);
}

void TestStructuralTrussDynamics::smallStepAtFiniteDisplacement()
{
    const auto system=axial();RStructuralTrussDynamics driver(system);
    const auto inputs=held(RStructuralTrussNoInput{});
    std::vector<double> q(6,0.0),v(6,0.0);q[3]=.05;v[3]=.07;
    auto state=driver.initialize(0,q,v,inputs);
    constexpr double dt=1e-6;
    for (int step=0;step<100;++step)
    {
        const auto trial=driver.trialStep(state,inputs,dt,tight());
        QVERIFY(trial.residualNorm<=tight().absForceTol+tight().relForceTol*2.1);
        state=trial.candidate;
    }
    const double omega=std::sqrt(20.0),time=100*dt;
    const double expectedQ=.05*std::cos(omega*time)+.07/omega*std::sin(omega*time);
    const double expectedV=-.05*omega*std::sin(omega*time)+.07*std::cos(omega*time);
    QVERIFY(std::abs(state.displacements[3]-expectedQ)<1e-13);
    QVERIFY(std::abs(state.velocities[3]-expectedV)<1e-12);
}

void TestStructuralTrussDynamics::acceptedStateCertificationReturnsValidatedAssembly()
{
    const auto system=axial();RStructuralTrussDynamics driver(system);
    std::vector<double> q(6,0.0),v(6,0.0);q[3]=.05;v[3]=.07;
    auto state=driver.initialize(0,q,v,held());
    state=advance(driver,state,held(),.001,20);
    const auto before=snapshot(state);
    const auto certified=driver.certifyAcceptedState(state);
    const auto expected=RStructuralTrussAssembly::evaluate(system.referenceNodes,state.displacements,
                                                           system.members,state.appliedInputs.memberInputs);
    QCOMPARE(certified.storedEnergy,expected.storedEnergy);
    QCOMPARE(certified.internalForces.size(),expected.internalForces.size());
    for(std::size_t i=0;i<expected.internalForces.size();++i)
        QCOMPARE(certified.internalForces[i],expected.internalForces[i]);
    QCOMPARE(snapshot(state),before);
    auto corrupt=state;corrupt.reactions[3]+=1.0;
    const auto corruptBefore=snapshot(corrupt);
    QVERIFY(throws([&]{driver.certifyAcceptedState(corrupt);}));
    QCOMPARE(snapshot(corrupt),corruptBefore);
}

QTEST_APPLESS_MAIN(TestStructuralTrussDynamics)
#include "tst_structural_truss_dynamics.moc"

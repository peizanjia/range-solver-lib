#include <QtTest>
#include <QDataStream>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <omp.h>
#include "rparallel_execution_context.h"
#include "rstructural_truss_dynamics.h"

namespace {
using System=RStructuralTrussDynamicsSystem;
using Inputs=RStructuralTrussDynamicsInputs;
using State=RStructuralTrussDynamicsState;
using Driver=RStructuralTrussDynamics;
int openMPThreadLimit() {
#if defined(_OPENMP) && _OPENMP >= 200805
    return omp_get_thread_limit();
#else
    return std::numeric_limits<int>::max();
#endif
}
bool near(double a,double b,double tol=2e-10) {return std::abs(a-b)<=tol*(1+std::abs(b));}
System axial(double c=.8) {
    System s; s.referenceNodes={{{0,0,0}},{{1,0,0}}};
    s.members={{1,0,1,RStructuralTrussPassiveLaw{1,40}}};
    s.nodalLumpedMass={3,2}; s.prescribedMask={true,true,true,false,true,true};
    s.fixedPrescribedDisplacements.assign(6,0);
    s.physicalDamping.assign(36,0); s.physicalDamping[3*6+3]=c; return s;
}
Inputs force(double f=0) {Inputs i; i.memberInputs={RStructuralTrussNoInput{}}; i.externalForces.assign(6,0); i.externalForces[3]=f; return i;}
RStructuralTrussDynamicsOptions tight() {RStructuralTrussDynamicsOptions o; o.absForceTol=1e-11; o.relForceTol=1e-11; o.linearRelativeTolerance=1e-12; return o;}
QByteArray snapshot(const State &s) {
    QByteArray b; QDataStream out(&b,QIODevice::WriteOnly);
    const auto vector=[&](const std::vector<double>&v) {out<<quint64(v.size()); for(double x:v) out<<x;};
    out<<s.time; vector(s.displacements); vector(s.velocities); vector(s.accelerations);
    vector(s.reactions); vector(s.appliedInputs.externalForces); vector(s.appliedInputs.rightExternalForces);
    out<<quint64(s.appliedInputs.memberInputs.size());
    for(const auto &i:s.appliedInputs.memberInputs) {out<<quint64(i.index()); std::visit([&](const auto &v){
        using T=std::decay_t<decltype(v)>; if constexpr(!std::is_same_v<T,RStructuralTrussNoInput>) out<<v.value;
    },i);}
    out<<s.storedEnergy<<s.kineticEnergy<<s.cumulativeExternalWork<<s.cumulativeSupportWork
       <<s.cumulativeElectricalWork<<s.cumulativeLengthWork<<s.cumulativeDampingDissipation
       <<s.initialTotalEnergy<<s.systemIdentity; return b;
}
// Independent closed-form m*x''+c*x'+k*x=f(t), including the particular
// solution's initial displacement/velocity in the homogeneous correction.
double load(int kind,double t) {return kind==0?1.2:kind==1?.3+1.7*t:1.1*std::sin(3*t);}
std::array<double,2> exact(int kind,double t) {
    constexpr double m=2,c=.8,k=40,omega=3;
    const double alpha=c/(2*m),wd=std::sqrt(k/m-alpha*alpha);
    double p=0,pv=0,p0=0,pv0=0;
    if(kind==0) p=p0=1.2/k;
    if(kind==1) {p=(.3+1.7*t)/k-c*1.7/(k*k); pv=pv0=1.7/k; p0=.3/k-c*1.7/(k*k);}
    if(kind==2) {
        const double a=k-m*omega*omega,b=c*omega,den=a*a+b*b;
        p=1.1*(a*std::sin(omega*t)-b*std::cos(omega*t))/den;
        pv=1.1*omega*(a*std::cos(omega*t)+b*std::sin(omega*t))/den;
        p0=-1.1*b/den; pv0=1.1*omega*a/den;
    }
    const double A=.035-p0,B=(.08-pv0+alpha*A)/wd,e=std::exp(-alpha*t);
    const double h=A*std::cos(wd*t)+B*std::sin(wd*t);
    return {p+e*h,pv+e*(-alpha*h+wd*(-A*std::sin(wd*t)+B*std::cos(wd*t)))};
}
}

class TestStructuralTrussDampedForcing:public QObject {
    Q_OBJECT
private slots:
    void oscillatorEndpointLoads_data();
    void oscillatorEndpointLoads();
    void fullCrossDampingAndMovingSupport();
    void localPsdAndIdentity();
    void endpointAndAcceptedStateValidation();
    void nonlinearHeldCommandJump();
    void canonicalTwentyThousandAndSplit();
    void failurePreservesAccepted();
    void coupledSolveIsBitwiseDeterministicAndRestoresCallerThreads();
};

void TestStructuralTrussDampedForcing::oscillatorEndpointLoads_data() {
    QTest::addColumn<int>("kind"); QTest::newRow("constant")<<0; QTest::newRow("linear")<<1; QTest::newRow("sine")<<2;
}
void TestStructuralTrussDampedForcing::oscillatorEndpointLoads() {
    QFETCH(int,kind); Driver d(axial()); std::vector<double> q(6,0),v(6,0); q[3]=.035; v[3]=.08;
    const auto initial=d.initialize(0,q,v,force(load(kind,0))); const auto reference=exact(kind,.4);
    std::array<double,2> error{};
    for(int level=0;level<2;++level) {
        const double dt=level?.005:.01; const int count=level?80:40; auto s=initial; double work=0,loss=0;
        for(int j=1;j<=count;++j) {
            auto i=force(load(kind,(j-1)*dt)); i.rightExternalForces=force(load(kind,j*dt)).externalForces;
            const auto r=d.trialStepOnGrid(s,i,0,std::uint64_t(j),dt,tight()); const auto &n=r.candidate;
            const double dq=n.displacements[3]-s.displacements[3],mean=.5*(n.velocities[3]+s.velocities[3]);
            work+=.5*(i.externalForces[3]+i.rightExternalForces[3])*dq;
            loss+=(n.time-s.time)*.8*mean*mean;
            QVERIFY(near(r.stepWork.externalWork,.5*(i.externalForces[3]+i.rightExternalForces[3])*dq));
            QVERIFY(near(r.stepWork.dampingDissipation,(n.time-s.time)*.8*mean*mean));
            QVERIFY(near(n.accelerations[3],(i.rightExternalForces[3]-.8*n.velocities[3]-40*n.displacements[3])/2));
            QVERIFY(n.appliedInputs.rightExternalForces.empty()); QCOMPARE(n.appliedInputs.externalForces,i.rightExternalForces);
            QVERIFY(std::abs(r.energyBalanceError)<1e-10); s=n;
        }
        QVERIFY(near(s.cumulativeExternalWork,work)); QVERIFY(near(s.cumulativeDampingDissipation,loss)); QVERIFY(loss>0);
        error[level]=std::hypot(s.displacements[3]-reference[0],s.velocities[3]-reference[1]);
    }
    QVERIFY(error[0]>1e-8); QVERIFY(error[1]<.28*error[0]); QVERIFY(error[1]<1e-5);
}

void TestStructuralTrussDampedForcing::fullCrossDampingAndMovingSupport() {
    auto s=axial(0); s.physicalDamping[0]=s.physicalDamping[21]=1; s.physicalDamping[3]=s.physicalDamping[18]=-1;
    RStructuralPrescribedMotion motion; motion.id="constant-speed-anchor"; motion.target.nodeId=0; motion.target.dof=RDofType::UX;
    motion.knots={{0,0,.1,0},{1,.1,.1,0}}; s.prescribedMotions={motion}; Driver d(s);
    std::vector<double> q(6,0),v(6,0); v[0]=.1; const auto initial=d.initialize(0,q,v,force());
    QVERIFY(near(initial.accelerations[3],.05)); QVERIFY(near(initial.reactions[0],.1));
    const double t=.4,alpha=.25,wd=std::sqrt(20-alpha*alpha),e=std::exp(-alpha*t);
    const double x=.1*t-.1/wd*e*std::sin(wd*t),velocity=.1-.1*e*(std::cos(wd*t)-alpha/wd*std::sin(wd*t));
    double ec=0,ef=0;
    for(int level=0;level<2;++level) {
        auto n=initial; const double dt=level?.005:.01; const int count=level?80:40;
        for(int j=1;j<=count;++j) {
            auto r=d.trialStepOnGrid(n,force(),0,std::uint64_t(j),dt,tight()); const auto &right=r.candidate;
            const double rel=.5*(n.velocities[0]+right.velocities[0]-n.velocities[3]-right.velocities[3]);
            QVERIFY(near(r.stepWork.dampingDissipation,dt*rel*rel));
            QVERIFY(near(right.reactions[0],40*(right.displacements[0]-right.displacements[3])+right.velocities[0]-right.velocities[3]));
            QVERIFY(near(2*right.accelerations[3]+40*(right.displacements[3]-right.displacements[0])+right.velocities[3]-right.velocities[0],0));
            QVERIFY(std::abs(r.energyBalanceError)<1e-10); n=right;
        }
        const double err=std::hypot(n.displacements[3]-x,n.velocities[3]-velocity); if(level) ef=err; else ec=err;
        QVERIFY(n.cumulativeDampingDissipation>0); QVERIFY(n.cumulativeSupportWork>0);
    }
    QVERIFY(ec>1e-8); QVERIFY(ef<.28*ec);
}

void TestStructuralTrussDampedForcing::localPsdAndIdentity() {
    auto s=axial(); s.physicalDamping.clear(); Driver legacy(s); auto zero=s; zero.physicalDamping.assign(36,0); Driver explicitZero(zero);
    QVERIFY(legacy.systemIdentity()!=explicitZero.systemIdentity());
    auto expectInvalid=[&](const System &bad){QVERIFY_EXCEPTION_THROWN(Driver{bad},std::invalid_argument);};
    auto bad=zero; bad.physicalDamping.pop_back(); expectInvalid(bad);
    bad=zero; bad.physicalDamping[0]=1e100; bad.physicalDamping[21]=-1e-100; expectInvalid(bad);
    bad=zero; bad.physicalDamping[3]=bad.physicalDamping[18]=1e-100; expectInvalid(bad);
    bad=zero; bad.physicalDamping[0]=bad.physicalDamping[21]=1e-100; bad.physicalDamping[3]=bad.physicalDamping[18]=2e-100; expectInvalid(bad);
    bad=zero; bad.physicalDamping[0]=bad.physicalDamping[21]=1; bad.physicalDamping[3]=.1; expectInvalid(bad);
    bad=zero; bad.physicalDamping[0]=std::numeric_limits<double>::infinity(); expectInvalid(bad);
    auto good=zero; good.physicalDamping[0]=1e100; good.physicalDamping[21]=1e-100; Driver scaled(good);
    QVERIFY(scaled.systemIdentity()!=explicitZero.systemIdentity());
    auto changed=good; changed.physicalDamping[21]=2e-100; Driver changedDriver(changed); QVERIFY(scaled.systemIdentity()!=changedDriver.systemIdentity());
}

void TestStructuralTrussDampedForcing::endpointAndAcceptedStateValidation() {
    Driver d(axial()); std::vector<double> z(6,0); const auto a=d.initialize(0,z,z,force()); auto inputs=force(1);
    inputs.rightExternalForces=force(2).externalForces;
    QVERIFY_EXCEPTION_THROWN(d.initialize(0,z,z,inputs),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(d.refresh(a,inputs),std::invalid_argument);
    auto bad=a; bad.appliedInputs.rightExternalForces=force(0).externalForces;
    QVERIFY_EXCEPTION_THROWN(d.trialStep(bad,force(),.01),std::invalid_argument);
    bad=a; bad.cumulativeDampingDissipation=-1; QVERIFY_EXCEPTION_THROWN(d.refresh(bad,force()),std::invalid_argument);
    bad=a; bad.cumulativeDampingDissipation=std::numeric_limits<double>::quiet_NaN(); QVERIFY_EXCEPTION_THROWN(d.refresh(bad,force()),std::invalid_argument);
    inputs.rightExternalForces.pop_back(); QVERIFY_EXCEPTION_THROWN(d.trialStep(a,inputs,.01),std::invalid_argument);
    inputs.rightExternalForces=force().externalForces; inputs.rightExternalForces[0]=std::numeric_limits<double>::quiet_NaN();
    QVERIFY_EXCEPTION_THROWN(d.trialStep(a,inputs,.01),std::invalid_argument);
    auto f=force(2); const auto refreshed=d.refresh(a,f); QVERIFY(near(refreshed.candidate.accelerations[3],1));
    auto both=f; both.rightExternalForces=force(-1).externalForces;
    const auto r=d.trialStep(a,both,.01,tight()); const auto &n=r.candidate;
    QVERIFY(near(n.displacements[3],.01*.01*.25*(1+n.accelerations[3])));
    QVERIFY(near(r.stepWork.externalWork,.5*(2-1)*n.displacements[3]));
}

void TestStructuralTrussDampedForcing::nonlinearHeldCommandJump() {
    System s; s.referenceNodes={{{-1,0,0}},{{1,0,0}},{{0,.5,0}}};
    s.members={{1,0,2,RStructuralTrussPassiveLaw{1,50}},{2,1,2,RStructuralTrussPiezoLaw{1,20,2,2}},{3,0,2,RStructuralTrussCableLaw{0,30}}};
    s.nodalLumpedMass={1,1,2}; s.prescribedMask.assign(9,true); s.prescribedMask[7]=false; s.fixedPrescribedDisplacements.assign(9,0);
    s.physicalDamping.assign(81,0); s.physicalDamping[70]=.6; Driver d(s);
    Inputs old; old.memberInputs={RStructuralTrussNoInput{},RStructuralTrussVoltage{.2},RStructuralTrussCableLength{3}}; old.externalForces.assign(9,0);
    std::vector<double> q(9,0),v(9,0); q[7]=.05; v[7]=.07; const auto a=d.initialize(0,q,v,old); const auto saved=snapshot(a);
    auto next=old; next.memberInputs[1]=RStructuralTrussVoltage{.8}; next.memberInputs[2]=RStructuralTrussCableLength{1}; next.externalForces[7]=1;
    const auto jump=d.refresh(a,next); QVERIFY(jump.jumpWork.electricalWork>0); QVERIFY(jump.jumpWork.lengthWork>0); QCOMPARE(jump.jumpWork.dampingDissipation,0.0);
    next.rightExternalForces=next.externalForces; next.rightExternalForces[7]=2;
    const auto r=d.trialStep(a,next,.01,tight()); const auto &n=r.candidate; const double y=.5+n.displacements[7],l=std::hypot(1.,y);
    const double restoring=(70*(l-1)+1.6+30*(l-1))*y/l;
    QVERIFY(near(2*n.accelerations[7]+.6*n.velocities[7]+restoring-2,0,2e-9));
    QVERIFY(near(r.jumpWork.electricalWork,jump.jumpWork.electricalWork)); QVERIFY(near(r.jumpWork.lengthWork,jump.jumpWork.lengthWork));
    QCOMPARE(snapshot(a),saved); QVERIFY(n.appliedInputs.rightExternalForces.empty());
    // Independently reduced oblique mixed rod/cable equation, integrated by
    // RK4 at two fine grids to certify the reference separately from Newmark.
    const auto rk4=[](double h,int count) {
        double y=.55,v=.07;
        const auto acceleration=[](double t,double y,double v) {
            const double l=std::hypot(1.,y);
            return (1+2.5*t-.6*v-(100*(l-1)+1.6)*y/l)/2;
        };
        for(int j=0;j<count;++j) {
            const double t=j*h,ky1=v,kv1=acceleration(t,y,v);
            const double ky2=v+.5*h*kv1,kv2=acceleration(t+.5*h,y+.5*h*ky1,v+.5*h*kv1);
            const double ky3=v+.5*h*kv2,kv3=acceleration(t+.5*h,y+.5*h*ky2,v+.5*h*kv2);
            const double ky4=v+h*kv3,kv4=acceleration(t+h,y+h*ky3,v+h*kv3);
            y+=h*(ky1+2*ky2+2*ky3+ky4)/6; v+=h*(kv1+2*kv2+2*kv3+kv4)/6;
        }
        return std::array<double,2>{y,v};
    };
    const auto reference=rk4(.0001,4000),check=rk4(.0002,2000);
    QVERIFY(near(reference[0],check[0],1e-12)); QVERIFY(near(reference[1],check[1],1e-12));
    std::array<double,2> errors{},balance{};
    for(int level=0;level<2;++level) {
        const double h=level?.005:.01; const int count=level?80:40; auto state=a;
        for(int j=1;j<=count;++j) {
            auto command=next; command.externalForces[7]=1+2.5*(j-1)*h;
            command.rightExternalForces[7]=1+2.5*j*h;
            const auto result=d.trialStepOnGrid(state,command,0,std::uint64_t(j),h,tight());
            state=result.candidate; balance[level]=std::abs(result.energyBalanceError);
        }
        errors[level]=std::hypot(.5+state.displacements[7]-reference[0],state.velocities[7]-reference[1]);
    }
    QVERIFY(errors[0]>1e-8); QVERIFY(errors[1]<.3*errors[0]); QVERIFY(errors[1]<1e-4);
    QVERIFY(balance[0]>1e-9); QVERIFY(balance[1]<.3*balance[0]);
}

void TestStructuralTrussDampedForcing::canonicalTwentyThousandAndSplit() {
    System s; s.referenceNodes={{{0,0,0}}}; s.nodalLumpedMass={2}; s.prescribedMask={false,true,true}; s.fixedPrescribedDisplacements.assign(3,0);
    s.physicalDamping.assign(9,0); s.physicalDamping[0]=.4; Driver d(s),restarted(s);
    Inputs i; i.externalForces={1,0,0}; std::vector<double> z(3,0); constexpr double t0=1e9,dt=1e-5;
    State current=d.initialize(t0,z,z,i),mid; double q=0,v=0,a=.5; bool sawDifferentH=false;
    for(std::uint64_t step=1;step<=20000;++step) {
        const double right=t0+double(step)*dt,h=right-current.time; if(h!=dt) sawDifferentH=true;
        const double qp=q+h*v+.25*h*h*a,vp=v+.5*h*a,an=(1-.4*vp)/(2+.5*h*.4);
        q=qp+.25*h*h*an; v=vp+.5*h*an; a=an;
        const auto r=d.trialStepOnGrid(current,i,t0,step,dt,tight()); current=r.candidate;
        QCOMPARE(current.time,right); QVERIFY(near(current.displacements[0],q,2e-12)); QVERIFY(near(current.velocities[0],v,2e-12));
        if(step==10000) mid=current;
    }
    QVERIFY(sawDifferentH); auto split=mid;
    for(std::uint64_t step=10001;step<=20000;++step) split=restarted.trialStepOnGrid(split,i,t0,step,dt,tight()).candidate;
    QCOMPARE(snapshot(split),snapshot(current));
    QVERIFY_EXCEPTION_THROWN(d.trialStepOnGrid(mid,i,t0,10002,dt),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(d.trialStepOnGrid(mid,i,t0,0,dt),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(d.trialStepOnGrid(mid,i,t0,(std::uint64_t(1)<<53)+1,dt),std::invalid_argument);
    const auto initial=d.initialize(t0,z,z,i);
    QVERIFY_EXCEPTION_THROWN(d.trialStepOnGrid(initial,i,t0,1,1e-20),std::invalid_argument);
}

void TestStructuralTrussDampedForcing::failurePreservesAccepted() {
    auto s=axial(); Driver d(s); std::vector<double> q(6,0),v(6,0); q[3]=.05; const auto a=d.initialize(0,q,v,force()); const auto bytes=snapshot(a);
    auto invalid=force(); invalid.rightExternalForces=force().externalForces; invalid.rightExternalForces[3]=std::numeric_limits<double>::infinity();
    QVERIFY_EXCEPTION_THROWN(d.trialStep(a,invalid,.01),std::invalid_argument); QCOMPARE(snapshot(a),bytes);
    QVERIFY_EXCEPTION_THROWN(d.trialStep(a,force(),std::numeric_limits<double>::denorm_min()),std::invalid_argument); QCOMPARE(snapshot(a),bytes);
    auto wrong=s; wrong.physicalDamping[21]*=2; Driver other(wrong); QVERIFY_EXCEPTION_THROWN(other.refresh(a,force()),std::invalid_argument);
    QCOMPARE(snapshot(a),bytes); const auto retry=d.trialStep(a,force(),.01,tight()); QVERIFY(retry.residualNorm<1e-9); QCOMPARE(snapshot(a),bytes);
}

void TestStructuralTrussDampedForcing::coupledSolveIsBitwiseDeterministicAndRestoresCallerThreads() {
    // Four coupled free physical DOFs exercise the real FE/GMRES reductions;
    // a scalar oscillator cannot reveal worker-order dot-product variation.
    System s; s.referenceNodes={{{0,0,0}},{{.25,0,0}},{{.5,0,0}},{{.75,0,0}},{{1,0,0}}};
    s.nodalLumpedMass={.125,.25,.25,.25,.125};s.prescribedMask.assign(15,true);
    for(size_t node=1;node<5;++node)s.prescribedMask[3*node]=false;
    s.fixedPrescribedDisplacements.assign(15,0);s.physicalDamping.assign(225,0);
    for(size_t node=0;node<5;++node)s.physicalDamping[(3*node)*15+3*node]=.02*s.nodalLumpedMass[node];
    for(size_t member=0;member<4;++member){
        s.members.push_back({member+1,member,member+1,member==0?RStructuralTrussLaw{RStructuralTrussPiezoLaw{.25,4000,2,.04}}:
            RStructuralTrussLaw{RStructuralTrussPassiveLaw{.25,1000}}});
        const size_t a=3*member,b=3*(member+1);const double c=.04;
        s.physicalDamping[a*15+a]+=c;s.physicalDamping[b*15+b]+=c;
        s.physicalDamping[a*15+b]-=c;s.physicalDamping[b*15+a]-=c;
    }
    const auto command=[](uint step){Inputs i;i.memberInputs={RStructuralTrussVoltage{2e-5*(step+1)},RStructuralTrussNoInput{},RStructuralTrussNoInput{},RStructuralTrussNoInput{}};
        i.externalForces.assign(15,0);for(size_t node=1;node<5;++node)i.externalForces[3*node]=1e-6*(node+1)+2e-5*double(step-1)*.00047;
        i.rightExternalForces=i.externalForces;for(size_t node=1;node<5;++node)i.rightExternalForces[3*node]+=2e-5*.00047;return i;};
    auto initialInputs=command(1);initialInputs.rightExternalForces.clear();std::get<RStructuralTrussVoltage>(initialInputs.memberInputs[0]).value=0;
    std::vector<double> q(15,0),v(15,0);for(size_t node=1;node<5;++node){q[3*node]=node*1e-8;v[3*node]=(node%2?-1:1)*double(node)*1e-5;}
    const int originalThreads=omp_get_max_threads(),originalDynamic=omp_get_dynamic();
    QByteArray initialBytes;std::vector<QByteArray> reference;
    {
        RParallelExecutionConfig serial;serial.setMode(RParallelExecutionMode::SerialDeterministic);const RParallelExecutionContext baseline(serial);
        Driver driver(s);auto state=driver.initialize(0,q,v,initialInputs);initialBytes=snapshot(state);
        for(uint step=1;step<=6;++step){const auto result=driver.trialStepOnGrid(state,command(step),0,step,.00047,tight());
            QVERIFY(result.iterations>0);state=result.candidate;reference.push_back(snapshot(state));}
    }
    for(uint requested:{1u,2u,4u}){
        RParallelExecutionConfig policy;policy.setMode(RParallelExecutionMode::FixedThreads);policy.setRequestedThreads(requested);
        const RParallelExecutionContext caller(policy);const int callerThreads=omp_get_max_threads(),callerDynamic=omp_get_dynamic();
        qInfo()<<"Axial deterministic regression caller workers"<<caller.getReport().actualThreads;
        if(requested>=2 && omp_get_num_procs()>=2 && openMPThreadLimit()>=2)QVERIFY(caller.getReport().actualThreads>=2);
        for(uint retry=0;retry<16;++retry){
            Driver fresh(s);auto state=fresh.initialize(0,q,v,initialInputs);QCOMPARE(snapshot(state),initialBytes);
            for(uint step=1;step<=6;++step){const auto before=snapshot(state);
                const auto result=fresh.trialStepOnGrid(state,command(step),0,step,.00047,tight());
                QCOMPARE(snapshot(state),before);QCOMPARE(snapshot(result.candidate),reference.at(step-1));
                QCOMPARE(omp_get_max_threads(),callerThreads);QCOMPARE(omp_get_dynamic(),callerDynamic);state=result.candidate;
            }
        }
        // The singular transverse tangent reaches the legacy solver and then
        // fails its fresh certificate. Failure must also restore caller ICVs.
        auto singular=axial(0);singular.members[0].law=RStructuralTrussPassiveLaw{2,100};
        singular.nodalLumpedMass[1]=3.125;singular.prescribedMask[4]=false;
        Driver failed(singular);auto input=force(-50);input.externalForces[4]=1;
        std::vector<double> zero(6,0),velocity(6,0);velocity[4]=-.04;
        const auto accepted=failed.initialize(0,zero,velocity,input);const auto saved=snapshot(accepted);
        QVERIFY_EXCEPTION_THROWN(failed.trialStep(accepted,input,.5,tight()),std::runtime_error);
        QCOMPARE(snapshot(accepted),saved);QCOMPARE(omp_get_max_threads(),callerThreads);QCOMPARE(omp_get_dynamic(),callerDynamic);
    }
    QCOMPARE(omp_get_max_threads(),originalThreads);QCOMPARE(omp_get_dynamic(),originalDynamic);
}

QTEST_APPLESS_MAIN(TestStructuralTrussDampedForcing)
#include "tst_structural_truss_damped_forcing.moc"

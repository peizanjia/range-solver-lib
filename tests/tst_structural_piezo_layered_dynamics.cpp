#include <QtTest>
#include <rml_model.h>
#include <rbl_error.h>
#include "rstructural_piezo_layered_dynamics.h"
#include "rparallel_execution_context.h"
#include <omp.h>
#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>

namespace {
using V=std::vector<double>;
using Setup=RStructuralPiezoLayeredSetup;
using Adapter=RStructuralPiezoLayeredModelAdapter;
using Driver=RStructuralPiezoLayeredDynamics;
using State=RStructuralPiezoLayeredDynamicsState;
using Dof=RStructuralPiezoPhysicalDof;
bool rejects(const std::function<void()> &f) {try{f();}catch(...){return true;}return false;}
bool near(double a,double b,double tolerance=2e-8) {return std::fabs(a-b)<=tolerance*std::max({1.0,std::fabs(a),std::fabs(b)});}
void vectorNear(const V &a,const V &b,double tolerance=2e-8) {
    QCOMPARE(a.size(),b.size()); for(size_t i=0;i<a.size();++i) QVERIFY2(near(a[i],b[i],tolerance),qPrintable(QString("component %1: %2 vs %3").arg(i).arg(a[i],0,'g',17).arg(b[i],0,'g',17)));
}
V mv(const V &a,const V &x,size_t rows) {V y(rows); for(size_t i=0;i<rows;++i)for(size_t j=0;j<x.size();++j)y[i]+=a[i*x.size()+j]*x[j];return y;}
V solve(V a,V b) {
    const size_t n=b.size();
    for(size_t k=0;k<n;++k) {
        size_t pivot=k;for(size_t i=k+1;i<n;++i)if(std::fabs(a[i*n+k])>std::fabs(a[pivot*n+k]))pivot=i;
        if(a[pivot*n+k]==0)throw std::runtime_error("Singular independent test reference");
        for(size_t j=k;j<n;++j)std::swap(a[k*n+j],a[pivot*n+j]);std::swap(b[k],b[pivot]);
        for(size_t i=k+1;i<n;++i){const double factor=a[i*n+k]/a[k*n+k];for(size_t j=k+1;j<n;++j)a[i*n+j]-=factor*a[k*n+j];b[i]-=factor*b[k];}
    }
    V x(n);for(size_t ii=n;ii>0;--ii){const size_t i=ii-1;x[i]=b[i];for(size_t j=i+1;j<n;++j)x[i]-=a[i*n+j]*x[j];x[i]/=a[i*n+i];}return x;
}
// MSVC long double is binary64. Recover product/sum roundoff explicitly rather
// than assuming an extended-precision type is available in this independent oracle.
struct AccurateSum {
    double high=0,low=0;
    void add(double value) {
        const double sum=high+value,tail=sum-high;
        low+=(high-(sum-tail))+(value-tail);high=sum;
    }
    void product(double a,double b) {const double value=a*b;add(value);add(std::fma(a,b,-value));}
    double value() const {return high+low;}
};
double gamma(size_t operations) {
    const double nu=double(operations)*std::numeric_limits<double>::epsilon()/2;
    return nu/(1-nu);
}
V accurateMassSolve(const V &mass,const std::vector<AccurateSum> &rhs) {
    const size_t n=rhs.size();V rounded(n);for(size_t i=0;i<n;++i)rounded[i]=rhs[i].value();
    auto answer=solve(mass,rounded);
    // Independent pivoted elimination followed by compensated ORIGINAL-matrix
    // residual corrections; this does not call the production linear solver.
    for(uint iteration=0;iteration<3;++iteration){
        V residual(n);for(size_t i=0;i<n;++i){auto sum=rhs[i];for(size_t j=0;j<n;++j)sum.product(-mass[i*n+j],answer[j]);residual[i]=sum.value();}
        const auto correction=solve(mass,residual);bool changed=false;
        for(size_t i=0;i<n;++i){const double next=answer[i]+correction[i];changed=changed||next!=answer[i];answer[i]=next;}
        if(!changed)break;
    }
    return answer;
}
struct MassOracle {V acceleration,roundoffBound;};
MassOracle newmarkMassOracle(const Adapter &adapter,const State &state,double h) {
    const auto &d=adapter.data();const auto &s=state.integrator.state;const size_t n=d.freeDofs;
    const double qFactor=.25*h*h,vFactor=.5*h,u=std::numeric_limits<double>::epsilon()/2;
    V force(n),predictorQ(n),predictorV(n),qRound(n),vRound(n);
    std::vector<AccurateSum> acceptedRhs(n);
    for(size_t i=0;i<n;++i){
        AccurateSum sum;sum.add(d.backgroundOffset[i]);
        for(size_t j=0;j<state.heldActuatorVoltages.size();++j)
            sum.product(d.input[i*state.heldActuatorVoltages.size()+j],state.heldActuatorVoltages[j]-d.actuatorVoltages0[j]);
        force[i]=sum.value();acceptedRhs[i]=sum;
        // Bound the unseen predictor from the rounded final update. Its
        // accumulation into q/v contributes a force error after multiplication
        // by K/C; this is why an accepted-q/v mass solve differs from Newmark a.
        predictorQ[i]=(std::abs(s.q[i])+std::abs(qFactor*s.a[i]))/(1-u);
        predictorV[i]=(std::abs(s.v[i])+std::abs(vFactor*s.a[i]))/(1-u);
        qRound[i]=gamma(2)*(predictorQ[i]+std::abs(qFactor*s.a[i]));
        vRound[i]=gamma(2)*(predictorV[i]+std::abs(vFactor*s.a[i]));
        for(size_t j=0;j<n;++j){acceptedRhs[i].product(-d.stiffness[i*n+j],s.q[j]);acceptedRhs[i].product(-d.damping[i*n+j],s.v[j]);}
    }
    V rowBudget(n),inverse(n*n);
    for(size_t i=0;i<n;++i){
        double predictorTerms=std::abs(force[i]),effectiveTerms=0,updateTerms=0;
        for(size_t j=0;j<n;++j){
            predictorTerms+=std::abs(d.stiffness[i*n+j])*predictorQ[j]+std::abs(d.damping[i*n+j])*predictorV[j];
            const double effective=std::abs(d.mass[i*n+j])+qFactor*std::abs(d.stiffness[i*n+j])+vFactor*std::abs(d.damping[i*n+j]);
            effectiveTerms+=effective*std::abs(s.a[j]);
            updateTerms+=std::abs(d.stiffness[i*n+j])*qRound[j]+std::abs(d.damping[i*n+j])*vRound[j];
        }
        // gamma(2n+4): separate products/sums and two RHS subtractions.
        // gamma(4): scaled effective-matrix entries and their assembly.
        // 1e-12 is the unchanged production ORIGINAL effective-row certificate.
        rowBudget[i]=gamma(2*n+4)*predictorTerms+gamma(4)*effectiveTerms+updateTerms
            +(1e-12+gamma(2*n+2))*(predictorTerms+effectiveTerms);
        std::vector<AccurateSum> unit(n);unit[i].add(1);const auto column=accurateMassSolve(d.mass,unit);
        for(size_t j=0;j<n;++j)inverse[j*n+i]=column[j];
    }
    const auto oracle=accurateMassSolve(d.mass,acceptedRhs);
    // A computed inverse is not assumed exact. The Neumann bound for
    // E=I-inverse*M supplies the remaining forward-error amplification.
    double defect=0;V propagated(n);
    for(size_t i=0;i<n;++i){
        double rowDefect=0;
        for(size_t j=0;j<n;++j){AccurateSum sum;sum.add(i==j?1.:0.);double scale=i==j?1.:0.;
            for(size_t k=0;k<n;++k){sum.product(-inverse[i*n+k],d.mass[k*n+j]);scale+=std::abs(inverse[i*n+k]*d.mass[k*n+j]);}
            rowDefect+=std::abs(sum.value())+gamma(2*n+1)*scale;
        }
        defect=std::max(defect,rowDefect);
        AccurateSum residual=acceptedRhs[i];for(size_t j=0;j<n;++j)residual.product(-d.mass[i*n+j],oracle[j]);
        rowBudget[i]+=std::abs(residual.value());
    }
    if(!(defect<1))throw std::runtime_error("Independent mass inverse has unresolved conditioning");
    for(size_t i=0;i<n;++i)for(size_t j=0;j<n;++j)propagated[i]+=std::abs(inverse[i*n+j])*rowBudget[j];
    const double maximum=*std::max_element(propagated.begin(),propagated.end());
    for(size_t i=0;i<n;++i)propagated[i]=propagated[i]/(1-gamma(2*n))+defect*maximum/(1-defect);
    return{oracle,propagated};
}
double quadratic(const V &a,const V &x) {const auto y=mv(a,x,x.size());double result=0;for(size_t i=0;i<x.size();++i)result+=x[i]*y[i];return result;}
RMatrixSolverConf solver() {RMatrixSolverConf c(RMatrixSolverConf::GMRES);c.setNInnerIterations(100);c.setNOuterIterations(100);c.setSolverCvgValue(1e-13);return c;}
RModel scene(bool plate=false) {
    RModel m;const double xy[4][2]={{0,0},{1,0},{1,1},{0,1}};
    for(uint i=0;i<(plate?4u:3u);++i){m.addNode(plate?RNode(xy[i][0],xy[i][1],0):RNode(double(i),0,0));m.setNodePhysicalId(i,101+11*i);}
    for(uint i=0;i<(plate?1u:2u);++i){RElement e(plate?R_ELEMENT_QUAD1:R_ELEMENT_BEAM1);for(uint j=0;j<(plate?4u:2u);++j)e.setNodeId(j,plate?j:i+j);m.addElement(e,true,0);m.setElementPhysicalId(i,501+13*i);}
    // This fixture declares its complete loading through the layered setup.
    // removeAllEnvironmentConditions() would reinstate these legacy defaults.
    for(uint i=0;i<(plate?m.getNSurfaces():m.getNLines());++i){
        RElementGroup &group=plate?static_cast<RElementGroup &>(m.getSurface(i)):static_cast<RElementGroup &>(m.getLine(i));
        group.removeEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
        group.removeEnvironmentCondition(R_ENVIRONMENT_CONDITION_TEMPERATURE);
    }
    m.setProblemTaskTree(RProblemTaskItem::generateDefaultTree(R_PROBLEM_STRESS));m.getProblemSetup().setMassMatrixType(RMassMatrixType::Consistent);
    auto &d=m.getProblemSetup().getStructuralDynamicsSetup();d.setEnabled(true);d.setTimeStep(.01);d.setNSteps(10);m.getTimeSolver().setInputStartTime(.01);return m;
}
Setup config(bool plate=false) {
    Setup s;s.enabled=true;if(plate)s.family=RStructuralPiezoLayeredModelFamily::FlatPlate;
    auto native=RStructuralPiezoContinuumMaterial::isotropic(30,.2,7);
    RStructuralPiezoLayeredMaterial host;host.materialId="host";host.elasticityOrCompliance=native.elasticity;host.density=7;host.sourceId="independent-example";host.sourceSha256=QString(64,'a');
    auto active=host;active.materialId="piezo";active.piezoelectric=true;active.coupling[12]=2;active.coupling[13]=2;active.permittivity={{4,0,0,0,4,0,0,0,4}};s.materials={host,active};
    RStructuralPiezoLayeredPort drive;drive.readout.portId="drive";drive.readout.role=RStructuralPiezoPortRole::ActuatorSensor;drive.readout.readout=RStructuralPiezoReadoutMode::ActuatorCharge;drive.initialValue=.12;
    RStructuralPiezoLayeredPort sense;sense.readout.portId="sense";sense.readout.role=RStructuralPiezoPortRole::Sensor;sense.readout.readout=RStructuralPiezoReadoutMode::OpenCircuitVoltage;sense.initialValue=.03;
    RStructuralPiezoLayeredPort shorted;shorted.readout.portId="short";shorted.readout.role=RStructuralPiezoPortRole::Sensor;shorted.readout.readout=RStructuralPiezoReadoutMode::ShortCircuitCharge;
    RStructuralPiezoLayeredPort amp;amp.readout.portId="amp";amp.readout.role=RStructuralPiezoPortRole::Sensor;amp.readout.readout=RStructuralPiezoReadoutMode::IdealChargeAmplifier;amp.readout.Cf=2;amp.readout.Qref=.2;amp.readout.amplifierPolarity=1;amp.readout.VoutBias=.7;
    s.ports={drive,sense,shorted,amp};
    for(uint i=0;i<(plate?1u:2u);++i){
        RStructuralPiezoLayeredModelLayer h;h.materialId="host";h.lowerZ=.1;h.upperZ=.3;
        RStructuralPiezoLayeredModelLayer p;p.materialId="piezo";p.lowerZ=.3;p.upperZ=.4;p.portId="drive";
        auto open=p;open.lowerZ=.4;open.upperZ=.45;open.portId="sense";
        auto sh=p;sh.lowerZ=.45;sh.upperZ=.5;sh.portId="short";
        auto am=p;am.lowerZ=.5;am.upperZ=.55;am.portId="amp";
        s.members.push_back({501+13*i,{h,p,open,sh,am}});
    }
    for(uint node:plate?std::vector<uint>{101,134}:std::vector<uint>{101})
        for(Dof d:plate?std::vector<Dof>{Dof::UX,Dof::UY,Dof::UZ,Dof::RX,Dof::RY}:std::vector<Dof>{Dof::UX,Dof::UZ,Dof::RY})s.supports.push_back({node,d,0});
    s.backgroundForces={{123,{{.2,0,-.1,0,0,0}}}};
    s.damping.kind=RStructuralPiezoDampingKind::Rayleigh;s.damping.rayleighAlpha=.06;s.damping.rayleighBeta=.005;
    return s;
}
void initialPerturbation(Setup &s,const Adapter &base) {
    const auto q=base.mesh().reconstructPhysical(base.data().q0);
    for(size_t i=0;i<base.mesh().nodeIds.size();++i){RStructuralPiezoLayeredNodeField field;field.physicalNodeId=base.mesh().nodeIds[i];std::copy(q.begin()+6*i,q.begin()+6*i+6,field.values.begin());if(field.physicalNodeId==123){field.values[0]+=.015;field.values[2]+=.02;}s.initialDisplacements.push_back(field);}
    s.initialVelocities={{123,{{.04,0,.08,0,0,0}}}};
}
V perturb(const Adapter &a,double t) {V f(6*a.mesh().nodeIds.size());const auto it=std::find(a.mesh().nodeIds.begin(),a.mesh().nodeIds.end(),123);const auto i=size_t(it-a.mesh().nodeIds.begin());f[6*i]=.2*t;f[6*i+2]=.1*std::sin(2*t);return f;}
V rhs(const Adapter &a,const V &physical,const V &va) {V result=a.data().backgroundOffset;const auto f=a.reduceForce(physical);for(size_t i=0;i<result.size();++i){result[i]+=f[i];for(size_t j=0;j<va.size();++j)result[i]+=a.data().input[i*va.size()+j]*(va[j]-a.data().actuatorVoltages0[j]);}return result;}
struct Reference {V q,v,a;};
Reference advance(const Adapter &adapter,const Reference &left,const V &rightForce,double h) {
    const auto &d=adapter.data();const size_t n=d.freeDofs;V q=left.q,v=left.v,operatorMatrix=d.mass;
    for(size_t i=0;i<n;++i){q[i]+=h*left.v[i]+.25*h*h*left.a[i];v[i]+=.5*h*left.a[i];}
    V f=rightForce,cv=mv(d.damping,v,n),kq=mv(d.stiffness,q,n);
    for(size_t i=0;i<n;++i)f[i]-=cv[i]+kq[i];
    for(size_t i=0;i<n*n;++i)operatorMatrix[i]+=.5*h*d.damping[i]+.25*h*h*d.stiffness[i];
    const auto a=solve(operatorMatrix,f);for(size_t i=0;i<n;++i){q[i]+=.25*h*h*a[i];v[i]+=.5*h*a[i];}return{q,v,a};
}
}

#define REPORT_LAYERED_EXCEPTIONS \
    catch(const RError &error){QFAIL(qPrintable(error.getMessage()));} \
    catch(const std::exception &error){QFAIL(error.what());}

class TestStructuralPiezoLayeredDynamics:public QObject {
    Q_OBJECT
private slots:
    void actualBeamPlateIndependentNewmark_data();
    void actualBeamPlateIndependentNewmark();
    void thinPlateVoltageJumpKeepsEveryFreeRowCertified();
    void rampEndpointForceAndVoltageJump();
    void scalarHarmonicTimeRefinement();
    void viscousContinuumTimeRefinement();
    void fullEnergyDampingReactionsAndReadouts();
    void splitRestoreBitEquivalence();
    void tamperedStateAndIdentityRejection();
    void forcingAndNonfiniteFailureRetry();
    void canonicalGridActualInterval_data();
    void canonicalGridActualInterval();
    void serialPolicyRestoresParallelHost();
};
void TestStructuralPiezoLayeredDynamics::actualBeamPlateIndependentNewmark_data(){QTest::addColumn<bool>("plate");QTest::newRow("beam")<<false;QTest::newRow("plate")<<true;}
void TestStructuralPiezoLayeredDynamics::actualBeamPlateIndependentNewmark() try {
    QFETCH(bool,plate);auto model=scene(plate);auto setup=config(plate);Adapter base(model,setup);initialPerturbation(setup,base);Adapter adapter(model,setup);
    Driver driver(adapter,solver(),.25,[&](double t){return perturb(adapter,t);},"ramp+sin/native-example/v1");
    auto state=driver.initialize();const auto &d=adapter.data();
    V z0(d.freeDofs),q0=d.q0;for(size_t i=0;i<q0.size();++i)q0[i]-=d.p[i];
    for(size_t i=0;i<d.fullDofs;++i)for(size_t j=0;j<d.freeDofs;++j)z0[j]+=d.T[i*d.freeDofs+j]*q0[i];
    Reference ref;ref.q=d.initialZ;for(size_t i=0;i<ref.q.size();++i)ref.q[i]-=z0[i];ref.v=d.initialVelocity;
    auto f=rhs(adapter,perturb(adapter,.25),d.actuatorVoltages0),cv=mv(d.damping,ref.v,d.freeDofs),kq=mv(d.stiffness,ref.q,d.freeDofs);
    for(size_t i=0;i<f.size();++i)f[i]-=cv[i]+kq[i];ref.a=solve(d.mass,f);
    vectorNear(state.integrator.state.q,ref.q);vectorNear(state.integrator.state.a,ref.a);
    bool offDiagonal=false;for(size_t i=0;i<d.freeDofs;++i)for(size_t j=0;j<d.freeDofs;++j)if(i!=j && d.mass[i*d.freeDofs+j]!=0)offDiagonal=true;QVERIFY(offDiagonal);
    for(uint step=0;step<24;++step){const double target=state.integrator.state.time+.002;ref=advance(adapter,ref,rhs(adapter,perturb(adapter,target),state.heldActuatorVoltages),.002);state=driver.trial(state,state.heldActuatorVoltages,.002);
        vectorNear(state.integrator.state.q,ref.q);vectorNear(state.integrator.state.v,ref.v);vectorNear(state.integrator.state.a,ref.a);
        V q=mv(d.T,ref.q,d.fullDofs);for(size_t i=0;i<q.size();++i)q[i]+=d.q0[i];vectorNear(state.physicalDisplacements,adapter.mesh().reconstructPhysical(q));
    }
} REPORT_LAYERED_EXCEPTIONS

void TestStructuralPiezoLayeredDynamics::thinPlateVoltageJumpKeepsEveryFreeRowCertified() try {
    // Same physical scales and asymmetric driven/open-circuit layers as the
    // workbench SI plate example. Its membrane/bending rows differ greatly in
    // scale; a norm-only solve can miss an almost-unloaded individual row.
    auto model=scene(true);auto setup=config(true);setup.backgroundForces.clear();setup.damping={};
    auto material=RStructuralPiezoContinuumMaterial::isotropic(7e10,.3,2700);
    auto &host=setup.materials[0];host.elasticityOrCompliance=material.elasticity;host.density=2700;
    auto &active=setup.materials[1];active.elasticityOrCompliance=material.elasticity;active.density=7500;
    active.coupling.fill(0);active.coupling[12]=.1;active.coupling[13]=.1;active.permittivity={{1e-8,0,0,0,1e-8,0,0,0,1e-8}};
    setup.ports.resize(2);for(auto &port:setup.ports)port.initialValue=0;
    auto &layers=setup.members[0].layers;layers.resize(3);
    layers[0].lowerZ=-.005;layers[0].upperZ=.005;
    layers[1].lowerZ=.005;layers[1].upperZ=.0055;
    layers[2].lowerZ=-.0055;layers[2].upperZ=-.005;
    Adapter adapter(model,setup);RMatrixSolverConf policy(RMatrixSolverConf::GMRES);
    Driver driver(adapter,policy,0,[](double){return V(24);},"unforced-controlled-plate-regression/v1");
    auto state=driver.initialize();auto voltage=state.heldActuatorVoltages;voltage[0]=.3;
    const auto initial=state;Reference reference{state.integrator.state.q,state.integrator.state.v,state.integrator.state.a};
    auto refreshed=rhs(adapter,V(24),voltage);const auto cv=mv(adapter.data().damping,reference.v,adapter.data().freeDofs),kq=mv(adapter.data().stiffness,reference.q,adapter.data().freeDofs);
    for(size_t i=0;i<refreshed.size();++i)refreshed[i]-=cv[i]+kq[i];reference.a=solve(adapter.data().mass,refreshed);
    for(std::uint64_t step=1;step<=3;++step){
        reference=advance(adapter,reference,rhs(adapter,V(24),voltage),.001);
        state=driver.trialOnGrid(state,voltage,step,.001);
        vectorNear(state.integrator.state.q,reference.q,1e-13);
        vectorNear(state.integrator.state.v,reference.v,1e-13);
        // Newmark a solves its effective operator BEFORE q/v updates round.
        // Compare an independent compensated ORIGINAL-mass oracle at accepted
        // q/v using the derived operation/conditioning contribution PLUS the
        // existing independent comparison tolerance, 1e-11. The contribution
        // accounts for the thin plate's K/M scales; it is not an FE certificate.
        const double leftTime=step==1?0.:double(step-1)*.001;
        const auto oracle=newmarkMassOracle(adapter,state,double(step)*.001-leftTime);
        for(size_t i=0;i<oracle.acceleration.size();++i){
            const double difference=std::abs(state.integrator.state.a[i]-oracle.acceleration[i]);
            const double comparisonTolerance=1e-11*std::max({1.,std::abs(oracle.acceleration[i]),std::abs(state.integrator.state.a[i])});
            QVERIFY2(difference<=oracle.roundoffBound[i]+comparisonTolerance,qPrintable(QString("step %1 row %2: independent mass discrepancy=%3 derived roundoff/conditioning bound=%4 existing comparison tolerance=%5")
                .arg(qulonglong(step)).arg(qulonglong(i)).arg(difference,0,'g',17).arg(oracle.roundoffBound[i],0,'g',17).arg(comparisonTolerance,0,'g',17)));
        }
        qInfo()<<"thin-plate independent mass oracle"<<"step"<<step<<"row3Discrepancy"
               <<std::abs(state.integrator.state.a[3]-oracle.acceleration[3])<<"row3DerivedRoundoffBound"<<oracle.roundoffBound[3]
               <<"row3ExistingComparisonTolerance"<<1e-11*std::max({1.,std::abs(oracle.acceleration[3]),std::abs(state.integrator.state.a[3])});
        QCOMPARE(driver.restore(state).integrityDigest,state.integrityDigest);
    }
    QVERIFY(state.electricalWork!=0);QVERIFY(state.freeForceResidual>=0);
    QCOMPARE(driver.restore(initial).integrityDigest,initial.integrityDigest);
} REPORT_LAYERED_EXCEPTIONS

void TestStructuralPiezoLayeredDynamics::viscousContinuumTimeRefinement() try {
    auto model=scene();auto setup=config();setup.backgroundForces={{123,{{.2,0,0,0,0,0}}}};
    for(auto dof:{Dof::UX,Dof::UZ,Dof::RY})setup.supports.push_back({112,dof,0});
    setup.supports.push_back({123,Dof::UZ,0});setup.supports.push_back({123,Dof::RY,0});
    setup.initialVelocities={{123,{{.3,0,0,0,0,0}}}};
    setup.damping.rayleighAlpha=.2;setup.damping.rayleighBeta=.005;Adapter adapter(model,setup);
    QCOMPARE(adapter.data().freeDofs,1u);
    const double m=adapter.data().mass[0],k=adapter.data().stiffness[0],c=adapter.data().damping[0],v0=.3,t=.5,alpha=c/(2*m),omega=std::sqrt(k/m-alpha*alpha);
    const double offset=adapter.data().backgroundOffset[0],shift=offset/k,e=std::exp(-alpha*t);
    const double q=shift*(1-e*(std::cos(omega*t)+alpha/omega*std::sin(omega*t)))+v0/omega*e*std::sin(omega*t);
    const double v=shift*e*(k/m)/omega*std::sin(omega*t)+v0*e*(std::cos(omega*t)-alpha/omega*std::sin(omega*t));
    const double exact=.5*m*v0*v0+offset*q-.5*m*v*v-.5*k*q*q;
    double errors[3]{};
    for(uint refinement=0;refinement<3;++refinement){const uint steps=20u<<refinement;const double h=t/steps;
        Driver driver(adapter,solver(),0,[](double){return V(18);},"scalar-free-viscous/v1");auto state=driver.initialize();
        for(uint step=0;step<steps;++step)state=driver.trial(state,state.heldActuatorVoltages,h);
        errors[refinement]=std::fabs(state.viscousDissipation-exact);
        // Both quadratures approximate continuum dissipation; midpoint is the
        // exact work paired with this Newmark trajectory, trapezoid is a rate diagnostic.
        QVERIFY(state.integrator.state.dampingDissipation>=state.viscousDissipation);
    }
    QVERIFY(errors[1]<errors[0]/3);QVERIFY(errors[2]<errors[1]/3);
} REPORT_LAYERED_EXCEPTIONS

void TestStructuralPiezoLayeredDynamics::rampEndpointForceAndVoltageJump() try {
    auto model=scene();auto setup=config();Adapter adapter(model,setup);
    auto ramp=[&](double t){V f(18);f[12]=3*t;return f;};Driver driver(adapter,solver(),0,ramp,"linear-ramp/v1");
    const auto left=driver.initialize();const double h=.05;auto va=left.heldActuatorVoltages;va[0]+=.3;
    auto force=rhs(adapter,ramp(0),va),kq=mv(adapter.data().stiffness,left.integrator.state.q,adapter.data().freeDofs),cv=mv(adapter.data().damping,left.integrator.state.v,adapter.data().freeDofs);
    for(size_t i=0;i<force.size();++i)force[i]-=kq[i]+cv[i];Reference ref{left.integrator.state.q,left.integrator.state.v,solve(adapter.data().mass,force)};
    const auto expected=advance(adapter,ref,rhs(adapter,ramp(h),va),h);
    const auto actual=driver.trial(left,va,h);
    vectorNear(actual.integrator.state.q,expected.q);vectorNear(actual.integrator.state.v,expected.v);
    auto endHold=rhs(adapter,ramp(h),va);for(size_t i=0;i<endHold.size();++i)endHold[i]-=kq[i]+cv[i];Reference wrong{ref.q,ref.v,solve(adapter.data().mass,endHold)};
    const auto incorrect=advance(adapter,wrong,rhs(adapter,ramp(h),va),h);double difference=0;for(size_t i=0;i<incorrect.q.size();++i)difference+=std::fabs(incorrect.q[i]-expected.q[i]);QVERIFY(difference>1e-7);
    double refreshDifference=0;for(size_t i=0;i<ref.a.size();++i)refreshDifference+=std::fabs(ref.a[i]-left.integrator.state.a[i]);QVERIFY(refreshDifference>1e-5);
    auto initialQ=adapter.mesh().projectPhysical(left.physicalDisplacements,RStructuralPiezoPhysicalDomain::Displacement);const auto jump=adapter.electrical().evaluate(initialQ,va);
    double work=0;for(size_t i=0;i<left.electrical.ports.size();++i)if(left.electrical.ports[i].role==RStructuralPiezoElectricalRole::VoltageActuator){work+=.5*(left.electrical.ports[i].voltage+jump.ports[i].voltage)*(jump.ports[i].charge-left.electrical.ports[i].charge);work+=jump.ports[i].voltage*(actual.electrical.ports[i].charge-jump.ports[i].charge);}
    QVERIFY(near(actual.electricalWork,work,1e-10));
    QVERIFY(std::fabs(actual.energyBalanceResidual)<1e-9*std::max(1.0,std::fabs(actual.initialEnergy)));
    QCOMPARE(left.integrityDigest,driver.initialize().integrityDigest);
} REPORT_LAYERED_EXCEPTIONS

void TestStructuralPiezoLayeredDynamics::scalarHarmonicTimeRefinement() try {
    auto model=scene();auto setup=config();setup.damping={};setup.backgroundForces={{123,{{.2,0,0,0,0,0}}}};
    // Actual layered mesh with one admitted axial free coordinate. All other
    // beam coordinates are held, retaining consistent inertia in that subspace.
    for(auto dof:{Dof::UX,Dof::UZ,Dof::RY})setup.supports.push_back({112,dof,0});
    setup.supports.push_back({123,Dof::UZ,0});setup.supports.push_back({123,Dof::RY,0});
    Adapter adapter(model,setup);QCOMPARE(adapter.data().freeDofs,1u);
    const double m=adapter.data().mass[0],k=adapter.data().stiffness[0],omega=.7,natural=std::sqrt(k/m),amplitude=.25,finalTime=.5;
    const double exact=amplitude/(k-m*omega*omega)*(std::sin(omega*finalTime)-omega/natural*std::sin(natural*finalTime))+adapter.data().backgroundOffset[0]/k*(1-std::cos(natural*finalTime));
    double errors[3]{};
    for(uint refinement=0;refinement<3;++refinement){const uint steps=20u<<refinement;const double h=finalTime/steps;
        Driver driver(adapter,solver(),0,[=](double t){V f(18);f[12]=amplitude*std::sin(omega*t);return f;},"analytic-scalar-sine/v1");auto state=driver.initialize();
        for(uint step=0;step<steps;++step)state=driver.trial(state,state.heldActuatorVoltages,h);errors[refinement]=std::fabs(state.integrator.state.q[0]-exact);
    }
    QVERIFY(errors[0]>0);QVERIFY(errors[1]<errors[0]/3);QVERIFY(errors[2]<errors[1]/3);
} REPORT_LAYERED_EXCEPTIONS

void TestStructuralPiezoLayeredDynamics::fullEnergyDampingReactionsAndReadouts() try {
    auto model=scene();auto setup=config();setup.supports[0].value=.01;setup.supports[1].value=.02;Adapter base(model,setup);initialPerturbation(setup,base);Adapter adapter(model,setup);
    Driver driver(adapter,solver(),0,[](double){return V(18);},"zero-perturbation/v1");auto state=driver.initialize();const auto initialEnergy=state.initialEnergy;
    for(uint step=0;step<100;++step){const auto left=state;state=driver.trial(left,left.heldActuatorVoltages,.002);
        const auto q=adapter.mesh().projectPhysical(state.physicalDisplacements,RStructuralPiezoPhysicalDomain::Displacement),v=adapter.mesh().projectPhysical(state.physicalVelocities,RStructuralPiezoPhysicalDomain::Velocity),a=adapter.mesh().projectPhysical(state.physicalAccelerations,RStructuralPiezoPhysicalDomain::Velocity);
        V voltage;for(const auto &p:state.electrical.ports)voltage.push_back(p.voltage);
        QVERIFY(near(state.kineticEnergy,.5*quadratic(adapter.mesh().operators.mass,v),1e-10));
        QVERIFY(near(state.storedEnergy,.5*quadratic(adapter.mesh().operators.stiffness,q)+.5*quadratic(adapter.mesh().operators.capacitance,voltage),1e-10));
        const auto ma=mv(adapter.mesh().operators.mass,a,a.size()),cv=mv(adapter.data().fullDamping,v,v.size());V reaction(a.size());for(size_t i=0;i<reaction.size();++i)reaction[i]=ma[i]+cv[i]+state.electrical.internalForce[i]-adapter.data().backgroundForce[i];vectorNear(state.physicalReactions,adapter.mesh().reconstructPhysical(reaction));
        auto mean=adapter.mesh().projectPhysical(left.physicalVelocities,RStructuralPiezoPhysicalDomain::Velocity);for(size_t i=0;i<mean.size();++i)mean[i]=.5*(mean[i]+v[i]);
        QVERIFY(near(state.viscousDissipation-left.viscousDissipation,.002*quadratic(adapter.data().fullDamping,mean),1e-11));
        QVERIFY(state.viscousDissipation>=left.viscousDissipation);QVERIFY(std::fabs(state.energyBalanceResidual)<1e-8*std::max(1.,initialEnergy));
        QCOMPARE(state.physicalDisplacements[0],.01);QCOMPARE(state.physicalDisplacements[2],.02);
        for(size_t i=0;i<state.electrical.ports.size();++i){const auto &p=state.electrical.ports[i];const auto &readout=state.readouts[i];
            if(p.portId=="sense"){QCOMPARE(p.charge,.03);QCOMPARE(readout.unit,QString("V"));}
            if(p.portId=="short"){QCOMPARE(p.voltage,0.0);QCOMPARE(readout.value,p.charge);}
            if(p.portId=="amp"){QCOMPARE(p.voltage,0.0);QVERIFY(near(readout.value,.7+(p.charge-.2)/2,1e-12));}
        }
    }
    QVERIFY(state.viscousDissipation>0);QVERIFY(state.integrator.state.dampingDissipation>=state.viscousDissipation);
} REPORT_LAYERED_EXCEPTIONS

void TestStructuralPiezoLayeredDynamics::splitRestoreBitEquivalence() try {
    auto model=scene(true);auto setup=config(true);Adapter adapter(model,setup);auto callback=[&](double t){return perturb(adapter,t);};Driver driver(adapter,solver(),0,callback,"ramp+sin/v1");
    auto continuous=driver.initialize(),split=continuous;
    for(uint i=0;i<18;++i){auto va=continuous.heldActuatorVoltages;va[0]=i<5?.12:.19;continuous=driver.trial(continuous,va,.003);}
    for(uint i=0;i<7;++i){auto va=split.heldActuatorVoltages;va[0]=i<5?.12:.19;split=driver.trial(split,va,.003);}
    const auto checkpoint=driver.checkpoint(split);Driver restarted(adapter,solver(),0,callback,"ramp+sin/v1");auto restored=restarted.restore(checkpoint);
    for(uint i=7;i<18;++i)restored=restarted.trial(restored,restored.heldActuatorVoltages,.003);
    QCOMPARE(restored.integrityDigest,continuous.integrityDigest);QCOMPARE(restored.physicalDisplacements,continuous.physicalDisplacements);QCOMPARE(restored.physicalAccelerations,continuous.physicalAccelerations);QCOMPARE(restored.electricalWork,continuous.electricalWork);
    QVERIFY(std::fabs(restored.energyBalanceResidual)<1e-9*std::max(1.0,std::fabs(restored.initialEnergy)));
} REPORT_LAYERED_EXCEPTIONS

void TestStructuralPiezoLayeredDynamics::tamperedStateAndIdentityRejection() try {
    auto model=scene();auto setup=config();Adapter adapter(model,setup);Driver driver(adapter,solver(),0,[](double){return V(18);},"zero/v1");
    const auto initial=driver.initialize();const auto state=driver.trial(initial,initial.heldActuatorVoltages,.002);
    const std::vector<std::function<void(State &)>> changes={
        [](State &s){s.sourceFingerprint="stale";},[](State &s){s.forcingIdentity="other";},[](State &s){s.integrator.state.q[0]+=.001;},
        [](State &s){s.integrator.state.a[0]+=.001;},[](State &s){s.integrator.fixedStep*=2;},[](State &s){s.integrator.initialMechanicalEnergy+=1;},
        [](State &s){s.physicalReactions[0]+=.1;},[](State &s){s.electrical.ports[0].charge+=.01;},[](State &s){s.readouts[0].value+=1;},
        [](State &s){s.heldActuatorVoltages[0]+=.01;},[](State &s){s.initialEnergy+=1;},[](State &s){s.electricalWork+=1;},
        [](State &s){s.perturbationForces[0]=1;},[](State &s){s.integrator.state.v.pop_back();},[](State &s){s.viscousDissipation=std::numeric_limits<double>::quiet_NaN();}};
    for(const auto &change:changes){auto bad=state;change(bad);QVERIFY(rejects([&]{driver.restore(bad);}));QVERIFY(rejects([&]{driver.trial(bad,state.heldActuatorVoltages,.002);}));}
    Driver other(adapter,solver(),0,[](double){return V(18);},"zero/v2");QVERIFY(rejects([&]{other.restore(state);}));
    auto changed=setup;changed.damping.rayleighAlpha+=.1;Adapter changedAdapter(model,changed);Driver changedDriver(changedAdapter,solver(),0,[](double){return V(18);},"zero/v1");QVERIFY(rejects([&]{changedDriver.restore(state);}));
    auto shortCommand=state.heldActuatorVoltages;shortCommand.clear();QVERIFY(rejects([&]{driver.trial(state,shortCommand,.002);}));
    QVERIFY(rejects([&]{driver.trial(state,state.heldActuatorVoltages,.003);}));QCOMPARE(driver.restore(state).integrityDigest,state.integrityDigest);
} REPORT_LAYERED_EXCEPTIONS

void TestStructuralPiezoLayeredDynamics::forcingAndNonfiniteFailureRetry() try {
    auto model=scene();auto setup=config();Adapter adapter(model,setup);bool throwAtRight=false,nonfiniteAtRight=false;
    Driver driver(adapter,solver(),0,[&](double t){if(t>0 && throwAtRight)throw std::runtime_error("Injected forcing failure");V f(18);if(t>0 && nonfiniteAtRight)f[12]=std::numeric_limits<double>::infinity();return f;},"fault-injection-zero/v1");
    const auto accepted=driver.initialize();throwAtRight=true;QVERIFY(rejects([&]{driver.trial(accepted,accepted.heldActuatorVoltages,.002);}));throwAtRight=false;
    nonfiniteAtRight=true;QVERIFY(rejects([&]{driver.trial(accepted,accepted.heldActuatorVoltages,.002);}));nonfiniteAtRight=false;
    const auto retry=driver.trial(accepted,accepted.heldActuatorVoltages,.002),again=driver.trial(accepted,accepted.heldActuatorVoltages,.002);QCOMPARE(retry.integrityDigest,again.integrityDigest);QCOMPARE(driver.restore(accepted).integrityDigest,accepted.integrityDigest);
    auto invalidCommand=accepted.heldActuatorVoltages;invalidCommand[0]=std::numeric_limits<double>::max();QVERIFY(rejects([&]{driver.trial(accepted,invalidCommand,.002);}));QCOMPARE(driver.trial(accepted,accepted.heldActuatorVoltages,.002).integrityDigest,retry.integrityDigest);
    Driver shape(adapter,solver(),0,[](double){return V(17);},"bad-shape/v1");QVERIFY(rejects([&]{shape.initialize();}));
    Driver unsupported(adapter,solver(),0,[](double){V f(18);f[13]=1;return f;},"unsupported-force/v1");QVERIFY(rejects([&]{unsupported.initialize();}));
    // Failure after the native FE solve, while deriving the final amplifier
    // observation, also leaves this driver's next attempt independent.
    auto amplifierSetup=config();amplifierSetup.backgroundForces.clear();
    for(auto &p:amplifierSetup.ports){p.initialValue=0;if(p.readout.portId=="amp"){p.readout.Cf=1e-320;p.readout.Qref=0;}}
    Adapter amplifierAdapter(model,amplifierSetup);bool excite=false;
    Driver amplifierDriver(amplifierAdapter,solver(),0,[&](double t){V f(18);if(excite && t>0)f[12]=.1;return f;},"final-readout-fault/v1");
    const auto beforeAmplifierFailure=amplifierDriver.initialize();excite=true;
    QVERIFY(rejects([&]{amplifierDriver.trial(beforeAmplifierFailure,beforeAmplifierFailure.heldActuatorVoltages,.002);}));excite=false;
    const auto recovered=amplifierDriver.trial(beforeAmplifierFailure,beforeAmplifierFailure.heldActuatorVoltages,.002);
    QCOMPARE(recovered.physicalDisplacements,V(18));QCOMPARE(amplifierDriver.restore(beforeAmplifierFailure).integrityDigest,beforeAmplifierFailure.integrityDigest);
} REPORT_LAYERED_EXCEPTIONS

void TestStructuralPiezoLayeredDynamics::canonicalGridActualInterval_data() {
    QTest::addColumn<double>("time0");
    QTest::newRow("decimal-origin")<<.37;
    QTest::newRow("large-origin")<<1.0e10;
}
void TestStructuralPiezoLayeredDynamics::canonicalGridActualInterval() try {
    QFETCH(double,time0);
    auto model=scene();auto setup=config();setup.backgroundForces={{123,{{.2,0,0,0,0,0}}}};
    for(auto d:{Dof::UX,Dof::UZ,Dof::RY})setup.supports.push_back({112,d,0});
    setup.supports.push_back({123,Dof::UZ,0});setup.supports.push_back({123,Dof::RY,0});
    setup.initialVelocities={{123,{{.3,0,0,0,0,0}}}};
    Adapter adapter(model,setup);QCOMPARE(adapter.data().freeDofs,1u);
    const auto callback=[time0](double t){V f(18);f[12]=.1*(t-time0);return f;};
    Driver driver(adapter,solver(),time0,callback,"canonical-endpoint-ramp/v1");
    auto state=driver.initialize();const auto initial=state;
    Reference reference{state.integrator.state.q,state.integrator.state.v,state.integrator.state.a};
    const double dt=.01;State saved;double expectedDamping=0;
    for(std::uint64_t index=1;index<=64;++index){
        const double target=time0+double(index)*dt,h=target-state.integrator.state.time;
        const auto next=advance(adapter,reference,rhs(adapter,callback(target),state.heldActuatorVoltages),h);
        const double meanVelocity=.5*(reference.v[0]+next.v[0]);
        expectedDamping+=h*adapter.data().damping[0]*meanVelocity*meanVelocity;
        reference=next;state=driver.trialOnGrid(state,state.heldActuatorVoltages,index,dt);
        QCOMPARE(state.integrator.state.time,target);QCOMPARE(state.integrator.fixedStep,dt);
        vectorNear(state.integrator.state.q,reference.q);vectorNear(state.integrator.state.v,reference.v);
        vectorNear(state.integrator.state.a,reference.a);QCOMPARE(state.perturbationForces,callback(target));
        QVERIFY(near(state.viscousDissipation,expectedDamping,1e-11));
        QVERIFY(std::fabs(state.energyBalanceResidual)<1e-9*std::max(1.0,std::fabs(state.initialEnergy)));
        if(index==23)saved=driver.checkpoint(state);
    }
    Driver restarted(adapter,solver(),time0,callback,"canonical-endpoint-ramp/v1");
    auto restored=restarted.restore(saved);
    for(std::uint64_t index=24;index<=64;++index)restored=restarted.trialOnGrid(restored,restored.heldActuatorVoltages,index,dt);
    QCOMPARE(restored.integrityDigest,state.integrityDigest);
    QVERIFY(rejects([&]{driver.trialOnGrid(initial,initial.heldActuatorVoltages,0,dt);}));
    QVERIFY(rejects([&]{driver.trialOnGrid(initial,initial.heldActuatorVoltages,2,dt);}));
    QVERIFY(rejects([&]{driver.trialOnGrid(saved,saved.heldActuatorVoltages,24,std::nextafter(dt,1.0));}));
    QVERIFY(rejects([&]{driver.trialOnGrid(saved,saved.heldActuatorVoltages,25,dt);}));
    QCOMPARE(driver.restore(saved).integrityDigest,saved.integrityDigest);
} REPORT_LAYERED_EXCEPTIONS

void TestStructuralPiezoLayeredDynamics::serialPolicyRestoresParallelHost() try {
    RParallelExecutionConfig hostPolicy;hostPolicy.setMode(RParallelExecutionMode::FixedThreads);
    hostPolicy.setRequestedThreads(4);hostPolicy.setDeterministicReduction(false);
    const RParallelExecutionContext host(hostPolicy);
    QVERIFY(host.getReport().availableThreads>=2);
    QVERIFY(host.getReport().actualThreads>=2);
    const int hostThreads=omp_get_max_threads(),hostDynamic=omp_get_dynamic();
    auto model=scene();auto setup=config();Adapter adapter(model,setup);
    bool sawSerial=false,fail=false;
    Driver driver(adapter,solver(),0,[&](double t){
        sawSerial=omp_get_max_threads()==1 && omp_get_dynamic()==0;
        if(fail && t>0)throw std::runtime_error("Host-policy restoration fault");
        return V(18);
    },"host-parallel-serial-native/v1");
    const auto first=driver.initialize(),repeat=driver.initialize();
    QVERIFY(sawSerial);QCOMPARE(first.integrityDigest,repeat.integrityDigest);
    QCOMPARE(omp_get_max_threads(),hostThreads);QCOMPARE(omp_get_dynamic(),hostDynamic);
    const auto step=driver.trialOnGrid(first,first.heldActuatorVoltages,1,.002);
    QCOMPARE(step.integrityDigest,driver.trialOnGrid(first,first.heldActuatorVoltages,1,.002).integrityDigest);
    QCOMPARE(omp_get_max_threads(),hostThreads);QCOMPARE(omp_get_dynamic(),hostDynamic);
    fail=true;QVERIFY(rejects([&]{driver.trialOnGrid(first,first.heldActuatorVoltages,1,.002);}));
    QCOMPARE(omp_get_max_threads(),hostThreads);QCOMPARE(omp_get_dynamic(),hostDynamic);
    fail=false;QCOMPARE(driver.trialOnGrid(first,first.heldActuatorVoltages,1,.002).integrityDigest,step.integrityDigest);
    QCOMPARE(omp_get_max_threads(),hostThreads);QCOMPARE(omp_get_dynamic(),hostDynamic);
} REPORT_LAYERED_EXCEPTIONS

QTEST_APPLESS_MAIN(TestStructuralPiezoLayeredDynamics)
#include "tst_structural_piezo_layered_dynamics.moc"

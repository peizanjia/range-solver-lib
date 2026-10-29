#include <QtTest>
#include <QTemporaryDir>
#include <QDataStream>
#include <rml_model.h>
#include <rbl_error.h>
#include "rstructural_piezo_axial_native_run.h"
#include "rstructural_piezo_axial_restart_codec.h"
#include "rstructural_piezo_element_operators.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
namespace {
using Run=RStructuralPiezoAxialNativeRun;
using Bundle=RStructuralPiezoAxialRuntimeCallbacks;
using Codec=RStructuralPiezoAxialRestartCodec;
using Role=RStructuralPiezoPortRole;
using Mode=RStructuralPiezoReadoutMode;
bool rejects(const std::function<void()> &f){try{f();}catch(...){return true;}return false;}
bool exact(double a,double b){return std::memcmp(&a,&b,8)==0;}
void component(RBoundaryCondition &bc,RVariableType type,double value,bool enabled=true)
{
    auto &c=bc.getComponent(bc.findComponentPosition(type)); c.setEnabled(enabled); c.clear(); if(enabled)c.add(0,value);
}
RBoundaryCondition support(uint axis,bool fixed,bool transverse=true)
{
    RBoundaryCondition bc(R_BOUNDARY_CONDITION_DISPLACEMENT);
    const RVariableType types[]={R_VARIABLE_DISPLACEMENT_X,R_VARIABLE_DISPLACEMENT_Y,R_VARIABLE_DISPLACEMENT_Z};
    for(uint a=0;a<3;++a) component(bc,types[a],0,fixed || (transverse && a!=axis));
    return bc;
}
RStructuralPiezoAxialMemberProperties material()
{
    auto m=RStructuralPiezoContinuumMaterial::isotropic(100,.25,7);
    m.piezoelectric=true; m.coupling[0]=2;
    m.permittivity[0]=m.permittivity[4]=m.permittivity[8]=.5;
    return rStructuralPiezoCompileAxialMember(m,1,1);
}
RModel fixture(uint axis=0,double q0=.04,bool transverse=true,bool explicitReference=true,
               double tension=1,double voltage=.7)
{
    RModel m; m.addNode(RNode(0,0,0)); m.addNode(axis==0?RNode(1,0,0):RNode(0,1,0));
    m.setNodePhysicalId(0,909); m.setNodePhysicalId(1,101);
    for(uint i=0;i<6;++i) {
        RElement e(R_ELEMENT_TRUSS1); e.setNodeId(0,0); e.setNodeId(1,1);
        m.addElement(e,true,i); m.setElementPhysicalId(i,901+i);
        RMaterial mat(RMaterial::Solid);
        RMaterialProperty young(RMaterialProperty::ModulusOfElasiticity); young.add(293.15,10); mat.add(young);
        RMaterialProperty density(RMaterialProperty::Density); density.add(293.15,1); mat.add(density);
        auto &line=m.getLine(i); line.setCrossArea(1); line.setMaterial(mat);
        if(transverse) line.addBoundaryCondition(support(axis,false));
    }
    for(uint node=0;node<2;++node) {
        RElement e(R_ELEMENT_POINT); e.setNodeId(0,node); m.addElement(e,true,node);
    }
    m.getPoint(0).addBoundaryCondition(support(axis,true));
    RInitialCondition initial(R_INITIAL_CONDITION_DISPLACEMENT);
    const RVariableType qs[]={R_VARIABLE_DISPLACEMENT_X,R_VARIABLE_DISPLACEMENT_Y,R_VARIABLE_DISPLACEMENT_Z};
    for(uint a=0;a<3;++a) {
        auto &c=initial.getComponent(initial.findComponentPosition(qs[a])); c.setEnabled(true); c.clear();
        c.add(0,a==axis?(explicitReference?.02:0):0);
    }
    m.getPoint(1).addInitialCondition(initial);
    const double charge=explicitReference?.03:0;
    const double force=140*q0+2*voltage-.4*(charge-.4*q0)/.5+std::max(0.0,tension+10*q0);
    RBoundaryCondition load(R_BOUNDARY_CONDITION_FORCE);
    const RVariableType fs[]={R_VARIABLE_FORCE_X,R_VARIABLE_FORCE_Y,R_VARIABLE_FORCE_Z};
    for(uint a=0;a<3;++a) component(load,fs[a],a==axis?force:0,a==axis);
    m.getPoint(1).addBoundaryCondition(load);
    for(uint i=0;i<m.getNElementGroups();++i) {
        auto *group=m.getElementGroupPtr(i);
        if(!group || !group->hasEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION)) continue;
        auto &g=group->getEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
        for(auto type:{R_VARIABLE_G_ACCELERATION_X,R_VARIABLE_G_ACCELERATION_Y,R_VARIABLE_G_ACCELERATION_Z}) {
            auto &c=g.getComponent(g.findComponentPosition(type)); c.setEnabled(true); c.clear(); c.add(0,0);
        }
    }
    m.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS)); m.getProblemSetup().setMassMatrixType(RMassMatrixType::Lumped);
    auto &d=m.getProblemSetup().getStructuralDynamicsSetup(); d.setEnabled(true); d.setTimeStep(.004); d.setNSteps(2); d.setOutputFrequency(1);
    auto &time=m.getTimeSolver(); time.setEnabled(true); time.setInputStartTime(.004); time.setInputTimeStepSize(.004);
    time.setInputNTimeSteps(2); time.setTimes({.004,.008});
    auto &truss=m.getProblemSetup().getStructuralTrussDynamicsSetup(); truss.enabled=true;
    truss.absoluteForceTolerance=1e-11; truss.relativeForceTolerance=1e-11; truss.linearRelativeTolerance=1e-12;
    truss.passiveMembers={{901,1}};
    const double coupling[]={2,-.4,.5,-.6}; const QString ids[]={"actuator","open","short","amplifier"};
    const auto property=material();
    for(uint i=0;i<4;++i) {
        RStructuralTrussActiveMember a; a.physicalElementId=902+i; a.portId=ids[i]; a.referenceLength=1;
        a.axialStiffness=i==0?property.k:10; a.coupling=coupling[i]; a.capacitance=.5;
        a.kind=i==1?RStructuralTrussElectricalKind::ConductorCharge:RStructuralTrussElectricalKind::Voltage;
        a.initialValue=i==0?voltage:(i==1?charge:0);
        if(i==0) { a.massSource=RStructuralTrussMassSource::FixedMemberMass; a.fixedMass=property.fixedMass; }
        truss.activeMembers.push_back(a);
    }
    RStructuralTrussCableMember cable; cable.physicalElementId=906; cable.initialReferenceLength=1; cable.initialTension=tension;
    truss.cableMembers.push_back(cable);
    auto &setup=m.getProblemSetup().getStructuralPiezoReadoutSetup();
    setup.ports={{"open",Role::Sensor,Mode::OpenCircuitVoltage},
                 {"amplifier",Role::Sensor,Mode::IdealChargeAmplifier,.25,-1,.012,1.3},
                 {"actuator",Role::ActuatorSensor,Mode::ActuatorCharge},
                 {"short",Role::Sensor,Mode::ShortCircuitCharge}};
    if(explicitReference) {
        setup.workpoint.available=true;
        RStructuralPiezoNodeWorkpoint node; node.physicalNodeId=101; node.q0[axis]=q0; node.f0[axis]=force;
        setup.workpoint.nodes={node,{909,{{0,0,0}},{{0,0,0}}}};
        setup.workpoint.ports={{"amplifier",0,.6*q0},{"short",0,-.5*q0},
                              {"open",(charge-.4*q0)/.5,charge},{"actuator",voltage,.5*voltage-2*q0}};
    }
    truss.bindToModel(m); return m;
}

RModel scene(bool configured=true,uint steps=12,double start=.004) {
    auto m=fixture();auto &d=m.getProblemSetup().getStructuralDynamicsSetup();d.setNSteps(steps);
    auto &t=m.getTimeSolver();t.setInputStartTime(start);t.setInputNTimeSteps(steps);
    std::vector<double> times;for(uint i=1;i<=steps;++i)times.push_back((start-.004)+i*.004);t.setTimes(times);
    m.getProblemSetup().getStructuralTrussDynamicsSetup().bindToModel(m);
    auto &s=m.getProblemSetup().getStructuralPiezoAxialSetup();s.enabled=true;
    s.recipeConfiguration=configured?QByteArray("{\"controller\":{\"test\":true},\"damping\":{\"test\":true},\"excitation\":{\"test\":true},\"family\":\"axialRod\"}"):QByteArray("{\"family\":\"axialRod\"}");s.bindToModel(m);return m;
}
QByteArray snapshot(const RStructuralTrussDynamicsState &p,const std::vector<RStructuralPiezoNativeReadoutSample> &y){
    QByteArray bytes;QDataStream a(&bytes,QIODevice::WriteOnly);a.setFloatingPointPrecision(QDataStream::DoublePrecision);
    a<<p.systemIdentity<<p.time;
    for(const auto *v:{&p.displacements,&p.velocities,&p.accelerations,&p.reactions}){a<<quint32(v->size());for(double x:*v)a<<x;}
    a<<quint32(y.size());for(const auto &v:y)a<<v.portId<<v.value<<v.deltaValue<<v.voltage<<v.charge;
    return bytes;
}
struct History {uint tick=0,completed=0;QByteArray witness;double applied=0,previous=0;};
QByteArray opaque(const History &h){QByteArray b;QDataStream a(&b,QIODevice::WriteOnly);a.setFloatingPointPrecision(QDataStream::DoublePrecision);a<<quint32(h.tick)<<quint32(h.completed)<<h.witness<<h.applied<<h.previous;return b;}
History parse(const QByteArray &b){QDataStream a(b);a.setFloatingPointPrecision(QDataStream::DoublePrecision);quint32 tick=0,completed=0;History h;a>>tick>>completed>>h.witness>>h.applied>>h.previous;h.tick=tick;h.completed=completed;if(a.status()!=QDataStream::Ok || !a.atEnd() || opaque(h)!=b)throw std::invalid_argument("bad history");return h;}
Bundle bundle(const RModel &m,uint stride=3){
    const auto &s=m.getProblemSetup().getStructuralPiezoAxialSetup();Bundle b;b.sourceFingerprint=s.computeSourceFingerprint(m);b.recipeConfiguration=s.recipeConfiguration;
    b.backendIdentity="native-axial-test-callbacks/v1";b.forcingIdentity="physical-time-force/v1";b.controllerIdentity="test-counter-observer/v1";b.controlEveryFESteps=stride;b.envelopes={{-1,1,2}};
    const size_t n=3*m.getNNodes();b.authoritativeFullC.assign(n*n,0);for(size_t i=0;i<n;++i)b.authoritativeFullC[i*n+i]=.2;
    const double t0=m.getTimeSolver().getInputStartTime()-.004;b.perturbationForce=[n,t0](double t){std::vector<double> f(n);f[3]=.03+.2*(t-t0);return f;};
    RModel definition(m);definition.getProblemSetup().setRestart(false);RStructuralTrussDynamicsModelAdapter adapter(definition);auto system=adapter.system();system.physicalDamping=b.authoritativeFullC;RStructuralTrussDynamics driver(system);const auto seed=adapter.sourceInitialState();auto inputs=adapter.initialCompleteInputs();auto perturb=b.perturbationForce(t0);for(size_t i=0;i<n;++i)inputs.externalForces[i]+=perturb[i];const auto initial=driver.initialize(t0,seed.displacements,seed.velocities,inputs);RStructuralPiezoNativeReadout readout(adapter,adapter.piezoReadoutSetup(),system);
    b.initialOpaqueControllerState=opaque({0,0,snapshot(initial,readout.observe(initial)),.7,.7});
    b.propose=[](const auto &physical,const auto &readouts,uint completed,uint tick,const QByteArray &old){auto h=parse(old);if(h.tick!=tick || h.completed!=completed || h.witness!=snapshot(physical,readouts))throw std::invalid_argument("proposal boundary mismatch");++h.tick;RStructuralPiezoAxialControllerProposal p;p.requestedVoltages={tick%2?-2.:2.};p.nextOpaqueState=opaque(h);p.lastFrameBytes=QByteArray::number(completed);return p;};
    b.postFE=[](const auto &physical,const auto &readouts,uint completed,uint tick,const auto &applied,const QByteArray &old){auto h=parse(old);if(h.tick!=tick || h.completed+1!=completed || applied.size()!=1)throw std::invalid_argument("observer phase mismatch");h.completed=completed;h.previous=h.applied;h.applied=applied[0];h.witness=snapshot(physical,readouts);return opaque(h);};
    b.validateOpaque=[](const QByteArray &old,uint completed,uint tick,uint){const auto h=parse(old);if(h.tick!=tick || h.completed!=completed)throw std::invalid_argument("controller counter mismatch");};
    b.validateAcceptedBoundary=[](const RStructuralPiezoAxialRunState &s){const auto h=parse(s.opaqueControllerState);if(h.completed!=s.completedFE || h.tick!=s.controlTick || h.witness!=snapshot(s.physical,s.readouts) || s.lastApplied.size()!=1 || !exact(h.applied,s.lastApplied[0]) || (s.completedFE && !exact(h.previous,s.actions[0].previousApplied)))throw std::invalid_argument("duplicate accepted history mismatch");};return b;
}
RModel advance(Run &run,const RModel &source){auto step=run.prepare(source);run.validatePrepared(*step);RModel result(step->candidateModel());run.commit(*step);return result;}
}
#define REPORT_EXCEPTIONS catch(const RError &e){QFAIL(qPrintable(e.getMessage()));}catch(const std::exception &e){QFAIL(e.what());}
class TestStructuralPiezoAxialNativeRun:public QObject {
    Q_OBJECT
private slots:
    void diskEveryPhase_data();
    void diskEveryPhase();
    void actualDampingForceAndSensors();
    void rollbackAndTamper();
    void callbackFailures();
    void forcingAliasAfterSnapshotRejected();
    void validateOpaqueAliasVisibleToBoundaryCallback();
    void declaredCallbacksRequired();
    void noControlPreviousAppliedCertified();
    void semanticRestoreAndCodec();
    void originalThousandStepGrid();
    void sensorOnly();
    void multiActuatorOrderAndBounds();
};
void TestStructuralPiezoAxialNativeRun::diskEveryPhase_data(){QTest::addColumn<QString>("extension");QTest::addColumn<uint>("stride");QTest::addColumn<double>("start");for(const auto &ext:{QString("ram"),QString("rbm")})for(uint n:{2u,3u,4u})for(double start:{.004,1000.004})QTest::newRow(qPrintable(QString("%1-N%2-start%3").arg(ext).arg(n).arg(start,0,'g',17)))<<ext<<n<<start;}
void TestStructuralPiezoAxialNativeRun::diskEveryPhase() try {
    QFETCH(QString,extension);QFETCH(uint,stride);QFETCH(double,start);const auto source=scene(true,12,start);auto cb=bundle(source,stride);Run full(source,cb);
    while(full.hasRemainingSteps())advance(full,source);const auto expected=Codec::encode(full.acceptedRuntime());QTemporaryDir dir;QVERIFY(dir.isValid());
    for(uint split=0;split<=12;++split){Run first(source,cb);for(uint i=0;i<split;++i)advance(first,source);auto model=first.acceptedModel(source);const auto path=dir.filePath(QString("split-%1.%2").arg(split).arg(extension));const auto written=model.nativeWritePaths(path,false).first();model.write(path,false);RModel disk;disk.read(written);
        QCOMPARE(disk.getProblemSetup().getStructuralPiezoAxialRestart().runtimeState,Codec::encode(first.acceptedRuntime()));QVERIFY(!disk.getProblemSetup().getStructuralTrussDynamicsRestartData().hasNonDefaultData());QVERIFY(!disk.getProblemSetup().getStructuralNativeControlState().hasNonDefaultData());QCOMPARE(disk.getTimeSolver().getInputStartTime(),start);
        disk.getProblemSetup().setRestart(true);Run restored(disk,bundle(disk,stride));QCOMPARE(Codec::encode(restored.acceptedRuntime()),Codec::encode(first.acceptedRuntime()));while(restored.hasRemainingSteps())advance(restored,disk);QCOMPARE(Codec::encode(restored.acceptedRuntime()),expected);
    }
} REPORT_EXCEPTIONS
void TestStructuralPiezoAxialNativeRun::actualDampingForceAndSensors() try {
    const auto source=scene();auto cb=bundle(source);Run run(source,cb);QCOMPARE(run.acceptedState().systemIdentity,RStructuralTrussDynamics(run.system()).systemIdentity());QVERIFY(run.acceptedState().systemIdentity!=run.adapter().driver().systemIdentity());
    QCOMPARE(run.acceptedState().appliedInputs.externalForces[3],run.adapter().initialCompleteInputs().externalForces[3]+.03);
    const auto initial=Codec::encode(run.acceptedRuntime());auto prepared=run.prepare(source);QCOMPARE(Codec::encode(run.acceptedRuntime()),initial);auto direct=run.adapter().initialCompleteInputs();direct.memberInputs[run.adapter().maps().memberByPortId.at("actuator")]=RStructuralTrussVoltage{prepared->runtime().lastApplied[0]};direct.externalForces=run.acceptedState().appliedInputs.externalForces;direct.rightExternalForces=direct.externalForces;direct.rightExternalForces[3]=run.adapter().initialCompleteInputs().externalForces[3]+(.03+.2*.004);
    const auto physical=RStructuralTrussDynamics(run.system()).trialStepOnGrid(run.acceptedState(),direct,0,1,.004,run.adapter().options()).candidate;QCOMPARE(snapshot(physical,prepared->runtime().readouts),snapshot(prepared->result(),prepared->runtime().readouts));
    QVERIFY(prepared->result().cumulativeDampingDissipation>0);QCOMPARE(prepared->runtime().actions.size(),size_t(1));QCOMPARE(prepared->runtime().actions[0].portId,QString("actuator"));QVERIFY(prepared->runtime().actions[0].amplitudeClipped);QVERIFY(prepared->runtime().actions[0].rateClipped);
    run.validatePrepared(*prepared);run.commit(*prepared);for(uint i=0;i<2;++i){advance(run,source);QCOMPARE(run.acceptedRuntime().controlTick,1u);QCOMPARE(run.acceptedRuntime().heldRequested[0],2.);}
    advance(run,source);QCOMPARE(run.acceptedRuntime().controlTick,2u);QCOMPARE(run.acceptedRuntime().heldRequested[0],-2.);
    for(const auto &p:run.acceptedRuntime().readouts){if(p.portId=="open")QCOMPARE(p.charge,.03);if(p.portId=="short" || p.portId=="amplifier")QCOMPARE(p.voltage,0.);}
    const auto out=run.acceptedModel(source);QVERIFY(out.findVariable(R_VARIABLE_STRESS_NORMAL)!=RConstants::eod);QVERIFY(out.findVariable(R_VARIABLE_VELOCITY)!=RConstants::eod);QVERIFY(out.findVariable(R_VARIABLE_ELECTRIC_POTENTIAL)==RConstants::eod);
} REPORT_EXCEPTIONS
void TestStructuralPiezoAxialNativeRun::rollbackAndTamper() try {
    const auto source=scene();Run run(source,bundle(source));auto a=run.prepare(source),b=run.prepare(source);QCOMPARE(Codec::encode(a->runtime()),Codec::encode(b->runtime()));QCOMPARE(run.completedSteps(),0u);
    const auto check=[&](const std::function<void(RModel &)> &mutate){RModel changed(a->candidateModel());mutate(changed);QVERIFY(rejects([&]{run.validatePrepared(*a,changed);}));};
    check([](RModel &m){m.getVariable(0).setValue(0,0,1);});check([](RModel &m){m.getVariable(0).getVariableData().setScale(2);});check([](RModel &m){m.getTimeSolver().setComputedTime(.1);});check([](RModel &m){m.getTimeSolver().getTimes()[0]=.002;});check([](RModel &m){auto &b=m.getProblemSetup().getStructuralPiezoAxialRestart().runtimeState;b[15]=char(uchar(b[15])^1);});check([](RModel &m){m.getProblemSetup().getStructuralDynamicsSetup().setTimeStep(.002);});
    run.validatePrepared(*b);run.commit(*b);QVERIFY(rejects([&]{run.validatePrepared(*a);}));QVERIFY(rejects([&]{run.validatePrepared(*b);}));QCOMPARE(run.completedSteps(),1u);auto next=run.prepare(source);run.markNonRetryable();QVERIFY(rejects([&]{run.validatePrepared(*next);}));QVERIFY(rejects([&]{run.prepare(source);}));
} REPORT_EXCEPTIONS
void TestStructuralPiezoAxialNativeRun::callbackFailures() try {
    const auto source=scene();for(uint which=0;which<4;++which){auto cb=bundle(source);bool fail=false;
        const auto propose=cb.propose;cb.propose=[&,propose](const auto &p,const auto &y,uint c,uint t,const QByteArray &o){if(fail && which==0)throw std::runtime_error("proposal failure");return propose(p,y,c,t,o);};
        const auto post=cb.postFE;cb.postFE=[&,post](const auto &p,const auto &y,uint c,uint t,const auto &v,const QByteArray &o){if(fail && which==1)throw std::runtime_error("observer failure");return post(p,y,c,t,v,o);};
        const auto forcing=cb.perturbationForce;cb.perturbationForce=[&,forcing](double t){if(fail && which==2 && t>0)throw std::runtime_error("forcing failure");return forcing(t);};
        const auto validate=cb.validateAcceptedBoundary;cb.validateAcceptedBoundary=[&,validate](const auto &s){if(fail && which==3 && s.completedFE)throw std::runtime_error("certificate failure");validate(s);};
        Run run(source,cb);const auto before=Codec::encode(run.acceptedRuntime());fail=true;QVERIFY(rejects([&]{run.prepare(source);}));QCOMPARE(Codec::encode(run.acceptedRuntime()),before);fail=false;advance(run,source);QCOMPARE(run.completedSteps(),1u);
    }
} REPORT_EXCEPTIONS
void TestStructuralPiezoAxialNativeRun::forcingAliasAfterSnapshotRejected() try {
    const auto source=scene();auto cb=bundle(source);RStructuralPiezoAxialRunState *aliasedState=nullptr;bool mutate=false;uint forceCalls=0;
    // certifyRuntime snapshots before inputs() invokes this callback; keep the factory path untouched.
    const auto forcing=cb.perturbationForce;cb.perturbationForce=[&,forcing](double t){++forceCalls;if(mutate && aliasedState)aliasedState->controllerIdentity="test-counter-observer/post-callback-v1";return forcing(t);};
    Run run(source,cb);auto unmutated=run.acceptedRuntime();aliasedState=&unmutated;forceCalls=0;run.certifyAcceptedRuntime(unmutated);QVERIFY(forceCalls>0);QCOMPARE(unmutated.controllerIdentity,cb.controllerIdentity);
    auto mutated=run.acceptedRuntime();aliasedState=&mutated;mutate=true;forceCalls=0;QVERIFY(rejects([&]{run.certifyAcceptedRuntime(mutated);}));QVERIFY(forceCalls>0);QCOMPARE(mutated.controllerIdentity,QString("test-counter-observer/post-callback-v1"));
} REPORT_EXCEPTIONS
void TestStructuralPiezoAxialNativeRun::validateOpaqueAliasVisibleToBoundaryCallback() try {
    const auto source=scene();auto cb=bundle(source);RStructuralPiezoAxialRunState *aliasedState=nullptr;bool mutateAfterOpaque=false,compareAcceptedBytes=false,opaqueInvoked=false,boundaryInvoked=false;QByteArray expectedAcceptedBytes;
    // The completed step gives this arbitrary controller-frame field valid accepted bytes.
    const auto validateOpaque=cb.validateOpaque;cb.validateOpaque=[&](const QByteArray &old,uint completed,uint tick,uint phase){validateOpaque(old,completed,tick,phase);opaqueInvoked=true;if(mutateAfterOpaque && aliasedState)aliasedState->lastFrameBytes.append("-opaque-mutated");};
    const auto validateBoundary=cb.validateAcceptedBoundary;cb.validateAcceptedBoundary=[&](const RStructuralPiezoAxialRunState &state){boundaryInvoked=true;if(compareAcceptedBytes && Codec::encode(state)!=expectedAcceptedBytes)throw std::invalid_argument("accepted state changed after opaque validation");validateBoundary(state);};
    Run run(source,cb);advance(run,source);auto unmutated=run.acceptedRuntime();QCOMPARE(unmutated.completedFE,1u);QVERIFY(!unmutated.lastFrameBytes.isEmpty());
    aliasedState=&unmutated;expectedAcceptedBytes=Codec::encode(unmutated);compareAcceptedBytes=true;opaqueInvoked=boundaryInvoked=false;run.certifyAcceptedRuntime(unmutated);QVERIFY(opaqueInvoked);QVERIFY(boundaryInvoked);QCOMPARE(Codec::encode(unmutated),expectedAcceptedBytes);
    auto mutated=run.acceptedRuntime();aliasedState=&mutated;expectedAcceptedBytes=Codec::encode(mutated);mutateAfterOpaque=true;opaqueInvoked=boundaryInvoked=false;const auto frameBefore=mutated.lastFrameBytes;
    QVERIFY(rejects([&]{run.certifyAcceptedRuntime(mutated);}));QVERIFY(opaqueInvoked);QVERIFY(boundaryInvoked);QVERIFY(mutated.lastFrameBytes!=frameBefore);QCOMPARE(mutated.lastFrameBytes,frameBefore+QByteArray("-opaque-mutated"));
} REPORT_EXCEPTIONS
void TestStructuralPiezoAxialNativeRun::declaredCallbacksRequired() try {
    const auto source=scene();auto control=bundle(source);
    control.propose={};control.postFE={};control.validateOpaque={};control.validateAcceptedBoundary={};
    control.initialOpaqueControllerState.clear();control.controllerIdentity.clear();control.envelopes.clear();
    QVERIFY(rejects([&]{Run invalid(source,control);}));
    auto force=bundle(source);force.perturbationForce={};force.forcingIdentity.clear();
    QVERIFY(rejects([&]{Run invalid(source,force);}));
} REPORT_EXCEPTIONS
void TestStructuralPiezoAxialNativeRun::noControlPreviousAppliedCertified() try {
    const auto source=scene(false);Run run(source);advance(run,source);advance(run,source);
    auto model=run.acceptedModel(source);model.getProblemSetup().setRestart(true);
    auto state=run.acceptedRuntime();QVERIFY(!state.actions.empty());
    state.actions[0].previousApplied+=.125;
    model.getProblemSetup().getStructuralPiezoAxialRestart().runtimeState=Codec::encode(state);
    QVERIFY(rejects([&]{Run invalid(model);}));
} REPORT_EXCEPTIONS
void TestStructuralPiezoAxialNativeRun::semanticRestoreAndCodec() try {
    const auto source=scene();const auto cb=bundle(source);Run run(source,cb);advance(run,source);auto model=run.acceptedModel(source);model.getProblemSetup().setRestart(true);
    const auto reject=[&](const std::function<void(RStructuralPiezoAxialRunState &)> &mutate){auto state=run.acceptedRuntime();mutate(state);model.getProblemSetup().getStructuralPiezoAxialRestart().runtimeState=Codec::encode(state);QVERIFY(rejects([&]{Run invalid(model,cb);}));};
    reject([](auto &s){s.opaqueControllerState="bad";});reject([](auto &s){s.initialPhysicalReference.displacements[3]+=.01;});reject([](auto &s){s.physical.initialTotalEnergy+=.1;});reject([](auto &s){s.physical.displacements[3]+=.01;});reject([](auto &s){s.readouts[0].value+=1;});reject([](auto &s){++s.controlTick;});
    reject([](auto &s){auto h=parse(s.opaqueControllerState);h.previous+=.001;s.opaqueControllerState=opaque(h);});
    model=run.acceptedModel(source);model.getProblemSetup().setRestart(true);auto changed=cb;changed.authoritativeFullC[0]+=.1;QVERIFY(rejects([&]{Run invalid(model,changed);}));changed=cb;changed.backendIdentity="wrong";QVERIFY(rejects([&]{Run invalid(model,changed);}));changed=cb;changed.recipeConfiguration+=" ";QVERIFY(rejects([&]{Run invalid(model,changed);}));
    model.getVariable(0).setValue(0,0,123);Run displayIgnored(model,cb);QCOMPARE(Codec::encode(displayIgnored.acceptedRuntime()),Codec::encode(run.acceptedRuntime()));
    auto state=run.acceptedRuntime();state.heldRequested[0]=-0.;const auto bytes=Codec::encode(state);QVERIFY(std::signbit(Codec::decode(bytes).heldRequested[0]));QCOMPARE(Codec::encode(Codec::decode(bytes)),bytes);QVERIFY(rejects([&]{Codec::decode(bytes.left(bytes.size()-1));}));QVERIFY(rejects([&]{Codec::decode(bytes+"x");}));state.dt=std::numeric_limits<double>::infinity();QVERIFY(rejects([&]{Codec::encode(state);}));
    changed=cb;changed.validateAcceptedBoundary={};QVERIFY(rejects([&]{Run invalid(source,changed);}));
} REPORT_EXCEPTIONS
void TestStructuralPiezoAxialNativeRun::originalThousandStepGrid() try {
    const auto source=scene(true,1000,1000.004);const auto cb=bundle(source,4);Run run(source,cb);while(run.hasRemainingSteps()){advance(run,source);QVERIFY(exact(run.acceptedState().time,1000.+double(run.completedSteps())*.004));const auto &a=run.acceptedRuntime().actions[0];const double h=run.acceptedState().time-(1000.+double(run.completedSteps()-1)*.004);QVERIFY(exact(a.applied,std::clamp(std::clamp(a.requested,-1.,1.),a.previousApplied-2*h,a.previousApplied+2*h)));}
    auto model=run.acceptedModel(source);model.getProblemSetup().setRestart(true);Run restored(model,bundle(model,4));QCOMPARE(Codec::encode(restored.acceptedRuntime()),Codec::encode(run.acceptedRuntime()));QVERIFY(!restored.hasRemainingSteps());
} REPORT_EXCEPTIONS
void TestStructuralPiezoAxialNativeRun::sensorOnly() try {
    auto source=scene(false);auto &catalog=source.getProblemSetup().getStructuralPiezoReadoutSetup();for(auto &port:catalog.ports)if(port.portId=="actuator"){port.role=Role::Sensor;port.readout=Mode::ShortCircuitCharge;}
    auto &truss=source.getProblemSetup().getStructuralTrussDynamicsSetup();for(auto &member:truss.activeMembers)if(member.portId=="actuator")member.initialValue=0;
    catalog.workpoint.available=false;catalog.workpoint.nodes.clear();catalog.workpoint.ports.clear();truss.bindToModel(source);source.getProblemSetup().getStructuralPiezoAxialSetup().bindToModel(source);Run run(source);QVERIFY(run.acceptedRuntime().lastApplied.empty());advance(run,source);QVERIFY(run.acceptedRuntime().actions.empty());
} REPORT_EXCEPTIONS
void TestStructuralPiezoAxialNativeRun::multiActuatorOrderAndBounds() try {
    auto source=scene();auto &catalog=source.getProblemSetup().getStructuralPiezoReadoutSetup();
    for(auto &p:catalog.ports)if(p.portId=="short"){p.role=Role::ActuatorSensor;p.readout=Mode::ActuatorCharge;}
    std::reverse(catalog.ports.begin(),catalog.ports.end());source.getProblemSetup().getStructuralTrussDynamicsSetup().bindToModel(source);source.getProblemSetup().getStructuralPiezoAxialSetup().bindToModel(source);
    Bundle cb;cb.sourceFingerprint=source.getProblemSetup().getStructuralPiezoAxialSetup().sourceFingerprint;cb.recipeConfiguration=source.getProblemSetup().getStructuralPiezoAxialSetup().recipeConfiguration;cb.backendIdentity="test-two-actuators/v1";cb.controllerIdentity="test-held-two-actions/v1";cb.controlEveryFESteps=3;cb.envelopes={{-.1,.8,3},{-1,1,5}};
    cb.forcingIdentity="test-declared-zero-time-force/v1";
    cb.perturbationForce=[n=size_t(source.getNNodes())*3](double){return std::vector<double>(n,0.);};
    cb.propose=[](const auto &,const auto &,uint,uint,const QByteArray &){RStructuralPiezoAxialControllerProposal p;p.requestedVoltages={2,-2};return p;};
    cb.validateAcceptedBoundary=[](const auto &s){
        if(!s.opaqueControllerState.isEmpty() || s.lastApplied.size()!=2)throw std::invalid_argument("invalid two-actuator history");
        std::vector<double> expected={.7,0},previous=expected;
        for(uint i=1;i<=s.completedFE;++i){previous=expected;const double h=(s.initialTime+i*s.dt)-(i==1?s.initialTime:s.initialTime+(i-1)*s.dt);expected[0]=std::clamp(.8,expected[0]-3*h,expected[0]+3*h);expected[1]=std::clamp(-1.,expected[1]-5*h,expected[1]+5*h);}
        for(size_t i=0;i<2;++i)if(!exact(expected[i],s.lastApplied[i]) || (s.completedFE && !exact(previous[i],s.actions[i].previousApplied)))throw std::invalid_argument("two-actuator previous-applied mismatch");
    };
    Run run(source,cb);while(run.hasRemainingSteps()){advance(run,source);const auto &r=run.acceptedRuntime();QCOMPARE(r.actions.size(),size_t(2));QCOMPARE(r.actions[0].portId,QString("actuator"));QCOMPARE(r.actions[1].portId,QString("short"));QCOMPARE(r.controlTick,(r.completedFE+2)/3);}
} REPORT_EXCEPTIONS
QTEST_APPLESS_MAIN(TestStructuralPiezoAxialNativeRun)
#include "tst_structural_piezo_axial_native_run.moc"

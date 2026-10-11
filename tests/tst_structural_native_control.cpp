#include <QtTest>
#include <QCryptographicHash>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <rbl_error.h>
#include <rml_model.h>
#include "rstructural_native_control.h"
#include "rstructural_truss_native_run.h"
#include "rstructural_truss_dynamics_restart_adapter.h"

namespace
{
using Adapter=RStructuralTrussDynamicsModelAdapter;
using Restart=RStructuralTrussDynamicsRestartAdapter;
using Run=RStructuralTrussNativeRun;
using Measurement=RStructuralNativeMeasurementQuantity;
using Action=RStructuralNativeActionQuantity;
template<class F> bool rejects(F f)
{ try { f(); } catch(const RError &) { return true; } catch(const std::exception &) { return true; } return false; }
bool near(double a,double b,double tolerance=2e-10)
{ return std::abs(a-b)<=tolerance*(1+std::abs(b)); }
void component(RBoundaryCondition &bc,RVariableType type,double value,bool enabled=true)
{
    auto &c=bc.getComponent(bc.findComponentPosition(type)); c.setEnabled(enabled); c.clear(); if(enabled)c.add(0,value);
}
RBoundaryCondition support(bool x,bool y,bool z)
{
    RBoundaryCondition bc(R_BOUNDARY_CONDITION_DISPLACEMENT);
    component(bc,R_VARIABLE_DISPLACEMENT_X,0,x); component(bc,R_VARIABLE_DISPLACEMENT_Y,0,y); component(bc,R_VARIABLE_DISPLACEMENT_Z,0,z);
    return bc;
}
void point(RModel &m,uint node,uint group)
{ RElement e(R_ELEMENT_POINT); e.setNodeId(0,node); m.addElement(e,true,group); }
void initial(RModel &m)
{
    RInitialCondition q(R_INITIAL_CONDITION_DISPLACEMENT),v(R_INITIAL_CONDITION_VELOCITY);
    const RVariableType qs[]={R_VARIABLE_DISPLACEMENT_X,R_VARIABLE_DISPLACEMENT_Y,R_VARIABLE_DISPLACEMENT_Z};
    const RVariableType vs[]={R_VARIABLE_VELOCITY_X,R_VARIABLE_VELOCITY_Y,R_VARIABLE_VELOCITY_Z};
    for(uint i=0;i<3;++i)
    {
        auto &a=q.getComponent(q.findComponentPosition(qs[i])); a.setEnabled(true); a.clear(); a.add(0,i==0?.035:0.0);
        auto &b=v.getComponent(v.findComponentPosition(vs[i])); b.setEnabled(true); b.clear(); b.add(0,i==0?.08:0.0);
    }
    m.getPoint(1).addInitialCondition(q); m.getPoint(1).addInitialCondition(v);
}
RStructuralNativeControlSensor nodeSensor(const QString &id,Measurement quantity,uint node,uint axis=0)
{ RStructuralNativeControlSensor s; s.sensorId=id; s.quantity=quantity; s.physicalNodeId=node; s.axis=axis; return s; }
RStructuralNativeControlSensor portSensor(const QString &id,Measurement quantity,const QString &port)
{ RStructuralNativeControlSensor s; s.sensorId=id; s.quantity=quantity; s.portId=port; return s; }
RStructuralNativeControlSensor memberSensor(const QString &id,uint member)
{ RStructuralNativeControlSensor s; s.sensorId=id; s.quantity=Measurement::MemberAxialForce; s.physicalMemberId=member; return s; }
RStructuralNativeControlAction forceAction()
{ RStructuralNativeControlAction a; a.actionId="force increment"; a.physicalNodeId=101; return a; }
void bind(RModel &m) { m.getProblemSetup().getStructuralTrussDynamicsSetup().bindToModel(m); }
RModel fixture(bool mixed=false,uint steps=24,double dt=.004)
{
    RModel m; m.addNode(RNode(0,0,0)); m.addNode(RNode(1,0,0));
    // Physical IDs deliberately disagree with sorted order and array indices.
    m.setNodePhysicalId(0,909); m.setNodePhysicalId(1,101);
    const uint ids[]={901,503,777,609};
    for(uint i=0;i<(mixed?4u:1u);++i)
    {
        RElement e(R_ELEMENT_TRUSS1); e.setNodeId(0,0); e.setNodeId(1,1); m.addElement(e,true,i); m.setElementPhysicalId(i,ids[i]);
        RMaterial material(RMaterial::Solid);
        RMaterialProperty young(RMaterialProperty::ModulusOfElasiticity); young.add(293.15,mixed?10:40); material.add(young);
        RMaterialProperty density(RMaterialProperty::Density); density.add(293.15,mixed?1:4); material.add(density);
        auto &line=m.getLine(i); line.setCrossArea(1); line.setMaterial(material); line.addBoundaryCondition(support(false,true,true));
    }
    point(m,0,0); m.getPoint(0).addBoundaryCondition(support(true,true,true)); point(m,1,1); initial(m);
    RBoundaryCondition load(R_BOUNDARY_CONDITION_FORCE); component(load,R_VARIABLE_FORCE_X,1.2);
    component(load,R_VARIABLE_FORCE_Y,0,false); component(load,R_VARIABLE_FORCE_Z,0,false); m.getPoint(1).addBoundaryCondition(load);
    for(uint i=0;i<m.getNElementGroups();++i)
    {
        auto *group=m.getElementGroupPtr(i);
        if(!group || !group->hasEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION)) continue;
        auto &gravity=group->getEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
        for(auto type:{R_VARIABLE_G_ACCELERATION_X,R_VARIABLE_G_ACCELERATION_Y,R_VARIABLE_G_ACCELERATION_Z})
        { auto &c=gravity.getComponent(gravity.findComponentPosition(type)); c.setEnabled(true); c.clear(); c.add(0,0); }
    }
    m.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS)); m.getProblemSetup().setMassMatrixType(RMassMatrixType::Lumped);
    auto &d=m.getProblemSetup().getStructuralDynamicsSetup(); d.setEnabled(true); d.setTimeStep(dt); d.setNSteps(steps); d.setOutputFrequency(1);
    auto &time=m.getTimeSolver(); time.setEnabled(true); time.setInputStartTime(dt); time.setInputTimeStepSize(dt); time.setInputNTimeSteps(steps);
    std::vector<double> times; for(uint i=1;i<=steps;++i)times.push_back(double(i)*dt); time.setTimes(times);
    auto &truss=m.getProblemSetup().getStructuralTrussDynamicsSetup(); truss.enabled=true;
    truss.absoluteForceTolerance=1e-11; truss.relativeForceTolerance=1e-11; truss.linearRelativeTolerance=1e-12;
    truss.passiveMembers={{901,1}};
    if(mixed)
    {
        for(uint i=0;i<2;++i)
        {
            RStructuralTrussActiveMember a; a.physicalElementId=ids[i+1]; a.portId=i==0?"voltage independent":"charge independent";
            a.referenceLength=1; a.axialStiffness=10; a.coupling=.2; a.capacitance=.5; a.initialValue=i==0?.7:.15;
            a.kind=i==0?RStructuralTrussElectricalKind::Voltage:RStructuralTrussElectricalKind::ConductorCharge;
            truss.activeMembers.push_back(a);
        }
        truss.cableMembers={{609,1,1}};
    }
    auto &spec=m.getProblemSetup().getStructuralNativeControlSpec(); spec.source=RStructuralNativeControlSource::LocalV3;
    RStructuralNativeControl::stampBuiltinIdentity(spec);
    spec.sensors={nodeSensor("displacement",Measurement::NodeDisplacement,101),
                  nodeSensor("left acceleration",Measurement::NodeAcceleration,101)};
    spec.actions={forceAction()}; spec.piLoops={{"force increment","displacement",.02,5,1.5,.3,.012}};
    bind(m); return m;
}
QString seal(const Run &run)
{
    const auto record=Restart::build(run.adapter(),run.acceptedState(),run.completedSteps(),run.nextCommandIndex(),run.acceptedControlState());
    return Restart::stateFingerprint(record,run.acceptedControlState());
}
void advance(Run &run,const RModel &source,uint until)
{
    while(run.completedSteps()<until)
    { auto step=run.prepare(source); run.validatePrepared(*step); run.commit(*step); }
}
}

#define REPORT_NATIVE_EXCEPTIONS \
    catch(const RError &error) { QFAIL(qPrintable(QString("RError: %1").arg(error.getMessage()))); } \
    catch(const std::exception &error) { QFAIL(error.what()); }

class TestStructuralNativeControl:public QObject
{
    Q_OBJECT
private slots:
    void independentScalarSampleHoldPI();
    void physicalMeasurementsAndOwnership();
    void movingSupportReactionIncludesInertia();
    void identityAndUnsupportedSpec_data();
    void identityAndUnsupportedSpec();
    void completeStateRestoreTamper_data();
    void completeStateRestoreTamper();
    void splitRestoreEveryPhase();
    void proposalAndTrialFailureRemainPure();
    void limiterPhysicalIntervalAndRetry();
    void defaultFingerprintAndOldRouteGuards();
};

void TestStructuralNativeControl::independentScalarSampleHoldPI() try
{
    const auto source=fixture(false,100); Run run(source);
    double q=.035,v=.08,a=(1.2-40*q)/2,integral=.012;
    for(uint n=0;n<100;++n)
    {
        const double leftA=a,error=.02-q,requested=.3+5*error+1.5*integral;
        integral+=.004*error;
        // Independent closed-form undamped scalar Newmark with the command
        // held over the interval and its left acceleration refreshed after sampling.
        const double f=1.2+requested,h=run.adapter().targetTime(n+1)-run.adapter().targetTime(n);
        const double refreshedA=(f-40*q)/2;
        const double predictor=q+h*v+.25*h*h*refreshedA;
        const double nextQ=(predictor+.25*h*h*f/2)/(1+.25*h*h*40/2);
        const double nextA=(f-40*nextQ)/2;
        const double nextV=v+.5*h*(refreshedA+nextA);
        auto step=run.prepare(source);
        const auto &state=step->candidateModel().getProblemSetup().getStructuralNativeControlState();
        QVERIFY(near(state.lastFrame.samples[0].value,q)); QVERIFY(near(state.lastFrame.samples[1].value,leftA));
        QVERIFY(near(state.lastActions[0].requested,requested)); QCOMPARE(state.lastActions[0].conditioned,state.lastActions[0].requested);
        if(n==0) QVERIFY(std::abs(state.lastFrame.samples[1].value-refreshedA)>1e-4);
        run.validatePrepared(*step); run.commit(*step);
        QVERIFY(near(run.acceptedState().displacements[3],nextQ)); QVERIFY(near(run.acceptedState().velocities[3],nextV));
        QVERIFY(near(run.acceptedState().accelerations[3],nextA));
        q=nextQ; v=nextV; a=nextA;
    }
    QCOMPARE(run.acceptedControlState().completedFE,quint64(100)); QCOMPARE(run.acceptedControlState().completedControlTicks,quint64(100));
} REPORT_NATIVE_EXCEPTIONS

void TestStructuralNativeControl::physicalMeasurementsAndOwnership() try
{
    auto source=fixture(true);
    auto &spec=source.getProblemSetup().getStructuralNativeControlSpec();
    spec.sensors={nodeSensor("u",Measurement::NodeDisplacement,101),nodeSensor("v",Measurement::NodeVelocity,101),
                  nodeSensor("a",Measurement::NodeAcceleration,101),memberSensor("rod N",901),memberSensor("cable N",609),
                  portSensor("V record echo",Measurement::PortVoltage,"voltage independent"),
                  portSensor("Q deformation",Measurement::PortCharge,"voltage independent"),
                  portSensor("V deformation",Measurement::PortVoltage,"charge independent"),
                  portSensor("Q record",Measurement::PortCharge,"charge independent"),
                  nodeSensor("support Rx",Measurement::SupportReaction,909)};
    auto voltage=RStructuralNativeControlAction{}; voltage.actionId="voltage"; voltage.quantity=Action::Voltage;
    voltage.portId="voltage independent"; voltage.initialCommand=.7;
    auto charge=RStructuralNativeControlAction{}; charge.actionId="charge"; charge.quantity=Action::Charge;
    charge.portId="charge independent"; charge.initialCommand=.15;
    auto cable=RStructuralNativeControlAction{}; cable.actionId="length"; cable.quantity=Action::CableReferenceLength;
    cable.physicalMemberId=609; cable.initialCommand=1;
    spec.actions={forceAction(),voltage,charge,cable};
    // PI envelope order is independent of action catalog order.
    spec.piLoops={{"length","u",0,0,0,.99,0},{"charge","u",0,0,0,.2,0},
                  {"force increment","u",0,0,0,.4,0},{"voltage","u",0,0,0,.8,0}};
    bind(source);
    const uint stale=source.addVariable(R_VARIABLE_DISPLACEMENT); auto &display=source.getVariable(stale);
    display.setApplyType(R_VARIABLE_APPLY_NODE); display.resize(3,2); display.setValue(0,1,123456);
    Adapter adapter(source); RStructuralNativeControl provider(adapter); const auto physical=adapter.initialize();
    const auto frame=provider.sample(physical,0);
    const auto fresh=RStructuralTrussAssembly::evaluate(adapter.system().referenceNodes,physical.displacements,adapter.system().members,physical.appliedInputs.memberInputs);
    const auto value=[&](const QString &id){for(const auto &s:frame.samples)if(s.sensorId==id)return s.value;return std::numeric_limits<double>::quiet_NaN();};
    QCOMPARE(value("u"),.035); QCOMPARE(value("v"),.08); QCOMPARE(value("a"),physical.accelerations[3]);
    QVERIFY(near(value("rod N"),.35)); QVERIFY(near(value("cable N"),1.35)); QVERIFY(value("cable N")>=0);
    QCOMPARE(value("V record echo"),.7); QVERIFY(near(value("Q deformation"),.5*.7-.2*.035));
    QVERIFY(near(value("V deformation"),(.15+.2*.035)/.5)); QCOMPARE(value("Q record"),.15);
    QCOMPARE(value("support Rx"),fresh.internalForces[0]-physical.appliedInputs.externalForces[0]);
    const auto proposal=provider.propose(provider.initialize(),physical,0);
    QCOMPARE(proposal.state.lastActions[0].requested,.4); QCOMPARE(proposal.completeInputs.externalForces[3],1.2+.4);
    QCOMPARE(std::get<RStructuralTrussVoltage>(proposal.completeInputs.memberInputs[1]).value,.8);
    QCOMPARE(std::get<RStructuralTrussCharge>(proposal.completeInputs.memberInputs[2]).value,.2);
    QCOMPARE(std::get<RStructuralTrussCableLength>(proposal.completeInputs.memberInputs[3]).value,.99);
    QCOMPARE(rStructuralNativeMeasurementUnit(Measurement::NodeVelocity),QString("m/s"));
    QCOMPARE(rStructuralNativeMeasurementUnit(Measurement::NodeAcceleration),QString("m/s^2"));
    QCOMPARE(rStructuralNativeMeasurementUnit(Measurement::SupportReaction),QString("N"));
    QCOMPARE(rStructuralNativeActionUnit(Action::Charge),QString("C"));
    QCOMPARE(rStructuralNativeActionUnit(Action::CableReferenceLength),QString("m"));
    QVERIFY(rejects([]{rStructuralNativeMeasurementUnit(Measurement(200));}));
    QVERIFY(rejects([]{rStructuralNativeActionUnit(Action(2));}));
    Run run(source); advance(run,source,2); const auto restored=Restart::restore(run.adapter(),
        Restart::build(run.adapter(),run.acceptedState(),2,0,run.acceptedControlState()),run.acceptedControlState());
    QCOMPARE(restored.controlState.fingerprintBytes(),run.acceptedControlState().fingerprintBytes());
} REPORT_NATIVE_EXCEPTIONS

void TestStructuralNativeControl::movingSupportReactionIncludesInertia() try
{
    auto source=fixture();
    source.getProblemSetup().getStructuralInteractionSetup().getPrescribedMotions().push_back(
        {"moving support",{0,RDofType::UX},{{0,0,0,.4},{.096,.0018432,.0384,.4}}});
    auto &spec=source.getProblemSetup().getStructuralNativeControlSpec();
    spec.sensors.push_back(nodeSensor("moving support Rx",Measurement::SupportReaction,909)); bind(source);
    Adapter adapter(source); RStructuralNativeControl provider(adapter); const auto physical=adapter.initialize();
    const auto fresh=RStructuralTrussAssembly::evaluate(adapter.system().referenceNodes,physical.displacements,
        adapter.system().members,physical.appliedInputs.memberInputs);
    const auto sample=provider.sample(physical,0).samples.back();
    QCOMPARE(physical.accelerations[0],.4);
    const double expected=2*.4+fresh.internalForces[0]-physical.appliedInputs.externalForces[0];
    QCOMPARE(sample.value,expected); QVERIFY(sample.value!=fresh.internalForces[0]);
} REPORT_NATIVE_EXCEPTIONS

void TestStructuralNativeControl::identityAndUnsupportedSpec_data()
{
    QTest::addColumn<int>("mutation");
    const char *names[]={"provider-id","revision","code-digest","environment-manifest-and-valid-digest","environment-digest", "N2","sample-phase","seed", "delay","hold","noise", "missing-node","missing-port","unsupported-measurement","free-reaction","sensor-alias","action-alias","port-mode","wrong-initial-V","nonzero-initial-force","wrong-initial-L0","cable-on-rod","held-schedule","target-fields","initial-outside-bounds","negative-rate","missing-loop-sensor","unsupported-action"};
    for(int i=0;i<int(sizeof(names)/sizeof(names[0]));++i)QTest::newRow(names[i])<<i;
}
void TestStructuralNativeControl::identityAndUnsupportedSpec() try
{
    QFETCH(int,mutation); auto source=fixture(true); auto &spec=source.getProblemSetup().getStructuralNativeControlSpec();
    switch(mutation)
    {
        case 0:spec.providerId="another-provider";break;
        case 1:spec.providerRevision="2";break;
        case 2:spec.implementationSha256=QString(64,'a');break;
        case 3:spec.dependencyEnvironmentManifest+="altered\n"; spec.dependencyEnvironmentSha256=QString::fromLatin1(QCryptographicHash::hash(spec.dependencyEnvironmentManifest,QCryptographicHash::Sha256).toHex());break;
        case 4:spec.dependencyEnvironmentSha256=QString(64,'b');break;
        case 5:spec.controlEveryFESteps=2;break;
        case 6:spec.samplePhase=1;break;
        case 7:spec.deterministicSeed=1;break;
        case 8:spec.actions[0].delayTicks=1;break;
        case 9:spec.actions[0].holdTicks=2;break;
        case 10:spec.actions[0].noiseAmplitude=.01;break;
        case 11:spec.sensors[0].physicalNodeId=999;break;
        case 12:spec.sensors[0]=portSensor("displacement",Measurement::PortVoltage,"missing");break;
        case 13:spec.sensors[0].quantity=Measurement(100);break;
        case 14:spec.sensors[0].quantity=Measurement::SupportReaction;break;
        case 15:{auto alias=spec.sensors[0];alias.sensorId="alias";spec.sensors.push_back(alias);break;}
        case 16:{auto alias=spec.actions[0];alias.actionId="alias";spec.actions.push_back(alias);spec.piLoops.push_back({"alias","displacement"});break;}
        case 17:spec.actions[0].quantity=Action::Charge;spec.actions[0].physicalNodeId=0;spec.actions[0].portId="voltage independent";spec.actions[0].initialCommand=.7;break;
        case 18:spec.actions[0].quantity=Action::Voltage;spec.actions[0].physicalNodeId=0;spec.actions[0].portId="voltage independent";spec.actions[0].initialCommand=.8;break;
        case 19:spec.actions[0].initialCommand=.01;break;
        case 20:spec.actions[0].quantity=Action::CableReferenceLength;spec.actions[0].physicalNodeId=0;spec.actions[0].physicalMemberId=609;spec.actions[0].initialCommand=.99;break;
        case 21:spec.actions[0].quantity=Action::CableReferenceLength;spec.actions[0].physicalNodeId=0;spec.actions[0].physicalMemberId=901;spec.actions[0].initialCommand=1;break;
        case 22:source.getProblemSetup().getStructuralTrussDynamicsSetup().heldCommands={{1,{},{},{{101,{{.4,0,0}}}}}};break;
        case 23:spec.sensors[0].physicalMemberId=901;break;
        case 24:spec.actions[0].minimum=.1;break;
        case 25:spec.actions[0].rateLimit=-1;break;
        case 26:spec.piLoops[0].sensorId="unknown";break;
        case 27:spec.actions[0].quantity=Action(2);break;
    }
    QVERIFY(rejects([&]{bind(source);Adapter adapter(source);}));
} REPORT_NATIVE_EXCEPTIONS

void TestStructuralNativeControl::completeStateRestoreTamper_data()
{
    QTest::addColumn<int>("mutation");
    const char *names[]={"after-envelope","before-envelope","left-sample","frame-time","frame-next-time","frame-physical-h","control-tick","FE-count","FE-phase","sample-phase","requested","conditioned","applied","held","last-applied","flags","physical-force","physical-port","source","provider-revision","sensor-quantity","action-quantity","previous-applied","unavailable"};
    for(int i=0;i<int(sizeof(names)/sizeof(names[0]));++i)QTest::newRow(names[i])<<i;
}
void TestStructuralNativeControl::completeStateRestoreTamper() try
{
    QFETCH(int,mutation); const auto source=fixture(true); Run run(source); advance(run,source,3);
    auto control=run.acceptedControlState(); auto record=Restart::build(run.adapter(),run.acceptedState(),3,0,control);
    const auto originalSeal=Restart::stateFingerprint(record,control);
    switch(mutation)
    {
        case 0:control.controllerState[control.controllerState.size()-1]^=char(1);break;
        case 1:control.lastControllerInputState[control.lastControllerInputState.size()-1]^=char(1);break;
        case 2:control.lastFrame.samples[0].value+=.001;break;
        case 3:control.lastFrame.sampleTime+=1e-12;break;
        case 4:control.lastFrame.nextFETime+=1e-12;break;
        case 5:control.lastFrame.physicalDt+=1e-12;break;
        case 6:++control.completedControlTicks;break;
        case 7:++control.completedFE;break;
        case 8:control.fePhase=1;break;
        case 9:control.samplePhase=1;break;
        case 10:control.lastActions[0].requested+=.01;break;
        case 11:control.lastActions[0].conditioned+=.01;break;
        case 12:control.lastActions[0].applied+=.01;break;
        case 13:control.heldValues[0]+=.01;break;
        case 14:control.lastAppliedValues[0]+=.01;break;
        case 15:control.lastActions[0].limitFlags=1;break;
        case 16:record.externalForce[3]+=.01;break;
        case 17:record.memberInputValues[1]+=.01;break;
        case 18:control.sourceFingerprint=QString(64,'a');break;
        case 19:control.providerRevision="2";break;
        case 20:control.lastFrame.samples[0].quantity=Measurement::NodeVelocity;break;
        case 21:control.lastActions[0].quantity=Action::Voltage;break;
        case 22:control.lastActions[0].previousApplied=std::numeric_limits<double>::infinity();break;
        case 23:control.available=false;break;
    }
    QVERIFY(rejects([&]{Restart::restore(run.adapter(),record,control);}));
    // Structural corruption may reject before fingerprint construction; valid
    // edits must change the complete seal even if they retain physical fields.
    if(!rejects([&]{Restart::stateFingerprint(record,control);}))
        QVERIFY(Restart::stateFingerprint(record,control)!=originalSeal);
} REPORT_NATIVE_EXCEPTIONS

void TestStructuralNativeControl::splitRestoreEveryPhase() try
{
    const auto source=fixture(false,12); Run continuous(source); advance(continuous,source,12); const auto expected=seal(continuous);
    for(uint split=0;split<=12;++split)
    {
        Run first(source); advance(first,source,split); auto checkpoint=first.acceptedModel(source);
        checkpoint.getProblemSetup().setRestart(true); Run resumed(checkpoint);
        QCOMPARE(seal(resumed),seal(first)); advance(resumed,checkpoint,12); QCOMPARE(seal(resumed),expected);
    }
} REPORT_NATIVE_EXCEPTIONS

void TestStructuralNativeControl::proposalAndTrialFailureRemainPure() try
{
    const auto source=fixture(true); Adapter adapter(source); RStructuralNativeControl provider(adapter);
    auto physical=adapter.initialize(); const auto control=provider.initialize(); const auto original=control.fingerprintBytes();
    const auto record=Restart::build(adapter,physical,0,0,control); const auto fingerprint=Restart::stateFingerprint(record,control);
    auto bad=physical; bad.accelerations[3]+=1;
    QVERIFY(rejects([&]{provider.propose(control,bad,0);})); QCOMPARE(control.fingerprintBytes(),original);
    bad=physical; bad.appliedInputs.externalForces[3]+=.1;
    QVERIFY(rejects([&]{provider.propose(control,bad,0);})); QCOMPARE(control.fingerprintBytes(),original);
    const auto proposal=provider.propose(control,physical,0); auto inputs=proposal.completeInputs;
    inputs.externalForces[3]=std::numeric_limits<double>::infinity();
    QVERIFY(rejects([&]{adapter.trialWithCompleteInputs(physical,0,inputs);}));
    QCOMPARE(Restart::stateFingerprint(Restart::build(adapter,physical,0,0,control),control),fingerprint);
    inputs=proposal.completeInputs; inputs.memberInputs[3]=RStructuralTrussCableLength{-1};
    QVERIFY(rejects([&]{adapter.trialWithCompleteInputs(physical,0,inputs);}));
    QCOMPARE(control.fingerprintBytes(),original);
    QVERIFY(rejects([&]{provider.completeInputs({std::numeric_limits<double>::quiet_NaN()});}));
} REPORT_NATIVE_EXCEPTIONS

void TestStructuralNativeControl::limiterPhysicalIntervalAndRetry() try
{
    auto source=fixture(false,4,.1); auto &spec=source.getProblemSetup().getStructuralNativeControlSpec();
    spec.actions[0].minimum=-.1; spec.actions[0].maximum=.1; spec.actions[0].rateLimit=.2;
    spec.piLoops[0].proportionalGain=0; spec.piLoops[0].integralGain=0; spec.piLoops[0].bias=2;
    // Large offset exposes a representable physical interval different from nominal Ts.
    source.getTimeSolver().setInputStartTime(10.1); bind(source); Run run(source);
    const auto before=seal(run); auto first=run.prepare(source); const auto candidate=first->candidateModel();
    const auto &control=candidate.getProblemSetup().getStructuralNativeControlState();
    const double h=run.adapter().targetTime(1)-run.adapter().targetTime(0);
    QVERIFY(h!=run.adapter().nominalTimeStep()); QCOMPARE(control.lastFrame.physicalDt,h);
    QCOMPARE(control.lastActions[0].requested,2.0); QCOMPARE(control.lastActions[0].conditioned,2.0);
    QCOMPARE(control.lastActions[0].applied,.2*h); QCOMPARE(control.lastActions[0].limitFlags,quint32(3));
    first.reset(); QCOMPARE(seal(run),before); auto retry=run.prepare(source);
    QCOMPARE(retry->candidateModel().getProblemSetup().getStructuralNativeControlState().fingerprintBytes(),control.fingerprintBytes());
    run.validatePrepared(*retry);
    auto &tamper=retry->candidateModel().getProblemSetup().getStructuralNativeControlState(); tamper.lastFrame.samples[0].value+=1e-4;
    QVERIFY(rejects([&]{run.validatePrepared(*retry);})); QCOMPARE(seal(run),before);
} REPORT_NATIVE_EXCEPTIONS

void TestStructuralNativeControl::defaultFingerprintAndOldRouteGuards() try
{
    auto source=fixture(); Run local(source); const auto control=local.acceptedControlState(); const auto physical=local.acceptedState();
    QVERIFY(rejects([&]{local.adapter().trial(physical,0,0);}));
    QVERIFY(rejects([&]{local.adapter().replayInputs(0);}));
    QVERIFY(rejects([&]{local.adapter().candidateModel(source,physical,0,0);}));
    QVERIFY(rejects([&]{Restart::build(local.adapter(),physical,0,0);}));
    auto record=Restart::build(local.adapter(),physical,0,0,control);
    QVERIFY(rejects([&]{Restart::restore(local.adapter(),record);}));
    source.getProblemSetup().setStructuralNativeControlSpec({}); source.getProblemSetup().setStructuralNativeControlState({}); bind(source);
    Adapter frozen(source); const auto initialized=frozen.initialize(); record=Restart::build(frozen,initialized,0,0);
    QCOMPARE(Restart::stateFingerprint(record),Restart::stateFingerprint(record,RStructuralNativeControlState{}));
    QVERIFY(rejects([&]{frozen.trialWithCompleteInputs(initialized,0,initialized.appliedInputs);}));
    const auto next=frozen.trial(initialized,0,0); QVERIFY(next.state.time>initialized.time);
} REPORT_NATIVE_EXCEPTIONS

QTEST_GUILESS_MAIN(TestStructuralNativeControl)
#include "tst_structural_native_control.moc"

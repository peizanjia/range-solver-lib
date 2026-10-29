#include <QtTest>
#include <QTemporaryDir>
#include <rml_model.h>
#include <rbl_error.h>
#include "rstructural_piezo_element_operators.h"
#include "rstructural_piezo_layered_native_run.h"
#include "rstructural_piezo_layered_restart_codec.h"
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <map>
#include <stdexcept>

namespace {
using Run=RStructuralPiezoLayeredNativeRun;
using Bundle=RStructuralPiezoLayeredRuntimeCallbacks;
using Codec=RStructuralPiezoLayeredRestartCodec;
bool rejects(const std::function<void()> &f){try{f();}catch(...){return true;}return false;}
bool exact(double a,double b){return std::memcmp(&a,&b,sizeof(double))==0;}
RModel scene(bool plate=false,bool configured=false) {
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
    auto &d=m.getProblemSetup().getStructuralDynamicsSetup();d.setEnabled(true);d.setTimeStep(.01);d.setNSteps(8);m.getTimeSolver().setInputStartTime(.01);
    auto &s=m.getProblemSetup().getStructuralPiezoLayeredSetup();s.enabled=true;if(plate)s.family=RStructuralPiezoLayeredModelFamily::FlatPlate;
    const auto native=RStructuralPiezoContinuumMaterial::isotropic(30,.2,7);
    RStructuralPiezoLayeredMaterial h;h.materialId="host";h.elasticityOrCompliance=native.elasticity;h.density=7;h.sourceId="independent-example";h.sourceSha256=QString(64,'a');
    auto p=h;p.materialId="piezo";p.piezoelectric=true;p.coupling[12]=2;p.coupling[13]=2;p.permittivity={{4,0,0,0,4,0,0,0,4}};s.materials={h,p};
    RStructuralPiezoLayeredPort drive;drive.readout.portId="drive";drive.readout.role=RStructuralPiezoPortRole::ActuatorSensor;drive.readout.readout=RStructuralPiezoReadoutMode::ActuatorCharge;drive.initialValue=.12;
    RStructuralPiezoLayeredPort sense;sense.readout.portId="sense";sense.readout.role=RStructuralPiezoPortRole::Sensor;sense.readout.readout=RStructuralPiezoReadoutMode::OpenCircuitVoltage;sense.initialValue=.03;s.ports={drive,sense};
    for(uint i=0;i<(plate?1u:2u);++i){RStructuralPiezoLayeredModelLayer host;host.materialId="host";host.lowerZ=.1;host.upperZ=.3;auto active=host;active.materialId="piezo";active.lowerZ=.3;active.upperZ=.4;active.portId="drive";auto sensor=active;sensor.lowerZ=.4;sensor.upperZ=.45;sensor.portId="sense";s.members.push_back({501+13*i,{host,active,sensor}});}
    using D=RStructuralPiezoPhysicalDof;
    for(uint node:plate?std::vector<uint>{101,134}:std::vector<uint>{101})for(D dof:plate?std::vector<D>{D::UX,D::UY,D::UZ,D::RX,D::RY}:std::vector<D>{D::UX,D::UZ,D::RY})s.supports.push_back({node,dof,0});
    s.backgroundForces={{123,{{.2,0,-.1,0,0,0}}}};s.initialVelocities={{123,{{.04,0,.08,0,0,0}}}};
    s.damping.kind=RStructuralPiezoDampingKind::Rayleigh;s.damping.rayleighAlpha=.06;s.damping.rayleighBeta=.005;
    if(configured){s.excitationConfiguration="ramp/v1";s.controllerConfiguration="counter+observer/v1";s.designConfiguration="native-experiment/v1";}
    s.bindToModel(m);return m;
}
RModel twoCellPlateScene() {
    RModel m=scene(true);
    m.addNode(RNode(2,0,0));m.setNodePhysicalId(4,145);
    m.addNode(RNode(2,1,0));m.setNodePhysicalId(5,156);
    RElement element(R_ELEMENT_QUAD1);const uint connectivity[4]={1,4,5,2};
    for(uint i=0;i<4;++i)element.setNodeId(i,connectivity[i]);
    m.addElement(element,true,0);m.setElementPhysicalId(1,514);
    for(uint i=0;i<m.getNSurfaces();++i) {
        auto &group=m.getSurface(i);
        if(group.hasEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION))
            group.removeEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
        if(group.hasEnvironmentCondition(R_ENVIRONMENT_CONDITION_TEMPERATURE))
            group.removeEnvironmentCondition(R_ENVIRONMENT_CONDITION_TEMPERATURE);
    }
    auto &setup=m.getProblemSetup().getStructuralPiezoLayeredSetup();
    auto second=setup.members.front();second.physicalElementId=514;setup.members.push_back(second);
    setup.bindToModel(m);
    return m;
}
QByteArray opaque(uint proposal,uint observer){return QByteArray::number(proposal)+":"+QByteArray::number(observer);}
std::pair<uint,uint> parse(const QByteArray &state){const auto fields=state.split(':');if(fields.size()!=2)throw std::invalid_argument("bad opaque");bool a=false,b=false;const uint x=fields[0].toUInt(&a),y=fields[1].toUInt(&b);if(!a||!b||opaque(x,y)!=state)throw std::invalid_argument("bad opaque");return {x,y};}
Bundle bundle(const RModel &m,uint stride=3){
    const auto &s=m.getProblemSetup().getStructuralPiezoLayeredSetup();Bundle b;b.sourceFingerprint=s.computeSourceFingerprint(m);b.excitationConfiguration=s.excitationConfiguration;b.controllerConfiguration=s.controllerConfiguration;b.designConfiguration=s.designConfiguration;
    b.backendIdentity="native-test-pure-callbacks/v1";b.forcingIdentity="physical-time-ramp/v1";b.controllerIdentity="counter-observer/v1";b.controlEveryFESteps=stride;b.initialOpaqueControllerState=opaque(0,0);b.envelopes={{-.5,.5,2}};
    const uint nodes=m.getNNodes();b.perturbationForce=[nodes](double t){std::vector<double> f(6*nodes);f[12]=.2*t;f[14]=-.1*t;return f;};
    b.propose=[](const RStructuralPiezoLayeredDynamicsState &physical,uint completed,uint tick,const QByteArray &old){const auto counters=parse(old);if(counters.first!=tick||counters.second!=completed||physical.physicalDisplacements.empty())throw std::invalid_argument("proposal physical/phase mismatch");RStructuralPiezoLayeredControllerProposal p;p.requestedVoltages={tick%2?-2.:2.};p.nextOpaqueState=opaque(tick+1,completed);p.lastFrameBytes=QByteArray::number(completed);return p;};
    b.postFE=[](const RStructuralPiezoLayeredDynamicsState &p,uint completed,uint tick,const std::vector<double> &applied,const QByteArray &old){const auto c=parse(old);if(c.first!=tick||c.second+1!=completed||p.heldActuatorVoltages!=applied)throw std::invalid_argument("observer applied mismatch");return opaque(tick,completed);};
    b.validateOpaque=[](const QByteArray &s,uint completed,uint tick,uint){if(parse(s)!=std::make_pair(tick,completed))throw std::invalid_argument("opaque semantic mismatch");};return b;
}
RModel advance(Run &run,const RModel &source){auto step=run.prepare(source);run.validatePrepared(*step);RModel result(step->candidateModel());run.commit(*step);return result;}
}
#define REPORT_EXCEPTIONS catch(const RError &e){QFAIL(qPrintable(e.getMessage()));}catch(const std::exception &e){QFAIL(e.what());}
class TestStructuralPiezoLayeredNativeRun:public QObject {
    Q_OBJECT
private slots:
    void diskEverySplit_data();
    void diskEverySplit();
    void preparedTamperAndSourceInvalidation();
    void heldConditionerAndRollback();
    void callbackFailureSameInstanceRetry();
    void canonicalCodecAndSemanticRestore();
    void initialAndRestoredVoltageEnvelopeRejection();
    void physicallySealedDifferentNominalStepRejected();
    void continuousBeamElectrodeChargeReadoutAndRestart();
    void continuousPlateElectrodeChargeReadoutAndRestart();
};
void TestStructuralPiezoLayeredNativeRun::diskEverySplit_data(){QTest::addColumn<bool>("plate");QTest::addColumn<QString>("extension");QTest::addColumn<uint>("stride");QTest::addColumn<double>("startTime");for(bool p:{false,true})for(const auto &ext:{QString("ram"),QString("rbm")})for(uint n:{0u,1u,3u})for(double start:{.01,1000.01})QTest::newRow(qPrintable(QString("%1-%2-N%3-start%4").arg(p?"plate":"beam",ext).arg(n).arg(start,0,'g',17)))<<p<<ext<<n<<start;}
void TestStructuralPiezoLayeredNativeRun::diskEverySplit() try {
    QFETCH(bool,plate);QFETCH(QString,extension);QFETCH(uint,stride);QFETCH(double,startTime);auto source=scene(plate,stride!=0);source.getTimeSolver().setInputStartTime(startTime);source.getProblemSetup().getStructuralPiezoLayeredSetup().bindToModel(source);auto cb=stride?bundle(source,stride):Bundle{};Run uninterrupted(source,cb);RModel final(source);
    const auto witness=[&](const Run &run){const auto &s=run.acceptedRuntime();QVERIFY(exact(s.physical.integrator.fixedStep,s.completedFE?.01:0.));QVERIFY(exact(s.nominalControlTs,double(s.controlEveryFESteps)*.01));QVERIFY(exact(s.physical.integrator.state.time,s.completedFE?(startTime-.01)+double(s.completedFE)*.01:startTime-.01));
        const double left=s.completedFE>1?s.initialTime+double(s.completedFE-1)*s.dt:s.initialTime;
        const double h=s.physical.integrator.state.time-left;
        for(const auto &a:s.actions){if(stride){const double amp=std::clamp(a.requested,-.5,.5);const double actual=std::clamp(amp,a.previousApplied-2*h,a.previousApplied+2*h);QVERIFY(exact(a.applied,actual));QVERIFY(a.applied>=-.5 && a.applied<=.5);}}
    };
    witness(uninterrupted);while(uninterrupted.hasRemainingSteps()){final=advance(uninterrupted,source);witness(uninterrupted);}const auto expected=Codec::encode(uninterrupted.acceptedRuntime());QTemporaryDir dir;QVERIFY(dir.isValid());
    for(uint split=0;split<=8;++split){Run first(source,cb);for(uint i=0;i<split;++i)advance(first,source);auto checkpoint=first.acceptedModel(source);const auto path=dir.filePath(QString("split-%1.%2").arg(split).arg(extension));const auto actualPath=checkpoint.nativeWritePaths(path,false).first();checkpoint.write(path,false);RModel disk;disk.read(actualPath);
        QCOMPARE(disk.getTimeSolver().getTimes().size(),size_t(8));QCOMPARE(disk.getProblemSetup().getStructuralDynamicsSetup().getNSteps(),8u);QCOMPARE(disk.getProblemSetup().getStructuralPiezoLayeredRestart().stateBytes,Codec::encode(first.acceptedRuntime()));
        disk.getProblemSetup().setRestart(true);Run restored(disk,cb);QCOMPARE(Codec::encode(restored.acceptedRuntime()),Codec::encode(first.acceptedRuntime()));witness(restored);while(restored.hasRemainingSteps()){advance(restored,disk);witness(restored);}QCOMPARE(Codec::encode(restored.acceptedRuntime()),expected);
        QCOMPARE(restored.acceptedState().electrical.ports.size(),size_t(2));QCOMPARE(final.findVariable(R_VARIABLE_ELECTRIC_POTENTIAL),RConstants::eod);QCOMPARE(final.findVariable(R_VARIABLE_ELECTRIC_CHARGE),RConstants::eod);
    }
} REPORT_EXCEPTIONS
void TestStructuralPiezoLayeredNativeRun::preparedTamperAndSourceInvalidation() try {
    auto source=scene();Run run(source);auto step=run.prepare(source,true);QVERIFY(step->identity().finalStep);QCOMPARE(step->runtime().targetSteps,8u);QCOMPARE(step->candidateModel().getNVariables(),8u);
    const auto check=[&](const std::function<void(RModel &)> &mutate){RModel altered(step->candidateModel());mutate(altered);QVERIFY(rejects([&]{run.validatePrepared(*step,altered);}));};
    check([](RModel &m){m.getVariable(0).setValue(0,0,1);});check([](RModel &m){m.getVariable(0).getVariableData().setScale(2);});check([](RModel &m){m.getTimeSolver().setComputedTime(.02);});check([](RModel &m){m.getTimeSolver().getTimes()[0]=.002;});
    check([](RModel &m){auto &b=m.getProblemSetup().getStructuralPiezoLayeredRestart().stateBytes;b[12]=char(uchar(b[12])^1);});
    check([](RModel &m){m.getProblemSetup().getStructuralPiezoLayeredSetup().materials[0].density+=1;});check([](RModel &m){m.getProblemSetup().getStructuralDynamicsSetup().setTimeStep(.02);});check([](RModel &m){m.getProblemSetup().getStructuralDynamicsSetup().setNSteps(4);});check([](RModel &m){m.getProblemSetup().getStructuralPiezoLayeredSetup().controllerConfiguration="changed";});
    run.validatePrepared(*step);run.commit(*step);QVERIFY(rejects([&]{run.validatePrepared(*step);}));QCOMPARE(run.completedSteps(),1u);
    auto retry=run.prepare(source);run.markNonRetryable();QVERIFY(!run.isRetryable());QVERIFY(rejects([&]{run.validatePrepared(*retry);}));QVERIFY(rejects([&]{run.prepare(source);}));QCOMPARE(run.completedSteps(),1u);
} REPORT_EXCEPTIONS
void TestStructuralPiezoLayeredNativeRun::heldConditionerAndRollback() try {
    auto source=scene(false,true);Run run(source,bundle(source));const auto original=Codec::encode(run.acceptedRuntime());auto abandoned=run.prepare(source);QCOMPARE(Codec::encode(run.acceptedRuntime()),original);auto retry=run.prepare(source);QCOMPARE(Codec::encode(retry->runtime()),Codec::encode(abandoned->runtime()));
    QCOMPARE(retry->runtime().actions[0].requested,2.);QCOMPARE(retry->runtime().actions[0].amplitudeLimited,.5);QVERIFY(retry->runtime().actions[0].amplitudeClipped);QVERIFY(retry->runtime().actions[0].rateClipped);QCOMPARE(retry->runtime().actions[0].previousApplied,.12);QCOMPARE(retry->runtime().lastApplied[0],.12+2*.01);
    run.validatePrepared(*retry);run.commit(*retry);QCOMPARE(run.acceptedRuntime().controlTick,1u);QCOMPARE(run.acceptedRuntime().controlPhase,1u);
    advance(run,source);QCOMPARE(run.acceptedRuntime().heldRequested[0],2.);QCOMPARE(run.acceptedRuntime().controlTick,1u);QCOMPARE(run.acceptedRuntime().controlPhase,2u);advance(run,source);QCOMPARE(run.acceptedRuntime().controlPhase,0u);advance(run,source);QCOMPARE(run.acceptedRuntime().heldRequested[0],-2.);QCOMPARE(run.acceptedRuntime().controlTick,2u);QCOMPARE(run.acceptedRuntime().opaqueControllerState,opaque(2,4));
} REPORT_EXCEPTIONS
void TestStructuralPiezoLayeredNativeRun::callbackFailureSameInstanceRetry() try {
    auto source=scene(false,true);auto cb=bundle(source);bool fail=false;const auto original=cb.propose;cb.propose=[&](const auto &state,uint completed,uint tick,const QByteArray &old){if(fail)throw std::runtime_error("injected pure proposal failure");return original(state,completed,tick,old);};Run run(source,cb);const auto before=Codec::encode(run.acceptedRuntime());fail=true;QVERIFY(rejects([&]{run.prepare(source);}));QCOMPARE(Codec::encode(run.acceptedRuntime()),before);fail=false;auto first=run.prepare(source);auto second=run.prepare(source);QCOMPARE(Codec::encode(first->runtime()),Codec::encode(second->runtime()));run.validatePrepared(*second);run.commit(*second);
    bool observerFail=false;cb=bundle(source);const auto post=cb.postFE;cb.postFE=[&](const auto &state,uint completed,uint tick,const auto &applied,const QByteArray &old){if(observerFail)throw std::runtime_error("injected pure observer failure");return post(state,completed,tick,applied,old);};Run observer(source,cb);observerFail=true;QVERIFY(rejects([&]{observer.prepare(source);}));QCOMPARE(observer.completedSteps(),0u);observerFail=false;advance(observer,source);
    cb=bundle(source);bool forceFail=false;const auto forcing=cb.perturbationForce;cb.perturbationForce=[&](double t){if(forceFail&&t>0)throw std::runtime_error("injected forcing failure");return forcing(t);};Run forced(source,cb);forceFail=true;QVERIFY(rejects([&]{forced.prepare(source);}));QCOMPARE(forced.completedSteps(),0u);forceFail=false;advance(forced,source);
} REPORT_EXCEPTIONS
void TestStructuralPiezoLayeredNativeRun::canonicalCodecAndSemanticRestore() try {
    auto source=scene(false,true);QVERIFY(rejects([&]{Run missing(source);}));auto cb=bundle(source);Run run(source,cb);advance(run,source);auto runtime=run.acceptedRuntime();runtime.heldRequested[0]=-0.;runtime.actions[0].requested=-0.;runtime.actions[0].amplitudeLimited=-0.;runtime.actions[0].applied=runtime.lastApplied[0];
    const auto bytes=Codec::encode(runtime);const auto decoded=Codec::decode(bytes);QVERIFY(std::signbit(decoded.heldRequested[0]));QCOMPARE(Codec::encode(decoded),bytes);QVERIFY(rejects([&]{Codec::decode(bytes.left(bytes.size()-1));}));QVERIFY(rejects([&]{Codec::decode(bytes+"x");}));runtime.dt=std::numeric_limits<double>::infinity();QVERIFY(rejects([&]{Codec::encode(runtime);}));
    auto model=run.acceptedModel(source);model.getProblemSetup().setRestart(true);auto altered=run.acceptedRuntime();altered.opaqueControllerState="invalid";model.getProblemSetup().getStructuralPiezoLayeredRestart().stateBytes=Codec::encode(altered);QVERIFY(rejects([&]{Run bad(model,cb);}));
    model=run.acceptedModel(source);model.getProblemSetup().setRestart(true);auto wrong=bundle(source);wrong.backendIdentity="changed-backend";QVERIFY(rejects([&]{Run bad(model,wrong);}));wrong=bundle(source);wrong.controlEveryFESteps=2;QVERIFY(rejects([&]{Run bad(model,wrong);}));
    altered=run.acceptedRuntime();altered.physical.physicalDisplacements.back()+=1;model.getProblemSetup().getStructuralPiezoLayeredRestart().stateBytes=Codec::encode(altered);QVERIFY(rejects([&]{Run bad(model,cb);}));
    model=run.acceptedModel(source);model.getVariable(0).setValue(0,0,1234);model.getProblemSetup().setRestart(true);Run displayNotAuthority(model,cb);QCOMPARE(Codec::encode(displayNotAuthority.acceptedRuntime()),Codec::encode(run.acceptedRuntime()));
} REPORT_EXCEPTIONS
void TestStructuralPiezoLayeredNativeRun::initialAndRestoredVoltageEnvelopeRejection() try {
    auto source=scene(false,true);auto &setup=source.getProblemSetup().getStructuralPiezoLayeredSetup();setup.ports[0].initialValue=10;setup.bindToModel(source);auto cb=bundle(source);cb.envelopes={{-1,1,2}};
    uint forceCalls=0,proposalCalls=0,observerCalls=0;
    const auto force=cb.perturbationForce;
    const auto proposal=cb.propose;
    const auto observer=cb.postFE;
    cb.perturbationForce=[&](double t){++forceCalls;return force(t);};cb.propose=[&](const auto &p,uint completed,uint tick,const QByteArray &old){++proposalCalls;return proposal(p,completed,tick,old);};cb.postFE=[&](const auto &p,uint completed,uint tick,const auto &applied,const QByteArray &old){++observerCalls;return observer(p,completed,tick,applied,old);};
    const auto sourceBytes=setup.identityBytes();const auto initialOpaque=cb.initialOpaqueControllerState;
    QVERIFY(rejects([&]{Run incompatible(source,cb);}));QCOMPARE(forceCalls,0u);QCOMPARE(proposalCalls,0u);QCOMPARE(observerCalls,0u);QCOMPARE(setup.identityBytes(),sourceBytes);QCOMPARE(cb.initialOpaqueControllerState,initialOpaque);QVERIFY(!source.getProblemSetup().getStructuralPiezoLayeredRestart().available);

    source=scene(false,true);cb=bundle(source);Run valid(source,cb);RStructuralPiezoLayeredDynamics physical(valid.adapter(),RMatrixSolverConf(RMatrixSolverConf::GMRES),0,cb.perturbationForce,cb.forcingIdentity);
    auto bad=valid.acceptedRuntime();bad.physical=physical.initialize({10});physical.restore(bad.physical);bad.heldRequested=bad.lastApplied={10};auto model=valid.acceptedModel(source);model.getProblemSetup().setRestart(true);model.getProblemSetup().getStructuralPiezoLayeredRestart().stateBytes=Codec::encode(bad);
    QVERIFY(rejects([&]{Run outsideApplied(model,cb);}));QCOMPARE(valid.completedSteps(),0u);

    advance(valid,source);bad=valid.acceptedRuntime();bad.physical=physical.trialOnGrid(physical.initialize(),{.5},1,.01);physical.restore(bad.physical);bad.lastApplied={.5};bad.actions[0].previousApplied=std::nextafter(.5,1.);bad.actions[0].applied=.5;bad.actions[0].rateClipped=false;
    model=valid.acceptedModel(source);model.getProblemSetup().setRestart(true);model.getProblemSetup().getStructuralPiezoLayeredRestart().stateBytes=Codec::encode(bad);QVERIFY(rejects([&]{Run outsidePrevious(model,cb);}));QCOMPARE(valid.completedSteps(),1u);
} REPORT_EXCEPTIONS
void TestStructuralPiezoLayeredNativeRun::physicallySealedDifferentNominalStepRejected() try {
    auto source=scene(false,true);source.getTimeSolver().setInputStartTime(1000.01);
    source.getProblemSetup().getStructuralPiezoLayeredSetup().bindToModel(source);
    auto cb=bundle(source);Run valid(source,cb);advance(valid,source);auto altered=valid.acceptedRuntime();
    RStructuralPiezoLayeredDynamics physical(valid.adapter(),RMatrixSolverConf(RMatrixSolverConf::GMRES),valid.acceptedRuntime().initialTime,cb.perturbationForce,cb.forcingIdentity);
    const double differentNominal=std::nextafter(.01,1.);altered.physical=physical.trialOnGrid(physical.initialize(),altered.lastApplied,1,differentNominal);
    physical.restore(altered.physical); // Production API supplied a valid new physical integrity seal.
    QVERIFY(exact(altered.physical.integrator.fixedStep,differentNominal));QVERIFY(!exact(altered.physical.integrator.fixedStep,altered.dt));
    QVERIFY(exact(altered.physical.integrator.state.time,valid.acceptedState().integrator.state.time));
    auto model=valid.acceptedModel(source);model.getProblemSetup().setRestart(true);model.getProblemSetup().getStructuralPiezoLayeredRestart().stateBytes=Codec::encode(altered);
    QVERIFY(rejects([&]{Run changedFixedStep(model,cb);}));QCOMPARE(valid.completedSteps(),1u);
} REPORT_EXCEPTIONS
void TestStructuralPiezoLayeredNativeRun::continuousBeamElectrodeChargeReadoutAndRestart() try {
    auto source=scene();
    auto &setup=source.getProblemSetup().getStructuralPiezoLayeredSetup();
    // Both physical members carry the same actuator and sensor electrode IDs.
    setup.members[1].layers[1].portId="drive";
    setup.bindToModel(source);

    Run run(source);
    QCOMPARE(run.adapter().mesh().memberIds,(std::vector<uint>{501,514}));
    QCOMPARE(run.adapter().mesh().operators.portIds,(std::vector<QString>{"drive","sense"}));
    QCOMPARE(run.acceptedState().electrical.ports.size(),size_t(2));
    QCOMPARE(run.acceptedState().electrical.ports[0].portId,QString("drive"));
    QCOMPARE(run.acceptedState().electrical.ports[0].voltage,.12);
    QCOMPARE(run.acceptedState().readouts.size(),size_t(2));
    QCOMPARE(run.acceptedState().readouts[0].portId,QString("drive"));
    QCOMPARE(run.acceptedState().readouts[0].unit,QString("C"));
    advance(run,source);

    // Independently evaluate each physical member's beam-layer operator at
    // its accepted local q, then sum the two conductor-charge contributions.
    const auto &mesh=run.adapter().mesh();
    const auto &physical=run.acceptedState().physicalDisplacements;
    double summedElementCharge=0;
    std::array<double,2> memberCharges{};
    const auto &sourceSetup=source.getProblemSetup().getStructuralPiezoLayeredSetup();
    const auto responsePort=[&](const QString &id){const auto &ports=run.acceptedState().electrical.ports;const auto it=std::find_if(ports.begin(),ports.end(),[&](const auto &port){return port.portId==id;});if(it==ports.end())throw std::invalid_argument("Missing electrical port in continuous-electrode test");return size_t(it-ports.begin());};
    const size_t driveIndex=responsePort("drive");
    for(size_t memberIndex=0;memberIndex<sourceSetup.members.size();++memberIndex) {
        const auto &member=sourceSetup.members[memberIndex];
        const uint elementPosition=source.findElementPositionByPhysicalId(member.physicalElementId);
        const auto &element=source.getElement(elementPosition);
        const uint node0=source.getNodePhysicalId(element.getNodeId(0));
        const uint node1=source.getNodePhysicalId(element.getNodeId(1));
        const auto nodeIndex=[&](uint id){const auto it=std::find(mesh.nodeIds.begin(),mesh.nodeIds.end(),id);if(it==mesh.nodeIds.end())throw std::invalid_argument("Missing physical node in continuous-electrode test");return size_t(it-mesh.nodeIds.begin());};
        std::array<double,6> localQ{};
        for(size_t end=0;end<2;++end) {
            const size_t ni=nodeIndex(end?node1:node0);
            localQ[3*end]=physical[6*ni];
            localQ[3*end+1]=physical[6*ni+2];
            localQ[3*end+2]=physical[6*ni+4];
        }
        std::map<QString,RStructuralPiezoContinuumMaterial> materials;
        for(const auto &material:sourceSetup.materials) {
            RStructuralPiezoContinuumMaterial native;
            native.elasticity=material.elasticityOrCompliance;native.coupling=material.coupling;
            native.permittivity=material.permittivity;native.density=material.density;native.piezoelectric=material.piezoelectric;
            materials.emplace(material.materialId,native);
        }
        std::vector<RStructuralPiezoLayer> layers;
        for(const auto &layer:member.layers) {
            RStructuralPiezoLayer native;native.material=materials.at(layer.materialId);native.lowerZ=layer.lowerZ;
            native.upperZ=layer.upperZ;native.width=layer.width;native.portId=layer.portId;native.polarity=layer.polarity;
            native.coverage=layer.coverage;layers.push_back(native);
        }
        const double length=std::hypot(source.getNode(element.getNodeId(1)).getX()-source.getNode(element.getNodeId(0)).getX(),
                                       source.getNode(element.getNodeId(1)).getZ()-source.getNode(element.getNodeId(0)).getZ());
        const auto elementOperators=RStructuralPiezoBeamLayerKernel::assemble(length,layers);
        QCOMPARE(elementOperators.portIds,(std::vector<QString>{"drive","sense"}));
        double contribution=0;const size_t nPorts=elementOperators.portIds.size();
        for(size_t p=0;p<nPorts;++p) contribution+=elementOperators.capacitance[driveIndex*nPorts+p]*run.acceptedState().electrical.ports[responsePort(elementOperators.portIds[p])].voltage;
        for(size_t i=0;i<6;++i) contribution-=elementOperators.coupling[i*nPorts+driveIndex]*localQ[i];
        memberCharges[memberIndex]=contribution;summedElementCharge+=contribution;
    }
    QVERIFY(std::abs(memberCharges[0])>1e-12);
    QVERIFY(std::abs(memberCharges[1])>1e-12);
    QVERIFY(std::abs(run.acceptedState().electrical.ports[driveIndex].charge-summedElementCharge)<
            1e-15+1e-7*std::abs(summedElementCharge));
    QCOMPARE(run.acceptedState().readouts[driveIndex].value,run.acceptedState().electrical.ports[driveIndex].charge);

    QTemporaryDir dir;QVERIFY(dir.isValid());
    auto checkpoint=run.acceptedModel(source);const auto path=dir.filePath("continuous-electrode.rbm");
    const auto actualPath=checkpoint.nativeWritePaths(path,false).first();checkpoint.write(path,false);
    RModel disk;disk.read(actualPath);disk.getProblemSetup().setRestart(true);Run restored(disk);
    QCOMPARE(restored.acceptedState().integrityDigest,run.acceptedState().integrityDigest);
    QCOMPARE(restored.acceptedState().electrical.ports.size(),size_t(2));
    QCOMPARE(restored.acceptedState().electrical.ports[0].portId,QString("drive"));
    QCOMPARE(restored.acceptedState().electrical.ports[0].voltage,run.acceptedState().electrical.ports[0].voltage);
    QCOMPARE(restored.acceptedState().electrical.ports[0].charge,run.acceptedState().electrical.ports[0].charge);
    QCOMPARE(restored.acceptedState().readouts[0].portId,QString("drive"));
    QCOMPARE(restored.acceptedState().readouts[0].value,run.acceptedState().readouts[0].value);

    auto broken=source;auto brokenSetup=sourceSetup;brokenSetup.members[1].layers[1].portId="unbound";
    QVERIFY(rejects([&]{brokenSetup.bindToModel(broken);broken.getProblemSetup().getStructuralPiezoLayeredSetup()=brokenSetup;Run invalid(broken);}));
} REPORT_EXCEPTIONS
void TestStructuralPiezoLayeredNativeRun::continuousPlateElectrodeChargeReadoutAndRestart() try {
    // Reuse MITC4 as-is; this qualifies two-cell aggregation, native readout, and restart only.
    auto source=twoCellPlateScene();
    const auto &setup=source.getProblemSetup().getStructuralPiezoLayeredSetup();
    QCOMPARE(source.getNElements(R_ELEMENT_QUAD1),2u);
    QCOMPARE(setup.members.size(),size_t(2));
    QCOMPARE(setup.members[0].physicalElementId,501u);
    QCOMPARE(setup.members[1].physicalElementId,514u);
    QCOMPARE(setup.ports.size(),size_t(2));
    QCOMPARE(setup.ports[0].readout.portId,QString("drive"));
    QCOMPARE(setup.ports[0].readout.readout,RStructuralPiezoReadoutMode::ActuatorCharge);
    QCOMPARE(setup.ports[0].initialValue,.12);
    QCOMPARE(setup.ports[1].readout.portId,QString("sense"));
    QCOMPARE(setup.ports[1].readout.readout,RStructuralPiezoReadoutMode::OpenCircuitVoltage);
    QCOMPARE(setup.ports[1].initialValue,.03);
    QVERIFY(std::any_of(setup.supports.begin(),setup.supports.end(),[](const auto &support){return support.physicalNodeId==101;}));
    QVERIFY(std::any_of(setup.supports.begin(),setup.supports.end(),[](const auto &support){return support.physicalNodeId==134;}));
    const auto force=std::find_if(setup.backgroundForces.begin(),setup.backgroundForces.end(),[](const auto &node){return node.physicalNodeId==123;});
    QVERIFY(force!=setup.backgroundForces.end());
    QCOMPARE(force->values[0],.2);
    QCOMPARE(force->values[2],-.1);
    const auto velocity=std::find_if(setup.initialVelocities.begin(),setup.initialVelocities.end(),[](const auto &node){return node.physicalNodeId==123;});
    QVERIFY(velocity!=setup.initialVelocities.end());
    QCOMPARE(velocity->values[0],.04);
    QCOMPARE(velocity->values[2],.08);
    for(const auto &member:setup.members) {
        QVERIFY(std::any_of(member.layers.begin(),member.layers.end(),[](const auto &layer){return layer.portId=="drive";}));
        QVERIFY(std::any_of(member.layers.begin(),member.layers.end(),[](const auto &layer){return layer.portId=="sense";}));
    }

    Run run(source);
    const auto &mesh=run.adapter().mesh();
    QCOMPARE(mesh.memberIds,(std::vector<uint>{501,514}));
    QCOMPARE(mesh.nodeIds,(std::vector<uint>{101,112,123,134,145,156}));
    QCOMPARE(mesh.operators.portIds,(std::vector<QString>{"drive","sense"}));
    QCOMPARE(run.acceptedState().electrical.ports.size(),size_t(2));
    QCOMPARE(run.acceptedState().readouts.size(),size_t(2));
    advance(run,source);

    const auto &accepted=run.acceptedState();
    const auto responsePort=[&](const QString &id) {
        const auto &ports=accepted.electrical.ports;
        const auto it=std::find_if(ports.begin(),ports.end(),[&](const auto &port){return port.portId==id;});
        if(it==ports.end())throw std::invalid_argument("Missing electrical port in continuous-plate test");
        return size_t(it-ports.begin());
    };
    const auto responseReadout=[&](const QString &id) {
        const auto &readouts=accepted.readouts;
        const auto it=std::find_if(readouts.begin(),readouts.end(),[&](const auto &readout){return readout.portId==id;});
        if(it==readouts.end())throw std::invalid_argument("Missing readout in continuous-plate test");
        return size_t(it-readouts.begin());
    };
    const size_t driveIndex=responsePort("drive"),senseIndex=responsePort("sense");
    const size_t driveReadoutIndex=responseReadout("drive"),senseReadoutIndex=responseReadout("sense");
    QCOMPARE(accepted.electrical.ports[driveIndex].portId,QString("drive"));
    QCOMPARE(accepted.electrical.ports[driveIndex].voltage,.12);
    QCOMPARE(accepted.electrical.ports[senseIndex].portId,QString("sense"));
    QCOMPARE(accepted.electrical.ports[senseIndex].charge,.03);
    QCOMPARE(accepted.readouts[driveReadoutIndex].mode,RStructuralPiezoReadoutMode::ActuatorCharge);
    QCOMPARE(accepted.readouts[driveReadoutIndex].unit,QString("C"));
    QCOMPARE(accepted.readouts[senseReadoutIndex].mode,RStructuralPiezoReadoutMode::OpenCircuitVoltage);
    QCOMPARE(accepted.readouts[senseReadoutIndex].unit,QString("V"));

    const auto &physical=accepted.physicalDisplacements;
    const auto &sourceSetup=source.getProblemSetup().getStructuralPiezoLayeredSetup();
    std::map<QString,RStructuralPiezoContinuumMaterial> materials;
    for(const auto &material:sourceSetup.materials) {
        RStructuralPiezoContinuumMaterial native;
        native.elasticity=material.elasticityOrCompliance;native.coupling=material.coupling;
        native.permittivity=material.permittivity;native.density=material.density;native.piezoelectric=material.piezoelectric;
        materials.emplace(material.materialId,native);
    }
    const auto meshNodeIndex=[&](uint physicalId) {
        const auto it=std::find(mesh.nodeIds.begin(),mesh.nodeIds.end(),physicalId);
        if(it==mesh.nodeIds.end())throw std::invalid_argument("Missing physical node in continuous-plate test");
        return size_t(it-mesh.nodeIds.begin());
    };
    const auto portIndex=[](const std::vector<QString> &ids,const QString &id) {
        const auto it=std::find(ids.begin(),ids.end(),id);
        if(it==ids.end())throw std::invalid_argument("Element kernel omitted a declared continuous-plate port");
        return size_t(it-ids.begin());
    };
    double summedDriveCharge=0.0,driveChargeScale=0.0;
    double summedSensorCoupling=0.0,summedSensorCouplingScale=0.0,summedSensorCapacitance=0.0;
    std::array<double,2> cellDriveCoupling{},cellSensorCoupling{},cellDriveCharge{};
    for(size_t memberIndex=0;memberIndex<sourceSetup.members.size();++memberIndex) {
        const auto &member=sourceSetup.members[memberIndex];
        const uint elementPosition=source.findElementPositionByPhysicalId(member.physicalElementId);
        const auto &element=source.getElement(elementPosition);
        std::array<std::array<double,2>,4> xy{};
        std::array<double,20> localQ{};
        for(size_t node=0;node<4;++node) {
            const uint modelNode=element.getNodeId(uint(node));
            xy[node]={{source.getNode(modelNode).getX(),source.getNode(modelNode).getY()}};
            const size_t globalNode=meshNodeIndex(source.getNodePhysicalId(modelNode));
            for(size_t dof=0;dof<5;++dof)localQ[5*node+dof]=physical[6*globalNode+dof];
        }
        std::vector<RStructuralPiezoLayer> layers;
        for(const auto &layer:member.layers) {
            RStructuralPiezoLayer native;native.material=materials.at(layer.materialId);native.lowerZ=layer.lowerZ;
            native.upperZ=layer.upperZ;native.width=layer.width;native.portId=layer.portId;native.polarity=layer.polarity;
            native.coverage=layer.coverage;layers.push_back(native);
        }
        const auto operators=RStructuralPiezoPlateMitc4Kernel::assemble(xy,layers);
        QCOMPARE(operators.dofs,20u);
        QCOMPARE(operators.portIds,(std::vector<QString>{"drive","sense"}));
        const size_t nPorts=operators.portIds.size();
        const size_t localDrive=portIndex(operators.portIds,"drive");
        const size_t localSense=portIndex(operators.portIds,"sense");
        double capacitiveDriveCharge=0.0,capacitiveDriveScale=0.0;
        double driveCoupling=0.0,driveCouplingScale=0.0;
        double sensorCoupling=0.0,sensorCouplingScale=0.0;
        for(size_t p=0;p<nPorts;++p) {
            const double term=operators.capacitance[localDrive*nPorts+p]
                *accepted.electrical.ports[responsePort(operators.portIds[p])].voltage;
            capacitiveDriveCharge+=term;capacitiveDriveScale+=std::abs(term);
        }
        for(size_t i=0;i<localQ.size();++i) {
            const double driveTerm=operators.coupling[i*nPorts+localDrive]*localQ[i];
            const double sensorTerm=operators.coupling[i*nPorts+localSense]*localQ[i];
            driveCoupling+=driveTerm;driveCouplingScale+=std::abs(driveTerm);
            sensorCoupling+=sensorTerm;sensorCouplingScale+=std::abs(sensorTerm);
        }
        cellDriveCoupling[memberIndex]=driveCoupling;
        cellSensorCoupling[memberIndex]=sensorCoupling;
        cellDriveCharge[memberIndex]=capacitiveDriveCharge-driveCoupling;
        summedDriveCharge+=cellDriveCharge[memberIndex];
        driveChargeScale+=capacitiveDriveScale+driveCouplingScale;
        summedSensorCoupling+=sensorCoupling;
        summedSensorCouplingScale+=sensorCouplingScale;
        const double sensorCapacitance=operators.capacitance[localSense*nPorts+localSense];
        QVERIFY(std::isfinite(sensorCapacitance));
        QVERIFY(sensorCapacitance>0.0);
        summedSensorCapacitance+=sensorCapacitance;
        QVERIFY(driveCouplingScale>0.0);
        QVERIFY(std::abs(driveCoupling)>1e-10*driveCouplingScale);
        QVERIFY(sensorCouplingScale>0.0);
        QVERIFY(std::abs(sensorCoupling)>1e-10*sensorCouplingScale);
        QVERIFY(std::abs(cellDriveCharge[memberIndex])>1e-10*(capacitiveDriveScale+driveCouplingScale));
    }
    QVERIFY(std::abs(cellDriveCoupling[0])>0.0);
    QVERIFY(std::abs(cellDriveCoupling[1])>0.0);
    QVERIFY(std::abs(cellSensorCoupling[0])>0.0);
    QVERIFY(std::abs(cellSensorCoupling[1])>0.0);
    QVERIFY(driveChargeScale>0.0);
    QVERIFY(std::abs(accepted.electrical.ports[driveIndex].charge-summedDriveCharge)<=1e-9*driveChargeScale);
    QCOMPARE(accepted.readouts[driveReadoutIndex].value,accepted.electrical.ports[driveIndex].charge);
    QVERIFY(summedSensorCapacitance>0.0);
    QVERIFY(summedSensorCouplingScale>0.0);
    const double fixedSensorCharge=setup.ports[1].initialValue;
    const double expectedSensorVoltage=(fixedSensorCharge+summedSensorCoupling)/summedSensorCapacitance;
    const double sensorVoltageScale=(std::abs(fixedSensorCharge)+summedSensorCouplingScale)/summedSensorCapacitance;
    QVERIFY(std::isfinite(expectedSensorVoltage));
    QVERIFY(sensorVoltageScale>0.0);
    QVERIFY(std::abs(accepted.electrical.ports[senseIndex].voltage-expectedSensorVoltage)<=1e-9*sensorVoltageScale);
    QCOMPARE(accepted.readouts[senseReadoutIndex].value,accepted.electrical.ports[senseIndex].voltage);

    QTemporaryDir dir;QVERIFY(dir.isValid());
    auto checkpoint=run.acceptedModel(source);const auto path=dir.filePath("continuous-plate-electrode.rbm");
    const auto actualPath=checkpoint.nativeWritePaths(path,false).first();checkpoint.write(path,false);
    RModel disk;disk.read(actualPath);disk.getProblemSetup().setRestart(true);Run restored(disk);
    const auto compareElectrical=[&](const Run &expected,const Run &actual) {
        const auto &expectedState=expected.acceptedState();const auto &actualState=actual.acceptedState();
        QCOMPARE(actualState.electrical.ports.size(),expectedState.electrical.ports.size());
        QCOMPARE(actualState.readouts.size(),expectedState.readouts.size());
        for(size_t i=0;i<expectedState.electrical.ports.size();++i) {
            QCOMPARE(actualState.electrical.ports[i].portId,expectedState.electrical.ports[i].portId);
            QCOMPARE(actualState.electrical.ports[i].voltage,expectedState.electrical.ports[i].voltage);
            QCOMPARE(actualState.electrical.ports[i].charge,expectedState.electrical.ports[i].charge);
        }
        for(size_t i=0;i<expectedState.readouts.size();++i) {
            QCOMPARE(actualState.readouts[i].portId,expectedState.readouts[i].portId);
            QCOMPARE(actualState.readouts[i].mode,expectedState.readouts[i].mode);
            QCOMPARE(actualState.readouts[i].unit,expectedState.readouts[i].unit);
            QCOMPARE(actualState.readouts[i].value,expectedState.readouts[i].value);
        }
        QCOMPARE(actualState.integrityDigest,expectedState.integrityDigest);
        QCOMPARE(Codec::encode(actual.acceptedRuntime()),Codec::encode(expected.acceptedRuntime()));
    };
    compareElectrical(run,restored);
    advance(run,source);advance(restored,source);
    compareElectrical(run,restored);

    auto broken=twoCellPlateScene();
    auto brokenSetup=broken.getProblemSetup().getStructuralPiezoLayeredSetup();
    brokenSetup.members[1].layers[1].portId="undeclared";
    QVERIFY(rejects([&]{brokenSetup.bindToModel(broken);}));
} REPORT_EXCEPTIONS
QTEST_APPLESS_MAIN(TestStructuralPiezoLayeredNativeRun)
#include "tst_structural_piezo_layered_native_run.moc"

#include "rstructural_piezo_axial_native_run.h"
#include "rstructural_piezo_axial_restart_codec.h"
#include <rml_model.h>
#include "rstructural_static_truss_model.h"
#include <QJsonDocument>
#include <QJsonObject>
#include <QCryptographicHash>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <stdexcept>

namespace {
void require(bool b,const char *message) {if(!b)throw std::invalid_argument(message);}
// The legacy adapter compiles a source definition, not axial restart authority.
// Its old restart gate is bypassed only in an independent source copy; every
// original physical, initial, declaration and clock field remains unchanged.
RModel sourceDefinition(const RModel &source) {
    RModel copy(source);
    if(source.getProblemSetup().getRestart()){
        const auto &r=source.getProblemSetup().getStructuralPiezoAxialRestart();r.validate();
        require(source.getProblemSetup().getStructuralPiezoAxialSetup().enabled && r.available,"Axial source extraction requires its complete restart authority");
        copy.getProblemSetup().setRestart(false);
    }
    return copy;
}
bool bits(double a,double b) {return std::memcmp(&a,&b,8)==0;}
bool sameVector(const std::vector<double> &a,const std::vector<double> &b) {return a.size()==b.size() && (a.empty() || std::memcmp(a.data(),b.data(),8*a.size())==0);}
void integer(QByteArray &b,quint64 n) {for(uint i=0;i<8;++i)b.append(char((n>>(8*i))&255));}
void blob(QByteArray &b,const QByteArray &v) {integer(b,quint64(v.size()));b.append(v);}
void real(QByteArray &b,double v) {require(std::isfinite(v),"Nonfinite axial runtime configuration");quint64 n;std::memcpy(&n,&v,8);integer(b,n);}
bool resultsEqual(const RResults &a,const RResults &b) {
    if(a.getNNodes()!=b.getNNodes() || a.getNElements()!=b.getNElements() || a.getNVariables()!=b.getNVariables())return false;
    for(uint i=0;i<a.getNVariables();++i){const auto &x=a.getVariable(i),&y=b.getVariable(i);
        if(x.getType()!=y.getType() || x.getApplyType()!=y.getApplyType() || x.getName()!=y.getName() || x.getUnits()!=y.getUnits() || x.getNVectors()!=y.getNVectors() || x.getNValues()!=y.getNValues())return false;
        const auto &xd=x.getVariableData(),&yd=y.getVariableData();
        if(!bits(xd.getMinDisplayValue(),yd.getMinDisplayValue()) || !bits(xd.getMaxDisplayValue(),yd.getMaxDisplayValue()) || xd.getDisplayValueRange()!=yd.getDisplayValueRange() || xd.getValueRangeName()!=yd.getValueRangeName() || !bits(xd.getScale(),yd.getScale()))return false;
        for(uint c=0;c<x.getNVectors();++c)for(uint n=0;n<x.getNValues();++n)if(!bits(x.getValue(c,n),y.getValue(c,n)))return false;
    }return true;
}
void publish(RModel &m,RVariableType type,const std::vector<double> &values,const QString &name) {
    require(values.size()==3*size_t(m.getNNodes()),"Axial result vector dimensions mismatch");
    auto &v=m.getVariable(m.addVariable(type));v.setName(name);v.setApplyType(R_VARIABLE_APPLY_NODE);v.resize(3,m.getNNodes());double maximum=0;
    for(uint n=0;n<m.getNNodes();++n){for(uint c=0;c<3;++c)v.setValue(c,n,values[3*n+c]);maximum=std::max(maximum,std::hypot(values[3*n],std::hypot(values[3*n+1],values[3*n+2])));}
    require(std::isfinite(maximum),"Axial display magnitude overflow");v.getVariableData().setMinMaxDisplayValue(0,maximum);
}

}
RStructuralPiezoAxialNativeRun::PreparedStep::PreparedStep()=default;
RStructuralPiezoAxialNativeRun::PreparedStep::~PreparedStep()=default;
RModel &RStructuralPiezoAxialNativeRun::PreparedStep::candidateModel() noexcept{return *model;}
const RModel &RStructuralPiezoAxialNativeRun::PreparedStep::candidateModel() const noexcept{return *model;}
const RAcceptedStepIdentity &RStructuralPiezoAxialNativeRun::PreparedStep::identity() const noexcept{return stepIdentity;}
const RStructuralTrussDynamicsState &RStructuralPiezoAxialNativeRun::PreparedStep::result() const noexcept{return sealedRuntime->physical;}
const RStructuralPiezoAxialRunState &RStructuralPiezoAxialNativeRun::PreparedStep::runtime() const noexcept{return *sealedRuntime;}
uint RStructuralPiezoAxialNativeRun::PreparedStep::completedSteps() const noexcept{return uint(stepIdentity.sequence);}

RStructuralPiezoAxialNativeRun::RStructuralPiezoAxialNativeRun(const RModel &source,RStructuralPiezoAxialRuntimeCallbacks bundle)
    :callbacks(std::move(bundle)),frozenAdapter(sourceDefinition(source)),actualSystem(frozenAdapter.system()) {
    const auto &setup=source.getProblemSetup().getStructuralPiezoAxialSetup();setup.validate(source);
    const auto &d=source.getProblemSetup().getStructuralDynamicsSetup();d.validate();
    require(setup.enabled && d.getEnabled() && d.getMethod()==RStructuralTimeIntegrationMethod::Newmark && d.getNewmarkGamma()==.5 && d.getNewmarkBeta()==.25,"Axial facade requires enabled average-acceleration Newmark source");
    dt=d.getTimeStep();target=d.getNSteps();initialTime=source.getTimeSolver().getInputStartTime()-dt;
    require(std::isfinite(initialTime) && target>0 && target<=1000000,"Invalid axial source clock");
    frozenSource=setup.computeSourceFingerprint(source);
    const bool custom=bool(callbacks.prepareFE)||bool(callbacks.bindForceState)||bool(callbacks.intervalPerturbationForce)||bool(callbacks.perturbationForce)||bool(callbacks.propose)||bool(callbacks.postFE)||bool(callbacks.validateOpaque)||bool(callbacks.validateAcceptedBoundary)||!callbacks.authoritativeFullC.empty()||!callbacks.initialActuatorVoltages.empty()||!callbacks.envelopes.empty()||!callbacks.backendIdentity.isEmpty()||!callbacks.sourceFingerprint.isEmpty()||!callbacks.recipeConfiguration.isEmpty();
    if(custom)require(callbacks.sourceFingerprint==frozenSource && !callbacks.backendIdentity.isEmpty() && callbacks.recipeConfiguration==setup.recipeConfiguration,"Axial runtime bundle must bind actual source, raw recipe and backend identity");
    else {
        const auto recipe=QJsonDocument::fromJson(setup.recipeConfiguration).object();
        for(const auto &key:{"controller","excitation","design","damping"})require(recipe.value(key).toObject().isEmpty(),"Configured axial recipe requires an application-certified callbacks bundle");
    }
    const auto declared=QJsonDocument::fromJson(setup.recipeConfiguration).object();
    require(declared.value("controller").toObject().isEmpty() || bool(callbacks.propose),"Declared axial controller has no proposal callback");
    require(declared.value("excitation").toObject().isEmpty() || bool(callbacks.perturbationForce),"Declared axial excitation has no time-only force callback");
    require(callbacks.controlEveryFESteps>0 && callbacks.controlEveryFESteps<=1000000,"Invalid axial controller FE stride");
    if(callbacks.propose)require(!callbacks.controllerIdentity.isEmpty() && bool(callbacks.validateAcceptedBoundary),"Axial controller requires semantic accepted-boundary validator and identity");
    else require(!callbacks.postFE && !callbacks.validateOpaque && !callbacks.validateAcceptedBoundary && callbacks.initialOpaqueControllerState.isEmpty() && callbacks.controllerIdentity.isEmpty(),"Controller runtime without proposal is unsupported");
    if(!callbacks.perturbationForce){const size_t n=actualSystem.prescribedMask.size();callbacks.perturbationForce=[n](double){return std::vector<double>(n,0);};callbacks.forcingIdentity="range-fea/axial-zero-time-force/v1";}
    require(!callbacks.forcingIdentity.isEmpty(),"Axial time callback needs declared identity");
    if(callbacks.backendIdentity.isEmpty())callbacks.backendIdentity="range-fea/native-axial-Newmark/no-control/v1";
    // Match the actual operator package's electrical input column order.
    for(const auto &id:frozenAdapter.maps().orderedPortIds){
        const auto &ports=frozenAdapter.piezoReadoutSetup().ports;
        const auto port=std::find_if(ports.begin(),ports.end(),[&](const auto &p){return p.portId==id;});
        require(port!=ports.end(),"Axial readout catalog is incomplete");if(port->role==RStructuralPiezoPortRole::Sensor)continue;
        const auto member=frozenAdapter.maps().memberByPortId.at(id);require(std::holds_alternative<RStructuralTrussVoltage>(frozenAdapter.initialCompleteInputs().memberInputs.at(member)),"Axial actuators require voltage boundaries");actuatorIds.push_back(id);actuatorMembers.push_back(member);
    }
    require(callbacks.envelopes.empty() || callbacks.envelopes.size()==actuatorIds.size(),"Axial envelope dimensions differ from actuator order");
    require(callbacks.initialActuatorVoltages.empty() || callbacks.initialActuatorVoltages.size()==actuatorIds.size(),"Axial initial voltage dimensions differ from actuator order");
    std::vector<double> initialVa;for(auto i:actuatorMembers)initialVa.push_back(std::get<RStructuralTrussVoltage>(frozenAdapter.initialCompleteInputs().memberInputs[i]).value);
    if(!callbacks.initialActuatorVoltages.empty())initialVa=callbacks.initialActuatorVoltages;
    actualSystem.physicalDamping=callbacks.authoritativeFullC;
    driver=std::make_unique<RStructuralTrussDynamics>(actualSystem);
    readout=std::make_unique<RStructuralPiezoNativeReadout>(frozenAdapter,frozenAdapter.piezoReadoutSetup(),actualSystem);
    QByteArray config;blob(config,frozenSource.toUtf8());blob(config,setup.identityBytes());blob(config,callbacks.backendIdentity.toUtf8());blob(config,callbacks.forcingIdentity.toUtf8());blob(config,callbacks.controllerIdentity.toUtf8());integer(config,callbacks.controlEveryFESteps);
    blob(config,callbacks.initialOpaqueControllerState);integer(config,bool(callbacks.postFE));integer(config,bool(callbacks.propose));integer(config,bool(callbacks.validateOpaque));integer(config,bool(callbacks.validateAcceptedBoundary));integer(config,callbacks.authoritativeFullC.size());for(double x:callbacks.authoritativeFullC)real(config,x);
    integer(config,initialVa.size());for(double x:initialVa)real(config,x);
    integer(config,callbacks.envelopes.size());for(const auto &e:callbacks.envelopes){require(e.minimum<=e.maximum && e.maximumRate>=0,"Invalid axial voltage envelope");real(config,e.minimum);real(config,e.maximum);real(config,e.maximumRate);}
    for(size_t i=0;i<callbacks.envelopes.size();++i)require(initialVa[i]>=callbacks.envelopes[i].minimum && initialVa[i]<=callbacks.envelopes[i].maximum,"Axial initial actuator voltage is outside its amplitude envelope");
    real(config,initialTime);real(config,dt);integer(config,target);blob(config,driver->systemIdentity());
    require(bool(callbacks.prepareFE)==bool(callbacks.bindForceState) && bool(callbacks.prepareFE)==bool(callbacks.intervalPerturbationForce),"Incomplete transactional force callbacks");
    require(!callbacks.prepareFE || (callbacks.bindForceState && callbacks.intervalPerturbationForce && callbacks.propose && callbacks.postFE),"Sampled forcing requires a complete transactional runtime bundle");
    if(callbacks.prepareFE){blob(config,"transactional-interval-force/v1");integer(config,1);}
    frozenConfiguration=QCryptographicHash::hash(config,QCryptographicHash::Sha256);
    auto state=std::make_unique<RStructuralPiezoAxialRunState>();
    if(source.getProblemSetup().getRestart()){
        const auto &r=source.getProblemSetup().getStructuralPiezoAxialRestart();r.validate();require(r.available && r.sourceFingerprint==frozenSource,"Axial restart needs matching complete runtime authority");*state=RStructuralPiezoAxialRestartCodec::decode(r.runtimeState);
    }else{
        state->opaqueControllerState=callbacks.initialOpaqueControllerState;
        if(callbacks.prepareFE)state->opaqueControllerState=callbacks.prepareFE(initialTime,0,state->opaqueControllerState);
    }
    if(callbacks.bindForceState)callbacks.bindForceState(state->opaqueControllerState);
    const auto seed=frozenAdapter.sourceInitialState();initialPhysical=driver->initialize(initialTime,seed.displacements,seed.velocities,inputs(initialTime,initialVa));
    if(!source.getProblemSetup().getRestart()){
        state->physical=state->initialPhysicalReference=initialPhysical;state->readouts=readout->observe(initialPhysical);state->sourceFingerprint=frozenSource;state->forcingIdentity=callbacks.forcingIdentity;state->controllerIdentity=callbacks.controllerIdentity;state->backendIdentity=callbacks.backendIdentity;state->configurationIdentity=frozenConfiguration;
        state->targetSteps=target;state->controlEveryFESteps=callbacks.controlEveryFESteps;state->initialTime=initialTime;state->dt=dt;state->nominalControlTs=double(callbacks.controlEveryFESteps)*dt;state->envelopes=callbacks.envelopes;state->heldRequested=state->lastApplied=initialVa;
    }
    certifyRuntime(*state);accepted=std::move(state);
}
RStructuralPiezoAxialNativeRun::~RStructuralPiezoAxialNativeRun()=default;
const RStructuralTrussDynamicsModelAdapter &RStructuralPiezoAxialNativeRun::adapter() const noexcept{return frozenAdapter;}
const RStructuralTrussDynamicsSystem &RStructuralPiezoAxialNativeRun::system() const noexcept{return actualSystem;}
const RStructuralTrussDynamicsState &RStructuralPiezoAxialNativeRun::acceptedState() const noexcept{return accepted->physical;}
const RStructuralPiezoAxialRunState &RStructuralPiezoAxialNativeRun::acceptedRuntime() const noexcept{return *accepted;}
uint RStructuralPiezoAxialNativeRun::completedSteps() const noexcept{return accepted->completedFE;}
uint RStructuralPiezoAxialNativeRun::targetSteps() const noexcept{return target;}
bool RStructuralPiezoAxialNativeRun::hasRemainingSteps() const noexcept{return accepted->completedFE<target;}
bool RStructuralPiezoAxialNativeRun::isRetryable() const noexcept{return retryable;}
void RStructuralPiezoAxialNativeRun::markNonRetryable() noexcept{retryable=false;}
double RStructuralPiezoAxialNativeRun::targetTime(uint step) const {const double t=step?initialTime+double(step)*dt:initialTime;require(std::isfinite(t),"Axial canonical clock overflow");return t;}
void RStructuralPiezoAxialNativeRun::validateCurrentSource(const RModel &source) const {
    const auto &s=source.getProblemSetup().getStructuralPiezoAxialSetup();s.validate(source);require(s.computeSourceFingerprint(source)==frozenSource,"Axial native source changed since freeze");
}
RStructuralPiezoAxialRunState RStructuralPiezoAxialNativeRun::certifyAcceptedRuntime(const RModel &model) const {
    validateCurrentSource(model);
    const auto &restart=model.getProblemSetup().getStructuralPiezoAxialRestart();restart.validate();
    require(restart.available && restart.sourceFingerprint==frozenSource,"Axial accepted MODEL requires complete matching restart authority");
    auto state=RStructuralPiezoAxialRestartCodec::decode(restart.runtimeState);certifyRuntime(state);
    const auto &time=model.getTimeSolver();
    require(time.getEnabled() && bits(time.getComputedTime(),state.physical.time)
        && time.getCurrentTimeStep()==(state.completedFE?state.completedFE-1:0)
        && time.getInputNTimeSteps()==target && bits(time.getInputTimeStepSize(),dt)
        && time.getTimes().size()==target
        && time.getOutputFrequency()==model.getProblemSetup().getStructuralDynamicsSetup().getOutputFrequency(),"Axial accepted MODEL clock changed");
    for(uint i=0;i<target;++i)require(bits(time.getTimes()[i],targetTime(i+1)),"Axial original accepted target grid changed");
    return state;
}
void RStructuralPiezoAxialNativeRun::certifyAcceptedRuntime(const RStructuralPiezoAxialRunState &state) const {
    certifyRuntime(state);
}
std::vector<double> RStructuralPiezoAxialNativeRun::totalForce(double time) const {
    auto f=callbacks.perturbationForce(time);const auto &background=frozenAdapter.initialCompleteInputs().externalForces;
    require(f.size()==background.size(),"Axial time force dimensions mismatch");for(size_t i=0;i<f.size();++i){require(std::isfinite(f[i]),"Nonfinite axial time force");f[i]+=background[i];require(std::isfinite(f[i]),"Axial total force overflow");}return f;
}
RStructuralTrussDynamicsInputs RStructuralPiezoAxialNativeRun::inputs(double time,const std::vector<double> &voltage) const {
    require(voltage.size()==actuatorMembers.size(),"Axial voltage dimensions mismatch");auto in=frozenAdapter.initialCompleteInputs();in.externalForces=totalForce(time);in.rightExternalForces.clear();
    for(size_t i=0;i<voltage.size();++i){require(std::isfinite(voltage[i]),"Nonfinite axial actuator voltage");in.memberInputs[actuatorMembers[i]]=RStructuralTrussVoltage{voltage[i]};}return in;
}
void RStructuralPiezoAxialNativeRun::certifyRuntime(const RStructuralPiezoAxialRunState &s) const {
    if(callbacks.bindForceState)callbacks.bindForceState(s.opaqueControllerState);
    require(s.sourceFingerprint==frozenSource && s.configurationIdentity==frozenConfiguration && s.forcingIdentity==callbacks.forcingIdentity && s.controllerIdentity==callbacks.controllerIdentity && s.backendIdentity==callbacks.backendIdentity,"Axial runtime source/backend/configuration mismatch");
    require(s.targetSteps==target && s.completedFE<=target && s.controlEveryFESteps==callbacks.controlEveryFESteps && s.controlPhase==s.completedFE%s.controlEveryFESteps && bits(s.initialTime,initialTime) && bits(s.dt,dt) && bits(s.nominalControlTs,double(s.controlEveryFESteps)*dt),"Axial runtime clock/control phase mismatch");
    require(s.controlTick==(callbacks.propose?(s.completedFE+s.controlEveryFESteps-1)/s.controlEveryFESteps:0),"Axial controller tick mismatch");
    require(bits(s.physical.time,targetTime(s.completedFE)),"Axial physical time differs from canonical original clock");
    const auto readoutAssembly=driver->certifyAcceptedState(s.physical);
    // Compare through the exact codec to preserve signed zero and every input/work field.
    auto expected=s;expected.initialPhysicalReference=initialPhysical;
    // The complete refresh above validated this unchanged state through the
    // same approved frozen system used to construct the native readout.
    expected.readouts=readout->observe(s.physical,readoutAssembly);
    expected.physical.appliedInputs=inputs(s.physical.time,s.lastApplied);
    require(RStructuralPiezoAxialRestartCodec::canonicalEquals(expected,s),"Axial physical boundary/reference/readout differs from frozen source and actual system");
    require(bits(s.physical.initialTotalEnergy,initialPhysical.initialTotalEnergy),"Axial initial mechanical energy changed");
    require(s.lastApplied.size()==actuatorIds.size() && s.heldRequested.size()==s.lastApplied.size() && s.envelopes.size()==callbacks.envelopes.size(),"Axial held/applied voltage dimensions mismatch");
    for(size_t i=0;i<s.envelopes.size();++i){const auto &a=s.envelopes[i],&b=callbacks.envelopes[i];require(bits(a.minimum,b.minimum)&&bits(a.maximum,b.maximum)&&bits(a.maximumRate,b.maximumRate),"Axial voltage envelope changed");require(s.lastApplied[i]>=a.minimum && s.lastApplied[i]<=a.maximum,"Axial accepted actuator voltage is outside its amplitude envelope");}
    for(double x:s.heldRequested)require(std::isfinite(x),"Nonfinite held axial voltage request");
    require(s.actions.size()==(s.completedFE?s.lastApplied.size():0),"Axial action table dimensions mismatch");
    const double h=s.completedFE?targetTime(s.completedFE)-targetTime(s.completedFE-1):0;
    for(size_t i=0;i<s.actions.size();++i){const auto &a=s.actions[i];require(a.portId==actuatorIds[i] && bits(a.requested,s.heldRequested[i]) && bits(a.applied,s.lastApplied[i]),"Axial action table differs from held/applied boundary");
        double amp=a.requested,applied=amp;if(!s.envelopes.empty()){const auto &e=s.envelopes[i];require(a.previousApplied>=e.minimum && a.previousApplied<=e.maximum,"Axial previous voltage outside amplitude envelope");amp=std::clamp(amp,e.minimum,e.maximum);const double delta=e.maximumRate*h;require(std::isfinite(delta),"Axial voltage rate overflow");applied=std::clamp(amp,a.previousApplied-delta,a.previousApplied+delta);}
        require(bits(amp,a.amplitudeLimited) && bits(applied,a.applied) && a.amplitudeClipped==(amp!=a.requested) && a.rateClipped==(applied!=amp),"Axial action conditioner record mismatch");}
    if(!s.completedFE){require((callbacks.prepareFE || s.opaqueControllerState==callbacks.initialOpaqueControllerState) && s.lastFrameBytes.isEmpty(),"Axial initial controller state changed");auto baseline=s;baseline.physical=initialPhysical;require(RStructuralPiezoAxialRestartCodec::canonicalEquals(baseline,s),"Axial zero-step state changed");}
    if(callbacks.validateOpaque)callbacks.validateOpaque(s.opaqueControllerState,s.completedFE,s.controlTick,s.controlPhase);
    if(callbacks.validateAcceptedBoundary)callbacks.validateAcceptedBoundary(s);
    std::vector<double> initialVa;for(auto member:actuatorMembers)initialVa.push_back(std::get<RStructuralTrussVoltage>(initialPhysical.appliedInputs.memberInputs[member]).value);
    if(!s.completedFE || !callbacks.propose)require(sameVector(s.heldRequested,initialVa) && sameVector(s.lastApplied,initialVa),"Axial initial/held no-control voltage changed");
    if(s.completedFE==1 || !callbacks.propose)for(size_t i=0;i<s.actions.size();++i)require(bits(s.actions[i].previousApplied,initialVa[i]),"Axial initial or no-control previous voltage changed");
    if(!callbacks.propose)require(s.opaqueControllerState.isEmpty() && s.lastFrameBytes.isEmpty(),"No-control runtime contains hidden controller data");
}
RModel RStructuralPiezoAxialNativeRun::candidate(const RModel &source,const RStructuralPiezoAxialRunState &s) const {
    validateCurrentSource(source);certifyRuntime(s);RModel m(source);RStructuralPiezoAxialRestart r;r.available=true;r.sourceFingerprint=frozenSource;r.runtimeState=RStructuralPiezoAxialRestartCodec::encode(s);m.getProblemSetup().setStructuralPiezoAxialRestart(r);
    auto &time=m.getTimeSolver();std::vector<double> times;times.reserve(target);for(uint i=1;i<=target;++i)times.push_back(targetTime(i));time.setTimes(times);time.setEnabled(true);time.setInputNTimeSteps(target);time.setInputTimeStepSize(dt);time.setOutputFrequency(source.getProblemSetup().getStructuralDynamicsSetup().getOutputFrequency());time.setCurrentTimeStep(s.completedFE?s.completedFE-1:0);time.setComputedTime(s.physical.time);
    m.removeAllVariables();const auto &p=s.physical;
    const auto assembly=RStructuralTrussAssembly::evaluate(actualSystem.referenceNodes,p.displacements,actualSystem.members,p.appliedInputs.memberInputs);
    std::vector<double> axial;for(const auto &member:assembly.memberResponses)axial.push_back(member.axialForce);
    RModel outputDefinition=sourceDefinition(m);RStructuralStaticTrussModel::publishMechanicalMixed(outputDefinition,p.displacements,p.reactions,axial);
    static_cast<RResults &>(m)=static_cast<const RResults &>(outputDefinition);
    m.getVariable(m.findVariable(R_VARIABLE_FORCE)).setName("Force reaction");
    publish(m,R_VARIABLE_VELOCITY,p.velocities,"Velocity");publish(m,R_VARIABLE_ACCELERATION,p.accelerations,"Acceleration");return m;
}
RModel RStructuralPiezoAxialNativeRun::acceptedModel(const RModel &source) const{return candidate(source,*accepted);}
RModel RStructuralPiezoAxialNativeRun::modelForAcceptedRuntime(const RModel &source,const RStructuralPiezoAxialRunState &state) const {
    auto model=candidate(source,state);
    certifyAcceptedRuntime(model);
    return model;
}
std::unique_ptr<RStructuralPiezoAxialNativeRun::PreparedStep> RStructuralPiezoAxialNativeRun::prepare(const RModel &source,bool boundary) const {
    require(retryable,"Axial native run cannot retry after incomplete publication rollback");validateCurrentSource(source);require(hasRemainingSteps(),"Axial native run reached original target");
    auto next=std::make_unique<RStructuralPiezoAxialRunState>(*accepted);
    if(callbacks.prepareFE)next->opaqueControllerState=callbacks.prepareFE(targetTime(accepted->completedFE),accepted->completedFE,accepted->opaqueControllerState);
    if(callbacks.bindForceState)callbacks.bindForceState(next->opaqueControllerState);
    if(callbacks.propose && accepted->completedFE%callbacks.controlEveryFESteps==0){const auto p=callbacks.propose(accepted->physical,accepted->readouts,accepted->completedFE,accepted->controlTick,next->opaqueControllerState);require(p.requestedVoltages.size()==accepted->lastApplied.size(),"Controller proposal voltage dimensions mismatch");next->heldRequested=p.requestedVoltages;next->opaqueControllerState=p.nextOpaqueState;next->lastFrameBytes=p.lastFrameBytes;++next->controlTick;validateCurrentSource(source);}
    const double endpoint=targetTime(accepted->completedFE+1);require(endpoint>accepted->physical.time,"Axial physical FE interval does not advance time");
    const double h=endpoint-accepted->physical.time;
    next->actions.clear();next->lastApplied.resize(next->heldRequested.size());const auto &ids=actuatorIds;
    for(size_t i=0;i<next->heldRequested.size();++i){RStructuralPiezoAxialVoltageAction a;a.portId=ids[i];a.requested=next->heldRequested[i];require(std::isfinite(a.requested),"Nonfinite controller voltage proposal");a.previousApplied=accepted->lastApplied[i];a.amplitudeLimited=a.requested;a.applied=a.requested;
        if(!next->envelopes.empty()){const auto &e=next->envelopes[i];a.amplitudeLimited=std::clamp(a.requested,e.minimum,e.maximum);const double delta=e.maximumRate*h;require(std::isfinite(delta),"Voltage rate overflow");a.applied=std::clamp(a.amplitudeLimited,a.previousApplied-delta,a.previousApplied+delta);}
        a.amplitudeClipped=a.amplitudeLimited!=a.requested;a.rateClipped=a.applied!=a.amplitudeLimited;next->lastApplied[i]=a.applied;next->actions.push_back(a);
    }
    auto applied=inputs(accepted->physical.time,next->lastApplied);applied.rightExternalForces=totalForce(endpoint);
    if(callbacks.intervalPerturbationForce){
        const auto &background=frozenAdapter.initialCompleteInputs().externalForces;
        auto interval=[&](double time){auto force=callbacks.intervalPerturbationForce(time);require(force.size()==background.size(),"Axial interval force dimensions mismatch");for(size_t i=0;i<force.size();++i){require(std::isfinite(force[i]),"Nonfinite axial interval force");force[i]+=background[i];require(std::isfinite(force[i]),"Axial interval total force overflow");}return force;};
        applied.externalForces=interval(accepted->physical.time);applied.rightExternalForces=interval(endpoint);
    }
    next->physical=driver->trialStepOnGrid(accepted->physical,applied,initialTime,std::uint64_t(accepted->completedFE)+1,dt,frozenAdapter.options()).candidate;next->readouts=readout->observe(next->physical);++next->completedFE;next->controlPhase=next->completedFE%next->controlEveryFESteps;
    if(callbacks.postFE)next->opaqueControllerState=callbacks.postFE(next->physical,next->readouts,next->completedFE,next->controlTick,next->lastApplied,next->opaqueControllerState);
    validateCurrentSource(source);certifyRuntime(*next);
    std::unique_ptr<PreparedStep> step(new PreparedStep);step->owner=this;step->base=accepted.get();step->next=std::move(next);step->sealedRuntime=std::make_unique<const RStructuralPiezoAxialRunState>(*step->next);step->model=std::make_unique<RModel>(candidate(source,*step->next));step->sealedResults=std::make_unique<RResults>(*step->model);
    step->stepIdentity.sequence=step->next->completedFE;step->stepIdentity.timeS=step->next->physical.time;step->stepIdentity.dtS=dt;step->stepIdentity.finalStep=boundary||step->next->completedFE==target;return step;
}
void RStructuralPiezoAxialNativeRun::validatePrepared(const PreparedStep &s) const {require(bool(s.model),"Axial candidate missing");validatePrepared(s,*s.model);}
void RStructuralPiezoAxialNativeRun::validatePrepared(const PreparedStep &s,const RModel &m) const {
    require(retryable && s.owner==this && s.base==accepted.get() && !s.committed && s.next && s.sealedRuntime && s.model && s.sealedResults && s.next->completedFE==accepted->completedFE+1,"Axial prepared step is stale or foreign");validateCurrentSource(m);certifyRuntime(*s.next);require(resultsEqual(m,*s.sealedResults),"Axial candidate result/display metadata changed");
    const auto &r=m.getProblemSetup().getStructuralPiezoAxialRestart();require(r.available && r.sourceFingerprint==frozenSource && r.runtimeState==RStructuralPiezoAxialRestartCodec::encode(*s.next),"Axial candidate complete runtime bytes changed");
    const auto &t=m.getTimeSolver();require(t.getEnabled() && bits(t.getComputedTime(),s.stepIdentity.timeS) && t.getCurrentTimeStep()==s.next->completedFE-1 && t.getInputNTimeSteps()==target && bits(t.getInputTimeStepSize(),dt) && t.getTimes().size()==target && t.getOutputFrequency()==m.getProblemSetup().getStructuralDynamicsSetup().getOutputFrequency(),"Axial candidate clock changed");
    for(uint i=0;i<target;++i)require(bits(t.getTimes()[i],targetTime(i+1)),"Axial original target grid changed");
}
void RStructuralPiezoAxialNativeRun::commit(PreparedStep &s) noexcept {
    assert(retryable && s.owner==this && s.base==accepted.get() && !s.committed && s.next);static_assert(noexcept(accepted.swap(s.next)),"Axial promotion must not throw");accepted.swap(s.next);s.committed=true;
}

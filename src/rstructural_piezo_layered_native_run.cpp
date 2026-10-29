#include "rstructural_piezo_layered_native_run.h"
#include "rstructural_piezo_layered_restart_codec.h"
#include <rml_model.h>
#include <QCryptographicHash>
#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstring>
#include <stdexcept>

namespace {
void require(bool b,const char *message) {if(!b)throw std::invalid_argument(message);}
bool bits(double a,double b) {return std::memcmp(&a,&b,8)==0;}
bool sameVector(const std::vector<double> &a,const std::vector<double> &b) {return a.size()==b.size() && (a.empty() || std::memcmp(a.data(),b.data(),8*a.size())==0);}
void integer(QByteArray &b,quint64 n) {for(uint i=0;i<8;++i)b.append(char((n>>(8*i))&255));}
void blob(QByteArray &b,const QByteArray &v) {integer(b,quint64(v.size()));b.append(v);}
void real(QByteArray &b,double v) {require(std::isfinite(v),"Nonfinite layered runtime configuration");quint64 n;std::memcpy(&n,&v,8);integer(b,n);}
bool resultsEqual(const RResults &a,const RResults &b) {
    if(a.getNNodes()!=b.getNNodes() || a.getNElements()!=b.getNElements() || a.getNVariables()!=b.getNVariables())return false;
    for(uint i=0;i<a.getNVariables();++i){const auto &x=a.getVariable(i),&y=b.getVariable(i);
        if(x.getType()!=y.getType() || x.getApplyType()!=y.getApplyType() || x.getName()!=y.getName() || x.getUnits()!=y.getUnits() || x.getNVectors()!=y.getNVectors() || x.getNValues()!=y.getNValues())return false;
        const auto &xd=x.getVariableData(),&yd=y.getVariableData();
        if(!bits(xd.getMinDisplayValue(),yd.getMinDisplayValue()) || !bits(xd.getMaxDisplayValue(),yd.getMaxDisplayValue()) || xd.getDisplayValueRange()!=yd.getDisplayValueRange() || xd.getValueRangeName()!=yd.getValueRangeName() || !bits(xd.getScale(),yd.getScale()))return false;
        for(uint c=0;c<x.getNVectors();++c)for(uint n=0;n<x.getNValues();++n)if(!bits(x.getValue(c,n),y.getValue(c,n)))return false;
    }return true;
}
void publish(RModel &m,const std::vector<uint> &ids,RVariableType type,const std::vector<double> &values,uint first,const QString &name) {
    auto &v=m.getVariable(m.addVariable(type));v.setName(name);v.setApplyType(R_VARIABLE_APPLY_NODE);v.resize(3,m.getNNodes());double maximum=0;
    for(uint n=0;n<m.getNNodes();++n){const auto it=std::lower_bound(ids.begin(),ids.end(),m.getNodePhysicalId(n));require(it!=ids.end() && *it==m.getNodePhysicalId(n),"Layered result node identity mismatch");const size_t k=size_t(it-ids.begin());
        for(uint c=0;c<3;++c)v.setValue(c,n,values[6*k+first+c]);maximum=std::max(maximum,std::hypot(values[6*k+first],std::hypot(values[6*k+first+1],values[6*k+first+2])));
    }require(std::isfinite(maximum),"Layered display magnitude overflow");v.getVariableData().setMinMaxDisplayValue(0,maximum);
}
}
RStructuralPiezoLayeredNativeRun::PreparedStep::PreparedStep()=default;
RStructuralPiezoLayeredNativeRun::PreparedStep::~PreparedStep()=default;
RModel &RStructuralPiezoLayeredNativeRun::PreparedStep::candidateModel() noexcept{return *model;}
const RModel &RStructuralPiezoLayeredNativeRun::PreparedStep::candidateModel() const noexcept{return *model;}
const RAcceptedStepIdentity &RStructuralPiezoLayeredNativeRun::PreparedStep::identity() const noexcept{return stepIdentity;}
const RStructuralPiezoLayeredDynamicsState &RStructuralPiezoLayeredNativeRun::PreparedStep::result() const noexcept{return sealedRuntime->physical;}
const RStructuralPiezoLayeredRunState &RStructuralPiezoLayeredNativeRun::PreparedStep::runtime() const noexcept{return *sealedRuntime;}
uint RStructuralPiezoLayeredNativeRun::PreparedStep::completedSteps() const noexcept{return uint(stepIdentity.sequence);}

RStructuralPiezoLayeredNativeRun::RStructuralPiezoLayeredNativeRun(const RModel &source,RStructuralPiezoLayeredRuntimeCallbacks bundle,RMatrixSolverConf conf)
    :callbacks(std::move(bundle)),frozenAdapter(source,source.getProblemSetup().getStructuralPiezoLayeredSetup(),callbacks.authoritativeFullC) {
    const auto &setup=source.getProblemSetup().getStructuralPiezoLayeredSetup();const auto &d=source.getProblemSetup().getStructuralDynamicsSetup();
    d.validate();require(setup.enabled && d.getEnabled() && d.getMethod()==RStructuralTimeIntegrationMethod::Newmark && d.getNewmarkGamma()==.5 && d.getNewmarkBeta()==.25,"Layered facade requires enabled average-acceleration Newmark source");
    dt=d.getTimeStep();target=d.getNSteps();initialTime=source.getTimeSolver().getInputStartTime()-dt;
    require(std::isfinite(initialTime) && target>0 && target<=1000000,"Invalid layered source clock");
    frozenSource=setup.computeSourceFingerprint(source);
    const bool custom=bool(callbacks.prepareFE)||bool(callbacks.bindForceState)||bool(callbacks.intervalPerturbationForce)||bool(callbacks.perturbationForce)||bool(callbacks.propose)||bool(callbacks.postFE)||bool(callbacks.validateOpaque)||bool(callbacks.validateAcceptedBoundary)||!callbacks.authoritativeFullC.empty()||!callbacks.backendIdentity.isEmpty()||!callbacks.sourceFingerprint.isEmpty();
    if(custom){require(callbacks.sourceFingerprint==frozenSource && !callbacks.backendIdentity.isEmpty(),"Layered runtime bundle must bind actual source and frozen backend identity");
        require(callbacks.excitationConfiguration==setup.excitationConfiguration && callbacks.controllerConfiguration==setup.controllerConfiguration && callbacks.designConfiguration==setup.designConfiguration,"Layered runtime declaration differs from source configuration");}
    else require(setup.excitationConfiguration.isEmpty() && setup.controllerConfiguration.isEmpty() && setup.designConfiguration.isEmpty(),"Nonempty layered runtime configuration requires a matching callbacks bundle");
    require(setup.excitationConfiguration.isEmpty() || bool(callbacks.perturbationForce),"Declared layered excitation has no time-only callback");
    require(setup.controllerConfiguration.isEmpty() || bool(callbacks.propose),"Declared layered controller has no proposal callback");
    require(callbacks.controlEveryFESteps>0 && callbacks.controlEveryFESteps<=1000000,"Invalid layered controller FE stride");
    if(callbacks.propose)require(!callbacks.controllerIdentity.isEmpty() && (bool(callbacks.validateOpaque)||bool(callbacks.validateAcceptedBoundary)),"Layered controller needs identity and semantic opaque validator");
    else require(!callbacks.postFE && !callbacks.validateOpaque && !callbacks.validateAcceptedBoundary && callbacks.initialOpaqueControllerState.isEmpty() && callbacks.controllerIdentity.isEmpty() && callbacks.envelopes.empty(),"Controller runtime without proposal is unsupported");
    if(!callbacks.perturbationForce){const size_t n=6*frozenAdapter.mesh().nodeIds.size();callbacks.perturbationForce=[n](double){return std::vector<double>(n,0);};callbacks.forcingIdentity="range-fea/layered-zero-time-force/v1";}
    require(!callbacks.forcingIdentity.isEmpty(),"Layered time callback needs declared identity");
    if(callbacks.backendIdentity.isEmpty())callbacks.backendIdentity="range-fea/native-layered-Newmark/no-control/v1";
    const size_t na=frozenAdapter.electrical().condensed().actuatorPortIds.size();require(callbacks.envelopes.empty() || callbacks.envelopes.size()==na,"Voltage envelope dimensions differ from actuator order");
    QByteArray config;blob(config,frozenSource.toUtf8());blob(config,setup.identityBytes());blob(config,callbacks.backendIdentity.toUtf8());blob(config,callbacks.forcingIdentity.toUtf8());blob(config,callbacks.controllerIdentity.toUtf8());integer(config,callbacks.controlEveryFESteps);
    blob(config,callbacks.initialOpaqueControllerState);integer(config,bool(callbacks.postFE));integer(config,bool(callbacks.propose));integer(config,bool(callbacks.validateAcceptedBoundary));integer(config,callbacks.authoritativeFullC.size());for(double x:callbacks.authoritativeFullC)real(config,x);
    integer(config,callbacks.envelopes.size());for(const auto &e:callbacks.envelopes){require(e.minimum<=e.maximum && e.maximumRate>=0,"Invalid voltage envelope");real(config,e.minimum);real(config,e.maximum);real(config,e.maximumRate);}
    for(size_t i=0;i<callbacks.envelopes.size();++i){const auto &e=callbacks.envelopes[i];const double initial=frozenAdapter.data().actuatorVoltages0[i];require(initial>=e.minimum && initial<=e.maximum,"Layered source initial actuator voltage is outside its amplitude envelope");}
    integer(config,conf.getType());integer(config,conf.getNInnerIterations());integer(config,conf.getNOuterIterations());real(config,conf.getSolverCvgValue());real(config,initialTime);real(config,dt);integer(config,target);
    require(bool(callbacks.prepareFE)==bool(callbacks.bindForceState) && bool(callbacks.prepareFE)==bool(callbacks.intervalPerturbationForce),"Incomplete transactional force callbacks");
    require(!callbacks.prepareFE || (callbacks.bindForceState && callbacks.intervalPerturbationForce && callbacks.propose && callbacks.postFE),"Sampled forcing requires a complete transactional runtime bundle");
    if(callbacks.prepareFE){blob(config,"transactional-interval-force/v1");integer(config,1);}
    frozenConfiguration=QCryptographicHash::hash(config,QCryptographicHash::Sha256);
    driver=std::make_unique<RStructuralPiezoLayeredDynamics>(frozenAdapter,conf,initialTime,callbacks.perturbationForce,callbacks.forcingIdentity);
    auto state=std::make_unique<RStructuralPiezoLayeredRunState>();
    if(source.getProblemSetup().getRestart()){
        const auto &r=source.getProblemSetup().getStructuralPiezoLayeredRestart();r.validate();require(r.available && r.sourceFingerprint==frozenSource,"Layered restart needs matching complete runtime authority");*state=RStructuralPiezoLayeredRestartCodec::decode(r.stateBytes);
    }else{
        state->opaqueControllerState=callbacks.initialOpaqueControllerState;
        if(callbacks.prepareFE)state->opaqueControllerState=callbacks.prepareFE(initialTime,0,state->opaqueControllerState);
        if(callbacks.bindForceState)callbacks.bindForceState(state->opaqueControllerState);
        state->physical=driver->initialize();state->sourceFingerprint=frozenSource;state->forcingIdentity=callbacks.forcingIdentity;state->controllerIdentity=callbacks.controllerIdentity;state->backendIdentity=callbacks.backendIdentity;state->configurationIdentity=frozenConfiguration;
        state->targetSteps=target;state->controlEveryFESteps=callbacks.controlEveryFESteps;state->initialTime=initialTime;state->dt=dt;state->nominalControlTs=double(callbacks.controlEveryFESteps)*dt;state->envelopes=callbacks.envelopes;
        state->heldRequested=state->lastApplied=state->physical.heldActuatorVoltages;
    }
    certifyRuntime(*state);accepted=std::move(state);
}
RStructuralPiezoLayeredNativeRun::~RStructuralPiezoLayeredNativeRun()=default;
const RStructuralPiezoLayeredModelAdapter &RStructuralPiezoLayeredNativeRun::adapter() const noexcept{return frozenAdapter;}
const RStructuralPiezoLayeredDynamicsState &RStructuralPiezoLayeredNativeRun::acceptedState() const noexcept{return accepted->physical;}
const RStructuralPiezoLayeredRunState &RStructuralPiezoLayeredNativeRun::acceptedRuntime() const noexcept{return *accepted;}
uint RStructuralPiezoLayeredNativeRun::completedSteps() const noexcept{return accepted->completedFE;}
uint RStructuralPiezoLayeredNativeRun::targetSteps() const noexcept{return target;}
bool RStructuralPiezoLayeredNativeRun::hasRemainingSteps() const noexcept{return accepted->completedFE<target;}
bool RStructuralPiezoLayeredNativeRun::isRetryable() const noexcept{return retryable;}
void RStructuralPiezoLayeredNativeRun::markNonRetryable() noexcept{retryable=false;}
double RStructuralPiezoLayeredNativeRun::targetTime(uint step) const {const double t=step?initialTime+double(step)*dt:initialTime;require(std::isfinite(t),"Layered canonical clock overflow");return t;}
void RStructuralPiezoLayeredNativeRun::validateCurrentSource(const RModel &source) const {
    const auto &s=source.getProblemSetup().getStructuralPiezoLayeredSetup();s.validate(source);require(s.computeSourceFingerprint(source)==frozenSource,"Layered native source changed since freeze");
}
RStructuralPiezoLayeredRunState RStructuralPiezoLayeredNativeRun::certifyAcceptedRuntime(const RModel &model) const {
    validateCurrentSource(model);
    const auto &restart=model.getProblemSetup().getStructuralPiezoLayeredRestart();restart.validate();
    require(restart.available && restart.sourceFingerprint==frozenSource,"Layered accepted MODEL requires complete matching restart authority");
    auto state=RStructuralPiezoLayeredRestartCodec::decode(restart.stateBytes);certifyRuntime(state);
    const auto &time=model.getTimeSolver();
    require(time.getEnabled() && bits(time.getComputedTime(),state.physical.integrator.state.time)
        && time.getCurrentTimeStep()==(state.completedFE?state.completedFE-1:0)
        && time.getInputNTimeSteps()==target && bits(time.getInputTimeStepSize(),dt)
        && time.getTimes().size()==target
        && time.getOutputFrequency()==model.getProblemSetup().getStructuralDynamicsSetup().getOutputFrequency(),"Layered accepted MODEL clock changed");
    for(uint i=0;i<target;++i)require(bits(time.getTimes()[i],targetTime(i+1)),"Layered original accepted target grid changed");
    return state;
}
void RStructuralPiezoLayeredNativeRun::certifyAcceptedRuntime(const RStructuralPiezoLayeredRunState &state) const {
    certifyRuntime(state);
}
void RStructuralPiezoLayeredNativeRun::certifyRuntime(const RStructuralPiezoLayeredRunState &s) const {
    if(callbacks.bindForceState)callbacks.bindForceState(s.opaqueControllerState);
    require(s.sourceFingerprint==frozenSource && s.configurationIdentity==frozenConfiguration && s.forcingIdentity==callbacks.forcingIdentity && s.controllerIdentity==callbacks.controllerIdentity && s.backendIdentity==callbacks.backendIdentity,"Layered runtime source/backend/configuration mismatch");
    require(s.targetSteps==target && s.completedFE<=target && s.controlEveryFESteps==callbacks.controlEveryFESteps && s.controlPhase==s.completedFE%s.controlEveryFESteps && bits(s.initialTime,initialTime) && bits(s.dt,dt) && bits(s.nominalControlTs,double(s.controlEveryFESteps)*dt),"Layered runtime clock/control phase mismatch");
    require(s.controlTick==(callbacks.propose?(s.completedFE+s.controlEveryFESteps-1)/s.controlEveryFESteps:0),"Layered controller tick mismatch");
    require(bits(s.physical.integrator.state.time,targetTime(s.completedFE)),"Layered physical time differs from canonical original clock");
    require(bits(s.physical.integrator.fixedStep,s.completedFE?dt:0.0),"Layered physical fixed step differs from frozen nominal FE interval");
    driver->restore(s.physical);
    require(sameVector(s.lastApplied,s.physical.heldActuatorVoltages) && s.heldRequested.size()==s.lastApplied.size() && s.envelopes.size()==callbacks.envelopes.size(),"Layered held/applied voltage dimensions mismatch");
    for(size_t i=0;i<s.envelopes.size();++i){const auto &a=s.envelopes[i],&b=callbacks.envelopes[i];require(bits(a.minimum,b.minimum)&&bits(a.maximum,b.maximum)&&bits(a.maximumRate,b.maximumRate),"Layered voltage envelope changed");require(s.lastApplied[i]>=a.minimum && s.lastApplied[i]<=a.maximum,"Layered accepted actuator voltage is outside its amplitude envelope");}
    for(double x:s.heldRequested)require(std::isfinite(x),"Nonfinite held layered voltage request");
    require(s.actions.size()==(s.completedFE?s.lastApplied.size():0),"Layered action table dimensions mismatch");
    const auto &ids=driver->actuatorPortIds();
    const double h=s.completedFE?targetTime(s.completedFE)-targetTime(s.completedFE-1):0;
    for(size_t i=0;i<s.actions.size();++i){const auto &a=s.actions[i];require(a.portId==ids[i] && bits(a.requested,s.heldRequested[i]) && bits(a.applied,s.lastApplied[i]),"Layered action table differs from held/applied boundary");
        double amp=a.requested,applied=amp;if(!s.envelopes.empty()){const auto &e=s.envelopes[i];require(a.previousApplied>=e.minimum && a.previousApplied<=e.maximum,"Layered previous applied actuator voltage is outside its amplitude envelope");amp=std::clamp(amp,e.minimum,e.maximum);const double delta=e.maximumRate*h;require(std::isfinite(delta),"Layered voltage rate overflow");applied=std::clamp(amp,a.previousApplied-delta,a.previousApplied+delta);}
        require(bits(amp,a.amplitudeLimited) && bits(applied,a.applied) && a.amplitudeClipped==(amp!=a.requested) && a.rateClipped==(applied!=amp),"Layered action conditioner record mismatch");}
    if(callbacks.validateOpaque)callbacks.validateOpaque(s.opaqueControllerState,s.completedFE,s.controlTick,s.controlPhase);
    if(callbacks.validateAcceptedBoundary)callbacks.validateAcceptedBoundary(s);
    if(!callbacks.validateOpaque && !callbacks.validateAcceptedBoundary)require(s.opaqueControllerState.isEmpty() && s.lastFrameBytes.isEmpty(),"No-control runtime contains hidden controller data");
    RStructuralPiezoLayeredRestartCodec::encode(s); // Bound every opaque/vector field before publication.
}
RModel RStructuralPiezoLayeredNativeRun::candidate(const RModel &source,const RStructuralPiezoLayeredRunState &s) const {
    validateCurrentSource(source);certifyRuntime(s);RModel m(source);RStructuralPiezoLayeredRestart r;r.available=true;r.sourceFingerprint=frozenSource;r.stateBytes=RStructuralPiezoLayeredRestartCodec::encode(s);m.getProblemSetup().setStructuralPiezoLayeredRestart(r);
    auto &time=m.getTimeSolver();std::vector<double> times;times.reserve(target);for(uint i=1;i<=target;++i)times.push_back(targetTime(i));time.setTimes(times);time.setEnabled(true);time.setInputNTimeSteps(target);time.setInputTimeStepSize(dt);time.setOutputFrequency(source.getProblemSetup().getStructuralDynamicsSetup().getOutputFrequency());time.setCurrentTimeStep(s.completedFE?s.completedFE-1:0);time.setComputedTime(s.physical.integrator.state.time);
    m.removeAllVariables();const auto &ids=frozenAdapter.mesh().nodeIds;const auto &p=s.physical;
    publish(m,ids,R_VARIABLE_DISPLACEMENT,p.physicalDisplacements,0,"Displacement");publish(m,ids,R_VARIABLE_ROTATION,p.physicalDisplacements,3,"Rotation");
    publish(m,ids,R_VARIABLE_VELOCITY,p.physicalVelocities,0,"Velocity");publish(m,ids,R_VARIABLE_ANGULAR_VELOCITY,p.physicalVelocities,3,"Angular velocity");
    publish(m,ids,R_VARIABLE_ACCELERATION,p.physicalAccelerations,0,"Acceleration");publish(m,ids,R_VARIABLE_ANGULAR_ACCELERATION,p.physicalAccelerations,3,"Angular acceleration");
    publish(m,ids,R_VARIABLE_FORCE,p.physicalReactions,0,"Force reaction");publish(m,ids,R_VARIABLE_MOMENT,p.physicalReactions,3,"Moment reaction");return m;
}
RModel RStructuralPiezoLayeredNativeRun::acceptedModel(const RModel &source) const{return candidate(source,*accepted);}
RModel RStructuralPiezoLayeredNativeRun::modelForAcceptedRuntime(const RModel &source,const RStructuralPiezoLayeredRunState &state) const {
    auto model=candidate(source,state);
    certifyAcceptedRuntime(model);
    return model;
}
std::unique_ptr<RStructuralPiezoLayeredNativeRun::PreparedStep> RStructuralPiezoLayeredNativeRun::prepare(const RModel &source,bool boundary) const {
    require(retryable,"Layered native run cannot retry after incomplete publication rollback");validateCurrentSource(source);require(hasRemainingSteps(),"Layered native run reached original target");
    auto next=std::make_unique<RStructuralPiezoLayeredRunState>(*accepted);
    if(callbacks.prepareFE)next->opaqueControllerState=callbacks.prepareFE(targetTime(accepted->completedFE),accepted->completedFE,accepted->opaqueControllerState);
    if(callbacks.bindForceState)callbacks.bindForceState(next->opaqueControllerState);
    if(callbacks.propose && accepted->completedFE%callbacks.controlEveryFESteps==0){const auto p=callbacks.propose(accepted->physical,accepted->completedFE,accepted->controlTick,next->opaqueControllerState);require(p.requestedVoltages.size()==accepted->lastApplied.size(),"Controller proposal voltage dimensions mismatch");next->heldRequested=p.requestedVoltages;next->opaqueControllerState=p.nextOpaqueState;next->lastFrameBytes=p.lastFrameBytes;++next->controlTick;validateCurrentSource(source);}
    const double endpoint=targetTime(accepted->completedFE+1);require(endpoint>accepted->physical.integrator.state.time,"Layered physical FE interval does not advance time");
    const double h=endpoint-accepted->physical.integrator.state.time;
    next->actions.clear();next->lastApplied.resize(next->heldRequested.size());const auto &ids=driver->actuatorPortIds();
    for(size_t i=0;i<next->heldRequested.size();++i){RStructuralPiezoLayeredVoltageAction a;a.portId=ids[i];a.requested=next->heldRequested[i];require(std::isfinite(a.requested),"Nonfinite controller voltage proposal");a.previousApplied=accepted->lastApplied[i];a.amplitudeLimited=a.requested;a.applied=a.requested;
        if(!next->envelopes.empty()){const auto &e=next->envelopes[i];a.amplitudeLimited=std::clamp(a.requested,e.minimum,e.maximum);const double delta=e.maximumRate*h;require(std::isfinite(delta),"Voltage rate overflow");a.applied=std::clamp(a.amplitudeLimited,a.previousApplied-delta,a.previousApplied+delta);}
        a.amplitudeClipped=a.amplitudeLimited!=a.requested;a.rateClipped=a.applied!=a.amplitudeLimited;next->lastApplied[i]=a.applied;next->actions.push_back(a);
    }
    next->physical=driver->trialOnGrid(accepted->physical,next->lastApplied,std::uint64_t(accepted->completedFE)+1,dt,callbacks.intervalPerturbationForce);++next->completedFE;next->controlPhase=next->completedFE%next->controlEveryFESteps;
    if(callbacks.postFE)next->opaqueControllerState=callbacks.postFE(next->physical,next->completedFE,next->controlTick,next->lastApplied,next->opaqueControllerState);
    validateCurrentSource(source);certifyRuntime(*next);
    std::unique_ptr<PreparedStep> step(new PreparedStep);step->owner=this;step->base=accepted.get();step->next=std::move(next);step->sealedRuntime=std::make_unique<const RStructuralPiezoLayeredRunState>(*step->next);step->model=std::make_unique<RModel>(candidate(source,*step->next));step->sealedResults=std::make_unique<RResults>(*step->model);
    step->stepIdentity.sequence=step->next->completedFE;step->stepIdentity.timeS=step->next->physical.integrator.state.time;step->stepIdentity.dtS=dt;step->stepIdentity.finalStep=boundary||step->next->completedFE==target;return step;
}
void RStructuralPiezoLayeredNativeRun::validatePrepared(const PreparedStep &s) const {require(bool(s.model),"Layered candidate missing");validatePrepared(s,*s.model);}
void RStructuralPiezoLayeredNativeRun::validatePrepared(const PreparedStep &s,const RModel &m) const {
    require(retryable && s.owner==this && s.base==accepted.get() && !s.committed && s.next && s.sealedRuntime && s.model && s.sealedResults && s.next->completedFE==accepted->completedFE+1,"Layered prepared step is stale or foreign");validateCurrentSource(m);certifyRuntime(*s.next);require(resultsEqual(m,*s.sealedResults),"Layered candidate result/display metadata changed");
    const auto &r=m.getProblemSetup().getStructuralPiezoLayeredRestart();require(r.available && r.sourceFingerprint==frozenSource && r.stateBytes==RStructuralPiezoLayeredRestartCodec::encode(*s.next),"Layered candidate complete runtime bytes changed");
    const auto &t=m.getTimeSolver();require(t.getEnabled() && bits(t.getComputedTime(),s.stepIdentity.timeS) && t.getCurrentTimeStep()==s.next->completedFE-1 && t.getInputNTimeSteps()==target && bits(t.getInputTimeStepSize(),dt) && t.getTimes().size()==target && t.getOutputFrequency()==m.getProblemSetup().getStructuralDynamicsSetup().getOutputFrequency(),"Layered candidate clock changed");
    for(uint i=0;i<target;++i)require(bits(t.getTimes()[i],targetTime(i+1)),"Layered original target grid changed");
}
void RStructuralPiezoLayeredNativeRun::commit(PreparedStep &s) noexcept {
    assert(retryable && s.owner==this && s.base==accepted.get() && !s.committed && s.next);static_assert(noexcept(accepted.swap(s.next)),"Layered promotion must not throw");accepted.swap(s.next);s.committed=true;
}

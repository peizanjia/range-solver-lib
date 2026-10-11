#include "rstructural_native_control.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>
#include <map>
#include <set>
#include <QCryptographicHash>
#include <QSysInfo>
#include <rbl_error.h>

#ifndef RANGE_NATIVE_CONTROL_IMPLEMENTATION_SHA256
#error "Native control provider requires the configured implementation SHA256."
#endif

namespace
{
[[noreturn]] void fail(const char *message)
{
    throw RError(RError::Type::InvalidInput,R_ERROR_REF,"%s",message);
}
bool sameBits(double a,double b) { return std::memcmp(&a,&b,sizeof(double)) == 0; }
bool sameInputs(const RStructuralTrussDynamicsInputs &a,const RStructuralTrussDynamicsInputs &b)
{
    if (a.externalForces.size()!=b.externalForces.size() || a.memberInputs.size()!=b.memberInputs.size()) return false;
    for (std::size_t i=0;i<a.externalForces.size();++i)
        if (!sameBits(a.externalForces[i],b.externalForces[i])) return false;
    for (std::size_t i=0;i<a.memberInputs.size();++i)
    {
        const auto &left=a.memberInputs[i],&right=b.memberInputs[i];
        if (left.index()!=right.index()) return false;
        if (const auto *v=std::get_if<RStructuralTrussVoltage>(&left))
        { if (!sameBits(v->value,std::get<RStructuralTrussVoltage>(right).value)) return false; }
        else if (const auto *q=std::get_if<RStructuralTrussCharge>(&left))
        { if (!sameBits(q->value,std::get<RStructuralTrussCharge>(right).value)) return false; }
        else if (const auto *l=std::get_if<RStructuralTrussCableLength>(&left))
        { if (!sameBits(l->value,std::get<RStructuralTrussCableLength>(right).value)) return false; }
    }
    return true;
}
QString sha(const QByteArray &bytes)
{ return QString::fromLatin1(QCryptographicHash::hash(bytes,QCryptographicHash::Sha256).toHex()); }
void appendInteger(QByteArray &bytes,std::uint64_t value)
{ for (unsigned i=0;i<8;++i) bytes.append(char((value>>(8*i))&255)); }
std::uint64_t readInteger(const QByteArray &bytes,qsizetype &position)
{
    if (position < 0 || bytes.size()-position < 8) fail("Native PI state envelope is truncated.");
    std::uint64_t value=0;
    for (unsigned i=0;i<8;++i) value |= std::uint64_t(uchar(bytes[position++]))<<(8*i);
    return value;
}
const QByteArray envelopeDomain("range-fea/native-pi/integral/v1\0",32);
QByteArray encode(const std::vector<double> &integrals)
{
    QByteArray bytes=envelopeDomain;
    appendInteger(bytes,integrals.size());
    for (double value:integrals)
    {
        if (!std::isfinite(value)) fail("Native PI integral is nonfinite.");
        std::uint64_t bits=0; std::memcpy(&bits,&value,sizeof(bits)); appendInteger(bytes,bits);
    }
    return bytes;
}
std::vector<double> decode(const QByteArray &bytes,std::size_t count)
{
    if (!bytes.startsWith(envelopeDomain)) fail("Native PI state envelope has a different schema.");
    qsizetype position=envelopeDomain.size();
    if (readInteger(bytes,position)!=count || count>std::size_t((bytes.size()-position)/8) ||
        bytes.size()-position != qsizetype(count*8)) fail("Native PI state envelope has an invalid integral count.");
    std::vector<double> values;
    for (std::size_t i=0;i<count;++i)
    {
        const auto bits=readInteger(bytes,position); double value=0.0; std::memcpy(&value,&bits,sizeof(value));
        if (!std::isfinite(value)) fail("Native PI state envelope contains a nonfinite integral.");
        values.push_back(value);
    }
    return values;
}
bool sameRecord(const RStructuralNativeControlActionRecord &a,const RStructuralNativeControlActionRecord &b)
{
    return a.actionId==b.actionId && a.quantity==b.quantity && a.limitFlags==b.limitFlags &&
        sameBits(a.previousApplied,b.previousApplied) && sameBits(a.requested,b.requested) &&
        sameBits(a.conditioned,b.conditioned) && sameBits(a.applied,b.applied);
}
}

QByteArray RStructuralNativeControl::builtinEnvironmentManifest()
{
    static_assert(sizeof(double)==8 && std::numeric_limits<double>::is_iec559 &&
                  std::numeric_limits<double>::digits==53,"Native PI requires IEEE binary64.");
    QByteArray bytes("range-fea/native-pi/environment/v1\n");
#if defined(__clang__)
    bytes += "compiler=clang:" __clang_version__ "\n";
#elif defined(_MSC_VER)
    bytes += "compiler=msvc:"; bytes += QByteArray::number(_MSC_VER); bytes += ":";
    bytes += QByteArray::number(_MSC_FULL_VER); bytes += "\n";
#elif defined(__GNUC__)
    bytes += "compiler=gcc:" __VERSION__ "\n";
#else
#error "Native PI environment manifest needs an explicit compiler identity."
#endif
    bytes += "qt.compile=" QT_VERSION_STR "\nqt.runtime=";
    bytes += qVersion(); bytes += "\nabi="; bytes += QSysInfo::buildAbi().toUtf8();
    bytes += "\ncpp="; bytes += QByteArray::number(qulonglong(__cplusplus));
    bytes += "\npointer.bytes="; bytes += QByteArray::number(qulonglong(sizeof(void*)));
    bytes += "\nsize_t.bytes="; bytes += QByteArray::number(qulonglong(sizeof(std::size_t)));
    bytes += "\ndouble=IEEE754-binary64;canonical=little-endian;preserve-signed-zero\n";
    bytes += "clock=canonical-FE-grid;N=1;left-accepted-sample;PI-integrates-nominal-Ts;rate-physical-h\n";
    bytes += "dependencies=QtCore;C++standard-library;no-imports;no-weights;no-device-I/O\n";
    bytes += "scope=declared-compiler-Qt-ABI-binary64-clock;not-full-process-environment\n";
    return bytes;
}

void RStructuralNativeControl::stampBuiltinIdentity(RStructuralNativeControlSpec &value)
{
    value.providerId=QStringLiteral("range-fea.native.pi");
    value.providerRevision=QStringLiteral("1");
    value.implementationSha256=QString::fromLatin1(RANGE_NATIVE_CONTROL_IMPLEMENTATION_SHA256);
    value.dependencyEnvironmentManifest=builtinEnvironmentManifest();
    value.dependencyEnvironmentSha256=sha(value.dependencyEnvironmentManifest);
}

RStructuralNativeControl::RStructuralNativeControl(const RStructuralTrussDynamicsModelAdapter &value)
    : adapter(value),spec(value.nativeControlSpec())
{
    spec.validate();
    RStructuralNativeControlSpec identity; stampBuiltinIdentity(identity);
    if (spec.source!=RStructuralNativeControlSource::LocalV3 || spec.providerId!=identity.providerId ||
        spec.providerRevision!=identity.providerRevision || spec.implementationSha256!=identity.implementationSha256 ||
        spec.dependencyEnvironmentManifest!=identity.dependencyEnvironmentManifest ||
        spec.dependencyEnvironmentSha256!=identity.dependencyEnvironmentSha256)
        fail("Native LocalV3 provider/code/declared environment identity differs from the compiled factory.");
    if (spec.controlEveryFESteps!=1 || spec.samplePhase!=0 || spec.deterministicSeed!=0)
        fail("Native built-in PI supports only N=1, phase zero and deterministic seed zero.");
    const auto &maps=adapter.maps(); const auto &system=adapter.system();
    const auto &inputs=adapter.initialCompleteInputs();
    std::set<QString> sensorNames,actionNames,sensorOwners,actionOwners;
    for (const auto &sensor:spec.sensors)
    {
        if (sensor.sensorId.isEmpty() || !sensorNames.insert(sensor.sensorId).second)
            fail("Native control sensor identities must be unique and nonempty.");
        std::size_t target=0; QString owner;
        switch (sensor.quantity)
        {
            case RStructuralNativeMeasurementQuantity::NodeDisplacement:
            case RStructuralNativeMeasurementQuantity::NodeVelocity:
            case RStructuralNativeMeasurementQuantity::NodeAcceleration:
            case RStructuralNativeMeasurementQuantity::SupportReaction:
                if (!maps.nodeByPhysicalId.count(sensor.physicalNodeId) || sensor.axis>2 ||
                    sensor.physicalMemberId!=0 || !sensor.portId.isEmpty())
                    fail("Native node sensor requires exactly one physical node and a global XYZ axis.");
                target=3*maps.nodeByPhysicalId.at(sensor.physicalNodeId)+sensor.axis;
                if (sensor.quantity==RStructuralNativeMeasurementQuantity::SupportReaction && !system.prescribedMask.at(target))
                    fail("Native support reaction sensor cannot address a free residual DOF.");
                owner=QString("%1:n:%2:%3").arg(uint(sensor.quantity)).arg(sensor.physicalNodeId).arg(sensor.axis); break;
            case RStructuralNativeMeasurementQuantity::MemberAxialForce:
                if (!maps.memberByPhysicalId.count(sensor.physicalMemberId) || sensor.physicalNodeId!=0 ||
                    !sensor.portId.isEmpty() || sensor.axis!=0)
                    fail("Native member force sensor requires exactly one physical member.");
                target=maps.memberByPhysicalId.at(sensor.physicalMemberId);
                owner=QString("m:%1").arg(sensor.physicalMemberId); break;
            case RStructuralNativeMeasurementQuantity::PortVoltage:
            case RStructuralNativeMeasurementQuantity::PortCharge:
                if (!maps.memberByPortId.count(sensor.portId) || sensor.physicalNodeId!=0 || sensor.physicalMemberId!=0 || sensor.axis!=0)
                    fail("Native electrical sensor requires exactly one independent physical port ID.");
                target=maps.memberByPortId.at(sensor.portId);
                owner=QString("%1:p:%2").arg(uint(sensor.quantity)).arg(sensor.portId); break;
            default: fail("Native sensor quantity is unsupported.");
        }
        if (!sensorOwners.insert(owner).second) fail("Native sensor aliases duplicate the same quantity/physical target.");
        sensorTargets.push_back(target);
    }
    for (const auto &action:spec.actions)
    {
        if (action.actionId.isEmpty() || !actionNames.insert(action.actionId).second || action.delayTicks!=0 ||
            action.holdTicks!=1 || action.noiseAmplitude!=0.0)
            fail("Native actions require unique IDs and delay=0, hold=1, noise=0.");
        std::size_t target=0; QString owner;
        switch (action.quantity)
        {
            case RStructuralNativeActionQuantity::Force:
                if (!maps.nodeByPhysicalId.count(action.physicalNodeId) || action.axis>2 ||
                    action.physicalMemberId!=0 || !action.portId.isEmpty() || action.initialCommand!=0.0)
                    fail("Native force action requires one global node DOF and an initial zero increment.");
                target=3*maps.nodeByPhysicalId.at(action.physicalNodeId)+action.axis;
                owner=QString("n:%1:%2").arg(action.physicalNodeId).arg(action.axis); break;
            case RStructuralNativeActionQuantity::Voltage:
            case RStructuralNativeActionQuantity::Charge:
                if (!maps.memberByPortId.count(action.portId) || action.physicalNodeId!=0 || action.physicalMemberId!=0 || action.axis!=0)
                    fail("Native electrical action requires exactly one physical port.");
                target=maps.memberByPortId.at(action.portId); owner=QString("m:%1").arg(qulonglong(target));
                if (action.quantity==RStructuralNativeActionQuantity::Voltage)
                {
                    const auto *input=std::get_if<RStructuralTrussVoltage>(&inputs.memberInputs.at(target));
                    if (!input || !sameBits(action.initialCommand,input->value)) fail("Native voltage action mode/initial value differs from its source port.");
                }
                else
                {
                    const auto *input=std::get_if<RStructuralTrussCharge>(&inputs.memberInputs.at(target));
                    if (!input || !sameBits(action.initialCommand,input->value)) fail("Native charge action mode/initial value differs from its source port.");
                }
                break;
            case RStructuralNativeActionQuantity::CableReferenceLength:
            {
                if (!maps.memberByPhysicalId.count(action.physicalMemberId) || action.physicalNodeId!=0 ||
                    !action.portId.isEmpty() || action.axis!=0)
                    fail("Native cable action requires exactly one physical member.");
                target=maps.memberByPhysicalId.at(action.physicalMemberId); owner=QString("m:%1").arg(qulonglong(target));
                const auto *input=std::get_if<RStructuralTrussCableLength>(&inputs.memberInputs.at(target));
                if (!input || !sameBits(action.initialCommand,input->value)) fail("Native cable action mode/initial value differs from its source cable.");
                break;
            }
            default: fail("Native action quantity is unsupported.");
        }
        if (!actionOwners.insert(owner).second) fail("Native action target has multiple owners.");
        if (!std::isfinite(action.initialCommand) || std::isnan(action.minimum) || std::isnan(action.maximum) ||
            action.minimum>action.maximum || action.initialCommand<action.minimum || action.initialCommand>action.maximum ||
            std::isnan(action.rateLimit) || action.rateLimit<0.0)
            fail("Native initial command or amplitude/rate bounds are invalid.");
        actionTargets.push_back(target); initialValues.push_back(action.initialCommand);
    }
    std::set<std::size_t> loopOwners;
    for (const auto &loop:spec.piLoops)
    {
        const auto a=std::find_if(spec.actions.begin(),spec.actions.end(),[&](const auto &x){return x.actionId==loop.actionId;});
        const auto s=std::find_if(spec.sensors.begin(),spec.sensors.end(),[&](const auto &x){return x.sensorId==loop.sensorId;});
        if (a==spec.actions.end() || s==spec.sensors.end()) fail("Native PI loop has an unresolved action or sensor.");
        const auto ai=std::size_t(a-spec.actions.begin());
        if (!loopOwners.insert(ai).second) fail("Native PI loops cannot share action ownership.");
        loopActions.push_back(ai); loopSensors.push_back(std::size_t(s-spec.sensors.begin()));
    }
    if (loopOwners.size()!=spec.actions.size() || spec.actions.empty()) fail("Native PI requires exactly one loop for every action.");
}

RStructuralTrussDynamicsInputs RStructuralNativeControl::completeInputs(const std::vector<double> &values) const
{
    if (values.size()!=spec.actions.size()) fail("Native applied action count differs from its frozen catalog.");
    auto inputs=adapter.initialCompleteInputs();
    for (std::size_t i=0;i<values.size();++i)
    {
        const auto &action=spec.actions[i]; const double value=values[i]; const auto target=actionTargets[i];
        if (!std::isfinite(value) || value<action.minimum || value>action.maximum)
            fail("Native applied action violates its finite amplitude interval.");
        switch (action.quantity)
        {
            case RStructuralNativeActionQuantity::Force:
                inputs.externalForces.at(target)=adapter.initialCompleteInputs().externalForces.at(target)+value;
                if (!std::isfinite(inputs.externalForces.at(target))) fail("Native force background plus increment is nonfinite.");
                break;
            case RStructuralNativeActionQuantity::Voltage: inputs.memberInputs.at(target)=RStructuralTrussVoltage{value}; break;
            case RStructuralNativeActionQuantity::Charge: inputs.memberInputs.at(target)=RStructuralTrussCharge{value}; break;
            case RStructuralNativeActionQuantity::CableReferenceLength:
                if (value<=0.0) fail("Native applied cable reference length must be positive.");
                inputs.memberInputs.at(target)=RStructuralTrussCableLength{value}; break;
            default: fail("Native applied action quantity is unsupported.");
        }
    }
    return inputs;
}

RStructuralNativeControlState RStructuralNativeControl::initialize() const
{
    RStructuralNativeControlState state;
    state.available=true; state.sourceFingerprint=adapter.sourceFingerprint(); state.providerRevision=spec.providerRevision;
    std::vector<double> integrals;
    for (const auto &loop:spec.piLoops) integrals.push_back(loop.initialIntegral);
    state.controllerState=encode(integrals);
    state.heldValues=initialValues; state.lastAppliedValues=initialValues;
    state.validate(spec);
    return state;
}

void RStructuralNativeControl::validateFrame(const RStructuralNativeControlFrame &frame,uint completed) const
{
    if (completed>=adapter.targetSteps() || frame.completedFE!=completed || frame.controlTick!=completed ||
        frame.fePhase!=0 || frame.samplePhase!=0 || !sameBits(frame.sampleTime,adapter.targetTime(completed)) ||
        !sameBits(frame.nextFETime,adapter.targetTime(completed+1)) || !sameBits(frame.nextControlTime,frame.nextFETime) ||
        !sameBits(frame.physicalDt,adapter.targetTime(completed+1)-adapter.targetTime(completed)) ||
        !(frame.physicalDt>0.0) || frame.samples.size()!=spec.sensors.size())
        fail("Native left sample frame has a different canonical clock/catalog or N=1 phase.");
    for (std::size_t i=0;i<spec.sensors.size();++i)
        if (frame.samples[i].sensorId!=spec.sensors[i].sensorId || frame.samples[i].quantity!=spec.sensors[i].quantity ||
            !std::isfinite(frame.samples[i].value)) fail("Native left sample has a different sensor identity/quantity or nonfinite value.");
}

RStructuralNativeControlFrame RStructuralNativeControl::sample(const RStructuralTrussDynamicsState &state,uint completed) const
{
    if (state.time!=adapter.targetTime(completed) || completed>=adapter.targetSteps())
        fail("Native sampling requires the canonical left accepted state before the target.");
    adapter.driver().refresh(state,state.appliedInputs); // validate only; sample the original left acceleration.
    const auto &system=adapter.system();
    const auto fresh=RStructuralTrussAssembly::evaluate(system.referenceNodes,state.displacements,system.members,state.appliedInputs.memberInputs);
    RStructuralNativeControlFrame frame;
    frame.completedFE=completed; frame.controlTick=completed; frame.sampleTime=state.time;
    frame.nextFETime=adapter.targetTime(completed+1); frame.nextControlTime=frame.nextFETime;
    frame.physicalDt=frame.nextFETime-frame.sampleTime;
    for (std::size_t i=0;i<spec.sensors.size();++i)
    {
        const auto &sensor=spec.sensors[i]; const auto target=sensorTargets[i]; double value=0.0;
        switch (sensor.quantity)
        {
            case RStructuralNativeMeasurementQuantity::NodeDisplacement: value=state.displacements.at(target); break;
            case RStructuralNativeMeasurementQuantity::NodeVelocity: value=state.velocities.at(target); break;
            case RStructuralNativeMeasurementQuantity::NodeAcceleration: value=state.accelerations.at(target); break;
            case RStructuralNativeMeasurementQuantity::MemberAxialForce: value=fresh.memberResponses.at(target).axialForce; break;
            case RStructuralNativeMeasurementQuantity::PortVoltage:
                // A voltage-driven V sample is an action-record echo, not independent deformation sensing.
                if (!fresh.memberResponses.at(target).electrical) fail("Native voltage sensor has no independent port response.");
                value=fresh.memberResponses.at(target).electrical->voltage; break;
            case RStructuralNativeMeasurementQuantity::PortCharge:
                if (!fresh.memberResponses.at(target).electrical) fail("Native charge sensor has no independent port response.");
                value=fresh.memberResponses.at(target).electrical->conductorCharge; break;
            case RStructuralNativeMeasurementQuantity::SupportReaction:
                value=system.nodalLumpedMass.at(target/3)*state.accelerations.at(target)+
                    (fresh.internalForces.at(target)-state.appliedInputs.externalForces.at(target)); break;
            default: fail("Native sensor quantity is unsupported.");
        }
        frame.samples.push_back({sensor.sensorId,sensor.quantity,value});
    }
    validateFrame(frame,completed);
    return frame;
}

RStructuralNativeControl::PureProposal RStructuralNativeControl::step(
    const RStructuralNativeControlFrame &frame,const QByteArray &opaque) const
{
    auto integrals=decode(opaque,spec.piLoops.size());
    PureProposal output; output.requested.resize(spec.actions.size());
    for (std::size_t i=0;i<spec.piLoops.size();++i)
    {
        const auto &loop=spec.piLoops[i];
        const double error=loop.referenceValue-frame.samples.at(loopSensors[i]).value;
        const double command=loop.bias+loop.proportionalGain*error+loop.integralGain*integrals[i];
        const double next=integrals[i]+adapter.nominalTimeStep()*error;
        if (!std::isfinite(error) || !std::isfinite(command) || !std::isfinite(next)) fail("Native PI proposal is nonfinite.");
        output.requested.at(loopActions[i])=command; integrals[i]=next;
    }
    output.nextState=encode(integrals);
    return output;
}

RStructuralNativeControlActionRecord RStructuralNativeControl::limit(
    std::size_t index,double previous,double requested,double h) const
{
    const auto &action=spec.actions.at(index);
    if (!std::isfinite(previous) || previous<action.minimum || previous>action.maximum ||
        !std::isfinite(requested) || !(h>0.0) || !std::isfinite(h)) fail("Native limiter received invalid previous/requested/clock values.");
    double low=action.minimum,high=action.maximum;
    double rateLow=-std::numeric_limits<double>::infinity(),rateHigh=std::numeric_limits<double>::infinity();
    if (std::isfinite(action.rateLimit))
    {
        const double delta=action.rateLimit*h;
        if (!std::isfinite(delta)) fail("Native physical rate interval is not representable.");
        rateLow=previous-delta; rateHigh=previous+delta;
        low=std::max(low,rateLow); high=std::min(high,rateHigh);
    }
    if (low>high) fail("Native amplitude/rate feasible interval is empty.");
    RStructuralNativeControlActionRecord record;
    record.actionId=action.actionId; record.quantity=action.quantity; record.previousApplied=previous;
    record.requested=requested; record.conditioned=requested;
    const double amplitude=std::max(action.minimum,std::min(action.maximum,requested));
    record.applied=std::max(low,std::min(high,amplitude));
    if (amplitude!=requested) record.limitFlags|=R_NATIVE_CONTROL_AMPLITUDE_LIMITED;
    if (record.applied!=amplitude) record.limitFlags|=R_NATIVE_CONTROL_RATE_LIMITED;
    return record;
}

void RStructuralNativeControl::validateAccepted(const RStructuralNativeControlState &state,uint completed) const
{
    state.validate(spec);
    if (!state.available || state.sourceFingerprint!=adapter.sourceFingerprint() || state.providerRevision!=spec.providerRevision ||
        state.completedFE!=completed || state.completedControlTicks!=completed || state.fePhase!=0 || state.samplePhase!=0 ||
        completed>adapter.targetSteps()) fail("Native complete control state differs from source/provider or accepted FE/control phase.");
    if (completed==0)
    {
        if (state.fingerprintBytes()!=initialize().fingerprintBytes()) fail("Native initial control state differs from frozen PI initialization.");
        return;
    }
    validateFrame(state.lastFrame,completed-1);
    const auto output=step(state.lastFrame,state.lastControllerInputState);
    if (output.nextState!=state.controllerState || state.lastActions.size()!=spec.actions.size() ||
        state.heldValues.size()!=spec.actions.size() || state.lastAppliedValues.size()!=spec.actions.size())
        fail("Native controller envelope/action catalog differs from the latest pure PI proposal.");
    if (completed==1 && state.lastControllerInputState!=initialize().controllerState)
        fail("Native first PI proposal did not consume the frozen initial envelope.");
    for (std::size_t i=0;i<spec.actions.size();++i)
    {
        const auto &record=state.lastActions[i];
        if (completed==1 && !sameBits(record.previousApplied,initialValues[i])) fail("Native first limit step has a different initial applied command.");
        const auto expected=limit(i,record.previousApplied,output.requested[i],state.lastFrame.physicalDt);
        if (!sameRecord(record,expected) || !sameBits(state.heldValues[i],expected.applied) ||
            !sameBits(state.lastAppliedValues[i],expected.applied)) fail("Native requested/conditioned/applied/limits differ from pure proposal reconstruction.");
    }
    completeInputs(state.lastAppliedValues);
}

RStructuralNativeControlProposal RStructuralNativeControl::propose(
    const RStructuralNativeControlState &accepted,const RStructuralTrussDynamicsState &physical,uint completed) const
{
    validateAccepted(accepted,completed);
    if (!sameInputs(physical.appliedInputs,completeInputs(accepted.lastAppliedValues)))
        fail("Native left physical inputs differ from its complete accepted control action mapping.");
    RStructuralNativeControlProposal proposal;
    proposal.state=accepted;
    auto &state=proposal.state;
    state.lastFrame=sample(physical,completed);
    const auto output=step(state.lastFrame,accepted.controllerState);
    state.lastControllerInputState=accepted.controllerState; state.controllerState=output.nextState;
    state.lastActions.clear(); state.heldValues.clear(); state.lastAppliedValues.clear();
    for (std::size_t i=0;i<spec.actions.size();++i)
    {
        const auto record=limit(i,accepted.lastAppliedValues.at(i),output.requested.at(i),state.lastFrame.physicalDt);
        state.lastActions.push_back(record); state.heldValues.push_back(record.applied); state.lastAppliedValues.push_back(record.applied);
    }
    ++state.completedFE; ++state.completedControlTicks;
    proposal.completeInputs=completeInputs(state.lastAppliedValues);
    validateAccepted(state,completed+1);
    return proposal;
}

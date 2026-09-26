#include "rstructural_typed_control.h"

#include <QCryptographicHash>
#include <QSet>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>
#include <utility>

namespace
{
constexpr qsizetype maxStateBytes = 1024*1024;
constexpr uint maxChannels = 4096;
constexpr uint maxDelaySlots = 16384;
constexpr quint64 maxCompletedSteps = quint64(1)<<32;
constexpr quint32 stateVersion = 1;

bool validQuantity(RStructuralTypedControlQuantity quantity)
{
    switch (quantity)
    {
        case RStructuralTypedControlQuantity::Force:
        case RStructuralTypedControlQuantity::CableReferenceLength:
        case RStructuralTypedControlQuantity::CableInitialTension:
        case RStructuralTypedControlQuantity::CableAnchorDisplacement:
        case RStructuralTypedControlQuantity::Voltage:
        case RStructuralTypedControlQuantity::Charge:
            return true;
    }
    return false;
}

bool validChannelId(const QString &id)
{
    if (id.isEmpty() || id.size() > 128 || id.trimmed() != id) return false;
    for (const QChar character : id)
    {
        if (character.isSpace() || character.isNull() || character == QLatin1Char('"'))
            return false;
    }
    return true;
}

void appendU8(QByteArray &out, quint8 value)
{
    out.append(char(value));
}

void appendU32(QByteArray &out, quint32 value)
{
    for (int i=0;i<4;i++) out.append(char((value >> (8*i)) & 0xffu));
}

void appendU64(QByteArray &out, quint64 value)
{
    for (int i=0;i<8;i++) out.append(char((value >> (8*i)) & 0xffu));
}

void appendDouble(QByteArray &out, double value)
{
    static_assert(sizeof(double) == sizeof(quint64),"Typed control state requires 64-bit doubles");
    static_assert(std::numeric_limits<double>::is_iec559,
                  "Typed control state requires IEEE-754 binary64 doubles");
    quint64 bits = 0;
    std::memcpy(&bits,&value,sizeof(bits));
    appendU64(out,bits);
}

class Reader
{
    public:
        explicit Reader(const QByteArray &bytes) : bytes(bytes) {}

        QByteArray take(qsizetype count)
        {
            if (count < 0 || this->position > this->bytes.size() ||
                count > this->bytes.size()-this->position)
                throw std::invalid_argument("Truncated typed structural control state");
            const QByteArray value = this->bytes.mid(this->position,count);
            this->position += count;
            return value;
        }

        quint8 u8()
        {
            return quint8(this->take(1)[0]);
        }

        quint32 u32()
        {
            const QByteArray field = this->take(4);
            quint32 value = 0;
            for (int i=0;i<4;i++) value |= quint32(quint8(field[i])) << (8*i);
            return value;
        }

        quint64 u64()
        {
            const QByteArray field = this->take(8);
            quint64 value = 0;
            for (int i=0;i<8;i++) value |= quint64(quint8(field[i])) << (8*i);
            return value;
        }

        double floating()
        {
            const quint64 bits = this->u64();
            double value = 0.0;
            std::memcpy(&value,&bits,sizeof(value));
            return value;
        }

        void finish() const
        {
            if (this->position != this->bytes.size())
                throw std::invalid_argument("Trailing bytes in typed structural control state");
        }

    private:
        const QByteArray &bytes;
        qsizetype position = 0;
};

quint64 splitMix64(quint64 value)
{
    value += 0x9e3779b97f4a7c15ull;
    value = (value ^ (value >> 30)) * 0xbf58476d1ce4e5b9ull;
    value = (value ^ (value >> 27)) * 0x94d049bb133111ebull;
    return value ^ (value >> 31);
}

quint64 stableIdHash(const QString &id)
{
    quint64 hash = 0xcbf29ce484222325ull;
    const QByteArray utf8 = id.toUtf8();
    for (const char byte : utf8)
    {
        hash ^= quint8(byte);
        hash *= 0x100000001b3ull;
    }
    return hash;
}

double uniformNoise(quint64 seed, const QString &id, quint64 step, double amplitude)
{
    if (amplitude == 0.0) return 0.0;
    const quint64 bits = splitMix64(seed ^ stableIdHash(id) ^ splitMix64(step));
    const double unit = double(bits >> 11)*(1.0/9007199254740992.0);
    return amplitude*(2.0*unit-1.0);
}

QByteArray hash(const QByteArray &bytes)
{
    return QCryptographicHash::hash(bytes,QCryptographicHash::Sha256);
}

void validateSensors(const std::vector<RStructuralControlSensorSample> &sensors)
{
    QSet<QString> ids;
    for (const auto &sensor : sensors)
    {
        if (!validChannelId(sensor.sensorId) || ids.contains(sensor.sensorId) ||
            sensor.unit.size() > 32 || !std::isfinite(sensor.value))
            throw std::invalid_argument("Invalid typed control observation sensor");
        for (const QChar character : sensor.unit)
            if (character.isNull())
                throw std::invalid_argument("Invalid typed control sensor unit");
        ids.insert(sensor.sensorId);
    }
}

bool timeMatchesStepCount(double initialTime,
                          double acceptedTime,
                          quint64 completedSteps,
                          double fixedDt)
{
    const double expected=std::fma(double(completedSteps),fixedDt,initialTime);
    if (!std::isfinite(expected) || acceptedTime!=expected) return false;
    return completedSteps==0 ||
        std::fma(double(completedSteps-1),fixedDt,initialTime)<expected;
}

void requireFiniteSample(const RStructuralTypedControlChannelSample &sample)
{
    if (!std::isfinite(sample.requestedValue) || !std::isfinite(sample.conditionedValue) ||
        !std::isfinite(sample.appliedValue) || !std::isfinite(sample.cumulativeWorkJ) ||
        !std::isfinite(sample.averagePowerW))
        throw std::invalid_argument("Non-finite typed structural control candidate sample");
}
}

QString rStructuralTypedControlUnit(RStructuralTypedControlQuantity quantity)
{
    switch (quantity)
    {
        case RStructuralTypedControlQuantity::Force: return QStringLiteral("N");
        case RStructuralTypedControlQuantity::CableReferenceLength: return QStringLiteral("m");
        case RStructuralTypedControlQuantity::CableInitialTension: return QStringLiteral("N");
        case RStructuralTypedControlQuantity::CableAnchorDisplacement: return QStringLiteral("m");
        case RStructuralTypedControlQuantity::Voltage: return QStringLiteral("V");
        case RStructuralTypedControlQuantity::Charge: return QStringLiteral("C");
    }
    throw std::invalid_argument("Unknown typed structural control quantity");
}

RStructuralTypedControlSession::RStructuralTypedControlSession(
    quint64 deterministicSeed,
    double fixedDt,
    std::vector<RStructuralTypedControlChannel> catalog,
    QString controllerIdentity)
    : seed(deterministicSeed), timeStep(fixedDt), channels(std::move(catalog)),
      providerIdentity(std::move(controllerIdentity))
{
    if (!std::isfinite(this->timeStep) || this->timeStep <= 0.0 ||
        this->channels.empty() || this->channels.size() > maxChannels ||
        !validChannelId(this->providerIdentity))
        throw std::invalid_argument("Invalid typed structural control step size or catalog size");

    QSet<QString> ids;
    uint delaySlots = 0;
    QByteArray identityBytes("RSTC-config-v1",14);
    appendU64(identityBytes,this->seed);
    appendDouble(identityBytes,this->timeStep);
    const QByteArray controllerUtf8 = this->providerIdentity.toUtf8();
    appendU32(identityBytes,quint32(controllerUtf8.size()));
    identityBytes.append(controllerUtf8);
    appendU32(identityBytes,quint32(this->channels.size()));
    for (const auto &channel : this->channels)
    {
        const bool validBounds = !std::isnan(channel.minimum) && !std::isnan(channel.maximum) &&
                                 channel.minimum <= channel.maximum;
        const bool validRate = !std::isnan(channel.rateLimit) && channel.rateLimit >= 0.0;
        if (!validChannelId(channel.channelId) || ids.contains(channel.channelId) ||
            !validQuantity(channel.quantity) || !validBounds || !validRate ||
            !std::isfinite(channel.initialCommand) ||
            channel.initialCommand < channel.minimum || channel.initialCommand > channel.maximum ||
            channel.holdSteps == 0 || channel.delaySteps > maxDelaySlots ||
            !std::isfinite(channel.noiseAmplitude) || channel.noiseAmplitude < 0.0 ||
            delaySlots > maxDelaySlots-channel.delaySteps)
            throw std::invalid_argument("Invalid typed structural control channel catalog");
        if (std::isfinite(channel.rateLimit) &&
            !std::isfinite(channel.rateLimit*this->timeStep))
            throw std::invalid_argument("Typed structural control rate limit overflows this time step");

        ids.insert(channel.channelId);
        delaySlots += channel.delaySteps;
        const QByteArray id = channel.channelId.toUtf8();
        appendU32(identityBytes,quint32(id.size()));
        identityBytes.append(id);
        appendU8(identityBytes,quint8(channel.quantity));
        appendDouble(identityBytes,channel.minimum);
        appendDouble(identityBytes,channel.maximum);
        appendDouble(identityBytes,channel.rateLimit);
        appendDouble(identityBytes,channel.initialCommand);
        appendU32(identityBytes,channel.delaySteps);
        appendU32(identityBytes,channel.holdSteps);
        appendDouble(identityBytes,channel.noiseAmplitude);
    }
    this->digest = hash(identityBytes);
}

RStructuralTypedControlSession::State RStructuralTypedControlSession::initialState(
    double initialTime,
    QByteArray controllerState) const
{
    if (!std::isfinite(initialTime) || controllerState.size() > maxStateBytes)
        throw std::invalid_argument("Invalid typed structural control initial state");

    State state;
    state.initialTime = initialTime;
    state.acceptedTime = initialTime;
    state.configurationSha256 = this->digest;
    state.controllerState = std::move(controllerState);
    state.lastAppliedValues.reserve(this->channels.size());
    state.heldValues.reserve(this->channels.size());
    state.delayedValues.reserve(this->channels.size());
    state.cumulativeWorkJ.resize(this->channels.size(),0.0);
    for (const auto &channel : this->channels)
    {
        state.lastAppliedValues.push_back(channel.initialCommand);
        state.heldValues.push_back(channel.initialCommand);
        state.delayedValues.emplace_back(channel.delaySteps,channel.initialCommand);
    }
    this->validateState(state);
    return state;
}

RStructuralControlObservation RStructuralTypedControlSession::makeObservation(
    const State &accepted,
    std::vector<RStructuralControlSensorSample> sensors) const
{
    this->validateState(accepted);
    validateSensors(sensors);
    RStructuralControlObservation observation;
    observation.acceptedTime = accepted.acceptedTime;
    if (accepted.completedSteps>=maxCompletedSteps)
        throw std::invalid_argument("Typed structural control step counter is exhausted");
    observation.nextTime = std::fma(double(accepted.completedSteps+1),this->timeStep,accepted.initialTime);
    if (!std::isfinite(observation.nextTime) || observation.nextTime <= accepted.acceptedTime ||
        std::fabs((observation.nextTime-accepted.acceptedTime)-this->timeStep)>1.0e-9*this->timeStep)
        throw std::invalid_argument("Typed structural control time step is not representable at this time");
    observation.deterministicSeed = this->seed;
    observation.sensors = std::move(sensors);
    observation.priorControllerState = accepted.controllerState;
    return observation;
}

RStructuralTypedControlSession::Candidate RStructuralTypedControlSession::advance(
    const State &accepted,
    const RStructuralControlObservation &observation,
    const RStructuralTypedControlDecision &decision) const
{
    this->validateState(accepted);
    validateSensors(observation.sensors);
    const double expectedNextTime = std::fma(double(accepted.completedSteps+1),this->timeStep,accepted.initialTime);
    if (!std::isfinite(expectedNextTime) || expectedNextTime <= accepted.acceptedTime ||
        std::fabs((expectedNextTime-accepted.acceptedTime)-this->timeStep)>1.0e-9*this->timeStep ||
        observation.acceptedTime != accepted.acceptedTime ||
        observation.nextTime != expectedNextTime ||
        observation.deterministicSeed != this->seed ||
        observation.priorControllerState != accepted.controllerState)
        throw std::invalid_argument("Typed structural control observation does not match accepted state");
    if (accepted.completedSteps >= maxCompletedSteps ||
        decision.nextControllerState.size() > maxStateBytes ||
        decision.requestedCommands.size() != this->channels.size())
        throw std::invalid_argument("Invalid typed structural control decision size or step counter");

    for (std::size_t i=0;i<this->channels.size();i++)
    {
        const auto &expected = this->channels[i];
        const auto &command = decision.requestedCommands[i];
        if (command.channelId != expected.channelId || command.quantity != expected.quantity ||
            !std::isfinite(command.value))
            throw std::invalid_argument("Typed control command ID, quantity, order, or value mismatch");
    }

    Candidate result;
    result.state = accepted;
    result.state.completedSteps++;
    result.state.acceptedTime = observation.nextTime;
    result.state.controllerState = decision.nextControllerState;
    result.channels.reserve(this->channels.size());

    for (std::size_t i=0;i<this->channels.size();i++)
    {
        const auto &channel = this->channels[i];
        const auto &requested = decision.requestedCommands[i];
        double held = accepted.heldValues[i];
        if (accepted.completedSteps % channel.holdSteps == 0)
        {
            held = requested.value+uniformNoise(this->seed,channel.channelId,
                                                accepted.completedSteps,channel.noiseAmplitude);
            if (!std::isfinite(held))
                throw std::invalid_argument("Typed control held command overflow");
        }

        double conditioned = held;
        if (channel.delaySteps != 0)
        {
            const std::size_t slot = std::size_t(accepted.completedSteps % channel.delaySteps);
            conditioned = accepted.delayedValues[i][slot];
            result.state.delayedValues[i][slot] = held;
        }
        if (!std::isfinite(conditioned))
            throw std::invalid_argument("Typed control conditioned command is non-finite");

        double lower = channel.minimum;
        double upper = channel.maximum;
        if (std::isfinite(channel.rateLimit))
        {
            const double rateDelta = channel.rateLimit*this->timeStep;
            const long double rateLower = static_cast<long double>(accepted.lastAppliedValues[i])-
                                          static_cast<long double>(rateDelta);
            const long double rateUpper = static_cast<long double>(accepted.lastAppliedValues[i])+
                                          static_cast<long double>(rateDelta);
            const long double feasibleLower = std::max(
                static_cast<long double>(channel.minimum),rateLower);
            const long double feasibleUpper = std::min(
                static_cast<long double>(channel.maximum),rateUpper);
            if (feasibleLower > feasibleUpper)
                throw std::invalid_argument("Typed control amplitude/rate feasible interval is empty");
            const long double maxDouble = std::numeric_limits<double>::max();
            lower = feasibleLower < -maxDouble ? -std::numeric_limits<double>::infinity()
                                               : double(feasibleLower);
            upper = feasibleUpper > maxDouble ? std::numeric_limits<double>::infinity()
                                              : double(feasibleUpper);
        }
        if (lower > upper)
            throw std::invalid_argument("Typed control amplitude/rate feasible interval is empty");

        double applied = conditioned;
        if (applied < lower) applied = lower;
        if (applied > upper) applied = upper;
        if (!std::isfinite(applied))
            throw std::invalid_argument("Typed control applied value is non-finite");

        result.state.heldValues[i] = held;
        result.state.lastAppliedValues[i] = applied;
        RStructuralTypedControlChannelSample sample;
        sample.channelId = channel.channelId;
        sample.quantity = channel.quantity;
        sample.requestedValue = requested.value;
        sample.conditionedValue = conditioned;
        sample.appliedValue = applied;
        sample.cumulativeWorkJ = result.state.cumulativeWorkJ[i];
        sample.averagePowerW = 0.0;
        result.channels.push_back(std::move(sample));
    }

    this->validateState(result.state);
    this->validateCandidate(result);
    return result;
}

RStructuralTypedControlSession::Candidate RStructuralTypedControlSession::trial(
    const State &accepted,
    std::vector<RStructuralControlSensorSample> sensors,
    const ControllerCallback &controller) const
{
    if (!controller) throw std::invalid_argument("Typed structural control callback is empty");
    const auto observation = this->makeObservation(accepted,std::move(sensors));
    const auto decision = controller(observation);
    return this->advance(accepted,observation,decision);
}

RStructuralTypedControlSession::Candidate RStructuralTypedControlSession::withWork(
    const Candidate &candidate,
    const std::vector<double> &workIncrementsJ) const
{
    this->validateCandidate(candidate);
    if (workIncrementsJ.size() != this->channels.size())
        throw std::invalid_argument("Typed control work increment count does not match catalog");

    Candidate result = candidate;
    for (std::size_t i=0;i<this->channels.size();i++)
    {
        const double increment = workIncrementsJ[i];
        if (!std::isfinite(increment))
            throw std::invalid_argument("Typed control work increment is non-finite");
        result.state.cumulativeWorkJ[i] += increment;
        result.channels[i].cumulativeWorkJ = result.state.cumulativeWorkJ[i];
        result.channels[i].averagePowerW = increment/this->timeStep;
        if (!std::isfinite(result.channels[i].cumulativeWorkJ) ||
            !std::isfinite(result.channels[i].averagePowerW))
            throw std::invalid_argument("Typed control work or average power overflow");
    }
    this->validateState(result.state);
    this->validateCandidate(result);
    return result;
}

QByteArray RStructuralTypedControlSession::serialize(const State &state) const
{
    this->validateState(state);
    QByteArray result("RTCS",4);
    appendU32(result,stateVersion);
    result.append(this->digest);
    appendDouble(result,state.initialTime);
    appendDouble(result,state.acceptedTime);
    appendU64(result,state.completedSteps);
    appendU32(result,quint32(state.controllerState.size()));
    result.append(state.controllerState);
    appendU32(result,quint32(this->channels.size()));
    for (std::size_t i=0;i<this->channels.size();i++)
    {
        appendDouble(result,state.lastAppliedValues[i]);
        appendDouble(result,state.heldValues[i]);
        appendDouble(result,state.cumulativeWorkJ[i]);
        appendU32(result,quint32(state.delayedValues[i].size()));
        for (double value : state.delayedValues[i]) appendDouble(result,value);
    }
    if (result.size()+32 > maxStateBytes)
        throw std::invalid_argument("Typed structural control serialized state exceeds 1 MiB");
    result.append(hash(result));
    return result;
}

RStructuralTypedControlSession::State RStructuralTypedControlSession::deserialize(
    const QByteArray &bytes) const
{
    if (bytes.size() > maxStateBytes || bytes.size() < 4+4+32+8+8+8+4+4+32)
        throw std::invalid_argument("Invalid typed structural control state size");
    const qsizetype payloadSize = bytes.size()-32;
    const QByteArray payload = bytes.left(payloadSize);
    if (bytes.mid(payloadSize) != hash(payload))
        throw std::invalid_argument("Typed structural control state checksum mismatch");

    Reader input(payload);
    if (input.take(4) != QByteArray("RTCS",4) || input.u32() != stateVersion ||
        input.take(32) != this->digest)
        throw std::invalid_argument("Typed structural control state format/configuration mismatch");

    State state;
    state.configurationSha256 = this->digest;
    state.initialTime = input.floating();
    state.acceptedTime = input.floating();
    state.completedSteps = input.u64();
    const quint32 controllerSize = input.u32();
    if (controllerSize > maxStateBytes)
        throw std::invalid_argument("Typed controller state exceeds 1 MiB");
    state.controllerState = input.take(controllerSize);

    const quint32 count = input.u32();
    if (count != this->channels.size())
        throw std::invalid_argument("Typed structural control state channel count mismatch");
    state.lastAppliedValues.reserve(count);
    state.heldValues.reserve(count);
    state.delayedValues.reserve(count);
    state.cumulativeWorkJ.reserve(count);
    for (quint32 i=0;i<count;i++)
    {
        state.lastAppliedValues.push_back(input.floating());
        state.heldValues.push_back(input.floating());
        state.cumulativeWorkJ.push_back(input.floating());
        const quint32 delayCount = input.u32();
        if (delayCount != this->channels[i].delaySteps)
            throw std::invalid_argument("Typed structural control delay state mismatch");
        auto &queue = state.delayedValues.emplace_back();
        queue.reserve(delayCount);
        for (quint32 j=0;j<delayCount;j++) queue.push_back(input.floating());
    }
    input.finish();
    this->validateState(state);
    return state;
}

void RStructuralTypedControlSession::validateState(const State &state) const
{
    if (state.configurationSha256 != this->digest || !std::isfinite(state.initialTime) ||
        !std::isfinite(state.acceptedTime) || state.acceptedTime < state.initialTime ||
        state.completedSteps > maxCompletedSteps ||
        !timeMatchesStepCount(state.initialTime,state.acceptedTime,
                              state.completedSteps,this->timeStep) ||
        (state.completedSteps == 0 && state.acceptedTime != state.initialTime) ||
        (state.completedSteps != 0 && state.acceptedTime <= state.initialTime) ||
        state.controllerState.size() > maxStateBytes ||
        state.lastAppliedValues.size() != this->channels.size() ||
        state.heldValues.size() != this->channels.size() ||
        state.delayedValues.size() != this->channels.size() ||
        state.cumulativeWorkJ.size() != this->channels.size())
        throw std::invalid_argument("Typed structural control state/configuration mismatch");

    for (std::size_t i=0;i<this->channels.size();i++)
    {
        const auto &channel = this->channels[i];
        if (state.delayedValues[i].size() != channel.delaySteps ||
            !std::isfinite(state.lastAppliedValues[i]) ||
            state.lastAppliedValues[i] < channel.minimum || state.lastAppliedValues[i] > channel.maximum ||
            !std::isfinite(state.heldValues[i]) || !std::isfinite(state.cumulativeWorkJ[i]))
            throw std::invalid_argument("Invalid typed structural control channel state");
        for (const double value : state.delayedValues[i])
        {
            if (!std::isfinite(value))
                throw std::invalid_argument("Non-finite delayed typed control value");
        }
    }
}

void RStructuralTypedControlSession::validateCandidate(const Candidate &candidate) const
{
    this->validateState(candidate.state);
    if (candidate.channels.size() != this->channels.size())
        throw std::invalid_argument("Typed structural control candidate channel count mismatch");
    for (std::size_t i=0;i<this->channels.size();i++)
    {
        const auto &channel = this->channels[i];
        const auto &sample = candidate.channels[i];
        if (sample.channelId != channel.channelId || sample.quantity != channel.quantity ||
            sample.appliedValue != candidate.state.lastAppliedValues[i] ||
            sample.cumulativeWorkJ != candidate.state.cumulativeWorkJ[i])
            throw std::invalid_argument("Typed structural control candidate/state mismatch");
        requireFiniteSample(sample);
    }
}

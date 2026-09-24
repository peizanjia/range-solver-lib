#include "rstructural_control_conditioner.h"

#include <QCryptographicHash>
#include <QJsonDocument>
#include <QJsonParseError>
#include <QSet>

#include <cmath>
#include <cstring>
#include <limits>
#include <stdexcept>

namespace
{
constexpr int maxStateBytes = 1048576;
constexpr uint maxChannels = 4096;
constexpr uint maxDelaySlots = 16384;

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
    static_assert(sizeof(double) == sizeof(quint64));
    quint64 bits = 0;
    std::memcpy(&bits,&value,sizeof(bits));
    appendU64(out,bits);
}

struct Reader
{
    const QByteArray &bytes;
    qsizetype position = 0;

    QByteArray take(qsizetype count)
    {
        if (count < 0 || position > bytes.size() || count > bytes.size()-position)
            throw std::invalid_argument("Truncated structural control state");
        QByteArray result = bytes.mid(position,count);
        position += count;
        return result;
    }

    quint32 u32()
    {
        const QByteArray field = take(4);
        quint32 value = 0;
        for (int i=0;i<4;i++) value |= quint32(quint8(field[i])) << (8*i);
        return value;
    }

    quint64 u64()
    {
        const QByteArray field = take(8);
        quint64 value = 0;
        for (int i=0;i<8;i++) value |= quint64(quint8(field[i])) << (8*i);
        return value;
    }

    double floating()
    {
        const quint64 bits = u64();
        double value = 0.0;
        std::memcpy(&value,&bits,sizeof(value));
        return value;
    }

    void finish() const
    {
        if (position != bytes.size())
            throw std::invalid_argument("Trailing bytes in structural control state");
    }
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
    const quint64 bits = splitMix64(seed ^ stableIdHash(id) ^
                                    splitMix64(step));
    const double unit = double(bits >> 11) * (1.0/9007199254740992.0);
    return amplitude*(2.0*unit-1.0);
}

QByteArray canonicalControllerJson(const QByteArray &bytes)
{
    QJsonParseError error;
    const QJsonDocument parsed = QJsonDocument::fromJson(bytes.isEmpty() ? QByteArray("{}") : bytes,&error);
    if (error.error != QJsonParseError::NoError || !parsed.isObject())
        throw std::invalid_argument("Controller-private state must be a JSON object");
    return parsed.toJson(QJsonDocument::Compact);
}
}

RStructuralControlConditioner::RStructuralControlConditioner(
    quint64 deterministicSeed,
    std::vector<RStructuralControlConditioningPolicy> policies)
    : seed(deterministicSeed), policies(std::move(policies))
{
    if (this->policies.empty() || this->policies.size() > maxChannels)
        throw std::invalid_argument("Structural control conditioner needs 1..4096 channels");
    QSet<QString> ids;
    uint delaySlots = 0;
    QByteArray identity("RCC-config-v1",13);
    appendU64(identity,this->seed);
    appendU32(identity,quint32(this->policies.size()));
    for (const auto &policy : this->policies)
    {
        if (policy.actuatorId.isEmpty() || policy.actuatorId.size() > 128 ||
            policy.actuatorId.trimmed() != policy.actuatorId ||
            ids.contains(policy.actuatorId) || policy.holdSteps == 0 ||
            policy.delaySteps > maxDelaySlots ||
            !std::isfinite(policy.noiseAmplitude) || policy.noiseAmplitude < 0.0 ||
            delaySlots > maxDelaySlots-policy.delaySteps)
            throw std::invalid_argument("Invalid structural control conditioning policy");
        for (const QChar character : policy.actuatorId)
        {
            if (character.isSpace() || character.isNull() || character == QLatin1Char('"'))
                throw std::invalid_argument("Invalid structural control conditioning channel ID");
        }
        ids.insert(policy.actuatorId);
        delaySlots += policy.delaySteps;
        const QByteArray id = policy.actuatorId.toUtf8();
        appendU32(identity,quint32(id.size()));
        identity.append(id);
        appendU32(identity,policy.delaySteps);
        appendU32(identity,policy.holdSteps);
        appendDouble(identity,policy.noiseAmplitude);
    }
    this->digest = QCryptographicHash::hash(identity,QCryptographicHash::Sha256);
}

RStructuralControlConditioner::State RStructuralControlConditioner::initialState() const
{
    State result;
    result.configurationSha256 = this->digest;
    result.heldCommands.resize(this->policies.size(),0.0);
    result.delayedCommands.reserve(this->policies.size());
    for (const auto &policy : this->policies)
        result.delayedCommands.emplace_back(policy.delaySteps,0.0);
    return result;
}

void RStructuralControlConditioner::validate(const State &state) const
{
    if (state.configurationSha256 != this->digest ||
        state.heldCommands.size() != this->policies.size() ||
        state.delayedCommands.size() != this->policies.size())
        throw std::invalid_argument("Structural control conditioner state/configuration mismatch");
    for (size_t i=0;i<this->policies.size();i++)
    {
        if (state.delayedCommands[i].size() != this->policies[i].delaySteps ||
            !std::isfinite(state.heldCommands[i]))
            throw std::invalid_argument("Invalid structural control conditioner channel state");
        for (const double value : state.delayedCommands[i])
            if (!std::isfinite(value))
                throw std::invalid_argument("Non-finite delayed structural control command");
    }
}

RStructuralControlConditioner::Result RStructuralControlConditioner::advance(
    const State &accepted,
    const std::vector<RStructuralControlCommand> &requested) const
{
    this->validate(accepted);
    if (accepted.completedSteps == std::numeric_limits<quint64>::max() ||
        requested.size() != this->policies.size())
        throw std::invalid_argument("Invalid structural control conditioner step");
    Result result;
    result.candidate = accepted;
    result.candidate.completedSteps++;
    result.appliedCommands.reserve(requested.size());
    for (size_t i=0;i<requested.size();i++)
    {
        const auto &policy = this->policies[i];
        if (requested[i].actuatorId != policy.actuatorId ||
            !std::isfinite(requested[i].value))
            throw std::invalid_argument("Structural control conditioner command order/value mismatch");
        double held = accepted.heldCommands[i];
        if (accepted.completedSteps % policy.holdSteps == 0)
        {
            held = requested[i].value + uniformNoise(this->seed,policy.actuatorId,
                                                       accepted.completedSteps,policy.noiseAmplitude);
            if (!std::isfinite(held))
                throw std::invalid_argument("Structural control conditioned command overflow");
        }
        result.candidate.heldCommands[i] = held;
        double applied = held;
        if (policy.delaySteps != 0)
        {
            const size_t slot = size_t(accepted.completedSteps % policy.delaySteps);
            applied = accepted.delayedCommands[i][slot];
            result.candidate.delayedCommands[i][slot] = held;
        }
        result.appliedCommands.push_back({policy.actuatorId,applied});
    }
    return result;
}

QByteArray RStructuralControlConditioner::serialize(const State &state) const
{
    this->validate(state);
    QByteArray result("RCC1",4);
    result.append(this->digest);
    appendU64(result,state.completedSteps);
    appendU32(result,quint32(this->policies.size()));
    for (size_t i=0;i<this->policies.size();i++)
    {
        appendDouble(result,state.heldCommands[i]);
        appendU32(result,quint32(state.delayedCommands[i].size()));
        for (const double value : state.delayedCommands[i]) appendDouble(result,value);
    }
    if (result.size() > maxStateBytes)
        throw std::invalid_argument("Structural control conditioner state exceeds 1 MiB");
    return result;
}

RStructuralControlConditioner::State RStructuralControlConditioner::deserialize(
    const QByteArray &bytes) const
{
    if (bytes.size() > maxStateBytes)
        throw std::invalid_argument("Structural control conditioner state exceeds 1 MiB");
    Reader input{bytes};
    if (input.take(4) != QByteArray("RCC1",4) || input.take(32) != this->digest)
        throw std::invalid_argument("Structural control conditioner format/configuration mismatch");
    State state;
    state.configurationSha256 = this->digest;
    state.completedSteps = input.u64();
    const uint count = input.u32();
    if (count != this->policies.size())
        throw std::invalid_argument("Structural control conditioner channel count mismatch");
    state.heldCommands.reserve(count);
    state.delayedCommands.reserve(count);
    for (uint i=0;i<count;i++)
    {
        state.heldCommands.push_back(input.floating());
        const uint delay = input.u32();
        if (delay != this->policies[i].delaySteps)
            throw std::invalid_argument("Structural control conditioner delay length mismatch");
        auto &queue = state.delayedCommands.emplace_back();
        queue.reserve(delay);
        for (uint j=0;j<delay;j++) queue.push_back(input.floating());
    }
    input.finish();
    this->validate(state);
    return state;
}

QByteArray RStructuralControlStateEnvelope::serialize() const
{
    const QByteArray json = canonicalControllerJson(this->controllerJson);
    if (json.size() > maxStateBytes || this->conditionerState.size() > maxStateBytes ||
        json.size() + this->conditionerState.size() > maxStateBytes-16)
        throw std::invalid_argument("Composite structural control state exceeds 1 MiB");
    QByteArray result("RCCE",4);
    appendU32(result,1);
    appendU32(result,quint32(json.size()));
    result.append(json);
    appendU32(result,quint32(this->conditionerState.size()));
    result.append(this->conditionerState);
    return result;
}

RStructuralControlStateEnvelope RStructuralControlStateEnvelope::deserialize(
    const QByteArray &bytes)
{
    if (bytes.size() > maxStateBytes)
        throw std::invalid_argument("Composite structural control state exceeds 1 MiB");
    Reader input{bytes};
    if (input.take(4) != QByteArray("RCCE",4) || input.u32() != 1)
        throw std::invalid_argument("Unsupported composite structural control state version");
    RStructuralControlStateEnvelope result;
    const uint jsonSize = input.u32();
    result.controllerJson = canonicalControllerJson(input.take(jsonSize));
    const uint conditionerSize = input.u32();
    result.conditionerState = input.take(conditionerSize);
    input.finish();
    return result;
}

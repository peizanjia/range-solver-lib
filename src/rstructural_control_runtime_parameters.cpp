#include "rstructural_control_runtime_parameters.h"

#include <QCryptographicHash>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonParseError>
#include <QJsonValue>
#include <QStringList>
#include <QUuid>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <utility>

namespace
{
using namespace RStructuralControlRuntimeParameters;

constexpr double minimumVelocityGainNsm = 1.0e-9;
constexpr double maximumVelocityGainNsm = 1.0e12;
constexpr int maximumActuatorIdLength = 128;
constexpr char stateSchema[] = "range.structural-runtime-parameters/v1";
constexpr char parameterHashMarker[] = "range.structural-runtime-parameters/hash/v1";

[[noreturn]] void invalid(const char *message)
{
    throw std::invalid_argument(message);
}

bool isAsciiAlphaNumeric(QChar character)
{
    const ushort value = character.unicode();
    return (value >= 'a' && value <= 'z') ||
           (value >= 'A' && value <= 'Z') ||
           (value >= '0' && value <= '9');
}

void validateActuatorId(const QString &id)
{
    if (id.isEmpty() || id.size() > maximumActuatorIdLength ||
        !isAsciiAlphaNumeric(id.front()))
        invalid("Invalid structural runtime actuator ID");
    for (const QChar character : id)
    {
        if (!isAsciiAlphaNumeric(character) && character != QLatin1Char('_') &&
            character != QLatin1Char('-') && character != QLatin1Char('.') &&
            character != QLatin1Char(':') && character != QLatin1Char('/'))
            invalid("Invalid structural runtime actuator ID");
    }
}

bool isLowerSha256(const QString &value)
{
    if (value.size() != 64) return false;
    for (const QChar character : value)
    {
        const ushort code = character.unicode();
        if (!((code >= '0' && code <= '9') || (code >= 'a' && code <= 'f')))
            return false;
    }
    return true;
}

bool isCanonicalUuid(const QString &value)
{
    if (value.size() != 36) return false;
    const QUuid uuid(value);
    return !uuid.isNull() && uuid.toString(QUuid::WithoutBraces) == value;
}

void validateLimit(const ActuatorLimit &limit)
{
    validateActuatorId(limit.actuatorId);
    if (!std::isfinite(limit.minN) || !std::isfinite(limit.maxN) ||
        limit.minN > limit.maxN || !std::isfinite(limit.rateNPerS) ||
        limit.rateNPerS < 0.0)
        invalid("Invalid structural runtime actuator limit");
}

double canonicalZero(double value)
{
    return value == 0.0 ? 0.0 : value;
}

std::vector<ActuatorLimit> sortedAndValidatedLimits(std::vector<ActuatorLimit> limits)
{
    if (limits.empty() || limits.size() > maximumActuatorCount)
        invalid("Structural runtime actuator count must be 1..256");
    for (auto &limit : limits)
    {
        validateLimit(limit);
        limit.minN = canonicalZero(limit.minN);
        limit.maxN = canonicalZero(limit.maxN);
        limit.rateNPerS = canonicalZero(limit.rateNPerS);
    }
    std::sort(limits.begin(), limits.end(), [](const auto &left, const auto &right) {
        return left.actuatorId < right.actuatorId;
    });
    for (std::size_t i = 1; i < limits.size(); ++i)
        if (limits[i - 1].actuatorId == limits[i].actuatorId)
            invalid("Duplicate structural runtime actuator ID");
    return limits;
}

void validateSortedLimits(const std::vector<ActuatorLimit> &limits)
{
    if (limits.empty() || limits.size() > maximumActuatorCount)
        invalid("Structural runtime actuator count must be 1..256");
    for (std::size_t i = 0; i < limits.size(); ++i)
    {
        validateLimit(limits[i]);
        if (i != 0 && !(limits[i - 1].actuatorId < limits[i].actuatorId))
            invalid("Structural runtime actuator IDs must be unique and sorted");
    }
}

bool sameActuatorIds(const std::vector<ActuatorLimit> &left,
                     const std::vector<ActuatorLimit> &right)
{
    if (left.size() != right.size()) return false;
    for (std::size_t i = 0; i < left.size(); ++i)
        if (left[i].actuatorId != right[i].actuatorId) return false;
    return true;
}

void appendU32BigEndian(QByteArray &bytes, quint32 value)
{
    bytes.append(char((value >> 24) & 0xffu));
    bytes.append(char((value >> 16) & 0xffu));
    bytes.append(char((value >> 8) & 0xffu));
    bytes.append(char(value & 0xffu));
}

void appendU64BigEndian(QByteArray &bytes, quint64 value)
{
    for (int shift = 56; shift >= 0; shift -= 8)
        bytes.append(char((value >> shift) & 0xffu));
}

void appendF64BigEndian(QByteArray &bytes, double value)
{
    static_assert(sizeof(double) == sizeof(quint64));
    static_assert(std::numeric_limits<double>::is_iec559,
                  "WF09 runtime parameter hashing requires IEEE-754 binary64");
    const double canonical = canonicalZero(value);
    quint64 bits = 0;
    std::memcpy(&bits, &canonical, sizeof(bits));
    appendU64BigEndian(bytes, bits);
}

void appendLengthDelimited(QByteArray &bytes, const QByteArray &field)
{
    if (quint64(field.size()) > std::numeric_limits<quint32>::max())
        invalid("Structural runtime revision field is too large");
    appendU32BigEndian(bytes, quint32(field.size()));
    bytes.append(field);
}

QJsonObject actuatorLimitObject(const ActuatorLimit &limit)
{
    return QJsonObject{
        {QStringLiteral("actuatorId"), limit.actuatorId},
        {QStringLiteral("minN"), canonicalZero(limit.minN)},
        {QStringLiteral("maxN"), canonicalZero(limit.maxN)},
        {QStringLiteral("rateNPerS"), canonicalZero(limit.rateNPerS)}};
}

QJsonObject semanticParametersObject(const Parameters &parameters)
{
    QJsonArray limits;
    for (const auto &limit : parameters.actuatorLimits)
        limits.append(actuatorLimitObject(limit));

    return QJsonObject{
        {QStringLiteral("revision"), QString::number(parameters.revision)},
        {QStringLiteral("velocityGainNsm"), canonicalZero(parameters.velocityGainNsm)},
        {QStringLiteral("actuatorLimits"), limits},
        {QStringLiteral("fixedBaseControllerRevision"), parameters.fixedBaseControllerRevision}};
}

QString sha256Hex(const QByteArray &bytes)
{
    return QString::fromLatin1(
        QCryptographicHash::hash(bytes, QCryptographicHash::Sha256).toHex());
}

QString computedParameterSha256(const Parameters &parameters)
{
    QByteArray identity(parameterHashMarker, qsizetype(sizeof(parameterHashMarker)));
    const QByteArray baseSha = QByteArray::fromHex(
        parameters.fixedBaseControllerRevision.toLatin1());
    if (baseSha.size() != 32 || parameters.actuatorLimits.size() > maximumActuatorCount)
        invalid("Invalid structural runtime parameter hash input");
    identity.append(baseSha);
    appendU64BigEndian(identity, parameters.revision);
    appendF64BigEndian(identity, parameters.velocityGainNsm);
    appendU32BigEndian(identity, quint32(parameters.actuatorLimits.size()));
    for (const auto &limit : parameters.actuatorLimits)
    {
        const QByteArray actuatorId = limit.actuatorId.toUtf8();
        if (quint64(actuatorId.size()) > std::numeric_limits<quint32>::max())
            invalid("Structural runtime actuator ID is too large to hash");
        appendU32BigEndian(identity, quint32(actuatorId.size()));
        identity.append(actuatorId);
        appendF64BigEndian(identity, limit.minN);
        appendF64BigEndian(identity, limit.maxN);
        appendF64BigEndian(identity, limit.rateNPerS);
    }
    return sha256Hex(identity);
}

QString computedRuntimeControllerRevision(const Parameters &parameters)
{
    const QByteArray base = parameters.fixedBaseControllerRevision.toLatin1();
    const QByteArray parameter = parameters.parameterSha256.toLatin1();
    QByteArray identity;
    appendLengthDelimited(identity, base);
    appendLengthDelimited(identity, parameter);
    return sha256Hex(identity);
}

void validateParameters(const Parameters &parameters)
{
    if (!std::isfinite(parameters.velocityGainNsm) ||
        parameters.velocityGainNsm < minimumVelocityGainNsm ||
        parameters.velocityGainNsm > maximumVelocityGainNsm)
        invalid("Structural runtime velocity gain must be finite and within [1e-9, 1e12] N*s/m");
    if (!isLowerSha256(parameters.fixedBaseControllerRevision))
        invalid("Fixed-base controller revision must be lowercase SHA-256 hex");
    validateSortedLimits(parameters.actuatorLimits);
    if (!isLowerSha256(parameters.parameterSha256) ||
        parameters.parameterSha256 != computedParameterSha256(parameters))
        invalid("Structural runtime parameter SHA-256 mismatch");
    if (!isLowerSha256(parameters.runtimeControllerRevision) ||
        parameters.runtimeControllerRevision != computedRuntimeControllerRevision(parameters))
        invalid("Structural runtime controller revision mismatch");
}

Parameters withComputedDigests(Parameters parameters)
{
    if (!std::isfinite(parameters.velocityGainNsm) ||
        parameters.velocityGainNsm < minimumVelocityGainNsm ||
        parameters.velocityGainNsm > maximumVelocityGainNsm)
        invalid("Structural runtime velocity gain must be finite and within [1e-9, 1e12] N*s/m");
    if (!isLowerSha256(parameters.fixedBaseControllerRevision))
        invalid("Fixed-base controller revision must be lowercase SHA-256 hex");
    validateSortedLimits(parameters.actuatorLimits);
    parameters.velocityGainNsm = canonicalZero(parameters.velocityGainNsm);
    parameters.parameterSha256 = computedParameterSha256(parameters);
    parameters.runtimeControllerRevision = computedRuntimeControllerRevision(parameters);
    return parameters;
}

QJsonObject parametersObject(const Parameters &parameters)
{
    QJsonObject result = semanticParametersObject(parameters);
    result.insert(QStringLiteral("parameterSha256"), parameters.parameterSha256);
    result.insert(QStringLiteral("runtimeControllerRevision"),
                  parameters.runtimeControllerRevision);
    return result;
}

QJsonObject changeObject(const Change &change)
{
    return QJsonObject{
        {QStringLiteral("proposalId"), change.proposalId},
        {QStringLiteral("previousRevision"), QString::number(change.previousRevision)},
        {QStringLiteral("revision"), QString::number(change.revision)},
        {QStringLiteral("appliedSeq"), QString::number(change.appliedSeq)},
        {QStringLiteral("acceptedTimeS"), canonicalZero(change.acceptedTimeS)},
        {QStringLiteral("parameterSha256"), change.parameterSha256}};
}

QJsonObject stateObject(const State &state)
{
    QJsonArray history;
    for (const auto &change : state.history) history.append(changeObject(change));
    return QJsonObject{
        {QStringLiteral("schema"), QString::fromLatin1(stateSchema)},
        {QStringLiteral("parameters"), parametersObject(state.parameters)},
        {QStringLiteral("history"), history}};
}

void validateState(const State &state)
{
    validateParameters(state.parameters);
    if (state.history.size() > maximumHistoryCount)
        invalid("Structural runtime accepted history exceeds 128 changes");
    if (state.history.size() != state.parameters.revision)
        invalid("Structural runtime accepted history does not match current revision");

    quint64 expectedRevision = 0;
    quint64 previousAppliedSeq = 0;
    double previousAcceptedTimeS = 0.0;
    QStringList proposalIds;
    for (std::size_t i = 0; i < state.history.size(); ++i)
    {
        const Change &change = state.history[i];
        if (!isCanonicalUuid(change.proposalId) || proposalIds.contains(change.proposalId))
            invalid("Structural runtime history proposal IDs must be unique canonical UUIDs");
        if (expectedRevision == std::numeric_limits<quint64>::max() ||
            change.previousRevision != expectedRevision ||
            change.revision != expectedRevision + 1)
            invalid("Structural runtime accepted history revisions are not contiguous");
        if (i != 0 && (change.appliedSeq <= previousAppliedSeq ||
                       change.acceptedTimeS <= previousAcceptedTimeS))
            invalid("Structural runtime accepted sequence and time must increase");
        if (change.appliedSeq == 0 || !std::isfinite(change.acceptedTimeS) || change.acceptedTimeS < 0.0 ||
            !isLowerSha256(change.parameterSha256))
            invalid("Invalid structural runtime accepted history record");
        proposalIds.append(change.proposalId);
        expectedRevision = change.revision;
        previousAppliedSeq = change.appliedSeq;
        previousAcceptedTimeS = change.acceptedTimeS;
    }
    if (!state.history.empty() &&
        state.history.back().parameterSha256 != state.parameters.parameterSha256)
        invalid("Structural runtime current parameter digest does not match accepted history");
}

void requireExactKeys(const QJsonObject &object,
                      std::initializer_list<QString> expectedKeys)
{
    if (object.size() != qsizetype(expectedKeys.size()))
        invalid("Unexpected or missing structural runtime JSON fields");
    for (const QString &key : expectedKeys)
        if (!object.contains(key)) invalid("Unexpected or missing structural runtime JSON fields");
}

QString readString(const QJsonObject &object, const QString &key)
{
    const QJsonValue value = object.value(key);
    if (!value.isString()) invalid("Structural runtime JSON field must be a string");
    return value.toString();
}

double readFiniteNumber(const QJsonObject &object, const QString &key)
{
    const QJsonValue value = object.value(key);
    if (!value.isDouble()) invalid("Structural runtime JSON field must be a number");
    const double number = value.toDouble();
    if (!std::isfinite(number)) invalid("Structural runtime JSON number must be finite");
    return number;
}

quint64 readDecimalString(const QJsonObject &object, const QString &key)
{
    const QString value = readString(object, key);
    bool ok = false;
    const quint64 parsed = value.toULongLong(&ok, 10);
    if (!ok || QString::number(parsed) != value)
        invalid("Structural runtime revision and sequence fields must be canonical decimal strings");
    return parsed;
}

ActuatorLimit parseActuatorLimit(const QJsonValue &value)
{
    if (!value.isObject()) invalid("Structural runtime actuator limit must be an object");
    const QJsonObject object = value.toObject();
    requireExactKeys(object, {QStringLiteral("actuatorId"), QStringLiteral("minN"),
                              QStringLiteral("maxN"), QStringLiteral("rateNPerS")});
    ActuatorLimit result;
    result.actuatorId = readString(object, QStringLiteral("actuatorId"));
    result.minN = readFiniteNumber(object, QStringLiteral("minN"));
    result.maxN = readFiniteNumber(object, QStringLiteral("maxN"));
    result.rateNPerS = readFiniteNumber(object, QStringLiteral("rateNPerS"));
    validateLimit(result);
    return result;
}

Parameters parseParameters(const QJsonValue &value)
{
    if (!value.isObject()) invalid("Structural runtime parameters must be an object");
    const QJsonObject object = value.toObject();
    requireExactKeys(object, {QStringLiteral("revision"), QStringLiteral("velocityGainNsm"),
                              QStringLiteral("actuatorLimits"),
                              QStringLiteral("fixedBaseControllerRevision"),
                              QStringLiteral("parameterSha256"),
                              QStringLiteral("runtimeControllerRevision")});

    Parameters result;
    result.revision = readDecimalString(object, QStringLiteral("revision"));
    result.velocityGainNsm = readFiniteNumber(object, QStringLiteral("velocityGainNsm"));
    result.fixedBaseControllerRevision =
        readString(object, QStringLiteral("fixedBaseControllerRevision"));
    result.parameterSha256 = readString(object, QStringLiteral("parameterSha256"));
    result.runtimeControllerRevision =
        readString(object, QStringLiteral("runtimeControllerRevision"));

    const QJsonValue limitsValue = object.value(QStringLiteral("actuatorLimits"));
    if (!limitsValue.isArray()) invalid("Structural runtime actuatorLimits must be an array");
    const QJsonArray limits = limitsValue.toArray();
    if (limits.size() > qsizetype(maximumActuatorCount))
        invalid("Structural runtime actuator count exceeds 256");
    result.actuatorLimits.reserve(std::size_t(limits.size()));
    for (const QJsonValue &limit : limits) result.actuatorLimits.push_back(parseActuatorLimit(limit));
    validateParameters(result);
    return result;
}

Change parseChange(const QJsonValue &value)
{
    if (!value.isObject()) invalid("Structural runtime history entry must be an object");
    const QJsonObject object = value.toObject();
    requireExactKeys(object, {QStringLiteral("proposalId"), QStringLiteral("previousRevision"),
                              QStringLiteral("revision"), QStringLiteral("appliedSeq"),
                              QStringLiteral("acceptedTimeS"),
                              QStringLiteral("parameterSha256")});
    Change result;
    result.proposalId = readString(object, QStringLiteral("proposalId"));
    result.previousRevision = readDecimalString(object, QStringLiteral("previousRevision"));
    result.revision = readDecimalString(object, QStringLiteral("revision"));
    result.appliedSeq = readDecimalString(object, QStringLiteral("appliedSeq"));
    result.acceptedTimeS = readFiniteNumber(object, QStringLiteral("acceptedTimeS"));
    result.parameterSha256 = readString(object, QStringLiteral("parameterSha256"));
    return result;
}

}

namespace RStructuralControlRuntimeParameters
{
State fresh(const QString &fixedBaseControllerRevision,
            double velocityGainNsm,
            std::vector<ActuatorLimit> actuatorLimits)
{
    Parameters parameters;
    parameters.revision = 0;
    parameters.velocityGainNsm = velocityGainNsm;
    parameters.actuatorLimits = sortedAndValidatedLimits(std::move(actuatorLimits));
    parameters.fixedBaseControllerRevision = fixedBaseControllerRevision;
    parameters = withComputedDigests(std::move(parameters));
    return State{std::move(parameters), {}};
}

Parameters candidate(const State &accepted,
                     double velocityGainNsm,
                     std::vector<ActuatorLimit> actuatorLimits)
{
    validateState(accepted);
    if (accepted.parameters.revision == std::numeric_limits<quint64>::max())
        invalid("Structural runtime revision overflow");
    actuatorLimits = sortedAndValidatedLimits(std::move(actuatorLimits));
    if (!sameActuatorIds(accepted.parameters.actuatorLimits, actuatorLimits))
        invalid("Structural runtime candidate actuator IDs must match the accepted state");

    Parameters result;
    result.revision = accepted.parameters.revision + 1;
    result.velocityGainNsm = velocityGainNsm;
    result.actuatorLimits = std::move(actuatorLimits);
    result.fixedBaseControllerRevision = accepted.parameters.fixedBaseControllerRevision;
    return withComputedDigests(std::move(result));
}

std::vector<AppliedForceInterval> validateCandidateAgainstBase(
    const Parameters &candidateParameters,
    const std::vector<ActuatorLimit> &immutableLimits,
    const QMap<QString, double> &lastAppliedById,
    double dt)
{
    validateParameters(candidateParameters);
    if (!std::isfinite(dt) || dt <= 0.0)
        invalid("Structural runtime limiter dt must be finite and positive");

    const std::vector<ActuatorLimit> base = sortedAndValidatedLimits(immutableLimits);
    if (!sameActuatorIds(base, candidateParameters.actuatorLimits) ||
        lastAppliedById.size() != qsizetype(base.size()))
        invalid("Structural runtime candidate, base limits and applied forces must use the same actuator IDs");

    std::vector<AppliedForceInterval> result;
    result.reserve(base.size());
    for (std::size_t i = 0; i < base.size(); ++i)
    {
        const ActuatorLimit &baseLimit = base[i];
        const ActuatorLimit &nextLimit = candidateParameters.actuatorLimits[i];
        if (nextLimit.minN < baseLimit.minN || nextLimit.maxN > baseLimit.maxN)
            invalid("Structural runtime candidate may only narrow base amplitude limits");
        if (baseLimit.rateNPerS > 0.0 &&
            (nextLimit.rateNPerS == 0.0 || nextLimit.rateNPerS > baseLimit.rateNPerS))
            invalid("Structural runtime candidate may only narrow a positive base rate limit");

        const auto appliedIt = lastAppliedById.constFind(baseLimit.actuatorId);
        if (appliedIt == lastAppliedById.cend() || !std::isfinite(appliedIt.value()))
            invalid("Structural runtime last-applied force is missing or non-finite");

        double minimum = nextLimit.minN;
        double maximum = nextLimit.maxN;
        if (nextLimit.rateNPerS > 0.0)
        {
            const double delta = nextLimit.rateNPerS * dt;
            const double rateMinimum = appliedIt.value() - delta;
            const double rateMaximum = appliedIt.value() + delta;
            if (!std::isfinite(delta) || !std::isfinite(rateMinimum) ||
                !std::isfinite(rateMaximum))
                invalid("Structural runtime amplitude/rate intersection overflow");
            minimum = std::max(minimum, rateMinimum);
            maximum = std::min(maximum, rateMaximum);
        }
        if (!std::isfinite(minimum) || !std::isfinite(maximum) || minimum > maximum)
            invalid("Structural runtime amplitude and rate limits have an empty intersection");
        result.push_back({baseLimit.actuatorId, minimum, maximum});
    }
    return result;
}

State acceptedChange(const State &accepted,
                     const QString &proposalId,
                     const Parameters &candidateParameters,
                     quint64 appliedSeq,
                     double acceptedTimeS)
{
    validateState(accepted);
    if (accepted.history.size() >= maximumHistoryCount)
        invalid("Structural runtime accepted history is full");
    if (!isCanonicalUuid(proposalId))
        invalid("Structural runtime proposal ID must be a canonical lowercase UUID");
    for (const auto &change : accepted.history)
        if (change.proposalId == proposalId)
            invalid("Structural runtime proposal ID was already accepted");
    if (accepted.parameters.revision == std::numeric_limits<quint64>::max() ||
        candidateParameters.revision != accepted.parameters.revision + 1)
        invalid("Structural runtime candidate revision is not current+1");
    validateParameters(candidateParameters);
    if (candidateParameters.fixedBaseControllerRevision !=
        accepted.parameters.fixedBaseControllerRevision)
        invalid("Structural runtime candidate changed the immutable base controller revision");
    if (!sameActuatorIds(accepted.parameters.actuatorLimits,
                         candidateParameters.actuatorLimits))
        invalid("Structural runtime candidate actuator IDs changed");
    if (appliedSeq == 0 || !std::isfinite(acceptedTimeS) || acceptedTimeS < 0.0)
        invalid("Structural runtime accepted time must be finite and non-negative");
    if (!accepted.history.empty() &&
        (appliedSeq <= accepted.history.back().appliedSeq ||
         acceptedTimeS <= accepted.history.back().acceptedTimeS))
        invalid("Structural runtime accepted sequence and time must increase");

    State result = accepted;
    Change change;
    change.proposalId = proposalId;
    change.previousRevision = accepted.parameters.revision;
    change.revision = candidateParameters.revision;
    change.appliedSeq = appliedSeq;
    change.acceptedTimeS = canonicalZero(acceptedTimeS);
    change.parameterSha256 = candidateParameters.parameterSha256;
    result.parameters = candidateParameters;
    result.history.push_back(std::move(change));
    // Return only a state that can actually be persisted in the bounded
    // checkpoint component; the caller can reject it before publication.
    serialize(result);
    return result;
}

QByteArray serialize(const State &state)
{
    validateState(state);
    const QByteArray bytes = QJsonDocument(stateObject(state)).toJson(QJsonDocument::Compact);
    if (bytes.size() > maximumStateBytes)
        invalid("Structural runtime parameter state exceeds 64 KiB");
    return bytes;
}

State deserialize(const QByteArray &bytes)
{
    if (bytes.isEmpty() || bytes.size() > maximumStateBytes)
        invalid("Structural runtime parameter state must be 1..65536 bytes");
    QJsonParseError parseError;
    const QJsonDocument document = QJsonDocument::fromJson(bytes, &parseError);
    if (parseError.error != QJsonParseError::NoError || !document.isObject())
        invalid("Malformed structural runtime parameter JSON");

    const QJsonObject object = document.object();
    requireExactKeys(object, {QStringLiteral("schema"), QStringLiteral("parameters"),
                              QStringLiteral("history")});
    if (readString(object, QStringLiteral("schema")) != QString::fromLatin1(stateSchema))
        invalid("Unsupported structural runtime parameter schema");

    State result;
    result.parameters = parseParameters(object.value(QStringLiteral("parameters")));
    const QJsonValue historyValue = object.value(QStringLiteral("history"));
    if (!historyValue.isArray()) invalid("Structural runtime history must be an array");
    const QJsonArray history = historyValue.toArray();
    if (history.size() > qsizetype(maximumHistoryCount))
        invalid("Structural runtime accepted history exceeds 128 changes");
    result.history.reserve(std::size_t(history.size()));
    for (const QJsonValue &change : history) result.history.push_back(parseChange(change));
    validateState(result);

    // This exact-byte check rejects duplicate fields (which Qt's JSON object
    // parser would otherwise collapse), non-compact variants, and any trailing
    // bytes. Only codec-produced state is accepted as persisted state.
    if (serialize(result) != bytes)
        invalid("Structural runtime parameter JSON is not canonical compact codec output");
    return result;
}
}

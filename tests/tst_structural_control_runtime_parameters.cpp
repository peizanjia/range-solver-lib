#include <QtTest>

#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include "rstructural_control_runtime_parameters.h"

using namespace RStructuralControlRuntimeParameters;

namespace
{
QString baseRevision()
{
    return QString(64, QLatin1Char('a'));
}

std::vector<ActuatorLimit> baseLimits()
{
    return {{QStringLiteral("force-a"), -10.0, 10.0, 2.0},
            {QStringLiteral("force-b"), 0.0, 100.0, 0.0}};
}

std::vector<ActuatorLimit> narrowedLimits()
{
    return {{QStringLiteral("force-a"), -5.0, 7.0, 1.0},
            {QStringLiteral("force-b"), 10.0, 90.0, 0.0}};
}

QString proposalId(quint64 index)
{
    return QStringLiteral("00000000-0000-4000-8000-%1")
        .arg(index, 12, 16, QLatin1Char('0'));
}

QByteArray compact(const QJsonObject &object)
{
    return QJsonDocument(object).toJson(QJsonDocument::Compact);
}
}

class TestStructuralControlRuntimeParameters : public QObject
{
    Q_OBJECT

private slots:
    void freshStateIsDistinctAndRoundTripsCanonicalBytes();
    void hashesAreDeterministicAndBindRevisionAndRuntimeBase();
    void portableHashesMatchGoldenVector();
    void signedZeroIsCanonicalForHashesAndStateBytes();
    void candidateValidationNarrowsAndIntersectsAmplitudeAndRate();
    void acceptedChangeKeepsOnlyAcceptedContiguousHistory();
    void rejectsInvalidLimitsAndNumericRanges();
    void rejectsMalformedNonCanonicalAndForgedJson();
    void rejectsHistoryOverflowWithoutMutation();
    void rejectsCombinedStateByteLimitBeforeReturningAcceptedState();
};

void TestStructuralControlRuntimeParameters::freshStateIsDistinctAndRoundTripsCanonicalBytes()
{
    const State state = fresh(baseRevision(), 125.0,
                              {{QStringLiteral("force-b"), 0.0, 100.0, 0.0},
                               {QStringLiteral("force-a"), -10.0, 10.0, 2.0}});
    const QByteArray bytes = serialize(state);
    QVERIFY(!bytes.isEmpty()); // A revision-zero parameter state is not an empty RCCE component.
    QVERIFY(bytes.contains(QByteArrayLiteral("range.structural-runtime-parameters/v1")));
    QVERIFY(bytes.contains(QByteArrayLiteral("\"history\":[]")));
    QCOMPARE(state.parameters.revision, quint64(0));
    QCOMPARE(state.parameters.actuatorLimits.front().actuatorId, QStringLiteral("force-a"));
    QCOMPARE(state.parameters.parameterSha256.size(), 64);
    QCOMPARE(state.parameters.runtimeControllerRevision.size(), 64);

    const State restored = deserialize(bytes);
    QCOMPARE(serialize(restored), bytes);
    QCOMPARE(restored.parameters.parameterSha256, state.parameters.parameterSha256);
    QCOMPARE(restored.parameters.runtimeControllerRevision,
             state.parameters.runtimeControllerRevision);
    QVERIFY(restored.history.empty());
}

void TestStructuralControlRuntimeParameters::hashesAreDeterministicAndBindRevisionAndRuntimeBase()
{
    const State first = fresh(baseRevision(), 125.0, baseLimits());
    const State sortedEquivalent = fresh(baseRevision(), 125.0,
        {{QStringLiteral("force-b"), 0.0, 100.0, 0.0},
         {QStringLiteral("force-a"), -10.0, 10.0, 2.0}});
    QCOMPARE(first.parameters.parameterSha256, sortedEquivalent.parameters.parameterSha256);
    QCOMPARE(first.parameters.runtimeControllerRevision,
             sortedEquivalent.parameters.runtimeControllerRevision);
    QCOMPARE(serialize(first), serialize(sortedEquivalent));

    const Parameters changed = candidate(first, 126.0, baseLimits());
    QCOMPARE(changed.revision, quint64(1));
    QVERIFY(changed.parameterSha256 != first.parameters.parameterSha256);
    QVERIFY(changed.runtimeControllerRevision != first.parameters.runtimeControllerRevision);
    QVERIFY(changed.fixedBaseControllerRevision == first.parameters.fixedBaseControllerRevision);
}

void TestStructuralControlRuntimeParameters::portableHashesMatchGoldenVector()
{
    const State state = fresh(QString(64, QLatin1Char('0')), 1.0,
        {{QStringLiteral("f"), -1.0, 1.0, 0.0}});
    QCOMPARE(state.parameters.parameterSha256,
             QStringLiteral("7f99575c73b6790a5748e6275a62127385d496df16aa7bb1336a425aa1ee54cc"));
    QCOMPARE(state.parameters.runtimeControllerRevision,
             QStringLiteral("a6546d489a5f102a1de13955556b011d89dbbef9e97be40a61048ea54e817669"));
}

void TestStructuralControlRuntimeParameters::signedZeroIsCanonicalForHashesAndStateBytes()
{
    const State negativeZero = fresh(QString(64, QLatin1Char('0')), 125.0,
        {{QStringLiteral("force"), -0.0, 1.0, -0.0}});
    const State positiveZero = fresh(QString(64, QLatin1Char('0')), 125.0,
        {{QStringLiteral("force"), 0.0, 1.0, 0.0}});
    QVERIFY(!std::signbit(negativeZero.parameters.actuatorLimits[0].minN));
    QVERIFY(!std::signbit(negativeZero.parameters.actuatorLimits[0].rateNPerS));
    QCOMPARE(negativeZero.parameters.parameterSha256, positiveZero.parameters.parameterSha256);
    QCOMPARE(serialize(negativeZero), serialize(positiveZero));

    const State initial = fresh(baseRevision(), 125.0, baseLimits());
    const State negativeTime = acceptedChange(initial, proposalId(1),
        candidate(initial, 130.0, baseLimits()), 1, -0.0);
    const State positiveTime = acceptedChange(initial, proposalId(1),
        candidate(initial, 130.0, baseLimits()), 1, 0.0);
    QVERIFY(!std::signbit(negativeTime.history[0].acceptedTimeS));
    QCOMPARE(serialize(negativeTime), serialize(positiveTime));
}

void TestStructuralControlRuntimeParameters::candidateValidationNarrowsAndIntersectsAmplitudeAndRate()
{
    const State state = fresh(baseRevision(), 125.0, baseLimits());
    const Parameters proposed = candidate(state, 130.0, narrowedLimits());
    const auto intervals = validateCandidateAgainstBase(
        proposed, baseLimits(),
        {{QStringLiteral("force-a"), 4.0}, {QStringLiteral("force-b"), 20.0}}, 2.0);
    QCOMPARE(intervals.size(), std::size_t(2));
    QCOMPARE(intervals[0].actuatorId, QStringLiteral("force-a"));
    QCOMPARE(intervals[0].minN, 2.0);
    QCOMPARE(intervals[0].maxN, 6.0);
    QCOMPARE(intervals[1].actuatorId, QStringLiteral("force-b"));
    QCOMPARE(intervals[1].minN, 10.0);
    QCOMPARE(intervals[1].maxN, 90.0);

    auto widenedAmplitude = narrowedLimits();
    widenedAmplitude[0].minN = -11.0;
    QVERIFY_EXCEPTION_THROWN(
        validateCandidateAgainstBase(candidate(state, 130.0, widenedAmplitude), baseLimits(),
            {{QStringLiteral("force-a"), 0.0}, {QStringLiteral("force-b"), 0.0}}, 0.5),
        std::invalid_argument);

    auto disabledPositiveBaseRate = narrowedLimits();
    disabledPositiveBaseRate[0].rateNPerS = 0.0;
    QVERIFY_EXCEPTION_THROWN(
        validateCandidateAgainstBase(candidate(state, 130.0, disabledPositiveBaseRate),
            baseLimits(), {{QStringLiteral("force-a"), 0.0},
                           {QStringLiteral("force-b"), 0.0}}, 0.5),
        std::invalid_argument);

    auto widenedRate = narrowedLimits();
    widenedRate[0].rateNPerS = 3.0;
    QVERIFY_EXCEPTION_THROWN(
        validateCandidateAgainstBase(candidate(state, 130.0, widenedRate), baseLimits(),
            {{QStringLiteral("force-a"), 0.0}, {QStringLiteral("force-b"), 0.0}}, 0.5),
        std::invalid_argument);

    QVERIFY_EXCEPTION_THROWN(
        validateCandidateAgainstBase(proposed, baseLimits(),
            {{QStringLiteral("force-a"), 50.0}, {QStringLiteral("force-b"), 20.0}}, 1.0),
        std::invalid_argument); // The new amplitude and rate interval do not overlap.

    const std::vector<ActuatorLimit> hugeRateBase = {
        {QStringLiteral("force"), -1.0e308, 1.0e308, 1.0e308}};
    const State hugeRateState = fresh(baseRevision(), 125.0, hugeRateBase);
    const Parameters hugeRateCandidate = candidate(hugeRateState, 130.0, hugeRateBase);
    QVERIFY_EXCEPTION_THROWN(validateCandidateAgainstBase(
        hugeRateCandidate, hugeRateBase, {{QStringLiteral("force"), 0.0}}, 2.0),
        std::invalid_argument); // rate * dt overflows.
    QVERIFY_EXCEPTION_THROWN(validateCandidateAgainstBase(
        hugeRateCandidate, hugeRateBase, {{QStringLiteral("force"), 1.0e308}}, 1.0),
        std::invalid_argument); // old force + rate * dt overflows.
}

void TestStructuralControlRuntimeParameters::acceptedChangeKeepsOnlyAcceptedContiguousHistory()
{
    const State initial = fresh(baseRevision(), 125.0, baseLimits());
    const QByteArray initialBytes = serialize(initial);
    const Parameters next = candidate(initial, 130.0, narrowedLimits());
    const State accepted = acceptedChange(initial, proposalId(1), next, 7, 0.25);

    QCOMPARE(serialize(initial), initialBytes);
    QCOMPARE(accepted.parameters.revision, quint64(1));
    QCOMPARE(accepted.history.size(), std::size_t(1));
    QCOMPARE(accepted.history[0].previousRevision, quint64(0));
    QCOMPARE(accepted.history[0].revision, quint64(1));
    QCOMPARE(accepted.history[0].appliedSeq, quint64(7));
    QCOMPARE(accepted.history[0].acceptedTimeS, 0.25);
    QCOMPARE(accepted.history[0].parameterSha256, accepted.parameters.parameterSha256);

    const State second = acceptedChange(accepted, proposalId(2),
        candidate(accepted, 140.0, narrowedLimits()), 9, 0.5);
    QCOMPARE(second.parameters.revision, quint64(2));
    QCOMPARE(second.history.size(), std::size_t(2));
    QCOMPARE(second.history[1].previousRevision, quint64(1));
    QCOMPARE(deserialize(serialize(second)).history.size(), std::size_t(2));

    QVERIFY_EXCEPTION_THROWN(acceptedChange(accepted, proposalId(1),
        candidate(accepted, 140.0, narrowedLimits()), 9, 0.5), std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(acceptedChange(accepted, proposalId(3),
        candidate(accepted, 140.0, narrowedLimits()), 7, 0.5), std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(acceptedChange(accepted, proposalId(3),
        candidate(accepted, 140.0, narrowedLimits()), 9, 0.25), std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(acceptedChange(initial, proposalId(4), initial.parameters, 1, 0.1),
                             std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(acceptedChange(initial, proposalId(4), next, 0, 0.1),
                             std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(validateCandidateAgainstBase(next, baseLimits(),
        QMap<QString, double>{{QStringLiteral("force-a"), 0.0}, {QStringLiteral("force-b"), 0.0}}, 0.0),
                             std::invalid_argument);
}

void TestStructuralControlRuntimeParameters::rejectsInvalidLimitsAndNumericRanges()
{
    QVERIFY_EXCEPTION_THROWN(fresh(baseRevision(), 125.0, {}), std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(fresh(QString(64, QLatin1Char('A')), 125.0, baseLimits()),
                             std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(fresh(baseRevision(), 0.0, baseLimits()), std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(fresh(baseRevision(), 1.0e-10, baseLimits()), std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(fresh(baseRevision(), 1.0e13, baseLimits()), std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(fresh(baseRevision(), std::numeric_limits<double>::infinity(),
                                   baseLimits()), std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(fresh(baseRevision(), std::numeric_limits<double>::quiet_NaN(),
                                   baseLimits()), std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(fresh(baseRevision(), 125.0,
        {{QStringLiteral("force-a"), 2.0, 1.0, 1.0}}), std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(fresh(baseRevision(), 125.0,
        {{QStringLiteral("force-a"), 0.0, 1.0, -1.0}}), std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(fresh(baseRevision(), 125.0,
        {{QStringLiteral("force-a"), 0.0, 1.0,
          std::numeric_limits<double>::infinity()}}), std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(fresh(baseRevision(), 125.0,
        {{QStringLiteral("bad id"), 0.0, 1.0, 1.0}}), std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(fresh(baseRevision(), 125.0,
        {{QStringLiteral("force-a"), 0.0, 1.0, 1.0},
         {QStringLiteral("force-a"), -1.0, 2.0, 1.0}}), std::invalid_argument);

    std::vector<ActuatorLimit> tooMany;
    tooMany.reserve(maximumActuatorCount + 1);
    for (std::size_t i = 0; i <= maximumActuatorCount; ++i)
        tooMany.push_back({QStringLiteral("actuator-%1").arg(qulonglong(i)), 0.0, 1.0, 0.0});
    QVERIFY_EXCEPTION_THROWN(fresh(baseRevision(), 125.0, std::move(tooMany)),
                             std::invalid_argument);

    const State state = fresh(baseRevision(), 125.0, baseLimits());
    const Parameters valid = candidate(state, 130.0, narrowedLimits());
    auto forged = valid;
    forged.fixedBaseControllerRevision = QString(64, QLatin1Char('b'));
    QVERIFY_EXCEPTION_THROWN(acceptedChange(state, proposalId(1), forged, 1, 0.1),
                             std::invalid_argument);
    forged = valid;
    forged.parameterSha256 = QString(64, QLatin1Char('f'));
    QVERIFY_EXCEPTION_THROWN(acceptedChange(state, proposalId(1), forged, 1, 0.1),
                             std::invalid_argument);
    forged = valid;
    forged.runtimeControllerRevision = QString(64, QLatin1Char('f'));
    QVERIFY_EXCEPTION_THROWN(acceptedChange(state, proposalId(1), forged, 1, 0.1),
                             std::invalid_argument);
}

void TestStructuralControlRuntimeParameters::rejectsMalformedNonCanonicalAndForgedJson()
{
    const QByteArray valid = serialize(fresh(baseRevision(), 125.0, baseLimits()));
    QVERIFY_EXCEPTION_THROWN(deserialize(valid + QByteArrayLiteral("x")), std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(deserialize(valid + QByteArrayLiteral(" ")), std::invalid_argument);

    QByteArray duplicateSchema = valid;
    duplicateSchema.replace(QByteArrayLiteral("\"schema\":"),
                            QByteArrayLiteral("\"schema\":\"range.structural-runtime-parameters/v1\",\"schema\":"));
    QVERIFY_EXCEPTION_THROWN(deserialize(duplicateSchema), std::invalid_argument);

    QJsonObject object = QJsonDocument::fromJson(valid).object();
    object.insert(QStringLiteral("unexpected"), true);
    QVERIFY_EXCEPTION_THROWN(deserialize(compact(object)), std::invalid_argument);

    object = QJsonDocument::fromJson(valid).object();
    object.insert(QStringLiteral("schema"), QStringLiteral("range.structural-runtime-parameters/v2"));
    QVERIFY_EXCEPTION_THROWN(deserialize(compact(object)), std::invalid_argument);

    object = QJsonDocument::fromJson(valid).object();
    QJsonObject parameters = object.value(QStringLiteral("parameters")).toObject();
    parameters.insert(QStringLiteral("revision"), 0);
    object.insert(QStringLiteral("parameters"), parameters);
    QVERIFY_EXCEPTION_THROWN(deserialize(compact(object)), std::invalid_argument);

    object = QJsonDocument::fromJson(valid).object();
    parameters = object.value(QStringLiteral("parameters")).toObject();
    parameters.insert(QStringLiteral("velocityGainNsm"), true);
    object.insert(QStringLiteral("parameters"), parameters);
    QVERIFY_EXCEPTION_THROWN(deserialize(compact(object)), std::invalid_argument);

    object = QJsonDocument::fromJson(valid).object();
    parameters = object.value(QStringLiteral("parameters")).toObject();
    parameters.insert(QStringLiteral("revision"), QStringLiteral("00"));
    object.insert(QStringLiteral("parameters"), parameters);
    QVERIFY_EXCEPTION_THROWN(deserialize(compact(object)), std::invalid_argument);

    object = QJsonDocument::fromJson(valid).object();
    parameters = object.value(QStringLiteral("parameters")).toObject();
    parameters.insert(QStringLiteral("parameterSha256"), QString(64, QLatin1Char('F')));
    object.insert(QStringLiteral("parameters"), parameters);
    QVERIFY_EXCEPTION_THROWN(deserialize(compact(object)), std::invalid_argument);

    object = QJsonDocument::fromJson(valid).object();
    parameters = object.value(QStringLiteral("parameters")).toObject();
    parameters.insert(QStringLiteral("runtimeControllerRevision"),
                      QString(64, QLatin1Char('f')));
    object.insert(QStringLiteral("parameters"), parameters);
    QVERIFY_EXCEPTION_THROWN(deserialize(compact(object)), std::invalid_argument);

    object = QJsonDocument::fromJson(valid).object();
    parameters = object.value(QStringLiteral("parameters")).toObject();
    QJsonArray limits = parameters.value(QStringLiteral("actuatorLimits")).toArray();
    QJsonObject firstLimit = limits[0].toObject();
    firstLimit.insert(QStringLiteral("extra"), 1.0);
    limits[0] = firstLimit;
    parameters.insert(QStringLiteral("actuatorLimits"), limits);
    object.insert(QStringLiteral("parameters"), parameters);
    QVERIFY_EXCEPTION_THROWN(deserialize(compact(object)), std::invalid_argument);

    const State accepted = acceptedChange(fresh(baseRevision(), 125.0, baseLimits()),
        proposalId(1), candidate(fresh(baseRevision(), 125.0, baseLimits()), 130.0, narrowedLimits()),
        1, 0.1);
    object = QJsonDocument::fromJson(serialize(accepted)).object();
    QJsonArray history = object.value(QStringLiteral("history")).toArray();
    QJsonObject change = history[0].toObject();
    change.insert(QStringLiteral("appliedSeq"), QStringLiteral("01"));
    history[0] = change;
    object.insert(QStringLiteral("history"), history);
    QVERIFY_EXCEPTION_THROWN(deserialize(compact(object)), std::invalid_argument);

    object = QJsonDocument::fromJson(serialize(accepted)).object();
    history = object.value(QStringLiteral("history")).toArray();
    change = history[0].toObject();
    change.insert(QStringLiteral("revision"), QStringLiteral("2"));
    history[0] = change;
    object.insert(QStringLiteral("history"), history);
    QVERIFY_EXCEPTION_THROWN(deserialize(compact(object)), std::invalid_argument);
}

void TestStructuralControlRuntimeParameters::rejectsHistoryOverflowWithoutMutation()
{
    State accepted = fresh(baseRevision(), 125.0, baseLimits());
    for (quint64 i = 0; i < maximumHistoryCount; ++i)
    {
        const Parameters next = candidate(accepted, 125.0 + double(i + 1), baseLimits());
        accepted = acceptedChange(accepted, proposalId(i + 1), next, i + 1, double(i + 1));
    }
    QCOMPARE(accepted.history.size(), maximumHistoryCount);
    QCOMPARE(accepted.parameters.revision, quint64(maximumHistoryCount));
    const QByteArray before = serialize(accepted);

    const Parameters beyondLimit = candidate(accepted, 999.0, baseLimits());
    QVERIFY_EXCEPTION_THROWN(acceptedChange(accepted, proposalId(999), beyondLimit,
                                            quint64(maximumHistoryCount + 1), 999.0),
                             std::invalid_argument);
    QCOMPARE(serialize(accepted), before);
}

void TestStructuralControlRuntimeParameters::rejectsCombinedStateByteLimitBeforeReturningAcceptedState()
{
    std::vector<ActuatorLimit> limits;
    for (std::size_t i = 0; i < maximumActuatorCount; ++i)
        limits.push_back({QStringLiteral("actuator-%1-").arg(qulonglong(i), 3, 10, QLatin1Char('0'))
                              .leftJustified(128, QLatin1Char('x')), -100.0, 100.0, 0.0});
    State accepted = fresh(baseRevision(), 125.0, limits);
    QVERIFY(serialize(accepted).size() <= maximumStateBytes);
    bool rejected = false;
    for (quint64 i = 0; i < maximumHistoryCount; ++i)
    {
        const QByteArray before = serialize(accepted);
        try
        {
            accepted = acceptedChange(accepted, proposalId(i + 1),
                candidate(accepted, 126.0 + double(i), limits), i + 1, double(i + 1));
        }
        catch (const std::invalid_argument &)
        {
            rejected = true;
            QCOMPARE(serialize(accepted), before);
            QVERIFY(accepted.history.size() < maximumHistoryCount);
            break;
        }
    }
    QVERIFY(rejected);
}

QTEST_APPLESS_MAIN(TestStructuralControlRuntimeParameters)

#include "tst_structural_control_runtime_parameters.moc"

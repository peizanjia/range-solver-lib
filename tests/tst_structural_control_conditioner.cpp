#include <QtTest>

#include <cmath>
#include <limits>
#include <stdexcept>

#include "rstructural_control_conditioner.h"
#include "rstructural_control_runtime_parameters.h"

class TestStructuralControlConditioner : public QObject
{
    Q_OBJECT

private slots:
    void defaultPolicyPreservesRequestedPhysicalCommandExactly();
    void delayAndZeroOrderHold();
    void deterministicNoiseAndRejectedTrialReplay();
    void multipleChannelsAndStateRoundTrip();
    void rejectsChangedConfigurationAndMalformedState();
    void compositeEnvelopeRoundTripAndValidation();
    void compositeV2RuntimeComponentPreservesLegacyBytesAndRejectsCorruption();
};

void TestStructuralControlConditioner::defaultPolicyPreservesRequestedPhysicalCommandExactly()
{
    RStructuralControlConditioner conditioner(31,{{QStringLiteral("force"),0,1,0.0}});
    auto state = conditioner.initialState();
    const double analyticCommands[] = {0.0,0.125,-0.5,1024.0,-1.0/3.0};
    for (size_t i=0;i<std::size(analyticCommands);++i)
    {
        const auto result = conditioner.advance(state,{{QStringLiteral("force"),analyticCommands[i]}});
        QCOMPARE(result.appliedCommands.size(),size_t(1));
        QCOMPARE(result.appliedCommands.front().value,analyticCommands[i]);
        QCOMPARE(result.candidate.heldCommands,std::vector<double>{analyticCommands[i]});
        QCOMPARE(result.candidate.completedSteps,quint64(i+1));
        QCOMPARE(result.candidate.delayedCommands.size(),size_t(1));
        QVERIFY(result.candidate.delayedCommands.front().empty());
        const auto roundTrip=conditioner.deserialize(conditioner.serialize(result.candidate));
        QCOMPARE(roundTrip.heldCommands,result.candidate.heldCommands);
        QCOMPARE(roundTrip.completedSteps,result.candidate.completedSteps);
        state=result.candidate;
    }
}

void TestStructuralControlConditioner::delayAndZeroOrderHold()
{
    RStructuralControlConditioner conditioner(7,{{QStringLiteral("a"),1,2,0.0}});
    auto accepted = conditioner.initialState();
    const double input[] = {1.0,2.0,3.0,4.0,5.0};
    const double expected[] = {0.0,1.0,1.0,3.0,3.0};
    for (int step=0;step<5;step++)
    {
        const auto candidate = conditioner.advance(accepted,{{QStringLiteral("a"),input[step]}});
        QCOMPARE(candidate.appliedCommands.size(),size_t(1));
        QCOMPARE(candidate.appliedCommands[0].actuatorId,QStringLiteral("a"));
        QCOMPARE(candidate.appliedCommands[0].value,expected[step]);
        accepted = candidate.candidate;
        QCOMPARE(accepted.completedSteps,quint64(step+1));
    }
}

void TestStructuralControlConditioner::deterministicNoiseAndRejectedTrialReplay()
{
    const std::vector<RStructuralControlConditioningPolicy> policy =
        {{QStringLiteral("force.x"),0,2,0.3}};
    RStructuralControlConditioner first(123456789,policy), second(123456789,policy);
    auto state = first.initialState();
    for (int step=0;step<8;step++)
    {
        const std::vector<RStructuralControlCommand> request =
            {{QStringLiteral("force.x"),double(step)}};
        const auto candidate = first.advance(state,request);
        // A rejected trial does not mutate the accepted state or RNG stream.
        const auto retry = first.advance(state,request);
        QCOMPARE(candidate.appliedCommands[0].value,retry.appliedCommands[0].value);
        QCOMPARE(first.serialize(candidate.candidate),first.serialize(retry.candidate));
        const auto recovered = second.deserialize(first.serialize(state));
        const auto replay = second.advance(recovered,request);
        QCOMPARE(candidate.appliedCommands[0].value,replay.appliedCommands[0].value);
        QCOMPARE(first.serialize(candidate.candidate),second.serialize(replay.candidate));
        QVERIFY(std::abs(candidate.appliedCommands[0].value-double(step)) <= 0.3 ||
                step % 2 == 1);
        state = candidate.candidate;
    }
    RStructuralControlConditioner otherSeed(123456790,policy);
    const auto different = otherSeed.advance(otherSeed.initialState(),
                                             {{QStringLiteral("force.x"),0.0}});
    const auto original = first.advance(first.initialState(),
                                         {{QStringLiteral("force.x"),0.0}});
    QVERIFY(different.appliedCommands[0].value != original.appliedCommands[0].value);
}

void TestStructuralControlConditioner::multipleChannelsAndStateRoundTrip()
{
    RStructuralControlConditioner conditioner(0,{
        {QStringLiteral("left"),2,1,0.0},
        {QStringLiteral("right"),0,3,0.0}});
    auto state = conditioner.initialState();
    const auto first = conditioner.advance(state,{
        {QStringLiteral("left"),10.0},{QStringLiteral("right"),20.0}});
    QCOMPARE(first.appliedCommands[0].value,0.0);
    QCOMPARE(first.appliedCommands[1].value,20.0);
    const QByteArray bytes = conditioner.serialize(first.candidate);
    const auto restored = conditioner.deserialize(bytes);
    QCOMPARE(conditioner.serialize(restored),bytes);
    const auto second = conditioner.advance(restored,{
        {QStringLiteral("left"),11.0},{QStringLiteral("right"),99.0}});
    QCOMPARE(second.appliedCommands[0].value,0.0);
    QCOMPARE(second.appliedCommands[1].value,20.0);
    const auto third = conditioner.advance(second.candidate,{
        {QStringLiteral("left"),12.0},{QStringLiteral("right"),99.0}});
    QCOMPARE(third.appliedCommands[0].value,10.0);
    QCOMPARE(third.appliedCommands[1].value,20.0);
}

void TestStructuralControlConditioner::rejectsChangedConfigurationAndMalformedState()
{
    RStructuralControlConditioner conditioner(42,{{QStringLiteral("a"),1,1,0.0}});
    RStructuralControlConditioner changed(43,{{QStringLiteral("a"),1,1,0.0}});
    const QByteArray bytes = conditioner.serialize(conditioner.initialState());
    QVERIFY_EXCEPTION_THROWN(changed.deserialize(bytes),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(conditioner.deserialize(bytes.left(bytes.size()-1)),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(conditioner.deserialize(bytes + 'x'),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(conditioner.advance(conditioner.initialState(),
                                                 {{QStringLiteral("wrong"),1.0}}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(conditioner.advance(conditioner.initialState(),
                                                 {{QStringLiteral("a"),std::numeric_limits<double>::quiet_NaN()}}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralControlConditioner(0,{{QStringLiteral("a"),0,0,0.0}}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralControlConditioner(0,{
        {QStringLiteral("a"),0,1,0.0},{QStringLiteral("a"),0,1,0.0}}),std::invalid_argument);
}

void TestStructuralControlConditioner::compositeEnvelopeRoundTripAndValidation()
{
    RStructuralControlConditioner conditioner(1,{{QStringLiteral("actuator"),1,1,0.0}});
    RStructuralControlStateEnvelope envelope;
    envelope.controllerState = QByteArray::fromHex("00343000feff");
    envelope.conditionerState = conditioner.serialize(conditioner.initialState());
    const QByteArray packed = envelope.serialize();
    const auto restored = RStructuralControlStateEnvelope::deserialize(packed);
    QCOMPARE(restored.controllerState,envelope.controllerState);
    QCOMPARE(restored.conditionerState,envelope.conditionerState);
    QCOMPARE(conditioner.serialize(conditioner.deserialize(restored.conditionerState)),
             envelope.conditionerState);
    QVERIFY_EXCEPTION_THROWN(RStructuralControlStateEnvelope::deserialize(packed.left(packed.size()-1)),
                             std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralControlStateEnvelope::deserialize(packed + 'x'),
                             std::invalid_argument);
    envelope.controllerState.clear();
    QCOMPARE(RStructuralControlStateEnvelope::deserialize(envelope.serialize()).controllerState,
             QByteArray());
    envelope.controllerState = QByteArray(1048576,'x');
    QVERIFY_EXCEPTION_THROWN(envelope.serialize(),std::invalid_argument);
}

void TestStructuralControlConditioner::compositeV2RuntimeComponentPreservesLegacyBytesAndRejectsCorruption()
{
    RStructuralControlStateEnvelope envelope{QByteArray("a"), QByteArray("bc")};
    const QByteArray v1 = QByteArray::fromHex("52434345010000000100000061020000006263");
    QCOMPARE(envelope.serialize(), v1);
    QCOMPARE(RStructuralControlStateEnvelope::deserialize(v1).serialize(), v1);
    namespace Parameters = RStructuralControlRuntimeParameters;
    envelope.runtimeParametersState = Parameters::serialize(Parameters::fresh(
        QString(64, QLatin1Char('a')), 125.0, {{QStringLiteral("force"), -100.0, 100.0, 0.0}}));
    const QByteArray v2 = envelope.serialize();
    QCOMPARE(v2.mid(4, 4), QByteArray::fromHex("02000000"));
    const auto restored = RStructuralControlStateEnvelope::deserialize(v2);
    QCOMPARE(restored.controllerState, envelope.controllerState);
    QCOMPARE(restored.conditionerState, envelope.conditionerState);
    QCOMPARE(restored.runtimeParametersState, envelope.runtimeParametersState);
    QCOMPARE(restored.serialize(), v2);
    QVERIFY_EXCEPTION_THROWN(RStructuralControlStateEnvelope::deserialize(v2.left(v2.size()-1)), std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralControlStateEnvelope::deserialize(v2 + 'x'), std::invalid_argument);
    QByteArray unsupported = v2;
    unsupported[4] = 3;
    QVERIFY_EXCEPTION_THROWN(RStructuralControlStateEnvelope::deserialize(unsupported), std::invalid_argument);
    QByteArray forged = v2;
    forged[forged.size()-2] = 'x';
    QVERIFY_EXCEPTION_THROWN(RStructuralControlStateEnvelope::deserialize(forged), std::invalid_argument);
    envelope.runtimeParametersState = QByteArray("{}");
    QVERIFY_EXCEPTION_THROWN(envelope.serialize(), std::invalid_argument);
    envelope.runtimeParametersState.clear();
    QCOMPARE(envelope.serialize(), v1);
}

QTEST_APPLESS_MAIN(TestStructuralControlConditioner)

#include "tst_structural_control_conditioner.moc"

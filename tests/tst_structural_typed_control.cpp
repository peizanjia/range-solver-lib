#include <QtTest>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include "rstructural_typed_control.h"

namespace
{
using Quantity = RStructuralTypedControlQuantity;
using Channel = RStructuralTypedControlChannel;
using Command = RStructuralTypedControlCommand;
using Decision = RStructuralTypedControlDecision;
using Session = RStructuralTypedControlSession;
using Sensor = RStructuralControlSensorSample;

Channel channel(QString id,
                Quantity quantity = Quantity::Force,
                double minimum = -100.0,
                double maximum = 100.0,
                double rateLimit = std::numeric_limits<double>::infinity(),
                double initialCommand = 0.0,
                uint delaySteps = 0,
                uint holdSteps = 1,
                double noiseAmplitude = 0.0)
{
    Channel result;
    result.channelId = std::move(id);
    result.quantity = quantity;
    result.minimum = minimum;
    result.maximum = maximum;
    result.rateLimit = rateLimit;
    result.initialCommand = initialCommand;
    result.delaySteps = delaySteps;
    result.holdSteps = holdSteps;
    result.noiseAmplitude = noiseAmplitude;
    return result;
}

Decision decision(const Session &session,
                  const std::vector<double> &values,
                  QByteArray nextControllerState = {})
{
    Decision result;
    result.nextControllerState = std::move(nextControllerState);
    const auto &catalog = session.catalog();
    for (std::size_t i=0;i<catalog.size();i++)
        result.requestedCommands.push_back(
            {catalog[i].channelId,catalog[i].quantity,values[i]});
    return result;
}

Session oneChannelSession(double dt = 0.5,
                          quint64 seed = 1234,
                          Channel channelDefinition = channel(QStringLiteral("actuator")),
                          QString controllerIdentity = QStringLiteral("cpp-typed-v1"))
{
    std::vector<Channel> channels;
    channels.push_back(std::move(channelDefinition));
    return Session(seed,dt,std::move(channels),std::move(controllerIdentity));
}

bool close(double actual, double expected, double tolerance = 1.0e-12)
{
    return std::isfinite(actual) && std::isfinite(expected) &&
           std::abs(actual-expected) <= tolerance*(1.0+std::abs(expected));
}

Session::State takeStep(const Session &session, const Session::State &accepted, double value)
{
    const auto candidate = session.trial(
        accepted,{{QStringLiteral("tip"),0.25,QStringLiteral("m")}},
        [&session,value](const RStructuralControlObservation &observation)
        {
            return decision(session,{value},
                QByteArray::number(observation.priorControllerState.toULongLong()+1));
        });
    return session.withWork(candidate,{0.25}).state;
}
}

class TestStructuralTypedControl : public QObject
{
    Q_OBJECT

private slots:
    void quantityUnitsAndCatalogValidation();
    void rejectWrongIdsQuantitiesDuplicatesAndNonFiniteCommands();
    void applyRateAndAmplitudeLimitsInChannelUnits();
    void applyDeterministicDelayHoldAndNoise();
    void retriesObserveOnlyPriorAcceptedState();
    void serializedSplitTrajectoryIsExactAndCorruptionIsRejected();
    void rejectForgedClockAndSensorMetadata();
    void accountForDriverSuppliedWorkPerChannel();
};

void TestStructuralTypedControl::quantityUnitsAndCatalogValidation()
{
    QCOMPARE(rStructuralTypedControlUnit(Quantity::Force),QStringLiteral("N"));
    QCOMPARE(rStructuralTypedControlUnit(Quantity::CableReferenceLength),QStringLiteral("m"));
    QCOMPARE(rStructuralTypedControlUnit(Quantity::CableInitialTension),QStringLiteral("N"));
    QCOMPARE(rStructuralTypedControlUnit(Quantity::CableAnchorDisplacement),QStringLiteral("m"));
    QCOMPARE(rStructuralTypedControlUnit(Quantity::Voltage),QStringLiteral("V"));
    QCOMPARE(rStructuralTypedControlUnit(Quantity::Charge),QStringLiteral("C"));
    QVERIFY_EXCEPTION_THROWN(
        rStructuralTypedControlUnit(static_cast<Quantity>(99)),std::invalid_argument);

    const auto unrestricted = oneChannelSession();
    QVERIFY(!unrestricted.configurationSha256().isEmpty());
    QVERIFY_EXCEPTION_THROWN(
        oneChannelSession(0.0),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(
        oneChannelSession(0.5,1,channel(QStringLiteral("bad"),Quantity::Voltage,2.0,1.0)),
        std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(
        oneChannelSession(0.5,1,channel(QStringLiteral("bad"),Quantity::Voltage,0.0,1.0,
                                         4.0,2.0)),
        std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(
        oneChannelSession(2.0,1,channel(QStringLiteral("bad"),Quantity::Voltage,0.0,1.0,
                                         std::numeric_limits<double>::max(),0.0)),
        std::invalid_argument);
}

void TestStructuralTypedControl::rejectWrongIdsQuantitiesDuplicatesAndNonFiniteCommands()
{
    Session session(4,0.25,{
        channel(QStringLiteral("force"),Quantity::Force),
        channel(QStringLiteral("voltage"),Quantity::Voltage)});
    const auto state = session.initialState();
    const auto observation = session.makeObservation(
        state,{{QStringLiteral("tip"),0.0,QStringLiteral("m")}});
    Decision valid = decision(session,{1.0,2.0},QByteArrayLiteral("candidate"));

    Decision wrongId = valid;
    wrongId.requestedCommands[0].channelId = QStringLiteral("unknown");
    QVERIFY_EXCEPTION_THROWN(session.advance(state,observation,wrongId),std::invalid_argument);

    Decision wrongQuantity = valid;
    wrongQuantity.requestedCommands[1].quantity = Quantity::Charge;
    QVERIFY_EXCEPTION_THROWN(session.advance(state,observation,wrongQuantity),std::invalid_argument);

    Decision duplicate = valid;
    duplicate.requestedCommands[1] = duplicate.requestedCommands[0];
    QVERIFY_EXCEPTION_THROWN(session.advance(state,observation,duplicate),std::invalid_argument);

    Decision nonFinite = valid;
    nonFinite.requestedCommands[0].value = std::numeric_limits<double>::quiet_NaN();
    QVERIFY_EXCEPTION_THROWN(session.advance(state,observation,nonFinite),std::invalid_argument);
}

void TestStructuralTypedControl::applyRateAndAmplitudeLimitsInChannelUnits()
{
    const auto session = oneChannelSession(
        0.5,7,channel(QStringLiteral("voltage"),Quantity::Voltage,0.0,3.0,1.0,0.0));
    auto state = session.initialState();
    auto result = session.trial(
        state,{},[&session](const RStructuralControlObservation &) {
            return decision(session,{20.0});
        });
    QVERIFY(close(result.channels[0].appliedValue,0.5));
    QVERIFY(result.channels[0].appliedValue >= 0.0 && result.channels[0].appliedValue <= 3.0);

    state = result.state;
    result = session.trial(
        state,{},[&session](const RStructuralControlObservation &) {
            return decision(session,{-20.0});
        });
    QVERIFY(close(result.channels[0].appliedValue,0.0));
    QVERIFY(result.channels[0].appliedValue >= 0.0 && result.channels[0].appliedValue <= 3.0);

    state = result.state;
    result = session.trial(
        state,{},[&session](const RStructuralControlObservation &) {
            return decision(session,{20.0});
        });
    QVERIFY(close(result.channels[0].appliedValue,0.5));
}

void TestStructuralTypedControl::applyDeterministicDelayHoldAndNoise()
{
    const auto session = oneChannelSession(
        1.0,991,channel(QStringLiteral("length"),Quantity::CableReferenceLength,
                        -100.0,100.0,std::numeric_limits<double>::infinity(),
                        0.0,2,2,0.25));
    auto first = session.initialState();
    const auto run = [&session](const Session::State &accepted, double value)
    {
        return session.trial(accepted,{},[&session,value](const RStructuralControlObservation &) {
            return decision(session,{value});
        });
    };

    const auto step0 = run(first,4.0);
    QVERIFY(close(step0.channels[0].conditionedValue,0.0));
    const double held0 = step0.state.heldValues[0];
    QVERIFY(std::abs(held0-4.0) <= 0.25);
    const auto step1 = run(step0.state,9.0);
    QVERIFY(close(step1.channels[0].conditionedValue,0.0));
    QVERIFY(close(step1.state.heldValues[0],held0));
    const auto step2 = run(step1.state,-8.0);
    QVERIFY(close(step2.channels[0].conditionedValue,held0));
    const double held2 = step2.state.heldValues[0];
    const auto step3 = run(step2.state,-7.0);
    QVERIFY(close(step3.channels[0].conditionedValue,held0));
    QVERIFY(close(step3.state.heldValues[0],held2));
    const auto step4 = run(step3.state,2.0);
    QVERIFY(close(step4.channels[0].conditionedValue,held2));

    auto twin = session.initialState();
    const auto twin0 = run(twin,4.0);
    QVERIFY(close(twin0.state.heldValues[0],held0));
    QVERIFY(close(twin0.channels[0].conditionedValue,step0.channels[0].conditionedValue));
}

void TestStructuralTypedControl::retriesObserveOnlyPriorAcceptedState()
{
    const auto session = oneChannelSession();
    auto accepted = session.initialState(3.0,QByteArrayLiteral("prior"));
    std::vector<QByteArray> observedStates;
    const auto callback = [&session,&observedStates](const RStructuralControlObservation &observation)
    {
        observedStates.push_back(observation.priorControllerState);
        return decision(session,{4.0},QByteArrayLiteral("candidate"));
    };

    const auto firstAttempt = session.trial(accepted,{},callback);
    const auto retryAttempt = session.trial(accepted,{},callback);
    QCOMPARE(observedStates.size(),std::size_t(2));
    QCOMPARE(observedStates[0],QByteArrayLiteral("prior"));
    QCOMPARE(observedStates[1],QByteArrayLiteral("prior"));
    QCOMPARE(accepted.controllerState,QByteArrayLiteral("prior"));
    QCOMPARE(session.serialize(firstAttempt.state),session.serialize(retryAttempt.state));

    accepted = firstAttempt.state; // The caller publishes only after acceptance.
    const auto nextAttempt = session.trial(accepted,{},callback);
    QCOMPARE(observedStates.back(),QByteArrayLiteral("candidate"));
    QCOMPARE(nextAttempt.state.controllerState,QByteArrayLiteral("candidate"));
}

void TestStructuralTypedControl::serializedSplitTrajectoryIsExactAndCorruptionIsRejected()
{
    const auto session = oneChannelSession(
        0.125,81,channel(QStringLiteral("charge"),Quantity::Charge,-10.0,10.0,2.0,0.0,
                         3,2,0.3));
    auto full = session.initialState(1.0,QByteArrayLiteral("0"));
    auto split = full;
    for (int i=0;i<9;i++) full = takeStep(session,full,2.0+double(i));
    for (int i=0;i<4;i++) split = takeStep(session,split,2.0+double(i));

    const QByteArray checkpoint = session.serialize(split);
    split = session.deserialize(checkpoint);
    for (int i=4;i<9;i++) split = takeStep(session,split,2.0+double(i));
    QCOMPARE(session.serialize(full),session.serialize(split));

    QByteArray corrupted = checkpoint;
    corrupted[20] = char(corrupted[20] ^ 0x40);
    QVERIFY_EXCEPTION_THROWN(session.deserialize(corrupted),std::invalid_argument);
    const auto wrongConfig = oneChannelSession(
        0.125,82,channel(QStringLiteral("charge"),Quantity::Charge,-10.0,10.0,2.0,0.0,
                         3,2,0.3));
    QVERIFY_EXCEPTION_THROWN(wrongConfig.deserialize(checkpoint),std::invalid_argument);
    const auto wrongProvider = oneChannelSession(
        0.125,81,channel(QStringLiteral("charge"),Quantity::Charge,-10.0,10.0,2.0,0.0,
                         3,2,0.3),QStringLiteral("other-controller-revision"));
    QVERIFY_EXCEPTION_THROWN(wrongProvider.deserialize(checkpoint),std::invalid_argument);

    QVERIFY_EXCEPTION_THROWN(
        session.initialState(0.0,QByteArray(1024*1024+1,'x')),std::invalid_argument);
}

void TestStructuralTypedControl::rejectForgedClockAndSensorMetadata()
{
    const auto session = oneChannelSession(0.125);
    const auto initial = session.initialState(4.0);
    const auto advanced = session.trial(
        initial,{},[&session](const RStructuralControlObservation &) {
            return decision(session,{1.0});
        });

    auto forgedTime = advanced.state;
    forgedTime.acceptedTime += 0.125;
    QVERIFY_EXCEPTION_THROWN(session.serialize(forgedTime),std::invalid_argument);
    auto forgedCount = advanced.state;
    forgedCount.completedSteps++;
    QVERIFY_EXCEPTION_THROWN(session.serialize(forgedCount),std::invalid_argument);

    // A one-ULP forgery must not fit an accumulated clock tolerance at scale.
    auto farState=initial;
    farState.completedSteps=1000000000;
    farState.acceptedTime=std::fma(double(farState.completedSteps),session.fixedDt(),farState.initialTime);
    QVERIFY(!session.serialize(farState).isEmpty());
    farState.acceptedTime=std::nextafter(farState.acceptedTime,std::numeric_limits<double>::infinity());
    QVERIFY_EXCEPTION_THROWN(session.serialize(farState),std::invalid_argument);
    const auto fineSession=oneChannelSession(0.1);
    const auto hugeInitial=fineSession.initialState(1.0e12);
    QVERIFY_EXCEPTION_THROWN(fineSession.makeObservation(hugeInitial,{}),std::invalid_argument);

    QVERIFY_EXCEPTION_THROWN(
        session.makeObservation(initial,{
            {QStringLiteral("tip"),0.0,QStringLiteral("m")},
            {QStringLiteral("tip"),1.0,QStringLiteral("m")}}),
        std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(
        session.makeObservation(initial,{
            {QStringLiteral("tip"),0.0,QString(33,QLatin1Char('m'))}}),
        std::invalid_argument);
}

void TestStructuralTypedControl::accountForDriverSuppliedWorkPerChannel()
{
    const auto session = oneChannelSession(0.25);
    const auto accepted = session.initialState();
    const auto candidate = session.trial(
        accepted,{},[&session](const RStructuralControlObservation &) {
            return decision(session,{1.0},QByteArrayLiteral("next"));
        });
    const auto worked = session.withWork(candidate,{2.0});
    QVERIFY(close(worked.channels[0].cumulativeWorkJ,2.0));
    QVERIFY(close(worked.channels[0].averagePowerW,8.0));
    QCOMPARE(candidate.state.cumulativeWorkJ[0],0.0);
    QCOMPARE(accepted.cumulativeWorkJ[0],0.0);

    const auto regenerated = session.withWork(worked,{-0.5});
    QVERIFY(close(regenerated.channels[0].cumulativeWorkJ,1.5));
    QVERIFY(close(regenerated.channels[0].averagePowerW,-2.0));
    QVERIFY_EXCEPTION_THROWN(session.withWork(candidate,{}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(
        session.withWork(candidate,{std::numeric_limits<double>::infinity()}),std::invalid_argument);
}

QTEST_APPLESS_MAIN(TestStructuralTypedControl)

#include "tst_structural_typed_control.moc"

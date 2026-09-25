#include <QtTest>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "rstructural_thermal_control.h"

namespace
{
RStructuralThermalControlConfig config(unsigned int thermalSubsteps = 1,
                                       unsigned int sampleEveryMacrosteps = 1,
                                       double minimumPowerWatts = 0.0,
                                       double maximumPowerSlewWattsPerSecond =
                                           std::numeric_limits<double>::infinity())
{
    RStructuralThermalControlConfig result;
    result.thermalSubstepSeconds = 1.0;
    result.thermalSubstepsPerMechanicalMacrostep = thermalSubsteps;
    result.mechanicalMacrostepsPerControllerSample = sampleEveryMacrosteps;
    result.channels = {{"heater-A","node-A","sensor-A","node-A",
                        302.999,303.0,10.0,minimumPowerWatts,10.0,
                        maximumPowerSlewWattsPerSecond}};
    result.sparseObservationNodeIds = {"node-A"};
    return result;
}

RStructuralThermalControlScheduler::Plant initialPlant()
{
    // A single lumped thermal state with C=10 J/K is the scheduler contract
    // fixture; it does not represent a finite-element thermal integration.
    return {{"node-A",300.0,0.0}};
}

RStructuralThermalControlScheduler::SubstepFunction lumpedBackend(
    unsigned int *callCount = nullptr,
    int failOnCall = -1)
{
    return [callCount,failOnCall](const auto &accepted,
                                  const auto &heaters,
                                  double dt,
                                  bool solveMechanical) {
        if (callCount) ++*callCount;
        if (callCount && int(*callCount) == failOnCall)
            throw std::runtime_error("injected thermal backend failure");
        if (heaters.size() != 1 || heaters.front().actuatorId != "heater-A")
            throw std::runtime_error("unexpected heater routing");
        auto candidate = accepted;
        candidate.front().temperatureKelvin += heaters.front().watts*dt/10.0;
        if (solveMechanical)
            candidate.front().axialDisplacement = (candidate.front().temperatureKelvin-300.0)*1.0e-6;
        return RStructuralThermalControlSubstepResult{candidate,heaters};
    };
}

bool sameState(const RStructuralThermalControlState &left,
               const RStructuralThermalControlState &right)
{
    if (left.configurationSha256 != right.configurationSha256 ||
        left.completedMechanicalMacrosteps != right.completedMechanicalMacrosteps ||
        left.completedControllerSamples != right.completedControllerSamples ||
        left.plant.size() != right.plant.size() || left.channels.size() != right.channels.size())
        return false;
    for (std::size_t i=0;i<left.plant.size();i++)
    {
        if (left.plant[i].nodeId != right.plant[i].nodeId ||
            left.plant[i].temperatureKelvin != right.plant[i].temperatureKelvin ||
            left.plant[i].axialDisplacement != right.plant[i].axialDisplacement)
            return false;
    }
    for (std::size_t i=0;i<left.channels.size();i++)
    {
        const auto &a = left.channels[i];
        const auto &b = right.channels[i];
        if (a.actuatorId != b.actuatorId || a.heaterOn != b.heaterOn ||
            a.requestedPowerWatts != b.requestedPowerWatts ||
            a.conditionedPowerWatts != b.conditionedPowerWatts ||
            a.appliedPowerWatts != b.appliedPowerWatts ||
            a.cumulativeRequestedEnergyJoules != b.cumulativeRequestedEnergyJoules ||
            a.cumulativeConditionedEnergyJoules != b.cumulativeConditionedEnergyJoules ||
            a.cumulativeAppliedEnergyJoules != b.cumulativeAppliedEnergyJoules)
            return false;
    }
    return true;
}
}

class TestStructuralThermalControl : public QObject
{
    Q_OBJECT

private slots:
    void bangBangAcceptsThreeJouleBalancedSteps();
    void positiveMinimumPowerAllowsOffAndRestart();
    void slewTransitionsMayPassBelowOnMinimum();
    void restartReplayIsDeterministic();
    void integerRatesHoldCommandsAndSolveMechanicalAtMacroBoundary();
    void failedSubstepPublishesNoCandidate();
    void checkpointRejectsTruncationCorruptionAndDifferentConfiguration();
    void observationsUseStableNodeAndSensorIds();
};

void TestStructuralThermalControl::bangBangAcceptsThreeJouleBalancedSteps()
{
    const RStructuralThermalControlScheduler scheduler(config());
    auto state = scheduler.initialState(initialPlant());
    double totalJoules = 0.0;
    for (int i=0;i<3;i++)
    {
        const auto result = scheduler.advance(state,lumpedBackend());
        QCOMPARE(result.candidate.completedMechanicalMacrosteps,quint64(i+1));
        QCOMPARE(result.controllerSampled,true);
        QCOMPARE(result.channels.size(),std::size_t(1));
        QCOMPARE(result.channels[0].requestedPowerWatts,10.0);
        QCOMPARE(result.channels[0].conditionedPowerWatts,10.0);
        QCOMPARE(result.channels[0].appliedAveragePowerWatts,10.0);
        QCOMPARE(result.channels[0].requestedEnergyJoules,10.0);
        QCOMPARE(result.channels[0].conditionedEnergyJoules,10.0);
        QCOMPARE(result.channels[0].appliedEnergyJoules,10.0);
        totalJoules += result.channels[0].appliedEnergyJoules;
        state = result.candidate;
    }
    QCOMPARE(state.plant[0].temperatureKelvin,303.0);
    QCOMPARE(totalJoules,30.0);
    QCOMPARE(state.channels[0].cumulativeRequestedEnergyJoules,30.0);
    QCOMPARE(state.channels[0].cumulativeConditionedEnergyJoules,30.0);
    QCOMPARE(state.channels[0].cumulativeAppliedEnergyJoules,30.0);

    const auto shutoff = scheduler.advance(state,lumpedBackend());
    QCOMPARE(shutoff.controllerSampled,true);
    QCOMPARE(shutoff.channels[0].requestedPowerWatts,0.0);
    QCOMPARE(shutoff.candidate.plant[0].temperatureKelvin,303.0);
}

void TestStructuralThermalControl::positiveMinimumPowerAllowsOffAndRestart()
{
    const RStructuralThermalControlScheduler scheduler(config(1,1,4.0));
    auto state = scheduler.initialState(initialPlant());
    QCOMPARE(state.channels[0].requestedPowerWatts,0.0);
    QCOMPARE(state.channels[0].conditionedPowerWatts,0.0);
    QCOMPARE(state.channels[0].appliedPowerWatts,0.0);

    for (int i=0;i<3;i++)
        state = scheduler.advance(state,lumpedBackend()).candidate;
    QCOMPARE(state.plant[0].temperatureKelvin,303.0);
    QCOMPARE(state.channels[0].conditionedPowerWatts,10.0);

    const auto restoredBeforeShutoff = scheduler.deserialize(scheduler.serialize(state));
    const auto shutoff = scheduler.advance(restoredBeforeShutoff,lumpedBackend());
    QCOMPARE(shutoff.channels[0].requestedPowerWatts,0.0);
    QCOMPARE(shutoff.channels[0].conditionedPowerWatts,0.0);
    QCOMPARE(shutoff.channels[0].appliedAveragePowerWatts,0.0);
    QCOMPARE(shutoff.channels[0].requestedEnergyJoules,0.0);
    QCOMPARE(shutoff.channels[0].conditionedEnergyJoules,0.0);
    QCOMPARE(shutoff.channels[0].appliedEnergyJoules,0.0);
    QCOMPARE(shutoff.candidate.plant[0].temperatureKelvin,303.0);

    const auto restoredOff = scheduler.deserialize(scheduler.serialize(shutoff.candidate));
    const auto heldOff = scheduler.advance(restoredOff,lumpedBackend());
    QCOMPARE(heldOff.channels[0].requestedPowerWatts,0.0);
    QCOMPARE(heldOff.channels[0].conditionedPowerWatts,0.0);
    QCOMPARE(heldOff.channels[0].appliedAveragePowerWatts,0.0);
    QCOMPARE(heldOff.candidate.channels[0].cumulativeAppliedEnergyJoules,30.0);
}

void TestStructuralThermalControl::slewTransitionsMayPassBelowOnMinimum()
{
    const RStructuralThermalControlScheduler scheduler(config(1,1,4.0,2.0));
    auto state = scheduler.initialState(initialPlant());
    const auto rampingOn = scheduler.advance(state,lumpedBackend());
    QCOMPARE(rampingOn.channels[0].requestedPowerWatts,10.0);
    QCOMPARE(rampingOn.channels[0].conditionedPowerWatts,2.0);
    QCOMPARE(rampingOn.channels[0].appliedAveragePowerWatts,2.0);
    state = rampingOn.candidate;
    state = scheduler.advance(state,lumpedBackend()).candidate;
    QCOMPARE(state.channels[0].conditionedPowerWatts,4.0);
    state = scheduler.advance(state,lumpedBackend()).candidate;
    QCOMPARE(state.channels[0].conditionedPowerWatts,6.0);

    // Move the sensed plant across the switch threshold while retaining a
    // valid accepted controller state, then verify an off slew can descend
    // through the positive on-target minimum.
    state.plant[0].temperatureKelvin = 303.0;
    const auto rampingOff = scheduler.advance(state,lumpedBackend());
    QCOMPARE(rampingOff.channels[0].requestedPowerWatts,0.0);
    QCOMPARE(rampingOff.channels[0].conditionedPowerWatts,4.0);
    QCOMPARE(rampingOff.channels[0].appliedAveragePowerWatts,4.0);
    const auto belowMinimumRampdown = scheduler.advance(rampingOff.candidate,lumpedBackend());
    QCOMPARE(belowMinimumRampdown.channels[0].requestedPowerWatts,0.0);
    QCOMPARE(belowMinimumRampdown.channels[0].conditionedPowerWatts,2.0);
    QCOMPARE(belowMinimumRampdown.channels[0].appliedAveragePowerWatts,2.0);
}

void TestStructuralThermalControl::restartReplayIsDeterministic()
{
    const RStructuralThermalControlScheduler scheduler(config(2,2));
    auto uninterrupted = scheduler.initialState(initialPlant());
    const auto first = scheduler.advance(uninterrupted,lumpedBackend());
    const QByteArray checkpoint = scheduler.serialize(first.candidate);
    const auto restored = scheduler.deserialize(checkpoint);
    QVERIFY(sameState(restored,first.candidate));

    const auto directSecond = scheduler.advance(first.candidate,lumpedBackend());
    const auto replaySecond = scheduler.advance(restored,lumpedBackend());
    QVERIFY(sameState(directSecond.candidate,replaySecond.candidate));
    QCOMPARE(directSecond.channels[0].requestedPowerWatts,
             replaySecond.channels[0].requestedPowerWatts);
    QCOMPARE(directSecond.channels[0].appliedEnergyJoules,
             replaySecond.channels[0].appliedEnergyJoules);
}

void TestStructuralThermalControl::integerRatesHoldCommandsAndSolveMechanicalAtMacroBoundary()
{
    const RStructuralThermalControlScheduler scheduler(config(3,2));
    auto state = scheduler.initialState(initialPlant());
    std::vector<bool> solveFlags;
    auto backend = [&solveFlags](const auto &accepted,const auto &heaters,double dt,bool solveMechanical) {
        solveFlags.push_back(solveMechanical);
        auto candidate = accepted;
        candidate.front().temperatureKelvin += heaters.front().watts*dt/10.0;
        if (solveMechanical)
            candidate.front().axialDisplacement = (candidate.front().temperatureKelvin-300.0)*1.0e-6;
        return RStructuralThermalControlSubstepResult{candidate,heaters};
    };

    const auto first = scheduler.advance(state,backend);
    QCOMPARE(first.controllerSampled,true);
    QCOMPARE(first.channels[0].requestedEnergyJoules,30.0);
    state = first.candidate;
    const auto second = scheduler.advance(state,backend);
    QCOMPARE(second.controllerSampled,false);
    QCOMPARE(second.channels[0].requestedPowerWatts,10.0);
    QCOMPARE(second.channels[0].requestedEnergyJoules,30.0);
    QCOMPARE(second.candidate.completedControllerSamples,quint64(1));
    QCOMPARE(solveFlags.size(),std::size_t(6));
    QCOMPARE(solveFlags[0],false);
    QCOMPARE(solveFlags[1],false);
    QCOMPARE(solveFlags[2],true);
    QCOMPARE(solveFlags[3],false);
    QCOMPARE(solveFlags[4],false);
    QCOMPARE(solveFlags[5],true);
    QCOMPARE(second.candidate.plant[0].temperatureKelvin,306.0);
}

void TestStructuralThermalControl::failedSubstepPublishesNoCandidate()
{
    const RStructuralThermalControlScheduler scheduler(config(3,1));
    const auto accepted = scheduler.initialState(initialPlant());
    const QByteArray before = scheduler.serialize(accepted);
    unsigned int callbackCount = 0;
    QVERIFY_EXCEPTION_THROWN(scheduler.advance(accepted,lumpedBackend(&callbackCount,2)),
                             std::runtime_error);
    QCOMPARE(callbackCount,2u);
    QCOMPARE(scheduler.serialize(accepted),before);
}

void TestStructuralThermalControl::checkpointRejectsTruncationCorruptionAndDifferentConfiguration()
{
    const RStructuralThermalControlScheduler scheduler(config());
    const auto state = scheduler.initialState(initialPlant());
    const QByteArray checkpoint = scheduler.serialize(state);
    QVERIFY_EXCEPTION_THROWN(scheduler.deserialize(checkpoint.left(checkpoint.size()-1)),
                             std::invalid_argument);
    QByteArray corrupted = checkpoint;
    corrupted[0] = char(corrupted[0] ^ 0x01);
    QVERIFY_EXCEPTION_THROWN(scheduler.deserialize(corrupted),std::invalid_argument);
    const RStructuralThermalControlScheduler otherScheduler(config(2,1));
    QVERIFY_EXCEPTION_THROWN(otherScheduler.deserialize(checkpoint),std::invalid_argument);
}

void TestStructuralThermalControl::observationsUseStableNodeAndSensorIds()
{
    const RStructuralThermalControlScheduler scheduler(config());
    auto reversed = initialPlant();
    reversed.push_back({"node-B",305.0,0.002});
    reversed[0].nodeId = "node-B";
    reversed[0].temperatureKelvin = 305.0;
    reversed[0].axialDisplacement = 0.002;
    reversed[1] = {"node-A",300.0,0.0};
    auto testConfig = config();
    testConfig.sparseObservationNodeIds = {"node-B"};
    RStructuralThermalControlScheduler twoNodeScheduler(testConfig);
    auto state = twoNodeScheduler.initialState(reversed);
    const auto result = twoNodeScheduler.advance(state,[](const auto &accepted,const auto &heaters,
                                                          double,bool) {
        return RStructuralThermalControlSubstepResult{accepted,heaters};
    });
    QCOMPARE(result.observation.fullField.size(),std::size_t(2));
    QCOMPARE(result.observation.fullField[0].nodeId,QString("node-A"));
    QCOMPARE(result.observation.fullField[1].nodeId,QString("node-B"));
    QCOMPARE(result.observation.sparseField.size(),std::size_t(1));
    QCOMPARE(result.observation.sparseField[0].nodeId,QString("node-B"));
    QCOMPARE(result.observation.sparseField[0].temperatureKelvin,305.0);
    QCOMPARE(result.observation.sparseField[0].axialDisplacementMetres,0.002);
    QCOMPARE(result.observation.temperatureSensors.size(),std::size_t(1));
    QCOMPARE(result.observation.temperatureSensors[0].sensorId,QString("sensor-A"));
    QCOMPARE(result.observation.temperatureSensors[0].nodeId,QString("node-A"));
}

QTEST_APPLESS_MAIN(TestStructuralThermalControl)
#include "tst_structural_thermal_control.moc"

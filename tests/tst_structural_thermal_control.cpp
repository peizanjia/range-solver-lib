#include <QtTest>
#include <QTemporaryDir>
#include <QFile>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <rbl_error.h>
#include <rml_model.h>
#include <rml_structural_interaction.h>

#include "rstructural_thermal_control.h"
#include "rstructural_thermal_fe_adapter.h"
#include "rsolver.h"

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
        left.acceptedTimeSeconds != right.acceptedTimeSeconds ||
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

void setComponent(RBoundaryCondition &condition,RVariableType type,double value)
{
    const uint position = condition.findComponentPosition(type);
    if (position == RConstants::eod) throw std::logic_error("missing boundary component");
    RConditionComponent &component = condition.getComponent(position);
    component.setEnabled(true);
    component.clear();
    component.add(0.0,value);
}

RModel thermalFeRod()
{
    RModel model;
    model.setNNodes(2);
    model.getNode(0).set(0.0,0.0,0.0);
    model.getNode(1).set(1.0,0.0,0.0);
    RElement rod(R_ELEMENT_TRUSS1);
    rod.setNodeId(0,0);
    rod.setNodeId(1,1);
    model.addElement(rod,true,0);
    RLine &line = model.getLine(0);
    line.setCrossArea(1.0);
    RMaterial material(RMaterial::Solid);
    const auto property = [&material](RMaterialProperty::Type type,double value) {
        RMaterialProperty entry(type);
        entry.add(293.15,value);
        material.add(entry);
    };
    property(RMaterialProperty::ModulusOfElasiticity,10.0);
    property(RMaterialProperty::PoissonRatio,0.0);
    property(RMaterialProperty::ThermalExpansionCoefficient,0.01);
    property(RMaterialProperty::ThermalConductivity,1.0);
    property(RMaterialProperty::HeatCapacity,3.0);
    property(RMaterialProperty::Density,1.0);
    property(RMaterialProperty::Emissivity,0.0);
    line.setMaterial(material);
    for (uint node=0;node<2;node++)
    {
        RElement point(R_ELEMENT_POINT);
        point.setNodeId(0,node);
        model.addElement(point,true,node);
        model.getPoint(node).setVolume(0.0);
        RBoundaryCondition support(R_BOUNDARY_CONDITION_DISPLACEMENT);
        if (node == 0) setComponent(support,R_VARIABLE_DISPLACEMENT_X,0.0);
        else support.getComponent(support.findComponentPosition(
                 R_VARIABLE_DISPLACEMENT_X)).setEnabled(false);
        setComponent(support,R_VARIABLE_DISPLACEMENT_Y,0.0);
        setComponent(support,R_VARIABLE_DISPLACEMENT_Z,0.0);
        model.getPoint(node).addBoundaryCondition(support);
    }
    RBoundaryCondition anchorTemperature(R_BOUNDARY_CONDITION_TEMPERATURE);
    setComponent(anchorTemperature,R_VARIABLE_TEMPERATURE,293.15);
    model.getPoint(0).addBoundaryCondition(anchorTemperature);
    RProblemTaskItem sequence;
    sequence.addChild(RProblemTaskItem(R_PROBLEM_HEAT));
    sequence.addChild(RProblemTaskItem(R_PROBLEM_STRESS));
    model.setProblemTaskTree(sequence);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setSolverCvgValue(1.0e-12);

    auto &catalog = model.getProblemSetup().getStructuralInteractionSetup();
    RStructuralActuator actuator;
    actuator.id = "heater-fe";
    actuator.type = RStructuralActuatorType::HeatSource;
    actuator.target.scope = RStructuralTargetScope::Node;
    actuator.target.nodeId = 1;
    actuator.variable = R_VARIABLE_HEAT;
    actuator.unit = "W";
    actuator.minimum = 0.0;
    actuator.maximum = 2.0;
    catalog.getActuators().push_back(actuator);
    RStructuralSensor sensor;
    sensor.id = "sensor-fe";
    sensor.target.scope = RStructuralTargetScope::Node;
    sensor.target.nodeId = 1;
    sensor.quantity = RStructuralSensorQuantity::Temperature;
    sensor.aggregation = RStructuralSensorAggregation::Point;
    sensor.unit = "K";
    catalog.getSensors().push_back(sensor);
    return model;
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
    void pointHeaterControlsTransientHeatAndStressFE();
    void legacyPointHeatBoundaryActuallyAppliesWatts();
    void pointThermalMassActuallyStoresHeat();
    void multirateFeSubstepsPreserveTemperature();
    void versionedModelSetupDrivesActualFe();
};

void TestStructuralThermalControl::versionedModelSetupDrivesActualFe()
{
    RModel model = thermalFeRod();
    const QString node = RStructuralThermalFEAdapter::nodeId(model,1);
    auto &setup = model.getProblemSetup().getStructuralThermalControlSetup();
    setup.enabled = true;
    setup.thermalSubstepSeconds = 1.0;
    setup.thermalSubstepsPerMechanicalMacrostep = 1;
    setup.mechanicalMacrostepsPerControllerSample = 1;
    RStructuralThermalControlChannelSetup channel;
    channel.actuatorId = "heater-fe";
    channel.heaterNodeId = node;
    channel.temperatureSensorId = "sensor-fe";
    channel.temperatureSensorNodeId = node;
    channel.turnOnBelowKelvin = 294.15;
    channel.turnOffAtOrAboveKelvin = 294.15;
    channel.onPowerWatts = 2.0;
    channel.maximumPowerWatts = 2.0;
    setup.channels = {channel};
    setup.sparseObservationNodeIds = {node};
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString file = directory.filePath(
        QStringLiteral("configured.") + RModel::getDefaultFileExtension(true));
    RModel loaded;
    try { model.write(file,false); loaded.read(file); }
    catch (const RError &error) { QFAIL(qPrintable(error.getMessage())); }
    catch (const std::exception &error) { QFAIL(error.what()); }
    try
    {
        const auto scheduler = RStructuralThermalFEAdapter::schedulerFromModel(loaded);
        const auto accepted = scheduler.initialState(RStructuralThermalFEAdapter::capture(loaded));
        const auto result = RStructuralThermalFEAdapter::advance(loaded,scheduler,accepted);
        QCOMPARE(result.channels[0].appliedEnergyJoules,2.0);
        QVERIFY(std::abs(loaded.getVariable(loaded.findVariable(R_VARIABLE_TEMPERATURE))
                             .getValue(0,1)-294.15) < 1.0e-8);
        const QString checkpoint = directory.filePath(QStringLiteral("configured.tfec"));
        RStructuralThermalFEAdapter::writeCheckpoint(
            checkpoint,loaded,scheduler,result.candidate);
        RModel restored;
        RStructuralThermalControlState restoredState;
        RStructuralThermalFEAdapter::readConfiguredCheckpoint(
            checkpoint,restored,restoredState);
        QVERIFY(sameState(restoredState,result.candidate));
        RModel conflicting(loaded);
        conflicting.getProblemSetup().getStructuralThermalControlSetup()
            .channels[0].turnOnBelowKelvin = 294.0;
        bool mismatchedSetupRejected = false;
        try { RStructuralThermalFEAdapter::writeCheckpoint(
                  directory.filePath(QStringLiteral("wrong.tfec")),
                  conflicting,scheduler,result.candidate); }
        catch (const RError &) { mismatchedSetupRejected = true; }
        catch (const std::exception &) { mismatchedSetupRejected = true; }
        QVERIFY(mismatchedSetupRejected);
        const auto restoredScheduler = RStructuralThermalFEAdapter::schedulerFromModel(restored);
        const auto uninterrupted = RStructuralThermalFEAdapter::advance(
            loaded,scheduler,result.candidate);
        const auto resumed = RStructuralThermalFEAdapter::advance(
            restored,restoredScheduler,restoredState);
        QVERIFY(sameState(uninterrupted.candidate,resumed.candidate));
    }
    catch (const RError &error) { QFAIL(qPrintable(error.getMessage())); }
    catch (const std::exception &error) { QFAIL(error.what()); }
}

void TestStructuralThermalControl::multirateFeSubstepsPreserveTemperature()
{
    RModel model = thermalFeRod();
    RStructuralThermalControlConfig settings;
    settings.thermalSubstepSeconds = 1.0;
    settings.thermalSubstepsPerMechanicalMacrostep = 2;
    settings.mechanicalMacrostepsPerControllerSample = 1;
    const QString node = RStructuralThermalFEAdapter::nodeId(model,1);
    settings.channels = {{"heater-fe",node,"sensor-fe",node,
                          400.0,401.0,2.0,0.0,2.0}};
    settings.sparseObservationNodeIds = {node};
    const RStructuralThermalControlScheduler scheduler(settings);
    const auto initial = scheduler.initialState(RStructuralThermalFEAdapter::capture(model));
    RStructuralThermalControlScheduler::Result result;
    try { result = RStructuralThermalFEAdapter::advance(model,scheduler,initial); }
    catch (const RError &error) { QFAIL(qPrintable(error.getMessage())); }
    catch (const std::exception &error) { QFAIL(error.what()); }
    // Backward Euler with the rod row C/dt=K=1 and 2 W input yields
    // theta_1=1 K, theta_2=1.5 K before the single mechanical solve.
    QVERIFY(std::abs(model.getVariable(model.findVariable(R_VARIABLE_TEMPERATURE)).getValue(0,1)-
                     294.65) < 1.0e-8);
    // The two-node rod uses mean element temperature, so its elastic free
    // elongation is alpha * (1.5 K / 2) * L = 0.0075 m.
    QVERIFY(std::abs(model.getVariable(model.findVariable(R_VARIABLE_DISPLACEMENT)).getValue(0,1)-
                     0.0075) < 1.0e-8);
    QCOMPARE(result.channels[0].appliedEnergyJoules,4.0);
    QCOMPARE(result.candidate.completedMechanicalMacrosteps,quint64(1));
    QCOMPARE(model.getTimeSolver().getComputedTime(),2.0);
}

void TestStructuralThermalControl::pointThermalMassActuallyStoresHeat()
{
    RModel model = thermalFeRod();
    model.getPoint(1).setVolume(1.0);
    model.getPoint(1).setMaterial(model.getLine(0).getMaterial());
    RBoundaryCondition heater(R_BOUNDARY_CONDITION_HEAT);
    setComponent(heater,R_VARIABLE_HEAT,2.0);
    model.getPoint(1).addBoundaryCondition(heater);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_HEAT));
    RTimeSolver &time = model.getTimeSolver();
    time.setEnabled(true);
    time.setInputNTimeSteps(1);
    time.setInputStartTime(0.0);
    time.setInputTimeStepSize(1.0);
    time.setTimeMarchApproximation(R_TIME_MARCH_BACKWARD);
    RSolver solver(model,QString(),QString());
    try { solver.run(); }
    catch (const RError &error) { QFAIL(qPrintable(error.getMessage())); }
    catch (const std::exception &error) { QFAIL(error.what()); }
    // The base rod row has C/dt + K = 2 W/K; the added point contributes
    // rho*c*V/dt = 3 W/K. Two watts therefore raise the free node by 0.4 K.
    QVERIFY(std::abs(model.getVariable(model.findVariable(R_VARIABLE_TEMPERATURE)).getValue(0,1)-
                     293.55) < 1.0e-8);
}

void TestStructuralThermalControl::legacyPointHeatBoundaryActuallyAppliesWatts()
{
    RModel model = thermalFeRod();
    RBoundaryCondition heater(R_BOUNDARY_CONDITION_HEAT);
    setComponent(heater,R_VARIABLE_HEAT,2.0);
    model.getPoint(1).addBoundaryCondition(heater);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_HEAT));
    RTimeSolver &time = model.getTimeSolver();
    time.setEnabled(true);
    time.setInputNTimeSteps(1);
    time.setInputStartTime(0.0);
    time.setInputTimeStepSize(1.0);
    time.setTimeMarchApproximation(R_TIME_MARCH_BACKWARD);
    RSolver solver(model,QString(),QString());
    try { solver.run(); }
    catch (const RError &error) { QFAIL(qPrintable(error.getMessage())); }
    catch (const std::exception &error) { QFAIL(error.what()); }
    QVERIFY(std::abs(model.getVariable(model.findVariable(R_VARIABLE_TEMPERATURE)).getValue(0,1)-
                     294.15) < 1.0e-8);
}

void TestStructuralThermalControl::pointHeaterControlsTransientHeatAndStressFE()
{
    RModel model = thermalFeRod();
    RStructuralThermalControlConfig settings;
    settings.thermalSubstepSeconds = 1.0;
    settings.thermalSubstepsPerMechanicalMacrostep = 1;
    settings.mechanicalMacrostepsPerControllerSample = 1;
    const QString heatedNode = RStructuralThermalFEAdapter::nodeId(model,1);
    settings.channels = {{"heater-fe",heatedNode,"sensor-fe",heatedNode,
                          294.15,294.15,2.0,0.0,2.0}};
    settings.sparseObservationNodeIds = {heatedNode};
    const RStructuralThermalControlScheduler scheduler(settings);
    auto accepted = scheduler.initialState(RStructuralThermalFEAdapter::capture(model));
    RStructuralThermalControlScheduler::Result first;
    try { first = RStructuralThermalFEAdapter::advance(model,scheduler,accepted); }
    catch (const RError &error) { QFAIL(qPrintable(error.getMessage())); }
    catch (const std::exception &error) { QFAIL(error.what()); }
    QCOMPARE(first.channels.size(),std::size_t(1));
    QCOMPARE(first.channels[0].appliedEnergyJoules,2.0);
    QVERIFY(std::abs(model.getVariable(model.findVariable(R_VARIABLE_TEMPERATURE)).getValue(0,1)-
                     294.15) < 1.0e-8);
    QVERIFY(std::abs(model.getVariable(model.findVariable(R_VARIABLE_DISPLACEMENT)).getValue(0,1)-
                     0.005) < 1.0e-8);
    accepted = first.candidate;
    RStructuralThermalControlScheduler::Result second;
    try { second = RStructuralThermalFEAdapter::advance(model,scheduler,accepted); }
    catch (const RError &error) { QFAIL(qPrintable(error.getMessage())); }
    catch (const std::exception &error) { QFAIL(error.what()); }
    QCOMPARE(second.channels[0].appliedEnergyJoules,0.0);
    // Conduction uses the accepted displaced length (1.005 m), while the
    // mass-conserving rod heat capacity remains at its reference value.
    const double deformedLength = 1.005;
    const double expectedCooldown = 293.15 +
        1.0/(1.0+1.0/deformedLength);
    QVERIFY(std::abs(model.getVariable(model.findVariable(R_VARIABLE_TEMPERATURE)).getValue(0,1)-
                     expectedCooldown) < 1.0e-6);
    QCOMPARE(second.candidate.channels[0].cumulativeAppliedEnergyJoules,2.0);
    auto wrongClock = second.candidate;
    wrongClock.acceptedTimeSeconds += 1.0;
    bool clockMismatchRejected = false;
    try { (void)RStructuralThermalFEAdapter::advance(model,scheduler,wrongClock); }
    catch (const std::exception &) { clockMismatchRejected = true; }
    QVERIFY(clockMismatchRejected);

    RModel resumedModel(model);
    const QByteArray serializedState = scheduler.serialize(second.candidate);
    QTemporaryDir checkpointDirectory;
    QVERIFY(checkpointDirectory.isValid());
    const QString checkpointFile = checkpointDirectory.filePath(
        QStringLiteral("thermal-fe.") + RModel::getDefaultFileExtension(true));
    RModel diskResumedModel;
    try
    {
        QCOMPARE(model.write(checkpointFile,false),checkpointFile);
        diskResumedModel.read(checkpointFile);
    }
    catch (const RError &error) { QFAIL(qPrintable(error.getMessage())); }
    catch (const std::exception &error) { QFAIL(error.what()); }
    QCOMPARE(diskResumedModel.getTimeSolver().getComputedTime(),
             model.getTimeSolver().getComputedTime());
    const auto diskState = scheduler.deserialize(serializedState);
    RStructuralThermalControlScheduler::Result diskContinuation;
    try { diskContinuation = RStructuralThermalFEAdapter::advance(
              diskResumedModel,scheduler,diskState); }
    catch (const RError &error) { QFAIL(qPrintable(error.getMessage())); }
    catch (const std::exception &error) { QFAIL(error.what()); }
    const QString bundleFile = checkpointDirectory.filePath(QStringLiteral("accepted.tfec"));
    RModel bundledModel;
    RStructuralThermalControlState bundledState;
    try
    {
        RStructuralThermalFEAdapter::writeCheckpoint(bundleFile,model,scheduler,second.candidate);
        RStructuralThermalFEAdapter::readCheckpoint(
            bundleFile,bundledModel,scheduler,bundledState);
    }
    catch (const RError &error) { QFAIL(qPrintable(error.getMessage())); }
    catch (const std::exception &error) { QFAIL(error.what()); }
    QVERIFY(sameState(bundledState,second.candidate));
    RStructuralThermalControlScheduler::Result bundledContinuation;
    try { bundledContinuation = RStructuralThermalFEAdapter::advance(
              bundledModel,scheduler,bundledState); }
    catch (const RError &error) { QFAIL(qPrintable(error.getMessage())); }
    catch (const std::exception &error) { QFAIL(error.what()); }

    const double acceptedTemperature =
        model.getVariable(model.findVariable(R_VARIABLE_TEMPERATURE)).getValue(0,1);
    const double acceptedDisplacement =
        model.getVariable(model.findVariable(R_VARIABLE_DISPLACEMENT)).getValue(0,1);
    const double acceptedTime = model.getTimeSolver().getComputedTime();
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setSolverCvgValue(-1.0);
    bool failed = false;
    try { (void)RStructuralThermalFEAdapter::advance(model,scheduler,second.candidate); }
    catch (const RError &) { failed = true; }
    catch (const std::exception &) { failed = true; }
    QVERIFY(failed);
    QCOMPARE(model.getVariable(model.findVariable(R_VARIABLE_TEMPERATURE)).getValue(0,1),
             acceptedTemperature);
    QCOMPARE(model.getVariable(model.findVariable(R_VARIABLE_DISPLACEMENT)).getValue(0,1),
             acceptedDisplacement);
    QCOMPARE(model.getTimeSolver().getComputedTime(),acceptedTime);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setSolverCvgValue(1.0e-12);
    RStructuralThermalControlScheduler::Result retry;
    try { retry = RStructuralThermalFEAdapter::advance(model,scheduler,second.candidate); }
    catch (const RError &error) { QFAIL(qPrintable(error.getMessage())); }
    catch (const std::exception &error) { QFAIL(error.what()); }
    QCOMPARE(retry.candidate.completedMechanicalMacrosteps,quint64(3));
    QCOMPARE(retry.channels[0].appliedEnergyJoules,2.0);
    QVERIFY(sameState(diskContinuation.candidate,retry.candidate));
    QVERIFY(sameState(bundledContinuation.candidate,retry.candidate));
    QFile corruptBundle(bundleFile);
    QVERIFY(corruptBundle.open(QIODevice::ReadWrite));
    QVERIFY(corruptBundle.seek(8));
    const QByteArray originalByte = corruptBundle.read(1);
    QCOMPARE(originalByte.size(),1);
    QVERIFY(corruptBundle.seek(8));
    QCOMPARE(corruptBundle.write(QByteArray(1,char(originalByte[0] ^ 1))),qint64(1));
    corruptBundle.close();
    const auto intactBundledState = bundledState;
    const double intactBundledTemperature = bundledModel.getVariable(
        bundledModel.findVariable(R_VARIABLE_TEMPERATURE)).getValue(0,1);
    bool rejectedCorruption = false;
    try { RStructuralThermalFEAdapter::readCheckpoint(
              bundleFile,bundledModel,scheduler,bundledState); }
    catch (const std::exception &) { rejectedCorruption = true; }
    QVERIFY(rejectedCorruption);
    QVERIFY(sameState(bundledState,intactBundledState));
    QCOMPARE(bundledModel.getVariable(
                 bundledModel.findVariable(R_VARIABLE_TEMPERATURE)).getValue(0,1),
             intactBundledTemperature);

    const auto restored = scheduler.deserialize(serializedState);
    RStructuralThermalControlScheduler::Result resumed;
    try { resumed = RStructuralThermalFEAdapter::advance(resumedModel,scheduler,restored); }
    catch (const RError &error) { QFAIL(qPrintable(error.getMessage())); }
    catch (const std::exception &error) { QFAIL(error.what()); }
    QCOMPARE(resumed.candidate.channels[0].cumulativeAppliedEnergyJoules,
             retry.candidate.channels[0].cumulativeAppliedEnergyJoules);
    QCOMPARE(resumedModel.getVariable(resumedModel.findVariable(R_VARIABLE_TEMPERATURE)).getValue(0,1),
             model.getVariable(model.findVariable(R_VARIABLE_TEMPERATURE)).getValue(0,1));
    QCOMPARE(resumedModel.getVariable(resumedModel.findVariable(R_VARIABLE_DISPLACEMENT)).getValue(0,1),
             model.getVariable(model.findVariable(R_VARIABLE_DISPLACEMENT)).getValue(0,1));

    RModel changedMesh(model);
    changedMesh.getNode(1).set(1.000001,0.0,0.0);
    bool rejectedMesh = false;
    try { (void)RStructuralThermalFEAdapter::advance(changedMesh,scheduler,retry.candidate); }
    catch (const RError &) { rejectedMesh = true; }
    catch (const std::exception &) { rejectedMesh = true; }
    QVERIFY(rejectedMesh);
}

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

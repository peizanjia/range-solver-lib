#include <QtTest>

#include <QTemporaryDir>

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>
#include <vector>

#include <rbl_error.h>

#include <rml_boundary_condition.h>
#include <rml_element_group.h>
#include <rml_element.h>
#include <rml_file_utils.h>
#include <rml_line.h>
#include <rml_material.h>
#include <rml_model.h>
#include <rml_problem_task_item.h>
#include <rml_structural_cable_setup.h>
#include <rml_structural_control_setup.h>
#include <rml_structural_dynamics_restart_data.h>
#include <rml_structural_dynamics_setup.h>
#include <rml_structural_interaction.h>
#include <rml_time_solver.h>
#include <rml_variable.h>

#include "rsolver.h"
#include "rstructural_control_conditioner.h"
#include "rstructural_control_protocol.h"

namespace
{
constexpr double cableLength = 2.0;
constexpr double cableArea = 0.01;
constexpr double youngModulus = 20000.0;
constexpr double density = 10.0;
constexpr double prestress = 1.0;
constexpr double timeStep = 0.005;
constexpr uint totalSteps = 6;
constexpr uint firstSegmentSteps = 3;
constexpr char actuatorId[] = "tip-force";
constexpr char controllerId[] = "cable-force-controller";

void setComponent(RBoundaryCondition &condition,
                  RVariableType type,
                  bool enabled,
                  double value = 0.0)
{
    const uint position = condition.findComponentPosition(type);
    if (position == RConstants::eod)
        throw std::logic_error("Boundary condition is missing a cable component.");
    RConditionComponent &component = condition.getComponent(position);
    component.setEnabled(enabled);
    component.clear();
    if (enabled) component.add(0.0,value);
}

RBoundaryCondition displacementCondition(const std::array<bool,3> &enabled)
{
    RBoundaryCondition condition(R_BOUNDARY_CONDITION_DISPLACEMENT);
    const RVariableType types[] = {
        R_VARIABLE_DISPLACEMENT_X,R_VARIABLE_DISPLACEMENT_Y,R_VARIABLE_DISPLACEMENT_Z
    };
    for (uint axis=0;axis<3;axis++) setComponent(condition,types[axis],enabled[axis]);
    return condition;
}

void addPoint(RModel &model,uint node)
{
    RElement point(R_ELEMENT_POINT);
    point.setNodeId(0,node);
    model.addElement(point,true,node);
}

void setZeroGravity(RModel &model)
{
    const RVariableType types[] = {
        R_VARIABLE_G_ACCELERATION_X,R_VARIABLE_G_ACCELERATION_Y,R_VARIABLE_G_ACCELERATION_Z
    };
    for (uint groupId=0;groupId<model.getNElementGroups();groupId++)
    {
        RElementGroup *group = model.getElementGroupPtr(groupId);
        if (!group) continue;
        for (RVariableType type : types)
        {
            REnvironmentCondition &condition =
                group->getEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
            const uint position = condition.findComponentPosition(type);
            if (position == RConstants::eod)
                throw std::logic_error("Environment condition is missing a gravity component.");
            RConditionComponent &component = condition.getComponent(position);
            component.setEnabled(true);
            component.clear();
            component.add(0.0,0.0);
        }
    }
}

RStructuralSensor makeSensor(const QString &id,RStructuralSensorQuantity quantity)
{
    RStructuralSensor sensor;
    sensor.id = id;
    sensor.target.scope = RStructuralTargetScope::Node;
    sensor.target.nodeId = 1;
    sensor.quantity = quantity;
    sensor.aggregation = RStructuralSensorAggregation::Point;
    sensor.frame = RStructuralFrame::Global;
    sensor.direction = RR3Vector(1.0,0.0,0.0);
    sensor.unit = RStructuralInteractionSetup::getSensorUnit(quantity);
    return sensor;
}

RModel makeCable(uint steps,bool enableControl = true)
{
    RModel model;
    model.setNNodes(2);
    model.getNode(0).set(0.0,0.0,0.0);
    model.getNode(1).set(cableLength,0.0,0.0);

    RElement cable(R_ELEMENT_TRUSS1);
    cable.setNodeId(0,0);
    cable.setNodeId(1,1);
    model.addElement(cable,true,0);

    RMaterial material(RMaterial::Solid);
    RMaterialProperty modulus(RMaterialProperty::ModulusOfElasiticity);
    modulus.add(293.15,youngModulus);
    material.add(modulus);
    RMaterialProperty poisson(RMaterialProperty::PoissonRatio);
    poisson.add(293.15,0.0);
    material.add(poisson);
    RMaterialProperty massDensity(RMaterialProperty::Density);
    massDensity.add(293.15,density);
    material.add(massDensity);
    RLine &line = model.getLine(0);
    line.setCrossArea(cableArea);
    line.setMaterial(material);

    RBoundaryCondition transverse(R_BOUNDARY_CONDITION_DISPLACEMENT);
    setComponent(transverse,R_VARIABLE_DISPLACEMENT_X,false);
    setComponent(transverse,R_VARIABLE_DISPLACEMENT_Y,true,0.0);
    setComponent(transverse,R_VARIABLE_DISPLACEMENT_Z,true,0.0);
    line.addBoundaryCondition(transverse);

    addPoint(model,0);
    addPoint(model,1);
    model.getPoint(0).addBoundaryCondition(displacementCondition({{true,true,true}}));
    model.getPoint(1).addBoundaryCondition(displacementCondition({{false,true,true}}));
    RBoundaryCondition deadForce(R_BOUNDARY_CONDITION_FORCE);
    setComponent(deadForce,R_VARIABLE_FORCE_X,true,prestress);
    setComponent(deadForce,R_VARIABLE_FORCE_Y,true,0.0);
    setComponent(deadForce,R_VARIABLE_FORCE_Z,true,0.0);
    model.getPoint(1).addBoundaryCondition(deadForce);
    setZeroGravity(model);

    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    RStructuralCableSetup &cableSetup = model.getProblemSetup().getStructuralCableSetup();
    cableSetup.enabled = true;
    cableSetup.absoluteForceTolerance = 1.0e-12;
    cableSetup.relativeForceTolerance = 1.0e-12;
    RStructuralCableMember member;
    member.elementId = 0;
    member.referenceLength = cableLength;
    member.initialTension = prestress;
    cableSetup.members.push_back(member);
    cableSetup.bindToModel(model);

    RStructuralDynamicsSetup &dynamics =
        model.getProblemSetup().getStructuralDynamicsSetup();
    dynamics.setEnabled(true);
    dynamics.setMethod(RStructuralTimeIntegrationMethod::Newmark);
    dynamics.setNewmarkGamma(0.5);
    dynamics.setNewmarkBeta(0.25);
    dynamics.setTimeStep(timeStep);
    dynamics.setNSteps(steps);
    dynamics.setOutputFrequency(1);

    if (enableControl)
    {
        RStructuralControlSetup &control = model.getProblemSetup().getStructuralControlSetup();
        control.enabled = true;
        control.controllerId = QString::fromLatin1(controllerId);
        control.deterministicSeed = 918273;
        control.channelConditioning.push_back({QString::fromLatin1(actuatorId),1,2,0.05});

        RStructuralActuator actuator;
        actuator.id = QString::fromLatin1(actuatorId);
        actuator.type = RStructuralActuatorType::NodalForce;
        actuator.target.scope = RStructuralTargetScope::Node;
        actuator.target.nodeId = 1;
        actuator.variable = R_VARIABLE_FORCE_X;
        actuator.frame = RStructuralFrame::Global;
        actuator.direction = RR3Vector(1.0,0.0,0.0);
        actuator.unit = RVariable::getUnits(R_VARIABLE_FORCE_X);
        actuator.minimum = -1.0;
        actuator.maximum = 1.0;
        actuator.rateLimit = 20.0;
        auto &interaction = model.getProblemSetup().getStructuralInteractionSetup();
        interaction.getActuators().push_back(actuator);
        interaction.getSensors().push_back(makeSensor(
            QStringLiteral("tip-displacement"),RStructuralSensorQuantity::Displacement));
        interaction.getSensors().push_back(makeSensor(
            QStringLiteral("tip-velocity"),RStructuralSensorQuantity::Velocity));
        interaction.getSensors().push_back(makeSensor(
            QStringLiteral("tip-acceleration"),RStructuralSensorQuantity::Acceleration));
    }
    return model;
}

struct Captures
{
    std::vector<RStructuralControlObservation> observations;
    std::vector<RStructuralControlStepRecord> acceptedSteps;
};

RStructuralControlDecision decide(const RStructuralControlObservation &observation)
{
    double displacement = 0.0;
    double velocity = 0.0;
    double acceleration = 0.0;
    for (const auto &sensor : observation.sensors)
    {
        if (sensor.sensorId == QStringLiteral("tip-displacement")) displacement = sensor.value;
        if (sensor.sensorId == QStringLiteral("tip-velocity")) velocity = sensor.value;
        if (sensor.sensorId == QStringLiteral("tip-acceleration")) acceleration = sensor.value;
    }
    RStructuralControlDecision decision;
    decision.requestedCommands.push_back({QString::fromLatin1(actuatorId),
        0.65-10.0*displacement-0.05*velocity-0.001*acceleration});
    decision.nextControllerState = QByteArray::number(observation.nextTime,'g',17);
    return decision;
}

bool solve(RModel &model,Captures &captures,const QString &outputPath,QString &error)
{
    error.clear();
    try
    {
        RSolver solver(model,outputPath,QString());
        solver.setStructuralControlCallback(
            QString::fromLatin1(controllerId),QStringLiteral("test-revision-1"),
            [&captures](const RStructuralControlObservation &observation)
            {
                captures.observations.push_back(observation);
                return decide(observation);
            });
        solver.setStructuralControlStepObserver(
            [&captures](const RStructuralControlStepRecord &record)
            {
                captures.acceptedSteps.push_back(record);
            });
        solver.run();
        return true;
    }
    catch (const RError &failure)
    {
        error = failure.getMessage();
    }
    catch (const std::exception &failure)
    {
        error = QString::fromUtf8(failure.what());
    }
    return false;
}

bool load(RModel &model,const QString &path,QString &error)
{
    error.clear();
    try
    {
        model.read(path);
        return true;
    }
    catch (const RError &failure)
    {
        error = failure.getMessage();
    }
    catch (const std::exception &failure)
    {
        error = QString::fromUtf8(failure.what());
    }
    return false;
}

bool close(double actual,double expected,double tolerance = 2.0e-10)
{
    return std::isfinite(actual) && std::isfinite(expected) &&
           std::fabs(actual-expected) <= tolerance*std::max({1.0,std::fabs(actual),std::fabs(expected)});
}

double resultValue(const RModel &model,RVariableType type,uint component,uint node)
{
    const uint position = model.findVariable(type);
    if (position == RConstants::eod)
        return std::numeric_limits<double>::quiet_NaN();
    return model.getVariable(position).getValue(component,node);
}

void compareCheckpoints(const RStructuralDynamicsRestartData &left,
                        const RStructuralDynamicsRestartData &right)
{
    QCOMPARE(left.available,right.available);
    QCOMPARE(left.completedSteps,right.completedSteps);
    QVERIFY(close(left.time,right.time));
    QCOMPARE(left.systemFingerprint,right.systemFingerprint);
    QCOMPARE(left.controlAvailable,right.controlAvailable);
    QCOMPARE(left.controlProtocolVersion,right.controlProtocolVersion);
    QCOMPARE(left.controlFingerprint,right.controlFingerprint);
    QCOMPARE(left.controlDeterministicSeed,right.controlDeterministicSeed);
    QCOMPARE(left.controlChannelIds.size(),right.controlChannelIds.size());
    for (size_t i=0;i<left.controlChannelIds.size();i++)
    {
        QCOMPARE(left.controlChannelIds[i],right.controlChannelIds[i]);
        QVERIFY(close(left.controlLastCommands[i],right.controlLastCommands[i]));
        QVERIFY(close(left.controlPowers[i],right.controlPowers[i]));
        QVERIFY(close(left.cumulativeControlWork[i],right.cumulativeControlWork[i]));
    }
    QCOMPARE(left.controlControllerState,right.controlControllerState);
    QCOMPARE(left.q.size(),right.q.size());
    QCOMPARE(left.v.size(),right.v.size());
    QCOMPARE(left.a.size(),right.a.size());
    QCOMPARE(left.externalForce.size(),right.externalForce.size());
    for (size_t i=0;i<left.q.size();i++)
    {
        QVERIFY(close(left.q[i],right.q[i]));
        QVERIFY(close(left.v[i],right.v[i]));
        QVERIFY(close(left.a[i],right.a[i]));
        QVERIFY(close(left.externalForce[i],right.externalForce[i]));
    }
    QVERIFY(close(left.cumulativeExternalWork,right.cumulativeExternalWork));
    QVERIFY(close(left.cumulativeDampingDissipation,right.cumulativeDampingDissipation));
}

void compareRecords(const RStructuralControlStepRecord &left,
                    const RStructuralControlStepRecord &right)
{
    QVERIFY(close(left.acceptedTime,right.acceptedTime));
    QVERIFY(close(left.nextTime,right.nextTime));
    QCOMPARE(left.sensors.size(),right.sensors.size());
    for (size_t i=0;i<left.sensors.size();i++)
    {
        QCOMPARE(left.sensors[i].sensorId,right.sensors[i].sensorId);
        QCOMPARE(left.sensors[i].unit,right.sensors[i].unit);
        QVERIFY(close(left.sensors[i].value,right.sensors[i].value));
    }
    QCOMPARE(left.requestedCommands.size(),right.requestedCommands.size());
    QCOMPARE(left.conditionedCommands.size(),right.conditionedCommands.size());
    for (size_t i=0;i<left.requestedCommands.size();i++)
    {
        QCOMPARE(left.requestedCommands[i].actuatorId,right.requestedCommands[i].actuatorId);
        QVERIFY(close(left.requestedCommands[i].value,right.requestedCommands[i].value));
        QCOMPARE(left.conditionedCommands[i].actuatorId,right.conditionedCommands[i].actuatorId);
        QVERIFY(close(left.conditionedCommands[i].value,right.conditionedCommands[i].value));
    }
    QCOMPARE(left.appliedChannels.size(),right.appliedChannels.size());
    for (size_t i=0;i<left.appliedChannels.size();i++)
    {
        QCOMPARE(left.appliedChannels[i].actuatorId,right.appliedChannels[i].actuatorId);
        QVERIFY(close(left.appliedChannels[i].command,right.appliedChannels[i].command));
        QVERIFY(close(left.appliedChannels[i].endpointPower,right.appliedChannels[i].endpointPower));
        QVERIFY(close(left.appliedChannels[i].cumulativeWork,right.appliedChannels[i].cumulativeWork));
    }
    QCOMPARE(left.controllerStateSha256,right.controllerStateSha256);
}
}

class TestStructuralCableControl : public QObject
{
    Q_OBJECT

private slots:
    void acceptedControlMatchesNewmarkWorkAndPower();
    void velocityFeedbackReducesActualCableVibration();
    void persistedSplitControlMatchesFullRun_data();
    void persistedSplitControlMatchesFullRun();
    void failedStepDoesNotCommitControllerCandidate();
    void controlDefinitionChangeRejectsRestart();
    void callbackCannotChangeAcceptedClock();
    void callbackCannotChangeConditioningSetup();
    void callbackCannotChangeUnrelatedResults();
    void unsupportedActuatorFailsClosed();
};

void TestStructuralCableControl::acceptedControlMatchesNewmarkWorkAndPower()
{
    RModel model = makeCable(1);
    auto &conditioning = model.getProblemSetup().getStructuralControlSetup().channelConditioning;
    conditioning.clear();
    auto &actuator = model.getProblemSetup().getStructuralInteractionSetup().getActuators().front();
    actuator.rateLimit = 1000.0;

    Captures captures;
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    QString error;
    QVERIFY2(solve(model,captures,directory.filePath(QStringLiteral("step.ram")),error),
             qPrintable(error));
    QCOMPARE(captures.observations.size(),std::size_t(1));
    QCOMPARE(captures.acceptedSteps.size(),std::size_t(1));
    const auto &record = captures.acceptedSteps.front();
    QCOMPARE(record.requestedCommands.size(),std::size_t(1));
    QCOMPARE(record.conditionedCommands.size(),std::size_t(1));
    QCOMPARE(record.appliedChannels.size(),std::size_t(1));
    const double force = record.appliedChannels[0].command;
    QVERIFY(close(force,record.requestedCommands[0].value));
    QVERIFY(force > 0.0);

    const auto &checkpoint = model.getProblemSetup().getStructuralDynamicsRestartData();
    QVERIFY(checkpoint.controlAvailable);
    QCOMPARE(checkpoint.controlChannelIds.size(),std::size_t(1));
    QVERIFY(close(checkpoint.controlLastCommands[0],force));

    const double mass = density*cableArea*cableLength/2.0;
    const double stiffness = youngModulus*cableArea/cableLength;
    const double effective = 4.0*mass/(timeStep*timeStep);
    const double expectedDisplacement = 2.0*force/(effective+stiffness);
    const double displacement = checkpoint.q[3];
    const double velocity = checkpoint.v[3];
    const double acceleration = checkpoint.a[3];
    QVERIFY2(close(displacement,expectedDisplacement,2.0e-7),
             qPrintable(QString("controlled Newmark displacement %1, expected %2")
                            .arg(displacement,0,'g',17).arg(expectedDisplacement,0,'g',17)));
    const double expectedEndpointPower = force*velocity;
    const double expectedControllerWork =
        0.5*force*(0.0+velocity)*timeStep;
    QVERIFY(close(checkpoint.controlPowers[0],expectedEndpointPower,2.0e-7));
    QVERIFY(close(checkpoint.cumulativeControlWork[0],expectedControllerWork,2.0e-7));
    QVERIFY(close(record.appliedChannels[0].endpointPower,expectedEndpointPower,2.0e-7));
    QVERIFY(close(record.appliedChannels[0].cumulativeWork,expectedControllerWork,2.0e-7));
    QVERIFY(std::isfinite(acceleration));
    QVERIFY(!checkpoint.controlControllerState.isEmpty());
}

void TestStructuralCableControl::velocityFeedbackReducesActualCableVibration()
{
    constexpr uint steps=200;
    constexpr double initialVelocity=0.01;
    constexpr double gain=1.0;
    const double mass=density*cableArea*cableLength/2.0;
    const double stiffness=youngModulus*cableArea/cableLength;
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    double passiveEnergy=0.0;
    for (bool controlled : {false,true})
    {
        RModel model=makeCable(steps,controlled);
        auto &initial=model.getVariable(model.addVariable(R_VARIABLE_VELOCITY));
        initial.setApplyType(R_VARIABLE_APPLY_NODE);
        initial.resize(3,model.getNNodes());
        initial.setValue(0,1,initialVelocity);
        if (controlled)
        {
            model.getProblemSetup().getStructuralControlSetup().channelConditioning.clear();
            model.getProblemSetup().getStructuralInteractionSetup().getActuators().front().rateLimit=1000.0;
        }
        RSolver solver(model,directory.filePath(controlled ? "feedback.rbm" : "passive.rbm"),{});
        if (controlled) solver.setStructuralControlCallback(controllerId,"velocity-feedback-v1",
            [](const RStructuralControlObservation &observation)
            {
                double velocity=0.0;
                for (const auto &sensor : observation.sensors)
                    if (sensor.sensorId=="tip-velocity") velocity=sensor.value;
                RStructuralControlDecision decision;
                decision.requestedCommands.push_back({actuatorId,-gain*velocity});
                return decision;
            });
        QString failure;
        try { solver.run(); }
        catch (const RError &error) { failure=error.getMessage(); }
        catch (const std::exception &error) { failure=QString::fromUtf8(error.what()); }
        QVERIFY2(failure.isEmpty(),qPrintable(failure));
        const auto &checkpoint=model.getProblemSetup().getStructuralDynamicsRestartData();
        double q=0.0,v=initialVelocity,work=0.0;
        // Independent scalar recurrence with the force held over each step
        // and right-limit acceleration refreshed at the accepted boundary.
        for (uint i=0;i<steps;i++)
        {
            const double force=controlled ? -gain*v : 0.0;
            const double effective=4.0*mass/(timeStep*timeStep);
            const double next=(2.0*force+(effective-stiffness)*q+4.0*mass/timeStep*v)/(effective+stiffness);
            const double nextVelocity=2.0*(next-q)/timeStep-v;
            work+=force*(next-q);
            q=next;
            v=nextVelocity;
        }
        QVERIFY2(close(checkpoint.q[3],q,2.0e-8),qPrintable(QString("controlled=%1 FE=%2 recurrence=%3")
            .arg(controlled).arg(checkpoint.q[3],0,'g',17).arg(q,0,'g',17)));
        QVERIFY(close(checkpoint.v[3],v,2.0e-8));
        const double energy=0.5*stiffness*checkpoint.q[3]*checkpoint.q[3]+0.5*mass*checkpoint.v[3]*checkpoint.v[3];
        if (!controlled)
        {
            passiveEnergy=energy;
            QVERIFY(close(energy,0.5*mass*initialVelocity*initialVelocity,1.0e-10));
        }
        else
        {
            QVERIFY(energy<0.001*passiveEnergy);
            QVERIFY(close(checkpoint.cumulativeControlWork[0],work,1.0e-10));
            QVERIFY(close(energy-0.5*mass*initialVelocity*initialVelocity,work,1.0e-10));
        }
    }
}

void TestStructuralCableControl::persistedSplitControlMatchesFullRun_data()
{
    QTest::addColumn<bool>("binary");
    QTest::newRow("ascii-ram") << false;
    QTest::newRow("binary-rbm") << true;
}

void TestStructuralCableControl::persistedSplitControlMatchesFullRun()
{
    QFETCH(bool,binary);
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString extension = binary ? QStringLiteral("rbm") : QStringLiteral("ram");
    const QString fullPath = directory.filePath(QStringLiteral("full.")+extension);
    const QString splitPath = directory.filePath(QStringLiteral("split.")+extension);

    Captures fullCaptures;
    Captures splitCaptures;
    RModel uninterrupted = makeCable(totalSteps);
    RModel first = makeCable(firstSegmentSteps);
    QString error;
    QVERIFY2(solve(uninterrupted,fullCaptures,fullPath,error),qPrintable(error));
    QVERIFY2(solve(first,splitCaptures,splitPath,error),qPrintable(error));
    QVERIFY(first.getProblemSetup().getStructuralDynamicsRestartData().controlAvailable);

    RModel resumed;
    QVERIFY2(load(resumed,splitPath,error),qPrintable(error));
    resumed.getProblemSetup().setRestart(true);
    resumed.getProblemSetup().getStructuralDynamicsSetup().setNSteps(
        totalSteps-firstSegmentSteps);
    QVERIFY2(solve(resumed,splitCaptures,splitPath,error),qPrintable(error));

    QCOMPARE(fullCaptures.observations.size(),std::size_t(totalSteps));
    QCOMPARE(splitCaptures.observations.size(),std::size_t(totalSteps));
    QCOMPARE(fullCaptures.acceptedSteps.size(),std::size_t(totalSteps));
    QCOMPARE(splitCaptures.acceptedSteps.size(),std::size_t(totalSteps));
    for (uint step=0;step<totalSteps;step++)
    {
        const auto &fullObservation = fullCaptures.observations[step];
        const auto &splitObservation = splitCaptures.observations[step];
        QVERIFY(close(fullObservation.acceptedTime,double(step)*timeStep));
        QVERIFY(close(fullObservation.nextTime,double(step+1)*timeStep));
        QCOMPARE(fullObservation.sensors.size(),std::size_t(3));
        QCOMPARE(splitObservation.sensors.size(),std::size_t(3));
        for (size_t sensor=0;sensor<3;sensor++)
        {
            QCOMPARE(fullObservation.sensors[sensor].sensorId,
                     splitObservation.sensors[sensor].sensorId);
            QCOMPARE(fullObservation.sensors[sensor].unit,
                     splitObservation.sensors[sensor].unit);
            QVERIFY(close(fullObservation.sensors[sensor].value,
                          splitObservation.sensors[sensor].value));
        }
        QCOMPARE(fullObservation.priorControllerState,splitObservation.priorControllerState);
        compareRecords(fullCaptures.acceptedSteps[step],splitCaptures.acceptedSteps[step]);
        if (step == 0)
        {
            const auto &record = fullCaptures.acceptedSteps[step];
            QVERIFY(close(record.requestedCommands[0].value,0.65));
            QVERIFY(close(record.conditionedCommands[0].value,0.0));
            QVERIFY(close(record.appliedChannels[0].command,0.0));
            QVERIFY(close(record.appliedChannels[0].endpointPower,0.0));
        }

        RModel fullState;
        RModel splitState;
        QVERIFY2(load(fullState,RFileUtils::getFileNameWithTimeStep(fullPath,step+1),error),
                 qPrintable(error));
        QVERIFY2(load(splitState,RFileUtils::getFileNameWithTimeStep(splitPath,step+1),error),
                 qPrintable(error));
        const auto &fullCheckpoint =
            fullState.getProblemSetup().getStructuralDynamicsRestartData();
        const auto &splitCheckpoint =
            splitState.getProblemSetup().getStructuralDynamicsRestartData();
        compareCheckpoints(fullCheckpoint,splitCheckpoint);
        const double q = fullCheckpoint.q[3];
        const double v = fullCheckpoint.v[3];
        const double tension = prestress+youngModulus*cableArea/cableLength*q;
        const double fullTension = resultValue(fullState,R_VARIABLE_FORCE,0,1);
        const double splitTension = resultValue(splitState,R_VARIABLE_FORCE,0,1);
        QVERIFY(close(fullTension,tension,binary ? 2.0e-10 : 1.0e-5));
        QVERIFY(close(splitTension,tension,binary ? 2.0e-10 : 1.0e-5));
        QVERIFY(close(fullTension,splitTension,binary ? 2.0e-10 : 1.0e-5));
        const double mechanicalEnergy = tension*tension/
            (2.0*youngModulus*cableArea/cableLength)+
            0.5*(density*cableArea*cableLength/2.0)*v*v;
        QVERIFY(close(mechanicalEnergy,
                      fullCheckpoint.initialMechanicalEnergy+
                          fullCheckpoint.cumulativeExternalWork,
                      2.0e-7));
    }

    const auto &fullFinal = uninterrupted.getProblemSetup().getStructuralDynamicsRestartData();
    const auto &splitFinal = resumed.getProblemSetup().getStructuralDynamicsRestartData();
    compareCheckpoints(fullFinal,splitFinal);
    QCOMPARE(splitFinal.completedSteps,totalSteps);
}

void TestStructuralCableControl::failedStepDoesNotCommitControllerCandidate()
{
    RModel model = makeCable(1);
    auto &conditioning = model.getProblemSetup().getStructuralControlSetup().channelConditioning;
    conditioning.clear();
    conditioning.push_back({QString::fromLatin1(actuatorId),0,2,0.0});
    auto &actuator = model.getProblemSetup().getStructuralInteractionSetup().getActuators().front();
    actuator.maximum = std::numeric_limits<double>::infinity();
    actuator.rateLimit = std::numeric_limits<double>::infinity();
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    RSolver solver(model,directory.filePath(QStringLiteral("retry.ram")),QString());
    std::vector<QByteArray> priorStates;
    std::vector<RStructuralControlStepRecord> accepted;
    uint callbackCount = 0;
    solver.setStructuralControlCallback(
        QString::fromLatin1(controllerId),QStringLiteral("retry-controller-v1"),
        [&priorStates,&callbackCount](const RStructuralControlObservation &observation)
        {
            priorStates.push_back(observation.priorControllerState);
            RStructuralControlDecision decision;
            decision.nextControllerState = callbackCount == 0
                ? QByteArrayLiteral("rejected-candidate") : QByteArrayLiteral("accepted-retry");
            decision.requestedCommands.push_back({QString::fromLatin1(actuatorId),
                callbackCount++ == 0 ? std::numeric_limits<double>::max() : 0.2});
            return decision;
        });
    solver.setStructuralControlStepObserver(
        [&accepted](const RStructuralControlStepRecord &record)
        {
            accepted.push_back(record);
        });

    QString rejectionError;
    try { solver.run(); }
    catch (const RError &failure) { rejectionError = failure.getMessage(); }
    catch (const std::exception &failure) { rejectionError = QString::fromUtf8(failure.what()); }
    QVERIFY2(!rejectionError.isEmpty(),"Expected the first huge-force candidate to be rejected.");
    QVERIFY2(rejectionError.contains(QStringLiteral("right-limit acceleration"),Qt::CaseInsensitive),
             qPrintable(rejectionError));
    QVERIFY(!model.getProblemSetup().getStructuralDynamicsRestartData().available);
    QVERIFY(!model.getProblemSetup().getStructuralDynamicsRestartData().controlAvailable);

    // The callback sees the same accepted controller state on retry; its
    // rejected candidate is not committed to either the model or the ledger.
    try { solver.run(); }
    catch (const RError &failure)
    {
        QFAIL(qPrintable(QString("Retry failed: %1").arg(failure.getMessage())));
    }
    QCOMPARE(priorStates.size(),std::size_t(2));
    QVERIFY(priorStates[0].isEmpty());
    QVERIFY(priorStates[1].isEmpty());
    QCOMPARE(accepted.size(),std::size_t(1));
    QCOMPARE(accepted[0].requestedCommands[0].value,0.2);
    const auto &checkpoint = model.getProblemSetup().getStructuralDynamicsRestartData();
    QVERIFY(checkpoint.available && checkpoint.controlAvailable);
    const auto envelope = RStructuralControlStateEnvelope::deserialize(
        checkpoint.controlControllerState);
    QCOMPARE(envelope.controllerState,QByteArrayLiteral("accepted-retry"));
    QCOMPARE(checkpoint.completedSteps,1u);
}

void TestStructuralCableControl::controlDefinitionChangeRejectsRestart()
{
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString path = directory.filePath(QStringLiteral("changed.ram"));
    Captures captures;
    RModel first = makeCable(firstSegmentSteps);
    QString error;
    QVERIFY2(solve(first,captures,path,error),qPrintable(error));

    RModel resumed;
    QVERIFY2(load(resumed,path,error),qPrintable(error));
    resumed.getProblemSetup().setRestart(true);
    resumed.getProblemSetup().getStructuralDynamicsSetup().setNSteps(
        totalSteps-firstSegmentSteps);
    resumed.getProblemSetup().getStructuralInteractionSetup().getActuators().front().maximum = 0.9;
    Captures resumedCaptures;
    QVERIFY(!solve(resumed,resumedCaptures,path,error));
    QVERIFY2(error.contains(QStringLiteral("control restart"),Qt::CaseInsensitive),
             qPrintable(error));
    QCOMPARE(resumed.getProblemSetup().getStructuralDynamicsRestartData().completedSteps,
             firstSegmentSteps);
}

void TestStructuralCableControl::callbackCannotChangeAcceptedClock()
{
    RModel model = makeCable(2);
    RSolver solver(model,QString(),QString());
    bool callbackCalled = false;
    solver.setStructuralControlCallback(
        QString::fromLatin1(controllerId),QStringLiteral("tamper-clock-v1"),
        [&model,&callbackCalled](const RStructuralControlObservation &)
        {
            callbackCalled = true;
            model.getTimeSolver().setCurrentTimeStep(1);
            RStructuralControlDecision decision;
            decision.requestedCommands.push_back({QString::fromLatin1(actuatorId),0.2});
            return decision;
        });
    QString error;
    try { solver.run(); }
    catch (const RError &failure) { error = failure.getMessage(); }
    QVERIFY(callbackCalled);
    QVERIFY2(error.contains(QStringLiteral("frozen time solver"),Qt::CaseInsensitive),
             qPrintable(error));
    QVERIFY(!model.getProblemSetup().getStructuralDynamicsRestartData().available);
    QCOMPARE(model.getTimeSolver().getCurrentTimeStep(),0u);
}

void TestStructuralCableControl::callbackCannotChangeConditioningSetup()
{
    RModel model = makeCable(1);
    const uint originalDelay = model.getProblemSetup().getStructuralControlSetup()
        .channelConditioning.front().delaySteps;
    RSolver solver(model,QString(),QString());
    bool callbackCalled = false;
    solver.setStructuralControlCallback(
        QString::fromLatin1(controllerId),QStringLiteral("tamper-conditioning-v1"),
        [&model,&callbackCalled](const RStructuralControlObservation &)
        {
            callbackCalled = true;
            ++model.getProblemSetup().getStructuralControlSetup()
                .channelConditioning.front().delaySteps;
            RStructuralControlDecision decision;
            decision.requestedCommands.push_back({QString::fromLatin1(actuatorId),0.2});
            return decision;
        });
    QString error;
    try { solver.run(); }
    catch (const RError &failure) { error = failure.getMessage(); }
    QVERIFY(callbackCalled);
    QVERIFY2(error.contains(QStringLiteral("control setup"),Qt::CaseInsensitive),
             qPrintable(error));
    QVERIFY(!model.getProblemSetup().getStructuralDynamicsRestartData().available);
    QCOMPARE(model.getProblemSetup().getStructuralControlSetup()
                 .channelConditioning.front().delaySteps,originalDelay);
}

void TestStructuralCableControl::callbackCannotChangeUnrelatedResults()
{
    for (uint change=0;change<3;change++)
    {
        RModel model = makeCable(1);
        RSolver solver(model,QString(),QString());
        bool first = true;
        std::vector<QByteArray> priorStates;
        solver.setStructuralControlCallback(
            QString::fromLatin1(controllerId),QStringLiteral("tamper-results-v1"),
            [&model,&first,&priorStates,change](const RStructuralControlObservation &observation)
            {
                priorStates.push_back(observation.priorControllerState);
                if (first)
                {
                    first = false;
                    RVariable field(change==2 ? R_VARIABLE_PATCH_SIZE : R_VARIABLE_PRESSURE);
                    field.setApplyType(R_VARIABLE_APPLY_NODE);
                    field.resize(1,model.getNNodes());
                    field.setValue(0,1,321.0);
                    if (change==1) field.setUnits("invalid-unit");
                    model.addVariable(field);
                }
                RStructuralControlDecision decision;
                decision.requestedCommands.push_back({QString::fromLatin1(actuatorId),0.2});
                decision.nextControllerState = QByteArrayLiteral("accepted-results-guard");
                return decision;
            });
        QString error;
        try { solver.run(); }
        catch (const RError &failure) { error = failure.getMessage(); }
        QVERIFY2(error.contains("frozen result fields",Qt::CaseInsensitive),qPrintable(error));
        QVERIFY(!model.getProblemSetup().getStructuralDynamicsRestartData().available);
        QCOMPARE(model.findVariable(R_VARIABLE_PRESSURE),RConstants::eod);
        QCOMPARE(model.findVariable(R_VARIABLE_PATCH_SIZE),RConstants::eod);
        try { solver.run(); }
        catch (const RError &failure) { QFAIL(qPrintable(failure.getMessage())); }
        QCOMPARE(priorStates.size(),std::size_t(2));
        QCOMPARE(priorStates[0],priorStates[1]);
        QVERIFY(model.getProblemSetup().getStructuralDynamicsRestartData().controlAvailable);
        QCOMPARE(model.findVariable(R_VARIABLE_PRESSURE),RConstants::eod);
    }
}

void TestStructuralCableControl::unsupportedActuatorFailsClosed()
{
    RModel model = makeCable(1);
    RStructuralActuator &actuator =
        model.getProblemSetup().getStructuralInteractionSetup().getActuators().front();
    actuator.type = RStructuralActuatorType::BaseMotion;
    actuator.variable = R_VARIABLE_VELOCITY;
    actuator.unit = RVariable::getUnits(R_VARIABLE_VELOCITY);

    Captures captures;
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    QString error;
    QVERIFY(!solve(model,captures,directory.filePath(QStringLiteral("unsupported.ram")),error));
    QVERIFY2(error.contains(QStringLiteral("NodalForce"),Qt::CaseInsensitive),
             qPrintable(error));
    QVERIFY(!model.getProblemSetup().getStructuralDynamicsRestartData().available);
}

QTEST_APPLESS_MAIN(TestStructuralCableControl)
#include "tst_structural_cable_control.moc"

#include <QtTest>

#include <QByteArray>

#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <rbl_error.h>

#include <rml_boundary_condition.h>
#include <rml_model.h>
#include <rml_problem_task_item.h>
#include <rml_structural_cable_setup.h>
#include <rml_structural_dynamics_setup.h>
#include <rml_variable.h>

#include "rsolver.h"
#include "rsolverstress.h"

namespace
{

void setComponent(RBoundaryCondition &condition,
                  RVariableType type,
                  bool enabled,
                  double value = 0.0)
{
    const uint position = condition.findComponentPosition(type);
    if (position == RConstants::eod)
    {
        throw std::logic_error("The boundary condition is missing a required component.");
    }
    RConditionComponent &component = condition.getComponent(position);
    component.setEnabled(enabled);
    component.clear();
    if (enabled)
    {
        component.add(0.0,value);
    }
}

RBoundaryCondition displacementCondition(const std::array<bool,3> &enabled,
                                         const std::array<double,3> &values)
{
    RBoundaryCondition condition(R_BOUNDARY_CONDITION_DISPLACEMENT);
    const RVariableType types[3] =
    {
        R_VARIABLE_DISPLACEMENT_X,
        R_VARIABLE_DISPLACEMENT_Y,
        R_VARIABLE_DISPLACEMENT_Z
    };
    for (uint axis=0;axis<3;axis++)
    {
        setComponent(condition,types[axis],enabled[axis],values[axis]);
    }
    return condition;
}

RBoundaryCondition forceCondition(const std::array<double,3> &values)
{
    RBoundaryCondition condition(R_BOUNDARY_CONDITION_FORCE);
    const RVariableType types[3] =
    {
        R_VARIABLE_FORCE_X,
        R_VARIABLE_FORCE_Y,
        R_VARIABLE_FORCE_Z
    };
    for (uint axis=0;axis<3;axis++)
    {
        setComponent(condition,types[axis],true,values[axis]);
    }
    return condition;
}

void addPointElement(RModel &model, uint node, uint pointGroup)
{
    RElement point(R_ELEMENT_POINT);
    point.setNodeId(0,node);
    model.addElement(point,true,pointGroup);
}

void assignCableMaterial(RLine &line, double area, double youngModulus,
                         double density = 0.0)
{
    line.setCrossArea(area);
    RMaterial material(RMaterial::Solid);
    RMaterialProperty modulus(RMaterialProperty::ModulusOfElasiticity);
    modulus.add(293.15,youngModulus);
    material.add(modulus);
    if (density > 0.0)
    {
        RMaterialProperty massDensity(RMaterialProperty::Density);
        massDensity.add(293.15,density);
        material.add(massDensity);
    }
    line.setMaterial(material);
}

void setEnvironmentValue(RElementGroup &group,
                         REnvironmentConditionType conditionType,
                         RVariableType variableType,
                         double value)
{
    REnvironmentCondition &condition = group.getEnvironmentCondition(conditionType);
    const uint position = condition.findComponentPosition(variableType);
    if (position == RConstants::eod)
    {
        throw std::logic_error("The environment condition is missing a required component.");
    }
    RConditionComponent &component = condition.getComponent(position);
    component.setEnabled(true);
    component.clear();
    component.add(0.0,value);
}

void setZeroGravity(RModel &model)
{
    const RVariableType gravityTypes[3] =
    {
        R_VARIABLE_G_ACCELERATION_X,
        R_VARIABLE_G_ACCELERATION_Y,
        R_VARIABLE_G_ACCELERATION_Z
    };
    for (uint groupId=0;groupId<model.getNElementGroups();groupId++)
    {
        RElementGroup *group = model.getElementGroupPtr(groupId);
        if (!group)
        {
            continue;
        }
        for (RVariableType type : gravityTypes)
        {
            setEnvironmentValue(*group,R_ENVIRONMENT_CONDITION_G_ACCELERATION,type,0.0);
        }
    }
}

void addCableMember(RStructuralCableSetup &setup,
                    uint elementId,
                    double referenceLength,
                    double initialTension)
{
    RStructuralCableMember member;
    member.elementId = elementId;
    member.referenceLength = referenceLength;
    member.initialTension = initialTension;
    setup.members.push_back(member);
}

RModel makeAxialCable(bool includeUnconfiguredSegment = false,
                      bool addContradictoryDisplacement = false,
                      double density = 0.0)
{
    RModel model;
    model.setNNodes(includeUnconfiguredSegment ? 3 : 2);
    model.getNode(0).set(0.0,0.0,0.0);
    model.getNode(1).set(2.0,0.0,0.0);
    if (includeUnconfiguredSegment)
    {
        model.getNode(2).set(4.0,0.0,0.0);
    }

    RElement first(R_ELEMENT_TRUSS1);
    first.setNodeId(0,0);
    first.setNodeId(1,1);
    model.addElement(first,true,0);
    if (includeUnconfiguredSegment)
    {
        RElement second(R_ELEMENT_TRUSS1);
        second.setNodeId(0,1);
        second.setNodeId(1,2);
        model.addElement(second,true,0);
    }
    assignCableMaterial(model.getLine(0),0.5,1000.0,density);

    addPointElement(model,0,0);
    addPointElement(model,1,1);
    model.getPoint(0).addBoundaryCondition(
        displacementCondition({{true,true,true}},{{0.0,0.0,0.0}}));
    model.getPoint(1).addBoundaryCondition(
        displacementCondition({{false,true,true}},{{0.0,0.0,0.0}}));
    model.getPoint(1).addBoundaryCondition(forceCondition({{10.0,0.0,0.0}}));

    if (addContradictoryDisplacement)
    {
        addPointElement(model,0,2);
        model.getPoint(2).addBoundaryCondition(
            displacementCondition({{true,false,false}},{{0.1,0.0,0.0}}));
    }

    setZeroGravity(model);

    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    RStructuralCableSetup &setup = model.getProblemSetup().getStructuralCableSetup();
    setup.enabled = true;
    addCableMember(setup,0,2.0,5.0);
    setup.bindToModel(model);
    return model;
}

RModel makeSymmetricSagCable()
{
    RModel model;
    model.setNNodes(3);
    model.getNode(0).set(-1.0,0.0,0.0);
    model.getNode(1).set(0.0,0.0,0.0);
    model.getNode(2).set(1.0,0.0,0.0);

    RElement left(R_ELEMENT_TRUSS1);
    left.setNodeId(0,0);
    left.setNodeId(1,1);
    model.addElement(left,true,0);
    RElement right(R_ELEMENT_TRUSS1);
    right.setNodeId(0,1);
    right.setNodeId(1,2);
    model.addElement(right,true,0);
    assignCableMaterial(model.getLine(0),0.1,10.0); // EA = 1 N.

    for (uint node=0;node<3;node++)
    {
        addPointElement(model,node,node);
    }
    model.getPoint(0).addBoundaryCondition(
        displacementCondition({{true,true,true}},{{0.0,0.0,0.0}}));
    model.getPoint(2).addBoundaryCondition(
        displacementCondition({{true,true,true}},{{0.0,0.0,0.0}}));
    model.getPoint(1).addBoundaryCondition(
        displacementCondition({{true,false,true}},{{0.0,0.0,0.0}}));
    model.getPoint(1).addBoundaryCondition(forceCondition({{0.0,-120.0,0.0}}));

    setZeroGravity(model);

    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    RStructuralCableSetup &setup = model.getProblemSetup().getStructuralCableSetup();
    setup.enabled = true;
    addCableMember(setup,0,1.25,100.0);
    addCableMember(setup,1,1.25,100.0);
    setup.bindToModel(model);
    return model;
}

RModel makeLoadedSlackCable()
{
    RModel model;
    model.setNNodes(2);
    model.getNode(0).set(0.0,0.0,0.0);
    model.getNode(1).set(2.0,0.0,0.0);
    RElement cable(R_ELEMENT_TRUSS1);
    cable.setNodeId(0,0);
    cable.setNodeId(1,1);
    model.addElement(cable,true,0);
    assignCableMaterial(model.getLine(0),0.5,1000.0);
    addPointElement(model,0,0);
    addPointElement(model,1,1);
    model.getPoint(0).addBoundaryCondition(
        displacementCondition({{true,true,true}},{{0.0,0.0,0.0}}));
    model.getPoint(1).addBoundaryCondition(
        displacementCondition({{true,false,true}},{{0.0,0.0,0.0}}));
    model.getPoint(1).addBoundaryCondition(forceCondition({{0.0,1.0,0.0}}));
    setZeroGravity(model);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    RStructuralCableSetup &setup = model.getProblemSetup().getStructuralCableSetup();
    setup.enabled = true;
    addCableMember(setup,0,2.0,0.0);
    setup.bindToModel(model);
    return model;
}

RModel makeFixedTwoMemberCableWithReversedSetupOrder()
{
    RModel model;
    model.setNNodes(4);
    model.getNode(0).set(0.0,0.0,0.0);
    model.getNode(1).set(1.0,0.0,0.0);
    model.getNode(2).set(3.0,0.0,0.0);
    model.getNode(3).set(5.0,0.0,0.0);

    RElement first(R_ELEMENT_TRUSS1);
    first.setNodeId(0,0);
    first.setNodeId(1,1);
    model.addElement(first,true,0);
    RElement second(R_ELEMENT_TRUSS1);
    second.setNodeId(0,2);
    second.setNodeId(1,3);
    model.addElement(second,true,0);
    assignCableMaterial(model.getLine(0),0.5,1000.0);

    for (uint node=0;node<4;node++)
    {
        addPointElement(model,node,node);
        model.getPoint(node).addBoundaryCondition(
            displacementCondition({{true,true,true}},{{0.0,0.0,0.0}}));
    }
    setZeroGravity(model);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    RStructuralCableSetup &setup = model.getProblemSetup().getStructuralCableSetup();
    setup.enabled = true;
    // The setup order is deliberately opposite to element order and each
    // member has a different initial tension to detect accidental index reuse.
    addCableMember(setup,1,2.0,13.0);
    addCableMember(setup,0,1.0,5.0);
    setup.bindToModel(model);
    return model;
}

double nodeValue(const RModel &model, RVariableType type, uint component, uint node)
{
    const uint position = model.findVariable(type);
    if (position == RConstants::eod)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return model.getVariable(position).getValue(component,node);
}

double elementValue(const RModel &model, RVariableType type, uint element)
{
    const uint position = model.findVariable(type);
    if (position == RConstants::eod)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return model.getVariable(position).getValue(0,element);
}

}

class TestStructuralCableFE : public QObject
{
    Q_OBJECT

private slots:
    void axialExtensionAndStressMatchAnalyticSolution();
    void symmetricPretensionedSagMatchesAnalyticEquilibrium();
    void staleMeshFingerprintRejectedBeforeCleanup();
    void unsupportedDynamicsOptionRejectedBeforeCleanup();
    void unconfiguredNonpointMemberRejected();
    void contradictoryPrescriptionsRejected();
    void loadedSlackCableFailsWithoutPublishingResults();
    void directModalCableSolverRejected();
    void cableModeCannotChangeAfterConstruction();
    void directSharedTemperatureMustBeReferenceTemperature();
    void nonReferenceEnvironmentTemperatureRejected();
    void nonzeroEnvironmentGravityRejected();
    void cableMembersFollowElementIdsWhenSetupOrderIsReversed();
    void transientSingleCableMatchesNewmarkAverageAcceleration();
    void failedSecondCableStepRollsBackToLastAcceptedState();
    void cableAndPiezoModesCannotBeEnabledTogether();
};

void TestStructuralCableFE::axialExtensionAndStressMatchAnalyticSolution()
{
    RModel model = makeAxialCable();
    try
    {
        RSolver solver(model,QString(),QString());
        solver.run();
    }
    catch (const RError &error)
    {
        const QByteArray message = error.getMessage().toUtf8();
        QFAIL(message.constData());
    }

    const double expectedExtension = (10.0-5.0)*2.0/(1000.0*0.5);
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_DISPLACEMENT,0,1)-expectedExtension) < 1.0e-9);
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_FORCE,0,0)+10.0) < 1.0e-8);
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_FORCE,0,1)-10.0) < 1.0e-8);
    QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_NORMAL,0)-20.0) < 1.0e-8);
    QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_VON_MISES,0)-20.0) < 1.0e-8);
    QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_X,0)-20.0) < 1.0e-8);
    QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_NORMAL,1)) < 1.0e-14);
    QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_VON_MISES,2)) < 1.0e-14);
}

void TestStructuralCableFE::symmetricPretensionedSagMatchesAnalyticEquilibrium()
{
    RModel model = makeSymmetricSagCable();
    try
    {
        RSolver solver(model,QString(),QString());
        solver.run();
    }
    catch (const RError &error)
    {
        const QByteArray message = error.getMessage().toUtf8();
        QFAIL(message.constData());
    }

    // The target chord length is the reference length, so each cable retains
    // exactly 100 N and the two vertical components balance the 120 N load.
    const double halfSpan = 1.0;
    const double sag = 0.75;
    const double chord = std::hypot(halfSpan,sag);
    const double expectedLoad = 2.0*100.0*sag/chord;
    QVERIFY(std::fabs(chord-1.25) < 1.0e-14);
    QVERIFY(std::fabs(expectedLoad-120.0) < 1.0e-14);
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_DISPLACEMENT,1,1)+sag) < 2.0e-7);
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_FORCE,0,1)) < 1.0e-7);
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_FORCE,1,1)+120.0) < 1.0e-7);
    QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_NORMAL,0)-1000.0) < 1.0e-6);
    QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_NORMAL,1)-1000.0) < 1.0e-6);
}

void TestStructuralCableFE::staleMeshFingerprintRejectedBeforeCleanup()
{
    RModel model = makeAxialCable();
    const uint displacementPosition = model.addVariable(R_VARIABLE_DISPLACEMENT);
    model.getVariable(displacementPosition).resize(3,model.getNNodes());
    model.getNode(1).setX(2.1);

    bool rejected = false;
    try
    {
        RSolver solver(model,QString(),QString());
    }
    catch (const RError &error)
    {
        rejected = error.getMessage().contains("different mesh",Qt::CaseInsensitive);
    }
    QVERIFY(rejected);
    QVERIFY(model.findVariable(R_VARIABLE_DISPLACEMENT) != RConstants::eod);
}

void TestStructuralCableFE::unsupportedDynamicsOptionRejectedBeforeCleanup()
{
    RModel model = makeAxialCable();
    const uint displacementPosition = model.addVariable(R_VARIABLE_DISPLACEMENT);
    model.getVariable(displacementPosition).resize(3,model.getNNodes());
    RStructuralDynamicsSetup &dynamics =
        model.getProblemSetup().getStructuralDynamicsSetup();
    dynamics.setEnabled(true);
    dynamics.setRayleighAlpha(0.1);

    bool rejected = false;
    try
    {
        RSolver solver(model,QString(),QString());
    }
    catch (const RError &error)
    {
        rejected = error.getMessage().contains("damping",Qt::CaseInsensitive);
    }
    QVERIFY(rejected);
    QVERIFY(model.findVariable(R_VARIABLE_DISPLACEMENT) != RConstants::eod);
}

void TestStructuralCableFE::unconfiguredNonpointMemberRejected()
{
    RModel model = makeAxialCable(true);
    bool rejected = false;
    try
    {
        RSolver solver(model,QString(),QString());
    }
    catch (const RError &error)
    {
        rejected = error.getMessage().contains("no explicit cable setup",Qt::CaseInsensitive);
    }
    QVERIFY(rejected);
}

void TestStructuralCableFE::contradictoryPrescriptionsRejected()
{
    RModel model = makeAxialCable(false,true);
    bool rejected = false;
    try
    {
        RSolver solver(model,QString(),QString());
    }
    catch (const RError &error)
    {
        rejected = error.getMessage().contains("contradictory",Qt::CaseInsensitive);
    }
    QVERIFY(rejected);
}

void TestStructuralCableFE::loadedSlackCableFailsWithoutPublishingResults()
{
    RModel model = makeLoadedSlackCable();
    RSolver solver(model,QString(),QString());
    bool failed = false;
    try
    {
        solver.run();
    }
    catch (const RError &)
    {
        failed = true;
    }
    QVERIFY(failed);
    QVERIFY(model.findVariable(R_VARIABLE_DISPLACEMENT) == RConstants::eod);
    QVERIFY(model.findVariable(R_VARIABLE_FORCE) == RConstants::eod);
    QVERIFY(model.findVariable(R_VARIABLE_STRESS_NORMAL) == RConstants::eod);
    QVERIFY(model.findVariable(R_VARIABLE_STRESS_VON_MISES) == RConstants::eod);
    QVERIFY(model.findVariable(R_VARIABLE_STRESS_X) == RConstants::eod);
}

void TestStructuralCableFE::directModalCableSolverRejected()
{
    RModel model = makeAxialCable();
    RSolverSharedData sharedData;
    bool rejected = false;
    try
    {
        RSolverStress solver(&model,QString(),QString(),sharedData,true);
    }
    catch (const RError &error)
    {
        rejected = error.getMessage().contains("non-modal",Qt::CaseInsensitive);
    }
    QVERIFY(rejected);
}

void TestStructuralCableFE::cableModeCannotChangeAfterConstruction()
{
    RModel model = makeAxialCable();
    RSolver solver(model,QString(),QString());
    model.getProblemSetup().getStructuralCableSetup().enabled = false;
    bool rejected = false;
    try
    {
        solver.run();
    }
    catch (const RError &error)
    {
        rejected = error.getMessage().contains("cannot be enabled, disabled or switched",Qt::CaseInsensitive);
    }
    QVERIFY(rejected);
    QVERIFY(model.findVariable(R_VARIABLE_DISPLACEMENT) == RConstants::eod);
}

void TestStructuralCableFE::directSharedTemperatureMustBeReferenceTemperature()
{
    RModel model = makeAxialCable();
    RSolverSharedData sharedData;
    sharedData.addData("element-temperature",RRVector(model.getNElements(),300.0));
    RSolverStress solver(&model,QString(),QString(),sharedData,false);

    bool rejected = false;
    try
    {
        solver.run(true,0);
    }
    catch (const RError &error)
    {
        rejected = error.getMessage().contains("isothermal at 293.15",Qt::CaseInsensitive);
    }
    QVERIFY(rejected);
    QVERIFY(model.findVariable(R_VARIABLE_DISPLACEMENT) == RConstants::eod);
}

void TestStructuralCableFE::nonReferenceEnvironmentTemperatureRejected()
{
    RModel model = makeAxialCable();
    setEnvironmentValue(model.getLine(0),R_ENVIRONMENT_CONDITION_TEMPERATURE,
                        R_VARIABLE_TEMPERATURE,300.0);

    bool rejected = false;
    try
    {
        RSolver solver(model,QString(),QString());
    }
    catch (const RError &error)
    {
        rejected = error.getMessage().contains("different temperature",Qt::CaseInsensitive);
    }
    QVERIFY(rejected);
    QVERIFY(model.findVariable(R_VARIABLE_DISPLACEMENT) == RConstants::eod);
}

void TestStructuralCableFE::nonzeroEnvironmentGravityRejected()
{
    RModel model = makeAxialCable();
    setEnvironmentValue(model.getLine(0),R_ENVIRONMENT_CONDITION_G_ACCELERATION,
                        R_VARIABLE_G_ACCELERATION_Z,-9.80665);

    bool rejected = false;
    try
    {
        RSolver solver(model,QString(),QString());
    }
    catch (const RError &error)
    {
        rejected = error.getMessage().contains("gravity components to be zero",Qt::CaseInsensitive);
    }
    QVERIFY(rejected);
    QVERIFY(model.findVariable(R_VARIABLE_DISPLACEMENT) == RConstants::eod);
}

void TestStructuralCableFE::cableMembersFollowElementIdsWhenSetupOrderIsReversed()
{
    RModel model = makeFixedTwoMemberCableWithReversedSetupOrder();
    try
    {
        RSolver solver(model,QString(),QString());
        solver.run();
    }
    catch (const RError &error)
    {
        const QByteArray message = error.getMessage().toUtf8();
        QFAIL(message.constData());
    }

    QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_NORMAL,0)-10.0) < 1.0e-9);
    QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_NORMAL,1)-26.0) < 1.0e-9);
    QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_X,0)-10.0) < 1.0e-9);
    QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_X,1)-26.0) < 1.0e-9);
}

void TestStructuralCableFE::transientSingleCableMatchesNewmarkAverageAcceleration()
{
    RModel model = makeAxialCable(false,false,2.0);
    model.getProblemSetup().getStructuralCableSetup().members[0].initialTension = 0.0;
    RStructuralDynamicsSetup &dynamics =
        model.getProblemSetup().getStructuralDynamicsSetup();
    dynamics.setEnabled(true);
    dynamics.setMethod(RStructuralTimeIntegrationMethod::Newmark);
    dynamics.setNewmarkGamma(0.5);
    dynamics.setNewmarkBeta(0.25);
    dynamics.setTimeStep(0.1);
    dynamics.setNSteps(1);
    dynamics.setOutputFrequency(1);

    try
    {
        RSolver solver(model,QString(),QString());
        solver.run();
    }
    catch (const RError &error)
    {
        const QByteArray message = error.getMessage().toUtf8();
        QFAIL(message.constData());
    }

    // rho*A*L = 2 kg, so the free endpoint receives one kilogram. For this
    // axial member k=EA/L0=250 N/m and the Newmark mass tangent is 400 N/m.
    const double expectedDisplacement = 10.0/(250.0+400.0);
    const double expectedAcceleration = 400.0*expectedDisplacement;
    const double expectedVelocity = 0.5*0.1*expectedAcceleration;
    const double expectedTension = 250.0*expectedDisplacement;
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_DISPLACEMENT,0,1)-expectedDisplacement) < 1.0e-11);
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_VELOCITY,0,1)-expectedVelocity) < 1.0e-10);
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_ACCELERATION,0,1)-expectedAcceleration) < 1.0e-9);
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_FORCE,0,1)-expectedTension) < 1.0e-8);
    QVERIFY(std::fabs(elementValue(model,R_VARIABLE_STRESS_NORMAL,0)-expectedTension/0.5) < 1.0e-7);
    QCOMPARE(model.getTimeSolver().getCurrentTimeStep(),uint(0));
    QVERIFY(std::fabs(model.getTimeSolver().getComputedTime()-0.1) < 1.0e-14);
}

void TestStructuralCableFE::failedSecondCableStepRollsBackToLastAcceptedState()
{
    RModel model = makeAxialCable(false,false,2.0);
    model.getProblemSetup().getStructuralCableSetup().members[0].initialTension = 0.0;
    RBoundaryCondition &force =
        model.getPoint(1).getBoundaryCondition(R_BOUNDARY_CONDITION_FORCE);
    setComponent(force,R_VARIABLE_FORCE_X,true,-1.6);
    RStructuralDynamicsSetup &dynamics =
        model.getProblemSetup().getStructuralDynamicsSetup();
    dynamics.setEnabled(true);
    dynamics.setTimeStep(1.0);
    dynamics.setNSteps(2);
    dynamics.setOutputFrequency(1);

    RSolver solver(model,QString(),QString());
    bool failed = false;
    QString failureMessage;
    try
    {
        solver.run();
    }
    catch (const RError &error)
    {
        failed = true;
        failureMessage = error.getMessage();
    }
    QVERIFY(failed);
    QVERIFY2(failureMessage.contains("current length must be positive",Qt::CaseInsensitive),
             qPrintable(failureMessage));

    // Step 1 reaches a valid slack state. The constant compression would
    // collapse the element at step 2, so the rejected candidate must restore
    // the already-published step-1 mechanics, kinematics, and accepted time.
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_DISPLACEMENT,0,1)+0.4) < 1.0e-10);
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_VELOCITY,0,1)+0.8) < 1.0e-10);
    QVERIFY(std::fabs(nodeValue(model,R_VARIABLE_ACCELERATION,0,1)+1.6) < 1.0e-10);
    QCOMPARE(model.getTimeSolver().getCurrentTimeStep(),uint(0));
    QVERIFY(std::fabs(model.getTimeSolver().getCurrentTime()-1.0) < 1.0e-14);
    QVERIFY(std::fabs(model.getTimeSolver().getComputedTime()-1.0) < 1.0e-14);
}

void TestStructuralCableFE::cableAndPiezoModesCannotBeEnabledTogether()
{
    RModel model = makeAxialCable();
    model.getProblemSetup().getStructuralPiezoSetup().enabled = true;

    bool rejected = false;
    try
    {
        RSolver solver(model,QString(),QString());
    }
    catch (const RError &error)
    {
        rejected = error.getMessage().contains("cannot be enabled together",Qt::CaseInsensitive);
    }
    QVERIFY(rejected);
    QVERIFY(model.findVariable(R_VARIABLE_DISPLACEMENT) == RConstants::eod);
}

QTEST_MAIN(TestStructuralCableFE)
#include "tst_structural_cable_fe.moc"

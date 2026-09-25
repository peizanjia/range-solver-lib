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

void assignCableMaterial(RLine &line, double area, double youngModulus)
{
    line.setCrossArea(area);
    RMaterial material(RMaterial::Solid);
    RMaterialProperty modulus(RMaterialProperty::ModulusOfElasiticity);
    modulus.add(293.15,youngModulus);
    material.add(modulus);
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
                      bool addContradictoryDisplacement = false)
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
    assignCableMaterial(model.getLine(0),0.5,1000.0);

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
    void unsupportedDynamicsRejectedBeforeCleanup();
    void unconfiguredNonpointMemberRejected();
    void contradictoryPrescriptionsRejected();
    void loadedSlackCableFailsWithoutPublishingResults();
    void directModalCableSolverRejected();
    void cableModeCannotChangeAfterConstruction();
    void directSharedTemperatureMustBeReferenceTemperature();
    void nonReferenceEnvironmentTemperatureRejected();
    void nonzeroEnvironmentGravityRejected();
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

void TestStructuralCableFE::unsupportedDynamicsRejectedBeforeCleanup()
{
    RModel model = makeAxialCable();
    const uint displacementPosition = model.addVariable(R_VARIABLE_DISPLACEMENT);
    model.getVariable(displacementPosition).resize(3,model.getNNodes());
    model.getProblemSetup().getStructuralDynamicsSetup().setEnabled(true);

    bool rejected = false;
    try
    {
        RSolver solver(model,QString(),QString());
    }
    catch (const RError &error)
    {
        rejected = error.getMessage().contains("dynamics",Qt::CaseInsensitive);
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
        rejected = error.getMessage().contains("cannot be enabled or disabled",Qt::CaseInsensitive);
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

QTEST_MAIN(TestStructuralCableFE)
#include "tst_structural_cable_fe.moc"

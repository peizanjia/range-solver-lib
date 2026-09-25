#include <QtTest>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

#include <rbl_error.h>

#include <rml_boundary_condition.h>
#include <rml_line.h>
#include <rml_model.h>
#include <rml_problem_setup.h>
#include <rml_problem_task_item.h>
#include <rml_structural_cable_setup.h>
#include <rml_structural_dynamics_restart_data.h>
#include <rml_structural_piezo_setup.h>
#include <rml_variable.h>

#include "rsolver.h"
#include "rstructural_piezo_model_adapter.h"

namespace
{
void setComponent(RBoundaryCondition &condition,
                  RVariableType type,
                  bool enabled,
                  double value = 0.0)
{
    const uint position = condition.findComponentPosition(type);
    if (position == RConstants::eod)
        throw std::logic_error("The boundary condition is missing a required component.");
    RConditionComponent &component = condition.getComponent(position);
    component.setEnabled(enabled);
    component.clear();
    if (enabled) component.add(0.0,value);
}

RBoundaryCondition displacementCondition(const std::array<bool,3> &enabled,
                                         const std::array<double,3> &values)
{
    RBoundaryCondition condition(R_BOUNDARY_CONDITION_DISPLACEMENT);
    const RVariableType types[3] = {
        R_VARIABLE_DISPLACEMENT_X,R_VARIABLE_DISPLACEMENT_Y,R_VARIABLE_DISPLACEMENT_Z
    };
    for (uint axis=0;axis<3;axis++) setComponent(condition,types[axis],enabled[axis],values[axis]);
    return condition;
}

RBoundaryCondition forceCondition(const std::array<double,3> &values)
{
    RBoundaryCondition condition(R_BOUNDARY_CONDITION_FORCE);
    const RVariableType types[3] = {R_VARIABLE_FORCE_X,R_VARIABLE_FORCE_Y,R_VARIABLE_FORCE_Z};
    for (uint axis=0;axis<3;axis++) setComponent(condition,types[axis],true,values[axis]);
    return condition;
}

void setEnvironmentValue(RElementGroup &group,
                         REnvironmentConditionType conditionType,
                         RVariableType variableType,
                         double value)
{
    REnvironmentCondition &condition = group.getEnvironmentCondition(conditionType);
    const uint position = condition.findComponentPosition(variableType);
    if (position == RConstants::eod)
        throw std::logic_error("The environment condition is missing a required component.");
    RConditionComponent &component = condition.getComponent(position);
    component.setEnabled(true);
    component.clear();
    component.add(0.0,value);
}

void setZeroGravity(RModel &model)
{
    const RVariableType gravityTypes[3] = {
        R_VARIABLE_G_ACCELERATION_X,R_VARIABLE_G_ACCELERATION_Y,R_VARIABLE_G_ACCELERATION_Z
    };
    for (uint groupId=0;groupId<model.getNElementGroups();groupId++)
    {
        RElementGroup *group = model.getElementGroupPtr(groupId);
        if (!group) continue;
        for (RVariableType type : gravityTypes)
            setEnvironmentValue(*group,R_ENVIRONMENT_CONDITION_G_ACCELERATION,type,0.0);
    }
}

void addPointElement(RModel &model,uint node,uint pointGroup)
{
    RElement point(R_ELEMENT_POINT);
    point.setNodeId(0,node);
    model.addElement(point,true,pointGroup);
}

void assignRodMaterial(RLine &line,double area,double youngModulus)
{
    line.setCrossArea(area);
    RMaterial material(RMaterial::Solid);
    RMaterialProperty modulus(RMaterialProperty::ModulusOfElasiticity);
    modulus.add(293.15,youngModulus);
    material.add(modulus);
    line.setMaterial(material);
}

RModel makeRod(double voltage,
               double coupling,
               double permittivity,
               bool reversedEndpoints = false,
               bool prescribeRightAxial = false,
               double rightAxialDisplacement = 0.0,
               double rightAxialForce = 0.0,
               uint segmentCount = 1,
               bool onlyFirstPort = false)
{
    constexpr double segmentLength = 2.0;
    RModel model;
    model.setNNodes(segmentCount+1);
    for (uint node=0;node<=segmentCount;node++)
        model.getNode(node).set(segmentLength*node,0.0,0.0);

    for (uint segment=0;segment<segmentCount;segment++)
    {
        RElement rod(R_ELEMENT_TRUSS1);
        const uint left = reversedEndpoints ? segment+1 : segment;
        const uint right = reversedEndpoints ? segment : segment+1;
        rod.setNodeId(0,left);
        rod.setNodeId(1,right);
        model.addElement(rod,true,0);
    }
    assignRodMaterial(model.getLine(0),0.5,10.0);
    for (uint node=0;node<=segmentCount;node++) addPointElement(model,node,node);

    model.getPoint(0).addBoundaryCondition(
        displacementCondition({{true,true,true}},{{0.0,0.0,0.0}}));
    for (uint node=1;node<segmentCount;node++)
        model.getPoint(node).addBoundaryCondition(
            displacementCondition({{false,true,true}},{{0.0,0.0,0.0}}));
    model.getPoint(segmentCount).addBoundaryCondition(displacementCondition(
        {{prescribeRightAxial,true,true}},
        {{rightAxialDisplacement,0.0,0.0}}));
    if (rightAxialForce != 0.0)
        model.getPoint(segmentCount).addBoundaryCondition(
            forceCondition({{rightAxialForce,0.0,0.0}}));

    setZeroGravity(model);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    RStructuralPiezoSetup &setup = model.getProblemSetup().getStructuralPiezoSetup();
    setup.enabled = true;
    const uint portCount = onlyFirstPort ? std::min(1u,segmentCount) : segmentCount;
    for (uint segment=0;segment<portCount;segment++)
    {
        RStructuralPiezoPort port;
        port.id = QString("port-%1").arg(segment);
        port.elementId = segment;
        port.coupling = reversedEndpoints ? -coupling : coupling;
        port.permittivity = permittivity;
        port.voltage = reversedEndpoints ? -voltage : voltage;
        setup.ports.push_back(port);
    }
    setup.bindToModel(model);
    return model;
}

double nodalValue(const RModel &model,RVariableType type,uint component,uint node)
{
    const uint position = model.findVariable(type);
    if (position == RConstants::eod) return std::numeric_limits<double>::quiet_NaN();
    return model.getVariable(position).getValue(component,node);
}

double elementValue(const RModel &model,RVariableType type,uint element)
{
    const uint position = model.findVariable(type);
    if (position == RConstants::eod) return std::numeric_limits<double>::quiet_NaN();
    return model.getVariable(position).getValue(0,element);
}

bool runSolver(RModel &model,QString &error)
{
    try
    {
        RSolver solver(model,QString(),QString());
        solver.run();
        return true;
    }
    catch (const RError &exception)
    {
        error = exception.getMessage();
    }
    catch (const std::exception &exception)
    {
        error = QString::fromUtf8(exception.what());
    }
    return false;
}

bool adapterRejects(const std::function<void()> &call,QString &message)
{
    try
    {
        call();
    }
    catch (const RError &error)
    {
        message = error.getMessage();
        return true;
    }
    catch (const std::exception &error)
    {
        message = QString::fromUtf8(error.what());
        return true;
    }
    return false;
}
}

class TestStructuralPiezoFE : public QObject
{
    Q_OBJECT

private slots:
    void freeVoltageActuationMatchesRodSolutionAndRightCharge();
    void clampedVoltagePublishesReactionAndCharge();
    void shortCircuitForceAndChargeMatchAnalyticRod();
    void polingReversalChangesTheExpectedSigns();
    void endpointReversalPreservesGlobalResponse();
    void rejectStaleMissingAndMixedSetups();
    void failedLinearSolveDoesNotPublishPiezoResults();
    void publishRejectsNonFiniteChargeBeforeMutation();
};

void TestStructuralPiezoFE::freeVoltageActuationMatchesRodSolutionAndRightCharge()
{
    constexpr double area = 0.5;
    constexpr double length = 2.0;
    constexpr double cE = 10.0;
    constexpr double e = 2.0;
    constexpr double epsS = 4.0;
    constexpr double voltage = 0.03;
    RModel model = makeRod(voltage,e,epsS);
    QString error;
    QVERIFY2(runSolver(model,error),error.toUtf8().constData());

    const double expectedExtension = -e*voltage/cE;
    const double expectedCharge = -(area/length)*(e*expectedExtension-epsS*voltage);
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_DISPLACEMENT,0,1)-expectedExtension) < 1.0e-10);
    QVERIFY(std::abs(elementValue(model,R_VARIABLE_ELECTRIC_CHARGE,0)-expectedCharge) < 1.0e-11);
    QVERIFY(std::abs(elementValue(model,R_VARIABLE_STRESS_NORMAL,0)) < 1.0e-8);
    QCOMPARE(model.getVariable(model.findVariable(R_VARIABLE_ELECTRIC_CHARGE)).getApplyType(),
             R_VARIABLE_APPLY_ELEMENT);
    QCOMPARE(model.getVariable(model.findVariable(R_VARIABLE_ELECTRIC_CHARGE)).getNValues(),
             model.getNElements());
    // The model has two point elements after its rod; charge is only assigned
    // to the explicit port element.
    QCOMPARE(elementValue(model,R_VARIABLE_ELECTRIC_CHARGE,1),0.0);
    QCOMPARE(elementValue(model,R_VARIABLE_ELECTRIC_CHARGE,2),0.0);
}

void TestStructuralPiezoFE::clampedVoltagePublishesReactionAndCharge()
{
    constexpr double area = 0.5;
    constexpr double length = 2.0;
    constexpr double e = 2.0;
    constexpr double epsS = 4.0;
    constexpr double voltage = 0.03;
    RModel model = makeRod(voltage,e,epsS,false,true,0.0);
    QString error;
    QVERIFY2(runSolver(model,error),error.toUtf8().constData());

    const double expectedForce = e*area/length*voltage;
    const double expectedRightCharge = epsS*area/length*voltage;
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_DISPLACEMENT,0,0)) < 1.0e-12);
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_DISPLACEMENT,0,1)) < 1.0e-12);
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_FORCE,0,0)+expectedForce) < 1.0e-10);
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_FORCE,0,1)-expectedForce) < 1.0e-10);
    QVERIFY(std::abs(elementValue(model,R_VARIABLE_ELECTRIC_CHARGE,0)-expectedRightCharge) < 1.0e-11);
    QVERIFY(std::abs(elementValue(model,R_VARIABLE_STRESS_NORMAL,0)-expectedForce/area) < 1.0e-10);
}

void TestStructuralPiezoFE::shortCircuitForceAndChargeMatchAnalyticRod()
{
    constexpr double area = 0.5;
    constexpr double length = 2.0;
    constexpr double cE = 10.0;
    constexpr double e = 2.0;
    constexpr double force = 0.025;
    RModel model = makeRod(0.0,e,4.0,false,false,0.0,force);
    QString error;
    QVERIFY2(runSolver(model,error),error.toUtf8().constData());

    const double expectedExtension = force*length/(cE*area);
    const double expectedRightCharge = -(area/length)*e*expectedExtension;
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_DISPLACEMENT,0,1)-expectedExtension) < 1.0e-10);
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_FORCE,0,0)+force) < 1.0e-10);
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_FORCE,0,1)-force) < 1.0e-10);
    QVERIFY(std::abs(elementValue(model,R_VARIABLE_ELECTRIC_CHARGE,0)-expectedRightCharge) < 1.0e-11);
    QVERIFY(std::abs(elementValue(model,R_VARIABLE_STRESS_NORMAL,0)-force/area) < 1.0e-10);
}

void TestStructuralPiezoFE::polingReversalChangesTheExpectedSigns()
{
    constexpr double imposedExtension = 0.01;
    const auto positiveFixed = RStructuralPiezoModelAdapter::solve(
        makeRod(0.0,2.0,4.0,false,true,imposedExtension));
    const auto negativeFixed = RStructuralPiezoModelAdapter::solve(
        makeRod(0.0,-2.0,4.0,false,true,imposedExtension));
    QVERIFY(std::abs(positiveFixed.memberRightCharges[0]+0.005) < 1.0e-12);
    QVERIFY(std::abs(negativeFixed.memberRightCharges[0]-0.005) < 1.0e-12);
    QVERIFY(std::abs(positiveFixed.memberAxialForces[0]-negativeFixed.memberAxialForces[0]) < 1.0e-12);

    constexpr double voltage = 0.03;
    const auto positiveFree = RStructuralPiezoModelAdapter::solve(
        makeRod(voltage,2.0,4.0));
    const auto negativeFree = RStructuralPiezoModelAdapter::solve(
        makeRod(voltage,-2.0,4.0));
    QVERIFY(std::abs(positiveFree.displacements[3]+0.006) < 1.0e-10);
    QVERIFY(std::abs(negativeFree.displacements[3]-0.006) < 1.0e-10);
    // At voltage-fixed mechanical free extension, the capacitive term and e*u
    // both retain sign, so the electrode charge does not simply reverse.
    QVERIFY(std::abs(positiveFree.memberRightCharges[0]-negativeFree.memberRightCharges[0]) < 1.0e-11);
}

void TestStructuralPiezoFE::endpointReversalPreservesGlobalResponse()
{
    constexpr double voltage = 0.03;
    const auto forward = RStructuralPiezoModelAdapter::solve(
        makeRod(voltage,2.0,4.0,false));
    const auto reversed = RStructuralPiezoModelAdapter::solve(
        makeRod(voltage,2.0,4.0,true));
    for (std::size_t dof=0;dof<forward.displacements.size();dof++)
        QVERIFY(std::abs(forward.displacements[dof]-reversed.displacements[dof]) < 1.0e-12);
    for (std::size_t dof=0;dof<forward.internalForces.size();dof++)
        QVERIFY(std::abs(forward.internalForces[dof]-reversed.internalForces[dof]) < 1.0e-12);
    QVERIFY(std::abs(forward.memberRightCharges[0]+reversed.memberRightCharges[0]) < 1.0e-11);
    QVERIFY(std::abs(forward.enthalpy-reversed.enthalpy) < 1.0e-12);
}

void TestStructuralPiezoFE::rejectStaleMissingAndMixedSetups()
{
    RModel stale = makeRod(0.03,2.0,4.0);
    stale.getNode(1).setX(2.1);
    QString message;
    QVERIFY(adapterRejects([&stale] { RStructuralPiezoModelAdapter::validate(stale); },message));
    QVERIFY(message.contains("different mesh",Qt::CaseInsensitive));

    RModel missing = makeRod(0.01,2.0,4.0,false,false,0.0,0.0,2,true);
    QVERIFY(adapterRejects([&missing] { RStructuralPiezoModelAdapter::validate(missing); },message));
    QVERIFY(message.contains("every TRUSS1 member",Qt::CaseInsensitive));

    RModel mixed = makeRod(0.01,2.0,4.0);
    mixed.getProblemSetup().getStructuralCableSetup().enabled = true;
    QVERIFY(adapterRejects([&mixed] { RStructuralPiezoModelAdapter::validate(mixed); },message));
    QVERIFY(message.contains("cannot be enabled together",Qt::CaseInsensitive));
}

void TestStructuralPiezoFE::failedLinearSolveDoesNotPublishPiezoResults()
{
    RModel model = makeRod(0.03,2.0,4.0);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setNOuterIterations(0);
    QString message;
    QVERIFY(adapterRejects([&model] { (void)RStructuralPiezoModelAdapter::solve(model); },message));
    QVERIFY(message.contains("residual certification",Qt::CaseInsensitive));
    QVERIFY(model.findVariable(R_VARIABLE_ELECTRIC_CHARGE) == RConstants::eod);
    QVERIFY(model.findVariable(R_VARIABLE_DISPLACEMENT) == RConstants::eod);

    QString runError;
    QVERIFY(!runSolver(model,runError));
    QVERIFY(model.findVariable(R_VARIABLE_ELECTRIC_CHARGE) == RConstants::eod);
    QVERIFY(model.findVariable(R_VARIABLE_DISPLACEMENT) == RConstants::eod);
}

void TestStructuralPiezoFE::publishRejectsNonFiniteChargeBeforeMutation()
{
    RModel model = makeRod(0.03,2.0,4.0);
    auto result = RStructuralPiezoModelAdapter::solve(model);
    result.memberRightCharges[0] = std::numeric_limits<double>::infinity();
    QString message;
    QVERIFY(adapterRejects([&model,&result] {
        RStructuralPiezoModelAdapter::publish(model,result);
    },message));
    QVERIFY(message.contains("non-finite",Qt::CaseInsensitive));
    QVERIFY(model.findVariable(R_VARIABLE_ELECTRIC_CHARGE) == RConstants::eod);
    QVERIFY(model.findVariable(R_VARIABLE_DISPLACEMENT) == RConstants::eod);
}

QTEST_MAIN(TestStructuralPiezoFE)
#include "tst_structural_piezo_fe.moc"

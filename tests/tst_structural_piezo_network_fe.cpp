#include <QtTest>

#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

#include <rbl_error.h>

#include <rml_boundary_condition.h>
#include <rml_line.h>
#include <rml_matrix_solver_conf.h>
#include <rml_model.h>
#include <rml_problem_setup.h>
#include <rml_problem_task_item.h>
#include <rml_structural_piezo_setup.h>
#include <rml_variable.h>

#include "rstructural_piezo_model_adapter.h"
#include "rsolver.h"

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
    if (enabled)
        component.add(0.0,value);
}

RBoundaryCondition displacementCondition(const std::array<bool,3> &enabled,
                                         const std::array<double,3> &values)
{
    RBoundaryCondition condition(R_BOUNDARY_CONDITION_DISPLACEMENT);
    const RVariableType types[3] = {
        R_VARIABLE_DISPLACEMENT_X,R_VARIABLE_DISPLACEMENT_Y,R_VARIABLE_DISPLACEMENT_Z
    };
    for (uint axis=0;axis<3;axis++)
        setComponent(condition,types[axis],enabled[axis],values[axis]);
    return condition;
}

RBoundaryCondition forceCondition(const std::array<double,3> &values)
{
    RBoundaryCondition condition(R_BOUNDARY_CONDITION_FORCE);
    const RVariableType types[3] = {R_VARIABLE_FORCE_X,R_VARIABLE_FORCE_Y,R_VARIABLE_FORCE_Z};
    for (uint axis=0;axis<3;axis++)
        setComponent(condition,types[axis],true,values[axis]);
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
        if (!group)
            continue;
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

RModel makeNodalPotentialRod(uint segmentCount,
                             bool shortCircuit,
                             double tipLoad,
                             double tipCharge = 0.0,
                             bool reverseMembers = false,
                             bool transverseTipLoad = false,
                             bool perturbMiddleNode = false,
                             bool transverseTipDisplacement = false)
{
    constexpr double segmentLength = 2.0;
    RModel model;
    model.setNNodes(segmentCount+1);
    for (uint node=0;node<=segmentCount;node++)
    {
        const double y = perturbMiddleNode && node == 1 ? 0.1 : 0.0;
        model.getNode(node).set(segmentLength*node,y,0.0);
    }

    for (uint segment=0;segment<segmentCount;segment++)
    {
        RElement rod(R_ELEMENT_TRUSS1);
        rod.setNodeId(0,reverseMembers ? segment+1 : segment);
        rod.setNodeId(1,reverseMembers ? segment : segment+1);
        model.addElement(rod,true,0);
    }
    assignRodMaterial(model.getLine(0),0.5,10.0);
    for (uint node=0;node<=segmentCount;node++)
        addPointElement(model,node,node);

    model.getPoint(0).addBoundaryCondition(
        displacementCondition({{true,true,true}},{{0.0,0.0,0.0}}));
    for (uint node=1;node<segmentCount;node++)
        model.getPoint(node).addBoundaryCondition(
            displacementCondition({{false,true,true}},{{0.0,0.0,0.0}}));
    model.getPoint(segmentCount).addBoundaryCondition(
        displacementCondition({{false,true,true}},
                              {{0.0,transverseTipDisplacement ? 0.1 : 0.0,0.0}}));
    if (tipLoad != 0.0)
        model.getPoint(segmentCount).addBoundaryCondition(
            forceCondition({{tipLoad,transverseTipLoad ? tipLoad : 0.0,0.0}}));
    else if (transverseTipLoad)
        model.getPoint(segmentCount).addBoundaryCondition(forceCondition({{0.0,1.0,0.0}}));

    setZeroGravity(model);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    RStructuralPiezoSetup &setup = model.getProblemSetup().getStructuralPiezoSetup();
    setup.enabled = true;
    setup.electricalMode = RStructuralPiezoElectricalMode::NodalPotential;
    for (uint segment=0;segment<segmentCount;segment++)
    {
        RStructuralPiezoPort port;
        port.id = QString("network-port-%1").arg(segment);
        port.elementId = segment;
        port.coupling = reverseMembers ? -2.0 : 2.0;
        port.permittivity = 4.0;
        port.voltage = 0.0;
        setup.ports.push_back(port);
    }
    setup.nodalBoundaryConditions.push_back(
        {0,RStructuralPiezoBoundaryType::Potential,0.0});
    if (shortCircuit)
    {
        for (uint node=1;node<=segmentCount;node++)
            setup.nodalBoundaryConditions.push_back(
                {node,RStructuralPiezoBoundaryType::Potential,0.0});
    }
    else
    {
        setup.nodalBoundaryConditions.push_back(
            {segmentCount,RStructuralPiezoBoundaryType::ConductorCharge,tipCharge});
    }
    setup.bindToModel(model);
    RMatrixSolverConf &gmres = model.getMatrixSolverConf(RMatrixSolverConf::GMRES);
    gmres.setNInnerIterations(30);
    gmres.setNOuterIterations(20);
    gmres.setSolverCvgValue(1.0e-12);
    return model;
}

double nodalValue(const RModel &model,RVariableType type,uint component,uint node)
{
    const uint position = model.findVariable(type);
    if (position == RConstants::eod)
        return std::numeric_limits<double>::quiet_NaN();
    return model.getVariable(position).getValue(component,node);
}

bool adapterSolve(RModel &model,RStructuralPiezoModelResult &result,QString &message)
{
    try
    {
        result = RStructuralPiezoModelAdapter::solve(model);
        return true;
    }
    catch (const RError &error)
    {
        message = error.getMessage();
    }
    catch (const std::exception &error)
    {
        message = QString::fromUtf8(error.what());
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

class TestStructuralPiezoNetworkFE : public QObject
{
    Q_OBJECT

private slots:
    void openCircuitRodMatchesAnalyticCoupledSolution();
    void shortCircuitRodMatchesAnalyticForceAndCharge();
    void prescribedConductorChargeMatchesAnalyticCoupledSolution();
    void openCircuitMultiElementNetworkMatchesAnalyticSolution();
    void reversedElementOrientationPreservesGlobalNodalFields();
    void rejectsMissingGaugeNonCollinearGeometryAndTransverseLoads();
    void failedSolveAndInvalidPublishDoNotLeaveNodalOutputs();
    void publishesNodalPotentialAndPhysicalChargeWithCorrectScope();
    void heatThenConstrainedPiezoRodMatchesForceAndCharge();
};

void TestStructuralPiezoNetworkFE::heatThenConstrainedPiezoRodMatchesForceAndCharge()
{
    RModel model = makeNodalPotentialRod(2,true,0.0);
    RMaterial material = model.getLine(0).getMaterial();
    auto addProperty = [&material](RMaterialProperty::Type type,double value) {
        RMaterialProperty property(type);
        property.add(293.15,value);
        material.add(property);
    };
    addProperty(RMaterialProperty::ThermalExpansionCoefficient,1.0e-3);
    addProperty(RMaterialProperty::ThermalConductivity,2.0);
    addProperty(RMaterialProperty::HeatCapacity,500.0);
    addProperty(RMaterialProperty::Density,8000.0);
    addProperty(RMaterialProperty::PoissonRatio,0.0);
    addProperty(RMaterialProperty::Emissivity,0.0);
    model.getLine(0).setMaterial(material);
    model.getPoint(2).addBoundaryCondition(
        displacementCondition({{true,false,false}},{{0.0,0.0,0.0}}));
    for (uint node : {0u,2u})
    {
        RBoundaryCondition temperature(R_BOUNDARY_CONDITION_TEMPERATURE);
        setComponent(temperature,R_VARIABLE_TEMPERATURE,true,393.15);
        model.getPoint(node).addBoundaryCondition(temperature);
    }
    RProblemTaskItem sequence;
    sequence.addChild(RProblemTaskItem(R_PROBLEM_HEAT));
    sequence.addChild(RProblemTaskItem(R_PROBLEM_STRESS));
    model.setProblemTaskTree(sequence);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setSolverCvgValue(1.0e-12);

    RSolver solver(model,QString(),QString());
    try { solver.run(); }
    catch (const RError &error) { QFAIL(qPrintable(error.getMessage())); }
    catch (const std::exception &error) { QFAIL(error.what()); }
    for (uint node=0;node<3;node++)
    {
        QVERIFY(std::abs(nodalValue(model,R_VARIABLE_TEMPERATURE,0,node)-393.15) < 1.0e-8);
        QVERIFY(std::abs(nodalValue(model,R_VARIABLE_DISPLACEMENT,0,node)) < 1.0e-10);
        QVERIFY(std::abs(nodalValue(model,R_VARIABLE_ELECTRIC_POTENTIAL,0,node)) < 1.0e-10);
    }
    // c=10, e=2, A=0.5, alpha*DeltaT=0.1: F=-0.5 N,
    // D=-0.2 C/m^2 and each directed right electrode carries +0.1 C.
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_ELECTRIC_CHARGE,0,0)+0.1) < 1.0e-9);
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_ELECTRIC_CHARGE,0,2)-0.1) < 1.0e-9);
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_STRESS_X,0,0)+1.0) < 1.0e-9);
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_STRESS_X,0,1)+1.0) < 1.0e-9);

    RModel openCircuit = makeNodalPotentialRod(2,false,0.0);
    openCircuit.getLine(0).setMaterial(material);
    for (uint node : {0u,2u})
    {
        RBoundaryCondition temperature(R_BOUNDARY_CONDITION_TEMPERATURE);
        setComponent(temperature,R_VARIABLE_TEMPERATURE,true,393.15);
        openCircuit.getPoint(node).addBoundaryCondition(temperature);
    }
    openCircuit.setProblemTaskTree(sequence);
    openCircuit.getMatrixSolverConf(RMatrixSolverConf::CG).setSolverCvgValue(1.0e-12);
    RSolver openSolver(openCircuit,QString(),QString());
    try { openSolver.run(); }
    catch (const RError &error) { QFAIL(qPrintable(error.getMessage())); }
    catch (const std::exception &error) { QFAIL(error.what()); }
    // Open circuit and traction-free: both the elastic strain and E vanish.
    const double expectedStrain = 0.1;
    QVERIFY(std::abs(nodalValue(openCircuit,R_VARIABLE_DISPLACEMENT,0,2)-
                     4.0*expectedStrain) < 1.0e-8);
    QVERIFY(std::abs(nodalValue(openCircuit,R_VARIABLE_ELECTRIC_POTENTIAL,0,2)) < 1.0e-8);
    QVERIFY(std::abs(nodalValue(openCircuit,R_VARIABLE_ELECTRIC_CHARGE,0,2)) < 1.0e-9);
    QVERIFY(std::abs(nodalValue(openCircuit,R_VARIABLE_STRESS_X,0,0)) < 1.0e-9);
}

void TestStructuralPiezoNetworkFE::openCircuitRodMatchesAnalyticCoupledSolution()
{
    constexpr double area = 0.5;
    constexpr double length = 2.0;
    constexpr double cE = 10.0;
    constexpr double e = 2.0;
    constexpr double epsS = 4.0;
    constexpr double force = 0.025;
    RModel model = makeNodalPotentialRod(1,false,force);

    RStructuralPiezoModelResult result;
    QString error;
    QVERIFY2(adapterSolve(model,result,error),error.toUtf8().constData());

    const double forceScale = force*length/area;
    const double determinant = cE*epsS+e*e;
    const double expectedExtension = forceScale*epsS/determinant;
    const double expectedVoltage = e*forceScale/determinant;
    QVERIFY(std::abs(result.displacements[3]-expectedExtension) < 1.0e-11);
    QVERIFY(std::abs(result.nodalPotentials[1]-expectedVoltage) < 1.0e-11);
    QVERIFY(std::abs(result.memberAxialForces[0]-force) < 1.0e-11);
    QVERIFY(std::abs(result.nodalConductorCharges[1]) < 1.0e-12);
    QVERIFY(std::abs(result.residualNorm) < 1.0e-11);
    QVERIFY(std::abs(result.electricResidualNorm) < 1.0e-12);
}

void TestStructuralPiezoNetworkFE::shortCircuitRodMatchesAnalyticForceAndCharge()
{
    constexpr double area = 0.5;
    constexpr double length = 2.0;
    constexpr double cE = 10.0;
    constexpr double e = 2.0;
    constexpr double force = 0.025;
    RModel model = makeNodalPotentialRod(1,true,force);

    RStructuralPiezoModelResult result;
    QString error;
    QVERIFY2(adapterSolve(model,result,error),error.toUtf8().constData());

    const double expectedExtension = force*length/(area*cE);
    const double expectedRightCharge = -(area/length)*e*expectedExtension;
    QVERIFY(std::abs(result.displacements[3]-expectedExtension) < 1.0e-11);
    QVERIFY(std::abs(result.nodalPotentials[0]) < 1.0e-14);
    QVERIFY(std::abs(result.nodalPotentials[1]) < 1.0e-14);
    QVERIFY(std::abs(result.memberAxialForces[0]-force) < 1.0e-11);
    QVERIFY(std::abs(result.memberRightCharges[0]-expectedRightCharge) < 1.0e-11);
    QVERIFY(std::abs(result.nodalConductorCharges[1]-expectedRightCharge) < 1.0e-11);
}

void TestStructuralPiezoNetworkFE::prescribedConductorChargeMatchesAnalyticCoupledSolution()
{
    constexpr double area = 0.5;
    constexpr double length = 2.0;
    constexpr double cE = 10.0;
    constexpr double e = 2.0;
    constexpr double epsS = 4.0;
    constexpr double force = 0.025;
    constexpr double charge = 1.0e-4;
    RModel model = makeNodalPotentialRod(1,false,force,charge);

    RStructuralPiezoModelResult result;
    QString error;
    QVERIFY2(adapterSolve(model,result,error),error.toUtf8().constData());

    const double forceScale = force*length/area;
    const double chargeScale = charge*length/area;
    const double determinant = cE*epsS+e*e;
    const double expectedExtension = (forceScale*epsS-e*chargeScale)/determinant;
    const double expectedVoltage = (cE*chargeScale+e*forceScale)/determinant;
    QVERIFY(std::abs(result.displacements[3]-expectedExtension) < 1.0e-11);
    QVERIFY(std::abs(result.nodalPotentials[1]-expectedVoltage) < 1.0e-11);
    QVERIFY(std::abs(result.nodalConductorCharges[1]-charge) < 1.0e-12);
}

void TestStructuralPiezoNetworkFE::openCircuitMultiElementNetworkMatchesAnalyticSolution()
{
    constexpr double area = 0.5;
    constexpr double segmentLength = 2.0;
    constexpr double cE = 10.0;
    constexpr double e = 2.0;
    constexpr double epsS = 4.0;
    constexpr double force = 0.025;
    RModel model = makeNodalPotentialRod(2,false,force);

    RStructuralPiezoModelResult result;
    QString error;
    QVERIFY2(adapterSolve(model,result,error),error.toUtf8().constData());

    const double forceScale = force*segmentLength/area;
    const double determinant = cE*epsS+e*e;
    const double segmentExtension = forceScale*epsS/determinant;
    const double segmentVoltage = e*forceScale/determinant;
    QVERIFY(std::abs(result.displacements[3]-segmentExtension) < 1.0e-11);
    QVERIFY(std::abs(result.displacements[6]-2.0*segmentExtension) < 1.0e-11);
    QVERIFY(std::abs(result.nodalPotentials[1]-segmentVoltage) < 1.0e-11);
    QVERIFY(std::abs(result.nodalPotentials[2]-2.0*segmentVoltage) < 1.0e-11);
    for (double value : result.nodalConductorCharges)
        QVERIFY(std::abs(value) < 1.0e-12);
    QVERIFY(std::abs(result.memberAxialForces[0]-force) < 1.0e-11);
    QVERIFY(std::abs(result.memberAxialForces[1]-force) < 1.0e-11);
}

void TestStructuralPiezoNetworkFE::reversedElementOrientationPreservesGlobalNodalFields()
{
    RModel forwardModel = makeNodalPotentialRod(2,false,0.025);
    RModel reversedModel = makeNodalPotentialRod(2,false,0.025,0.0,true);
    RStructuralPiezoModelResult forward;
    RStructuralPiezoModelResult reversed;
    QString error;
    QVERIFY2(adapterSolve(forwardModel,forward,error),error.toUtf8().constData());
    QVERIFY2(adapterSolve(reversedModel,reversed,error),error.toUtf8().constData());
    QCOMPARE(reversed.displacements.size(),forward.displacements.size());
    QCOMPARE(reversed.nodalPotentials.size(),forward.nodalPotentials.size());
    for (std::size_t dof=0;dof<forward.displacements.size();dof++)
        QVERIFY(std::abs(forward.displacements[dof]-reversed.displacements[dof]) < 1.0e-11);
    for (std::size_t node=0;node<forward.nodalPotentials.size();node++)
    {
        QVERIFY(std::abs(forward.nodalPotentials[node]-reversed.nodalPotentials[node]) < 1.0e-11);
        QVERIFY(std::abs(forward.nodalConductorCharges[node]-reversed.nodalConductorCharges[node]) < 1.0e-11);
    }
}

void TestStructuralPiezoNetworkFE::rejectsMissingGaugeNonCollinearGeometryAndTransverseLoads()
{
    RModel noGauge = makeNodalPotentialRod(1,false,0.0);
    noGauge.getProblemSetup().getStructuralPiezoSetup().nodalBoundaryConditions.clear();
    QString message;
    QVERIFY(adapterRejects([&noGauge] { RStructuralPiezoModelAdapter::validate(noGauge); },message));
    QVERIFY(message.contains("potential",Qt::CaseInsensitive) ||
            message.contains("gauge",Qt::CaseInsensitive));

    RModel nonCollinear = makeNodalPotentialRod(2,false,0.0,0.0,false,false,true);
    QVERIFY(adapterRejects([&nonCollinear] { RStructuralPiezoModelAdapter::validate(nonCollinear); },message));
    QVERIFY(message.contains("collinear",Qt::CaseInsensitive));

    RModel transverse = makeNodalPotentialRod(1,false,0.025,0.0,false,true);
    QVERIFY(adapterRejects([&transverse] { RStructuralPiezoModelAdapter::validate(transverse); },message));
    QVERIFY(message.contains("transverse",Qt::CaseInsensitive));

    RModel transverseDisplacement =
        makeNodalPotentialRod(1,false,0.0,0.0,false,false,false,true);
    QVERIFY(adapterRejects([&transverseDisplacement] {
        RStructuralPiezoModelAdapter::validate(transverseDisplacement);
    },message));
    QVERIFY(message.contains("transverse",Qt::CaseInsensitive));
}

void TestStructuralPiezoNetworkFE::failedSolveAndInvalidPublishDoNotLeaveNodalOutputs()
{
    RModel failed = makeNodalPotentialRod(1,false,0.025);
    RStructuralPiezoModelResult accepted;
    QString message;
    QVERIFY2(adapterSolve(failed,accepted,message),message.toUtf8().constData());
    RStructuralPiezoModelAdapter::publish(failed,accepted);
    QVERIFY(failed.findVariable(R_VARIABLE_ELECTRIC_POTENTIAL) != RConstants::eod);
    QVERIFY(failed.findVariable(R_VARIABLE_ELECTRIC_CHARGE) != RConstants::eod);
    RStructuralPiezoModelAdapter::clearResults(failed);
    failed.getMatrixSolverConf(RMatrixSolverConf::GMRES).setNOuterIterations(0);
    RStructuralPiezoModelResult candidate;
    QVERIFY(!adapterSolve(failed,candidate,message));
    QVERIFY(message.contains("GMRES",Qt::CaseInsensitive));
    QVERIFY(failed.findVariable(R_VARIABLE_ELECTRIC_POTENTIAL) == RConstants::eod);
    QVERIFY(failed.findVariable(R_VARIABLE_ELECTRIC_CHARGE) == RConstants::eod);

    RModel invalidPublish = makeNodalPotentialRod(1,false,0.025);
    RStructuralPiezoModelResult result;
    QVERIFY2(adapterSolve(invalidPublish,result,message),message.toUtf8().constData());
    result.nodalPotentials[0] = std::numeric_limits<double>::quiet_NaN();
    QVERIFY(adapterRejects([&invalidPublish,&result] {
        RStructuralPiezoModelAdapter::publish(invalidPublish,result);
    },message));
    QVERIFY(invalidPublish.findVariable(R_VARIABLE_ELECTRIC_POTENTIAL) == RConstants::eod);
    QVERIFY(invalidPublish.findVariable(R_VARIABLE_ELECTRIC_CHARGE) == RConstants::eod);

    RModel changedMode = makeNodalPotentialRod(1,false,0.025);
    QVERIFY2(adapterSolve(changedMode,result,message),message.toUtf8().constData());
    RStructuralPiezoSetup &setup = changedMode.getProblemSetup().getStructuralPiezoSetup();
    setup.electricalMode = RStructuralPiezoElectricalMode::LegacyTerminalVoltage;
    setup.nodalBoundaryConditions.clear();
    setup.bindToModel(changedMode);
    QVERIFY(adapterRejects([&changedMode,&result] {
        RStructuralPiezoModelAdapter::publish(changedMode,result);
    },message));
    QVERIFY(message.contains("mode",Qt::CaseInsensitive));
    QVERIFY(changedMode.findVariable(R_VARIABLE_ELECTRIC_POTENTIAL) == RConstants::eod);
    QVERIFY(changedMode.findVariable(R_VARIABLE_ELECTRIC_CHARGE) == RConstants::eod);
}

void TestStructuralPiezoNetworkFE::publishesNodalPotentialAndPhysicalChargeWithCorrectScope()
{
    RModel model = makeNodalPotentialRod(1,false,0.025);
    RStructuralPiezoModelResult result;
    QString error;
    QVERIFY2(adapterSolve(model,result,error),error.toUtf8().constData());
    RStructuralPiezoModelAdapter::publish(model,result);

    const uint potentialPosition = model.findVariable(R_VARIABLE_ELECTRIC_POTENTIAL);
    const uint chargePosition = model.findVariable(R_VARIABLE_ELECTRIC_CHARGE);
    QVERIFY(potentialPosition != RConstants::eod);
    QVERIFY(chargePosition != RConstants::eod);
    const RVariable &potential = model.getVariable(potentialPosition);
    const RVariable &charge = model.getVariable(chargePosition);
    QCOMPARE(potential.getApplyType(),R_VARIABLE_APPLY_NODE);
    QCOMPARE(charge.getApplyType(),R_VARIABLE_APPLY_NODE);
    QCOMPARE(potential.getNValues(),model.getNNodes());
    QCOMPARE(charge.getNValues(),model.getNNodes());
    QCOMPARE(RVariable::getUnits(R_VARIABLE_ELECTRIC_POTENTIAL),QString("V"));
    QCOMPARE(RVariable::getUnits(R_VARIABLE_ELECTRIC_CHARGE),QString("C"));
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_ELECTRIC_POTENTIAL,0,0)) < 1.0e-14);
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_ELECTRIC_POTENTIAL,0,1)-
                     result.nodalPotentials[1]) < 1.0e-14);
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_ELECTRIC_CHARGE,0,1)-
                     result.nodalConductorCharges[1]) < 1.0e-14);

    RStructuralPiezoModelAdapter::clearResults(model);
    QVERIFY(model.findVariable(R_VARIABLE_ELECTRIC_POTENTIAL) == RConstants::eod);
    QVERIFY(model.findVariable(R_VARIABLE_ELECTRIC_CHARGE) == RConstants::eod);
}

QTEST_APPLESS_MAIN(TestStructuralPiezoNetworkFE)
#include "tst_structural_piezo_network_fe.moc"

#include <QtTest>

#include <cmath>
#include <limits>
#include <stdexcept>

#include <rbl_error.h>
#include <rml_line.h>
#include <rml_boundary_condition.h>
#include <rml_environment_condition.h>
#include <rml_model.h>
#include <rml_problem_task_item.h>
#include <rml_variable.h>

#include "rstructural_thermo_coupled_fe_adapter.h"

namespace
{
RModel makeRodModel(bool addMaterialDynamicViscosity = false)
{
    RModel model;
    model.setNNodes(2);
    model.getNode(0).set(0.0,0.0,0.0);
    model.getNode(1).set(1.0,0.0,0.0);

    RElement member(R_ELEMENT_TRUSS1);
    member.setNodeId(0,0);
    member.setNodeId(1,1);
    model.addElement(member,true,0);

    RMaterial material(RMaterial::Solid);
    material.setName("adapter test solid");
    const auto addProperty = [&material](RMaterialProperty::Type type,double value) {
        RMaterialProperty property(type);
        property.clear();
        property.add(RVariable::getInitValue(R_VARIABLE_TEMPERATURE),value);
        material.add(property);
    };
    addProperty(RMaterialProperty::Density,1.0);
    addProperty(RMaterialProperty::HeatCapacity,5.0);
    addProperty(RMaterialProperty::ThermalConductivity,1.0);
    addProperty(RMaterialProperty::ModulusOfElasiticity,100.0);
    addProperty(RMaterialProperty::ThermalExpansionCoefficient,0.0);
    if (addMaterialDynamicViscosity)
        addProperty(RMaterialProperty::DynamicViscosity,9.0);

    model.getLine(0).setCrossArea(1.0);
    model.getLine(0).setMaterial(material);
    REnvironmentCondition &gravity = model.getLine(0).getEnvironmentCondition(
        R_ENVIRONMENT_CONDITION_G_ACCELERATION);
    for (RVariableType type : {R_VARIABLE_G_ACCELERATION_X,
                               R_VARIABLE_G_ACCELERATION_Y,
                               R_VARIABLE_G_ACCELERATION_Z})
    {
        const uint position = gravity.findComponentPosition(type);
        gravity.getComponent(position).clear();
        gravity.getComponent(position).add(0.0,0.0);
    }
    REnvironmentCondition &ambient = model.getLine(0).getEnvironmentCondition(
        R_ENVIRONMENT_CONDITION_TEMPERATURE);
    const uint ambientPosition = ambient.findComponentPosition(R_VARIABLE_TEMPERATURE);
    ambient.getComponent(ambientPosition).clear();
    ambient.getComponent(ambientPosition).add(0.0,300.0);

    RProblemTaskItem heatStress;
    heatStress.addChild(RProblemTaskItem(R_PROBLEM_HEAT));
    heatStress.addChild(RProblemTaskItem(R_PROBLEM_STRESS));
    model.setProblemTaskTree(heatStress);

    RVariable temperature(R_VARIABLE_TEMPERATURE,R_VARIABLE_APPLY_NODE);
    temperature.resize(1,2,true);
    temperature.setValue(0,0,300.0);
    temperature.setValue(0,1,300.0);
    model.addVariable(temperature);

    RVariable displacement(R_VARIABLE_DISPLACEMENT,R_VARIABLE_APPLY_NODE);
    displacement.resize(3,2,true);
    model.addVariable(displacement);
    for (RVariableType type : {R_VARIABLE_STRESS_VON_MISES,
                               R_VARIABLE_STRESS_NORMAL,
                               R_VARIABLE_STRESS_SHEAR})
    {
        RVariable stale(type,R_VARIABLE_APPLY_ELEMENT);
        stale.resize(1,1,true);
        stale.setValue(0,0,77.0);
        model.addVariable(stale);
    }
    return model;
}

RStructuralThermoCoupled1DConditions fixedExtension(double rightMetres = 0.1)
{
    RStructuralThermoCoupled1DConditions conditions;
    conditions.prescribedDisplacements = {{"0",0.0},{"1",rightMetres}};
    return conditions;
}

RStructuralThermoCoupled1DOptions twoWayOptions()
{
    RStructuralThermoCoupled1DOptions options;
    options.dt = 1.0;
    options.twoWayDissipationHeat = true;
    return options;
}

double value(const RModel &model,RVariableType type,uint vector,uint entity)
{
    const uint position = model.findVariable(type);
    if (position == RConstants::eod)
        return std::numeric_limits<double>::quiet_NaN();
    return model.getVariable(position).getValue(vector,entity);
}
}

class TestStructuralThermoCoupledFEAdapter : public QObject
{
    Q_OBJECT

private slots:
    void explicitKelvinVoigtWorkFeedsHeatAndPublishesAxialFields();
    void failedCertifiedStepLeavesOriginalRModelUntouched();
    void materialDynamicViscosityCannotReplaceExplicitStructuralInput();
    void activeModelConditionsAreRejected();
    void persistedPhysicalIdsDriveDissipationAndRejectMissingInputs();
};

void TestStructuralThermoCoupledFEAdapter::persistedPhysicalIdsDriveDissipationAndRejectMissingInputs()
{
    RModel model = makeRodModel();
    model.setNodePhysicalId(0,41);
    model.setNodePhysicalId(1,97);
    model.setElementPhysicalId(0,501);
    RStructuralThermoCoupled1DConditions conditions;
    conditions.prescribedDisplacements = {{"41",0.0},{"97",0.1}};
    const RStructuralThermoFEViscosityByElement physicalViscosity{{501U,2.0}};
    RStructuralThermoCoupled1DResult result;
    try
    {
        result = RStructuralThermoCoupledFEAdapter::advanceByPhysicalId(
            model,physicalViscosity,conditions,twoWayOptions(),300.0);
    }
    catch (const std::exception &error) { QFAIL(error.what()); }
    catch (const RError &error) { QFAIL(qPrintable(error.getMessage())); }
    QCOMPARE(result.nextState[0].nodeId,QStringLiteral("41"));
    QCOMPARE(result.nextState[1].nodeId,QStringLiteral("97"));
    QCOMPARE(result.nodes[1].nodeId,QStringLiteral("97"));
    QVERIFY(std::abs(result.dissipationHeatEnergyJoules-0.02) < 1.0e-12);
    QCOMPARE(model.getNodePhysicalId(0),41U);
    QCOMPARE(model.getElementPhysicalId(0),501U);
    QCOMPARE(model.getTimeSolver().getComputedTime(),1.0);

    const double acceptedTemperature = value(model,R_VARIABLE_TEMPERATURE,0,1);
    conditions.prescribedDisplacements[1].nodeId = "98";
    QVERIFY_EXCEPTION_THROWN(RStructuralThermoCoupledFEAdapter::advanceByPhysicalId(
        model,physicalViscosity,conditions,twoWayOptions(),300.0),std::invalid_argument);
    QCOMPARE(value(model,R_VARIABLE_TEMPERATURE,0,1),acceptedTemperature);
    QCOMPARE(model.getTimeSolver().getComputedTime(),1.0);
}

void TestStructuralThermoCoupledFEAdapter::explicitKelvinVoigtWorkFeedsHeatAndPublishesAxialFields()
{
    RModel model = makeRodModel();
    const RStructuralThermoFEViscosityByElement viscosity{{0U,2.0}};
    RStructuralThermoCoupled1DResult result;
    try
    {
        result = RStructuralThermoCoupledFEAdapter::advance(
            model,viscosity,fixedExtension(),twoWayOptions(),300.0);
    }
    catch (const std::exception &error)
    {
        QFAIL(error.what());
    }
    catch (const RError &error)
    {
        QFAIL(qPrintable(error.getMessage()));
    }

    QCOMPARE(result.memberDissipatedEnergyJoules.size(),std::size_t(1));
    QVERIFY(std::abs(result.memberDissipatedEnergyJoules[0]-0.02) < 1.0e-12);
    QVERIFY(std::abs(result.dissipationHeatEnergyJoules-0.02) < 1.0e-12);
    QVERIFY(std::abs(result.thermalEnergyChangeJoules-0.02) < 1.0e-10);
    QVERIFY(std::abs(value(model,R_VARIABLE_TEMPERATURE,0,0)-300.004) < 1.0e-10);
    QVERIFY(std::abs(value(model,R_VARIABLE_TEMPERATURE,0,1)-300.004) < 1.0e-10);
    QCOMPARE(model.getVariable(model.findVariable(R_VARIABLE_DISPLACEMENT)).getNVectors(),3U);
    QVERIFY(std::abs(value(model,R_VARIABLE_DISPLACEMENT,0,1)-0.1) < 1.0e-12);
    QVERIFY(std::abs(value(model,R_VARIABLE_DISPLACEMENT,1,1)) < 1.0e-12);
    QVERIFY(std::abs(value(model,R_VARIABLE_DISPLACEMENT,2,1)) < 1.0e-12);

    const uint sigmaX = model.findVariable(R_VARIABLE_STRESS_X);
    QVERIFY(sigmaX != RConstants::eod);
    QCOMPARE(model.getVariable(sigmaX).getApplyType(),R_VARIABLE_APPLY_ELEMENT);
    QCOMPARE(model.getVariable(sigmaX).getNValues(),1U);
    QVERIFY(std::abs(value(model,R_VARIABLE_STRESS_X,0,0)-10.2) < 1.0e-10);
    QVERIFY(value(model,R_VARIABLE_STRESS_Y,0,0) == 0.0);
    QVERIFY(value(model,R_VARIABLE_STRESS_Z,0,0) == 0.0);
    QVERIFY(value(model,R_VARIABLE_STRESS_XY,0,0) == 0.0);
    QVERIFY(value(model,R_VARIABLE_STRESS_XZ,0,0) == 0.0);
    QVERIFY(value(model,R_VARIABLE_STRESS_YZ,0,0) == 0.0);
    QVERIFY(std::abs(value(model,R_VARIABLE_STRESS_VON_MISES,0,0)-10.2) < 1.0e-10);
    QVERIFY(std::abs(value(model,R_VARIABLE_STRESS_NORMAL,0,0)-10.2) < 1.0e-10);
    QVERIFY(value(model,R_VARIABLE_STRESS_SHEAR,0,0) == 0.0);

    const uint variableCountAfterFirstStep = model.getNVariables();
    QCOMPARE(model.getTimeSolver().getCurrentTimeStep(),0U);
    QCOMPARE(model.getTimeSolver().getComputedTime(),1.0);
    RStructuralThermoCoupled1DResult second;
    try
    {
        second = RStructuralThermoCoupledFEAdapter::advance(
            model,viscosity,fixedExtension(),twoWayOptions(),300.0);
    }
    catch (const std::exception &error)
    {
        QFAIL(error.what());
    }
    catch (const RError &error)
    {
        QFAIL(qPrintable(error.getMessage()));
    }
    QVERIFY(std::abs(second.viscousDissipationEnergyJoules) < 1.0e-12);
    QVERIFY(std::abs(value(model,R_VARIABLE_TEMPERATURE,0,0)-300.004) < 1.0e-10);
    QVERIFY(std::abs(value(model,R_VARIABLE_DISPLACEMENT,0,1)-0.1) < 1.0e-12);
    QVERIFY(std::abs(value(model,R_VARIABLE_STRESS_X,0,0)-10.0) < 1.0e-10);
    QVERIFY(std::abs(value(model,R_VARIABLE_STRESS_VON_MISES,0,0)-10.0) < 1.0e-10);
    QCOMPARE(model.getNVariables(),variableCountAfterFirstStep);
    QCOMPARE(model.getTimeSolver().getCurrentTimeStep(),0U);
    QCOMPARE(model.getTimeSolver().getComputedTime(),2.0);
}

void TestStructuralThermoCoupledFEAdapter::failedCertifiedStepLeavesOriginalRModelUntouched()
{
    RModel model = makeRodModel();
    const double oldTemperature0 = value(model,R_VARIABLE_TEMPERATURE,0,0);
    const double oldTemperature1 = value(model,R_VARIABLE_TEMPERATURE,0,1);
    const double oldDisplacement1 = value(model,R_VARIABLE_DISPLACEMENT,0,1);
    const uint oldVariableCount = model.getNVariables();
    RStructuralThermoCoupled1DConditions invalid = fixedExtension();
    invalid.heaters = {{"0",10.01,10.0}};

    QVERIFY_EXCEPTION_THROWN(
        RStructuralThermoCoupledFEAdapter::advance(
            model,{{0U,2.0}},invalid,twoWayOptions(),300.0),
        std::invalid_argument);

    QCOMPARE(value(model,R_VARIABLE_TEMPERATURE,0,0),oldTemperature0);
    QCOMPARE(value(model,R_VARIABLE_TEMPERATURE,0,1),oldTemperature1);
    QCOMPARE(value(model,R_VARIABLE_DISPLACEMENT,0,1),oldDisplacement1);
    QCOMPARE(model.getTimeSolver().getComputedTime(),0.0);
    QCOMPARE(model.getNVariables(),oldVariableCount);
    QCOMPARE(model.findVariable(R_VARIABLE_STRESS_X),RConstants::eod);
    QCOMPARE(value(model,R_VARIABLE_STRESS_VON_MISES,0,0),77.0);
}

void TestStructuralThermoCoupledFEAdapter::materialDynamicViscosityCannotReplaceExplicitStructuralInput()
{
    RModel model = makeRodModel(true);
    const uint oldVariableCount = model.getNVariables();
    QVERIFY_EXCEPTION_THROWN(
        RStructuralThermoCoupledFEAdapter::advance(
            model,{},fixedExtension(),twoWayOptions(),300.0),
        std::invalid_argument);
    QCOMPARE(model.getNVariables(),oldVariableCount);
    QCOMPARE(model.findVariable(R_VARIABLE_STRESS_X),RConstants::eod);
}

void TestStructuralThermoCoupledFEAdapter::activeModelConditionsAreRejected()
{
    RModel modelWithBoundaryCondition = makeRodModel();
    RBoundaryCondition fixed(R_BOUNDARY_CONDITION_DISPLACEMENT);
    const uint xPosition = fixed.findComponentPosition(R_VARIABLE_DISPLACEMENT_X);
    fixed.getComponent(xPosition).setEnabled(true);
    fixed.getComponent(xPosition).clear();
    fixed.getComponent(xPosition).add(0.0,0.0);
    modelWithBoundaryCondition.getLine(0).addBoundaryCondition(fixed);
    const uint boundaryVariableCount = modelWithBoundaryCondition.getNVariables();
    QVERIFY_EXCEPTION_THROWN(
        RStructuralThermoCoupledFEAdapter::advance(
            modelWithBoundaryCondition,{{0U,2.0}},fixedExtension(),twoWayOptions(),300.0),
        std::invalid_argument);
    QCOMPARE(modelWithBoundaryCondition.getNVariables(),boundaryVariableCount);
    QCOMPARE(modelWithBoundaryCondition.getTimeSolver().getComputedTime(),0.0);

    RModel modelWithEnvironment = makeRodModel();
    REnvironmentCondition &ambient = modelWithEnvironment.getLine(0).getEnvironmentCondition(
        R_ENVIRONMENT_CONDITION_TEMPERATURE);
    const uint temperaturePosition = ambient.findComponentPosition(R_VARIABLE_TEMPERATURE);
    ambient.getComponent(temperaturePosition).clear();
    ambient.getComponent(temperaturePosition).add(0.0,301.0);
    const uint environmentVariableCount = modelWithEnvironment.getNVariables();
    QVERIFY_EXCEPTION_THROWN(
        RStructuralThermoCoupledFEAdapter::advance(
            modelWithEnvironment,{{0U,2.0}},fixedExtension(),twoWayOptions(),300.0),
        std::invalid_argument);
    QCOMPARE(modelWithEnvironment.getNVariables(),environmentVariableCount);
    QCOMPARE(modelWithEnvironment.getTimeSolver().getComputedTime(),0.0);
}

QTEST_APPLESS_MAIN(TestStructuralThermoCoupledFEAdapter)
#include "tst_structural_thermo_coupled_fe_adapter.moc"

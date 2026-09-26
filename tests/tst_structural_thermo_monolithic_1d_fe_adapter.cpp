#include <QtTest>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include <rml_boundary_condition.h>
#include <rml_condition_component.h>
#include <rml_environment_condition.h>
#include <rml_line.h>
#include <rml_material.h>
#include <rml_material_property.h>
#include <rml_model.h>
#include <rml_problem_task_item.h>
#include <rml_time_solver.h>
#include <rml_variable.h>

#include "rstructural_thermo_monolithic_1d_fe_adapter.h"

namespace
{
using Vector3 = std::array<double,3>;

RModel makeTwoRodModel(bool temperatureDependentYoungModulus = false)
{
    RModel model;
    model.setNNodes(3);
    model.getNode(0).set(0.0,0.0,0.0);
    model.getNode(1).set(1.0,0.0,0.0);
    model.getNode(2).set(2.0,0.0,0.0);

    const std::pair<uint,uint> connectivity[] = {{0U,1U},{1U,2U}};
    for (const auto &nodes : connectivity)
    {
        RElement member(R_ELEMENT_TRUSS1);
        member.setNodeId(0,nodes.first);
        member.setNodeId(1,nodes.second);
        model.addElement(member,true,0U);
    }

    RMaterial material(RMaterial::Solid);
    material.setName("monolithic adapter test solid");
    const auto addProperty = [&material](RMaterialProperty::Type type,
                                         double value,
                                         double secondValue = std::numeric_limits<double>::quiet_NaN()) {
        RMaterialProperty property(type);
        property.clear();
        property.add(300.0,value);
        if (std::isfinite(secondValue))
        {
            property.add(400.0,secondValue);
        }
        material.add(property);
    };
    addProperty(RMaterialProperty::Density,1.0);
    addProperty(RMaterialProperty::HeatCapacity,10.0);
    addProperty(RMaterialProperty::ThermalConductivity,0.0);
    addProperty(RMaterialProperty::ModulusOfElasiticity,1000.0,
                 temperatureDependentYoungModulus ? 1200.0
                                                   : std::numeric_limits<double>::quiet_NaN());
    addProperty(RMaterialProperty::ThermalExpansionCoefficient,0.001);
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
    temperature.resize(1U,3U,true);
    for (uint node=0;node<3U;node++)
    {
        temperature.setValue(0U,node,300.0);
    }
    model.addVariable(temperature);

    RVariable displacement(R_VARIABLE_DISPLACEMENT,R_VARIABLE_APPLY_NODE);
    displacement.resize(3U,3U,true);
    model.addVariable(displacement);
    for (RVariableType type : {R_VARIABLE_STRESS_VON_MISES,
                               R_VARIABLE_STRESS_NORMAL,
                               R_VARIABLE_STRESS_SHEAR})
    {
        RVariable stale(type,R_VARIABLE_APPLY_ELEMENT);
        stale.resize(1U,2U,true);
        stale.setValue(0U,0U,77.0);
        stale.setValue(0U,1U,77.0);
        model.addVariable(stale);
    }
    return model;
}

RStructuralThermoMonolithic1DConditions endLoadConditions()
{
    RStructuralThermoMonolithic1DConditions conditions;
    conditions.mechanicalForceNewtons = {0.0,0.0,10.0};
    conditions.displacementPrescribed = {true,false,false};
    conditions.prescribedDisplacementMetres = {0.0,0.0,0.0};
    conditions.heatPowerWatts = {0.0,0.0,0.0};
    conditions.temperaturePrescribed = {false,false,false};
    conditions.prescribedTemperatureKelvin = {300.0,300.0,300.0};
    return conditions;
}

RStructuralThermoMonolithic1DOptions coupledOptions()
{
    RStructuralThermoMonolithic1DOptions options;
    options.timeStepSeconds = 2.0;
    options.heatKelvinVoigtDissipation = true;
    options.mechanicalAbsoluteResidualToleranceNewtons = 1.0e-12;
    options.mechanicalRelativeResidualTolerance = 1.0e-13;
    options.thermalAbsoluteResidualToleranceWatts = 1.0e-12;
    options.thermalRelativeResidualTolerance = 1.0e-13;
    return options;
}

RStructuralThermoMonolithicFEViscosityByElement twoRodViscosity()
{
    return {{0U,2000.0},{1U,2000.0}};
}

void orientModel(RModel &model,const Vector3 &direction,bool reverseConnectivity)
{
    for (uint node=0;node<model.getNNodes();node++)
    {
        const double axialCoordinate = model.getNode(node).getX();
        model.getNode(node).set(4.0+axialCoordinate*direction[0],
                                -3.0+axialCoordinate*direction[1],
                                8.0+axialCoordinate*direction[2]);
    }
    if (reverseConnectivity)
    {
        for (uint elementIndex=0;elementIndex<model.getNElements();elementIndex++)
        {
            RElement &element = model.getElement(elementIndex);
            const uint firstNode = element.getNodeId(0U);
            element.setNodeId(0U,element.getNodeId(1U));
            element.setNodeId(1U,firstNode);
        }
    }
}

double value(const RModel &model,RVariableType type,uint vector,uint entity)
{
    const uint position = model.findVariable(type);
    if (position == RConstants::eod)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return model.getVariable(position).getValue(vector,entity);
}

bool near(double actual,double expected,double tolerance = 1.0e-9)
{
    return std::isfinite(actual) &&
        std::abs(actual-expected) <= tolerance*std::max(1.0,std::abs(expected));
}
}

class TestStructuralThermoMonolithic1DFEAdapter : public QObject
{
    Q_OBJECT

private slots:
    void twoRodRModelMatchesIndependentMonolithicResult();
    void arbitraryOrientationAndReversedConnectivityMatchScalarKernel();
    void tiedDominantAxisUsesDeterministicPositiveSign();
    void malformedConditionsFailWithoutChangingAcceptedModel();
    void everyElementRequiresExplicitCurrentIndexViscosity();
    void unsupportedGeometryOrTemperatureDependentMaterialIsRejected();
    void transverseAcceptedDisplacementIsRejectedAtomically();
    void activeModelBoundaryConditionIsRejected();
};

void TestStructuralThermoMonolithic1DFEAdapter::twoRodRModelMatchesIndependentMonolithicResult()
{
    RModel model = makeTwoRodModel();
    const uint firstPhysicalElementId = model.getElementPhysicalId(0U);
    const uint secondPhysicalElementId = model.getElementPhysicalId(1U);
    const uint initialVariableCount = model.getNVariables();
    const auto conditions = endLoadConditions();
    const auto options = coupledOptions();
    RStructuralThermoMonolithic1DResult result;
    try
    {
        result = RStructuralThermoMonolithic1DFEAdapter::advance(
            model,twoRodViscosity(),conditions,options,300.0);
    }
    catch (const std::exception &error)
    {
        QFAIL(error.what());
    }

    const double young = 1000.0;
    const double expansion = 0.001;
    const double densityCapacity = 10.0;
    const double viscosity = 2000.0;
    const double referenceTemperature = 300.0;
    const double length = 1.0;
    const double area = 1.0;
    const double timeStep = options.timeStepSeconds;
    const double forcePerArea = conditions.mechanicalForceNewtons.back()/area;
    const double quadratic = young*expansion*viscosity/
        (densityCapacity*length*length*timeStep);
    const double linear = (young+viscosity/timeStep)/length+
        referenceTemperature*young*young*expansion*expansion/
            (densityCapacity*length);
    const double discriminant = linear*linear-4.0*quadratic*forcePerArea;
    QVERIFY(discriminant > 0.0);
    const double memberExtension = 2.0*forcePerArea/(linear+std::sqrt(discriminant));
    const double temperatureChange =
        (viscosity*memberExtension*memberExtension/(length*timeStep)-
         referenceTemperature*young*expansion*memberExtension)/
        (densityCapacity*length);
    const double memberDissipation = viscosity*area*length*
        std::pow(memberExtension/(length*timeStep),2)*timeStep;

    QCOMPARE(result.nextState.temperatureKelvin.size(),std::size_t(3));
    QCOMPARE(result.nextState.axialDisplacementMetres.size(),std::size_t(3));
    for (uint node=0;node<3U;node++)
    {
        const double expectedDisplacement = double(node)*memberExtension;
        QVERIFY(near(result.nextState.temperatureKelvin[node],
                     referenceTemperature+temperatureChange));
        QVERIFY(near(result.nextState.axialDisplacementMetres[node],expectedDisplacement));
        QVERIFY(near(value(model,R_VARIABLE_TEMPERATURE,0U,node),
                     referenceTemperature+temperatureChange));
        QVERIFY(near(value(model,R_VARIABLE_DISPLACEMENT,0U,node),expectedDisplacement));
        QCOMPARE(value(model,R_VARIABLE_DISPLACEMENT,1U,node),0.0);
        QCOMPARE(value(model,R_VARIABLE_DISPLACEMENT,2U,node),0.0);
    }
    for (uint element=0;element<2U;element++)
    {
        QVERIFY(near(result.memberAxialStressPascals[element],10.0));
        QVERIFY(near(result.memberDissipatedEnergyJoules[element],memberDissipation));
        QVERIFY(near(value(model,R_VARIABLE_STRESS_X,0U,element),10.0));
        QCOMPARE(value(model,R_VARIABLE_STRESS_Y,0U,element),0.0);
        QCOMPARE(value(model,R_VARIABLE_STRESS_Z,0U,element),0.0);
        QCOMPARE(value(model,R_VARIABLE_STRESS_XY,0U,element),0.0);
        QCOMPARE(value(model,R_VARIABLE_STRESS_XZ,0U,element),0.0);
        QCOMPARE(value(model,R_VARIABLE_STRESS_YZ,0U,element),0.0);
        QVERIFY(near(value(model,R_VARIABLE_STRESS_VON_MISES,0U,element),10.0));
        QVERIFY(near(value(model,R_VARIABLE_STRESS_NORMAL,0U,element),10.0));
        QCOMPARE(value(model,R_VARIABLE_STRESS_SHEAR,0U,element),0.0);
    }
    QVERIFY(near(result.thermalEnergyChangeJoules,
                 densityCapacity*2.0*length*area*temperatureChange));
    QVERIFY(near(result.reversibleThermoelasticEnergyJoules,
                 2.0*referenceTemperature*young*expansion*area*memberExtension));
    QVERIFY(near(result.viscousDissipationEnergyJoules,2.0*memberDissipation));
    QVERIFY(near(result.dissipationHeatEnergyJoules,2.0*memberDissipation));
    QVERIFY(std::abs(result.thermalEnergyBalanceResidualJoules) <= 1.0e-9);
    QCOMPARE(model.getTimeSolver().getComputedTime(),timeStep);
    QCOMPARE(model.getTimeSolver().getCurrentTimeStep(),0U);
    QCOMPARE(model.getElementPhysicalId(0U),firstPhysicalElementId);
    QCOMPARE(model.getElementPhysicalId(1U),secondPhysicalElementId);
    QCOMPARE(model.getNVariables(),initialVariableCount+6U);
}

void TestStructuralThermoMonolithic1DFEAdapter::arbitraryOrientationAndReversedConnectivityMatchScalarKernel()
{
    RModel axialModel = makeTwoRodModel();
    RModel rotatedModel = makeTwoRodModel();
    const Vector3 physicalDirection = {-2.0/3.0,1.0/3.0,2.0/3.0};
    orientModel(rotatedModel,physicalDirection,true);

    RStructuralThermoMonolithic1DResult axialResult;
    RStructuralThermoMonolithic1DResult rotatedResult;
    auto axialConditions = endLoadConditions();
    axialConditions.prescribedDisplacementMetres[0U] = 0.02;
    auto rotatedConditions = endLoadConditions();
    rotatedConditions.mechanicalForceNewtons.back() = -10.0;
    rotatedConditions.prescribedDisplacementMetres[0U] = -0.02;
    // Second advance reprojects the nonzero accepted global displacement.
    for (unsigned int acceptedStep = 0; acceptedStep < 2U; ++acceptedStep)
    {
    try
    {
        axialResult = RStructuralThermoMonolithic1DFEAdapter::advance(
            axialModel,twoRodViscosity(),axialConditions,coupledOptions(),300.0);
        rotatedResult = RStructuralThermoMonolithic1DFEAdapter::advance(
            rotatedModel,twoRodViscosity(),rotatedConditions,coupledOptions(),300.0);
    }
    catch (const std::exception &error)
    {
        QFAIL(error.what());
    }

    for (uint node=0;node<axialModel.getNNodes();node++)
    {
        QVERIFY(near(rotatedResult.nextState.temperatureKelvin[node],
                     axialResult.nextState.temperatureKelvin[node]));
        QVERIFY(near(rotatedResult.nextState.axialDisplacementMetres[node],
                     -axialResult.nextState.axialDisplacementMetres[node]));
        QVERIFY(near(value(rotatedModel,R_VARIABLE_TEMPERATURE,0U,node),
                     value(axialModel,R_VARIABLE_TEMPERATURE,0U,node)));
        for (uint component=0;component<3U;component++)
        {
            QVERIFY(near(value(rotatedModel,R_VARIABLE_DISPLACEMENT,component,node),
                         physicalDirection[component]*
                             axialResult.nextState.axialDisplacementMetres[node]));
        }
    }
    for (uint element=0;element<axialModel.getNElements();element++)
    {
        QVERIFY(near(rotatedResult.memberAxialStressPascals[element],
                     axialResult.memberAxialStressPascals[element]));
        QVERIFY(near(rotatedResult.memberDissipatedEnergyJoules[element],
                     axialResult.memberDissipatedEnergyJoules[element]));
        QVERIFY(near(value(rotatedModel,R_VARIABLE_STRESS_X,0U,element),
                     value(axialModel,R_VARIABLE_STRESS_X,0U,element)));
    }
    QVERIFY(near(rotatedResult.thermalEnergyChangeJoules,
                 axialResult.thermalEnergyChangeJoules));
    QVERIFY(near(rotatedResult.heaterEnergyJoules,axialResult.heaterEnergyJoules));
    QVERIFY(near(rotatedResult.reversibleThermoelasticEnergyJoules,
                 axialResult.reversibleThermoelasticEnergyJoules));
    QVERIFY(near(rotatedResult.viscousDissipationEnergyJoules,
                 axialResult.viscousDissipationEnergyJoules));
    QVERIFY(near(rotatedResult.dissipationHeatEnergyJoules,
                 axialResult.dissipationHeatEnergyJoules));
    QVERIFY(near(rotatedResult.thermalBoundaryEnergyJoules,
                 axialResult.thermalBoundaryEnergyJoules));
    QVERIFY(near(rotatedResult.thermalEnergyBalanceResidualJoules,
                 axialResult.thermalEnergyBalanceResidualJoules));
    }
}

void TestStructuralThermoMonolithic1DFEAdapter::tiedDominantAxisUsesDeterministicPositiveSign()
{
    RModel model = makeTwoRodModel();
    const double diagonal = std::sqrt(0.5);
    const Vector3 physicalDirection = {-diagonal,diagonal,0.0};
    orientModel(model,physicalDirection,false);
    auto conditions = endLoadConditions();
    conditions.mechanicalForceNewtons.back() = -10.0;

    RStructuralThermoMonolithic1DResult result;
    try
    {
        result = RStructuralThermoMonolithic1DFEAdapter::advance(
            model,twoRodViscosity(),conditions,coupledOptions(),300.0);
    }
    catch (const std::exception &error)
    {
        QFAIL(error.what());
    }

    // X and Y tie for the dominant component. X wins the tie and is made
    // positive, so the scalar state has the opposite sign to this geometry.
    QVERIFY(result.nextState.axialDisplacementMetres[2U] < 0.0);
    for (uint component=0;component<3U;component++)
    {
        QVERIFY(near(value(model,R_VARIABLE_DISPLACEMENT,component,2U),
                     physicalDirection[component]*
                         std::abs(result.nextState.axialDisplacementMetres[2U])));
    }
}

void TestStructuralThermoMonolithic1DFEAdapter::malformedConditionsFailWithoutChangingAcceptedModel()
{
    RModel model = makeTwoRodModel();
    const double oldTemperature = value(model,R_VARIABLE_TEMPERATURE,0U,2U);
    const double oldDisplacement = value(model,R_VARIABLE_DISPLACEMENT,0U,2U);
    const double oldStress = value(model,R_VARIABLE_STRESS_VON_MISES,0U,0U);
    const uint oldVariableCount = model.getNVariables();
    auto conditions = endLoadConditions();
    conditions.heatPowerWatts.pop_back();

    QVERIFY_EXCEPTION_THROWN(RStructuralThermoMonolithic1DFEAdapter::advance(
        model,twoRodViscosity(),conditions,coupledOptions(),300.0),std::invalid_argument);
    QCOMPARE(value(model,R_VARIABLE_TEMPERATURE,0U,2U),oldTemperature);
    QCOMPARE(value(model,R_VARIABLE_DISPLACEMENT,0U,2U),oldDisplacement);
    QCOMPARE(value(model,R_VARIABLE_STRESS_VON_MISES,0U,0U),oldStress);
    QCOMPARE(model.getNVariables(),oldVariableCount);
    QCOMPARE(model.getTimeSolver().getComputedTime(),0.0);

    auto limitedOptions = coupledOptions();
    limitedOptions.maximumNewtonIterations = 1U;
    QVERIFY_EXCEPTION_THROWN(RStructuralThermoMonolithic1DFEAdapter::advance(
        model,twoRodViscosity(),endLoadConditions(),limitedOptions,300.0),std::runtime_error);
    QCOMPARE(value(model,R_VARIABLE_TEMPERATURE,0U,2U),oldTemperature);
    QCOMPARE(value(model,R_VARIABLE_DISPLACEMENT,0U,2U),oldDisplacement);
    QCOMPARE(value(model,R_VARIABLE_STRESS_VON_MISES,0U,0U),oldStress);
    QCOMPARE(model.getNVariables(),oldVariableCount);
    QCOMPARE(model.getTimeSolver().getComputedTime(),0.0);
}

void TestStructuralThermoMonolithic1DFEAdapter::everyElementRequiresExplicitCurrentIndexViscosity()
{
    RModel model = makeTwoRodModel();
    const double acceptedTemperature = value(model,R_VARIABLE_TEMPERATURE,0U,1U);
    const double acceptedDisplacement = value(model,R_VARIABLE_DISPLACEMENT,0U,2U);
    const uint variableCount = model.getNVariables();
    const RStructuralThermoMonolithicFEViscosityByElement incomplete{{0U,2000.0}};

    QVERIFY_EXCEPTION_THROWN(RStructuralThermoMonolithic1DFEAdapter::advance(
        model,incomplete,endLoadConditions(),coupledOptions(),300.0),std::invalid_argument);
    QCOMPARE(value(model,R_VARIABLE_TEMPERATURE,0U,1U),acceptedTemperature);
    QCOMPARE(value(model,R_VARIABLE_DISPLACEMENT,0U,2U),acceptedDisplacement);
    QCOMPARE(model.getNVariables(),variableCount);
    QCOMPARE(model.getTimeSolver().getComputedTime(),0.0);
}

void TestStructuralThermoMonolithic1DFEAdapter::unsupportedGeometryOrTemperatureDependentMaterialIsRejected()
{
    RModel misaligned = makeTwoRodModel();
    misaligned.getNode(2U).set(2.0,0.01,0.0);
    const double acceptedTemperature = value(misaligned,R_VARIABLE_TEMPERATURE,0U,2U);
    const double acceptedDisplacement = value(misaligned,R_VARIABLE_DISPLACEMENT,0U,2U);
    const double acceptedStress = value(misaligned,R_VARIABLE_STRESS_VON_MISES,0U,0U);
    const uint variableCount = misaligned.getNVariables();
    QVERIFY_EXCEPTION_THROWN(RStructuralThermoMonolithic1DFEAdapter::advance(
        misaligned,twoRodViscosity(),endLoadConditions(),coupledOptions(),300.0),
        std::invalid_argument);
    QCOMPARE(value(misaligned,R_VARIABLE_TEMPERATURE,0U,2U),acceptedTemperature);
    QCOMPARE(value(misaligned,R_VARIABLE_DISPLACEMENT,0U,2U),acceptedDisplacement);
    QCOMPARE(value(misaligned,R_VARIABLE_STRESS_VON_MISES,0U,0U),acceptedStress);
    QCOMPARE(misaligned.getNVariables(),variableCount);
    QCOMPARE(misaligned.getTimeSolver().getComputedTime(),0.0);

    RModel temperatureDependent = makeTwoRodModel(true);
    QVERIFY_EXCEPTION_THROWN(RStructuralThermoMonolithic1DFEAdapter::advance(
        temperatureDependent,twoRodViscosity(),endLoadConditions(),coupledOptions(),300.0),
        std::invalid_argument);
    QCOMPARE(value(temperatureDependent,R_VARIABLE_TEMPERATURE,0U,2U),300.0);
    QCOMPARE(temperatureDependent.getTimeSolver().getComputedTime(),0.0);
}

void TestStructuralThermoMonolithic1DFEAdapter::transverseAcceptedDisplacementIsRejectedAtomically()
{
    RModel model = makeTwoRodModel();
    const uint displacementPosition = model.findVariable(R_VARIABLE_DISPLACEMENT);
    model.getVariable(displacementPosition).setValue(1U,1U,0.01);
    const double acceptedTemperature = value(model,R_VARIABLE_TEMPERATURE,0U,2U);
    const double acceptedTransverseDisplacement = value(model,R_VARIABLE_DISPLACEMENT,1U,1U);
    const double acceptedStress = value(model,R_VARIABLE_STRESS_VON_MISES,0U,0U);
    const uint variableCount = model.getNVariables();

    QVERIFY_EXCEPTION_THROWN(RStructuralThermoMonolithic1DFEAdapter::advance(
        model,twoRodViscosity(),endLoadConditions(),coupledOptions(),300.0),
        std::invalid_argument);
    QCOMPARE(value(model,R_VARIABLE_TEMPERATURE,0U,2U),acceptedTemperature);
    QCOMPARE(value(model,R_VARIABLE_DISPLACEMENT,1U,1U),acceptedTransverseDisplacement);
    QCOMPARE(value(model,R_VARIABLE_STRESS_VON_MISES,0U,0U),acceptedStress);
    QCOMPARE(model.getNVariables(),variableCount);
    QCOMPARE(model.getTimeSolver().getComputedTime(),0.0);
}

void TestStructuralThermoMonolithic1DFEAdapter::activeModelBoundaryConditionIsRejected()
{
    RModel model = makeTwoRodModel();
    RBoundaryCondition fixed(R_BOUNDARY_CONDITION_DISPLACEMENT);
    const uint xPosition = fixed.findComponentPosition(R_VARIABLE_DISPLACEMENT_X);
    fixed.getComponent(xPosition).setEnabled(true);
    fixed.getComponent(xPosition).clear();
    fixed.getComponent(xPosition).add(0.0,0.0);
    model.getLine(0U).addBoundaryCondition(fixed);

    QVERIFY_EXCEPTION_THROWN(RStructuralThermoMonolithic1DFEAdapter::advance(
        model,twoRodViscosity(),endLoadConditions(),coupledOptions(),300.0),
        std::invalid_argument);
    QCOMPARE(value(model,R_VARIABLE_TEMPERATURE,0U,2U),300.0);
    QCOMPARE(value(model,R_VARIABLE_STRESS_VON_MISES,0U,0U),77.0);
    QCOMPARE(model.getTimeSolver().getComputedTime(),0.0);
}

QTEST_APPLESS_MAIN(TestStructuralThermoMonolithic1DFEAdapter)
#include "tst_structural_thermo_monolithic_1d_fe_adapter.moc"

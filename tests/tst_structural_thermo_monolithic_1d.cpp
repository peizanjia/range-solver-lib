#include <QtTest>

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "rstructural_thermo_monolithic_1d.h"

namespace
{
RStructuralThermoMonolithic1DModel makeModel(double length = 1.0,
                                             double area = 1.0,
                                             double youngModulus = 1000.0,
                                             double expansion = 0.01,
                                             double density = 1.0,
                                             double specificHeat = 10.0,
                                             double conductivity = 0.0,
                                             double viscosity = 0.0)
{
    RStructuralThermoMonolithic1DModel model;
    model.nodeCoordinatesMetres = {0.0,length};
    model.referenceTemperatureKelvin = 300.0;
    RStructuralThermoMonolithic1DMember member;
    member.node0 = 0;
    member.node1 = 1;
    member.areaSquareMetres = area;
    member.youngModulusPascals = youngModulus;
    member.expansionPerKelvin = expansion;
    member.densityKilogramsPerCubicMetre = density;
    member.specificHeatJoulesPerKilogramKelvin = specificHeat;
    member.conductivityWattsPerMetreKelvin = conductivity;
    member.kelvinVoigtViscosityPascalSeconds = viscosity;
    model.members = {member};
    return model;
}

RStructuralThermoMonolithic1DModel makeTwoMemberSeriesModel()
{
    auto model = makeModel(1.0,1.0,1000.0,0.001,1.0,10.0,0.0,2000.0);
    model.nodeCoordinatesMetres = {0.0,1.0,2.0};
    RStructuralThermoMonolithic1DMember left = model.members.front();
    RStructuralThermoMonolithic1DMember right = left;
    left.node0 = 0;
    left.node1 = 1;
    right.node0 = 1;
    right.node1 = 2;
    model.members = {left,right};
    return model;
}

RStructuralThermoMonolithic1DState uniformState(double temperature = 300.0)
{
    return {{temperature,temperature},{0.0,0.0}};
}

RStructuralThermoMonolithic1DConditions conditions()
{
    RStructuralThermoMonolithic1DConditions value;
    value.mechanicalForceNewtons = {0.0,0.0};
    value.displacementPrescribed = {false,false};
    value.prescribedDisplacementMetres = {0.0,0.0};
    value.heatPowerWatts = {0.0,0.0};
    value.temperaturePrescribed = {false,false};
    value.prescribedTemperatureKelvin = {300.0,300.0};
    return value;
}

RStructuralThermoMonolithic1DOptions options(double timeStep = 1.0)
{
    RStructuralThermoMonolithic1DOptions value;
    value.timeStepSeconds = timeStep;
    return value;
}

bool near(double actual, double expected, double tolerance = 1.0e-9)
{
    return std::isfinite(actual) &&
        std::abs(actual-expected) <= tolerance*std::max(1.0,std::abs(expected));
}
}

class TestStructuralThermoMonolithic1D : public QObject
{
    Q_OBJECT
private slots:
    void adiabaticTensionMatchesThermoelasticCoolingFormula();
    void adiabaticCompressionWarms();
    void zeroExpansionDecouplesThermalAndMechanicalFields();
    void prescribedTemperatureReportsBoundaryHeat();
    void optionalKelvinVoigtDissipationConservesHeat();
    void freeTwoMemberThermoelasticDissipationMatchesIndependentRoot();
    void overflowingNewtonTrialBacktracksToFiniteCandidate();
    void failedNewtonPublishesNoCandidate();
    void unattainableCoolingNeverPublishesNonpositiveTemperature();
};

void TestStructuralThermoMonolithic1D::adiabaticTensionMatchesThermoelasticCoolingFormula()
{
    const auto model = makeModel();
    const auto accepted = uniformState();
    auto load = conditions();
    load.displacementPrescribed[0] = true;
    load.mechanicalForceNewtons[1] = 1.0;

    const auto result = RStructuralThermoMonolithic1D::step(
        model,accepted,load,options());
    const double coupling = model.referenceTemperatureKelvin*
        model.members[0].youngModulusPascals*
        model.members[0].expansionPerKelvin*
        model.members[0].expansionPerKelvin/
        (model.members[0].densityKilogramsPerCubicMetre*
         model.members[0].specificHeatJoulesPerKilogramKelvin);
    const double expectedDisplacement =
        (1.0*1.0/(model.members[0].youngModulusPascals*
                  model.members[0].areaSquareMetres))/(1.0+coupling);
    const double expectedTemperatureChange =
        -model.referenceTemperatureKelvin*model.members[0].youngModulusPascals*
         model.members[0].expansionPerKelvin*expectedDisplacement/
         (model.members[0].densityKilogramsPerCubicMetre*
          model.members[0].specificHeatJoulesPerKilogramKelvin);

    QVERIFY(near(result.nextState.axialDisplacementMetres[1],expectedDisplacement));
    QVERIFY(near(result.nextState.temperatureKelvin[0],300.0+expectedTemperatureChange));
    QVERIFY(near(result.nextState.temperatureKelvin[1],300.0+expectedTemperatureChange));
    QVERIFY(near(result.memberAxialForceNewtons[0],1.0));
    QVERIFY(result.nextState.temperatureKelvin[0] > 0.0);
    QVERIFY(result.mechanicalResidualNormNewtons <= 1.0e-8);
    QVERIFY(result.thermalResidualNormWatts <= 1.0e-8);
    QVERIFY(std::abs(result.thermalEnergyBalanceResidualJoules) <= 1.0e-8);
    QVERIFY(near(result.reversibleThermoelasticEnergyJoules,
                 -result.thermalEnergyChangeJoules));
}

void TestStructuralThermoMonolithic1D::adiabaticCompressionWarms()
{
    const auto model = makeModel();
    const auto accepted = uniformState();
    auto load = conditions();
    load.displacementPrescribed[0] = true;
    load.mechanicalForceNewtons[1] = -1.0;

    const auto result = RStructuralThermoMonolithic1D::step(
        model,accepted,load,options());
    QVERIFY(near(result.nextState.axialDisplacementMetres[1],-0.00025));
    QVERIFY(near(result.nextState.temperatureKelvin[0],300.075));
    QVERIFY(near(result.nextState.temperatureKelvin[1],300.075));
    QVERIFY(result.nextState.temperatureKelvin[0] > accepted.temperatureKelvin[0]);
    QVERIFY(near(result.memberAxialForceNewtons[0],-1.0));
    QVERIFY(std::abs(result.thermalEnergyBalanceResidualJoules) <= 1.0e-8);
}

void TestStructuralThermoMonolithic1D::zeroExpansionDecouplesThermalAndMechanicalFields()
{
    const auto model = makeModel(2.0,1.0,1000.0,0.0,1.0,10.0,1.0,0.0);
    const auto accepted = uniformState();
    auto load = conditions();
    load.displacementPrescribed[0] = true;
    load.mechanicalForceNewtons[1] = 5.0;
    load.heatPowerWatts[0] = 10.0;

    const auto result = RStructuralThermoMonolithic1D::step(
        model,accepted,load,options());
    QVERIFY(near(result.nextState.axialDisplacementMetres[1],0.01));
    QVERIFY(near(result.memberAxialForceNewtons[0],5.0));
    QVERIFY(near(result.heaterEnergyJoules,10.0));
    QVERIFY(near(result.thermalEnergyChangeJoules,10.0));
    QVERIFY(near(result.reversibleThermoelasticEnergyJoules,0.0));
    QVERIFY(std::abs(result.thermalEnergyBalanceResidualJoules) <= 1.0e-8);
}

void TestStructuralThermoMonolithic1D::prescribedTemperatureReportsBoundaryHeat()
{
    const auto model = makeModel(1.0,1.0,1000.0,0.0,1.0,10.0,10.0,0.0);
    const auto accepted = uniformState();
    auto boundary = conditions();
    boundary.displacementPrescribed = {true,true};
    boundary.temperaturePrescribed[0] = true;
    boundary.prescribedTemperatureKelvin[0] = 310.0;

    const auto result = RStructuralThermoMonolithic1D::step(
        model,accepted,boundary,options());
    QVERIFY(near(result.nextState.temperatureKelvin[0],310.0));
    QVERIFY(near(result.nextState.temperatureKelvin[1],306.25));
    QVERIFY(near(result.thermalEnergyChangeJoules,81.25));
    QVERIFY(near(result.thermalBoundaryEnergyJoules,81.25));
    QVERIFY(std::abs(result.thermalEnergyBalanceResidualJoules) <= 1.0e-8);
}

void TestStructuralThermoMonolithic1D::optionalKelvinVoigtDissipationConservesHeat()
{
    const auto model = makeModel(2.0,1.0,1000.0,0.0,1.0,10.0,0.0,4.0);
    const auto accepted = uniformState();
    auto load = conditions();
    load.displacementPrescribed = {true,true};
    load.prescribedDisplacementMetres = {0.0,0.2};

    auto oneWay = options(2.0);
    const auto withoutHeatFeedback = RStructuralThermoMonolithic1D::step(
        model,accepted,load,oneWay);
    QVERIFY(near(withoutHeatFeedback.viscousDissipationEnergyJoules,0.04));
    QVERIFY(near(withoutHeatFeedback.dissipationHeatEnergyJoules,0.0));
    QVERIFY(near(withoutHeatFeedback.thermalEnergyChangeJoules,0.0));

    oneWay.heatKelvinVoigtDissipation = true;
    const auto withHeatFeedback = RStructuralThermoMonolithic1D::step(
        model,accepted,load,oneWay);
    QVERIFY(near(withHeatFeedback.viscousDissipationEnergyJoules,0.04));
    QVERIFY(near(withHeatFeedback.dissipationHeatEnergyJoules,0.04));
    QVERIFY(near(withHeatFeedback.memberDissipatedEnergyJoules[0],0.04));
    QVERIFY(near(withHeatFeedback.thermalEnergyChangeJoules,0.04));
    QVERIFY(near(withHeatFeedback.nextState.temperatureKelvin[0],300.002));
    QVERIFY(near(withHeatFeedback.nextState.temperatureKelvin[1],300.002));
    QVERIFY(std::abs(withHeatFeedback.thermalEnergyBalanceResidualJoules) <= 1.0e-10);
}

void TestStructuralThermoMonolithic1D::freeTwoMemberThermoelasticDissipationMatchesIndependentRoot()
{
    const auto model = makeTwoMemberSeriesModel();
    const RStructuralThermoMonolithic1DState accepted{
        {300.0,300.0,300.0},{0.0,0.0,0.0}};
    RStructuralThermoMonolithic1DConditions load;
    load.mechanicalForceNewtons = {0.0,0.0,10.0};
    load.displacementPrescribed = {true,false,false};
    load.prescribedDisplacementMetres = {0.0,0.0,0.0};
    load.heatPowerWatts = {0.0,0.0,0.0};
    load.temperaturePrescribed = {false,false,false};
    load.prescribedTemperatureKelvin = {300.0,300.0,300.0};
    auto stepOptions = options(2.0);
    stepOptions.heatKelvinVoigtDissipation = true;
    stepOptions.mechanicalAbsoluteResidualToleranceNewtons = 1.0e-12;
    stepOptions.mechanicalRelativeResidualTolerance = 1.0e-13;
    stepOptions.thermalAbsoluteResidualToleranceWatts = 1.0e-12;
    stepOptions.thermalRelativeResidualTolerance = 1.0e-13;

    const auto result = RStructuralThermoMonolithic1D::step(
        model,accepted,load,stepOptions);

    // An independent scalar reduction follows from equal member force/strain
    // and uniform temperature in the two-member series chain. It checks the
    // shared-node assembly and the thermal-to-free-displacement Newton block.
    const auto &member = model.members.front();
    const double length = model.nodeCoordinatesMetres[1]-model.nodeCoordinatesMetres[0];
    const double forcePerArea = load.mechanicalForceNewtons.back()/member.areaSquareMetres;
    const double heatQuadratic = member.youngModulusPascals*member.expansionPerKelvin*
        member.kelvinVoigtViscosityPascalSeconds/
        (member.densityKilogramsPerCubicMetre*
         member.specificHeatJoulesPerKilogramKelvin*length*length*stepOptions.timeStepSeconds);
    const double linear = (member.youngModulusPascals+
        member.kelvinVoigtViscosityPascalSeconds/stepOptions.timeStepSeconds)/length+
        model.referenceTemperatureKelvin*member.youngModulusPascals*
        member.youngModulusPascals*member.expansionPerKelvin*member.expansionPerKelvin/
        (member.densityKilogramsPerCubicMetre*
         member.specificHeatJoulesPerKilogramKelvin*length);
    const double discriminant = linear*linear-4.0*heatQuadratic*forcePerArea;
    QVERIFY(discriminant > 0.0);
    const double memberExtension = 2.0*forcePerArea/
        (linear+std::sqrt(discriminant));
    const double expectedTemperatureChange =
        (member.kelvinVoigtViscosityPascalSeconds*memberExtension*memberExtension/
             (length*stepOptions.timeStepSeconds)-
         model.referenceTemperatureKelvin*member.youngModulusPascals*
             member.expansionPerKelvin*memberExtension)/
        (member.densityKilogramsPerCubicMetre*
         member.specificHeatJoulesPerKilogramKelvin*length);
    const double memberDissipatedEnergy = member.kelvinVoigtViscosityPascalSeconds*
        member.areaSquareMetres*length*
        std::pow(memberExtension/(length*stepOptions.timeStepSeconds),2)*
        stepOptions.timeStepSeconds;

    QVERIFY(near(result.nextState.axialDisplacementMetres[1],memberExtension));
    QVERIFY(near(result.nextState.axialDisplacementMetres[2],2.0*memberExtension));
    for (double temperature : result.nextState.temperatureKelvin)
    {
        QVERIFY(near(temperature,300.0+expectedTemperatureChange));
    }
    for (double force : result.memberAxialForceNewtons)
    {
        QVERIFY(near(force,load.mechanicalForceNewtons.back()));
    }
    for (double energy : result.memberDissipatedEnergyJoules)
    {
        QVERIFY(near(energy,memberDissipatedEnergy));
    }
    QVERIFY(near(result.viscousDissipationEnergyJoules,2.0*memberDissipatedEnergy));
    QVERIFY(near(result.dissipationHeatEnergyJoules,2.0*memberDissipatedEnergy));
    QVERIFY(near(result.thermalEnergyChangeJoules,
                 member.densityKilogramsPerCubicMetre*
                     member.specificHeatJoulesPerKilogramKelvin*
                     member.areaSquareMetres*2.0*length*expectedTemperatureChange));
    QVERIFY(std::abs(result.thermalEnergyBalanceResidualJoules) <= 1.0e-9);
    QVERIFY(result.newtonIterations > 0);
    QVERIFY(result.newtonIterations <= 8);
}

void TestStructuralThermoMonolithic1D::overflowingNewtonTrialBacktracksToFiniteCandidate()
{
    const auto model = makeModel(1.0,1.0,1000.0,0.0,1.0,10.0,0.0,1.0);
    const auto accepted = uniformState();
    auto load = conditions();
    load.displacementPrescribed[0] = true;
    load.mechanicalForceNewtons[1] = 1.5e157;
    auto bounded = options();
    bounded.maximumNewtonIterations = 1;
    bounded.maximumLineSearchSteps = 2;
    bounded.mechanicalAbsoluteResidualToleranceNewtons = 1.0e157;

    const auto result = RStructuralThermoMonolithic1D::step(
        model,accepted,load,bounded);
    const double fullNewtonDisplacement = load.mechanicalForceNewtons[1]/
        (model.members[0].areaSquareMetres*(model.members[0].youngModulusPascals+
         model.members[0].kelvinVoigtViscosityPascalSeconds/bounded.timeStepSeconds));
    QVERIFY(near(result.nextState.axialDisplacementMetres[1]/fullNewtonDisplacement,0.5));
    QCOMPARE(result.newtonIterations,1U);
    QVERIFY(std::isfinite(result.memberDissipatedEnergyJoules[0]));
    QVERIFY(result.mechanicalResidualNormNewtons <=
            bounded.mechanicalAbsoluteResidualToleranceNewtons);
}

void TestStructuralThermoMonolithic1D::failedNewtonPublishesNoCandidate()
{
    const auto model = makeModel(2.0,1.0,1000.0,0.0,1.0,10.0,0.0,4.0);
    const auto accepted = uniformState();
    const auto before = accepted;
    auto load = conditions();
    load.displacementPrescribed[0] = true;
    load.mechanicalForceNewtons[1] = 100.0;
    auto limited = options(2.0);
    limited.heatKelvinVoigtDissipation = true;
    limited.maximumNewtonIterations = 1;

    QVERIFY_EXCEPTION_THROWN(RStructuralThermoMonolithic1D::step(
        model,accepted,load,limited),std::runtime_error);
    QCOMPARE(accepted.temperatureKelvin.size(),before.temperatureKelvin.size());
    QCOMPARE(accepted.axialDisplacementMetres.size(),before.axialDisplacementMetres.size());
    for (std::size_t i=0;i<accepted.temperatureKelvin.size();i++)
    {
        QCOMPARE(accepted.temperatureKelvin[i],before.temperatureKelvin[i]);
        QCOMPARE(accepted.axialDisplacementMetres[i],before.axialDisplacementMetres[i]);
    }
}

void TestStructuralThermoMonolithic1D::unattainableCoolingNeverPublishesNonpositiveTemperature()
{
    const auto model = makeModel();
    const auto accepted = uniformState();
    const auto before = accepted;
    auto load = conditions();
    load.displacementPrescribed[0] = true;
    // The linear adiabatic solution would cool below 0 K. The bounded Newton
    // line search must fail without exposing that nonphysical candidate.
    load.mechanicalForceNewtons[1] = 10000.0;

    QVERIFY_EXCEPTION_THROWN(RStructuralThermoMonolithic1D::step(
        model,accepted,load,options()),std::runtime_error);
    QCOMPARE(accepted.temperatureKelvin,before.temperatureKelvin);
    QCOMPARE(accepted.axialDisplacementMetres,before.axialDisplacementMetres);
}

QTEST_APPLESS_MAIN(TestStructuralThermoMonolithic1D)
#include "tst_structural_thermo_monolithic_1d.moc"

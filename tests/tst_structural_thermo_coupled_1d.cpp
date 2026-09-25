#include <QtTest>

#include <cmath>
#include <stdexcept>

#include "rstructural_thermo_coupled_1d.h"

namespace
{
RStructuralThermoCoupled1DMember rod(double heatCapacity,
                                    double youngModulus,
                                    std::function<double(double)> expansion,
                                    double viscosity = 0.0)
{
    RStructuralThermoCoupled1DMember member;
    member.elementId = "rod-0";
    member.node0Id = "left";
    member.node1Id = "right";
    member.area = 1.0;
    member.length = 1.0;
    member.conductivity = 1.0;
    member.density = 1.0;
    member.specificHeat = heatCapacity;
    member.internalViscosity = viscosity;
    member.youngModulusAtKelvin = [youngModulus](double) { return youngModulus; };
    member.expansionAtKelvin = std::move(expansion);
    return member;
}

RStructuralThermoCoupled1DModel modelWith(const RStructuralThermoCoupled1DMember &member)
{
    RStructuralThermoCoupled1DModel model;
    model.nodeIds = {"left","right"};
    model.members = {member};
    model.referenceTemperatureKelvin = 300.0;
    return model;
}

std::vector<RStructuralThermoCoupled1DNodeState> uniformState(double temperature)
{
    return {{"left",temperature,0.0},{"right",temperature,0.0}};
}

RStructuralThermoCoupled1DOptions options(double dt = 1.0)
{
    RStructuralThermoCoupled1DOptions result;
    result.dt = dt;
    return result;
}

bool near(double actual,double expected,double tolerance = 1.0e-10)
{
    return std::abs(actual-expected) <= tolerance*std::max(1.0,std::abs(expected));
}
}

class TestStructuralThermoCoupled1D : public QObject
{
    Q_OBJECT

private slots:
    void freeHeatedRodIntegratesTemperatureDependentExpansion();
    void fixedHeatedRodReportsThermalSupportReactions();
    void boundedHeaterConservesTenJoules();
    void temperatureDependentModulusChangesMechanicalResponse();
    void twoWayModeFeedsOnlyPositiveViscousWorkAsHeat();
    void partitionedTwoWayMatchesIndependentSimultaneousRoot();
    void heaterBoundFailurePublishesNoCandidate();
};

void TestStructuralThermoCoupled1D::freeHeatedRodIntegratesTemperatureDependentExpansion()
{
    auto member = rod(1.0,1000.0,[](double temperature) {
        return 1.0e-5+2.0e-7*(temperature-300.0);
    });
    const auto model = modelWith(member);
    const auto accepted = uniformState(300.0);
    RStructuralThermoCoupled1DConditions conditions;
    conditions.prescribedDisplacements = {{"left",0.0}}; // remove rigid translation only
    conditions.heaters = {{"left",5.0,5.0},{"right",5.0,5.0}};

    const auto result = RStructuralThermoCoupled1D::step(model,accepted,conditions,options());
    QVERIFY(near(result.nextState[0].temperatureKelvin,310.0));
    QVERIFY(near(result.nextState[1].temperatureKelvin,310.0));
    const double integratedStrain = 1.0e-5*10.0+0.5*2.0e-7*10.0*10.0;
    QVERIFY(near(result.nextState[1].axialDisplacement,integratedStrain));
    QVERIFY(std::abs(result.nextState[1].axialDisplacement-1.2e-4) > 1.0e-6);
    QVERIFY(std::abs(result.mechanicalResidualNormNewtons) < 1.0e-10);
    QVERIFY(std::abs(result.thermalEnergyBalanceResidualJoules) < 1.0e-10);
}

void TestStructuralThermoCoupled1D::fixedHeatedRodReportsThermalSupportReactions()
{
    const auto model = modelWith(rod(1.0,1000.0,[](double temperature) {
        return 1.0e-5+2.0e-7*(temperature-300.0);
    }));
    const auto accepted = uniformState(300.0);
    RStructuralThermoCoupled1DConditions conditions;
    conditions.prescribedDisplacements = {{"left",0.0},{"right",0.0}};
    conditions.heaters = {{"left",5.0,5.0},{"right",5.0,5.0}};

    const auto result = RStructuralThermoCoupled1D::step(model,accepted,conditions,options());
    const double expectedReaction = 1000.0*1.1e-4;
    QVERIFY(near(result.memberAxialForcesNewtons[0],-expectedReaction));
    QVERIFY(near(result.nodes[0].supportReactionNewtons,expectedReaction));
    QVERIFY(near(result.nodes[1].supportReactionNewtons,-expectedReaction));
    QCOMPARE(result.mechanicalResidualNormNewtons,0.0);
}

void TestStructuralThermoCoupled1D::boundedHeaterConservesTenJoules()
{
    const auto model = modelWith(rod(10.0,100.0,[](double) { return 0.0; }));
    const auto accepted = uniformState(300.0);
    RStructuralThermoCoupled1DConditions conditions;
    conditions.prescribedDisplacements = {{"left",0.0},{"right",0.0}};
    conditions.heaters = {{"left",10.0,10.0}};

    const auto result = RStructuralThermoCoupled1D::step(model,accepted,conditions,options());
    QVERIFY(near(result.heaterEnergyJoules,10.0));
    QVERIFY(near(result.thermalEnergyChangeJoules,10.0));
    QVERIFY(near(result.nextState[0].temperatureKelvin,301.7142857142857));
    QVERIFY(near(result.nextState[1].temperatureKelvin,300.2857142857143));
    QVERIFY(std::abs(result.thermalEnergyBalanceResidualJoules) < 1.0e-9);
}

void TestStructuralThermoCoupled1D::temperatureDependentModulusChangesMechanicalResponse()
{
    auto member = rod(1.0,100.0,[](double) { return 0.0; });
    member.youngModulusAtKelvin = [](double temperature) {
        return 100.0+2.0*(temperature-300.0);
    };
    const auto model = modelWith(member);
    const auto accepted = uniformState(300.0);
    RStructuralThermoCoupled1DConditions conditions;
    conditions.prescribedDisplacements = {{"left",0.0}};
    conditions.mechanicalForces = {{"right",14.0}};
    conditions.heaters = {{"left",10.0,10.0},{"right",10.0,10.0}};

    const auto result = RStructuralThermoCoupled1D::step(model,accepted,conditions,options());
    QVERIFY(near(result.nextState[0].temperatureKelvin,320.0));
    QVERIFY(near(result.nextState[1].temperatureKelvin,320.0));
    QVERIFY(near(result.nextState[1].axialDisplacement,0.1));
    QVERIFY(near(result.memberAxialForcesNewtons[0],14.0));
    QVERIFY(result.mechanicalResidualNormNewtons < 1.0e-9);
}

void TestStructuralThermoCoupled1D::twoWayModeFeedsOnlyPositiveViscousWorkAsHeat()
{
    const auto model = modelWith(rod(5.0,100.0,[](double) { return 0.0; },2.0));
    const auto accepted = uniformState(300.0);
    RStructuralThermoCoupled1DConditions conditions;
    conditions.prescribedDisplacements = {{"left",0.0},{"right",0.1}};
    auto stepOptions = options();
    stepOptions.twoWayDissipationHeat = true;

    const auto result = RStructuralThermoCoupled1D::step(model,accepted,conditions,stepOptions);
    QVERIFY(near(result.memberDissipatedEnergyJoules[0],0.02));
    QVERIFY(near(result.viscousDissipationEnergyJoules,0.02));
    QVERIFY(near(result.dissipationHeatEnergyJoules,0.02));
    QVERIFY(near(result.thermalEnergyChangeJoules,0.02));
    QVERIFY(near(result.nextState[0].temperatureKelvin,300.004));
    QVERIFY(near(result.nextState[1].temperatureKelvin,300.004));
    QVERIFY(result.memberDissipatedEnergyJoules[0] >= 0.0);
    QVERIFY(result.couplingPowerResidualWatts <= 1.0e-10);
    QVERIFY(std::abs(result.thermalEnergyBalanceResidualJoules) < 1.0e-10);
    QVERIFY(near(result.elasticStrainEnergyChangeJoules,0.5));
    QVERIFY(near(result.mechanicalExternalWorkJoules,0.0));
    QVERIFY(near(result.mechanicalSupportWorkJoules,0.52));
    QVERIFY(std::abs(result.mechanicalWorkBalanceResidualJoules) < 1.0e-10);
}

void TestStructuralThermoCoupled1D::partitionedTwoWayMatchesIndependentSimultaneousRoot()
{
    auto member = rod(5.0,100.0,[](double) { return 0.0; },20.0);
    member.youngModulusAtKelvin = [](double kelvin) {
        return 100.0+2.0*(kelvin-300.0);
    };
    const auto model = modelWith(member);
    const auto accepted = uniformState(300.0);
    RStructuralThermoCoupled1DConditions conditions;
    conditions.prescribedDisplacements = {{"left",0.0}};
    conditions.mechanicalForces = {{"right",20.0}};
    auto stepOptions = options();
    stepOptions.twoWayDissipationHeat = true;

    // Independently solve the two simultaneously coupled scalar equations:
    // (E(T)+eta/dt)u=F, C(T-Tn)/dt=eta*(u/dt)^2.
    const auto simultaneousResidual = [](double kelvin) {
        const double modulus = 100.0+2.0*(kelvin-300.0);
        const double displacement = 20.0/(modulus+20.0);
        return 5.0*(kelvin-300.0)-20.0*displacement*displacement;
    };
    double lower = 300.0;
    double upper = 301.0;
    QVERIFY(simultaneousResidual(lower) < 0.0);
    QVERIFY(simultaneousResidual(upper) > 0.0);
    for (unsigned int iteration=0;iteration<80;iteration++)
    {
        const double middle = 0.5*(lower+upper);
        if (simultaneousResidual(middle) > 0.0) upper = middle;
        else lower = middle;
    }
    const double simultaneousTemperature = 0.5*(lower+upper);
    const double simultaneousDisplacement =
        20.0/(100.0+2.0*(simultaneousTemperature-300.0)+20.0);

    const auto partitioned = RStructuralThermoCoupled1D::step(
        model,accepted,conditions,stepOptions);
    QVERIFY(partitioned.couplingIterations > 1);
    QVERIFY(near(partitioned.nextState[0].temperatureKelvin,simultaneousTemperature,1.0e-9));
    QVERIFY(near(partitioned.nextState[1].temperatureKelvin,simultaneousTemperature,1.0e-9));
    QVERIFY(near(partitioned.nextState[1].axialDisplacement,simultaneousDisplacement,1.0e-9));
    QVERIFY(near(partitioned.dissipationHeatEnergyJoules,
                 20.0*simultaneousDisplacement*simultaneousDisplacement,1.0e-9));
    QVERIFY(std::abs(partitioned.thermalEnergyBalanceResidualJoules) < 1.0e-9);
}

void TestStructuralThermoCoupled1D::heaterBoundFailurePublishesNoCandidate()
{
    const auto model = modelWith(rod(1.0,100.0,[](double) { return 0.0; }));
    const auto accepted = uniformState(300.0);
    RStructuralThermoCoupled1DConditions conditions;
    conditions.prescribedDisplacements = {{"left",0.0},{"right",0.0}};
    conditions.heaters = {{"left",10.01,10.0}};
    const auto before = accepted;

    QVERIFY_EXCEPTION_THROWN(RStructuralThermoCoupled1D::step(model,accepted,conditions,options()),
                             std::invalid_argument);
    QCOMPARE(accepted.size(),before.size());
    for (std::size_t i=0;i<accepted.size();i++)
    {
        QCOMPARE(accepted[i].nodeId,before[i].nodeId);
        QCOMPARE(accepted[i].temperatureKelvin,before[i].temperatureKelvin);
        QCOMPARE(accepted[i].axialDisplacement,before[i].axialDisplacement);
    }
}

QTEST_APPLESS_MAIN(TestStructuralThermoCoupled1D)
#include "tst_structural_thermo_coupled_1d.moc"

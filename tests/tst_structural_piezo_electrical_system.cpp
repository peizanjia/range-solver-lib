#include <QtTest>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "rstructural_piezo_electrical_system.h"

namespace
{
using Role = RStructuralPiezoElectricalRole;
using Boundary = RStructuralPiezoElectricalBoundary;

RStructuralPiezoElementOperators operators()
{
    RStructuralPiezoElementOperators result;
    result.dofs = 2;
    result.portIds = {"actuator","open","short","amplifier"};
    result.mass = {2.0,0.25,0.25,3.0};
    result.stiffness = {10.0,2.0,2.0,20.0};
    result.coupling = {2e-3,3e-3,-1e-3,4e-3,-1e-3,2e-3,2e-3,-3e-3};
    result.capacitance = {2e-6,0,0,0,0,5e-6,0,0,0,0,7e-6,0,0,0,0,3e-6};
    result.totalMass = 5.0;
    return result;
}

std::vector<Boundary> boundaries()
{
    // Intentionally different from physical operator order: exact IDs own roles.
    return {{"amplifier",Role::IdealChargeAmplifier,0,2e-6,3e-6},
            {"short",Role::ShortCircuitChargeSensor,0,0,0},
            {"open",Role::FixedChargeOpenCircuit,1e-6,0,0},
            {"actuator",Role::VoltageActuator,0,0,0}};
}

bool close(double actual,double expected,double relativeTolerance=2e-12)
{
    if (expected == 0.0) return actual == 0.0;
    return std::abs(actual-expected) <= relativeTolerance*std::abs(expected);
}
}

class TestStructuralPiezoElectricalSystem : public QObject
{
    Q_OBJECT
private slots:
    void mixedPortsAnalyticalReference();
    void openSubsetStiffeningIsPositiveSemidefinite();
    void multipleActuatorsAndOpenSensors();
    void finiteDifferenceReciprocityAndActuatorDirectTerm();
    void fixedChargeBiasAndReferenceIncrements();
    void shortCircuitAndAmplifierReadout();
    void scaleAwareCapacitance();
    void representableFinalValuesWithExtremeIntermediates();
    void malformedBoundariesAndOperators();
    void failedEvaluationPreservesWorkspace();
    void nonrepresentableCondensationAndReadout();
};

void TestStructuralPiezoElectricalSystem::mixedPortsAnalyticalReference()
{
    const RStructuralPiezoElectricalSystem system(operators(),boundaries(),{0.002,0.001},{3.0});
    const auto &condensed = system.condensed();
    QVERIFY(condensed.mass == operators().mass);
    QVERIFY(condensed.actuatorPortIds == std::vector<QString>{"actuator"});
    QVERIFY(close(condensed.stiffness[0],11.8));
    QVERIFY(close(condensed.stiffness[1],3.2));
    QVERIFY(close(condensed.stiffness[2],3.2));
    QVERIFY(close(condensed.stiffness[3],20.8));
    QVERIFY(close(condensed.coupling[0],2e-3));
    QVERIFY(close(condensed.coupling[1],-1e-3));
    QVERIFY(close(condensed.fixedChargeForce[0],6e-4));
    QVERIFY(close(condensed.fixedChargeForce[1],4e-4));
    const auto response = system.evaluate({0.01,-0.02},{12});
    QVERIFY(close(response.internalForce[0],0.0786));
    QVERIFY(close(response.internalForce[1],-0.3956));
    QVERIFY(close(response.deltaInternalForce[0],0.0452));
    QVERIFY(close(response.deltaInternalForce[1],-0.4202));
    QVERIFY(close(response.ports[0].voltage,12));
    QVERIFY(close(response.ports[0].charge,-16e-6));
    QVERIFY(close(response.ports[1].voltage,-1.8));
    QVERIFY(close(response.ports[1].charge,1e-6));
    QVERIFY(close(response.ports[2].charge,50e-6));
    QVERIFY(close(response.ports[3].charge,-100e-6));
    QVERIFY(close(response.ports[3].amplifierVoltage,51.5));
}

void TestStructuralPiezoElectricalSystem::openSubsetStiffeningIsPositiveSemidefinite()
{
    const auto source = operators();
    const RStructuralPiezoElectricalSystem system(source,boundaries(),{0,0},{0});
    const auto &k = system.condensed().stiffness;
    // Independent analytic rank-one difference: [1.8,1.2;1.2,0.8].
    QVERIFY(close(k[0]-source.stiffness[0],1.8));
    QVERIFY(close(k[1]-source.stiffness[1],1.2));
    QVERIFY(close(k[3]-source.stiffness[3],0.8));
    for (const auto &q : std::vector<std::vector<double>>{{1,0},{0,1},{1,-1},{2,-3},{-7,2}})
    {
        const double energy = q[0]*(1.8*q[0]+1.2*q[1])+q[1]*(1.2*q[0]+0.8*q[1]);
        QVERIFY(energy >= -1e-14);
        const double expected = std::pow(3e-3*q[0]+2e-3*q[1],2)/5e-6;
        QVERIFY(std::abs(energy-expected) < 2e-14*std::max(1.0,expected));
    }
    auto shortSettings = boundaries();
    shortSettings[2] = {"open",Role::ShortCircuitChargeSensor,0,0,0};
    const RStructuralPiezoElectricalSystem shortSystem(source,shortSettings,{0,0},{0});
    QVERIFY(shortSystem.condensed().stiffness == source.stiffness);
}

void TestStructuralPiezoElectricalSystem::multipleActuatorsAndOpenSensors()
{
    auto settings = boundaries();
    settings[0] = {"amplifier",Role::FixedChargeOpenCircuit,-2e-6,0,0};
    settings[1] = {"short",Role::VoltageActuator,0,0,0};
    const RStructuralPiezoElectricalSystem system(operators(),settings,{0,0},{0,0});
    const auto &condensed = system.condensed();
    const std::vector<QString> expectedIds{"actuator","short"};
    QVERIFY(condensed.actuatorPortIds == expectedIds);
    QVERIFY(close(condensed.stiffness[0],11.8+16.0/3.0));
    QVERIFY(close(condensed.stiffness[1],-0.8));
    QVERIFY(close(condensed.stiffness[3],23.8));
    const std::vector<double> expectedCoupling{2e-3,-1e-3,-1e-3,2e-3};
    QVERIFY(condensed.coupling == expectedCoupling);
    const auto response = system.evaluate({0.01,-0.02},{12,-4});
    QVERIFY(close(response.ports[2].charge,22e-6));
    QVERIFY(close(response.ports[3].voltage,98.0/3.0));
    QVERIFY(close(response.ports[3].charge,-2e-6));
}

void TestStructuralPiezoElectricalSystem::finiteDifferenceReciprocityAndActuatorDirectTerm()
{
    const RStructuralPiezoElectricalSystem system(operators(),boundaries(),{0,0},{0});
    const double h = 1e-4;
    const auto plusV = system.evaluate({0.01,-0.02},{12+h});
    const auto minusV = system.evaluate({0.01,-0.02},{12-h});
    const double g[] = {2e-3,-1e-3};
    for (size_t i=0;i<2;++i)
    {
        QVERIFY(close((plusV.internalForce[i]-minusV.internalForce[i])/(2*h),g[i],1e-8));
        std::vector<double> qp{0.01,-0.02}, qm = qp;
        qp[i] += h; qm[i] -= h;
        const auto plusQ = system.evaluate(qp,{12});
        const auto minusQ = system.evaluate(qm,{12});
        QVERIFY(close((plusQ.ports[0].charge-minusQ.ports[0].charge)/(2*h),-g[i],1e-8));
    }
    QVERIFY(close((plusV.ports[0].charge-minusV.ports[0].charge)/(2*h),2e-6,1e-8));
    QVERIFY(close(system.evaluate({0,0},{12}).ports[0].charge,24e-6));
}

void TestStructuralPiezoElectricalSystem::fixedChargeBiasAndReferenceIncrements()
{
    auto q0 = std::vector<double>{0.002,0.001};
    auto va0 = std::vector<double>{3};
    auto settings = boundaries();
    auto source = operators();
    const RStructuralPiezoElectricalSystem system(source,settings,q0,va0);
    // Input objects may change later; reference/operator snapshots are frozen.
    q0[0]=999; va0[0]=999; settings[2].fixedCharge=999; source.coupling[0]=999;
    QVERIFY(close(system.referenceDisplacement()[0],0.002));
    QVERIFY(close(system.referenceActuatorVoltages()[0],3));
    QVERIFY(close(system.referenceResponse().ports[1].voltage,1.8));
    const auto response = system.evaluate({0.01,-0.02},{12});
    QVERIFY(close(response.ports[1].deltaVoltage,-3.6));
    QCOMPARE(response.ports[1].deltaCharge,0.0);
    QVERIFY(close(response.ports[0].deltaCharge,-19e-6));
    QVERIFY(close(response.ports[3].deltaAmplifierVoltage,47.5));
    const auto atReference = system.evaluate({0.002,0.001},{3});
    for (const auto &port : atReference.ports)
    {
        QCOMPARE(port.deltaVoltage,0.0);
        QCOMPARE(port.deltaCharge,0.0);
        QCOMPARE(port.deltaAmplifierVoltage,0.0);
    }
    const RStructuralPiezoElectricalSystem zero(operators(),boundaries(),{0,0},{0});
    const auto origin = zero.evaluate({0,0},{0});
    QVERIFY(close(origin.ports[1].voltage,0.2));
    QVERIFY(close(origin.internalForce[0],6e-4));
    QVERIFY(close(origin.internalForce[1],4e-4));
    // Large fixed/reference charge biases can swallow differences of absolutes;
    // increments still follow the independent linear increment equations.
    settings = boundaries();
    settings[0].referenceCharge = 1e10;
    settings[2].fixedCharge = 1e10;
    const RStructuralPiezoElectricalSystem biased(operators(),settings,{0,0},{0});
    const auto signal = biased.evaluate({1e-6,0},{0});
    QVERIFY(close(signal.ports[1].deltaVoltage,6e-4));
    QVERIFY(close(signal.ports[3].deltaAmplifierVoltage,2e-3));
    QVERIFY(close(signal.deltaInternalForce[0],11.8e-6));
}

void TestStructuralPiezoElectricalSystem::shortCircuitAndAmplifierReadout()
{
    const RStructuralPiezoElectricalSystem system(operators(),boundaries(),{0,0},{0});
    const auto response = system.evaluate({0.01,-0.02},{0});
    QCOMPARE(response.ports[2].voltage,0.0);
    QCOMPARE(response.ports[3].voltage,0.0);
    QVERIFY(close(response.ports[2].charge,50e-6));
    QVERIFY(close(response.ports[3].amplifierVoltage,51.5));
    QVERIFY(close(response.ports[3].deltaAmplifierVoltage,50));
    // Positive electrode charge produces negative amplifier voltage relative to Qref.
    const auto reverse = system.evaluate({-0.01,0.02},{0});
    QVERIFY(close(reverse.ports[3].amplifierVoltage,-48.5));
}

void TestStructuralPiezoElectricalSystem::scaleAwareCapacitance()
{
    auto source = operators();
    source.coupling.assign(8,0.0);
    source.coupling[1] = 1e-20;
    source.capacitance[5] = 1e-30;
    auto settings = boundaries();
    settings[2].fixedCharge = 0;
    const RStructuralPiezoElectricalSystem system(source,settings,{0,0},{0});
    const auto response = system.evaluate({1,0},{0});
    QVERIFY(close(response.ports[1].voltage,1e10));
    QVERIFY(close(system.condensed().stiffness[0]-10,1e-10,1e-5));
    // Widely separated dielectric and coupling magnitudes, including on MSVC.
    source.coupling[1] = 1e-150;
    source.capacitance[5] = 1e-300;
    const RStructuralPiezoElectricalSystem scaled(source,settings,{0,0},{0});
    QVERIFY(close(scaled.condensed().stiffness[0],11));
    QVERIFY(close(scaled.evaluate({1,0},{0}).ports[1].voltage,1e150));
    // No electrical ports is a valid purely mechanical operator.
    source.portIds.clear(); source.coupling.clear(); source.capacitance.clear();
    const RStructuralPiezoElectricalSystem passive(source,{}, {0,0},{});
    QVERIFY(passive.evaluate({1,2},{}).ports.empty());
    QVERIFY(close(passive.evaluate({1,2},{}).internalForce[1],42));
}

void TestStructuralPiezoElectricalSystem::representableFinalValuesWithExtremeIntermediates()
{
    RStructuralPiezoElementOperators source;
    source.dofs = 1;
    source.mass = {1}; source.totalMass = 1;
    source.stiffness = {1};
    source.portIds = {"open"}; source.coupling = {1e-200}; source.capacitance = {1e-300};
    const std::vector<Boundary> settings{{"open",Role::FixedChargeOpenCircuit,0,0,0}};
    const RStructuralPiezoElectricalSystem system(source,settings,{0},{});
    const auto response = system.evaluate({1e-130},{});
    // G*q = 1e-330 is below double range, but it is NOT a returned quantity.
    QVERIFY(close(response.ports[0].voltage,1e-30));
    QVERIFY(close(response.ports[0].deltaVoltage,1e-30));
    QCOMPARE(response.ports[0].charge,0.0);
    QVERIFY(close(response.internalForce[0],1e-130));
    QVERIFY(close(response.deltaInternalForce[0],1e-130));
    // G*G = 1e-400 is below double range; its stiffness correction is 1e-100.
    source.stiffness = {1e-100};
    const RStructuralPiezoElectricalSystem soft(source,settings,{0},{});
    QVERIFY(close(soft.condensed().stiffness[0],2e-100));
    const auto softResponse = soft.evaluate({1e-130},{});
    QVERIFY(close(softResponse.internalForce[0],2e-230));
    QVERIFY(close(softResponse.deltaInternalForce[0],2e-230));
    // With fixed Q, absolute bias is retained but the tiny deformation increment
    // is calculated independently of subtraction of large absolute voltages.
    auto biasedSettings = settings; biasedSettings[0].fixedCharge = 1e-200;
    const RStructuralPiezoElectricalSystem biased(source,biasedSettings,{0},{});
    const auto biasedResponse = biased.evaluate({1e-130},{});
    QVERIFY(close(biasedResponse.ports[0].voltage,1e100));
    QVERIFY(close(biasedResponse.ports[0].deltaVoltage,1e-30));
    QCOMPARE(biasedResponse.ports[0].charge,1e-200);

    source.dofs = 3;
    source.mass = {1,0,0,0,1,0,0,0,1}; source.totalMass=3;
    source.portIds.clear(); source.coupling.clear(); source.capacitance.clear();
    source.stiffness.assign(9,1e200);
    const RStructuralPiezoElectricalSystem cancel(source,{}, {0,0,0},{});
    // Individual force products are +/-1e400. Their cancellation must happen
    // before conversion to double, preserving the remaining 1e-100 force.
    const auto cancelled = cancel.evaluate({1e200,1e-300,-1e200},{});
    for (size_t i=0;i<3;++i)
    {
        QVERIFY(close(cancelled.internalForce[i],1e-100));
        QVERIFY(close(cancelled.deltaInternalForce[i],1e-100));
    }
    source.stiffness.assign(9,1e-200);
    const RStructuralPiezoElectricalSystem subnormal(source,{}, {0,0,0},{});
    // Each term is about 2e-324 and rounds to zero separately; their sum is
    // about 6e-324 and rounds to the first representable subnormal value.
    const auto tiny = subnormal.evaluate({2e-124,2e-124,2e-124},{});
    for (size_t i=0;i<3;++i)
    {
        QCOMPARE(tiny.internalForce[i],std::numeric_limits<double>::denorm_min());
        QCOMPARE(tiny.deltaInternalForce[i],std::numeric_limits<double>::denorm_min());
    }
    auto workspace = tiny;
    QVERIFY_EXCEPTION_THROWN(subnormal.evaluateInto({2e-124,0,0},{},workspace),std::invalid_argument);
    QVERIFY(workspace.internalForce == tiny.internalForce);

    source.dofs=1; source.mass={1}; source.totalMass=1; source.stiffness={1.5};
    const double adjacent = std::nextafter(2.0,0.0);
    const RStructuralPiezoElectricalSystem adjacentReference(source,{}, {adjacent},{});
    const auto adjacentResponse = adjacentReference.evaluate({2},{});
    // Difference must precede multiplication: subtracting rounded K*q products
    // gives 2^-51 instead of the independent exact reference 1.5*2^-52.
    QCOMPARE(adjacentResponse.deltaInternalForce[0],1.5*std::ldexp(1.0,-52));
    source.stiffness={1e-200};
    const RStructuralPiezoElectricalSystem wideReference(source,{}, {-1e308},{});
    const auto wideResponse = wideReference.evaluate({1e308},{});
    QVERIFY(close(wideResponse.internalForce[0],1e108));
    QVERIFY(close(wideResponse.deltaInternalForce[0],2e108));
    // Check the same adjacent-reference semantics for open and short readouts.
    source.stiffness={1}; source.portIds={"sensor"}; source.coupling={1.5}; source.capacitance={1};
    const RStructuralPiezoElectricalSystem adjacentOpen(source,
        {{"sensor",Role::FixedChargeOpenCircuit,0,0,0}}, {adjacent},{});
    QCOMPARE(adjacentOpen.evaluate({2},{}).ports[0].deltaVoltage,1.5*std::ldexp(1.0,-52));
    const RStructuralPiezoElectricalSystem adjacentShort(source,
        {{"sensor",Role::ShortCircuitChargeSensor,0,0,0}}, {adjacent},{});
    QCOMPARE(adjacentShort.evaluate({2},{}).ports[0].deltaCharge,-1.5*std::ldexp(1.0,-52));
}

void TestStructuralPiezoElectricalSystem::malformedBoundariesAndOperators()
{
    const auto source = operators();
    auto settings = boundaries();
    settings.pop_back();
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoElectricalSystem(source,settings,{0,0},{0}),std::invalid_argument);
    settings=boundaries(); settings[0].portId="open";
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoElectricalSystem(source,settings,{0,0},{0}),std::invalid_argument);
    settings=boundaries(); settings[0].portId="amplifier ";
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoElectricalSystem(source,settings,{0,0},{0}),std::invalid_argument);
    settings=boundaries(); settings[0].role=static_cast<Role>(999);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoElectricalSystem(source,settings,{0,0},{0}),std::invalid_argument);
    for (double cf : {0.0,-1.0,std::numeric_limits<double>::infinity()})
    {
        settings=boundaries(); settings[0].feedbackCapacitance=cf;
        QVERIFY_EXCEPTION_THROWN(RStructuralPiezoElectricalSystem(source,settings,{0,0},{0}),std::invalid_argument);
    }
    settings=boundaries(); settings[1].fixedCharge=1;
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoElectricalSystem(source,settings,{0,0},{0}),std::invalid_argument);
    settings=boundaries(); settings[2].referenceCharge=1;
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoElectricalSystem(source,settings,{0,0},{0}),std::invalid_argument);
    auto invalid = source; invalid.capacitance[1]=invalid.capacitance[4]=1e-20;
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoElectricalSystem(invalid,boundaries(),{0,0},{0}),std::invalid_argument);
    invalid=source; invalid.capacitance[5]=0;
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoElectricalSystem(invalid,boundaries(),{0,0},{0}),std::invalid_argument);
    invalid=source; invalid.coupling[0]=std::numeric_limits<double>::quiet_NaN();
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoElectricalSystem(invalid,boundaries(),{0,0},{0}),std::invalid_argument);
    invalid=source; invalid.mass.pop_back();
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoElectricalSystem(invalid,boundaries(),{0,0},{0}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoElectricalSystem(source,boundaries(),{0},{0}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoElectricalSystem(source,boundaries(),{0,0},{}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoElectricalSystem(source,boundaries(),
                             {0,std::numeric_limits<double>::quiet_NaN()},{0}),std::invalid_argument);
    settings=boundaries(); settings[0].referenceCharge=std::numeric_limits<double>::quiet_NaN();
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoElectricalSystem(source,settings,{0,0},{0}),std::invalid_argument);
}

void TestStructuralPiezoElectricalSystem::failedEvaluationPreservesWorkspace()
{
    const RStructuralPiezoElectricalSystem system(operators(),boundaries(),{0,0},{0});
    auto workspace = system.evaluate({0.01,-0.02},{12});
    const auto original = workspace;
    for (const auto &q : std::vector<std::vector<double>>{{0},{0,std::numeric_limits<double>::quiet_NaN()},
                                                       {std::numeric_limits<double>::max(),0}})
    {
        QVERIFY_EXCEPTION_THROWN(system.evaluateInto(q,{12},workspace),std::invalid_argument);
        QVERIFY(workspace.internalForce == original.internalForce);
        QVERIFY(workspace.deltaInternalForce == original.deltaInternalForce);
        QCOMPARE(workspace.ports.size(),original.ports.size());
        for (size_t i=0;i<workspace.ports.size();++i)
        {
            QVERIFY(workspace.ports[i].portId == original.ports[i].portId);
            QCOMPARE(workspace.ports[i].voltage,original.ports[i].voltage);
            QCOMPARE(workspace.ports[i].charge,original.ports[i].charge);
            QCOMPARE(workspace.ports[i].amplifierVoltage,original.ports[i].amplifierVoltage);
            QCOMPARE(workspace.ports[i].deltaVoltage,original.ports[i].deltaVoltage);
            QCOMPARE(workspace.ports[i].deltaCharge,original.ports[i].deltaCharge);
            QCOMPARE(workspace.ports[i].deltaAmplifierVoltage,original.ports[i].deltaAmplifierVoltage);
        }
    }
    QVERIFY_EXCEPTION_THROWN(system.evaluateInto({0,0},{},workspace),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(system.evaluateInto({0,0},{std::numeric_limits<double>::infinity()},workspace),std::invalid_argument);
    QVERIFY(workspace.internalForce == original.internalForce);
}

void TestStructuralPiezoElectricalSystem::nonrepresentableCondensationAndReadout()
{
    auto source = operators();
    source.coupling[1]=1e200;
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoElectricalSystem(source,boundaries(),{0,0},{0}),std::invalid_argument);
    auto settings=boundaries(); settings[0].feedbackCapacitance=1e-310;
    const RStructuralPiezoElectricalSystem amplifier(operators(),settings,{0,0},{0});
    auto workspace = amplifier.evaluate({0,0},{0});
    const auto before = workspace;
    QVERIFY_EXCEPTION_THROWN(amplifier.evaluateInto({1e10,0},{0},workspace),std::invalid_argument);
    QVERIFY(workspace.internalForce == before.internalForce);
    QCOMPARE(workspace.ports[3].amplifierVoltage,before.ports[3].amplifierVoltage);
    settings=boundaries(); settings[0].feedbackCapacitance=1e300; settings[0].referenceCharge=0;
    const RStructuralPiezoElectricalSystem underflow(operators(),settings,{0,0},{0});
    workspace = underflow.evaluate({0,0},{0});
    QVERIFY_EXCEPTION_THROWN(underflow.evaluateInto({1e-100,0},{0},workspace),std::invalid_argument);
    QCOMPARE(workspace.ports[3].amplifierVoltage,0.0);
}

QTEST_APPLESS_MAIN(TestStructuralPiezoElectricalSystem)
#include "tst_structural_piezo_electrical_system.moc"

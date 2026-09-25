#include <QtTest>

#include <cmath>
#include <limits>
#include <stdexcept>

#include "rstructural_piezo_coupled_1d.h"

namespace
{
RStructuralPiezoCoupled1DInput oneRod()
{
    RStructuralPiezoCoupled1DInput input;
    input.coordinates = {0.0,1.0};
    input.members = {{0,1,2.0,120.0,18.0,6.0}};
    input.mechanicalLoads = {0.0,0.0};
    input.prescribedConductorCharges = {0.0,0.0};
    input.displacementPrescribed = {true,false};
    input.prescribedDisplacements = {0.0,0.0};
    input.potentialPrescribed = {true,true};
    input.prescribedPotentials = {0.0,0.0};
    input.mechanicalAbsoluteTolerance = 1.0e-10;
    input.mechanicalRelativeTolerance = 1.0e-11;
    input.electricAbsoluteTolerance = 1.0e-12;
    input.electricRelativeTolerance = 1.0e-11;
    return input;
}

RStructuralPiezoCoupled1DResult solveOpenCircuitRod(double force)
{
    RStructuralPiezoCoupled1DInput input = oneRod();
    input.mechanicalLoads[1] = force;
    input.potentialPrescribed[1] = false;
    return RStructuralPiezoCoupled1D::solve(input);
}

RStructuralPiezoCoupled1DInput twoRodSeries()
{
    RStructuralPiezoCoupled1DInput input;
    input.coordinates = {0.0,1.0,2.0};
    input.members = {{0,1,1.0,50.0,5.0,10.0},
                     {1,2,1.0,50.0,5.0,10.0}};
    input.mechanicalLoads = {0.0,0.0,0.0};
    input.prescribedConductorCharges = {0.0,0.0,0.0};
    input.displacementPrescribed = {true,false,false};
    input.prescribedDisplacements = {0.0,0.0,0.0};
    input.potentialPrescribed = {true,false,true};
    input.prescribedPotentials = {0.0,0.0,1.0};
    input.mechanicalAbsoluteTolerance = 1.0e-10;
    input.mechanicalRelativeTolerance = 1.0e-11;
    input.electricAbsoluteTolerance = 1.0e-12;
    input.electricRelativeTolerance = 1.0e-11;
    return input;
}
}

class TestStructuralPiezoCoupled1D : public QObject
{
    Q_OBJECT

private slots:
    void freeActuatorHasAnalyticExtensionChargeAndEnthalpy();
    void shortCircuitAndOpenCircuitStiffnesses();
    void seriesRodRecoversInternalDisplacementAndPotential();
    void prescribedChargeAndElectrodeReactionBalance();
    void condensedEnergySatisfiesElectromechanicalReciprocity();
    void reversingConnectivityPreservesPhysicalResponseAndMapsTerminalCharge();
    void rejectMissingMechanicalOrElectricalPins();
    void failClosedOnUnderconvergenceAndInvalidInputs();
    void highStiffnessSelfEquilibratedThermalRod();
};

void TestStructuralPiezoCoupled1D::highStiffnessSelfEquilibratedThermalRod()
{
    RStructuralPiezoCoupled1DInput input;
    input.coordinates = {0.0,1.0,2.0};
    constexpr double young = 2.0e11;
    constexpr double area = 0.01;
    constexpr double alpha = 1.2e-5;
    constexpr double referenceKelvin = 293.15;
    const double firstStrain = alpha*(325.0-referenceKelvin);
    const double secondStrain = alpha*(375.0-referenceKelvin);
    input.members = {{0,1,area,young,2.0,4.0,firstStrain},
                     {1,2,area,young,2.0,4.0,secondStrain}};
    input.mechanicalLoads = {0.0,0.0,0.0};
    input.prescribedConductorCharges = {0.0,0.0,0.0};
    input.displacementPrescribed = {true,false,true};
    input.prescribedDisplacements = {0.0,0.0,0.0};
    input.potentialPrescribed = {true,true,true};
    input.prescribedPotentials = {0.0,0.0,0.0};

    RStructuralPiezoCoupled1DResult result;
    try
    {
        result = RStructuralPiezoCoupled1D::solve(input);
    }
    catch (const std::exception &error)
    {
        QFAIL(error.what());
    }
    const double expectedStress = -young*0.5*(firstStrain+secondStrain);
    const double expectedMiddle = 0.5*(firstStrain-secondStrain);
    QVERIFY(std::abs(result.displacements[1]-expectedMiddle) < 1.0e-12);
    for (double force : result.memberAxialForces)
        QVERIFY(std::abs(force-area*expectedStress) < 1.0e-5);
    QVERIFY(std::abs(result.internalForces[1]) < 1.0e-5);
    QVERIFY(std::isfinite(result.mechanicalResidualNorm));
}

void TestStructuralPiezoCoupled1D::freeActuatorHasAnalyticExtensionChargeAndEnthalpy()
{
    constexpr double voltage = 2.0;
    RStructuralPiezoCoupled1DInput input = oneRod();
    input.potentialPrescribed[1] = true;
    input.prescribedPotentials[1] = voltage;
    const auto result = RStructuralPiezoCoupled1D::solve(input);

    const double expectedDisplacement = -18.0*voltage/120.0;
    const double expectedCharge = 2.0*(6.0+18.0*18.0/120.0)*voltage;
    QVERIFY(std::abs(result.displacements[1]-expectedDisplacement) < 1.0e-12);
    QVERIFY(std::abs(result.internalForces[1]) < 1.0e-10);
    QVERIFY(std::abs(result.memberAxialForces[0]) < 1.0e-10);
    QVERIFY(std::abs(result.conductorCharges[1]-expectedCharge) < 1.0e-10);
    QVERIFY(std::abs(result.memberRightCharges[0]-expectedCharge) < 1.0e-10);
    QVERIFY(std::abs(result.enthalpy+0.5*expectedCharge*voltage) < 1.0e-10);
    QVERIFY(result.mechanicalResidualNorm < 1.0e-10);
    QCOMPARE(result.electricResidualNorm,0.0);
}

void TestStructuralPiezoCoupled1D::shortCircuitAndOpenCircuitStiffnesses()
{
    constexpr double force = 15.0;
    RStructuralPiezoCoupled1DInput shortCircuit = oneRod();
    shortCircuit.mechanicalLoads[1] = force;
    const auto shortResult = RStructuralPiezoCoupled1D::solve(shortCircuit);
    const double shortStiffness = 120.0*2.0;
    QVERIFY(std::abs(shortResult.displacements[1]-force/shortStiffness) < 1.0e-12);

    const auto openResult = solveOpenCircuitRod(force);
    const double openStiffness = (120.0+18.0*18.0/6.0)*2.0;
    QVERIFY(std::abs(openResult.displacements[1]-force/openStiffness) < 1.0e-12);
    QVERIFY(std::abs(openResult.potentials[1]-3.0*openResult.displacements[1]) < 1.0e-12);
    QVERIFY(std::abs(openResult.conductorCharges[1]) < 1.0e-11);
    QVERIFY(openStiffness > shortStiffness);
}

void TestStructuralPiezoCoupled1D::seriesRodRecoversInternalDisplacementAndPotential()
{
    const auto result = RStructuralPiezoCoupled1D::solve(twoRodSeries());
    QVERIFY(std::abs(result.displacements[2]+0.1) < 1.0e-12);
    QVERIFY(std::abs(result.displacements[1]+0.05) < 1.0e-12);
    QVERIFY(std::abs(result.potentials[1]-0.5) < 1.0e-12);
    QVERIFY(std::abs(result.conductorCharges[1]) < 1.0e-11);
    QVERIFY(std::abs(result.memberAxialForces[0]) < 1.0e-10);
    QVERIFY(std::abs(result.memberAxialForces[1]) < 1.0e-10);
}

void TestStructuralPiezoCoupled1D::prescribedChargeAndElectrodeReactionBalance()
{
    RStructuralPiezoCoupled1DInput input = oneRod();
    input.potentialPrescribed[1] = false;
    input.prescribedConductorCharges[1] = 3.0;
    const auto result = RStructuralPiezoCoupled1D::solve(input);

    const double effectiveCapacitance = 2.0*(6.0+18.0*18.0/120.0);
    const double expectedPotential = 3.0/effectiveCapacitance;
    const double expectedDisplacement = -18.0*expectedPotential/120.0;
    QVERIFY(std::abs(result.potentials[1]-expectedPotential) < 1.0e-12);
    QVERIFY(std::abs(result.displacements[1]-expectedDisplacement) < 1.0e-12);
    QVERIFY(std::abs(result.conductorCharges[1]-3.0) < 1.0e-10);
    QVERIFY(std::abs(result.conductorCharges[0]+3.0) < 1.0e-10);
    QVERIFY(std::abs(result.conductorCharges[0]+result.conductorCharges[1]) < 1.0e-12);
}

void TestStructuralPiezoCoupled1D::condensedEnergySatisfiesElectromechanicalReciprocity()
{
    const auto solveWithBoundaryValues = [](double tipDisplacement,double tipPotential)
    {
        RStructuralPiezoCoupled1DInput input = twoRodSeries();
        input.displacementPrescribed[2] = true;
        input.prescribedDisplacements[2] = tipDisplacement;
        input.prescribedPotentials[2] = tipPotential;
        return RStructuralPiezoCoupled1D::solve(input);
    };

    constexpr double displacement = 0.04;
    constexpr double potential = 0.7;
    constexpr double dStep = 1.0e-6;
    constexpr double vStep = 1.0e-5;
    const auto voltagePlus = solveWithBoundaryValues(displacement,potential+vStep);
    const auto voltageMinus = solveWithBoundaryValues(displacement,potential-vStep);
    const double forceByVoltage =
        (voltagePlus.internalForces[2]-voltageMinus.internalForces[2])/(2.0*vStep);

    const auto displacementPlus = solveWithBoundaryValues(displacement+dStep,potential);
    const auto displacementMinus = solveWithBoundaryValues(displacement-dStep,potential);
    const double chargeByDisplacement =
        (displacementPlus.conductorCharges[2]-displacementMinus.conductorCharges[2])/(2.0*dStep);

    // The electric-enthalpy Hessian is symmetric: dF/dV = -d(q_right)/du
    // under the kernel's physical right-electrode charge convention.
    QVERIFY(std::abs(forceByVoltage+chargeByDisplacement) < 2.0e-7);
}

void TestStructuralPiezoCoupled1D::reversingConnectivityPreservesPhysicalResponseAndMapsTerminalCharge()
{
    RStructuralPiezoCoupled1DInput forward = oneRod();
    forward.prescribedPotentials[1] = 1.25;
    const auto forwardResult = RStructuralPiezoCoupled1D::solve(forward);

    RStructuralPiezoCoupled1DInput reversed = forward;
    reversed.members[0].node0 = 1;
    reversed.members[0].node1 = 0;
    reversed.members[0].coupling = -reversed.members[0].coupling;
    const auto reversedResult = RStructuralPiezoCoupled1D::solve(reversed);

    for (std::size_t node=0;node<2;node++)
    {
        QVERIFY(std::abs(forwardResult.displacements[node]-reversedResult.displacements[node]) < 1.0e-12);
        QVERIFY(std::abs(forwardResult.internalForces[node]-reversedResult.internalForces[node]) < 1.0e-10);
        QVERIFY(std::abs(forwardResult.conductorCharges[node]-reversedResult.conductorCharges[node]) < 1.0e-10);
    }
    QVERIFY(std::abs(forwardResult.memberRightCharges[0]+reversedResult.memberRightCharges[0]) < 1.0e-10);
}

void TestStructuralPiezoCoupled1D::rejectMissingMechanicalOrElectricalPins()
{
    RStructuralPiezoCoupled1DInput noMechanicalAnchor = oneRod();
    noMechanicalAnchor.displacementPrescribed = {false,false};
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoCoupled1D::solve(noMechanicalAnchor),std::invalid_argument);

    RStructuralPiezoCoupled1DInput noElectricGauge = oneRod();
    noElectricGauge.potentialPrescribed = {false,false};
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoCoupled1D::solve(noElectricGauge),std::invalid_argument);
}

void TestStructuralPiezoCoupled1D::failClosedOnUnderconvergenceAndInvalidInputs()
{
    RStructuralPiezoCoupled1DInput underconverged = oneRod();
    underconverged.potentialPrescribed[1] = false;
    underconverged.prescribedConductorCharges[1] = 3.0;
    underconverged.gmresRestart = 1;
    underconverged.gmresCycles = 1;
    underconverged.gmresRelativeTolerance = 1.0e-14;
    underconverged.mechanicalAbsoluteTolerance = 1.0e-30;
    underconverged.mechanicalRelativeTolerance = 1.0e-30;
    underconverged.electricAbsoluteTolerance = 1.0e-30;
    underconverged.electricRelativeTolerance = 1.0e-30;
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoCoupled1D::solve(underconverged),std::runtime_error);

    RStructuralPiezoCoupled1DInput nonfinite = oneRod();
    nonfinite.members[0].coupling = std::numeric_limits<double>::infinity();
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoCoupled1D::solve(nonfinite),std::invalid_argument);

    RStructuralPiezoCoupled1DInput mixedBoundary = oneRod();
    mixedBoundary.potentialPrescribed[1] = true;
    mixedBoundary.prescribedConductorCharges[1] = 1.0;
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoCoupled1D::solve(mixedBoundary),std::invalid_argument);
}

QTEST_APPLESS_MAIN(TestStructuralPiezoCoupled1D)

#include "tst_structural_piezo_coupled_1d.moc"

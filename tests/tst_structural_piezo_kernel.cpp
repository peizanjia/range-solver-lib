#include <QtTest>

#include <cmath>
#include <limits>
#include <stdexcept>

#include "rstructural_piezo_kernel.h"

class TestStructuralPiezoKernel : public QObject
{
    Q_OBJECT

private slots:
    void analyticEnthalpyResponseAndStiffnesses();
    void enthalpyGradientAndMixedDerivatives();
    void fixedVoltageFreeExtensionAndPolingReversal();
    void dimensionalScaling();
    void rejectInvalidAndNonfiniteInputsAndOutputs();
};

void TestStructuralPiezoKernel::analyticEnthalpyResponseAndStiffnesses()
{
    const auto result = RStructuralPiezoKernel::evaluate(0.3,-0.7,10.0,2.0,4.0,1.0,1.0);
    QVERIFY(std::abs(result.enthalpy-(-0.95)) < 1.0e-14);
    QVERIFY(std::abs(result.mechanicalForce-1.6) < 1.0e-14);
    QVERIFY(std::abs(result.dielectricDisplacementFlux-3.4) < 1.0e-14);
    QCOMPARE(result.strain,0.3);
    QCOMPARE(result.electricField,0.7);
    QVERIFY(std::abs(result.stress-1.6) < 1.0e-14);
    QVERIFY(std::abs(result.electricDisplacement-3.4) < 1.0e-14);
    QCOMPARE(result.hessian[0][0],10.0);
    QCOMPARE(result.hessian[0][1],2.0);
    QCOMPARE(result.hessian[1][0],2.0);
    QCOMPARE(result.hessian[1][1],-4.0);
    QCOMPARE(result.voltageFixedStiffness,10.0);
    QCOMPARE(result.openCircuitStiffness,11.0);
    QVERIFY(result.openCircuitStiffness >= result.voltageFixedStiffness);
    const auto openCircuit = RStructuralPiezoKernel::evaluate(0.25,2.0*0.25/4.0,
        10.0,2.0,4.0,1.0,1.0);
    QVERIFY(std::abs(openCircuit.dielectricDisplacementFlux) < 1.0e-14);
    QVERIFY(std::abs(openCircuit.mechanicalForce-11.0*0.25) < 1.0e-14);
    // dielectricDisplacementFlux is A*D at the right dielectric face; a
    // conductor electrode's free charge has the opposite sign.
    QCOMPARE(result.dielectricDisplacementFlux,1.0*result.electricDisplacement);
}

void TestStructuralPiezoKernel::enthalpyGradientAndMixedDerivatives()
{
    constexpr double cE = 10.0;
    constexpr double e = -1.7;
    constexpr double epsS = 3.2;
    constexpr double area = 2.4;
    constexpr double length = 0.8;
    constexpr double u = 0.12;
    constexpr double voltage = 1.1;
    const auto base = RStructuralPiezoKernel::evaluate(u,voltage,cE,e,epsS,area,length);
    const double h = 1.0e-6;
    const auto uPlus = RStructuralPiezoKernel::evaluate(u+h,voltage,cE,e,epsS,area,length);
    const auto uMinus = RStructuralPiezoKernel::evaluate(u-h,voltage,cE,e,epsS,area,length);
    const auto vPlus = RStructuralPiezoKernel::evaluate(u,voltage+h,cE,e,epsS,area,length);
    const auto vMinus = RStructuralPiezoKernel::evaluate(u,voltage-h,cE,e,epsS,area,length);

    const double dHdu = (uPlus.enthalpy-uMinus.enthalpy)/(2.0*h);
    const double dHdV = (vPlus.enthalpy-vMinus.enthalpy)/(2.0*h);
    QVERIFY(std::abs(dHdu-base.mechanicalForce) < 2.0e-8);
    QVERIFY(std::abs(dHdV-base.dielectricDisplacementFlux) < 2.0e-8);

    const double dFdu = (uPlus.mechanicalForce-uMinus.mechanicalForce)/(2.0*h);
    const double dFdV = (vPlus.mechanicalForce-vMinus.mechanicalForce)/(2.0*h);
    const double dQdu = (uPlus.dielectricDisplacementFlux-uMinus.dielectricDisplacementFlux)/(2.0*h);
    const double dQdV = (vPlus.dielectricDisplacementFlux-vMinus.dielectricDisplacementFlux)/(2.0*h);
    QVERIFY(std::abs(dFdu-base.hessian[0][0]) < 2.0e-8);
    QVERIFY(std::abs(dFdV-base.hessian[0][1]) < 2.0e-8);
    QVERIFY(std::abs(dQdu-base.hessian[1][0]) < 2.0e-8);
    QVERIFY(std::abs(dQdV-base.hessian[1][1]) < 2.0e-8);
    QCOMPARE(base.hessian[0][1],base.hessian[1][0]);
}

void TestStructuralPiezoKernel::fixedVoltageFreeExtensionAndPolingReversal()
{
    constexpr double cE = 10.0;
    constexpr double e = 2.0;
    constexpr double epsS = 4.0;
    constexpr double area = 1.0;
    constexpr double length = 1.0;
    constexpr double voltage = 1.5;
    const double freeExtension = -e*voltage/cE;
    const auto free = RStructuralPiezoKernel::evaluate(
        freeExtension,voltage,cE,e,epsS,area,length);
    QVERIFY(std::abs(free.mechanicalForce) < 1.0e-14);

    // Reversing poling reverses the coupling terms while preserving both
    // stiffnesses, including the open/short-circuit stiffness difference.
    const auto positive = RStructuralPiezoKernel::evaluate(0.3,-0.7,cE,e,epsS,area,length);
    const auto reversed = RStructuralPiezoKernel::evaluate(0.3,-0.7,cE,-e,epsS,area,length);
    QCOMPARE(positive.voltageFixedStiffness,reversed.voltageFixedStiffness);
    QCOMPARE(positive.openCircuitStiffness,reversed.openCircuitStiffness);
    QCOMPARE(positive.hessian[0][1],-reversed.hessian[0][1]);
    QCOMPARE(positive.hessian[1][0],-reversed.hessian[1][0]);
    QVERIFY(std::abs((positive.mechanicalForce+reversed.mechanicalForce)-
                     2.0*cE*0.3) < 1.0e-14);
    QVERIFY(std::abs((positive.dielectricDisplacementFlux-
                      reversed.dielectricDisplacementFlux)-2.0*e*0.3) < 1.0e-14);
}

void TestStructuralPiezoKernel::dimensionalScaling()
{
    const auto base = RStructuralPiezoKernel::evaluate(0.12,-0.8,10.0,2.0,4.0,2.0,4.0);
    const auto scaled = RStructuralPiezoKernel::evaluate(0.12,-0.8,10.0,2.0,4.0,6.0,8.0);
    // A/L changes by 3/2; strain and field change by 1/2.
    QVERIFY(std::abs(scaled.enthalpy-1.5*base.enthalpy) < 1.0e-14);
    QVERIFY(std::abs(scaled.mechanicalForce-1.5*base.mechanicalForce) < 1.0e-14);
    QVERIFY(std::abs(scaled.dielectricDisplacementFlux-1.5*base.dielectricDisplacementFlux) < 1.0e-14);
    QVERIFY(std::abs(scaled.hessian[0][1]-1.5*base.hessian[0][1]) < 1.0e-14);
    QVERIFY(std::abs(scaled.voltageFixedStiffness-1.5*base.voltageFixedStiffness) < 1.0e-14);
    QVERIFY(std::abs(scaled.openCircuitStiffness-1.5*base.openCircuitStiffness) < 1.0e-14);
    QVERIFY(std::abs(scaled.strain-0.5*base.strain) < 1.0e-15);
    QVERIFY(std::abs(scaled.electricField-0.5*base.electricField) < 1.0e-15);
    QVERIFY(std::abs(scaled.stress-0.5*base.stress) < 1.0e-14);
    QVERIFY(std::abs(scaled.electricDisplacement-0.5*base.electricDisplacement) < 1.0e-14);
}

void TestStructuralPiezoKernel::rejectInvalidAndNonfiniteInputsAndOutputs()
{
    constexpr double validC = 10.0;
    constexpr double validE = 2.0;
    constexpr double validEps = 4.0;
    constexpr double validArea = 1.0;
    constexpr double validLength = 1.0;
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoKernel::evaluate(0.0,0.0,0.0,validE,validEps,
        validArea,validLength),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoKernel::evaluate(0.0,0.0,-validC,validE,validEps,
        validArea,validLength),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoKernel::evaluate(0.0,0.0,validC,validE,0.0,
        validArea,validLength),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoKernel::evaluate(0.0,0.0,validC,validE,-validEps,
        validArea,validLength),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoKernel::evaluate(0.0,0.0,validC,validE,validEps,
        0.0,validLength),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoKernel::evaluate(0.0,0.0,validC,validE,validEps,
        -validArea,validLength),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoKernel::evaluate(0.0,0.0,validC,validE,validEps,
        validArea,0.0),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoKernel::evaluate(0.0,0.0,validC,validE,validEps,
        validArea,-validLength),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoKernel::evaluate(0.0,0.0,validC,
        std::numeric_limits<double>::infinity(),validEps,validArea,validLength),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoKernel::evaluate(
        std::numeric_limits<double>::quiet_NaN(),0.0,validC,validE,validEps,validArea,validLength),
        std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoKernel::evaluate(0.0,
        std::numeric_limits<double>::infinity(),validC,validE,validEps,validArea,validLength),
        std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoKernel::evaluate(0.0,0.0,validC,validE,validEps,
        std::numeric_limits<double>::max(),std::numeric_limits<double>::min()),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoKernel::evaluate(
        std::numeric_limits<double>::max(),0.0,std::numeric_limits<double>::max(),validE,
        validEps,validArea,validLength),std::invalid_argument);
}

QTEST_APPLESS_MAIN(TestStructuralPiezoKernel)

#include "tst_structural_piezo_kernel.moc"

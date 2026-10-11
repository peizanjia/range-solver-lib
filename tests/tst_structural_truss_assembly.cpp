#include <QtTest>
#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include "rstructural_truss_assembly.h"

namespace
{
using V3 = RStructuralTrussVector3;
using V6 = RStructuralTrussVector6;
using Law = RStructuralTrussLaw;
using Input = RStructuralTrussInput;

bool near(double actual, double expected, double tolerance = 2.0e-9)
{
    return std::abs(actual-expected) <= tolerance*(1.0+std::abs(expected));
}
RStructuralTrussElementResult element(const V6 &x, const Law &law, const Input &input)
{
    return RStructuralTrussAssembly::evaluateElement(
        V3{x[0],x[1],x[2]},V3{x[3],x[4],x[5]},law,input);
}

// Independent scalar potentials, differentiated numerically below. These do
// not use element responses, force/tangent formulas or the cable kernel.
double potential(const V6 &x, int kind)
{
    const double length = std::sqrt(std::pow(x[3]-x[0],2)+
                                   std::pow(x[4]-x[1],2)+std::pow(x[5]-x[2],2));
    const double delta = length-2.0;
    if (kind == 0) return 0.5*40.0*delta*delta;
    if (kind == 1) return 0.5*40.0*delta*delta+3.0*delta*1.7-0.5*2.0*1.7*1.7;
    if (kind == 2)
    {
        const double voltage = (0.8+3.0*delta)/2.0;
        return 0.5*40.0*delta*delta+0.5*2.0*voltage*voltage;
    }
    const double tension = std::max(0.0,4.0+40.0*(length-2.0));
    return tension*tension/80.0;
}
V6 rotate(const V6 &x)
{
    // Proper orthogonal rotation, with an oblique axis (cyclic permutation).
    return V6{x[2],x[0],x[1],x[5],x[3],x[4]};
}
}

class TestStructuralTrussAssembly : public QObject
{
    Q_OBJECT
private slots:
    void passiveCompressionAndVoltageBoundaries();
    void chargeBoundaryAndRepresentableCoupling();
    void cableBranchesAndLengthWork();
    void representableEnergyAndLengthDerivative();
    void independentPotentialGradientsAndTangents();
    void endpointReversalAndRigidMotion();
    void commonMixedAssembly();
    void rejectMalformedInputsWithoutMutation();
};

void TestStructuralTrussAssembly::passiveCompressionAndVoltageBoundaries()
{
    const auto compressed = element(V6{0,0,0,0.8,0,0},RStructuralTrussPassiveLaw{1,100},RStructuralTrussNoInput{});
    QVERIFY(near(compressed.axialForce,-20));
    QVERIFY(near(compressed.storedEnergy,2));
    QVERIFY(near(compressed.tangent[1][1],-25)); // Negative transverse geometric stiffness.
    QCOMPARE(compressed.axialTangent,100.0);
    QVERIFY(!compressed.electrical);
    QVERIFY(!compressed.cableTaut);
    QVERIFY(!compressed.storedEnergyLengthDerivative);

    const Law piezo = RStructuralTrussPiezoLaw{1,100,2,0.5};
    const Input voltage = RStructuralTrussVoltage{3};
    const auto free = element(V6{0,0,0,0.94,0,0},piezo,voltage);
    QVERIFY(near(free.axialForce,0)); // delta_free = -g V/k.
    QVERIFY(near(free.electrical->conductorCharge,1.62));
    const auto blocked = element(V6{0,0,0,1,0,0},piezo,voltage);
    QCOMPARE(blocked.axialForce,6.0);
    QCOMPARE(blocked.axialTangent,100.0);
    QCOMPARE(blocked.storedEnergy,2.25);
    QCOMPARE(blocked.mechanicalPotential,-2.25);
    QCOMPARE(blocked.electrical->voltage,3.0);
    QCOMPARE(blocked.electrical->conductorCharge,1.5);
    QCOMPARE(blocked.electrical->enthalpy,-2.25);
    const auto shorted = element(V6{0,0,0,1.1,0,0},piezo,RStructuralTrussVoltage{0});
    QVERIFY(near(shorted.axialForce,10));
    QVERIFY(near(shorted.electrical->conductorCharge,-0.2));
}

void TestStructuralTrussAssembly::representableEnergyAndLengthDerivative()
{
    // Each complete result is finite, although squaring N in double overflows.
    // This must also work on MSVC where long double has no wider exponent.
    const auto cable = element(V6{0,0,0,1,0,0},
        RStructuralTrussCableLaw{1.0e200,1.0e200},RStructuralTrussCableLength{1});
    QVERIFY(near(cable.axialForce/1.0e200,1.0));
    QVERIFY(near(cable.storedEnergy/1.0e200,0.5));
    QVERIFY(near(*cable.storedEnergyLengthDerivative/1.0e200,-0.5));
    QVERIFY(near(cable.tangent[0][0]/1.0e200,1.0));
    QVERIFY(near(cable.tangent[1][1]/1.0e200,1.0));

    // The two derivative terms overflow separately and cancel to finite dU/dL0.
    const auto cancellation = element(V6{0,0,0,1,0,0},
        RStructuralTrussCableLaw{1.25e308,1.0e308},RStructuralTrussCableLength{0.75});
    QVERIFY(near(cancellation.storedEnergy/1.0e308,361.0/384.0));
    QVERIFY(near(*cancellation.storedEnergyLengthDerivative/1.0e308,-247.0/288.0));

    // Do not lose a representable subnormal energy by rounding 0.5*N first.
    const double tiny = std::numeric_limits<double>::denorm_min();
    const auto passive = element(V6{0,0,0,101.25,0,0},
        RStructuralTrussPassiveLaw{100,100*tiny},RStructuralTrussNoInput{});
    QCOMPARE(passive.storedEnergy,tiny);
}

void TestStructuralTrussAssembly::chargeBoundaryAndRepresentableCoupling()
{
    const Law law = RStructuralTrussPiezoLaw{1,100,2,0.5};
    const auto open = element(V6{0,0,0,1.1,0,0},law,RStructuralTrussCharge{0});
    QVERIFY(near(open.axialForce,10.8));
    QCOMPARE(open.axialTangent,108.0);
    QVERIFY(near(open.electrical->voltage,0.4));
    QCOMPARE(open.electrical->conductorCharge,0.0);
    QVERIFY(near(open.storedEnergy,0.54));
    QCOMPARE(open.mechanicalPotential,open.storedEnergy);
    QVERIFY(near(open.electrical->enthalpy,open.storedEnergy));
    const auto charged = element(V6{0,0,0,1,0,0},law,RStructuralTrussCharge{1});
    QCOMPARE(charged.electrical->voltage,2.0);
    QCOMPARE(charged.axialForce,4.0);
    QCOMPARE(charged.storedEnergy,1.0);

    // g*g overflows double, but g^2/Ce = 1e200 is representable.
    const auto large = element(V6{0,0,0,1,0,0},RStructuralTrussPiezoLaw{1,1,1e200,1e200},RStructuralTrussCharge{0});
    QVERIFY(near(large.axialTangent,1e200));
    QVERIFY_EXCEPTION_THROWN(element(V6{0,0,0,1,0,0},
        RStructuralTrussPiezoLaw{1,1,1e200,1},RStructuralTrussCharge{0}),std::invalid_argument);

    // Independent references are dimensionless arithmetic followed by one
    // multiplication by 1e308. No overflowing product/sum appears in them.
    const Law huge = RStructuralTrussPiezoLaw{1,1,1e308,1.7e308};
    const double delta = 1.01-1.0;
    const double v = (1.79+delta)/1.7;
    const auto chargeDriven = element(V6{0,0,0,1.01,0,0},huge,RStructuralTrussCharge{1.79e308});
    QVERIFY(near(chargeDriven.electrical->voltage,v));
    QVERIFY(near(chargeDriven.axialForce,v*1e308));
    QVERIFY(near(chargeDriven.axialTangent,(1.0/1.7)*1e308));
    QVERIFY(near(chargeDriven.storedEnergy,(0.5*1.7*v*v)*1e308));
    QVERIFY(near(chargeDriven.electrical->enthalpy,(delta*v-0.5*1.7*v*v)*1e308));
    const auto voltageDriven = element(V6{0,0,0,1.01,0,0},huge,RStructuralTrussVoltage{1.06});
    QVERIFY(near(voltageDriven.electrical->conductorCharge,(1.7*1.06-delta)*1e308));
    QVERIFY(near(voltageDriven.axialForce,1.06e308));
    QVERIFY(near(voltageDriven.storedEnergy,(0.5*1.7*1.06*1.06)*1e308));
    QVERIFY(near(voltageDriven.electrical->enthalpy,(delta*1.06-0.5*1.7*1.06*1.06)*1e308));
    const auto cancelled = element(V6{0,0,0,2.5,0,0},RStructuralTrussPiezoLaw{1,1,1.4e308,1.7e308},RStructuralTrussCharge{-1.7e308});
    const double cv = (-1.7+1.4*1.5)/1.7;
    QVERIFY(near(cancelled.electrical->voltage,cv));
    QVERIFY(near(cancelled.axialForce,(1.4*cv)*1e308));
    QVERIFY(near(cancelled.axialTangent,(1.4*1.4/1.7)*1e308));
    QVERIFY(near(cancelled.storedEnergy,(0.5*1.7*cv*cv)*1e308));
    QVERIFY(near(cancelled.electrical->enthalpy,(1.4*1.5*cv-0.5*1.7*cv*cv)*1e308));
    // Same cancellation under voltage control: Ce*V and g*delta both exceed
    // double max, while Q2, N, U, H and the full tangent all remain finite.
    const auto cancelledVoltage = element(V6{0,0,0,2.5,0,0},RStructuralTrussPiezoLaw{1,1,1.4e308,1.7e308},RStructuralTrussVoltage{1.25});
    QVERIFY(near(cancelledVoltage.electrical->conductorCharge,(1.7*1.25-1.4*1.5)*1e308));
    QVERIFY(near(cancelledVoltage.axialForce,1.75e308));
    QVERIFY(near(cancelledVoltage.storedEnergy,(0.5*1.7*1.25*1.25)*1e308));
    const double tinyCe=std::numeric_limits<double>::denorm_min();
    const auto tiny = element(V6{0,0,0,1,0,0},RStructuralTrussPiezoLaw{1,1,0,tinyCe},RStructuralTrussVoltage{1e154});
    QVERIFY(near(tiny.storedEnergy,((tinyCe*1e154)*1e154)*0.5,1e-25));
    QVERIFY(tiny.storedEnergy>0.0);
    // Public outputs that truly exceed double max must still be rejected.
    QVERIFY_EXCEPTION_THROWN(element(V6{0,0,0,1.8,0,0},huge,RStructuralTrussCharge{1.79e308}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(element(V6{0,0,0,1,0,0},huge,RStructuralTrussVoltage{1.5}),std::invalid_argument);
}

void TestStructuralTrussAssembly::cableBranchesAndLengthWork()
{
    const Law law = RStructuralTrussCableLaw{10,100};
    for (double length : {0.8,0.9})
    {
        // 0.9 is not exactly the kink in binary arithmetic; use exact data below.
        const auto slack = element(V6{0,0,0,length,0,0},law,RStructuralTrussCableLength{1});
        QVERIFY(slack.cableTaut.has_value());
        if (length == 0.8)
        {
            QVERIFY(!*slack.cableTaut);
            QCOMPARE(slack.axialForce,0.0);
            QCOMPARE(slack.axialTangent,0.0);
            QCOMPARE(*slack.storedEnergyLengthDerivative,0.0);
        }
        QVERIFY(!slack.electrical);
    }
    const auto kink = element(V6{0,0,0,0.75,0,0},RStructuralTrussCableLaw{25,100},RStructuralTrussCableLength{1});
    QVERIFY(!*kink.cableTaut);
    QCOMPARE(kink.storedEnergy,0.0);
    QCOMPARE(*kink.storedEnergyLengthDerivative,0.0);
    for (const auto &row : kink.tangent) for (double value : row) QCOMPARE(value,0.0);
    const V6 x{0,0,0,1.2,0,0};
    const auto taut = element(x,law,RStructuralTrussCableLength{1});
    QVERIFY(*taut.cableTaut);
    QVERIFY(near(taut.axialForce,30));
    QVERIFY(near(taut.storedEnergy,4.5));
    QVERIFY(near(*taut.storedEnergyLengthDerivative,-31.5));
    const auto scalarEnergy = [](double l0) {
        const double n = 10+100*(1.2/l0-1);
        return n*n*l0/200;
    };
    const double h = 1e-6;
    QVERIFY(near(*taut.storedEnergyLengthDerivative,(scalarEnergy(1+h)-scalarEnergy(1-h))/(2*h),1e-7));
}

void TestStructuralTrussAssembly::independentPotentialGradientsAndTangents()
{
    const std::array<Law,4> laws{RStructuralTrussPassiveLaw{2,80},
        RStructuralTrussPiezoLaw{2,40,3,2},RStructuralTrussPiezoLaw{2,40,3,2},
        RStructuralTrussCableLaw{4,80}};
    const std::array<Input,4> inputs{RStructuralTrussNoInput{},RStructuralTrussVoltage{1.7},
        RStructuralTrussCharge{0.8},RStructuralTrussCableLength{2}};
    // Both tensile and compressive oblique states; taut cable case kept away from kink.
    for (const V6 x : {V6{0.1,-0.2,0.3,2.4,0.8,0.9},V6{0.1,-0.2,0.3,1.1,0.5,0.8}})
        for (int kind=0;kind<4;++kind)
        {
            const auto base = element(x,laws[kind],inputs[kind]);
            QVERIFY(near(base.mechanicalPotential,potential(x,kind)));
            const double h = 1e-5;
            for (std::size_t i=0;i<6;++i)
            {
                V6 plus=x,minus=x;
                plus[i]+=h; minus[i]-=h;
                const double gradient = (potential(plus,kind)-potential(minus,kind))/(2*h);
                QVERIFY(near(base.internalForce[i],gradient,2e-7));
                const auto fp=element(plus,laws[kind],inputs[kind]).internalForce;
                const auto fm=element(minus,laws[kind],inputs[kind]).internalForce;
                for (std::size_t j=0;j<6;++j)
                {
                    QVERIFY(near(base.tangent[j][i],(fp[j]-fm[j])/(2*h),2e-7));
                    QCOMPARE(base.tangent[i][j],base.tangent[j][i]);
                    V6 pp=x,pm=x,mp=x,mm=x;
                    pp[i]+=h; pp[j]+=h;
                    pm[i]+=h; pm[j]-=h;
                    mp[i]-=h; mp[j]+=h;
                    mm[i]-=h; mm[j]-=h;
                    const double hessian=(potential(pp,kind)-potential(pm,kind)-potential(mp,kind)+potential(mm,kind))/(4*h*h);
                    QVERIFY(near(base.tangent[i][j],hessian,1e-4));
                }
            }
        }
}

void TestStructuralTrussAssembly::endpointReversalAndRigidMotion()
{
    const V6 x{0.1,-0.2,0.3,1.6,0.8,0.9};
    const std::array<Law,4> laws{RStructuralTrussPassiveLaw{1,70},
        RStructuralTrussPiezoLaw{1,70,-3,2},RStructuralTrussCableLaw{4,70},RStructuralTrussPiezoLaw{1,70,-3,2}};
    const std::array<Input,4> inputs{RStructuralTrussNoInput{},RStructuralTrussVoltage{1.7},RStructuralTrussCableLength{1},RStructuralTrussCharge{0.4}};
    for (std::size_t kind=0;kind<4;++kind)
    {
        const auto original = element(x,laws[kind],inputs[kind]);
        const auto reversed = element(V6{x[3],x[4],x[5],x[0],x[1],x[2]},laws[kind],inputs[kind]);
        QCOMPARE(original.axialForce,reversed.axialForce);
        QCOMPARE(original.storedEnergy,reversed.storedEnergy);
        QCOMPARE(original.mechanicalPotential,reversed.mechanicalPotential);
        if (original.electrical)
        {
            QCOMPARE(original.electrical->voltage,reversed.electrical->voltage);
            QCOMPARE(original.electrical->conductorCharge,reversed.electrical->conductorCharge);
            QCOMPARE(original.electrical->enthalpy,reversed.electrical->enthalpy);
        }
        V6 translated=x;
        for (std::size_t i=0;i<6;++i) translated[i] += i%3+2;
        const auto shifted=element(translated,laws[kind],inputs[kind]);
        const auto rotated=element(rotate(x),laws[kind],inputs[kind]);
        const V6 rotatedForce=rotate(original.internalForce);
        QVERIFY(near(original.storedEnergy,rotated.storedEnergy));
        const std::array<std::size_t,6> p{2,0,1,5,3,4};
        for (std::size_t i=0;i<6;++i)
        {
            QVERIFY(near(original.internalForce[i],reversed.internalForce[(i+3)%6]));
            QVERIFY(near(original.internalForce[i],shifted.internalForce[i]));
            QVERIFY(near(rotatedForce[i],rotated.internalForce[i]));
            for (std::size_t j=0;j<6;++j)
            {
                QVERIFY(near(original.tangent[i][j],reversed.tangent[(i+3)%6][(j+3)%6]));
                QVERIFY(near(rotated.tangent[i][j],original.tangent[p[i]][p[j]]));
            }
        }
    }
    // Electrode relabeling is a different operation from endpoint reversal:
    // g and prescribed V/Q2 change sign together, while mechanical responses
    // and physical U remain invariant. Electrical V/Q2 change sign.
    for (const Input input : {Input{RStructuralTrussVoltage{1.7}},Input{RStructuralTrussCharge{0.4}}})
    {
        const Input invertedInput = std::holds_alternative<RStructuralTrussVoltage>(input)
            ? Input{RStructuralTrussVoltage{-1.7}} : Input{RStructuralTrussCharge{-0.4}};
        const auto original=element(x,RStructuralTrussPiezoLaw{1,70,3,2},input);
        const auto inverted=element(x,RStructuralTrussPiezoLaw{1,70,-3,2},invertedInput);
        QCOMPARE(original.axialForce,inverted.axialForce);
        QCOMPARE(original.axialTangent,inverted.axialTangent);
        QCOMPARE(original.storedEnergy,inverted.storedEnergy);
        QCOMPARE(original.mechanicalPotential,inverted.mechanicalPotential);
        QCOMPARE(original.electrical->voltage,-inverted.electrical->voltage);
        QCOMPARE(original.electrical->conductorCharge,-inverted.electrical->conductorCharge);
        QCOMPARE(original.electrical->enthalpy,inverted.electrical->enthalpy);
        for (std::size_t i=0;i<6;++i)
        {
            QCOMPARE(original.internalForce[i],inverted.internalForce[i]);
            for (std::size_t j=0;j<6;++j) QCOMPARE(original.tangent[i][j],inverted.tangent[i][j]);
        }
    }
}

void TestStructuralTrussAssembly::commonMixedAssembly()
{
    const std::vector<V3> nodes{{0,0,0},{1,0,0},{2,0,0},{3,0,0},{4,0,0}};
    std::vector<double> u(15,0.0);
    u[6]=0.1; u[9]=0.2; // Current segment lengths 1,1.1,1.1.
    const std::vector<RStructuralTrussMember> members{
        {19,0,1,RStructuralTrussPassiveLaw{0.9,90}},
        {7,1,2,RStructuralTrussPiezoLaw{1,50,2,0.5}},
        {42,2,3,RStructuralTrussCableLaw{0,100}}};
    const std::vector<Input> inputs{RStructuralTrussNoInput{},RStructuralTrussVoltage{2.5},RStructuralTrussCableLength{1}};
    const auto result=RStructuralTrussAssembly::evaluate(nodes,u,members,inputs);
    QCOMPARE(result.tangent.getNRows(),uint(15));
    QCOMPARE(result.memberResponses.size(),std::size_t(3));
    // All three members carry 10 N: common interior nodes are in equilibrium.
    const std::array<double,15> force{-10,0,0,0,0,0,0,0,0,10,0,0,0,0,0};
    for (std::size_t i=0;i<15;++i) QVERIFY(near(result.internalForces[i],force[i]));
    double reference[15][15]{};
    // Independent analytic 3x3 blocks on the x-axis; no element-kernel calls.
    const double longitudinal[3]{100,50,100};
    const double transverse[3]{10,10/1.1,10/1.1};
    for (std::size_t e=0;e<3;++e)
        for (std::size_t i=0;i<3;++i)
        {
            const double a=i==0 ? longitudinal[e] : transverse[e];
            const std::size_t left=3*e+i,right=left+3;
            reference[left][left]+=a; reference[right][right]+=a;
            reference[left][right]-=a; reference[right][left]-=a;
        }
    for (uint i=0;i<15;++i) for (uint j=0;j<15;++j)
        QVERIFY(near(result.tangent.findValue(i,j),reference[i][j]));
    QVERIFY(near(result.storedEnergy,0.5+0.25+1.5625+0.5));
    QVERIFY(near(result.mechanicalPotential,0.5+0.25+0.5-1.5625+0.5));
    QVERIFY(!result.memberResponses[0].electrical);
    QVERIFY(result.memberResponses[1].electrical);
    QVERIFY(!result.memberResponses[2].electrical);
    // Null rigid-translation direction and action/reaction on each axis.
    for (uint row=0;row<15;++row) for (uint axis=0;axis<3;++axis)
    {
        double sum=0;
        for (uint node=0;node<5;++node) sum+=result.tangent.findValue(row,3*node+axis);
        QVERIFY(near(sum,0));
    }

    // Non-collinear network: oblique compressed passive member, slack cable,
    // and charge-driven piezo share all global DOFs in a single assembly.
    const std::vector<V3> skewNodes{{0,0,0},{0.6,0.8,0},{0.6,0.8,1}};
    const std::vector<RStructuralTrussMember> skewMembers{
        {101,2,0,RStructuralTrussCableLaw{0,100}},
        {102,0,1,RStructuralTrussPassiveLaw{1.2,120}},
        {103,1,2,RStructuralTrussPiezoLaw{0.9,50,2,0.5}}};
    const auto skew=RStructuralTrussAssembly::evaluate(skewNodes,std::vector<double>(9,0.0),skewMembers,
        {RStructuralTrussCableLength{2},RStructuralTrussNoInput{},RStructuralTrussCharge{0.1}});
    const double skewForces[9]{12,16,0,-12,-16,-6.2,0,0,6.2};
    const double passiveBlock[3][3]{{23.2,57.6,0},{57.6,56.8,0},{0,0,-20}};
    const double piezoBlock[3][3]{{6.2,0,0},{0,6.2,0},{0,0,58}};
    double skewReference[9][9]{};
    for (uint i=0;i<3;++i) for (uint j=0;j<3;++j)
    {
        skewReference[i][j]+=passiveBlock[i][j];
        skewReference[i+3][j+3]+=passiveBlock[i][j];
        skewReference[i][j+3]-=passiveBlock[i][j];
        skewReference[i+3][j]-=passiveBlock[i][j];
        skewReference[i+3][j+3]+=piezoBlock[i][j];
        skewReference[i+6][j+6]+=piezoBlock[i][j];
        skewReference[i+3][j+6]-=piezoBlock[i][j];
        skewReference[i+6][j+3]-=piezoBlock[i][j];
    }
    for (uint i=0;i<9;++i)
    {
        QVERIFY(near(skew.internalForces[i],skewForces[i]));
        for (uint j=0;j<9;++j) QVERIFY(near(skew.tangent.findValue(i,j),skewReference[i][j]));
    }
    QVERIFY(near(skew.storedEnergy,2.34));
    QVERIFY(near(skew.mechanicalPotential,2.34));
    QVERIFY(!*skew.memberResponses[0].cableTaut);
    QVERIFY(near(skew.memberResponses[1].axialForce,-20));
    QVERIFY(near(skew.memberResponses[2].electrical->voltage,0.6));
    QVERIFY(near(skew.memberResponses[2].electrical->conductorCharge,0.1));
    for (uint axis=0;axis<3;++axis)
        QVERIFY(near(skew.internalForces[axis]+skew.internalForces[axis+3]+skew.internalForces[axis+6],0));
}

void TestStructuralTrussAssembly::rejectMalformedInputsWithoutMutation()
{
    const V6 x{0,0,0,1,0,0};
    const double nan=std::numeric_limits<double>::quiet_NaN();
    const double max=std::numeric_limits<double>::max();
    const std::array<Law,3> laws{RStructuralTrussPassiveLaw{1,100},RStructuralTrussPiezoLaw{1,100,2,0.5},RStructuralTrussCableLaw{0,100}};
    const std::array<Input,4> inputs{RStructuralTrussNoInput{},RStructuralTrussVoltage{1},RStructuralTrussCharge{1},RStructuralTrussCableLength{1}};
    for (std::size_t l=0;l<3;++l) for (std::size_t i=0;i<4;++i)
        if (!(l==0 && i==0) && !(l==1 && (i==1 || i==2)) && !(l==2 && i==3))
            QVERIFY_EXCEPTION_THROWN(element(x,laws[l],inputs[i]),std::invalid_argument);
    for (const Law law : {Law{RStructuralTrussPassiveLaw{0,100}},Law{RStructuralTrussPassiveLaw{1,nan}},
        Law{RStructuralTrussPassiveLaw{1e-300,max}},Law{RStructuralTrussPiezoLaw{1,0,2,1}},
        Law{RStructuralTrussPiezoLaw{1,100,nan,1}},Law{RStructuralTrussPiezoLaw{1,100,2,0}}})
    {
        const Input input=std::holds_alternative<RStructuralTrussPassiveLaw>(law) ? Input{RStructuralTrussNoInput{}} : Input{RStructuralTrussVoltage{1}};
        QVERIFY_EXCEPTION_THROWN(element(x,law,input),std::invalid_argument);
    }
    QVERIFY_EXCEPTION_THROWN(element(x,laws[1],RStructuralTrussVoltage{nan}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(element(x,laws[1],RStructuralTrussCharge{nan}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(element(x,laws[1],RStructuralTrussVoltage{max}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(element(x,laws[2],RStructuralTrussCableLength{0}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(element(x,laws[2],RStructuralTrussCableLength{nan}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(element(x,RStructuralTrussCableLaw{-1,100},RStructuralTrussCableLength{1}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(element(x,RStructuralTrussCableLaw{nan,100},RStructuralTrussCableLength{1}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(element(x,RStructuralTrussCableLaw{0,0},RStructuralTrussCableLength{1}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(element(x,RStructuralTrussCableLaw{0,max},RStructuralTrussCableLength{1e-300}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(element(V6{},laws[0],inputs[0]),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(element(V6{nan,0,0,1,0,0},laws[0],inputs[0]),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(element(V6{-max,0,0,max,0,0},laws[0],inputs[0]),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(element(V6{0,0,0,max,max,0},laws[0],inputs[0]),std::invalid_argument);

    const std::vector<V3> nodes{{0,0,0},{1,0,0}};
    const std::vector<double> u(6,0.0);
    const std::vector<RStructuralTrussMember> members{{1,0,1,laws[0]}};
    const std::vector<Input> applied{RStructuralTrussNoInput{}};
    auto candidate=RStructuralTrussAssembly::evaluate(nodes,u,members,applied);
    const auto oldForce=candidate.internalForces;
    const double oldEnergy=candidate.storedEnergy;
    auto bad=members;
    bad.push_back(bad.front());
    QVERIFY_EXCEPTION_THROWN(candidate=RStructuralTrussAssembly::evaluate(nodes,u,bad,
        std::vector<Input>{applied[0],applied[0]}),std::invalid_argument);
    QVERIFY(candidate.internalForces==oldForce);
    QCOMPARE(candidate.storedEnergy,oldEnergy);
    for (auto malformed : {RStructuralTrussMember{0,0,1,laws[0]},RStructuralTrussMember{1,0,0,laws[0]},RStructuralTrussMember{1,0,2,laws[0]}})
        QVERIFY_EXCEPTION_THROWN(RStructuralTrussAssembly::evaluate(nodes,u,{malformed},applied),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralTrussAssembly::evaluate(nodes,std::vector<double>(5),members,applied),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralTrussAssembly::evaluate(nodes,u,members,{}),std::invalid_argument);
    auto badU=u; badU[0]=nan;
    QVERIFY_EXCEPTION_THROWN(RStructuralTrussAssembly::evaluate(nodes,badU,members,applied),std::invalid_argument);
    auto badNodes=nodes; badNodes[1][0]=nan;
    QVERIFY_EXCEPTION_THROWN(RStructuralTrussAssembly::evaluate(badNodes,u,members,applied),std::invalid_argument);
    badNodes=nodes; badNodes[0][0]=max; badU=u; badU[0]=max;
    QVERIFY_EXCEPTION_THROWN(RStructuralTrussAssembly::evaluate(badNodes,badU,members,applied),std::invalid_argument);
    // Each huge compressed member remains representable, but their common
    // longitudinal tangent sum cannot be represented. Reject the whole candidate.
    const Law huge=RStructuralTrussPassiveLaw{1,1e308};
    const std::vector<V3> shortNodes{{0,0,0},{0.5,0,0}};
    const std::vector<RStructuralTrussMember> hugeMembers{{1,0,1,huge},{2,0,1,huge}};
    QVERIFY_EXCEPTION_THROWN(candidate=RStructuralTrussAssembly::evaluate(shortNodes,u,hugeMembers,
        std::vector<Input>{applied[0],applied[0]}),std::invalid_argument);
    QVERIFY(candidate.internalForces==oldForce);
    QCOMPARE(candidate.storedEnergy,oldEnergy);
    // Caller arrays are still the original immutable inputs.
    QCOMPARE(nodes[1][0],1.0);
    QVERIFY(u==std::vector<double>(6,0.0));
    QCOMPARE(members[0].physicalElementId,std::uint64_t(1));
    QVERIFY(std::holds_alternative<RStructuralTrussNoInput>(applied[0]));
}

QTEST_APPLESS_MAIN(TestStructuralTrussAssembly)
#include "tst_structural_truss_assembly.moc"

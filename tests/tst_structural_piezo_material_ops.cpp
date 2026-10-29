#include <QtTest>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "rstructural_piezo_element_operators.h"

namespace
{
bool closeRelative(double actual,double expected,double tolerance=1.0e-12)
{
    return std::abs(actual-expected) <= tolerance*std::max(std::abs(actual),std::abs(expected));
}

RStructuralPiezoContinuumMaterial coupledReference()
{
    RStructuralPiezoContinuumMaterial material;
    material.density = 7500.0;
    material.piezoelectric = true;
    // A small independently solvable anisotropic constitutive example.
    const std::array<double,6> diagonal{{10.0,8.0,4.0,3.0,2.0,1.0}};
    for (size_t i=0;i<6;++i) material.elasticity[i*6+i] = diagonal[i];
    material.elasticity[2] = material.elasticity[12] = 2.0;
    material.elasticity[8] = material.elasticity[13] = 1.0;
    material.coupling[12] = 3.0;
    material.coupling[13] = 2.0;
    material.coupling[14] = 6.0;
    material.coupling[16] = -1.0;
    material.permittivity[0] = material.permittivity[4] = material.permittivity[8] = 2.0;
    return material;
}

RStructuralPiezoElementOperators operatorReference()
{
    RStructuralPiezoElementOperators result;
    result.dofs = 2;
    result.portIds = {QStringLiteral("electrode")};
    result.mass = {2.0,0.25,0.25,1.0};
    result.stiffness = {3.0,-3.0,-3.0,3.0}; // One rigid translation.
    result.coupling = {2.0,-2.0};
    result.capacitance = {1.0e-10};
    result.totalMass = 3.5;
    return result;
}
}

class TestStructuralPiezoMaterialOps : public QObject
{
    Q_OBJECT
private slots:
    void isotropicEngineeringVoigtAndAxialReference();
    void isotropicPlaneStressReference();
    void coupledPlaneStressAndReciprocityReference();
    void componentOrderingAndNoElimination();
    void legitimateSmallPhysicalScales();
    void rejectInvalidMaterial();
    void rejectInvalidRetainedComponents();
    void layerGeometryAndElectrodeValidation();
    void elementOperatorValidation();
    void strainChargeConversionAxesAndReplacement();
    void anisotropicStrainStressEquivalence();
    void axialFiniteScales();
    void rejectCompilerInputs();
};

void TestStructuralPiezoMaterialOps::isotropicEngineeringVoigtAndAxialReference()
{
    // E=120, nu=1/4 gives lambda=48 and mu=48 exactly.
    const auto material = RStructuralPiezoContinuumMaterial::isotropic(120.0,0.25,2700.0);
    QCOMPARE(material.elasticity[0],144.0);
    QCOMPARE(material.elasticity[1],48.0);
    QCOMPARE(material.elasticity[7],144.0);
    for (uint i=3;i<6;++i) QCOMPARE(material.elasticity[i*6+i],48.0);
    QCOMPARE(material.density,2700.0);
    QVERIFY(!material.piezoelectric);
    for (double e : material.coupling) QCOMPARE(e,0.0);
    for (double eps : material.permittivity) QCOMPARE(eps,0.0);
    const auto axial = rStructuralPiezoReduceFreeStress(material,{0});
    QCOMPARE(axial.components.size(),size_t(1));
    QCOMPARE(axial.components[0],uint(0));
    QVERIFY(closeRelative(axial.elasticity[0],120.0));
    QCOMPARE(axial.coupling[0],0.0);
    QCOMPARE(axial.permittivity,0.0);
    // Longitudinal isotropic reduction is E over the entire admissible nu range.
    for (double nu : {-0.75,0.0,0.3,0.49})
    {
        const auto m = RStructuralPiezoContinuumMaterial::isotropic(70.0e9,nu,2700.0);
        QVERIFY(closeRelative(rStructuralPiezoReduceFreeStress(m,{0}).elasticity[0],70.0e9));
    }
}

void TestStructuralPiezoMaterialOps::isotropicPlaneStressReference()
{
    constexpr double young = 70.0e9, nu = 0.3;
    const auto material = RStructuralPiezoContinuumMaterial::isotropic(young,nu,2700.0);
    const auto plate = rStructuralPiezoReduceFreeStress(material,{0,1,5,4,3});
    const double normal = young/(1.0-nu*nu);
    const double shear = young/(2.0*(1.0+nu));
    for (size_t i=0;i<5;++i)
        for (size_t j=0;j<5;++j)
        {
            double expected = 0.0;
            if (i == j) expected = i < 2 ? normal : shear;
            else if (i < 2 && j < 2) expected = nu*normal;
            QVERIFY(closeRelative(plate.elasticity[i*5+j],expected));
        }
}

void TestStructuralPiezoMaterialOps::coupledPlaneStressAndReciprocityReference()
{
    const auto material = coupledReference();
    const auto plate = rStructuralPiezoReduceFreeStress(material,{0,1,5,4,3});
    // sigma_zz=2*Sxx+Syy+4*Szz-6*Ez=0, hence
    // Szz=1.5*Ez-0.5*Sxx-0.25*Syy, substituted into sigma and Dz.
    const std::array<double,25> expected{{
        9.0,-0.5,0.0,0.0,0.0,
        -0.5,7.75,0.0,0.0,0.0,
        0.0,0.0,1.0,0.0,0.0,
        0.0,0.0,0.0,2.0,0.0,
        0.0,0.0,0.0,0.0,3.0}};
    for (size_t i=0;i<25;++i) QVERIFY(closeRelative(plate.elasticity[i],expected[i]));
    const std::array<double,5> expectedCoupling{{0.0,0.5,0.0,-1.0,0.0}};
    for (size_t i=0;i<5;++i) QVERIFY(closeRelative(plate.coupling[i],expectedCoupling[i]));
    QVERIFY(closeRelative(plate.permittivity,11.0));
    QVERIFY(plate.permittivity >= material.permittivity[8]);

    const double strain = 0.2, voltage = -0.4;
    // Unit thickness/area: Ez=-V; force=K*S+e*V, Q=eps*V-e*S.
    const double force = 7.75*strain+0.5*voltage;
    const double charge = 11.0*voltage-0.5*strain;
    QVERIFY(closeRelative(plate.elasticity[6]*strain+plate.coupling[1]*voltage,force));
    QVERIFY(closeRelative(plate.permittivity*voltage-plate.coupling[1]*strain,charge));
    const double openVoltage = plate.coupling[1]*strain/plate.permittivity;
    const double openStiffness = 7.75+0.25/11.0;
    QVERIFY(closeRelative(7.75*strain+0.5*openVoltage,openStiffness*strain));
    QVERIFY(openStiffness > 7.75);
    auto reversed = material;
    for (double &e : reversed.coupling) e = -e;
    const auto opposite = rStructuralPiezoReduceFreeStress(reversed,{0,1,5,4,3});
    for (size_t i=0;i<5;++i) QVERIFY(closeRelative(opposite.coupling[i],-plate.coupling[i]));
    QVERIFY(closeRelative(opposite.permittivity,plate.permittivity));
    for (size_t i=0;i<25;++i) QVERIFY(closeRelative(opposite.elasticity[i],plate.elasticity[i]));
}

void TestStructuralPiezoMaterialOps::componentOrderingAndNoElimination()
{
    const auto material = coupledReference();
    const auto full = rStructuralPiezoReduceFreeStress(material,{5,3,1,0,4,2});
    const std::array<uint,6> order{{5,3,1,0,4,2}};
    for (size_t i=0;i<6;++i)
    {
        QCOMPARE(full.components[i],order[i]);
        QCOMPARE(full.coupling[i],material.coupling[12+order[i]]);
        for (size_t j=0;j<6;++j)
            QCOMPARE(full.elasticity[i*6+j],material.elasticity[order[i]*6+order[j]]);
    }
    QCOMPARE(full.permittivity,2.0);
    const auto ordered = rStructuralPiezoReduceFreeStress(material,{1,0});
    // Eliminating zz and xz leaves the known normal 2x2 block, independent
    // of which order the caller chooses for the retained coordinates.
    QVERIFY(closeRelative(ordered.elasticity[0],7.75));
    QVERIFY(closeRelative(ordered.elasticity[1],-0.5));
    QVERIFY(closeRelative(ordered.elasticity[2],-0.5));
    QVERIFY(closeRelative(ordered.elasticity[3],9.0));
    QVERIFY(closeRelative(ordered.permittivity,11.5));
}

void TestStructuralPiezoMaterialOps::legitimateSmallPhysicalScales()
{
    auto tiny = RStructuralPiezoContinuumMaterial::isotropic(1.0e-200,0.25,1.0e-200);
    const auto axial = rStructuralPiezoReduceFreeStress(tiny,{0});
    QVERIFY(closeRelative(axial.elasticity[0],1.0e-200));
    tiny.piezoelectric = true;
    tiny.permittivity[0] = tiny.permittivity[4] = tiny.permittivity[8] = 1.0e-250;
    tiny.validate();
    auto operators = operatorReference();
    operators.mass = {1.0e-200,0.0,0.0,1.0e-300};
    operators.stiffness = {1.0e-200,0.0,0.0,0.0};
    operators.capacitance = {1.0e-250};
    operators.totalMass = 1.0e-200;
    operators.validate();
}

void TestStructuralPiezoMaterialOps::rejectInvalidMaterial()
{
    auto m = coupledReference();
    m.validate();
    m.density = 0.0;
    QVERIFY_EXCEPTION_THROWN(m.validate(),std::invalid_argument);
    m = coupledReference(); m.density = -1.0;
    QVERIFY_EXCEPTION_THROWN(m.validate(),std::invalid_argument);
    m = coupledReference(); m.coupling[0] = std::numeric_limits<double>::quiet_NaN();
    QVERIFY_EXCEPTION_THROWN(m.validate(),std::invalid_argument);
    m = coupledReference(); m.elasticity[0] = std::numeric_limits<double>::infinity();
    QVERIFY_EXCEPTION_THROWN(m.validate(),std::invalid_argument);
    m = coupledReference(); m.elasticity[1] = 1.0;
    QVERIFY_EXCEPTION_THROWN(m.validate(),std::invalid_argument);
    m = coupledReference(); m.elasticity[2] = m.elasticity[12] = 10.0;
    QVERIFY_EXCEPTION_THROWN(m.validate(),std::invalid_argument);
    m = coupledReference(); m.permittivity[8] = 0.0;
    QVERIFY_EXCEPTION_THROWN(m.validate(),std::invalid_argument);
    m = coupledReference(); m.permittivity[1] = 0.5;
    QVERIFY_EXCEPTION_THROWN(m.validate(),std::invalid_argument);
    m = coupledReference(); m.permittivity[1] = m.permittivity[3] = 3.0;
    QVERIFY_EXCEPTION_THROWN(m.validate(),std::invalid_argument);
    m = coupledReference(); m.piezoelectric = false;
    QVERIFY_EXCEPTION_THROWN(m.validate(),std::invalid_argument);
    m = RStructuralPiezoContinuumMaterial::isotropic(10.0,0.2,1.0);
    m.permittivity[8] = -1.0;
    QVERIFY_EXCEPTION_THROWN(m.validate(),std::invalid_argument);
    m.permittivity.fill(0.0); m.permittivity[0] = std::numeric_limits<double>::quiet_NaN();
    QVERIFY_EXCEPTION_THROWN(m.validate(),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoContinuumMaterial::isotropic(10.0,0.5,1.0),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoContinuumMaterial::isotropic(10.0,-1.0,1.0),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoContinuumMaterial::isotropic(-10.0,0.2,1.0),std::invalid_argument);
}

void TestStructuralPiezoMaterialOps::rejectInvalidRetainedComponents()
{
    const auto m = coupledReference();
    QVERIFY_EXCEPTION_THROWN(rStructuralPiezoReduceFreeStress(m,{}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(rStructuralPiezoReduceFreeStress(m,{0,0}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(rStructuralPiezoReduceFreeStress(m,{6}),std::invalid_argument);
}

void TestStructuralPiezoMaterialOps::layerGeometryAndElectrodeValidation()
{
    RStructuralPiezoLayer base;
    base.material = coupledReference();
    base.lowerZ = -0.001; base.upperZ = 0.001;
    base.portId = QStringLiteral("patch");
    base.validate();
    auto layer = base; layer.polarity = -1; layer.validate();
    layer = base; layer.polarity = 0;
    QVERIFY_EXCEPTION_THROWN(layer.validate(),std::invalid_argument);
    layer = base; layer.upperZ = layer.lowerZ;
    QVERIFY_EXCEPTION_THROWN(layer.validate(),std::invalid_argument);
    layer = base; layer.lowerZ = std::numeric_limits<double>::quiet_NaN();
    QVERIFY_EXCEPTION_THROWN(layer.validate(),std::invalid_argument);
    layer = base; layer.width = -1.0;
    QVERIFY_EXCEPTION_THROWN(layer.validate(),std::invalid_argument);
    layer = base; layer.coverage[0] = -1.01;
    QVERIFY_EXCEPTION_THROWN(layer.validate(),std::invalid_argument);
    layer = base; layer.coverage[3] = layer.coverage[2];
    QVERIFY_EXCEPTION_THROWN(layer.validate(),std::invalid_argument);
    layer = base; layer.coverage[1] = std::numeric_limits<double>::infinity();
    QVERIFY_EXCEPTION_THROWN(layer.validate(),std::invalid_argument);
    layer = base; layer.portId.clear();
    QVERIFY_EXCEPTION_THROWN(layer.validate(),std::invalid_argument);
    layer.portId = QStringLiteral("  ");
    QVERIFY_EXCEPTION_THROWN(layer.validate(),std::invalid_argument);
    layer.portId = QStringLiteral("bad")+QChar(0);
    QVERIFY_EXCEPTION_THROWN(layer.validate(),std::invalid_argument);
    layer = base; layer.material = RStructuralPiezoContinuumMaterial::isotropic(10.0,0.2,1.0);
    QVERIFY_EXCEPTION_THROWN(layer.validate(),std::invalid_argument);
    layer.portId.clear(); layer.validate();
}

void TestStructuralPiezoMaterialOps::elementOperatorValidation()
{
    auto base = operatorReference(); base.validate();
    auto op = base; op.mass.pop_back();
    QVERIFY_EXCEPTION_THROWN(op.validate(),std::invalid_argument);
    op = base; op.stiffness.pop_back();
    QVERIFY_EXCEPTION_THROWN(op.validate(),std::invalid_argument);
    op = base; op.coupling.pop_back();
    QVERIFY_EXCEPTION_THROWN(op.validate(),std::invalid_argument);
    op = base; op.capacitance.clear();
    QVERIFY_EXCEPTION_THROWN(op.validate(),std::invalid_argument);
    op = base; op.dofs = 0;
    QVERIFY_EXCEPTION_THROWN(op.validate(),std::invalid_argument);
    op = base; op.totalMass = -1.0;
    QVERIFY_EXCEPTION_THROWN(op.validate(),std::invalid_argument);
    op = base; op.mass[1] = 0.0;
    QVERIFY_EXCEPTION_THROWN(op.validate(),std::invalid_argument);
    op = base; op.mass = {1.0,2.0,2.0,1.0};
    QVERIFY_EXCEPTION_THROWN(op.validate(),std::invalid_argument);
    op = base; op.stiffness[0] = -1.0;
    QVERIFY_EXCEPTION_THROWN(op.validate(),std::invalid_argument);
    op = base; op.stiffness[1] = 0.0;
    QVERIFY_EXCEPTION_THROWN(op.validate(),std::invalid_argument);
    op = base; op.coupling[0] = std::numeric_limits<double>::quiet_NaN();
    QVERIFY_EXCEPTION_THROWN(op.validate(),std::invalid_argument);
    op = base; op.capacitance[0] = -1.0e-10;
    QVERIFY_EXCEPTION_THROWN(op.validate(),std::invalid_argument);
    op = base; op.portIds[0] = QStringLiteral("bad")+QChar(0);
    QVERIFY_EXCEPTION_THROWN(op.validate(),std::invalid_argument);
    op = base; op.portIds.push_back(op.portIds[0]);
    op.coupling = {2.0,1.0,-2.0,-1.0}; op.capacitance = {2.0,0.0,0.0,1.0};
    QVERIFY_EXCEPTION_THROWN(op.validate(),std::invalid_argument);
    op.portIds[1] = QStringLiteral("other"); op.validate();
    op.capacitance = {1.0,2.0,2.0,1.0};
    QVERIFY_EXCEPTION_THROWN(op.validate(),std::invalid_argument);
    op = base; op.portIds.clear(); op.coupling.clear(); op.capacitance.clear(); op.validate();
}

void TestStructuralPiezoMaterialOps::strainChargeConversionAxesAndReplacement()
{
    std::array<double,36> s{}; for(size_t i=0;i<6;++i) s[6*i+i]=.1;
    s[1]=s[6]=-.02;
    std::array<double,18> d{}; d[0]=.2; d[1]=-.06; d[12]=-.03; d[13]=.08; d[14]=.25;
    const std::array<double,9> epsT{{10,0,0,0,10,0,0,0,10}};
    const auto m=RStructuralPiezoContinuumMaterial::fromStrainCharge(s,d,epsT,7);
    const auto axial=rStructuralPiezoReduceFreeStress(m,{0},0);
    QVERIFY(closeRelative(axial.elasticity[0],10)); QVERIFY(closeRelative(axial.coupling[0],2));
    QVERIFY(closeRelative(axial.permittivity,9.6));
    const auto transverse=rStructuralPiezoReduceFreeStress(m,{0},2);
    QVERIFY(closeRelative(transverse.coupling[0],-.3)); QVERIFY(closeRelative(transverse.permittivity,9.991));
    // d33 must be oriented into local x/x; z/x is transverse field (d31).
    QVERIFY(!closeRelative(axial.coupling[0],transverse.coupling[0]));
    QCOMPARE(rStructuralPiezoReduceFreeStress(m,{0}).coupling,transverse.coupling);
    const auto member=rStructuralPiezoCompileAxialMember(m,2,3);
    QVERIFY(closeRelative(member.k,20./3)); QVERIFY(closeRelative(member.g,4./3));
    QVERIFY(closeRelative(member.Ce,19.2/3)); QCOMPARE(member.fixedMass,42.);
    // One replacement mass; the compiler does not add the original host mass.
    QCOMPARE(member.fixedMass,7.*2.*3.);
    s[1]=s[6]=0; d.fill(0); d[0]=.2;
    const auto one=RStructuralPiezoContinuumMaterial::fromStrainCharge(s,d,epsT,7);
    QVERIFY(closeRelative(one.elasticity[0],10)); QVERIFY(closeRelative(one.coupling[0],2));
    QVERIFY(closeRelative(one.permittivity[0],9.6));
}

void TestStructuralPiezoMaterialOps::anisotropicStrainStressEquivalence()
{
    std::array<double,36> lower{},s{};
    for(size_t i=0;i<6;++i) for(size_t j=0;j<=i;++j) lower[6*i+j]=i==j?.3+.04*i:.01*(i+j+1);
    for(size_t i=0;i<6;++i) for(size_t j=0;j<6;++j) for(size_t k=0;k<6;++k) s[6*i+j]+=lower[6*i+k]*lower[6*j+k];
    std::array<double,18> d{}; for(size_t i=0;i<18;++i) d[i]=.003*(double(i)-7);
    const std::array<double,9> epsT{{4,.2,-.1,.2,5,.3,-.1,.3,6}};
    const auto m=RStructuralPiezoContinuumMaterial::fromStrainCharge(s,d,epsT,7500);
    const std::array<double,6> strain{{.07,-.04,.02,.03,-.01,.05}};
    const std::array<double,3> electric{{.8,-.6,.3}};
    std::array<double,6> sigma{};
    for(size_t i=0;i<6;++i)
    {
        for(size_t j=0;j<6;++j) sigma[i]+=m.elasticity[6*i+j]*strain[j];
        for(size_t j=0;j<3;++j) sigma[i]-=m.coupling[6*j+i]*electric[j];
    }
    for(size_t i=0;i<6;++i)
    {
        double reconstructed=0; for(size_t j=0;j<6;++j) reconstructed+=s[6*i+j]*sigma[j];
        for(size_t j=0;j<3;++j) reconstructed+=d[6*j+i]*electric[j];
        QVERIFY(closeRelative(reconstructed,strain[i],2e-12));
    }
    for(size_t i=0;i<3;++i)
    {
        double stressForm=0,strainForm=0;
        for(size_t j=0;j<6;++j) {stressForm+=m.coupling[6*i+j]*strain[j]; strainForm+=d[6*i+j]*sigma[j];}
        for(size_t j=0;j<3;++j) {stressForm+=m.permittivity[3*i+j]*electric[j]; strainForm+=epsT[3*i+j]*electric[j];}
        QVERIFY(closeRelative(stressForm,strainForm,2e-12));
        const auto axial=rStructuralPiezoReduceFreeStress(m,{0},uint(i));
        QVERIFY(closeRelative(axial.elasticity[0],1/s[0]));
        QVERIFY(closeRelative(axial.coupling[0],d[6*i]/s[0]));
        QVERIFY(closeRelative(axial.permittivity,epsT[4*i]-d[6*i]*d[6*i]/s[0]));
    }
}

void TestStructuralPiezoMaterialOps::axialFiniteScales()
{
    // Intermediate area*c overflows, while all final SI properties are finite.
    auto huge=RStructuralPiezoContinuumMaterial::isotropic(1e200,0,1e-300);
    huge.piezoelectric=true; huge.permittivity={{1e200,0,0,0,1e200,0,0,0,1e200}};
    const auto robust=rStructuralPiezoCompileAxialMember(huge,1e200,1e200);
    QVERIFY(closeRelative(robust.k,1e200)); QVERIFY(closeRelative(robust.Ce,1e200));
    QVERIFY(closeRelative(robust.fixedMass,1e100)); QCOMPARE(robust.g,0.);
    auto tiny=RStructuralPiezoContinuumMaterial::isotropic(1e-200,0,1e-200);
    tiny.piezoelectric=true; tiny.permittivity={{1e-250,0,0,0,1e-250,0,0,0,1e-250}};
    const auto small=rStructuralPiezoCompileAxialMember(tiny,1,1);
    QVERIFY(closeRelative(small.k,1e-200)); QVERIFY(closeRelative(small.Ce,1e-250));
}

void TestStructuralPiezoMaterialOps::rejectCompilerInputs()
{
    std::array<double,36> s{}; for(size_t i=0;i<6;++i) s[6*i+i]=.1;
    std::array<double,18> d{}; d[0]=.2;
    std::array<double,9> epsT{{10,0,0,0,10,0,0,0,10}};
    const auto convert=[&](const std::array<double,36> &ss,const std::array<double,18> &dd,const std::array<double,9> &ee,double rho=7.) {return RStructuralPiezoContinuumMaterial::fromStrainCharge(ss,dd,ee,rho);};
    auto badS=s; badS[0]=-1; QVERIFY_EXCEPTION_THROWN(convert(badS,d,epsT),std::invalid_argument);
    badS=s; badS[1]=badS[6]=1; QVERIFY_EXCEPTION_THROWN(convert(badS,d,epsT),std::invalid_argument);
    badS=s; badS[1]=.01; QVERIFY_EXCEPTION_THROWN(convert(badS,d,epsT),std::invalid_argument);
    badS=s; badS[0]=std::numeric_limits<double>::quiet_NaN(); QVERIFY_EXCEPTION_THROWN(convert(badS,d,epsT),std::invalid_argument);
    auto badD=d; badD[0]=std::numeric_limits<double>::infinity(); QVERIFY_EXCEPTION_THROWN(convert(s,badD,epsT),std::invalid_argument);
    auto badE=epsT; badE[0]=.1; // epsilonT SPD, corrected epsilonS indefinite.
    QVERIFY_EXCEPTION_THROWN(convert(s,d,badE),std::invalid_argument);
    badE=epsT; badE[0]=0; QVERIFY_EXCEPTION_THROWN(convert(s,d,badE),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(convert(s,d,epsT,0),std::invalid_argument);
    const auto m=convert(s,d,epsT);
    QVERIFY_EXCEPTION_THROWN(rStructuralPiezoReduceFreeStress(m,{0},3),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(rStructuralPiezoReduceFreeStress(m,{0},uint(-1)),std::invalid_argument);
    for(double invalid:{0.,-1.,std::numeric_limits<double>::infinity(),std::numeric_limits<double>::quiet_NaN()})
    {
        QVERIFY_EXCEPTION_THROWN(rStructuralPiezoCompileAxialMember(m,invalid,1),std::invalid_argument);
        QVERIFY_EXCEPTION_THROWN(rStructuralPiezoCompileAxialMember(m,1,invalid),std::invalid_argument);
    }
    QVERIFY_EXCEPTION_THROWN(rStructuralPiezoCompileAxialMember(m,1e308,1e-308),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(rStructuralPiezoCompileAxialMember(m,1e-308,1e308),std::invalid_argument);
    const auto passive=RStructuralPiezoContinuumMaterial::isotropic(10,0,7);
    QVERIFY_EXCEPTION_THROWN(rStructuralPiezoCompileAxialMember(passive,1,1),std::invalid_argument);
}

QTEST_APPLESS_MAIN(TestStructuralPiezoMaterialOps)
#include "tst_structural_piezo_material_ops.moc"

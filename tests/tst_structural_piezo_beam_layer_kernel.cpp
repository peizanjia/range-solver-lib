#include <QtTest>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "rstructural_piezo_element_operators.h"

namespace
{
using Vector = std::array<double,6>;
using Polynomial = std::array<long double,4>;

RStructuralPiezoLayer layer(double lower,double upper,bool active=false)
{
    RStructuralPiezoLayer result;
    result.material = RStructuralPiezoContinuumMaterial::isotropic(70.0e9,0.0,2700.0);
    result.lowerZ = lower;
    result.upperZ = upper;
    result.width = 0.03;
    if (active)
    {
        result.material.piezoelectric = true;
        result.material.coupling[12] = -5.2; // e_z,xx, in element axes.
        result.material.permittivity[0] = 1.3e-8;
        result.material.permittivity[4] = 1.3e-8;
        result.material.permittivity[8] = 1.3e-8;
        result.portId = QStringLiteral("electrode");
    }
    return result;
}

bool close(double actual,double expected,double relative=3.0e-12,double absolute=1.0e-12)
{
    return std::isfinite(actual) && std::isfinite(expected) &&
           std::abs(actual-expected) <= absolute+relative*std::abs(expected);
}

Vector multiply(const std::vector<double> &matrix,const Vector &q)
{
    Vector result{};
    for (size_t i=0;i<6;++i)
        for (size_t j=0;j<6;++j) result[i] += matrix[i*6+j]*q[j];
    return result;
}

double quadratic(const std::vector<double> &matrix,const Vector &q)
{
    const Vector product = multiply(matrix,q);
    double result = 0.0;
    for (size_t i=0;i<6;++i) result += q[i]*product[i];
    return result;
}

double dotCoupling(const RStructuralPiezoElementOperators &operators,const Vector &q)
{
    double result = 0.0;
    for (size_t i=0;i<6;++i) result += operators.coupling[i]*q[i];
    return result;
}

// Independent power-polynomial antiderivatives, rather than quadrature, for
// reference integrals on a partially covered element.
long double productIntegral(const Polynomial &a,const Polynomial &b,long double lo,long double hi)
{
    long double result = 0.0L;
    for (size_t i=0;i<4;++i)
        for (size_t j=0;j<4;++j)
            result += a[i]*b[j]*(std::pow(hi,int(i+j+1))-std::pow(lo,int(i+j+1)))/
                      static_cast<long double>(i+j+1);
    return result;
}

std::vector<double> analyticMass(double length,const RStructuralPiezoLayer &item)
{
    const long double l = length;
    const std::array<Polynomial,6> axial{{{{1,-1,0,0}},{{0,0,0,0}},{{0,0,0,0}},
                                         {{0,1,0,0}},{{0,0,0,0}},{{0,0,0,0}}}};
    const std::array<Polynomial,6> transverse{{{{0,0,0,0}},{{1,0,-3,2}},{{0,-l,2*l,-l}},
                                              {{0,0,0,0}},{{0,0,3,-2}},{{0,0,l,-l}}}};
    const std::array<Polynomial,6> rotation{{{{0,0,0,0}},{{0,6/l,-6/l,0}},{{1,-4,3,0}},
                                            {{0,0,0,0}},{{0,-6/l,6/l,0}},{{0,-2,3,0}}}};
    const long double lo = (item.coverage[0]+1.0L)/2.0L;
    const long double hi = (item.coverage[1]+1.0L)/2.0L;
    const long double t = item.upperZ-item.lowerZ;
    const long double c = (item.upperZ+item.lowerZ)/2.0L;
    const long double i0 = item.width*t,i1 = i0*c,i2 = i0*(c*c+t*t/12.0L);
    std::vector<double> result(36,0.0);
    for (size_t i=0;i<6;++i)
        for (size_t j=0;j<6;++j)
            result[i*6+j] = static_cast<double>(item.material.density*l*(
                i0*(productIntegral(axial[i],axial[j],lo,hi)+
                    productIntegral(transverse[i],transverse[j],lo,hi))+
                i1*(productIntegral(axial[i],rotation[j],lo,hi)+
                    productIntegral(rotation[i],axial[j],lo,hi))+
                i2*productIntegral(rotation[i],rotation[j],lo,hi)));
    return result;
}

Vector integratedStrain(double length,double center,double xiLower,double xiUpper)
{
    const double a = (xiLower+1.0)/2.0,b = (xiUpper+1.0)/2.0;
    // Integral B dx from endpoint differences in the displacement derivatives.
    return {{a-b,center*(6.0*(b-a)-6.0*(b*b-a*a))/length,
             center*(-4.0*(b-a)+3.0*(b*b-a*a)),b-a,
             center*(-6.0*(b-a)+6.0*(b*b-a*a))/length,
             center*(-2.0*(b-a)+3.0*(b*b-a*a))}};
}
}

class TestStructuralPiezoBeamLayerKernel : public QObject
{
    Q_OBJECT
private slots:
    void homogeneousClosedStiffnessAndConsistentMass();
    void rigidMotionMassCentroidAndZeroStrain();
    void symmetricAndEccentricLayers();
    void patchOwnSecondMoment();
    void electrodePolarityPortOrderAndParallelCapacitance();
    void partialCoverageExactMassForceAndCharge();
    void enthalpyDerivativesAndOpenShortCircuit();
    void passiveAndZeroCouplingLimits();
    void cantileverStaticCompliance();
    void tinyRepresentableGeometry();
    void representableOperatorsAfterVolumeUnderflow();
    void rejectsInvalidGeometryMaterialAndOverlap();
};

void TestStructuralPiezoBeamLayerKernel::homogeneousClosedStiffnessAndConsistentMass()
{
    const double l = 0.8,t = 0.012;
    const auto item = layer(-t/2.0,t/2.0);
    const auto result = RStructuralPiezoBeamLayerKernel::assemble(l,{item});
    QCOMPARE(result.dofs,6u);
    QVERIFY(result.portIds.empty());
    QVERIFY(result.coupling.empty());
    QVERIFY(result.capacitance.empty());
    const double area = item.width*t,inertia = item.width*t*t*t/12.0;
    const double e = 70.0e9,rho = 2700.0;
    std::vector<double> k(36,0.0),m(36,0.0);
    k[0] = k[21] = e*area/l;
    k[3] = k[18] = -e*area/l;
    m[0] = m[21] = rho*area*l/3.0;
    m[3] = m[18] = rho*area*l/6.0;
    const std::array<size_t,4> bend{{1,2,4,5}};
    const std::array<double,16> kb{{12,-6*l,-12,-6*l,
        -6*l,4*l*l,6*l,2*l*l,-12,6*l,12,6*l,-6*l,2*l*l,6*l,4*l*l}};
    const std::array<double,16> mt{{156,-22*l,54,13*l,
        -22*l,4*l*l,-13*l,-3*l*l,54,-13*l,156,22*l,13*l,-3*l*l,22*l,4*l*l}};
    const std::array<double,16> mr{{36,-3*l,-36,-3*l,
        -3*l,4*l*l,3*l,-l*l,-36,3*l,36,3*l,-3*l,-l*l,3*l,4*l*l}};
    for (size_t i=0;i<4;++i)
        for (size_t j=0;j<4;++j)
        {
            k[bend[i]*6+bend[j]] = e*inertia/(l*l*l)*kb[i*4+j];
            m[bend[i]*6+bend[j]] = rho*area*l/420.0*mt[i*4+j]+
                                  rho*inertia/(30.0*l)*mr[i*4+j];
        }
    for (size_t i=0;i<36;++i)
    {
        QVERIFY2(close(result.stiffness[i],k[i],3.0e-12,1.0e-9),qPrintable(QString::number(i)));
        QVERIFY2(close(result.mass[i],m[i]),qPrintable(QString::number(i)));
    }
    QVERIFY(close(result.totalMass,rho*area*l));
    // All six supported DOFs carry physical inertia; no unrelated torsion DOF.
    for (size_t i=0;i<6;++i) QVERIFY(result.mass[i*6+i] > 0.0);
}

void TestStructuralPiezoBeamLayerKernel::rigidMotionMassCentroidAndZeroStrain()
{
    const double l = 0.7;
    const auto item = layer(0.015,0.025);
    const auto result = RStructuralPiezoBeamLayerKernel::assemble(l,{item});
    const double u = 0.31,w = -0.23,theta = 0.17;
    const Vector motion{{u,w,theta,u,w-l*theta,theta}};
    const double t = item.upperZ-item.lowerZ,z = (item.upperZ+item.lowerZ)/2.0;
    const double expected = result.totalMass*(u*u+2*z*u*theta+
        (z*z+t*t/12.0)*theta*theta+w*w-l*w*theta+l*l*theta*theta/3.0);
    QVERIFY(close(quadratic(result.mass,motion),expected));
    const std::array<Vector,3> rigid{{{{1,0,0,1,0,0}},{{0,1,0,0,1,0}},{{0,0,1,0,-l,1}}}};
    const double scale = *std::max_element(result.stiffness.begin(),result.stiffness.end());
    for (const auto &q : rigid)
        for (double force : multiply(result.stiffness,q)) QVERIFY(std::abs(force) < scale*3.0e-14);
    const auto reference = analyticMass(l,item);
    for (size_t i=0;i<36;++i) QVERIFY(close(result.mass[i],reference[i]));
}

void TestStructuralPiezoBeamLayerKernel::symmetricAndEccentricLayers()
{
    const double l = 0.9;
    const auto bottom = layer(-0.012,-0.006),top = layer(0.006,0.012);
    const auto symmetric = RStructuralPiezoBeamLayerKernel::assemble(l,{bottom,top});
    for (size_t axial : {size_t(0),size_t(3)})
        for (size_t bend : {size_t(1),size_t(2),size_t(4),size_t(5)})
        {
            QVERIFY(std::abs(symmetric.mass[axial*6+bend]) < 1.0e-14);
            QVERIFY(std::abs(symmetric.stiffness[axial*6+bend]) < 1.0e-9);
        }
    const auto offset = RStructuralPiezoBeamLayerKernel::assemble(l,{top});
    const double i1 = top.width*(top.upperZ-top.lowerZ)*(top.upperZ+top.lowerZ)/2.0;
    const std::array<size_t,4> bend{{1,2,4,5}};
    const std::array<double,8> factors{{0.5,l/12.0,-0.5,-l/12.0,
                                       0.5,-l/12.0,-0.5,l/12.0}};
    for (size_t a=0;a<2;++a)
        for (size_t j=0;j<4;++j)
            QVERIFY(close(offset.mass[a*3*6+bend[j]],2700.0*i1*factors[a*4+j]));
    QVERIFY(close(offset.stiffness[2],70.0e9*i1/l,3.0e-12,1.0e-8));
    QVERIFY(close(offset.stiffness[5],-70.0e9*i1/l,3.0e-12,1.0e-8));
    QVERIFY(offset.mass[1] != 0.0);
    QVERIFY(offset.stiffness[2] != 0.0);
}

void TestStructuralPiezoBeamLayerKernel::patchOwnSecondMoment()
{
    const double l = 0.6;
    const auto host = layer(-0.005,0.005);
    const auto patch = layer(0.005,0.013,true);
    const auto base = RStructuralPiezoBeamLayerKernel::assemble(l,{host});
    const auto both = RStructuralPiezoBeamLayerKernel::assemble(l,{host,patch});
    const double t = patch.upperZ-patch.lowerZ,z = (patch.upperZ+patch.lowerZ)/2.0;
    const double inertia = patch.width*t*(z*z+t*t/12.0);
    QVERIFY(close(both.stiffness[7]-base.stiffness[7],12.0*70.0e9*inertia/(l*l*l),3.0e-12,1.0e-8));
    const double centerOnly = 12.0*70.0e9*patch.width*t*z*z/(l*l*l);
    QVERIFY((both.stiffness[7]-base.stiffness[7])-centerOnly > 0.01*centerOnly);
    const auto patchMass = analyticMass(l,patch);
    for (size_t i=0;i<36;++i) QVERIFY(close(both.mass[i]-base.mass[i],patchMass[i]));
    QVERIFY(close(both.totalMass-base.totalMass,2700.0*patch.width*t*l));
}

void TestStructuralPiezoBeamLayerKernel::electrodePolarityPortOrderAndParallelCapacitance()
{
    const double l = 0.75;
    auto top = layer(0.004,0.006,true);
    auto bottom = layer(-0.006,-0.004,true);
    const auto positive = RStructuralPiezoBeamLayerKernel::assemble(l,{top});
    auto reversed = top;
    reversed.polarity = -1;
    const auto negative = RStructuralPiezoBeamLayerKernel::assemble(l,{reversed});
    for (size_t i=0;i<6;++i) QVERIFY(close(positive.coupling[i],-negative.coupling[i]));
    QCOMPARE(positive.capacitance,negative.capacitance);
    const auto lower = RStructuralPiezoBeamLayerKernel::assemble(l,{bottom});
    QVERIFY(close(positive.coupling[0],lower.coupling[0]));
    QVERIFY(close(positive.coupling[2],-lower.coupling[2]));
    QVERIFY(close(positive.capacitance[0],lower.capacitance[0]));
    const auto parallel = RStructuralPiezoBeamLayerKernel::assemble(l,{top,bottom});
    QCOMPARE(parallel.portIds.size(),size_t(1));
    QVERIFY(close(parallel.capacitance[0],2.0*positive.capacitance[0]));
    for (size_t i=0;i<6;++i) QVERIFY(close(parallel.coupling[i],positive.coupling[i]+lower.coupling[i]));
    QVERIFY(std::abs(parallel.coupling[2]) < 1.0e-15);
    top.portId = QStringLiteral("z_first");
    bottom.portId = QStringLiteral("a_second");
    const auto separate = RStructuralPiezoBeamLayerKernel::assemble(l,{top,bottom});
    QCOMPARE(separate.portIds[0],top.portId);
    QCOMPARE(separate.portIds[1],bottom.portId);
    QCOMPARE(separate.capacitance[1],0.0);
    QCOMPARE(separate.capacitance[2],0.0);
    QVERIFY(close(separate.capacitance[0],positive.capacitance[0]));
    QVERIFY(close(separate.capacitance[3],lower.capacitance[0]));
}

void TestStructuralPiezoBeamLayerKernel::partialCoverageExactMassForceAndCharge()
{
    const double l = 1.1;
    const auto host = layer(-0.006,0.006);
    auto patch = layer(0.006,0.009,true);
    patch.coverage = {{-0.6,0.4,-1.0,1.0}};
    const auto result = RStructuralPiezoBeamLayerKernel::assemble(l,{host,patch});
    const auto hostMass = analyticMass(l,host),patchMass = analyticMass(l,patch);
    for (size_t i=0;i<36;++i) QVERIFY(close(result.mass[i],hostMass[i]+patchMass[i]));
    const double z = (patch.lowerZ+patch.upperZ)/2.0;
    const Vector integrated = integratedStrain(l,z,patch.coverage[0],patch.coverage[1]);
    for (size_t i=0;i<6;++i)
        QVERIFY(close(result.coupling[i],patch.width*patch.material.coupling[12]*integrated[i],3.0e-12,1.0e-15));
    const double coveredLength = l*(patch.coverage[1]-patch.coverage[0])/2.0;
    const double t = patch.upperZ-patch.lowerZ;
    QVERIFY(close(result.capacitance[0],1.3e-8*patch.width*coveredLength/t,3.0e-12,1.0e-18));
    QVERIFY(close(result.totalMass,2700.0*(host.width*(host.upperZ-host.lowerZ)*l+
                  patch.width*t*coveredLength)));
    // A linear axial strain plus constant curvature is represented exactly.
    const double axialStrain = 2.0e-5,curvature = -3.0e-4;
    const Vector q{{0,0,0,axialStrain*l,-curvature*l*l/2.0,curvature*l}};
    const double expectedCharge = -patch.material.coupling[12]*patch.width*coveredLength*
                                  (axialStrain+z*curvature);
    QVERIFY(close(-dotCoupling(result,q),expectedCharge,3.0e-12,1.0e-16));
    Vector expectedMechanicalForce{};
    for (const auto &item : {host,patch})
    {
        const double thickness = item.upperZ-item.lowerZ;
        const double center = (item.upperZ+item.lowerZ)/2.0;
        const double i0 = item.width*thickness,i1 = i0*center;
        const double i2 = i0*(center*center+thickness*thickness/12.0);
        const Vector axialIntegral = integratedStrain(l,0.0,item.coverage[0],item.coverage[1]);
        const Vector combinedIntegral = integratedStrain(l,1.0,item.coverage[0],item.coverage[1]);
        for (size_t i=0;i<6;++i)
            expectedMechanicalForce[i] += 70.0e9*(
                axialIntegral[i]*(i0*axialStrain+i1*curvature)+
                (combinedIntegral[i]-axialIntegral[i])*(i1*axialStrain+i2*curvature));
    }
    const Vector mechanicalForce = multiply(result.stiffness,q);
    for (size_t i=0;i<6;++i)
        QVERIFY(close(mechanicalForce[i],expectedMechanicalForce[i],3.0e-11,1.0e-8));
    const double voltage = 2.3;
    for (size_t i=0;i<6;++i)
        QVERIFY(close(result.coupling[i]*voltage,patch.width*patch.material.coupling[12]*
                      integrated[i]*voltage,3.0e-12,1.0e-15));
    // Separately meshed electrode pieces share the same voltage/charge port.
    auto left = patch,right = patch;
    left.coverage[1] = -0.1;
    right.coverage[0] = -0.1;
    const auto split = RStructuralPiezoBeamLayerKernel::assemble(l,{host,left,right});
    for (size_t i=0;i<36;++i)
    {
        QVERIFY(close(split.mass[i],result.mass[i]));
        QVERIFY(close(split.stiffness[i],result.stiffness[i],3.0e-12,1.0e-8));
    }
    for (size_t i=0;i<6;++i) QVERIFY(close(split.coupling[i],result.coupling[i],3.0e-12,1.0e-15));
    QVERIFY(close(split.capacitance[0],result.capacitance[0],3.0e-12,1.0e-18));
}

void TestStructuralPiezoBeamLayerKernel::enthalpyDerivativesAndOpenShortCircuit()
{
    const auto result = RStructuralPiezoBeamLayerKernel::assemble(0.8,
        {layer(-0.004,0.004),layer(0.004,0.006,true)});
    const Vector q{{1.0e-5,2.0e-5,-3.0e-5,-2.0e-5,4.0e-5,5.0e-5}};
    const double voltage = 2.7;
    const auto enthalpy = [&](const Vector &displacement,double v) {
        return 0.5*quadratic(result.stiffness,displacement)+dotCoupling(result,displacement)*v-
               0.5*result.capacitance[0]*v*v;
    };
    const auto force = [&](const Vector &displacement,double v) {
        Vector f = multiply(result.stiffness,displacement);
        for (size_t i=0;i<6;++i) f[i] += result.coupling[i]*v;
        return f;
    };
    const auto charge = [&](const Vector &displacement,double v) {
        return result.capacitance[0]*v-dotCoupling(result,displacement);
    };
    const double hq = 1.0e-7,hv = 0.01;
    const Vector baseForce = force(q,voltage),fp = force(q,voltage+hv),fm = force(q,voltage-hv);
    for (size_t i=0;i<6;++i)
    {
        auto qp = q,qm = q;
        qp[i] += hq;
        qm[i] -= hq;
        const double gradient = (enthalpy(qp,voltage)-enthalpy(qm,voltage))/(2.0*hq);
        QVERIFY(close(gradient,baseForce[i],1.0e-7,1.0e-7));
        QVERIFY(close((fp[i]-fm[i])/(2.0*hv),result.coupling[i],1.0e-7,1.0e-10));
        QVERIFY(close((charge(qp,voltage)-charge(qm,voltage))/(2.0*hq),-result.coupling[i],1.0e-8,1.0e-11));
    }
    QVERIFY(close((enthalpy(q,voltage+hv)-enthalpy(q,voltage-hv))/(2.0*hv),
                  -charge(q,voltage),1.0e-7,1.0e-13));
    const double openVoltage = dotCoupling(result,q)/result.capacitance[0];
    QVERIFY(std::abs(charge(q,openVoltage)) < 1.0e-20);
    auto openStiffness = result.stiffness;
    std::vector<double> difference(36,0.0);
    for (size_t i=0;i<6;++i)
        for (size_t j=0;j<6;++j)
        {
            difference[i*6+j] = result.coupling[i]*result.coupling[j]/result.capacitance[0];
            openStiffness[i*6+j] += difference[i*6+j];
        }
    const auto openForce = force(q,openVoltage),condensed = multiply(openStiffness,q);
    for (size_t i=0;i<6;++i) QVERIFY(close(openForce[i],condensed[i],3.0e-12,1.0e-8));
    QVERIFY(close(quadratic(difference,q),std::pow(dotCoupling(result,q),2)/result.capacitance[0]));
    // Deterministic vectors probe the rank-one, nonnegative open/short difference.
    for (int k=1;k<=12;++k)
    {
        Vector probe{};
        for (size_t i=0;i<6;++i) probe[i] = std::sin(double(k)*double(i+1));
        QVERIFY(quadratic(difference,probe) >= -1.0e-10);
    }
}

void TestStructuralPiezoBeamLayerKernel::passiveAndZeroCouplingLimits()
{
    auto active = layer(-0.006,0.006,true);
    active.material.coupling.fill(0.0);
    auto passive = active;
    passive.material.piezoelectric = false;
    passive.material.permittivity.fill(0.0);
    passive.portId.clear();
    const auto a = RStructuralPiezoBeamLayerKernel::assemble(0.8,{active});
    const auto p = RStructuralPiezoBeamLayerKernel::assemble(0.8,{passive});
    QCOMPARE(a.mass,p.mass);
    QCOMPARE(a.stiffness,p.stiffness);
    for (double g : a.coupling) QCOMPARE(g,0.0);
    QVERIFY(a.capacitance[0] > 0.0);
    QVERIFY(p.coupling.empty());
    QVERIFY(p.capacitance.empty());
    // Nonzero Poisson ratio must use uniaxial free-stress E, not c_xxxx.
    passive.material = RStructuralPiezoContinuumMaterial::isotropic(70.0e9,0.3,2700.0);
    const auto freeStress = RStructuralPiezoBeamLayerKernel::assemble(0.8,{passive});
    for (size_t i=0;i<36;++i) QVERIFY(close(freeStress.stiffness[i],p.stiffness[i],3.0e-12,1.0e-9));
}

void TestStructuralPiezoBeamLayerKernel::cantileverStaticCompliance()
{
    const double l = 0.8,t = 0.012,force = 1.2;
    const auto item = layer(-t/2,t/2);
    const auto result = RStructuralPiezoBeamLayerKernel::assemble(l,{item});
    // Clamp node one, solve the exact supported tip w/theta block.
    const double a = result.stiffness[28],b = result.stiffness[29],d = result.stiffness[35];
    const double determinant = a*d-b*b;
    const double w = d*force/determinant,theta = -b*force/determinant;
    const double ei = 70.0e9*item.width*t*t*t/12.0;
    QVERIFY(close(w,force*l*l*l/(3.0*ei)));
    QVERIFY(close(theta,-force*l*l/(2.0*ei)));
}

void TestStructuralPiezoBeamLayerKernel::tinyRepresentableGeometry()
{
    auto tiny = layer(-5.0e-61,5.0e-61,true);
    tiny.width = 1.0e-60;
    RStructuralPiezoElementOperators result;
    try
    {
        result = RStructuralPiezoBeamLayerKernel::assemble(1.0e-60,{tiny});
    }
    catch (const std::exception &exception)
    {
        QFAIL(qPrintable(QStringLiteral("Tiny beam assembly: %1").arg(exception.what())));
    }
    QVERIFY(close(result.totalMass,2700.0e-180,3.0e-12,0.0));
    QVERIFY(result.stiffness[0] > 0.0);
    QVERIFY(result.mass[35] > 0.0);
    QVERIFY(close(result.capacitance[0],1.3e-68,3.0e-12,0.0));
    for (double value : result.mass) QVERIFY(std::isfinite(value));
    // Offset squared overflows in double, while the volume-weighted inertia
    // remains representable. The quadrature must retain that inertia.
    auto offset = layer(1.0e155,2.0e155);
    offset.width = 1.0e-300;
    RStructuralPiezoElementOperators largeOffset;
    try
    {
        // Resolve translation and rotation together; a unit-length member
        // loses its rigid-translation inertia at this enormous eccentricity.
        offset.material = RStructuralPiezoContinuumMaterial::isotropic(1.0e100,0.0,1.0e-100);
        largeOffset = RStructuralPiezoBeamLayerKernel::assemble(1.0e155,{offset});
    }
    catch (const std::exception &exception)
    {
        QFAIL(qPrintable(QStringLiteral("Large-offset beam assembly: %1").arg(exception.what())));
    }
    QVERIFY(largeOffset.mass[35] > 0.0);
    QVERIFY(std::isfinite(largeOffset.mass[35]));
    QVERIFY(close(largeOffset.totalMass,1.0e-90,3.0e-12,0.0));
    // Independently scaled homogeneous Hermite and rotary-inertia integrals.
    QVERIFY(close(largeOffset.mass[35],1.0e220*(1.0/105.0+
                  2.0/15.0*(1.5*1.5+1.0/12.0)),3.0e-12,0.0));
    offset.material = RStructuralPiezoContinuumMaterial::isotropic(1.0e-100,0.0,1.0e-100);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoBeamLayerKernel::assemble(1.0,{offset}),std::invalid_argument);
}

void TestStructuralPiezoBeamLayerKernel::representableOperatorsAfterVolumeUnderflow()
{
    auto item = layer(-5.0e-51,5.0e-51);
    item.width = 1.0e-300;
    item.material = RStructuralPiezoContinuumMaterial::isotropic(1.0e200,0.0,1.0e200);
    // The geometric volume is 1e-400, but material scaling makes every
    // nonzero operator entry representable. A precomputed volume loses it.
    const auto result = RStructuralPiezoBeamLayerKernel::assemble(1.0e-50,{item});
    QVERIFY(close(result.totalMass,1.0e-200,3.0e-12,0.0));
    QVERIFY(close(result.stiffness[0],1.0e-100,3.0e-12,0.0));
    // Independent consistent Hermite mass plus rotary inertia:
    // rho*A*L^3/105 + rho*I*L*(2/15).
    QVERIFY(close(result.mass[35],(1.0/105.0+1.0/90.0)*1.0e-300,3.0e-12,0.0));
    for (size_t i=0;i<6;++i)
    {
        QVERIFY(result.mass[i*6+i] > 0.0);
        QVERIFY(result.stiffness[i*6+i] > 0.0);
    }
}

void TestStructuralPiezoBeamLayerKernel::rejectsInvalidGeometryMaterialAndOverlap()
{
    const auto good = layer(-0.005,0.005);
    for (double length : {0.0,-1.0,std::numeric_limits<double>::infinity(),
                          std::numeric_limits<double>::quiet_NaN()})
        QVERIFY_EXCEPTION_THROWN(RStructuralPiezoBeamLayerKernel::assemble(length,{good}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoBeamLayerKernel::assemble(1.0,{}),std::invalid_argument);
    auto invalid = good;
    invalid.width = -1.0;
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoBeamLayerKernel::assemble(1.0,{invalid}),std::invalid_argument);
    invalid = good;
    invalid.upperZ = invalid.lowerZ;
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoBeamLayerKernel::assemble(1.0,{invalid}),std::invalid_argument);
    invalid = good;
    std::swap(invalid.lowerZ,invalid.upperZ);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoBeamLayerKernel::assemble(1.0,{invalid}),std::invalid_argument);
    invalid = good;
    invalid.coverage[0] = 0.5;
    invalid.coverage[1] = -0.5;
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoBeamLayerKernel::assemble(1.0,{invalid}),std::invalid_argument);
    invalid = good;
    invalid.coverage[0] = -1.1;
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoBeamLayerKernel::assemble(1.0,{invalid}),std::invalid_argument);
    invalid = good;
    invalid.coverage[2] = -0.5;
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoBeamLayerKernel::assemble(1.0,{invalid}),std::invalid_argument);
    invalid = good;
    invalid.material.density = -1.0;
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoBeamLayerKernel::assemble(1.0,{invalid}),std::invalid_argument);
    invalid = good;
    invalid.material.elasticity[0] = -1.0;
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoBeamLayerKernel::assemble(1.0,{invalid}),std::invalid_argument);
    invalid = layer(0.005,0.007,true);
    invalid.material.permittivity[8] = 0.0;
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoBeamLayerKernel::assemble(1.0,{invalid}),std::invalid_argument);
    invalid = good;
    invalid.lowerZ = std::numeric_limits<double>::quiet_NaN();
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoBeamLayerKernel::assemble(1.0,{invalid}),std::invalid_argument);
    invalid = good;
    invalid.portId = QStringLiteral("illegal_passive_port");
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoBeamLayerKernel::assemble(1.0,{invalid}),std::invalid_argument);
    invalid = layer(0.005,0.007,true);
    invalid.portId.clear();
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoBeamLayerKernel::assemble(1.0,{invalid}),std::invalid_argument);
    invalid = layer(0.005,0.007,true);
    invalid.polarity = 0;
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoBeamLayerKernel::assemble(1.0,{invalid}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoBeamLayerKernel::assemble(1.0,{good,good}),std::invalid_argument);
    auto overlaps = good;
    overlaps.lowerZ = 0.0;
    overlaps.upperZ = 0.01;
    overlaps.width = 0.001;
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoBeamLayerKernel::assemble(1.0,{good,overlaps}),std::invalid_argument);
    // An adjoining z boundary, or disjoint x intervals, never overlaps volume.
    const auto adjoining = RStructuralPiezoBeamLayerKernel::assemble(1.0,{good,layer(0.005,0.01)});
    QVERIFY(adjoining.totalMass > 0.0);
    auto left = good,right = good;
    left.coverage[1] = 0.0;
    right.coverage[0] = 0.0;
    const auto partitioned = RStructuralPiezoBeamLayerKernel::assemble(1.0,{left,right});
    const auto whole = RStructuralPiezoBeamLayerKernel::assemble(1.0,{good});
    QVERIFY(close(partitioned.totalMass,whole.totalMass));
    for (size_t i=0;i<36;++i) QVERIFY(close(partitioned.mass[i],whole.mass[i]));
    invalid = good;
    invalid.width = std::numeric_limits<double>::max();
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoBeamLayerKernel::assemble(
        std::numeric_limits<double>::max(),{invalid}),std::invalid_argument);
}

QTEST_APPLESS_MAIN(TestStructuralPiezoBeamLayerKernel)

#include "tst_structural_piezo_beam_layer_kernel.moc"

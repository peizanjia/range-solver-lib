#include <QtTest>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "rstructural_piezo_element_operators.h"

namespace
{
using Nodes=std::array<std::array<double,2>,4>;
using Vector=std::vector<double>;
const Nodes square{{{{0.0,0.0}},{{2.0,0.0}},{{2.0,2.0}},{{0.0,2.0}}}};
const Nodes distorted{{{{-0.2,0.1}},{{2.1,-0.1}},{{2.5,1.8}},{{0.1,2.3}}}};

RStructuralPiezoLayer passive(double lo=-0.1,double hi=0.1)
{
    RStructuralPiezoLayer layer;
    layer.material=RStructuralPiezoContinuumMaterial::isotropic(30.0,0.2,7.0);
    layer.lowerZ=lo; layer.upperZ=hi;
    return layer;
}

RStructuralPiezoLayer active(double lo,double hi,const QString &port=QStringLiteral("p"))
{
    auto layer=passive(lo,hi);
    layer.material.piezoelectric=true;
    layer.material.coupling[12]=2.0;
    layer.material.coupling[13]=2.0;
    layer.material.permittivity={{4.0,0.0,0.0,0.0,4.0,0.0,0.0,0.0,4.0}};
    layer.portId=port;
    return layer;
}

double dot(const Vector &a,const Vector &b)
{
    double value=0.0;
    for (std::size_t i=0;i<a.size();++i) value+=a[i]*b[i];
    return value;
}

Vector product(const Vector &matrix,const Vector &q)
{
    Vector r(q.size(),0.0);
    for (std::size_t i=0;i<q.size();++i)
        for (std::size_t j=0;j<q.size();++j) r[i]+=matrix[i*q.size()+j]*q[j];
    return r;
}

double quadratic(const Vector &matrix,const Vector &q) { return dot(q,product(matrix,q)); }

bool close(double a,double b,double relative=1.0e-11,double absolute=1.0e-12)
{
    return std::abs(a-b)<=absolute+relative*std::max(std::abs(a),std::abs(b));
}

double area(const Nodes &nodes)
{
    double a=0.0;
    for (unsigned i=0;i<4;++i)
    {
        const auto &p=nodes[i],&q=nodes[(i+1)%4];
        a+=p[0]*q[1]-p[1]*q[0];
    }
    return 0.5*a;
}

// Nodal values of an exact physical quadratic Kirchhoff field. Its tangential
// derivative at each straight edge midpoint is exact, independently of mapping.
Vector field(const Nodes &nodes,double ex,double ey,double gxy,
             double kx,double ky,double kxy)
{
    Vector q(20,0.0);
    for (unsigned i=0;i<4;++i)
    {
        const double x=nodes[i][0],y=nodes[i][1];
        q[5*i]=ex*x+0.5*gxy*y;
        q[5*i+1]=ey*y+0.5*gxy*x;
        q[5*i+2]=0.5*kx*x*x+0.5*ky*y*y+kxy*x*y;
        q[5*i+3]=ky*y+kxy*x;
        q[5*i+4]=-kx*x-kxy*y;
    }
    return q;
}

// Independent cyclic Jacobi spectrum, solely for the small element rank test.
Vector eigenvalues(Vector a,unsigned n)
{
    for (unsigned sweep=0;sweep<100;++sweep)
    {
        double off=0.0,scale=0.0;
        for (unsigned i=0;i<n;++i) scale=std::max(scale,std::abs(a[i*n+i]));
        for (unsigned p=0;p<n;++p)
            for (unsigned q=p+1;q<n;++q)
            {
                const double apq=a[p*n+q];
                off=std::max(off,std::abs(apq));
                if (std::abs(apq)<=1.0e-15*scale) continue;
                const double angle=0.5*std::atan2(2.0*apq,a[q*n+q]-a[p*n+p]);
                const double c=std::cos(angle),s=std::sin(angle);
                const double app=a[p*n+p],aqq=a[q*n+q];
                a[p*n+p]=c*c*app-2.0*s*c*apq+s*s*aqq;
                a[q*n+q]=s*s*app+2.0*s*c*apq+c*c*aqq;
                a[p*n+q]=a[q*n+p]=0.0;
                for (unsigned k=0;k<n;++k)
                    if (k!=p && k!=q)
                    {
                        const double akp=a[k*n+p],akq=a[k*n+q];
                        a[k*n+p]=a[p*n+k]=c*akp-s*akq;
                        a[k*n+q]=a[q*n+k]=s*akp+c*akq;
                    }
            }
        if (off<=1.0e-13*scale) break;
    }
    Vector values(n);
    for (unsigned i=0;i<n;++i) values[i]=a[i*n+i];
    std::sort(values.begin(),values.end());
    return values;
}

Vector solve(Vector a,Vector rhs)
{
    const std::size_t n=rhs.size();
    // Cholesky uses the positive supported plate, without a numerical shift.
    for (std::size_t i=0;i<n;++i)
        for (std::size_t j=0;j<=i;++j)
        {
            double v=a[i*n+j];
            for (std::size_t k=0;k<j;++k) v-=a[i*n+k]*a[j*n+k];
            if (i==j)
            {
                if (!(v>0.0)) throw std::runtime_error("Supported plate is not positive definite");
                a[i*n+i]=std::sqrt(v);
            }
            else a[i*n+j]=v/a[j*n+j];
        }
    for (std::size_t i=0;i<n;++i)
    {
        for (std::size_t j=0;j<i;++j) rhs[i]-=a[i*n+j]*rhs[j];
        rhs[i]/=a[i*n+i];
    }
    for (std::size_t i=n;i-->0;)
    {
        for (std::size_t j=i+1;j<n;++j) rhs[i]-=a[j*n+i]*rhs[j];
        rhs[i]/=a[i*n+i];
    }
    return rhs;
}

// Actual multi-element assembly of the source 20-DOF kernel, unit square,
// sin(pi*x)sin(pi*y) pressure, hard simply supported edges. Membrane is fixed
// here and decoupled for this symmetric isotropic section; no surrogate beam.
double supportedPlate(unsigned divisions,double thickness)
{
    const unsigned side=divisions+1,nodeCount=side*side;
    std::vector<int> equation(3*nodeCount,-1);
    unsigned equations=0;
    for (unsigned y=0;y<=divisions;++y)
        for (unsigned x=0;x<=divisions;++x)
        {
            const unsigned node=y*side+x;
            if (x>0 && x<divisions && y>0 && y<divisions) equation[3*node]=int(equations++);
            if (x>0 && x<divisions) equation[3*node+1]=int(equations++);
            if (y>0 && y<divisions) equation[3*node+2]=int(equations++);
        }
    Vector stiffness(equations*equations,0.0),force(equations,0.0);
    auto layer=passive(-0.5*thickness,0.5*thickness);
    layer.material=RStructuralPiezoContinuumMaterial::isotropic(1.0,0.3,1.0);
    const double h=1.0/divisions,pi=std::acos(-1.0);
    const double gp[3]={-std::sqrt(3.0/5.0),0.0,std::sqrt(3.0/5.0)};
    const double gw[3]={5.0/9.0,8.0/9.0,5.0/9.0};
    for (unsigned y=0;y<divisions;++y)
        for (unsigned x=0;x<divisions;++x)
        {
            const std::array<unsigned,4> ids{{y*side+x,y*side+x+1,(y+1)*side+x+1,(y+1)*side+x}};
            const Nodes nodes{{{{x*h,y*h}},{{(x+1)*h,y*h}},{{(x+1)*h,(y+1)*h}},{{x*h,(y+1)*h}}}};
            const auto op=RStructuralPiezoPlateMitc4Kernel::assemble(nodes,{layer});
            for (unsigned i=0;i<12;++i)
            {
                const int row=equation[3*ids[i/3]+i%3];
                if (row<0) continue;
                for (unsigned j=0;j<12;++j)
                {
                    const int col=equation[3*ids[j/3]+j%3];
                    if (col>=0) stiffness[std::size_t(row)*equations+std::size_t(col)]
                        +=op.stiffness[(5*(i/3)+2+i%3)*20+5*(j/3)+2+j%3];
                }
            }
            for (unsigned ix=0;ix<3;++ix)
                for (unsigned iy=0;iy<3;++iy)
                {
                    const double xi=gp[ix],eta=gp[iy];
                    const double pressure=std::sin(pi*(x+0.5*(1.0+xi))*h)
                                          *std::sin(pi*(y+0.5*(1.0+eta))*h);
                    const double n[4]={0.25*(1-xi)*(1-eta),0.25*(1+xi)*(1-eta),
                                       0.25*(1+xi)*(1+eta),0.25*(1-xi)*(1+eta)};
                    for (unsigned i=0;i<4;++i)
                    {
                        const int row=equation[3*ids[i]];
                        if (row>=0) force[std::size_t(row)]+=pressure*n[i]*h*h*0.25*gw[ix]*gw[iy];
                    }
                }
        }
    const auto q=solve(stiffness,force);
    // The independent Navier single Fourier-mode Reissner-Mindlin solution.
    const double rigidity=thickness*thickness*thickness/(12.0*(1.0-0.3*0.3));
    const double shear=(5.0/6.0)*thickness/(2.0*(1.0+0.3));
    const double exact=1.0/(4.0*std::pow(pi,4)*rigidity)+1.0/(2.0*pi*pi*shear);
    const unsigned center=(divisions/2)*side+divisions/2;
    return q[std::size_t(equation[3*center])]/exact;
}
}

class TestStructuralPiezoPlateMitc4Kernel:public QObject
{
    Q_OBJECT
private slots:
    void exactConsistentMassCentroidAndRotaryInertia();
    void membraneAndConstantCurvaturePatchOnDistortedQ4();
    void sixRigidModesAndNoAdditionalZeroEnergyModes();
    void symmetricAndEccentricMembraneBendingCoupling();
    void anisotropicNormalShearCouplingAndCorrection();
    void isolatedPortsPolarityReciprocityAndOpenShortDifference();
    void partialCoverageUsesHostKinematicsAndActualArea();
    void inPlaneRotationAndTranslationInvariance();
    void simplySupportedMeshAndThicknessConvergence();
    void rejectInvalidGeometryLayersAndShearCorrection();
};

void TestStructuralPiezoPlateMitc4Kernel::exactConsistentMassCentroidAndRotaryInertia()
{
    const auto layer=passive(0.2,0.5);
    const auto op=RStructuralPiezoPlateMitc4Kernel::assemble(square,{layer});
    QCOMPARE(op.dofs,20u);
    QVERIFY(op.portIds.empty()); QVERIFY(op.capacitance.empty());
    const double t=0.3,rho=7.0,z=0.35,i0=rho*t,i1=i0*z,i2=i0*(z*z+t*t/12.0);
    QVERIFY(close(op.totalMass,4.0*i0));
    const double pattern[4][4]={{4,2,1,2},{2,4,2,1},{1,2,4,2},{2,1,2,4}};
    for (unsigned i=0;i<4;++i)
        for (unsigned j=0;j<4;++j)
        {
            const double nn=4.0*pattern[i][j]/36.0;
            QVERIFY(close(op.mass[(5*i)*20+5*j],nn*i0));
            QVERIFY(close(op.mass[(5*i+4)*20+5*j+4],nn*i2));
            QVERIFY(close(op.mass[(5*i)*20+5*j+4],nn*i1));
            QVERIFY(close(op.mass[(5*i+1)*20+5*j+3],-nn*i1));
        }
    Vector q(20,0.0);
    for (unsigned i=0;i<4;++i) {q[5*i]=0.7; q[5*i+4]=-0.4; q[5*i+1]=0.2; q[5*i+3]=0.6;}
    const double expected=4.0*(i0*(0.7*0.7+0.2*0.2)
        +2.0*i1*(0.7*(-0.4)-0.2*0.6)+i2*(0.4*0.4+0.6*0.6));
    QVERIFY(close(quadratic(op.mass,q),expected));
    const auto spectrum=eigenvalues(op.mass,20);
    QVERIFY(spectrum.front()>0.0);
}

void TestStructuralPiezoPlateMitc4Kernel::membraneAndConstantCurvaturePatchOnDistortedQ4()
{
    for (const Nodes &nodes:{square,distorted})
        for (const auto &layer:{passive(),passive(0.1,0.35)})
        {
            const auto op=RStructuralPiezoPlateMitc4Kernel::assemble(nodes,{layer});
            const double t=layer.upperZ-layer.lowerZ,m=0.5*(layer.upperZ+layer.lowerZ);
            const double e=30.0,nu=0.2,c=e/(1.0-nu*nu),g=e/(2.0*(1.0+nu));
            const double ex=0.07,ey=-0.03,xy=0.04,kx=0.13,ky=-0.11,kxy=0.08;
            const auto q=field(nodes,ex,ey,xy,kx,ky,kxy);
            const double membrane=c*(ex*ex+ey*ey+2*nu*ex*ey)+g*xy*xy;
            const double bending=c*(kx*kx+ky*ky+2*nu*kx*ky)+g*4*kxy*kxy;
            const double cross=-c*(ex*kx+ey*ky+nu*(ex*ky+ey*kx))-g*xy*2*kxy;
            const double expected=area(nodes)*t*(membrane+2*m*cross+(m*m+t*t/12)*bending);
            QVERIFY(close(quadratic(op.stiffness,q),expected));
            // Kirchhoff constant curvature has identically zero tied shear,
            // therefore changing kappa cannot change its energy.
            const auto other=RStructuralPiezoPlateMitc4Kernel::assemble(nodes,{layer},0.23);
            QVERIFY(close(quadratic(other.stiffness,q),expected));
        }
}

void TestStructuralPiezoPlateMitc4Kernel::sixRigidModesAndNoAdditionalZeroEnergyModes()
{
    for (const Nodes &nodes:{square,distorted})
    {
        const auto op=RStructuralPiezoPlateMitc4Kernel::assemble(nodes,{passive()});
        for (unsigned mode=0;mode<6;++mode)
        {
            Vector q(20,0.0);
            for (unsigned i=0;i<4;++i)
            {
                if (mode<3) q[5*i+mode]=1.0;
                if (mode==3) {q[5*i+2]=nodes[i][1]; q[5*i+3]=1.0;}
                if (mode==4) {q[5*i+2]=-nodes[i][0]; q[5*i+4]=1.0;}
                if (mode==5) {q[5*i]=-nodes[i][1]; q[5*i+1]=nodes[i][0];}
            }
            for (double residual:product(op.stiffness,q)) QVERIFY(std::abs(residual)<2.0e-13);
        }
        const auto spectrum=eigenvalues(op.stiffness,20);
        for (unsigned i=0;i<6;++i) QVERIFY(std::abs(spectrum[i])<1.0e-10*spectrum.back());
        QVERIFY(spectrum[6]>1.0e-5*spectrum.back());
    }
}

void TestStructuralPiezoPlateMitc4Kernel::symmetricAndEccentricMembraneBendingCoupling()
{
    const auto membrane=field(square,1,0,0,0,0,0),bend=field(square,0,0,0,1,0,0);
    const auto symmetric=RStructuralPiezoPlateMitc4Kernel::assemble(square,{passive(-.3,-.1),passive(.1,.3)});
    QVERIFY(std::abs(dot(membrane,product(symmetric.stiffness,bend)))<1.0e-13);
    const auto eccentric=RStructuralPiezoPlateMitc4Kernel::assemble(square,{passive(.1,.3)});
    const double expected=-area(square)*30.0/(1.0-.2*.2)*(.3*.3-.1*.1)/2.0;
    QVERIFY(close(dot(membrane,product(eccentric.stiffness,bend)),expected));
    const auto top=RStructuralPiezoPlateMitc4Kernel::assemble(square,{active(.1,.3)});
    const auto bottom=RStructuralPiezoPlateMitc4Kernel::assemble(square,{active(-.3,-.1)});
    QVERIFY(close(dot(top.coupling,bend),-dot(bottom.coupling,bend)));
    QVERIFY(close(dot(top.coupling,membrane),dot(bottom.coupling,membrane)));
    QVERIFY(close(top.capacitance[0],bottom.capacitance[0]));
}

void TestStructuralPiezoPlateMitc4Kernel::anisotropicNormalShearCouplingAndCorrection()
{
    auto layer=active(-.12,.18);
    layer.material.elasticity.fill(0.0);
    const double diagonal[6]={30,40,50,7,9,11};
    for (unsigned i=0;i<6;++i) layer.material.elasticity[6*i+i]=diagonal[i];
    for (const auto &entry:std::array<std::array<unsigned,2>,4>{{{{0,4}},{{1,3}},{{0,2}},{{2,4}}}})
    {
        const double v=(entry[0]==0 && entry[1]==4)?2.0:(entry[0]==1?3.0:(entry[0]==0?4.0:1.0));
        layer.material.elasticity[6*entry[0]+entry[1]]=v;
        layer.material.elasticity[6*entry[1]+entry[0]]=v;
    }
    const double ez[6]={2,-1,.6,.4,-.3,.2};
    for (unsigned i=0;i<6;++i) layer.material.coupling[12+i]=ez[i];
    const unsigned retained[5]={0,1,5,4,3};
    Vector q=field(distorted,.2,-.3,.15,0,0,0);
    for (unsigned i=0;i<4;++i) q[5*i+2]=.12*distorted[i][0]-.09*distorted[i][1];
    for (double correction:{1.0,0.37})
    {
        const auto op=RStructuralPiezoPlateMitc4Kernel::assemble(distorted,{layer},correction);
        const double strain[5]={.2,-.3,.15,.12*std::sqrt(correction),-.09*std::sqrt(correction)};
        double stiffness=0,coupling=0;
        for (unsigned i=0;i<5;++i)
        {
            const unsigned a=retained[i];
            coupling+=strain[i]*(ez[a]-ez[2]*layer.material.elasticity[12+a]/50.0);
            for (unsigned j=0;j<5;++j)
            {
                const unsigned b=retained[j];
                const double reduced=layer.material.elasticity[6*a+b]
                    -layer.material.elasticity[6*a+2]*layer.material.elasticity[12+b]/50.0;
                stiffness+=strain[i]*reduced*strain[j];
            }
        }
        QVERIFY(close(quadratic(op.stiffness,q),area(distorted)*.3*stiffness));
        QVERIFY(close(dot(op.coupling,q),area(distorted)*coupling));
        QVERIFY(close(op.capacitance[0],area(distorted)*(4.0+.6*.6/50.0)/.3));
    }
}

void TestStructuralPiezoPlateMitc4Kernel::isolatedPortsPolarityReciprocityAndOpenShortDifference()
{
    auto upper=active(.1,.2,QStringLiteral("b")),lower=active(-.2,-.1,QStringLiteral("a"));
    const auto op=RStructuralPiezoPlateMitc4Kernel::assemble(square,{upper,lower});
    QCOMPARE(op.portIds.size(),std::size_t(2));
    QCOMPARE(op.portIds[0],QStringLiteral("b")); QCOMPARE(op.portIds[1],QStringLiteral("a"));
    QCOMPARE(op.capacitance[1],0.0); QCOMPARE(op.capacitance[2],0.0);
    upper.polarity=-1;
    const auto reversed=RStructuralPiezoPlateMitc4Kernel::assemble(square,{upper,lower});
    for (unsigned i=0;i<20;++i)
    {
        QVERIFY(close(reversed.coupling[2*i],-op.coupling[2*i]));
        QVERIFY(close(reversed.coupling[2*i+1],op.coupling[2*i+1]));
    }
    QCOMPARE(reversed.capacitance,op.capacitance);
    lower.portId=QStringLiteral("b"); upper.polarity=1;
    const auto continuous=RStructuralPiezoPlateMitc4Kernel::assemble(square,{upper,lower});
    QCOMPARE(continuous.portIds.size(),std::size_t(1));
    for (unsigned i=0;i<20;++i) QVERIFY(close(continuous.coupling[i],op.coupling[2*i]+op.coupling[2*i+1]));
    QVERIFY(close(continuous.capacitance[0],op.capacitance[0]+op.capacitance[3]));
    Vector q=field(square,.07,-.03,.04,.13,-.11,.08),voltage{.4,-.2};
    const auto force=[&](const Vector &state,const Vector &v) {
        auto f=product(op.stiffness,state);
        for (unsigned i=0;i<20;++i) for (unsigned p=0;p<2;++p) f[i]+=op.coupling[2*i+p]*v[p];
        return f;
    };
    const auto charge=[&](const Vector &state,const Vector &v) {
        Vector charges(2,0.0);
        for (unsigned p=0;p<2;++p)
        {
            charges[p]=op.capacitance[3*p]*v[p];
            for (unsigned i=0;i<20;++i) charges[p]-=op.coupling[2*i+p]*state[i];
        }
        return charges;
    };
    const double delta=1.0e-5;
    for (unsigned p=0;p<2;++p)
    {
        Vector plus=voltage,minus=voltage; plus[p]+=delta; minus[p]-=delta;
        const auto fp=force(q,plus),fm=force(q,minus);
        for (unsigned i=0;i<20;++i)
        {
            QVERIFY(close((fp[i]-fm[i])/(2*delta),op.coupling[2*i+p],1.0e-8));
            auto qp=q,qm=q; qp[i]+=delta; qm[i]-=delta;
            QVERIFY(close((charge(qp,voltage)[p]-charge(qm,voltage)[p])/(2*delta),-op.coupling[2*i+p],1.0e-8,1.0e-9));
        }
    }
    Vector difference(400,0.0);
    for (unsigned i=0;i<20;++i) for (unsigned j=0;j<20;++j) for (unsigned p=0;p<2;++p)
        difference[20*i+j]+=op.coupling[2*i+p]*op.coupling[2*j+p]/op.capacitance[3*p];
    const auto spectrum=eigenvalues(difference,20);
    QVERIFY(spectrum.front()>-1.0e-12*spectrum.back());
    QVERIFY(spectrum.back()>0.0);
    Vector openVoltage(2,0.0);
    for (unsigned p=0;p<2;++p) for (unsigned i=0;i<20;++i)
        openVoltage[p]+=op.coupling[2*i+p]*q[i]/op.capacitance[3*p];
    for (double value:charge(q,openVoltage)) QVERIFY(std::abs(value)<1.0e-13);
    const auto openForce=force(q,openVoltage),shortForce=force(q,{0.0,0.0}),increment=product(difference,q);
    for (unsigned i=0;i<20;++i) QVERIFY(close(openForce[i]-shortForce[i],increment[i]));
}

void TestStructuralPiezoPlateMitc4Kernel::partialCoverageUsesHostKinematicsAndActualArea()
{
    auto layer=active(.1,.3);
    layer.coverage={{-1,0,-1,0}};
    const auto patch=RStructuralPiezoPlateMitc4Kernel::assemble(square,{layer});
    QVERIFY(close(patch.totalMass,7.0*.2));
    QVERIFY(close(patch.capacitance[0],4.0/.2));
    const auto curvature=field(square,0,0,0,.2,-.1,.05);
    layer.coverage={{-1,1,-1,1}};
    const auto full=RStructuralPiezoPlateMitc4Kernel::assemble(square,{layer});
    QVERIFY(close(quadratic(patch.stiffness,curvature),.25*quadratic(full.stiffness,curvature)));
    QVERIFY(close(dot(patch.coupling,curvature),.25*dot(full.coupling,curvature)));
    // A partial layer retains all host DOFs, rather than assembling a smaller
    // fictitious plate element or clipping the tying kinematics.
    QVERIFY(patch.mass[19*20+19]>0.0);
    std::vector<RStructuralPiezoLayer> tiles;
    for (unsigned x=0;x<2;++x) for (unsigned y=0;y<2;++y)
    {
        auto tile=layer;
        tile.coverage={{double(x)-1.0,double(x),double(y)-1.0,double(y)}};
        tiles.push_back(tile);
    }
    const auto tiled=RStructuralPiezoPlateMitc4Kernel::assemble(distorted,tiles);
    const auto distortedFull=RStructuralPiezoPlateMitc4Kernel::assemble(distorted,{layer});
    QVERIFY(close(tiled.totalMass,distortedFull.totalMass));
    for (unsigned i=0;i<400;++i) QVERIFY(close(tiled.mass[i],distortedFull.mass[i]));
    QVERIFY(close(tiled.capacitance[0],distortedFull.capacitance[0]));
    // Stiffness of a distorted map is rational and is quadrature qualified;
    // the constant physical patch energy is nevertheless exactly integrated.
    const auto q=field(distorted,.03,.04,.05,.07,.02,.06);
    QVERIFY(close(quadratic(tiled.stiffness,q),quadratic(distortedFull.stiffness,q)));
}

void TestStructuralPiezoPlateMitc4Kernel::inPlaneRotationAndTranslationInvariance()
{
    const double angle=.71,c=std::cos(angle),s=std::sin(angle);
    Nodes rotated=distorted;
    for (unsigned i=0;i<4;++i)
    {
        rotated[i][0]=c*distorted[i][0]-s*distorted[i][1]+17.0;
        rotated[i][1]=s*distorted[i][0]+c*distorted[i][1]-11.0;
    }
    const auto layer=active(.1,.3);
    const auto original=RStructuralPiezoPlateMitc4Kernel::assemble(distorted,{layer});
    const auto changed=RStructuralPiezoPlateMitc4Kernel::assemble(rotated,{layer});
    Vector transform(400,0.0);
    for (unsigned i=0;i<4;++i)
    {
        transform[(5*i)*20+5*i]=c; transform[(5*i)*20+5*i+1]=-s;
        transform[(5*i+1)*20+5*i]=s; transform[(5*i+1)*20+5*i+1]=c;
        transform[(5*i+2)*20+5*i+2]=1;
        transform[(5*i+3)*20+5*i+3]=c; transform[(5*i+3)*20+5*i+4]=-s;
        transform[(5*i+4)*20+5*i+3]=s; transform[(5*i+4)*20+5*i+4]=c;
    }
    for (unsigned i=0;i<20;++i)
    {
        double g=0.0;
        for (unsigned k=0;k<20;++k) g+=transform[k*20+i]*changed.coupling[k];
        QVERIFY(close(g,original.coupling[i]));
        for (unsigned j=0;j<20;++j)
        {
            double mass=0,stiffness=0;
            for (unsigned k=0;k<20;++k) for (unsigned l=0;l<20;++l)
            {
                const double factor=transform[k*20+i]*transform[l*20+j];
                mass+=factor*changed.mass[k*20+l]; stiffness+=factor*changed.stiffness[k*20+l];
            }
            QVERIFY(close(mass,original.mass[i*20+j]));
            QVERIFY(close(stiffness,original.stiffness[i*20+j]));
        }
    }
    QVERIFY(close(changed.capacitance[0],original.capacitance[0]));
}

void TestStructuralPiezoPlateMitc4Kernel::simplySupportedMeshAndThicknessConvergence()
{
    // Compare the same physical sine pressure with an independently derived
    // continuum Fourier solution; the FE solve uses the assembled native kernel.
    for (double thickness:{.1,.001,.0001})
    {
        const double coarse=supportedPlate(2,thickness);
        const double middle=supportedPlate(4,thickness);
        const double fine=supportedPlate(8,thickness);
        qInfo("MITC4 simply-supported h=%g: wFE/wExact n=2 %.9g n=4 %.9g n=8 %.9g",thickness,coarse,middle,fine);
        QVERIFY(std::abs(fine-1.0)<std::abs(middle-1.0));
        QVERIFY(std::abs(middle-1.0)<std::abs(coarse-1.0));
        QVERIFY(std::abs(fine-1.0)<.06);
    }
}

void TestStructuralPiezoPlateMitc4Kernel::rejectInvalidGeometryLayersAndShearCorrection()
{
    auto nodes=square;
    std::swap(nodes[1],nodes[3]);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoPlateMitc4Kernel::assemble(nodes,{passive()}),std::invalid_argument);
    nodes=square; std::swap(nodes[1],nodes[2]);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoPlateMitc4Kernel::assemble(nodes,{passive()}),std::invalid_argument);
    nodes=square; nodes[2]={{.5,.5}};
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoPlateMitc4Kernel::assemble(nodes,{passive()}),std::invalid_argument);
    nodes=square; nodes[1]=nodes[0];
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoPlateMitc4Kernel::assemble(nodes,{passive()}),std::invalid_argument);
    nodes=square; nodes[0][0]=std::numeric_limits<double>::quiet_NaN();
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoPlateMitc4Kernel::assemble(nodes,{passive()}),std::invalid_argument);
    for (double kappa:{0.0,-.1,1.01,std::numeric_limits<double>::infinity()})
        QVERIFY_EXCEPTION_THROWN(RStructuralPiezoPlateMitc4Kernel::assemble(square,{passive()},kappa),std::invalid_argument);
    auto layer=passive(); layer.width=2;
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoPlateMitc4Kernel::assemble(square,{layer}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoPlateMitc4Kernel::assemble(square,{}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoPlateMitc4Kernel::assemble(square,{passive(),passive()}),std::invalid_argument);
    layer=passive(); layer.coverage={{-.2,.3,-.4,.5}};
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoPlateMitc4Kernel::assemble(square,{passive(),layer}),std::invalid_argument);
    layer=passive(); layer.upperZ=layer.lowerZ;
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoPlateMitc4Kernel::assemble(square,{layer}),std::invalid_argument);
    layer=passive(); layer.coverage={{-1,1.1,-1,1}};
    QVERIFY_EXCEPTION_THROWN(RStructuralPiezoPlateMitc4Kernel::assemble(square,{layer}),std::invalid_argument);
    // A geometrically valid narrow plate is accepted without an arbitrary
    // aspect-ratio cutoff, mass shift, or Jacobian clamp.
    nodes=square; nodes[2][1]=nodes[3][1]=1.0e-6;
    const auto narrow=RStructuralPiezoPlateMitc4Kernel::assemble(nodes,{passive()});
    QVERIFY(close(narrow.totalMass,2.0e-6*.2*7.0));
}

QTEST_APPLESS_MAIN(TestStructuralPiezoPlateMitc4Kernel)

#include "tst_structural_piezo_plate_mitc4_kernel.moc"

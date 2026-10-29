#include <QtTest>

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <vector>

#include <rbl_error.h>
#include <rml_model.h>
#include <rml_boundary_condition.h>
#include <rml_environment_condition.h>
#include "rstructural_piezo_layered_model_adapter.h"

namespace
{
constexpr double side=0.1;
constexpr double area=side*side;
constexpr double e31=5.0;
constexpr double epsilon33=1.5e-8;
constexpr double shearCorrection=5.0/6.0;
constexpr size_t dofs=20;
using Vector=std::vector<double>;
using Row=std::array<double,dofs>;
using StrainRows=std::array<Row,3>;
using Dof=RStructuralPiezoPhysicalDof;

struct LayerData
{
    double lower,upper,young,density,permittivity;
};

constexpr std::array<LayerData,2> layers{{
    {-0.001,0.001,70.0e9,2700.0,0.0},
    { 0.001,0.0015,60.0e9,7500.0,epsilon33}
}};

struct Section
{
    double a=0,b=0,d=0,as=0,bs=0,ds=0;
    double i0=0,i1=0,i2=0;
};

Section section()
{
    Section result;
    for(const auto &layer:layers)
    {
        const double t=layer.upper-layer.lower;
        const double z1=layer.lower,z2=layer.upper;
        const double first=(z2*z2-z1*z1)/2.0;
        const double second=(z2*z2*z2-z1*z1*z1)/3.0;
        result.a+=layer.young*t;
        result.b+=layer.young*first;
        result.d+=layer.young*second;
        result.as+=0.5*layer.young*t;
        result.bs+=0.5*layer.young*first;
        result.ds+=0.5*layer.young*second;
        result.i0+=layer.density*t;
        result.i1+=layer.density*first;
        result.i2+=layer.density*second;
    }
    return result;
}

bool near(double actual,double expected,double relative,double absolute)
{
    return std::isfinite(actual) && std::isfinite(expected) &&
        std::abs(actual-expected)<=absolute+relative*std::max(std::abs(actual),std::abs(expected));
}

QString difference(const char *name,size_t row,size_t column,double actual,double expected)
{
    return QString("%1(%2,%3): actual=%4 expected=%5 error=%6")
        .arg(QString::fromLatin1(name)).arg(qulonglong(row)).arg(qulonglong(column))
        .arg(actual,0,'g',17).arg(expected,0,'g',17).arg(actual-expected,0,'g',17);
}

struct PlateOracle
{
    Vector mass=Vector(dofs*dofs,0.0),stiffness=Vector(dofs*dofs,0.0);
    Vector coupling=Vector(dofs,0.0);
    double capacitance=0.0;
    double activeCenter=0.00125;
    Section moments=section();
};

struct Shape
{
    std::array<double,4> n{},dx{},dy{},dxi{},deta{};
};

Shape shape(double xi,double eta)
{
    constexpr std::array<double,4> sx{{-1,1,1,-1}};
    constexpr std::array<double,4> sy{{-1,-1,1,1}};
    Shape result;
    for(size_t i=0;i<4;++i)
    {
        result.n[i]=0.25*(1+sx[i]*xi)*(1+sy[i]*eta);
        result.dxi[i]=0.25*sx[i]*(1+sy[i]*eta);
        result.deta[i]=0.25*sy[i]*(1+sx[i]*xi);
        result.dx[i]=2.0*result.dxi[i]/side;
        result.dy[i]=2.0*result.deta[i]/side;
    }
    return result;
}

Row covariantXiAtEdge(double eta)
{
    const Shape s=shape(0.0,eta);
    Row row{};
    for(size_t i=0;i<4;++i)
    {
        row[5*i+2]=s.dxi[i];
        row[5*i+4]=0.5*side*s.n[i];
    }
    return row;
}

Row covariantEtaAtEdge(double xi)
{
    const Shape s=shape(xi,0.0);
    Row row{};
    for(size_t i=0;i<4;++i)
    {
        row[5*i+2]=s.deta[i];
        row[5*i+3]=-0.5*side*s.n[i];
    }
    return row;
}

PlateOracle analyticPlate()
{
    PlateOracle result;
    const double consistentShapeMass[4][4]={{4,2,1,2},{2,4,2,1},{1,2,4,2},{2,1,2,4}};
    for(size_t i=0;i<4;++i)
        for(size_t j=0;j<4;++j)
        {
            const double nij=area*consistentShapeMass[i][j]/36.0;
            for(size_t component=0;component<3;++component)
                result.mass[(5*i+component)*dofs+5*j+component]=nij*result.moments.i0;
            result.mass[(5*i+3)*dofs+5*j+3]=nij*result.moments.i2;
            result.mass[(5*i+4)*dofs+5*j+4]=nij*result.moments.i2;
            result.mass[(5*i+0)*dofs+5*j+4]=nij*result.moments.i1;
            result.mass[(5*i+4)*dofs+5*j+0]=nij*result.moments.i1;
            result.mass[(5*i+1)*dofs+5*j+3]=-nij*result.moments.i1;
            result.mass[(5*i+3)*dofs+5*j+1]=-nij*result.moments.i1;
        }

    const double g=1.0/std::sqrt(3.0);
    const Row xiMinus=covariantXiAtEdge(-1.0),xiPlus=covariantXiAtEdge(1.0);
    const Row etaMinus=covariantEtaAtEdge(-1.0),etaPlus=covariantEtaAtEdge(1.0);
    const double jacobian=area/4.0;
    const double shearScale=std::sqrt(shearCorrection);
    for(double xi:{-g,g})
        for(double eta:{-g,g})
        {
            const Shape s=shape(xi,eta);
            StrainRows membrane{},curvature{};
            Row covXi{},covEta{};
            for(size_t node=0;node<4;++node)
            {
                const size_t ux=5*node,uy=ux+1,rx=ux+3,ry=ux+4;
                membrane[0][ux]=s.dx[node];
                membrane[1][uy]=s.dy[node];
                membrane[2][ux]=s.dy[node]; membrane[2][uy]=s.dx[node];
                curvature[0][ry]=s.dx[node];
                curvature[1][rx]=-s.dy[node];
                curvature[2][rx]=-s.dx[node]; curvature[2][ry]=s.dy[node];
            }
            for(size_t i=0;i<dofs;++i)
            {
                covXi[i]=0.5*((1.0-eta)*xiMinus[i]+(1.0+eta)*xiPlus[i]);
                covEta[i]=0.5*((1.0-xi)*etaMinus[i]+(1.0+xi)*etaPlus[i]);
            }
            Row shearX{},shearY{};
            for(size_t i=0;i<dofs;++i)
            {
                shearX[i]=shearScale*(2.0/side)*covXi[i];
                shearY[i]=shearScale*(2.0/side)*covEta[i];
            }
            for(size_t i=0;i<dofs;++i)
                for(size_t j=0;j<dofs;++j)
                {
                    double value=0.0;
                    for(size_t strain=0;strain<3;++strain)
                    {
                        const double a=strain==2?result.moments.as:result.moments.a;
                        const double b=strain==2?result.moments.bs:result.moments.b;
                        const double d=strain==2?result.moments.ds:result.moments.d;
                        value+=a*membrane[strain][i]*membrane[strain][j];
                        value+=b*(membrane[strain][i]*curvature[strain][j]+
                                  curvature[strain][i]*membrane[strain][j]);
                        value+=d*curvature[strain][i]*curvature[strain][j];
                    }
                    value+=result.moments.as*(shearX[i]*shearX[j]+shearY[i]*shearY[j]);
                    result.stiffness[i*dofs+j]+=jacobian*value;
                }
        }

    const double gradIntegral[4]={-side/2.0,side/2.0,side/2.0,-side/2.0};
    for(size_t i=0;i<4;++i)
    {
        result.coupling[5*i]=e31*gradIntegral[i];
        result.coupling[5*i+4]=e31*result.activeCenter*gradIntegral[i];
    }
    result.capacitance=epsilon33*area/(layers[1].upper-layers[1].lower);
    return result;
}

Vector patchDisplacement(double epsilon,double curvature)
{
    Vector q(dofs,0.0);
    for(size_t node=0;node<4;++node)
    {
        const double x=(node==1 || node==2)?side:0.0;
        q[5*node]=epsilon*x;
        q[5*node+2]=0.5*curvature*x*x;
        q[5*node+4]=-curvature*x;
    }
    return q;
}

Vector multiply(const Vector &matrix,const Vector &x,size_t n)
{
    Vector result(n,0.0);
    for(size_t i=0;i<n;++i) for(size_t j=0;j<n;++j) result[i]+=matrix[i*n+j]*x[j];
    return result;
}

double dot(const Vector &a,const Vector &b)
{
    double result=0.0;
    for(size_t i=0;i<a.size();++i) result+=a[i]*b[i];
    return result;
}

RStructuralPiezoLayeredMaterial material(const char *id,double young,double density,bool active)
{
    RStructuralPiezoLayeredMaterial result;
    result.materialId=QString::fromLatin1(id);
    result.sourceId="independent-analytic-example";
    result.sourceSha256=QString(64,'a');
    for(size_t i=0;i<3;++i) result.elasticityOrCompliance[i*6+i]=young;
    for(size_t i=3;i<6;++i) result.elasticityOrCompliance[i*6+i]=0.5*young;
    result.density=density;
    if(active)
    {
        result.piezoelectric=true;
        result.coupling[12]=e31;
        result.permittivity[0]=epsilon33;
        result.permittivity[4]=epsilon33;
        result.permittivity[8]=epsilon33;
    }
    return result;
}

RModel plateScene(bool openCircuit)
{
    RModel model;
    const double xy[4][2]={{0,0},{side,0},{side,side},{0,side}};
    for(uint i=0;i<4;++i)
    {
        model.addNode(RNode(xy[i][0],xy[i][1],0.0));
        model.setNodePhysicalId(i,101+11*i);
    }
    RElement element(R_ELEMENT_QUAD1);
    for(uint i=0;i<4;++i) element.setNodeId(i,i);
    model.addElement(element,true,0);
    model.setElementPhysicalId(0,701);
    for(uint i=0;i<model.getNSurfaces();++i)
    {
        auto &surface=model.getSurface(i);
        surface.removeEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
        surface.removeEnvironmentCondition(R_ENVIRONMENT_CONDITION_TEMPERATURE);
    }
    model.setProblemTaskTree(RProblemTaskItem::generateDefaultTree(R_PROBLEM_STRESS));
    model.getProblemSetup().setMassMatrixType(RMassMatrixType::Consistent);
    auto &dynamics=model.getProblemSetup().getStructuralDynamicsSetup();
    dynamics.setEnabled(true);
    dynamics.setMethod(RStructuralTimeIntegrationMethod::Newmark);
    dynamics.setTimeStep(1.0e-5);
    dynamics.setNSteps(8);
    model.getTimeSolver().setInputStartTime(1.0e-5);

    auto &setup=model.getProblemSetup().getStructuralPiezoLayeredSetup();
    setup.enabled=true; setup.family=RStructuralPiezoLayeredModelFamily::FlatPlate;
    setup.materials={material("host",70.0e9,2700.0,false),material("active-pzt",60.0e9,7500.0,true)};
    RStructuralPiezoLayeredModelLayer host; host.materialId="host"; host.lowerZ=-0.001; host.upperZ=0.001; host.width=1.0;
    RStructuralPiezoLayeredModelLayer patch; patch.materialId="active-pzt"; patch.lowerZ=0.001; patch.upperZ=0.0015; patch.width=1.0; patch.portId="electrode";
    setup.members={{701,{host,patch}}};
    RStructuralPiezoLayeredPort port;
    port.readout.portId="electrode";
    port.readout.role=openCircuit?RStructuralPiezoPortRole::Sensor:RStructuralPiezoPortRole::ActuatorSensor;
    port.readout.readout=openCircuit?RStructuralPiezoReadoutMode::OpenCircuitVoltage:RStructuralPiezoReadoutMode::ActuatorCharge;
    port.initialValue=openCircuit?3.75e-6:25.0;
    setup.ports={port};
    for(uint node:{101u,134u})
        for(Dof dof:{Dof::UX,Dof::UY,Dof::UZ,Dof::RX,Dof::RY}) setup.supports.push_back({node,dof,0.0});

    const PlateOracle oracle=analyticPlate();
    const Vector q=patchDisplacement(1.0e-4,0.02);
    Vector internal=multiply(oracle.stiffness,q,dofs);
    for(size_t i=0;i<dofs;++i) internal[i]+=oracle.coupling[i]*25.0;
    setup.workpoint.available=true;
    for(uint node=0;node<4;++node)
    {
        RStructuralPiezoLayeredNodeWorkpoint field;
        field.physicalNodeId=101+11*node;
        for(size_t component=0;component<5;++component) field.q0[component]=q[5*node+component];
        // The two x=side nodes are free. Apply the independently integrated
        // section-edge work there and let the fixed left edge carry reaction.
        if(node==1 || node==2)
        {
            for(size_t component=0;component<5;++component) field.f0[component]=internal[5*node+component];
        }
        setup.workpoint.nodes.push_back(field);
        RStructuralPiezoLayeredNodeField applied;
        applied.physicalNodeId=field.physicalNodeId;
        applied.values=field.f0;
        setup.backgroundForces.push_back(applied);
    }
    setup.workpoint.ports={{"electrode",25.0,3.75e-6}};
    setup.bindToModel(model);
    return model;
}

double charge(const RStructuralPiezoElectricalResponse &response)
{
    if(response.ports.size()!=1) throw std::runtime_error("Plate oracle expected exactly one electrode");
    return response.ports[0].charge;
}
double voltage(const RStructuralPiezoElectricalResponse &response)
{
    if(response.ports.size()!=1) throw std::runtime_error("Plate oracle expected exactly one electrode");
    return response.ports[0].voltage;
}
}

class TestStructuralPiezoLayeredPlateOracle : public QObject
{
    Q_OBJECT
private slots:
    void independentSectionMatricesAndPatchWorkpoint();
    void openCircuitDrivenReciprocityAndTangent();
};

void TestStructuralPiezoLayeredPlateOracle::independentSectionMatricesAndPatchWorkpoint() try
{
    const PlateOracle expected=analyticPlate();
    QVERIFY(near(expected.moments.a,170.0e6,2e-14,2e-8));
    QVERIFY(near(expected.moments.b,37500.0,2e-14,1e-11));
    QVERIFY(near(expected.moments.d,94.16666666666667,2e-14,2e-12));
    QVERIFY(near(expected.moments.i0,9.15,2e-14,2e-14));
    QVERIFY(near(expected.moments.i1,0.0046875,2e-14,2e-17));
    QVERIFY(near(expected.moments.i2,7.7375e-6,2e-14,2e-18));
    QVERIFY(near(expected.capacitance,3.0e-7,2e-14,2e-20));

    RModel model=plateScene(false);
    const auto &setup=model.getProblemSetup().getStructuralPiezoLayeredSetup();
    const RStructuralPiezoLayeredModelAdapter adapter(model,setup);
    const auto &actual=adapter.mesh().operators;
    QCOMPARE(actual.dofs,uint(dofs));
    QCOMPARE(actual.portIds,(std::vector<QString>{"electrode"}));
    for(size_t i=0;i<dofs*dofs;++i)
    {
        const size_t row=i/dofs,column=i%dofs;
        QVERIFY2(near(actual.mass[i],expected.mass[i],4e-12,2e-15),qPrintable(difference("M",row,column,actual.mass[i],expected.mass[i])));
        QVERIFY2(near(actual.stiffness[i],expected.stiffness[i],5e-11,4e-7),qPrintable(difference("K",row,column,actual.stiffness[i],expected.stiffness[i])));
    }
    for(size_t i=0;i<dofs;++i)
        QVERIFY2(near(actual.coupling[i],expected.coupling[i],3e-12,2e-14),
            qPrintable(difference("G",i,0,actual.coupling[i],expected.coupling[i])));
    QVERIFY2(near(actual.capacitance[0],expected.capacitance,3e-12,2e-20),
        qPrintable(difference("Ce",0,0,actual.capacitance[0],expected.capacitance)));

    const Vector q=patchDisplacement(1.0e-4,0.02);
    const double qkq=dot(q,multiply(expected.stiffness,q,dofs));
    const double analyticEnergy=0.5*area*(expected.moments.a*1.0e-8-
        2.0*expected.moments.b*1.0e-4*0.02+expected.moments.d*0.0004);
    QVERIFY(near(analyticEnergy,0.007938333333333334,3e-13,1e-15));
    QVERIFY2(near(0.5*qkq,analyticEnergy,4e-12,2e-13),
        qPrintable(QString("plate patch energy actual=%1 expected=%2").arg(0.5*qkq,0,'g',17).arg(analyticEnergy,0,'g',17)));
    QVERIFY(near(qkq,0.015876666666666668,4e-12,4e-13));
    QVERIFY(near(dot(expected.coupling,q),3.75e-6,3e-13,2e-18));

    const auto &response=adapter.electrical().referenceResponse();
    QVERIFY(near(voltage(response),25.0,2e-13,2e-12));
    QVERIFY(near(charge(response),3.75e-6,3e-12,2e-18));
    Vector expectedForce=multiply(expected.stiffness,q,dofs);
    for(size_t i=0;i<dofs;++i) expectedForce[i]+=expected.coupling[i]*25.0;
    for(size_t i=0;i<dofs;++i)
    {
        QVERIFY2(near(response.internalForce[i],expectedForce[i],5e-11,3e-7),
            qPrintable(QString("plate workpoint internal force %1 actual=%2 expected=%3")
                .arg(qulonglong(i)).arg(response.internalForce[i],0,'g',17).arg(expectedForce[i],0,'g',17)));
        const double applied=(i/5==1 || i/5==2)?expectedForce[i]:0.0;
        QVERIFY2(near(adapter.data().supportReactions[i],expectedForce[i]-applied,5e-11,3e-7),
            qPrintable(QString("plate workpoint reaction %1 actual=%2 expected=%3")
                .arg(qulonglong(i)).arg(adapter.data().supportReactions[i],0,'g',17).arg(expectedForce[i]-applied,0,'g',17)));
    }
    QVERIFY(adapter.data().forceResidual<1e-10);
} catch(const RError &error) { QFAIL(qPrintable(error.getMessage())); }
  catch(const std::exception &error) { QFAIL(error.what()); }

void TestStructuralPiezoLayeredPlateOracle::openCircuitDrivenReciprocityAndTangent() try
{
    const PlateOracle expected=analyticPlate();
    const Vector q=patchDisplacement(1.0e-4,0.02);
    RModel drivenModel=plateScene(false);
    const RStructuralPiezoLayeredModelAdapter driven(drivenModel,drivenModel.getProblemSetup().getStructuralPiezoLayeredSetup());
    RModel openModel=plateScene(true);
    const RStructuralPiezoLayeredModelAdapter open(openModel,openModel.getProblemSetup().getStructuralPiezoLayeredSetup());
    const auto &drivenResponse=driven.electrical().referenceResponse();
    const auto &openResponse=open.electrical().referenceResponse();
    QVERIFY(near(voltage(drivenResponse),25.0,2e-13,2e-12));
    QVERIFY(near(charge(drivenResponse),3.75e-6,3e-12,2e-18));
    QVERIFY(near(voltage(openResponse),25.0,2e-13,2e-12));
    QVERIFY(near(charge(openResponse),3.75e-6,3e-12,2e-18));
    const auto &actualDriven=driven.electrical().condensed().stiffness;
    for(size_t i=0;i<dofs*dofs;++i)
        QVERIFY2(near(actualDriven[i],expected.stiffness[i],5e-11,4e-7),
            qPrintable(difference("Kdriven",i/dofs,i%dofs,actualDriven[i],expected.stiffness[i])));

    Vector expectedOpen=expected.stiffness;
    for(size_t i=0;i<dofs;++i) for(size_t j=0;j<dofs;++j)
        expectedOpen[i*dofs+j]+=expected.coupling[i]*expected.coupling[j]/expected.capacitance;
    const auto &actualOpen=open.electrical().condensed().stiffness;
    for(size_t i=0;i<dofs*dofs;++i)
        QVERIFY2(near(actualOpen[i],expectedOpen[i],5e-11,4e-7),
            qPrintable(difference("Kopen",i/dofs,i%dofs,actualOpen[i],expectedOpen[i])));

    // The same displacement and voltage yield the same equilibrium force in
    // driven and open-circuit forms; their tangents differ by GG^T/Ce.
    for(size_t i=0;i<dofs;++i)
        QVERIFY2(near(openResponse.internalForce[i],drivenResponse.internalForce[i],5e-11,4e-7),
            qPrintable(QString("plate workpoint force mismatch at %1").arg(qulonglong(i))));

    Vector direction(dofs,0.0);
    direction[0]=0.13; direction[2]=-0.17; direction[9]=0.21; direction[14]=-0.09; direction[19]=0.11;
    const double hq=2.0e-7,hv=1.0e-2;
    Vector qPlus=q,qMinus=q;
    for(size_t i=0;i<dofs;++i) {qPlus[i]+=hq*direction[i];qMinus[i]-=hq*direction[i];}
    const auto drivenPlus=driven.electrical().evaluate(q,{25.0+hv});
    const auto drivenMinus=driven.electrical().evaluate(q,{25.0-hv});
    for(size_t i=0;i<dofs;++i)
    {
        const double derivative=(drivenPlus.internalForce[i]-drivenMinus.internalForce[i])/(2.0*hv);
        QVERIFY2(near(derivative,expected.coupling[i],2e-8,2e-9),
            qPrintable(QString("dF/dV row %1 actual=%2 expected=%3").arg(qulonglong(i)).arg(derivative,0,'g',17).arg(expected.coupling[i],0,'g',17)));
    }
    const auto drivenQPlus=driven.electrical().evaluate(qPlus,{25.0});
    const auto drivenQMinus=driven.electrical().evaluate(qMinus,{25.0});
    const double expectedChargeDerivative=-dot(expected.coupling,direction);
    const double actualChargeDerivative=(charge(drivenQPlus)-charge(drivenQMinus))/(2.0*hq);
    QVERIFY2(near(actualChargeDerivative,expectedChargeDerivative,3e-7,2e-9),
        qPrintable(QString("dQ/dq virtual-work derivative actual=%1 expected=%2")
            .arg(actualChargeDerivative,0,'g',17).arg(expectedChargeDerivative,0,'g',17)));

    const auto openPlus=open.electrical().evaluate(qPlus,{});
    const auto openMinus=open.electrical().evaluate(qMinus,{});
    const Vector tangentExpected=multiply(expectedOpen,direction,dofs);
    for(size_t i=0;i<dofs;++i)
    {
        const double tangent=(openPlus.internalForce[i]-openMinus.internalForce[i])/(2.0*hq);
        QVERIFY2(near(tangent,tangentExpected[i],5e-7,3e-7),
            qPrintable(QString("open-circuit tangent row %1 actual=%2 expected=%3")
                .arg(qulonglong(i)).arg(tangent,0,'g',17).arg(tangentExpected[i],0,'g',17)));
        QVERIFY2(near(charge(openPlus),charge(openMinus),2e-10,2e-18),
            qPrintable(QString("open-circuit charge changed under virtual displacement: %1 vs %2")
                .arg(charge(openPlus),0,'g',17).arg(charge(openMinus),0,'g',17)));
    }
} catch(const RError &error) { QFAIL(qPrintable(error.getMessage())); }
  catch(const std::exception &error) { QFAIL(error.what()); }

QTEST_APPLESS_MAIN(TestStructuralPiezoLayeredPlateOracle)
#include "tst_structural_piezo_layered_plate_oracle.moc"

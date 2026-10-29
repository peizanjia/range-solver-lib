#include <QtTest>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <stdexcept>
#include <vector>

#include <rml_model.h>
#include <rml_boundary_condition.h>
#include <rml_environment_condition.h>
#include "rstructural_piezo_layered_model_adapter.h"
#include "rstructural_piezo_layered_native_run.h"

namespace
{
using Vector = std::vector<double>;
using Poly = std::array<long double,4>;
using ShapeRow = std::array<Poly,6>;
using Dof = RStructuralPiezoPhysicalDof;

struct LayerData
{
    double lower,upper,width,young,density,permittivity;
    const char *material;
    const char *port;
};

constexpr double length = 0.5;
constexpr double e31 = 5.0;
constexpr double epsilon33 = 1.5e-8;
constexpr std::array<LayerData,3> layers{{
    {-0.005, 0.005, 0.02, 70.0e9, 2700.0, 0.0, "host", ""},
    { 0.005, 0.007, 0.02, 60.0e9, 7500.0, epsilon33, "drive-pzt", "drive"},
    {-0.006,-0.005, 0.02, 60.0e9, 7500.0, epsilon33, "sense-pzt", "sense"}
}};

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

long double integrateProduct(const Poly &a,const Poly &b)
{
    long double result=0.0L;
    for(size_t i=0;i<a.size();++i)
        for(size_t j=0;j<b.size();++j)
            result+=a[i]*b[j]/static_cast<long double>(i+j+1);
    return result;
}

struct Section
{
    long double a=0.0L,b=0.0L,d=0.0L;
    long double i0=0.0L,i1=0.0L,i2=0.0L;
};

Section section()
{
    Section result;
    for(const auto &layer:layers)
    {
        const long double lo=layer.lower,hi=layer.upper;
        const long double t=hi-lo;
        result.a+=layer.width*layer.young*t;
        result.b+=layer.width*layer.young*(hi*hi-lo*lo)/2.0L;
        result.d+=layer.width*layer.young*(hi*hi*hi-lo*lo*lo)/3.0L;
        result.i0+=layer.width*layer.density*t;
        result.i1+=layer.width*layer.density*(hi*hi-lo*lo)/2.0L;
        result.i2+=layer.width*layer.density*(hi*hi*hi-lo*lo*lo)/3.0L;
    }
    return result;
}

struct BeamOracle
{
    Vector mass=Vector(36,0.0),stiffness=Vector(36,0.0);
    Vector coupling=Vector(12,0.0),capacitance=Vector(4,0.0);
    Section moments=section();
};

BeamOracle analyticBeam()
{
    BeamOracle result;
    ShapeRow axial{},transverse{},rotation{},axialDerivative{},rotationDerivative{};
    axial[0]={{1.0L,-1.0L,0.0L,0.0L}};
    axial[3]={{0.0L,1.0L,0.0L,0.0L}};
    transverse[1]={{1.0L,0.0L,-3.0L,2.0L}};
    transverse[2]={{0.0L,-length,2.0L*length,-length}};
    transverse[4]={{0.0L,0.0L,3.0L,-2.0L}};
    transverse[5]={{0.0L,0.0L,length,-length}};
    rotation[1]={{0.0L,6.0L/length,-6.0L/length,0.0L}};
    rotation[2]={{1.0L,-4.0L,3.0L,0.0L}};
    rotation[4]={{0.0L,-6.0L/length,6.0L/length,0.0L}};
    rotation[5]={{0.0L,-2.0L,3.0L,0.0L}};
    axialDerivative[0]={{-1.0L/length,0.0L,0.0L,0.0L}};
    axialDerivative[3]={{ 1.0L/length,0.0L,0.0L,0.0L}};
    for(size_t dof=0;dof<6;++dof)
        for(size_t power=1;power<4;++power)
            rotationDerivative[dof][power-1]=static_cast<long double>(power)*rotation[dof][power]/length;

    for(size_t i=0;i<6;++i)
        for(size_t j=0;j<6;++j)
        {
            const long double m=length*(result.moments.i0*(integrateProduct(axial[i],axial[j])+integrateProduct(transverse[i],transverse[j]))+
                result.moments.i1*(integrateProduct(axial[i],rotation[j])+integrateProduct(rotation[i],axial[j]))+
                result.moments.i2*integrateProduct(rotation[i],rotation[j]));
            const long double k=length*(result.moments.a*integrateProduct(axialDerivative[i],axialDerivative[j])+
                result.moments.b*(integrateProduct(axialDerivative[i],rotationDerivative[j])+integrateProduct(rotationDerivative[i],axialDerivative[j]))+
                result.moments.d*integrateProduct(rotationDerivative[i],rotationDerivative[j]));
            result.mass[i*6+j]=static_cast<double>(m);
            result.stiffness[i*6+j]=static_cast<double>(k);
        }

    // Integrating the independently derived Euler-Bernoulli strain row over
    // the element gives e*b*[delta(u)+zbar*delta(theta)] for each electrode.
    size_t port=0;
    for(const auto &layer:layers)
    {
        if(layer.port[0]=='\0') continue;
        const double center=0.5*(layer.lower+layer.upper);
        const double scale=e31*layer.width;
        result.coupling[0*2+port]=-scale;
        result.coupling[1*2+port]=0.0;
        result.coupling[2*2+port]=-scale*center;
        result.coupling[3*2+port]=scale;
        result.coupling[4*2+port]=0.0;
        result.coupling[5*2+port]=scale*center;
        const double capacitance=layer.permittivity*layer.width*length/(layer.upper-layer.lower);
        result.capacitance[port*2+port]=capacitance;
        ++port;
    }
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

RModel beamScene(double timeStep=1.0e-5,uint steps=8)
{
    RModel model;
    model.addNode(RNode(0.0,0.0,0.0));
    model.setNodePhysicalId(0,101);
    model.addNode(RNode(length,0.0,0.0));
    model.setNodePhysicalId(1,112);
    RElement element(R_ELEMENT_BEAM1);
    element.setNodeId(0,0); element.setNodeId(1,1);
    model.addElement(element,true,0);
    model.setElementPhysicalId(0,501);
    model.getLine(0).removeEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
    model.getLine(0).removeEnvironmentCondition(R_ENVIRONMENT_CONDITION_TEMPERATURE);
    model.setProblemTaskTree(RProblemTaskItem::generateDefaultTree(R_PROBLEM_STRESS));
    model.getProblemSetup().setMassMatrixType(RMassMatrixType::Consistent);
    auto &dynamics=model.getProblemSetup().getStructuralDynamicsSetup();
    dynamics.setEnabled(true); dynamics.setMethod(RStructuralTimeIntegrationMethod::Newmark);
    dynamics.setTimeStep(timeStep); dynamics.setNSteps(steps);
    model.getTimeSolver().setInputStartTime(timeStep);

    auto &setup=model.getProblemSetup().getStructuralPiezoLayeredSetup();
    setup.enabled=true;
    setup.family=RStructuralPiezoLayeredModelFamily::PlanarBeam;
    setup.materials={material("host",70.0e9,2700.0,false),
                     material("drive-pzt",60.0e9,7500.0,true),
                     material("sense-pzt",60.0e9,7500.0,true)};
    RStructuralPiezoLayeredModelLayer host; host.materialId="host"; host.lowerZ=-0.005; host.upperZ=0.005; host.width=0.02;
    RStructuralPiezoLayeredModelLayer drive; drive.materialId="drive-pzt"; drive.lowerZ=0.005; drive.upperZ=0.007; drive.width=0.02; drive.portId="drive";
    RStructuralPiezoLayeredModelLayer sense; sense.materialId="sense-pzt"; sense.lowerZ=-0.006; sense.upperZ=-0.005; sense.width=0.02; sense.portId="sense";
    setup.members={{501,{host,drive,sense}}};
    RStructuralPiezoLayeredPort actuator; actuator.readout.portId="drive";
    actuator.readout.role=RStructuralPiezoPortRole::ActuatorSensor;
    actuator.readout.readout=RStructuralPiezoReadoutMode::ActuatorCharge;
    actuator.initialValue=25.0;
    RStructuralPiezoLayeredPort sensor; sensor.readout.portId="sense";
    sensor.readout.role=RStructuralPiezoPortRole::Sensor;
    sensor.readout.readout=RStructuralPiezoReadoutMode::OpenCircuitVoltage;
    sensor.initialValue=-3.25e-6;
    setup.ports={actuator,sensor};
    setup.supports={{101,Dof::UX,0.001},{101,Dof::UZ,0.0},{101,Dof::RY,0.0},
                    {112,Dof::UZ,0.0},{112,Dof::RY,0.0}};
    setup.backgroundForces={{112,{{353.0,0.0,0.0,0.0,0.0,0.0}}}};
    setup.workpoint.available=true;
    RStructuralPiezoLayeredNodeWorkpoint base; base.physicalNodeId=101; base.q0[0]=0.001;
    RStructuralPiezoLayeredNodeWorkpoint tip; tip.physicalNodeId=112; tip.q0[0]=0.00101; tip.f0[0]=353.0;
    setup.workpoint.nodes={base,tip};
    setup.workpoint.ports={{"drive",25.0,0.875e-6},{"sense",-15.0,-3.25e-6}};
    setup.bindToModel(model);
    return model;
}

RMatrixSolverConf accurateSolver()
{
    RMatrixSolverConf conf(RMatrixSolverConf::GMRES);
    conf.setNInnerIterations(100); conf.setNOuterIterations(100); conf.setSolverCvgValue(1.0e-13);
    return conf;
}

struct ScalarState { double q=0.0,v=0.0,a=0.0; };

ScalarState scalarNewmarkStep(ScalarState state,double force,double mass,double stiffness,double dt)
{
    const double qPredict=state.q+dt*state.v+0.25*dt*dt*state.a;
    const double vPredict=state.v+0.5*dt*state.a;
    const double acceleration=(force-stiffness*qPredict)/(mass+0.25*dt*dt*stiffness);
    return {qPredict+0.25*dt*dt*acceleration,vPredict+0.5*dt*acceleration,acceleration};
}

struct ByPort
{
    double voltage=0.0,charge=0.0;
};

ByPort port(const RStructuralPiezoElectricalResponse &response,const QString &id)
{
    const auto found=std::find_if(response.ports.begin(),response.ports.end(),
        [&](const auto &candidate){return candidate.portId==id;});
    if(found==response.ports.end()) throw std::runtime_error("Native beam response omitted an expected electrical port");
    return {found->voltage,found->charge};
}
}

class TestStructuralPiezoLayeredBeamOracle : public QObject
{
    Q_OBJECT
private slots:
    void independentSectionMatricesAndCertifiedWorkpoint();
    void acceptedVoltageStepMatchesScalarNewmark();
};

void TestStructuralPiezoLayeredBeamOracle::independentSectionMatricesAndCertifiedWorkpoint()
{
    const BeamOracle expected=analyticBeam();
    QVERIFY(near(static_cast<double>(expected.moments.a),17.6e6,2e-14,1e-8));
    QVERIFY(near(static_cast<double>(expected.moments.b),7800.0,2e-14,1e-12));
    QVERIFY(near(static_cast<double>(expected.moments.d),240.26666666666665,2e-14,1e-12));
    QVERIFY(near(static_cast<double>(expected.moments.i0),0.99,2e-14,1e-14));
    QVERIFY(near(static_cast<double>(expected.moments.i1),0.000975,2e-14,1e-16));
    QVERIFY(near(static_cast<double>(expected.moments.i2),1.995e-5,2e-14,1e-18));

    const RModel model=beamScene();
    const auto &setup=model.getProblemSetup().getStructuralPiezoLayeredSetup();
    const RStructuralPiezoLayeredModelAdapter adapter(model,setup);
    const auto &actual=adapter.mesh().operators;
    QCOMPARE(actual.dofs,6u);
    QCOMPARE(actual.portIds,(std::vector<QString>{"drive","sense"}));
    for(size_t i=0;i<36;++i)
    {
        const size_t row=i/6,column=i%6;
        QVERIFY2(near(actual.mass[i],expected.mass[i],3e-12,2e-15),qPrintable(difference("M",row,column,actual.mass[i],expected.mass[i])));
        QVERIFY2(near(actual.stiffness[i],expected.stiffness[i],3e-12,2e-8),qPrintable(difference("K",row,column,actual.stiffness[i],expected.stiffness[i])));
    }
    for(size_t i=0;i<12;++i)
    {
        const double value=actual.coupling[i];
        QVERIFY2(near(value,expected.coupling[i],3e-12,2e-14),qPrintable(difference("G",i/2,i%2,value,expected.coupling[i])));
    }
    for(size_t i=0;i<4;++i)
        QVERIFY2(near(actual.capacitance[i],expected.capacitance[i],3e-12,2e-20),
            qPrintable(difference("Ce",i/2,i%2,actual.capacitance[i],expected.capacitance[i])));

    // A rigid translation plus section rotation must project through the
    // first and second mass moments with the declared positive UX-RY sign.
    const double u=0.17,w=-0.23,theta=0.31;
    const std::array<double,6> rigid{{u,w,theta,u,w-length*theta,theta}};
    long double projected=0.0L;
    for(size_t i=0;i<6;++i) for(size_t j=0;j<6;++j)
        projected+=static_cast<long double>(rigid[i])*expected.mass[i*6+j]*rigid[j];
    const long double zbar=expected.moments.i1/expected.moments.i0;
    const long double zsecond=expected.moments.i2/expected.moments.i0;
    const long double rigidExpected=expected.moments.i0*length*
        (u*u+2.0L*zbar*u*theta+zsecond*theta*theta+w*w-length*w*theta+length*length*theta*theta/3.0L);
    QVERIFY(near(static_cast<double>(projected),static_cast<double>(rigidExpected),2e-12,2e-15));
    QVERIFY(expected.mass[0*6+2]>0.0); // positive I1 UX-RY coupling
    QVERIFY(expected.mass[1*6+2]<0.0); // transverse rigid rotation is w'=-RY

    // A constant axial strain and curvature provide an independent patch
    // force from the section resultants N=A*eps+B*kappa and M=B*eps+D*kappa.
    const double epsilon=2.0e-5,curvature=3.0e-4;
    const std::array<double,6> patch{{0,0,0,epsilon*length,-curvature*length*length/2.0,curvature*length}};
    const double resultant=static_cast<double>(expected.moments.a*epsilon+expected.moments.b*curvature);
    const double moment=static_cast<double>(expected.moments.b*epsilon+expected.moments.d*curvature);
    for(size_t i=0;i<6;++i)
    {
        double force=0.0;
        for(size_t j=0;j<6;++j) force+=expected.stiffness[i*6+j]*patch[j];
        const double analytic=i==0?-resultant:(i==2?-moment:(i==3?resultant:(i==5?moment:0.0)));
        QVERIFY2(near(force,analytic,4e-12,2e-8),qPrintable(QString("beam patch force row %1 actual=%2 expected=%3")
            .arg(qulonglong(i)).arg(force,0,'g',17).arg(analytic,0,'g',17)));
    }

    Vector openStiffness=expected.stiffness;
    const size_t senseColumn=1;
    for(size_t i=0;i<6;++i) for(size_t j=0;j<6;++j)
        openStiffness[i*6+j]+=expected.coupling[i*2+senseColumn]*expected.coupling[j*2+senseColumn]/expected.capacitance[3];
    const auto &actualOpen=adapter.electrical().condensed().stiffness;
    for(size_t i=0;i<36;++i)
        QVERIFY2(near(actualOpen[i],openStiffness[i],4e-12,2e-8),
            qPrintable(difference("Kopen",i/6,i%6,actualOpen[i],openStiffness[i])));

    const auto &response=adapter.electrical().referenceResponse();
    const auto drive=port(response,"drive");
    const auto sense=port(response,"sense");
    QVERIFY(near(drive.voltage,25.0,1e-13,1e-12));
    QVERIFY(near(drive.charge,0.875e-6,2e-12,2e-18));
    QVERIFY(near(sense.voltage,-15.0,1e-13,1e-12));
    QVERIFY(near(sense.charge,-3.25e-6,2e-12,2e-18));
    const std::array<double,6> q0{{0.001,0,0,0.00101,0,0}};
    for(size_t i=0;i<6;++i)
    {
        double force=0.0;
        for(size_t j=0;j<6;++j) force+=expected.stiffness[i*6+j]*q0[j];
        for(size_t p=0;p<2;++p) force+=expected.coupling[i*2+p]*(p==0?25.0:-15.0);
        const double applied=i==3?353.0:0.0;
        QVERIFY2(near(response.internalForce[i],force,3e-12,2e-8),
            qPrintable(QString("workpoint internal force row %1 actual=%2 expected=%3")
                .arg(qulonglong(i)).arg(response.internalForce[i],0,'g',17).arg(force,0,'g',17)));
        QVERIFY2(near(adapter.data().supportReactions[i],force-applied,3e-12,2e-8),
            qPrintable(QString("workpoint reaction row %1 actual=%2 expected=%3")
                .arg(qulonglong(i)).arg(adapter.data().supportReactions[i],0,'g',17).arg(force-applied,0,'g',17)));
    }
    QVERIFY(near(response.internalForce[3],353.0,2e-12,2e-8));
    QVERIFY(near(adapter.data().supportReactions[0],-353.0,2e-12,2e-8));

    const auto &compiled=adapter.data();
    QCOMPARE(compiled.freeDofs,1u);
    const size_t free=compiled.freeDofs;
    size_t tipUxColumn=free;
    for(size_t j=0;j<free;++j)
        if(near(compiled.T[3*free+j],1.0,0.0,1e-14)) tipUxColumn=j;
    QVERIFY(tipUxColumn<free);
    const double sensorG=expected.coupling[3*2+senseColumn];
    const double expectedK=static_cast<double>(expected.moments.a/length)+sensorG*sensorG/expected.capacitance[3];
    const double expectedM=static_cast<double>(expected.moments.i0*length/3.0L);
    QVERIFY2(near(compiled.mass[tipUxColumn*free+tipUxColumn],expectedM,3e-12,2e-15),
        qPrintable(QString("supported axial mass actual=%1 expected=%2").arg(compiled.mass[tipUxColumn*free+tipUxColumn],0,'g',17).arg(expectedM,0,'g',17)));
    QVERIFY2(near(compiled.stiffness[tipUxColumn*free+tipUxColumn],expectedK,3e-12,2e-8),
        qPrintable(QString("supported open-circuit axial stiffness actual=%1 expected=%2").arg(compiled.stiffness[tipUxColumn*free+tipUxColumn],0,'g',17).arg(expectedK,0,'g',17)));
    QVERIFY(near(expectedM,0.165,2e-14,1e-15));
    QVERIFY(near(expectedK,35266666.666666664,2e-14,1e-7));
}

void TestStructuralPiezoLayeredBeamOracle::acceptedVoltageStepMatchesScalarNewmark()
{
    constexpr double dt=1.0e-5;
    RModel source=beamScene(dt,8);
    auto &setup=source.getProblemSetup().getStructuralPiezoLayeredSetup();
    setup.controllerConfiguration="independent-fixed-voltage-step/v1";
    setup.bindToModel(source);
    RStructuralPiezoLayeredRuntimeCallbacks callbacks;
    callbacks.sourceFingerprint=setup.computeSourceFingerprint(source);
    callbacks.controllerConfiguration=setup.controllerConfiguration;
    callbacks.backendIdentity="independent-analytic-oracle/v1";
    callbacks.controllerIdentity="step-to-125V/v1";
    callbacks.controlEveryFESteps=1;
    callbacks.initialOpaqueControllerState="initial";
    callbacks.propose=[](const RStructuralPiezoLayeredDynamicsState &,uint,uint,const QByteArray &)
    {
        RStructuralPiezoLayeredControllerProposal proposal;
        proposal.requestedVoltages={125.0};
        proposal.nextOpaqueState="accepted-voltage-step";
        return proposal;
    };
    callbacks.validateOpaque=[](const QByteArray &,uint,uint,uint) {};
    RStructuralPiezoLayeredNativeRun run(source,callbacks,accurateSolver());
    const BeamOracle expected=analyticBeam();
    const double mass=static_cast<double>(expected.moments.i0*length/3.0L);
    const double senseG=expected.coupling[3*2+1];
    const double stiffness=static_cast<double>(expected.moments.a/length)+senseG*senseG/expected.capacitance[3];
    constexpr double voltageIncrement=100.0;
    const double force=-expected.coupling[3*2]*voltageIncrement;
    ScalarState oracle; oracle.a=force/mass; // refresh at the left endpoint after the held-voltage step.
    QCOMPARE(run.completedSteps(),0u);
    for(uint stepIndex=1;stepIndex<=8;++stepIndex)
    {
        oracle=scalarNewmarkStep(oracle,force,mass,stiffness,dt);
        auto prepared=run.prepare(source);
        QCOMPARE(prepared->completedSteps(),stepIndex);
        QCOMPARE(prepared->result().heldActuatorVoltages,(std::vector<double>{125.0}));
        run.validatePrepared(*prepared);
        run.commit(*prepared);
        const auto &state=run.acceptedState();
        QCOMPARE(run.completedSteps(),stepIndex);
        QVERIFY2(near(state.physicalDisplacements[0],0.001,1e-12,2e-15),"base UX support moved during the native step");
        QVERIFY2(near(state.physicalDisplacements[6],0.00101+oracle.q,3e-8,2e-14),
            qPrintable(QString("step %1 tip UX actual=%2 expected=%3").arg(stepIndex).arg(state.physicalDisplacements[6],0,'g',17).arg(0.00101+oracle.q,0,'g',17)));
        QVERIFY2(near(state.physicalVelocities[6],oracle.v,3e-8,2e-12),
            qPrintable(QString("step %1 tip velocity actual=%2 expected=%3").arg(stepIndex).arg(state.physicalVelocities[6],0,'g',17).arg(oracle.v,0,'g',17)));
        QVERIFY2(near(state.physicalAccelerations[6],oracle.a,3e-8,2e-9),
            qPrintable(QString("step %1 tip acceleration actual=%2 expected=%3").arg(stepIndex).arg(state.physicalAccelerations[6],0,'g',17).arg(oracle.a,0,'g',17)));
        for(size_t i=0;i<state.physicalDisplacements.size();++i)
            if(i!=0 && i!=6) QVERIFY2(std::abs(state.physicalDisplacements[i])<2e-12,
                qPrintable(QString("unexpected non-axial displacement component %1=%2").arg(qulonglong(i)).arg(state.physicalDisplacements[i],0,'g',17)));

        const double extension=1.0e-5+oracle.q;
        const double expectedDriveCharge=expected.capacitance[0]*125.0-expected.coupling[3*2]*extension;
        const double expectedSenseVoltage=(-3.25e-6+expected.coupling[3*2+1]*extension)/expected.capacitance[3];
        const auto drive=port(state.electrical,"drive");
        const auto sense=port(state.electrical,"sense");
        QVERIFY2(near(drive.voltage,125.0,1e-13,1e-12),qPrintable(QString("step %1 drive voltage=%2").arg(stepIndex).arg(drive.voltage,0,'g',17)));
        QVERIFY2(near(drive.charge,expectedDriveCharge,3e-8,2e-17),
            qPrintable(QString("step %1 drive charge actual=%2 expected=%3").arg(stepIndex).arg(drive.charge,0,'g',17).arg(expectedDriveCharge,0,'g',17)));
        QVERIFY2(near(sense.charge,-3.25e-6,2e-12,2e-18),qPrintable(QString("step %1 open-sensor charge=%2").arg(stepIndex).arg(sense.charge,0,'g',17)));
        QVERIFY2(near(sense.voltage,expectedSenseVoltage,3e-8,2e-11),
            qPrintable(QString("step %1 open-sensor voltage actual=%2 expected=%3").arg(stepIndex).arg(sense.voltage,0,'g',17).arg(expectedSenseVoltage,0,'g',17)));

        const double expectedBaseReaction=-353.0+(expected.moments.i0*length/6.0)*oracle.a-stiffness*oracle.q+force;
        QVERIFY2(near(state.physicalReactions[0],expectedBaseReaction,4e-8,2e-7),
            qPrintable(QString("step %1 base UX reaction actual=%2 expected=%3").arg(stepIndex).arg(state.physicalReactions[0],0,'g',17).arg(expectedBaseReaction,0,'g',17)));
        QVERIFY2(state.freeForceResidual<1e-7,
            qPrintable(QString("step %1 free-force residual=%2").arg(stepIndex).arg(state.freeForceResidual,0,'g',17)));
    }
}

QTEST_APPLESS_MAIN(TestStructuralPiezoLayeredBeamOracle)
#include "tst_structural_piezo_layered_beam_oracle.moc"

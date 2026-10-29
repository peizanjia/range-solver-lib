#include <QtTest>
#include <QTemporaryDir>
#include <rml_model.h>
#include <rbl_error.h>
#include "rstructural_piezo_solid_kernel.h"
#include "rstructural_piezo_layered_native_run.h"
#include "rstructural_piezo_layered_restart_codec.h"
#include <algorithm>
#include <cmath>
#include <functional>

namespace {
using V=std::vector<double>;using V3=std::array<double,3>;
bool close(double a,double b,double tol=2e-9){return std::abs(a-b)<=tol*std::max({1.,std::abs(a),std::abs(b)});}
bool rejects(const std::function<void()>&f){try{f();}catch(...){return true;}return false;}
double quadratic(const V&a,const V&x){double v=0;for(size_t i=0;i<x.size();++i)for(size_t j=0;j<x.size();++j)v+=x[i]*a[i*x.size()+j]*x[j];return v;}
V multiply(const V&a,const V&x){V y(x.size());for(size_t i=0;i<x.size();++i)for(size_t j=0;j<x.size();++j)y[i]+=a[i*x.size()+j]*x[j];return y;}
RStructuralPiezoContinuumMaterial material(bool active=true){auto m=RStructuralPiezoContinuumMaterial::isotropic(120,.23,7);if(active){m.piezoelectric=true;m.coupling[12]=2;m.coupling[13]=1;m.coupling[14]=3;m.permittivity={{4,0,0,0,5,0,0,0,6}};}return m;}
struct Grid{std::vector<RStructuralPiezoLayeredNode>nodes;std::vector<RStructuralPiezoLayeredMember>members;};
Grid grid(bool tet,uint nx=1,uint ny=1,uint nz=1,bool active=true){
    Grid g;const auto id=[&](uint i,uint j,uint k){return 101+i+(nx+1)*(j+(ny+1)*k);};
    for(uint k=0;k<=nz;++k)for(uint j=0;j<=ny;++j)for(uint i=0;i<=nx;++i)g.nodes.push_back({id(i,j,k),{{double(i)/nx,double(j)/ny,double(k)/nz}}});
    RStructuralPiezoLayer layer;layer.material=material(active);layer.lowerZ=0;layer.upperZ=1;if(active)layer.portId="drive";
    constexpr uint pattern[6][4]={{0,1,2,6},{0,2,3,6},{0,3,7,6},{0,7,4,6},{0,4,5,6},{0,5,1,6}};
    for(uint k=0;k<nz;++k)for(uint j=0;j<ny;++j)for(uint i=0;i<nx;++i){
        const uint c[8]={id(i,j,k),id(i+1,j,k),id(i+1,j+1,k),id(i,j+1,k),id(i,j,k+1),id(i+1,j,k+1),id(i+1,j+1,k+1),id(i,j+1,k+1)};
        if(tet)for(const auto&p:pattern){RStructuralPiezoLayeredMember member;member.physicalId=501+uint(g.members.size());for(uint n:p)member.nodeIds.push_back(c[n]);member.layers={layer};g.members.push_back(member);}
        else g.members.push_back({501+uint(g.members.size()),std::vector<uint>(c,c+8),{layer}});
    }return g;
}
RStructuralPiezoLayeredMeshOperators assemble(const Grid&g){return RStructuralPiezoLayeredMesh::assemble(RStructuralPiezoLayeredFamily::Solid,{},g.nodes,g.members);}
RModel scene(bool tet,RStructuralPiezoReadoutMode mode=RStructuralPiezoReadoutMode::ActuatorCharge){
    auto g=grid(tet,1,1,2);RModel m;
    for(const auto&n:g.nodes){m.addNode(RNode(n.position[0],n.position[1],n.position[2]));m.setNodePhysicalId(m.getNNodes()-1,n.physicalId);}
    for(const auto&member:g.members){RElement e(tet?R_ELEMENT_TETRA1:R_ELEMENT_HEXA1);for(uint j=0;j<member.nodeIds.size();++j)e.setNodeId(j,m.findNodePositionByPhysicalId(member.nodeIds[j]));m.addElement(e,true,0);m.setElementPhysicalId(m.getNElements()-1,member.physicalId);}
    for(uint i=0;i<m.getNVolumes();++i){m.getVolume(i).removeEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION);m.getVolume(i).removeEnvironmentCondition(R_ENVIRONMENT_CONDITION_TEMPERATURE);}
    m.setProblemTaskTree(RProblemTaskItem::generateDefaultTree(R_PROBLEM_STRESS));m.getProblemSetup().setMassMatrixType(RMassMatrixType::Consistent);
    auto&d=m.getProblemSetup().getStructuralDynamicsSetup();d.setEnabled(true);d.setTimeStep(.01);d.setNSteps(6);m.getTimeSolver().setInputStartTime(.01);
    auto&s=m.getProblemSetup().getStructuralPiezoLayeredSetup();s.enabled=true;s.family=RStructuralPiezoLayeredModelFamily::Solid3D;
    auto native=material();RStructuralPiezoLayeredMaterial mat;mat.materialId="pzt";mat.elasticityOrCompliance=native.elasticity;mat.coupling=native.coupling;mat.permittivity=native.permittivity;mat.density=native.density;mat.piezoelectric=true;mat.sourceId="independent-solid-test";mat.sourceSha256=QString(64,'b');s.materials={mat};
    for(const auto&member:g.members){RStructuralPiezoLayeredModelLayer l;l.materialId=mat.materialId;l.lowerZ=0;l.upperZ=1;l.portId="drive";s.members.push_back({member.physicalId,{l}});}
    RStructuralPiezoLayeredPort port;port.readout.portId="drive";port.readout.readout=mode;port.readout.role=mode==RStructuralPiezoReadoutMode::ActuatorCharge?RStructuralPiezoPortRole::ActuatorSensor:RStructuralPiezoPortRole::Sensor;
    if(mode==RStructuralPiezoReadoutMode::IdealChargeAmplifier)port.readout.Cf=.5;s.ports={port};
    for(const auto&node:g.nodes)if(node.position[0]==0)for(auto dof:{RStructuralPiezoPhysicalDof::UX,RStructuralPiezoPhysicalDof::UY,RStructuralPiezoPhysicalDof::UZ})s.supports.push_back({node.physicalId,dof,0});
    s.initialVelocities={{g.nodes.back().physicalId,{{.05,.02,-.03,0,0,0}}}};
    s.damping.kind=RStructuralPiezoDampingKind::Rayleigh;s.damping.rayleighAlpha=.02;s.damping.rayleighBeta=.001;s.bindToModel(m);return m;
}
}
#define CATCH_FAILURE catch(const RError&e){QFAIL(qPrintable(e.getMessage()));}catch(const std::exception&e){QFAIL(e.what());}
class TestStructuralPiezoSolid:public QObject{
    Q_OBJECT
private slots:
    void mechanicalPatchAndRigidBody_data();void mechanicalPatchAndRigidBody();
    void gaugeReciprocityAndCondensedElectrodes_data();void gaugeReciprocityAndCondensedElectrodes();
    void meshAdmission();
    void nativeDynamicsReadoutsAndRestart_data();void nativeDynamicsReadoutsAndRestart();
    void rotatedMaterialFrameAndPolarity();
};
void TestStructuralPiezoSolid::mechanicalPatchAndRigidBody_data(){QTest::addColumn<bool>("tet");QTest::newRow("TET4")<<true;QTest::newRow("HEX8")<<false;}
void TestStructuralPiezoSolid::mechanicalPatchAndRigidBody()try{
    QFETCH(bool,tet);auto g=grid(tet,2,2,2,false);auto mesh=assemble(g);const auto&op=mesh.operators;QCOMPARE(op.dofs,uint(g.nodes.size()*3));QVERIFY(close(op.totalMass,7));
    const auto refined=assemble(grid(tet,4,1,2,true));for(size_t i=0;i<refined.operators.dofs;++i)for(size_t j=i+1;j<refined.operators.dofs;++j){QCOMPARE(refined.operators.mass[i*refined.operators.dofs+j],refined.operators.mass[j*refined.operators.dofs+i]);QCOMPARE(refined.operators.stiffness[i*refined.operators.dofs+j],refined.operators.stiffness[j*refined.operators.dofs+i]);}
    const double a[3][3]={{.03,.02,-.01},{-.01,.04,.05},{.02,-.02,-.06}};
    V q(op.dofs);for(size_t i=0;i<g.nodes.size();++i)for(size_t j=0;j<3;++j)for(size_t k=0;k<3;++k)q[3*i+j]+=a[j][k]*g.nodes[i].position[k];
    const V strain={a[0][0],a[1][1],a[2][2],a[1][2]+a[2][1],a[0][2]+a[2][0],a[0][1]+a[1][0]};const auto m=material(false);QVERIFY(close(quadratic(op.stiffness,q),quadratic(V(m.elasticity.begin(),m.elasticity.end()),strain)));
    const auto f=multiply(op.stiffness,q);for(size_t i=0;i<g.nodes.size();++i)if(g.nodes[i].position==V3{{.5,.5,.5}})for(size_t k=0;k<3;++k)QVERIFY(std::abs(f[3*i+k])<1e-12);
    for(uint rigid=0;rigid<6;++rigid){V v(op.dofs);for(size_t i=0;i<g.nodes.size();++i){const auto&x=g.nodes[i].position;if(rigid<3)v[3*i+rigid]=1;else{const size_t axis=rigid-3;v[3*i+(axis+1)%3]=-x[(axis+2)%3];v[3*i+(axis+2)%3]=x[(axis+1)%3];}}
        const auto r=multiply(op.stiffness,v);for(double x:r)QVERIFY(std::abs(x)<1e-11);if(rigid<3)QVERIFY(close(quadratic(op.mass,v),7));}
    // Affine TET4 and HEX8 have consistent, non-diagonal positive mass.
    QVERIFY(op.mass[3]>0);const auto physical=mesh.reconstructPhysical(q);QCOMPARE(mesh.projectPhysical(physical,RStructuralPiezoPhysicalDomain::Displacement),q);
    // Independent non-diagonal Jacobian check: shear/skew the actual element,
    // prescribe the same GLOBAL affine strain, and use det(L)*reference volume.
    const double map[3][3]={{1.2,.3,.1},{.2,.9,-.15},{-.1,.25,1.4}};
    const double det=map[0][0]*(map[1][1]*map[2][2]-map[1][2]*map[2][1])-map[0][1]*(map[1][0]*map[2][2]-map[1][2]*map[2][0])+map[0][2]*(map[1][0]*map[2][1]-map[1][1]*map[2][0]);
    std::vector<V3> skew;for(uint id:g.members.front().nodeIds){V3 x{};const auto&source=g.nodes[id-101].position;for(size_t i=0;i<3;++i)for(size_t j=0;j<3;++j)x[i]+=map[i][j]*source[j];skew.push_back(x);}
    const auto raw=RStructuralPiezoSolidKernel::assemble(skew,m);V affine(raw.mechanical.dofs);for(size_t i=0;i<skew.size();++i)for(size_t j=0;j<3;++j)for(size_t k=0;k<3;++k)affine[3*i+j]+=a[j][k]*skew[i][k];
    const double volume=det/(tet?48.:8.);QVERIFY(close(raw.mechanical.totalMass,7*volume));QVERIFY(close(quadratic(raw.mechanical.stiffness,affine),volume*quadratic(V(m.elasticity.begin(),m.elasticity.end()),strain)));
}CATCH_FAILURE
void TestStructuralPiezoSolid::gaugeReciprocityAndCondensedElectrodes_data(){mechanicalPatchAndRigidBody_data();}
void TestStructuralPiezoSolid::gaugeReciprocityAndCondensedElectrodes()try{
    QFETCH(bool,tet);auto g=grid(tet,1,1,2);auto mesh=assemble(g);const auto&op=mesh.operators;QCOMPARE(op.portIds,std::vector<QString>{"drive"});QVERIFY(close(op.capacitance[0],6));
    V q(op.dofs);for(size_t i=0;i<g.nodes.size();++i)q[3*i+2]=.1*g.nodes[i].position[2];
    double gtq=0;for(size_t i=0;i<q.size();++i)gtq+=op.coupling[i]*q[i];QVERIFY(close(gtq,.3));
    const auto&member=g.members.front();std::vector<V3>x;for(uint id:member.nodeIds)x.push_back(g.nodes[id-101].position);const auto raw=RStructuralPiezoSolidKernel::assemble(x,material());const size_t nn=x.size();
    for(size_t i=0;i<nn;++i){double sum=0;for(size_t j=0;j<nn;++j)sum+=raw.dielectric[i*nn+j];QVERIFY(std::abs(sum)<1e-12);}
    for(size_t i=0;i<3*nn;++i){double sum=0;for(size_t j=0;j<nn;++j)sum+=raw.coupling[i*nn+j];QVERIFY(std::abs(sum)<1e-12);}
    RStructuralPiezoElectricalSystem electric(op,{{"drive",RStructuralPiezoElectricalRole::VoltageActuator}},V(op.dofs),{0});
    const double v=.4;const auto response=electric.evaluate(q,{v});QVERIFY(close(response.ports[0].charge,6*v-.3));
    // Mixed second derivatives of enthalpy carry opposite signs in force/Q.
    const double step=1e-5;V dq=q;dq.back()+=step;auto shifted=electric.evaluate(dq,{v});auto voltage=electric.evaluate(q,{v+step});
    QVERIFY(close((voltage.internalForce.back()-response.internalForce.back())/step,-(shifted.ports[0].charge-response.ports[0].charge)/step));
    for(size_t i=0;i<op.dofs;++i)for(size_t j=0;j<op.dofs;++j)QVERIFY(close(op.stiffness[i*op.dofs+j],op.stiffness[j*op.dofs+i]));
    // Condensing a free middle electrode plane must stiffen nonuniform strain.
    auto passive=g;for(auto&member:passive.members){member.layers[0].material=material(false);member.layers[0].portId.clear();}const auto bare=assemble(passive);
    V nonuniform(op.dofs);for(size_t i=0;i<g.nodes.size();++i)if(g.nodes[i].position[2]==.5)nonuniform[3*i+2]=.1;
    QVERIFY(quadratic(op.stiffness,nonuniform)>quadratic(bare.operators.stiffness,nonuniform)+1e-5);
}CATCH_FAILURE
void TestStructuralPiezoSolid::meshAdmission()try{
    auto g=grid(false);auto bad=g;bad.members.push_back(bad.members[0]);bad.members.back().physicalId+=10;QVERIFY(rejects([&]{assemble(bad);}));
    bad=g;bad.nodes[6].position[2]+=.1;QVERIFY(rejects([&]{assemble(bad);}));
    bad=grid(true);std::swap(bad.members[0].nodeIds[1],bad.members[0].nodeIds[2]);QVERIFY(rejects([&]{assemble(bad);}));
    bad=g;bad.members[0].layers[0].upperZ=2;QVERIFY(rejects([&]{assemble(bad);}));
    bad=g;bad.members[0].layers[0].coverage[1]=0;QVERIFY(rejects([&]{assemble(bad);}));
    bad=g;bad.members[0].layers.push_back(bad.members[0].layers[0]);QVERIFY(rejects([&]{assemble(bad);}));
    auto mesh=assemble(g);V unsupported(g.nodes.size()*6);unsupported[3]=1;QVERIFY(rejects([&]{mesh.projectPhysical(unsupported,RStructuralPiezoPhysicalDomain::ForceAndMoment);}));
    auto model=scene(false);auto&s=model.getProblemSetup().getStructuralPiezoLayeredSetup();s.supports.clear();s.sourceFingerprint.clear();QVERIFY(rejects([&]{RStructuralPiezoLayeredModelAdapter invalid(model,s);}));
}CATCH_FAILURE
void TestStructuralPiezoSolid::nativeDynamicsReadoutsAndRestart_data(){QTest::addColumn<bool>("tet");QTest::addColumn<int>("mode");for(bool tet:{false,true})for(auto mode:{RStructuralPiezoReadoutMode::ActuatorCharge,RStructuralPiezoReadoutMode::OpenCircuitVoltage,RStructuralPiezoReadoutMode::ShortCircuitCharge,RStructuralPiezoReadoutMode::IdealChargeAmplifier})QTest::newRow(qPrintable(QString("%1-mode%2").arg(tet?"tet":"hex").arg(int(mode))))<<tet<<int(mode);}
void TestStructuralPiezoSolid::nativeDynamicsReadoutsAndRestart()try{
    QFETCH(bool,tet);QFETCH(int,mode);auto model=scene(tet,RStructuralPiezoReadoutMode(mode));RStructuralPiezoLayeredNativeRun run(model);QCOMPARE(run.adapter().data().fullDofs,3*model.getNNodes());
    for(const auto&unit:run.adapter().data().generalizedCoordinateUnits)QCOMPARE(unit,QString("m"));const auto initial=run.acceptedState().physicalDisplacements;
    auto prepared=run.prepare(model);run.validatePrepared(*prepared);QCOMPARE(run.completedSteps(),0u);run.commit(*prepared);QCOMPARE(run.completedSteps(),1u);QVERIFY(run.acceptedState().physicalDisplacements!=initial);
    const auto&s=run.acceptedState();QVERIFY(std::abs(s.energyBalanceResidual)<1e-8);QVERIFY(s.freeForceResidual<1e-8);QVERIFY(!s.readouts.empty());
    for(size_t i=0;i<model.getNNodes();++i)for(size_t k=3;k<6;++k){QCOMPARE(s.physicalDisplacements[6*i+k],0.);QCOMPARE(s.physicalVelocities[6*i+k],0.);}
    const auto&e=s.electrical.ports[0];if(mode==int(RStructuralPiezoReadoutMode::OpenCircuitVoltage))QCOMPARE(e.charge,0.);else if(mode!=int(RStructuralPiezoReadoutMode::ActuatorCharge))QCOMPARE(e.voltage,0.);
    if(mode==int(RStructuralPiezoReadoutMode::OpenCircuitVoltage))QVERIFY(std::abs(e.voltage)>1e-12);else QVERIFY(std::abs(e.charge)>1e-12);
    if(mode==int(RStructuralPiezoReadoutMode::IdealChargeAmplifier))QVERIFY(close(s.readouts[0].value,-e.charge/.5));
    if(mode==int(RStructuralPiezoReadoutMode::ActuatorCharge)){
        RStructuralPiezoLayeredDynamics dynamics(run.adapter(),RMatrixSolverConf(RMatrixSolverConf::GMRES),0,[count=model.getNNodes()](double){return V(6*count);},"solid-test/zero-force");auto zero=dynamics.initialize();auto step=dynamics.trial(zero,{.2},.01);auto unforced=dynamics.trial(zero,{0},.01);QVERIFY(step.physicalDisplacements!=unforced.physicalDisplacements);QVERIFY(std::abs(step.energyBalanceResidual)<1e-8);
    }
    auto checkpoint=run.acceptedModel(model);const auto firstState=RStructuralPiezoLayeredRestartCodec::encode(run.acceptedRuntime());
    while(run.hasRemainingSteps()){auto next=run.prepare(model);run.validatePrepared(*next);run.commit(*next);QVERIFY(std::abs(run.acceptedState().energyBalanceResidual)<1e-8);QVERIFY(run.acceptedState().freeForceResidual<1e-8);}
    const auto finalState=RStructuralPiezoLayeredRestartCodec::encode(run.acceptedRuntime());QTemporaryDir dir;QVERIFY(dir.isValid());for(const auto&ext:{QString("ram"),QString("rbm")}){
        auto path=dir.filePath("solid."+ext);const auto actual=checkpoint.nativeWritePaths(path,false).first();checkpoint.write(path,false);RModel disk;disk.read(actual);QCOMPARE(disk.getProblemSetup().getStructuralPiezoLayeredSetup().family,RStructuralPiezoLayeredModelFamily::Solid3D);disk.getProblemSetup().setRestart(true);RStructuralPiezoLayeredNativeRun restored(disk);QCOMPARE(RStructuralPiezoLayeredRestartCodec::encode(restored.acceptedRuntime()),firstState);
        while(restored.hasRemainingSteps()){auto next=restored.prepare(disk);restored.validatePrepared(*next);restored.commit(*next);}QCOMPARE(RStructuralPiezoLayeredRestartCodec::encode(restored.acceptedRuntime()),finalState);
    }
}CATCH_FAILURE
void TestStructuralPiezoSolid::rotatedMaterialFrameAndPolarity()try{
    auto g=grid(false,1,1,2);const auto plain=assemble(g);auto reverse=g;for(auto&m:reverse.members)m.layers[0].polarity=-1;auto reversed=assemble(reverse);QVERIFY(close(plain.operators.capacitance[0],reversed.operators.capacitance[0]));for(size_t i=0;i<plain.operators.dofs;++i)QVERIFY(close(plain.operators.coupling[i],-reversed.operators.coupling[i]));
    const double c=std::sqrt(.5);RStructuralPiezoPlanarFrame frame;frame.ex={{c,c,0}};frame.ey={{-c,c,0}};frame.ez={{0,0,1}};frame.origin={{11,-7,5}};
    for(auto&n:g.nodes){auto x=n.position;for(size_t i=0;i<3;++i)n.position[i]=frame.origin[i]+frame.ex[i]*x[0]+frame.ey[i]*x[1]+frame.ez[i]*x[2];}
    const auto rotated=RStructuralPiezoLayeredMesh::assemble(RStructuralPiezoLayeredFamily::Solid,frame,g.nodes,g.members);QVERIFY(close(rotated.operators.capacitance[0],plain.operators.capacitance[0]));QVERIFY(close(rotated.operators.totalMass,7));
    V q(plain.operators.dofs),r(q.size());for(size_t i=0;i<g.nodes.size();++i){q[3*i]=double(i)/100;q[3*i+1]=double(i*i)/1000;q[3*i+2]=-.02;for(size_t k=0;k<3;++k)r[3*i+k]=frame.ex[k]*q[3*i]+frame.ey[k]*q[3*i+1]+frame.ez[k]*q[3*i+2];}QVERIFY(close(quadratic(plain.operators.stiffness,q),quadratic(rotated.operators.stiffness,r)));
}CATCH_FAILURE
QTEST_APPLESS_MAIN(TestStructuralPiezoSolid)
#include "tst_structural_piezo_solid.moc"

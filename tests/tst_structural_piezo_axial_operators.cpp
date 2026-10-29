#include <QtTest>
#include <rbl_error.h>
#include <rml_model.h>
#include "rstructural_piezo_axial_operators.h"
#include "rstructural_piezo_element_operators.h"
#include <algorithm>
#include <cmath>

namespace {
using Operators=RStructuralPiezoAxialOperators;
using Mode=RStructuralPiezoReadoutMode;
using Role=RStructuralPiezoPortRole;
bool near(double a,double b,double tol=2e-9)
{ return std::fabs(a-b)<=tol*std::max({1.0,std::fabs(a),std::fabs(b)}); }
template<class F> bool rejects(F f)
{ try { f(); } catch(const RError &) { return true; } catch(const std::exception &) { return true; } return false; }
template<class F> QString rejectionMessage(F f)
{ try { f(); } catch(const RError &error) { return error.getMessage(); } catch(const std::exception &error) { return QString::fromUtf8(error.what()); } return {}; }
void component(RBoundaryCondition &bc,RVariableType type,double value,bool enabled=true)
{
    auto &c=bc.getComponent(bc.findComponentPosition(type)); c.setEnabled(enabled); c.clear(); if(enabled)c.add(0,value);
}
RBoundaryCondition support(uint axis,bool fixed,bool transverse=true)
{
    RBoundaryCondition bc(R_BOUNDARY_CONDITION_DISPLACEMENT);
    const RVariableType types[]={R_VARIABLE_DISPLACEMENT_X,R_VARIABLE_DISPLACEMENT_Y,R_VARIABLE_DISPLACEMENT_Z};
    for(uint a=0;a<3;++a) component(bc,types[a],0,fixed || (transverse && a!=axis));
    return bc;
}
RStructuralPiezoAxialMemberProperties material()
{
    auto m=RStructuralPiezoContinuumMaterial::isotropic(100,.25,7);
    m.piezoelectric=true; m.coupling[0]=2;
    m.permittivity[0]=m.permittivity[4]=m.permittivity[8]=.5;
    return rStructuralPiezoCompileAxialMember(m,1,1);
}
RModel fixture(uint axis=0,double q0=.04,bool transverse=true,bool explicitReference=true,
               double tension=1,double voltage=.7)
{
    RModel m; m.addNode(RNode(0,0,0)); m.addNode(axis==0?RNode(1,0,0):RNode(0,1,0));
    m.setNodePhysicalId(0,909); m.setNodePhysicalId(1,101);
    for(uint i=0;i<6;++i) {
        RElement e(R_ELEMENT_TRUSS1); e.setNodeId(0,0); e.setNodeId(1,1);
        m.addElement(e,true,i); m.setElementPhysicalId(i,901+i);
        RMaterial mat(RMaterial::Solid);
        RMaterialProperty young(RMaterialProperty::ModulusOfElasiticity); young.add(293.15,10); mat.add(young);
        RMaterialProperty density(RMaterialProperty::Density); density.add(293.15,1); mat.add(density);
        auto &line=m.getLine(i); line.setCrossArea(1); line.setMaterial(mat);
        if(transverse) line.addBoundaryCondition(support(axis,false));
    }
    for(uint node=0;node<2;++node) {
        RElement e(R_ELEMENT_POINT); e.setNodeId(0,node); m.addElement(e,true,node);
    }
    m.getPoint(0).addBoundaryCondition(support(axis,true));
    RInitialCondition initial(R_INITIAL_CONDITION_DISPLACEMENT);
    const RVariableType qs[]={R_VARIABLE_DISPLACEMENT_X,R_VARIABLE_DISPLACEMENT_Y,R_VARIABLE_DISPLACEMENT_Z};
    for(uint a=0;a<3;++a) {
        auto &c=initial.getComponent(initial.findComponentPosition(qs[a])); c.setEnabled(true); c.clear();
        c.add(0,a==axis?(explicitReference?.02:0):0);
    }
    m.getPoint(1).addInitialCondition(initial);
    const double charge=explicitReference?.03:0;
    const double force=140*q0+2*voltage-.4*(charge-.4*q0)/.5+std::max(0.0,tension+10*q0);
    RBoundaryCondition load(R_BOUNDARY_CONDITION_FORCE);
    const RVariableType fs[]={R_VARIABLE_FORCE_X,R_VARIABLE_FORCE_Y,R_VARIABLE_FORCE_Z};
    for(uint a=0;a<3;++a) component(load,fs[a],a==axis?force:0,a==axis);
    m.getPoint(1).addBoundaryCondition(load);
    for(uint i=0;i<m.getNElementGroups();++i) {
        auto *group=m.getElementGroupPtr(i);
        if(!group || !group->hasEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION)) continue;
        auto &g=group->getEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
        for(auto type:{R_VARIABLE_G_ACCELERATION_X,R_VARIABLE_G_ACCELERATION_Y,R_VARIABLE_G_ACCELERATION_Z}) {
            auto &c=g.getComponent(g.findComponentPosition(type)); c.setEnabled(true); c.clear(); c.add(0,0);
        }
    }
    m.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS)); m.getProblemSetup().setMassMatrixType(RMassMatrixType::Lumped);
    auto &d=m.getProblemSetup().getStructuralDynamicsSetup(); d.setEnabled(true); d.setTimeStep(.004); d.setNSteps(2); d.setOutputFrequency(1);
    auto &time=m.getTimeSolver(); time.setEnabled(true); time.setInputStartTime(.004); time.setInputTimeStepSize(.004);
    time.setInputNTimeSteps(2); time.setTimes({.004,.008});
    auto &truss=m.getProblemSetup().getStructuralTrussDynamicsSetup(); truss.enabled=true;
    truss.absoluteForceTolerance=1e-11; truss.relativeForceTolerance=1e-11; truss.linearRelativeTolerance=1e-12;
    truss.passiveMembers={{901,1}};
    const double coupling[]={2,-.4,.5,-.6}; const QString ids[]={"actuator","open","short","amplifier"};
    const auto property=material();
    for(uint i=0;i<4;++i) {
        RStructuralTrussActiveMember a; a.physicalElementId=902+i; a.portId=ids[i]; a.referenceLength=1;
        a.axialStiffness=i==0?property.k:10; a.coupling=coupling[i]; a.capacitance=.5;
        a.kind=i==1?RStructuralTrussElectricalKind::ConductorCharge:RStructuralTrussElectricalKind::Voltage;
        a.initialValue=i==0?voltage:(i==1?charge:0);
        if(i==0) { a.massSource=RStructuralTrussMassSource::FixedMemberMass; a.fixedMass=property.fixedMass; }
        truss.activeMembers.push_back(a);
    }
    RStructuralTrussCableMember cable; cable.physicalElementId=906; cable.initialReferenceLength=1; cable.initialTension=tension;
    truss.cableMembers.push_back(cable);
    auto &setup=m.getProblemSetup().getStructuralPiezoReadoutSetup();
    setup.ports={{"open",Role::Sensor,Mode::OpenCircuitVoltage},
                 {"amplifier",Role::Sensor,Mode::IdealChargeAmplifier,.25,-1,.012,1.3},
                 {"actuator",Role::ActuatorSensor,Mode::ActuatorCharge},
                 {"short",Role::Sensor,Mode::ShortCircuitCharge}};
    if(explicitReference) {
        setup.workpoint.available=true;
        RStructuralPiezoNodeWorkpoint node; node.physicalNodeId=101; node.q0[axis]=q0; node.f0[axis]=force;
        setup.workpoint.nodes={node,{909,{{0,0,0}},{{0,0,0}}}};
        setup.workpoint.ports={{"amplifier",0,.6*q0},{"short",0,-.5*q0},
                              {"open",(charge-.4*q0)/.5,charge},{"actuator",voltage,.5*voltage-2*q0}};
    }
    truss.bindToModel(m); return m;
}
std::size_t outputRow(const Operators &o,const QString &id)
{ for(std::size_t i=0;i<o.outputs.size();++i) if(o.outputs[i].portId==id) return i; throw std::invalid_argument("Missing output"); }
}

#define REPORT_AXIAL_EXCEPTIONS \
    catch(const RError &error) { QFAIL(qPrintable(error.getMessage())); } \
    catch(const std::exception &error) { QFAIL(error.what()); }

class TestStructuralPiezoAxialOperators:public QObject
{
    Q_OBJECT
private slots:
    void actualMixedMassTangentAndPreload();
    void nonparallelMultiNodeNetwork();
    void independentElectricalFiniteDifferences();
    void readoutBiasAndDirectTerms();
    void rotationAndConstraintMapping();
    void rejectsUnsupportedChargeActuator();
    void implicitZeroReferenceCertificate();
    void scaledCondensationAndExtremeScaleGates();
    void rejectsInvalidReferenceAndSourceAtomically();
};

void TestStructuralPiezoAxialOperators::actualMixedMassTangentAndPreload() try
{
    const auto m=fixture(); const auto o=Operators::compile(m);
    QCOMPARE(o.nativeSystem.members.size(),std::size_t(6)); QCOMPARE(o.free.dofs,std::size_t(1));
    QCOMPARE(o.freeToFull.front(),std::size_t(3)); QCOMPARE(o.freeDofs.front().physicalNodeId,uint(101));
    QCOMPARE(o.fullDofs.front().physicalNodeId,uint(909)); QCOMPARE(o.orderedMemberPhysicalIds.front(),uint(901));
    QVERIFY(near(o.free.M.front(),6)); // compiled replacement mass 7 + five independent unit masses, split in half.
    for(std::size_t i=0;i<6;++i) QCOMPARE(o.full.M[i*6+i],o.nativeSystem.nodalLumpedMass[i/3]);
    QVERIFY(near(o.free.Kmaterial.front(),150)); QVERIFY(near(o.free.KfixedCharge.front(),.32));
    QVERIFY(near(o.free.K.front(),150.32)); QVERIFY(near(o.free.Kshort.front(),150));
    QVERIFY(near(o.full.Kgeometric[4*6+4],o.f0[3]/1.04));
    QVERIFY(o.full.Kgeometric[4*6+4]>0); QCOMPARE(o.free.C.front(),0.0);
    QVERIFY(o.reference.equilibriumCertified); QVERIFY(o.reference.explicitWorkpoint);
    QVERIFY(o.reference.freeForceResidualNorm<=o.reference.forceTolerance);
    QVERIFY(near(o.freeInitialDeltaQ.front(),-.02)); QVERIFY(near(o.q0[3],.04));
    QCOMPARE(o.inputs.size(),std::size_t(2)); QCOMPARE(o.inputs[0].unit,QString("V"));
    QCOMPARE(o.inputs[1].physicalNodeId,uint(101)); QCOMPARE(o.inputs[1].axis,uint(0));
    QVERIFY(near(o.free.B[0],-2)); QCOMPARE(o.free.B[1],1.0);
    // Fixed-charge minus short-circuit tangent is PSD, including full incidence cross blocks.
    for(const std::vector<double> v:std::vector<std::vector<double>>{{1,2,3,4,5,6},{-2,0,1,3,0,-1},{1,1,1,1,1,1}}) {
        double energy=0;
        for(std::size_t i=0;i<6;++i) for(std::size_t j=0;j<6;++j) energy+=v[i]*(o.full.K[i*6+j]-o.full.Kshort[i*6+j])*v[j];
        QVERIFY(energy>=-1e-12);
    }
} REPORT_AXIAL_EXCEPTIONS

void TestStructuralPiezoAxialOperators::nonparallelMultiNodeNetwork() try
{
    auto m=fixture(); m.addNode(RNode(1,1,0)); m.setNodePhysicalId(2,707);
    RElement e(R_ELEMENT_TRUSS1); e.setNodeId(0,1); e.setNodeId(1,2);
    // The fixture's POINT elements already own physical IDs 907 and 908.
    const uint memberPhysicalId=910;
    const auto index=m.getNElements(); m.addElement(e,true,6); m.setElementPhysicalId(index,memberPhysicalId);
    m.getLine(6).setCrossArea(1); m.getLine(6).setMaterial(m.getLine(0).getMaterial());
    RElement p(R_ELEMENT_POINT); p.setNodeId(0,2); m.addElement(p,true,2);
    m.getPoint(2).addBoundaryCondition(support(0,true));
    const double length=std::hypot(.04,1.0);
    auto &truss=m.getProblemSetup().getStructuralTrussDynamicsSetup(); truss.passiveMembers.push_back({memberPhysicalId,length});
    m.getProblemSetup().getStructuralPiezoReadoutSetup().workpoint.nodes.push_back({707,{{0,0,0}},{{0,0,0}}});
    for(uint i=0;i<m.getNElementGroups();++i) {
        auto *group=m.getElementGroupPtr(i);
        if(!group || !group->hasEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION)) continue;
        auto &g=group->getEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
        for(auto type:{R_VARIABLE_G_ACCELERATION_X,R_VARIABLE_G_ACCELERATION_Y,R_VARIABLE_G_ACCELERATION_Z}) {
            auto &c=g.getComponent(g.findComponentPosition(type)); c.setEnabled(true); c.clear(); c.add(0,0);
        }
    }
    truss.bindToModel(m); const auto o=Operators::compile(m);
    QCOMPARE(o.nativeSystem.members.size(),std::size_t(7)); QCOMPARE(o.full.dofs,std::size_t(9));
    QCOMPARE(o.fullDofs[6].physicalNodeId,uint(707)); QVERIFY(near(o.free.M.front(),6.5));
    QVERIFY(near(o.free.K.front(),150.32+10*.04*.04/(length*length*length)));
    QVERIFY(near(o.full.Kmaterial[3*9+4],-10*.04/(length*length*length)));
    QCOMPARE(o.freeToFull,std::vector<std::size_t>{3});
} REPORT_AXIAL_EXCEPTIONS

void TestStructuralPiezoAxialOperators::independentElectricalFiniteDifferences() try
{
    const auto o=Operators::compile(fixture()); const auto &s=o.nativeSystem;
    const double h=1e-6;
    for(std::size_t j=0;j<o.orderedPortIds.size();++j) {
        const auto member=std::find_if(s.members.begin(),s.members.end(),[&](const auto &m){
            return m.physicalElementId==uint(902+j);});
        QVERIFY(member!=s.members.end()); const auto index=std::size_t(member-s.members.begin());
        auto plus=o.reference.inputs.memberInputs,minus=plus;
        plus[index]=RStructuralTrussVoltage{o.V0[j]+h}; minus[index]=RStructuralTrussVoltage{o.V0[j]-h};
        const auto fp=RStructuralTrussAssembly::evaluate(s.referenceNodes,o.q0,s.members,plus);
        const auto fm=RStructuralTrussAssembly::evaluate(s.referenceNodes,o.q0,s.members,minus);
        for(std::size_t d=0;d<6;++d) {
            QVERIFY(near((fp.internalForces[d]-fm.internalForces[d])/(2*h),o.full.G[d*4+j],2e-8));
            auto qp=o.q0,qm=o.q0; qp[d]+=h; qm[d]-=h;
            plus[index]=RStructuralTrussVoltage{o.V0[j]};
            const auto ep=RStructuralTrussAssembly::evaluate(s.referenceNodes,qp,s.members,plus);
            const auto em=RStructuralTrussAssembly::evaluate(s.referenceNodes,qm,s.members,plus);
            QVERIFY(near((ep.memberResponses[index].electrical->conductorCharge-em.memberResponses[index].electrical->conductorCharge)/(2*h),
                         -o.full.G[d*4+j],2e-8));
        }
    }
    // Updated geometry tangent also agrees with native force differences under the actual mixed boundaries.
    for(std::size_t j=0;j<6;++j) {
        auto qp=o.q0,qm=o.q0; qp[j]+=h; qm[j]-=h;
        const auto fp=RStructuralTrussAssembly::evaluate(s.referenceNodes,qp,s.members,o.reference.inputs.memberInputs);
        const auto fm=RStructuralTrussAssembly::evaluate(s.referenceNodes,qm,s.members,o.reference.inputs.memberInputs);
        for(std::size_t i=0;i<6;++i) QVERIFY(near((fp.internalForces[i]-fm.internalForces[i])/(2*h),o.full.K[i*6+j],2e-8));
    }
} REPORT_AXIAL_EXCEPTIONS

void TestStructuralPiezoAxialOperators::readoutBiasAndDirectTerms() try
{
    const auto model=fixture(); const auto o=Operators::compile(model); RStructuralTrussDynamicsModelAdapter a(model);
    const auto open=outputRow(o,"open"),amp=outputRow(o,"amplifier"),act=outputRow(o,"actuator"),shorted=outputRow(o,"short");
    QVERIFY(near(o.free.Lq[open],-.8)); QVERIFY(near(o.free.Lq[amp],-.6/.25));
    QVERIFY(near(o.free.Lq[act],-2)); QVERIFY(near(o.free.Lq[shorted],-.5));
    QVERIFY(near(o.free.D[act*2],.5)); QCOMPARE(o.free.D[open*2],0.0);
    QVERIFY(near(o.outputs[amp].value0,1.3-(.6*.04-.012)/.25));
    auto q=o.q0; q[3]+=.001; auto inputs=o.reference.inputs; inputs.memberInputs[1]=RStructuralTrussVoltage{.703};
    const auto state=a.driver().initialize(0,q,std::vector<double>(6,0),inputs);
    const auto samples=a.piezoReadout().observe(state);
    for(std::size_t row=0;row<samples.size();++row) {
        const double predicted=o.free.Lq[row]*.001+o.free.D[row*2]*.003;
        QVERIFY(near(predicted,samples[row].deltaValue)); QVERIFY(near(o.outputs[row].value0+predicted,samples[row].value));
    }
    for(double value:o.full.Lv) QCOMPARE(value,0.0);
    for(double value:o.full.La) QCOMPARE(value,0.0);
    auto echo=model; auto &setup=echo.getProblemSetup().getStructuralPiezoReadoutSetup();
    setup.ports[2].role=Role::Actuator; setup.ports[2].readout=Mode::None;
    echo.getProblemSetup().getStructuralTrussDynamicsSetup().bindToModel(echo);
    const auto e=Operators::compile(echo); const auto row=outputRow(e,"actuator");
    QVERIFY(e.outputs[row].actionRecord); QCOMPARE(e.free.Lq[row],0.0); QCOMPARE(e.free.D[row*2],1.0);
} REPORT_AXIAL_EXCEPTIONS

void TestStructuralPiezoAxialOperators::rotationAndConstraintMapping() try
{
    const auto x=Operators::compile(fixture(0)),y=Operators::compile(fixture(1));
    QCOMPARE(y.freeToFull.front(),std::size_t(4)); QCOMPARE(y.freeDofs.front().axis,uint(1));
    QCOMPARE(y.inputs.back().axis,uint(1)); QCOMPARE(x.free.M,y.free.M);
    QCOMPARE(x.free.K,y.free.K); QCOMPARE(x.free.G,y.free.G); QCOMPARE(x.free.B,y.free.B); QCOMPARE(x.free.Lq,y.free.Lq);
    // With sufficient actual tensile preload the transverse DOFs are stable and remain in the free map.
    const auto all=Operators::compile(fixture(0,.04,false)); QCOMPARE(all.free.dofs,std::size_t(3));
    QVERIFY(near(all.free.K[4],all.f0[3]/1.04)); QVERIFY(all.minimumNormalizedStiffnessPivot>0);
    auto reversed=fixture(); for(uint i=0;i<6;++i) { reversed.getElement(i).setNodeId(0,1); reversed.getElement(i).setNodeId(1,0); }
    reversed.getProblemSetup().getStructuralTrussDynamicsSetup().bindToModel(reversed);
    const auto r=Operators::compile(reversed); QCOMPARE(x.free.G,r.free.G); QCOMPARE(x.free.K,r.free.K);
    // An arbitrary 3D rotation retains actual off-diagonal physical operators.
    auto diagonal=fixture(0,.04,false); const double direction=1/std::sqrt(3.0);
    diagonal.getNode(1).setX(direction); diagonal.getNode(1).setY(direction); diagonal.getNode(1).setZ(direction);
    const RVariableType qs[]={R_VARIABLE_DISPLACEMENT_X,R_VARIABLE_DISPLACEMENT_Y,R_VARIABLE_DISPLACEMENT_Z};
    auto &initial=diagonal.getPoint(1).getInitialCondition(R_INITIAL_CONDITION_DISPLACEMENT);
    auto &load=diagonal.getPoint(1).getBoundaryCondition(R_BOUNDARY_CONDITION_FORCE);
    const RVariableType fs[]={R_VARIABLE_FORCE_X,R_VARIABLE_FORCE_Y,R_VARIABLE_FORCE_Z};
    auto &node=diagonal.getProblemSetup().getStructuralPiezoReadoutSetup().workpoint.nodes[0];
    for(uint axis=0;axis<3;++axis) {
        auto &c=initial.getComponent(initial.findComponentPosition(qs[axis])); c.clear(); c.add(0,.02*direction);
        node.q0[axis]=.04*direction; node.f0[axis]=all.f0[3]*direction;
        component(load,fs[axis],node.f0[axis]);
    }
    diagonal.getProblemSetup().getStructuralTrussDynamicsSetup().bindToModel(diagonal);
    const auto rotated=Operators::compile(diagonal); QCOMPARE(rotated.free.dofs,std::size_t(3));
    for(std::size_t i=0;i<3;++i) {
        QVERIFY(near(rotated.free.G[i*4],2*direction));
        for(std::size_t j=0;j<3;++j) {
            const double nn=direction*direction;
            QVERIFY(near(rotated.free.K[i*3+j],150.32*nn+(all.f0[3]/1.04)*((i==j?1:0)-nn)));
        }
    }
} REPORT_AXIAL_EXCEPTIONS

void TestStructuralPiezoAxialOperators::rejectsUnsupportedChargeActuator() try
{
    auto m=fixture(); auto &t=m.getProblemSetup().getStructuralTrussDynamicsSetup();
    auto &p=t.activeMembers.front(); p.kind=RStructuralTrussElectricalKind::ConductorCharge; p.initialValue=.5*.7-2*.04;
    QVERIFY(rejects([&]{t.bindToModel(m); Operators::compile(m);}));
} REPORT_AXIAL_EXCEPTIONS

void TestStructuralPiezoAxialOperators::implicitZeroReferenceCertificate() try
{
    const auto o=Operators::compile(fixture(0,0,true,false,1,0));
    QVERIFY(!o.reference.explicitWorkpoint); QVERIFY(o.reference.equilibriumCertified);
    QCOMPARE(o.certificateAuthority,QString("native-source-reference-with-operator-compiler-equilibrium-check"));
    QCOMPARE(o.reference.freeForceResidualNorm,0.0); QCOMPARE(o.freeInitialDeltaQ.front(),0.0);
    auto slack=fixture(0,0,true,false,0,0);
    auto &t=slack.getProblemSetup().getStructuralTrussDynamicsSetup(); t.cableMembers.front().initialReferenceLength=1.1;
    t.bindToModel(slack); const auto compiled=Operators::compile(slack);
    QVERIFY(near(compiled.free.K.front(),140.32));
    QVERIFY(!*compiled.reference.assembly.memberResponses.back().cableTaut);
} REPORT_AXIAL_EXCEPTIONS

void TestStructuralPiezoAxialOperators::rejectsInvalidReferenceAndSourceAtomically() try
{
    const auto pristine=fixture(); auto model=pristine; const auto before=Operators::compile(pristine);
    model.getNode(1).setX(1.01); QVERIFY(rejects([&]{Operators::compile(model);})); // stale source binding
    QCOMPARE(before.q0[3],.04); QVERIFY(before.reference.equilibriumCertified);
    auto invalid=fixture(); invalid.getProblemSetup().getStructuralPiezoReadoutSetup().workpoint.ports.back().Q0+=1e-4;
    invalid.getProblemSetup().getStructuralTrussDynamicsSetup().bindToModel(invalid);
    QVERIFY(rejects([&]{Operators::compile(invalid);}));
    QVERIFY(rejects([&]{Operators::compile(fixture(0,0,true,true,0));})); // exact cable kink
    QVERIFY(rejects([&]{Operators::compile(fixture(0,-.04,false));})); // compressive transverse instability
    auto mechanism=fixture(0,0,false,false,0,0);
    auto &t=mechanism.getProblemSetup().getStructuralTrussDynamicsSetup(); t.cableMembers.front().initialReferenceLength=1.1;
    t.bindToModel(mechanism); QVERIFY(rejects([&]{Operators::compile(mechanism);})); // slack cable away from kink, free transverse mechanism
    auto nonzero=fixture(); nonzero.getProblemSetup().getStructuralPiezoReadoutSetup().workpoint={};
    nonzero.getProblemSetup().getStructuralTrussDynamicsSetup().bindToModel(nonzero);
    QVERIFY(rejects([&]{Operators::compile(nonzero);}));
    auto moment=fixture(); RBoundaryCondition torque(R_BOUNDARY_CONDITION_MOMENT);
    component(torque,R_VARIABLE_MOMENT_Z,1); moment.getPoint(1).addBoundaryCondition(torque);
    QVERIFY(rejects([&]{moment.getProblemSetup().getStructuralTrussDynamicsSetup().bindToModel(moment); Operators::compile(moment);}));
    auto damped=fixture(); damped.getProblemSetup().getStructuralDynamicsSetup().setRayleighAlpha(.01);
    QVERIFY(rejects([&]{damped.getProblemSetup().getStructuralTrussDynamicsSetup().bindToModel(damped); Operators::compile(damped);}));
    const auto after=Operators::compile(pristine); QCOMPARE(before.sourceFingerprint,after.sourceFingerprint);
    QCOMPARE(before.full.K,after.full.K); QCOMPARE(before.nativeSystemIdentity,after.nativeSystemIdentity);
} REPORT_AXIAL_EXCEPTIONS

void TestStructuralPiezoAxialOperators::scaledCondensationAndExtremeScaleGates() try
{
    auto model=fixture(0,0,true,false,1,0);
    model.getNode(1).setY(1e-100);
    auto &truss=model.getProblemSetup().getStructuralTrussDynamicsSetup();
    truss.activeMembers[1].coupling=1e100; truss.activeMembers[1].capacitance=1e300;
    truss.bindToModel(model); const auto o=Operators::compile(model);
    // These expectations deliberately use a relative ratio, with no unit-sized
    // absolute tolerance that could accept erased tiny off-diagonal entries.
    const double xy=o.full.KfixedCharge[3*6+4],yx=o.full.KfixedCharge[4*6+3];
    QVERIFY(xy>0); QCOMPARE(xy,yx); QVERIFY(std::fabs(xy/1e-200-1)<2e-14);
    QVERIFY(std::fabs(o.full.KfixedCharge[4*6+4]/1e-300-1)<2e-14);
    const auto row=outputRow(o,"open");
    QVERIFY(std::fabs(o.full.Lq[row*6+4]/1e-300-1)<2e-14);
    for(std::size_t i=0;i<6;++i) for(std::size_t j=0;j<6;++j)
        QCOMPARE(o.full.KfixedCharge[i*6+j],o.full.KfixedCharge[j*6+i]);

    // Real source chord (1,1e-200,0), g=1e100, Ce=1e300: xy would
    // be 1e-300, but yy contains nonzero terms below binary64's range.
    // The compiler explicitly rejects this unsupported scale rather than
    // silently returning asymmetric condensation or claiming completeness.
    auto tiny=model; tiny.getNode(1).setY(1e-200);
    tiny.getProblemSetup().getStructuralTrussDynamicsSetup().bindToModel(tiny);
    const auto tinyMessage=rejectionMessage([&]{Operators::compile(tiny);});
    QVERIFY2(tinyMessage.contains("outside binary64 representability"),qPrintable(tinyMessage));

    // The source assembler forms n_y*n_y first. With this huge material
    // stiffness it erases a representable yy contribution (1e-200). A large
    // unrelated xx entry (1e200) must not excuse the local discrepancy.
    auto unresolved=fixture(0,0,true,false,0,0); unresolved.getNode(1).setY(1e-200);
    auto &setup=unresolved.getProblemSetup().getStructuralTrussDynamicsSetup();
    setup.cableMembers.front().initialReferenceLength=1.1;
    for(auto &member:setup.activeMembers) { member.axialStiffness=1e200; member.coupling=0; }
    auto &mat=unresolved.getLine(0).getMaterial();
    mat.get(mat.findPosition(RMaterialProperty::ModulusOfElasiticity)).add(293.15,1e200);
    setup.bindToModel(unresolved);
    const auto localMessage=rejectionMessage([&]{Operators::compile(unresolved);});
    QVERIFY2(localMessage.contains("cannot resolve the material/geometric entry at its local scale"),qPrintable(localMessage));
    QCOMPARE(Operators::compile(model).full.KfixedCharge,o.full.KfixedCharge);
} REPORT_AXIAL_EXCEPTIONS

QTEST_APPLESS_MAIN(TestStructuralPiezoAxialOperators)
#include "tst_structural_piezo_axial_operators.moc"

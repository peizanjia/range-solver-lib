#include <QtTest>
#include <rbl_error.h>
#include <rml_model.h>
#include "rstructural_piezo_native_readout.h"
#include "rstructural_truss_dynamics_restart_adapter.h"
#include "rstructural_piezo_element_operators.h"
#include <cmath>
#include <limits>

namespace {
using Adapter=RStructuralTrussDynamicsModelAdapter;
using Readout=RStructuralPiezoNativeReadout;
using Mode=RStructuralPiezoReadoutMode;
using Role=RStructuralPiezoPortRole;
template<class F> bool rejects(F f)
{ try { f(); } catch(const RError &) { return true; } catch(const std::exception &) { return true; } return false; }
bool near(double a,double b,double relative=2e-11)
{ return a==b || (b!=0.0 && std::fabs((a-b)/b)<=relative); }
void component(RBoundaryCondition &bc,RVariableType type,double value,bool enabled=true)
{
    auto &c=bc.getComponent(bc.findComponentPosition(type)); c.setEnabled(enabled); c.clear(); if(enabled)c.add(0,value);
}
RBoundaryCondition support(uint axis,bool fixed)
{
    RBoundaryCondition bc(R_BOUNDARY_CONDITION_DISPLACEMENT);
    const RVariableType types[]={R_VARIABLE_DISPLACEMENT_X,R_VARIABLE_DISPLACEMENT_Y,R_VARIABLE_DISPLACEMENT_Z};
    for(uint a=0;a<3;++a) component(bc,types[a],0.0,fixed || a!=axis);
    return bc;
}
void point(RModel &m,uint node,uint group)
{ RElement e(R_ELEMENT_POINT); e.setNodeId(0,node); m.addElement(e,true,group); }
RModel fixture(uint axis=0,double force=.7,double displacement=.02,double voltage=.7)
{
    RModel m; m.addNode(RNode(0,0,0)); m.addNode(axis==0?RNode(1,0,0):RNode(0,1,0));
    m.setNodePhysicalId(0,909); m.setNodePhysicalId(1,101);
    for(uint i=0;i<5;++i) {
        RElement e(R_ELEMENT_TRUSS1); e.setNodeId(0,0); e.setNodeId(1,1);
        m.addElement(e,true,i); m.setElementPhysicalId(i,901+i);
        RMaterial material(RMaterial::Solid);
        RMaterialProperty young(RMaterialProperty::ModulusOfElasiticity); young.add(293.15,10); material.add(young);
        RMaterialProperty density(RMaterialProperty::Density); density.add(293.15,1); material.add(density);
        auto &line=m.getLine(i); line.setCrossArea(1); line.setMaterial(material); line.addBoundaryCondition(support(axis,false));
    }
    point(m,0,0); m.getPoint(0).addBoundaryCondition(support(axis,true)); point(m,1,1);
    RInitialCondition q(R_INITIAL_CONDITION_DISPLACEMENT);
    const RVariableType qs[]={R_VARIABLE_DISPLACEMENT_X,R_VARIABLE_DISPLACEMENT_Y,R_VARIABLE_DISPLACEMENT_Z};
    for(uint a=0;a<3;++a) { auto &c=q.getComponent(q.findComponentPosition(qs[a])); c.setEnabled(true); c.clear(); c.add(0,a==axis?displacement:0); }
    m.getPoint(1).addInitialCondition(q);
    RBoundaryCondition load(R_BOUNDARY_CONDITION_FORCE);
    const RVariableType fs[]={R_VARIABLE_FORCE_X,R_VARIABLE_FORCE_Y,R_VARIABLE_FORCE_Z};
    for(uint a=0;a<3;++a) component(load,fs[a],a==axis?force:0.0,a==axis);
    m.getPoint(1).addBoundaryCondition(load);
    for(uint i=0;i<m.getNElementGroups();++i) {
        auto *group=m.getElementGroupPtr(i);
        if(!group || !group->hasEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION)) continue;
        auto &gravity=group->getEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
        for(auto type:{R_VARIABLE_G_ACCELERATION_X,R_VARIABLE_G_ACCELERATION_Y,R_VARIABLE_G_ACCELERATION_Z}) {
            auto &c=gravity.getComponent(gravity.findComponentPosition(type)); c.setEnabled(true); c.clear(); c.add(0,0);
        }
    }
    m.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS)); m.getProblemSetup().setMassMatrixType(RMassMatrixType::Lumped);
    auto &d=m.getProblemSetup().getStructuralDynamicsSetup(); d.setEnabled(true); d.setTimeStep(.004); d.setNSteps(24); d.setOutputFrequency(1);
    auto &time=m.getTimeSolver(); time.setEnabled(true); time.setInputStartTime(.004); time.setInputTimeStepSize(.004); time.setInputNTimeSteps(24);
    std::vector<double> times; for(uint i=1;i<=24;++i) times.push_back(i*.004); time.setTimes(times);
    auto &truss=m.getProblemSetup().getStructuralTrussDynamicsSetup(); truss.enabled=true;
    truss.absoluteForceTolerance=1e-11; truss.relativeForceTolerance=1e-11; truss.linearRelativeTolerance=1e-12;
    truss.passiveMembers={{901,1}};
    const double coupling[]={.3,-.4,.5,-.6};
    const QString ports[]={"actuator","open","short","amplifier"};
    for(uint i=0;i<4;++i) {
        RStructuralTrussActiveMember a; a.physicalElementId=902+i; a.portId=ports[i];
        a.referenceLength=1; a.axialStiffness=10; a.coupling=coupling[i]; a.capacitance=.5;
        a.kind=i==1?RStructuralTrussElectricalKind::ConductorCharge:RStructuralTrussElectricalKind::Voltage;
        a.initialValue=i==0?voltage:(i==1?.03:0.0); truss.activeMembers.push_back(a);
    }
    truss.bindToModel(m); return m;
}
RStructuralPiezoReadoutSetup catalog()
{
    RStructuralPiezoReadoutSetup s;
    // Catalog order deliberately differs from physical member order.
    s.ports={{"open",Role::Sensor,Mode::OpenCircuitVoltage},
             {"amplifier",Role::Sensor,Mode::IdealChargeAmplifier,.25,-1,.012,1.3},
             {"actuator",Role::ActuatorSensor,Mode::ActuatorCharge},
             {"short",Role::Sensor,Mode::ShortCircuitCharge}};
    return s;
}
RStructuralPiezoReadoutSetup explicitReference(const Adapter &a,double q0=.04)
{
    auto s=catalog(); s.workpoint.available=true;
    s.workpoint.nodes={{101,{{q0,0,0}},{{a.initialCompleteInputs().externalForces[3],0,0}}},
                       {909,{{0,0,0}},{{0,0,0}}}};
    s.workpoint.ports={{"amplifier",0,.6*q0},{"short",0,-.5*q0},
                       {"open",(.03-.4*q0)/.5,.03},{"actuator",.7,.5*.7-.3*q0}};
    return s;
}
const RStructuralPiezoNativeReadoutSample &sample(const std::vector<RStructuralPiezoNativeReadoutSample> &s,const QString &id)
{ for(const auto &p:s) if(p.portId==id) return p; throw std::invalid_argument("Test port missing"); }
}

#define REPORT_PIEZO_EXCEPTIONS \
    catch(const RError &error) { QFAIL(qPrintable(error.getMessage())); } \
    catch(const std::exception &error) { QFAIL(error.what()); }

class TestStructuralPiezoNativeReadout:public QObject
{
    Q_OBJECT
private slots:
    void mixedNativeTrialAndDirectTerm();
    void compiledMaterialNativeRodOpenShortFrequency();
    void explicitPreloadCertificateAndRejection();
    void frozenSensorBoundariesAndAtomicFailure();
    void rotationAndIndependentElectrodeSign();
    void smallIncrementSurvivesLargeBias();
    void splitRestorePreservesReference();
};

void TestStructuralPiezoNativeReadout::mixedNativeTrialAndDirectTerm() try
{
    auto model=fixture(); Adapter a(model); Readout r(a,catalog());
    QVERIFY(!r.reference().equilibriumCertified); QVERIFY(!r.reference().explicitWorkpoint);
    auto state=a.initialize();
    const auto initial=r.observe(state); QCOMPARE(initial.size(),std::size_t(4));
    QVERIFY(near(sample(initial,"open").voltage,(.03-.4*.02)/.5));
    QCOMPARE(sample(initial,"open").charge,.03);
    QVERIFY(near(sample(initial,"short").value,-.5*.02));
    QVERIFY(near(sample(initial,"amplifier").value,1.3-(.6*.02-.012)/.25));
    QVERIFY(near(sample(initial,"actuator").charge,.5*.7-.3*.02));
    for(const auto &s:initial) { QCOMPARE(s.deltaVoltage,0.0); QCOMPARE(s.deltaCharge,0.0); QCOMPARE(s.deltaValue,0.0); }
    auto inputs=state.appliedInputs; const auto index=a.maps().memberByPortId.at("actuator");
    inputs.memberInputs[index]=RStructuralTrussVoltage{.9};
    const auto trial=a.driver().trialStep(state,inputs,.004,a.options());
    const auto read=r.observe(trial.candidate); const double dq=trial.candidate.displacements[3]-.02;
    QVERIFY(near(sample(read,"actuator").deltaCharge,.5*(.9-.7)-.3*dq));
    QVERIFY(near(sample(read,"open").deltaVoltage,-.4*dq/.5));
    QCOMPARE(sample(read,"open").deltaCharge,0.0);
    QVERIFY(near(sample(read,"short").deltaCharge,-.5*dq));
    QVERIFY(near(sample(read,"amplifier").deltaValue,-.6*dq/.25));
    auto echoCatalog=catalog(); echoCatalog.ports[2].role=Role::Actuator; echoCatalog.ports[2].readout=Mode::None;
    Readout echo(a,echoCatalog); const auto value=sample(echo.observe(trial.candidate),"actuator");
    QVERIFY(value.actionRecord); QCOMPARE(value.value,.9); QCOMPARE(value.unit,QString("V"));
} REPORT_PIEZO_EXCEPTIONS

void TestStructuralPiezoNativeReadout::compiledMaterialNativeRodOpenShortFrequency() try
{
    auto material=RStructuralPiezoContinuumMaterial::isotropic(100,.25,7);
    material.piezoelectric=true;
    material.coupling[0]=2; // e_x,xx; other stress-electric entries are zero.
    material.permittivity[0]=material.permittivity[4]=material.permittivity[8]=.5;
    const auto property=rStructuralPiezoCompileAxialMember(material,1,1);
    QVERIFY(near(property.k,100)); QVERIFY(near(property.g,2));
    QVERIFY(near(property.Ce,.5)); QVERIFY(near(property.fixedMass,7));
    double frequencies[2]{};
    for(uint open=0;open<2;++open) {
        auto model=fixture(0,0,0,0);
        auto &truss=model.getProblemSetup().getStructuralTrussDynamicsSetup();
        truss.passiveMembers={{901,1},{903,1},{904,1},{905,1}};
        truss.activeMembers.resize(1);
        auto &rod=truss.activeMembers.front();
        rod.axialStiffness=property.k; rod.coupling=property.g; rod.capacitance=property.Ce;
        rod.massSource=RStructuralTrussMassSource::FixedMemberMass; rod.fixedMass=property.fixedMass;
        rod.kind=open?RStructuralTrussElectricalKind::ConductorCharge:RStructuralTrussElectricalKind::Voltage;
        rod.initialValue=0;
        truss.bindToModel(model);
        RStructuralPiezoReadoutSetup setup;
        setup.ports={{"actuator",Role::Sensor,open?Mode::OpenCircuitVoltage:Mode::ShortCircuitCharge}};
        Adapter adapter(model); Readout readout(adapter,setup);
        // The replacement rod contributes PZT mass exactly once. Other four
        // mechanical rods retain their independent unit member masses.
        QVERIFY(near(adapter.system().nodalLumpedMass[1],.5*(7+4)));
        const double tangent=readout.reference().assembly.tangent.findValue(3,3);
        const double expected=140+(open?8:0);
        QVERIFY(near(tangent,expected));
        frequencies[open]=std::sqrt(tangent/adapter.system().nodalLumpedMass[1]);
        QVERIFY(near(frequencies[open],std::sqrt(expected/5.5)));
        const auto state=adapter.initialize();
        const auto trial=adapter.trial(state,0,0);
        QCOMPARE(readout.observe(trial.state).front().value,0.0);
    }
    QVERIFY(frequencies[1]>frequencies[0]);
} REPORT_PIEZO_EXCEPTIONS

void TestStructuralPiezoNativeReadout::explicitPreloadCertificateAndRejection() try
{
    const double q0=.04,force=50*q0+.3*.7-.4*(.03-.4*q0)/.5;
    auto model=fixture(0,force); Adapter a(model); auto setup=explicitReference(a,q0); Readout r(a,setup);
    QVERIFY(r.reference().equilibriumCertified); QVERIFY(r.reference().explicitWorkpoint);
    QVERIFY(r.reference().freeForceResidualNorm<=r.reference().forceTolerance);
    QCOMPARE(r.reference().electricalResiduals.size(),std::size_t(4));
    // Preload contributes the real transverse geometric tangent N/L.
    QVERIFY(near(r.reference().assembly.tangent.findValue(4,4),force/(1+q0)));
    QVERIFY(near(r.reference().assembly.tangent.findValue(3,3),50+.4*.4/.5));
    auto bad=setup; bad.workpoint.nodes[0].f0[0]+=.01; QVERIFY(rejects([&]{Readout invalid(a,bad);}));
    bad=setup; bad.workpoint.nodes[0].q0[0]+=.01; QVERIFY(rejects([&]{Readout invalid(a,bad);}));
    bad=setup; bad.workpoint.nodes[1].q0[0]=.001; QVERIFY(rejects([&]{Readout invalid(a,bad);}));
    bad=setup; bad.workpoint.ports[2].Q0+=1e-8; QVERIFY(rejects([&]{Readout invalid(a,bad);}));
    bad=setup; bad.workpoint.ports[3].Q0+=1e-8; QVERIFY(rejects([&]{Readout invalid(a,bad);}));
    bad=setup; bad.workpoint.ports[3].V0+=.1; QVERIFY(rejects([&]{Readout invalid(a,bad);}));
    bad=setup; bad.workpoint.nodes.pop_back(); QVERIFY(rejects([&]{Readout invalid(a,bad);}));
    auto state=a.driver().initialize(0,r.reference().displacements,std::vector<double>(6,0),r.reference().inputs);
    for(const auto &s:r.observe(state)) QCOMPARE(s.deltaValue,0.0);
    // A sub-picocoulomb error is not accepted by a one-coulomb tolerance floor.
    bad=setup; bad.workpoint.ports[3].Q0+=1e-14; QVERIFY(rejects([&]{Readout invalid(a,bad);}));
} REPORT_PIEZO_EXCEPTIONS

void TestStructuralPiezoNativeReadout::frozenSensorBoundariesAndAtomicFailure() try
{
    auto model=fixture(); Adapter a(model); Readout r(a,catalog()); const auto original=a.initialize();
    for(const QString id:{QString("open"),QString("short"),QString("amplifier")}) {
        auto inputs=original.appliedInputs; const auto i=a.maps().memberByPortId.at(id);
        if(id=="open") inputs.memberInputs[i]=RStructuralTrussCharge{.031}; else inputs.memberInputs[i]=RStructuralTrussVoltage{.001};
        // Even a driver-certified physically consistent replacement cannot
        // change the frozen sensor-only boundary.
        const auto altered=a.driver().initialize(original.time,original.displacements,original.velocities,inputs);
        QVERIFY(rejects([&]{r.observe(altered);}));
    }
    auto bad=catalog(); bad.ports[1].Cf=1e-320; Readout overflow(a,bad);
    auto inputs=original.appliedInputs;
    auto q=original.displacements; q[3]=.03;
    const auto changed=a.driver().initialize(0,q,original.velocities,inputs);
    QVERIFY(rejects([&]{overflow.observe(changed);}));
    QCOMPARE(original.displacements,a.initialize().displacements);
    QCOMPARE(r.observe(original).size(),std::size_t(4));
    bad=catalog(); bad.ports[0].portId="missing"; QVERIFY(rejects([&]{Readout invalid(a,bad);}));
} REPORT_PIEZO_EXCEPTIONS

void TestStructuralPiezoNativeReadout::rotationAndIndependentElectrodeSign() try
{
    auto xModel=fixture(0),yModel=fixture(1); Adapter x(xModel),y(yModel); Readout xr(x,catalog()),yr(y,catalog());
    const auto xs=x.trial(x.initialize(),0,0).state,ys=y.trial(y.initialize(),0,0).state;
    const auto xv=xr.observe(xs),yv=yr.observe(ys);
    for(std::size_t i=0;i<xv.size();++i) {
        QVERIFY(near(xv[i].voltage,yv[i].voltage)); QVERIFY(near(xv[i].charge,yv[i].charge));
        QVERIFY(near(xv[i].deltaValue,yv[i].deltaValue));
    }
    // Reverse every mechanical endpoint while retaining signed electrode laws.
    auto reversed=fixture();
    for(uint i=0;i<5;++i) { reversed.getElement(i).setNodeId(0,1); reversed.getElement(i).setNodeId(1,0); }
    reversed.getProblemSetup().getStructuralTrussDynamicsSetup().bindToModel(reversed);
    Adapter z(reversed); Readout zr(z,catalog()); const auto zv=zr.observe(z.trial(z.initialize(),0,0).state);
    for(std::size_t i=0;i<xv.size();++i) QVERIFY(near(xv[i].value,zv[i].value));
} REPORT_PIEZO_EXCEPTIONS

void TestStructuralPiezoNativeReadout::smallIncrementSurvivesLargeBias() try
{
    auto model=fixture(0,.7,0,1e10); Adapter a(model); Readout r(a,catalog());
    auto q=a.initialize().displacements; q[3]=1e-12;
    const auto state=a.driver().initialize(0,q,std::vector<double>(6,0),a.initialCompleteInputs());
    const auto observed=r.observe(state); const auto &actuator=sample(observed,"actuator");
    // Absolute rounded charge loses the tiny signal; the displacement-based
    // increment retains it without subtracting two huge biased readings.
    QCOMPARE(actuator.charge,r.reference().assembly.memberResponses[a.maps().memberByPortId.at("actuator")].electrical->conductorCharge);
    QVERIFY(near(actuator.deltaCharge,-.3e-12));
    QVERIFY(near(sample(observed,"open").deltaVoltage,-.4e-12/.5));
    QVERIFY(near(sample(observed,"amplifier").deltaValue,-.6e-12/.25));
} REPORT_PIEZO_EXCEPTIONS

void TestStructuralPiezoNativeReadout::splitRestorePreservesReference() try
{
    auto model=fixture(); Adapter a(model); Readout r(a,catalog()); auto continuous=a.initialize();
    for(uint n=0;n<6;++n) continuous=a.trial(continuous,n,0).state;
    auto split=a.initialize(); for(uint n=0;n<3;++n) split=a.trial(split,n,0).state;
    const auto record=RStructuralTrussDynamicsRestartAdapter::build(a,split,3,0);
    Adapter other(model); Readout restoredReadout(other,catalog());
    auto restored=RStructuralTrussDynamicsRestartAdapter::restore(other,record).state;
    for(uint n=3;n<6;++n) restored=other.trial(restored,n,0).state;
    const auto full=r.observe(continuous),resumed=restoredReadout.observe(restored);
    for(std::size_t i=0;i<full.size();++i) {
        QCOMPARE(full[i].value,resumed[i].value); QCOMPARE(full[i].deltaValue,resumed[i].deltaValue);
    }
    QCOMPARE(r.reference().displacements,restoredReadout.reference().displacements);
} REPORT_PIEZO_EXCEPTIONS

QTEST_APPLESS_MAIN(TestStructuralPiezoNativeReadout)
#include "tst_structural_piezo_native_readout.moc"

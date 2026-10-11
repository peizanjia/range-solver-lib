#include <QtTest>
#include <QTemporaryDir>
#include <QDataStream>
#include <algorithm>
#include <cmath>
#include <cstring>
#include <limits>
#include <type_traits>
#include <rbl_error.h>
#include <rml_model.h>
#include "rstructural_truss_dynamics_model_adapter.h"
#include "rstructural_truss_dynamics_restart_adapter.h"
#include "rstructural_static_truss_model.h"
#include "rstructural_cable_model_adapter.h"
#include "rstructural_piezo_model_adapter.h"

namespace {
using Adapter=RStructuralTrussDynamicsModelAdapter;
using Restart=RStructuralTrussDynamicsRestartAdapter;
using State=RStructuralTrussDynamicsState;
bool near(double a,double b,double tolerance=2e-10)
{return std::abs(a-b)<=tolerance*(1+std::abs(b));}
template<class F> bool rejects(F f)
{try{f();}catch(const RError &){return true;}catch(const std::exception &){return true;}return false;}
void component(RBoundaryCondition &bc,RVariableType type,double value,bool enabled=true)
{
    auto &c=bc.getComponent(bc.findComponentPosition(type));c.setEnabled(enabled);c.clear();if(enabled)c.add(0,value);
}
RBoundaryCondition translations(bool x,bool y,bool z)
{
    RBoundaryCondition bc(R_BOUNDARY_CONDITION_DISPLACEMENT);
    component(bc,R_VARIABLE_DISPLACEMENT_X,0,x);component(bc,R_VARIABLE_DISPLACEMENT_Y,0,y);component(bc,R_VARIABLE_DISPLACEMENT_Z,0,z);
    return bc;
}
void point(RModel &m,uint node,uint group)
{
    RElement e(R_ELEMENT_POINT);e.setNodeId(0,node);m.addElement(e,true,group);
    auto &g=m.getPoint(group);if(g.hasEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION)){
        auto &gravity=g.getEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
        for(auto type:{R_VARIABLE_G_ACCELERATION_X,R_VARIABLE_G_ACCELERATION_Y,R_VARIABLE_G_ACCELERATION_Z}){
            auto &c=gravity.getComponent(gravity.findComponentPosition(type));c.setEnabled(true);c.clear();c.add(0,0);
        }
    }
}
void initial(RModel &m,uint group,double q,double v,uint axis=0)
{
    RInitialCondition displacement(R_INITIAL_CONDITION_DISPLACEMENT),velocity(R_INITIAL_CONDITION_VELOCITY);
    const RVariableType qs[]={R_VARIABLE_DISPLACEMENT_X,R_VARIABLE_DISPLACEMENT_Y,R_VARIABLE_DISPLACEMENT_Z};
    const RVariableType vs[]={R_VARIABLE_VELOCITY_X,R_VARIABLE_VELOCITY_Y,R_VARIABLE_VELOCITY_Z};
    for(uint i=0;i<3;++i){
        auto &a=displacement.getComponent(displacement.findComponentPosition(qs[i]));a.setEnabled(true);a.clear();a.add(0,i==axis?q:0);
        auto &b=velocity.getComponent(velocity.findComponentPosition(vs[i]));b.setEnabled(true);b.clear();b.add(0,i==axis?v:0);
    }
    m.getPoint(group).addInitialCondition(displacement);m.getPoint(group).addInitialCondition(velocity);
}
void force(RModel &m,uint group,double f)
{
    RBoundaryCondition bc(R_BOUNDARY_CONDITION_FORCE);component(bc,R_VARIABLE_FORCE_X,f);
    component(bc,R_VARIABLE_FORCE_Y,0,false);component(bc,R_VARIABLE_FORCE_Z,0,false);m.getPoint(group).addBoundaryCondition(bc);
}
void material(RLine &line,double young=40,double density=4)
{
    RMaterial mat(RMaterial::Solid);
    RMaterialProperty e(RMaterialProperty::ModulusOfElasiticity);e.add(293.15,young);mat.add(e);
    RMaterialProperty rho(RMaterialProperty::Density);rho.add(293.15,density);mat.add(rho);
    line.setCrossArea(1);line.setMaterial(mat);
}
void grid(RModel &m,double dt,uint steps)
{
    m.getProblemSetup().setMassMatrixType(RMassMatrixType::Lumped);
    // These reference problems have no gravity load; RLine defaults otherwise
    // contain a nonzero environment acceleration which this route rejects.
    for(uint i=0;i<m.getNElementGroups();++i){auto *group=m.getElementGroupPtr(i);
        if(!group || !group->hasEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION))continue;
        auto &gravity=group->getEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
        for(auto type:{R_VARIABLE_G_ACCELERATION_X,R_VARIABLE_G_ACCELERATION_Y,R_VARIABLE_G_ACCELERATION_Z}){
            auto &c=gravity.getComponent(gravity.findComponentPosition(type));c.setEnabled(true);c.clear();c.add(0,0);
        }
    }
    m.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));auto &d=m.getProblemSetup().getStructuralDynamicsSetup();
    d.setEnabled(true);d.setTimeStep(dt);d.setNSteps(steps);d.setOutputFrequency(1);
    auto &t=m.getTimeSolver();t.setEnabled(true);t.setInputStartTime(dt);t.setInputTimeStepSize(dt);t.setInputNTimeSteps(steps);
    std::vector<double> times;for(uint n=1;n<=steps;++n)times.push_back(double(n)*dt);t.setTimes(times);
}
void bind(RModel &m)
{m.getProblemSetup().getStructuralTrussDynamicsSetup().bindToModel(m);}
RModel axial(int kind=0,double dt=.004,uint steps=100,bool fixedActive=false)
{
    RModel m;m.addNode(RNode(0,0,0));m.addNode(RNode(1,0,0));m.setNodePhysicalId(0,101);m.setNodePhysicalId(1,202);
    RElement e(R_ELEMENT_TRUSS1);e.setNodeId(0,0);e.setNodeId(1,1);m.addElement(e,true,0);m.setElementPhysicalId(0,501);
    if(!fixedActive)material(m.getLine(0));else m.getLine(0).setCrossArea(1);
    m.getLine(0).addBoundaryCondition(translations(false,true,true));
    point(m,0,0);m.getPoint(0).addBoundaryCondition(translations(true,true,true));
    point(m,1,1);initial(m,1,.035,.08);if(kind!=3)force(m,1,1.2);
    grid(m,dt,steps);auto &s=m.getProblemSetup().getStructuralTrussDynamicsSetup();s.enabled=true;
    s.absoluteForceTolerance=1e-11;s.relativeForceTolerance=1e-11;s.linearRelativeTolerance=1e-12;
    if(kind==0)s.passiveMembers={{501,1}};
    else if(kind==3)s.cableMembers={{501,3,0}};
    else {
        RStructuralTrussActiveMember a;a.physicalElementId=501;a.portId="electrode 19";a.referenceLength=1;
        a.axialStiffness=40;a.coupling=2;a.capacitance=.5;a.initialValue=kind==1?.7:.15;
        a.kind=kind==1?RStructuralTrussElectricalKind::Voltage:RStructuralTrussElectricalKind::ConductorCharge;
        if(fixedActive){a.massSource=RStructuralTrussMassSource::FixedMemberMass;a.fixedMass=4;}
        s.activeMembers={a};
    }
    bind(m);return m;
}
RModel mixed(double dt=.002,uint steps=12)
{
    RModel m;m.addNode(RNode(0,0,0));m.addNode(RNode(1,0,0));m.setNodePhysicalId(0,101);m.setNodePhysicalId(1,202);
    // Deliberately nonmonotone physical member IDs: array positions cannot be used as commands.
    const uint ids[]={901,503,777,609};
    for(uint i=0;i<4;++i){RElement e(R_ELEMENT_TRUSS1);e.setNodeId(0,0);e.setNodeId(1,1);m.addElement(e,true,i);
        m.setElementPhysicalId(i,ids[i]);material(m.getLine(i),10,1);m.getLine(i).addBoundaryCondition(translations(false,true,true));}
    point(m,0,0);m.getPoint(0).addBoundaryCondition(translations(true,true,true));point(m,1,1);initial(m,1,.035,.08);force(m,1,1.2);
    grid(m,dt,steps);auto &s=m.getProblemSetup().getStructuralTrussDynamicsSetup();s.enabled=true;
    s.absoluteForceTolerance=1e-11;s.relativeForceTolerance=1e-11;s.linearRelativeTolerance=1e-12;
    s.passiveMembers={{901,1}};
    for(int i=0;i<2;++i){RStructuralTrussActiveMember a;a.physicalElementId=ids[i+1];a.portId=i==0?"voltage independent":"charge independent";
        a.referenceLength=1;a.axialStiffness=10;a.coupling=.2;a.capacitance=.5;a.initialValue=i==0?.7:.15;
        a.kind=i==0?RStructuralTrussElectricalKind::Voltage:RStructuralTrussElectricalKind::ConductorCharge;s.activeMembers.push_back(a);}
    s.cableMembers={{609,1,1}};
    s.heldCommands={{1,{{"charge independent",.2},{"voltage independent",.8}},{{609,.99}},{{202,{{.4,0,0}}}}},
                    {3,{},{},{{202,{{2,0,0}}}}},
                    {5,{{"voltage independent",-.3},{"charge independent",-.1}},{{609,1.02}},{} }};
    bind(m);return m;
}
struct End {State state;uint completed=0,cursor=0;};
End run(const Adapter &a,uint target,End e)
{
    while(e.completed<target){auto t=a.trial(e.state,e.completed,e.cursor);e.state=t.state;e.cursor=t.nextCommandIndex;++e.completed;}return e;
}
End fresh(const Adapter &a){return {a.initialize(),0,0};}
QByteArray stateBytes(const State &s)
{
    QByteArray b;QDataStream out(&b,QIODevice::WriteOnly);out.setFloatingPointPrecision(QDataStream::DoublePrecision);
    const auto real=[&](double v){quint64 bits;std::memcpy(&bits,&v,sizeof(bits));out<<bits;};
    const auto vec=[&](const std::vector<double>&v){out<<quint64(v.size());for(double x:v)real(x);};
    real(s.time);vec(s.displacements);vec(s.velocities);vec(s.accelerations);vec(s.reactions);vec(s.appliedInputs.externalForces);
    out<<quint64(s.appliedInputs.memberInputs.size());
    for(const auto &i:s.appliedInputs.memberInputs){out<<quint64(i.index());std::visit([&](const auto &v){
        if constexpr(!std::is_same_v<std::decay_t<decltype(v)>,RStructuralTrussNoInput>)real(v.value);},i);}
    real(s.storedEnergy);real(s.kineticEnergy);real(s.initialTotalEnergy);real(s.cumulativeExternalWork);
    real(s.cumulativeSupportWork);real(s.cumulativeElectricalWork);real(s.cumulativeLengthWork);out<<s.systemIdentity;return b;
}
void exactRecord(const RStructuralTrussDynamicsRestartData &a,const RStructuralTrussDynamicsRestartData &b)
{
#define EQ(f) QVERIFY(a.f==b.f)
    EQ(available);EQ(completedSteps);EQ(runTargetSteps);EQ(nextCommandIndex);EQ(orderedPhysicalNodeIds);
    EQ(orderedMemberPhysicalIds);EQ(orderedPortIds);EQ(sourceFingerprint);EQ(systemIdentity);EQ(memberInputKinds);
#undef EQ
#define BITS(f) QVERIFY(std::memcmp(&a.f,&b.f,sizeof(double))==0)
    BITS(initialTime);BITS(time);BITS(fixedTimeStep);BITS(storedEnergy);BITS(kineticEnergy);BITS(initialTotalEnergy);
    BITS(cumulativeExternalWork);BITS(cumulativeSupportWork);BITS(cumulativeElectricalWork);BITS(cumulativeRestLengthWork);
#undef BITS
    const auto same=[](const std::vector<double>&x,const std::vector<double>&y){if(x.size()!=y.size())return false;
        for(size_t i=0;i<x.size();++i)if(std::memcmp(&x[i],&y[i],sizeof(double)))return false;return true;};
#define VEC(f) QVERIFY(same(a.f,b.f))
    VEC(q);VEC(v);VEC(a);VEC(reactions);VEC(externalForce);VEC(memberInputValues);VEC(portVoltages);VEC(portCharges);VEC(portEnthalpies);
#undef VEC
}
}
#define REPORT_FE_EXCEPTIONS \
    catch(const RError &error) { QFAIL(qPrintable(QString("RError: %1").arg(error.getMessage()))); } \
    catch(const std::exception &error) { QFAIL(error.what()); } \
    catch(...) { QFAIL("Unknown exception during native mixed FE behavior"); }
class TestStructuralTrussDynamicsFE:public QObject {
    Q_OBJECT
private slots:
    void oscillators_data();void oscillators();
    void stiffPiezoObjectivity_data();void stiffPiezoObjectivity();
    void massSourcesAndSourceInitialConditions();
    void absoluteCommandsAndIndependentPorts();
    void continuousSplitDisk_data();void continuousSplitDisk();
    void canonicalMovingAnchor_data();void canonicalMovingAnchor();
    void obliquePhysicalModel();
    void staleSource_data();void staleSource();
    void corruptRestore_data();void corruptRestore();
    void legacyAndUnsupportedGates();
    void tinySourceContradictions_data();void tinySourceContradictions();
    void tinySlackSupportRestartCorruption();
};
void TestStructuralTrussDynamicsFE::oscillators_data()
{QTest::addColumn<int>("kind");QTest::newRow("passive")<<0;QTest::newRow("voltage")<<1;QTest::newRow("conductor-charge")<<2;}
void TestStructuralTrussDynamicsFE::oscillators() try
{
    QFETCH(int,kind);auto m=axial(kind);Adapter a(m);const auto initial=a.initialize();
    QCOMPARE(initial.displacements[3],.035);QCOMPARE(initial.velocities[3],.08);QCOMPARE(a.system().nodalLumpedMass[1],2.0);
    const double k=kind==2?48:40,offset=kind==1?1.4:(kind==2?.6:0),eq=(1.2-offset)/k,w=std::sqrt(k/2);
    const auto end=run(a,100,fresh(a));const double t=.4;
    const double q=eq+(.035-eq)*std::cos(w*t)+.08/w*std::sin(w*t);
    const double v=-(.035-eq)*w*std::sin(w*t)+.08*std::cos(w*t);
    QVERIFY(near(end.state.displacements[3],q,5e-6));QVERIFY(near(end.state.velocities[3],v,2e-5));
    QVERIFY(near(end.state.accelerations[3],(1.2-offset-k*end.state.displacements[3])/2));
    const double delta=end.state.displacements[3],voltage=kind==1?.7:(kind==2?(.15+2*delta)/.5:0);
    const double energy=20*delta*delta+(kind?.25*voltage*voltage:0);
    QVERIFY(near(end.state.storedEnergy,energy));
    auto candidate=a.candidateModel(m,end.state,100,end.cursor);
    QCOMPARE(candidate.getTimeSolver().getCurrentTimeStep(),99u);QCOMPARE(candidate.getTimeSolver().getComputedTime(),a.targetTime(100));
    const auto &forceVar=candidate.getVariable(candidate.findVariable(R_VARIABLE_FORCE));
    QVERIFY(near(forceVar.getValue(0,1),k*delta+offset));
    QVERIFY(std::abs(forceVar.getValue(0,1)-end.state.reactions[3])>.01);
    if(kind){const auto r=Restart::build(a,end.state,100,end.cursor);QCOMPARE(r.orderedPortIds[0],QString("electrode 19"));
        QVERIFY(near(r.portVoltages[0],voltage));QVERIFY(near(r.portCharges[0],.5*voltage-2*delta));
        QVERIFY(near(r.portEnthalpies[0],20*delta*delta+2*delta*voltage-.25*voltage*voltage));}
} REPORT_FE_EXCEPTIONS
void TestStructuralTrussDynamicsFE::stiffPiezoObjectivity_data()
{
    QTest::addColumn<double>("origin");
    QTest::newRow("origin-0-m")<<0.0;QTest::newRow("origin-3-m")<<3.0;QTest::newRow("origin-1000-m")<<1000.0;
}
void TestStructuralTrussDynamicsFE::stiffPiezoObjectivity() try
{
    QFETCH(double,origin);
    constexpr double stiffness=1e9,coupling=-.1,capacitance=1e-9,voltage=10,memberMass=.1;
    constexpr double freeMass=memberMass/2,equilibrium=1e-9,coarseDt=1e-6,fineDt=.5e-6;
    constexpr double duration=20*coarseDt,objectivityMeters=3e-15;
    const auto fixture=[&](double shift,double dt,uint steps) {
        auto m=axial(1,dt,steps,true);
        m.getNode(0).setX(shift);m.getNode(1).setX(shift+1);
        initial(m,1,0,0);force(m,1,0);
        auto &active=m.getProblemSetup().getStructuralTrussDynamicsSetup().activeMembers[0];
        active.referenceLength=1;active.axialStiffness=stiffness;active.coupling=coupling;
        active.capacitance=capacitance;active.initialValue=voltage;active.fixedMass=memberMass;
        bind(m);return m;
    };
    auto coarseModel=fixture(origin,coarseDt,20),fineModel=fixture(origin,fineDt,40);
    Adapter coarse(coarseModel),fine(fineModel);
    QCOMPARE(fine.system().nodalLumpedMass[1],freeMass);
    const double omega=std::sqrt(stiffness/freeMass);
    const double exactQ=equilibrium*(1-std::cos(omega*duration));
    const double exactV=equilibrium*omega*std::sin(omega*duration);
    const auto coarseEnd=run(coarse,20,fresh(coarse)),fineEnd=run(fine,40,fresh(fine));
    const double coarseError=std::abs(coarseEnd.state.displacements[3]-exactQ);
    const double fineError=std::abs(fineEnd.state.displacements[3]-exactQ);
    // Metres and metres/second: no one-metre normalization floor.
    QVERIFY2(coarseError<2e-12,qPrintable(QString("coarse displacement error %1 m").arg(coarseError,0,'g',17)));
    QVERIFY2(fineError<5e-13,qPrintable(QString("fine displacement error %1 m").arg(fineError,0,'g',17)));
    QVERIFY2(fineError<.30*coarseError,"Expected second-order Newmark refinement in physical displacement units");
    QVERIFY(std::abs(fineEnd.state.velocities[3]-exactV)<2e-7);
    auto zeroModel=fixture(0,fineDt,40);Adapter zero(zeroModel);
    auto translated=fresh(fine),baseline=fresh(zero);
    for(uint step=1;step<=40;++step) {
        translated=run(fine,step,translated);baseline=run(zero,step,baseline);
        const double difference=std::abs(translated.state.displacements[3]-baseline.state.displacements[3]);
        QVERIFY2(difference<=objectivityMeters,qPrintable(QString("origin %1 m step %2: translation changed q by %3 m")
            .arg(origin,0,'g',17).arg(step).arg(difference,0,'g',17)));
        const double extension=translated.state.displacements[3]-translated.state.displacements[0];
        const double expectedForce=stiffness*extension+coupling*voltage;
        const double expectedCharge=capacitance*voltage-coupling*extension;
        auto candidate=fine.candidateModel(fineModel,translated.state,step,translated.cursor);
        const auto &internal=candidate.getVariable(candidate.findVariable(R_VARIABLE_FORCE));
        const double forceError=std::abs(internal.getValue(0,1)-expectedForce);
        QVERIFY2(forceError<=stiffness*objectivityMeters+1e-12,
            qPrintable(QString("constitutive internal force error %1 N").arg(forceError,0,'g',17)));
        const auto record=Restart::build(fine,translated.state,step,translated.cursor);
        const double chargeError=std::abs(record.portCharges[0]-expectedCharge);
        QVERIFY2(chargeError<=std::abs(coupling)*objectivityMeters+1e-23,
            qPrintable(QString("constitutive port charge error %1 C").arg(chargeError,0,'g',17)));
        QCOMPARE(record.portVoltages[0],voltage);
    }
} REPORT_FE_EXCEPTIONS
void TestStructuralTrussDynamicsFE::massSourcesAndSourceInitialConditions() try
{
    auto m=axial();m.getProblemSetup().getStructuralTrussDynamicsSetup().passiveMembers[0].referenceLength=.8;bind(m);Adapter a(m);
    QCOMPARE(a.system().nodalLumpedMass[0],2.0);QCOMPARE(a.system().nodalLumpedMass[1],2.0); // rho*A*mesh chord / 2, not Lrest.
    auto fixed=axial(1,.004,100,true);QCOMPARE(fixed.getLine(0).getMaterial().size(),0u);Adapter calibrated(fixed);
    QCOMPARE(calibrated.system().nodalLumpedMass[1],2.0);QCOMPARE(calibrated.initialize().displacements[3],.035);
    auto source=axial();for(auto type:{R_VARIABLE_DISPLACEMENT,R_VARIABLE_VELOCITY,R_VARIABLE_ELECTRIC_POTENTIAL}) {
        auto p=source.addVariable(RVariable(type));auto &v=source.getVariable(p);v.setApplyType(R_VARIABLE_APPLY_NODE);v.resize(3,source.getNNodes());v.setValue(0,1,999);
    }
    const RVariableType oldElectrical[]={R_VARIABLE_ELECTRIC_FIELD,R_VARIABLE_ELECTRIC_ENERGY,R_VARIABLE_CHARGE_DENSITY,R_VARIABLE_CURRENT_DENSITY};
    for(auto type:oldElectrical){const auto p=source.addVariable(RVariable(type));auto &v=source.getVariable(p);
        v.setApplyType(R_VARIABLE_APPLY_NODE);v.resize(3,source.getNNodes());v.setValue(0,1,123);}
    Adapter stale(source);QCOMPARE(stale.initialize().displacements[3],.035);QCOMPARE(stale.initialize().velocities[3],.08);
    const auto first=stale.trial(stale.initialize(),0,0);const auto candidate=stale.candidateModel(source,first.state,1,first.nextCommandIndex);
    QCOMPARE(candidate.findVariable(R_VARIABLE_ELECTRIC_POTENTIAL),RConstants::eod);
    for(auto type:oldElectrical){QCOMPARE(candidate.findVariable(type),RConstants::eod);
        QCOMPARE(source.getVariable(source.findVariable(type)).getValue(0,1),123.0);}
    QCOMPARE(source.getVariable(source.findVariable(R_VARIABLE_DISPLACEMENT)).getValue(0,1),999.0);
} REPORT_FE_EXCEPTIONS
void TestStructuralTrussDynamicsFE::absoluteCommandsAndIndependentPorts() try
{
    auto m=mixed();Adapter a(m);auto left=a.initialize();QCOMPARE(a.maps().orderedMemberPhysicalIds,std::vector<uint>({901,503,777,609}));
    QCOMPARE(a.maps().memberByPortId.at("charge independent"),size_t(2));
    const auto action=a.actionForStep(left.appliedInputs,0,1);QCOMPARE(action.nextCommandIndex,1u);
    QCOMPARE(std::get<RStructuralTrussVoltage>(action.completeInputs.memberInputs[1]).value,.8);
    QCOMPARE(std::get<RStructuralTrussCharge>(action.completeInputs.memberInputs[2]).value,.2);
    QCOMPARE(std::get<RStructuralTrussCableLength>(action.completeInputs.memberInputs[3]).value,.99);
    QCOMPARE(action.completeInputs.externalForces[3],.4);const auto oldMass=a.system().nodalLumpedMass;
    const auto trial=a.trial(left,0,0);
    const double displacement=.035,oldQVoltage=(.15+.2*displacement)/.5,newQVoltage=(.2+.2*displacement)/.5;
    const double electrical=.25*(.8*.8-.7*.7)+.25*(newQVoltage*newQVoltage-oldQVoltage*oldQVoltage);
    const auto cableEnergy=[](double l0){const double stiffness=10/l0,n=std::max(0.0,1+stiffness*(1.035-l0));return n*n/(2*stiffness);};
    QVERIFY(near(trial.result.jumpWork.electricalWork,electrical));
    QVERIFY(near(trial.result.jumpWork.lengthWork,cableEnergy(.99)-cableEnergy(1)));
    QVERIFY(a.system().nodalLumpedMass==oldMass);
    auto step2=a.actionForStep(trial.state.appliedInputs,1,2);QCOMPARE(step2.nextCommandIndex,1u);QCOMPARE(step2.completeInputs.externalForces[3],.4);
    auto step3=a.actionForStep(step2.completeInputs,1,3);QCOMPARE(step3.completeInputs.externalForces[3],2.0); // absolute replacement, not .4+2.
    QCOMPARE(std::get<RStructuralTrussVoltage>(step3.completeInputs.memberInputs[1]).value,.8);
    const auto r=Restart::build(a,trial.state,1,trial.nextCommandIndex);QCOMPARE(r.orderedPortIds,std::vector<QString>({"voltage independent","charge independent"}));
    const double q=trial.state.displacements[3];QVERIFY(near(r.portVoltages[0],.8));QVERIFY(near(r.portCharges[0],.4-.2*q));
    QVERIFY(near(r.portCharges[1],.2));QVERIFY(near(r.portVoltages[1],(.2+.2*q)/.5));
    for(size_t p=0;p<2;++p){const double v=r.portVoltages[p];QVERIFY(near(r.portEnthalpies[p],5*q*q+.2*q*v-.25*v*v));}
} REPORT_FE_EXCEPTIONS
void TestStructuralTrussDynamicsFE::continuousSplitDisk_data()
{
    QTest::addColumn<bool>("binary");QTest::addColumn<int>("kind");
    for(bool binary:{false,true})for(int kind:{0,1,2,4})QTest::newRow(qPrintable(QString("%1-kind-%2").arg(binary?"rbm":"ram").arg(kind)))<<binary<<kind;
}
void TestStructuralTrussDynamicsFE::continuousSplitDisk() try
{
    QFETCH(bool,binary);QFETCH(int,kind);QTemporaryDir dir;QVERIFY(dir.isValid());
    auto source=kind==4?mixed(.0034567890123456789,24):axial(kind,.0034567890123456789,24);Adapter continuous(source);
    const auto expected=run(continuous,24,fresh(continuous));const auto split=run(continuous,9,fresh(continuous));
    const auto record=Restart::build(continuous,split.state,9,split.cursor);
    auto checkpoint=continuous.candidateModel(source,split.state,9,split.cursor);
    QCOMPARE(checkpoint.getTimeSolver().getCurrentTimeStep(),8u);QCOMPARE(checkpoint.getTimeSolver().getComputedTime(),record.time);
    const auto path=dir.filePath("split."+RModel::getDefaultFileExtension(binary));checkpoint.write(path,false);
    RModel disk;disk.read(path);exactRecord(disk.getProblemSetup().getStructuralTrussDynamicsRestartData(),record);
    disk.getProblemSetup().setRestart(true);Adapter resumed(disk);QCOMPARE(resumed.targetSteps(),24u);QCOMPARE(resumed.nominalTimeStep(),.0034567890123456789);
    const auto restored=Restart::restore(resumed,disk.getProblemSetup().getStructuralTrussDynamicsRestartData());
    QCOMPARE(stateBytes(restored.state),stateBytes(split.state));
    const auto actual=run(resumed,24,{restored.state,restored.completedSteps,restored.nextCommandIndex});
    QCOMPARE(stateBytes(actual.state),stateBytes(expected.state));QCOMPARE(actual.cursor,expected.cursor);
    exactRecord(Restart::build(resumed,actual.state,24,actual.cursor),Restart::build(continuous,expected.state,24,expected.cursor));
    // Display fields are never authoritative restore state.
    disk.getVariable(disk.findVariable(R_VARIABLE_DISPLACEMENT)).setValue(0,1,9999);
    resumed.validateCurrentSource(disk);QCOMPARE(stateBytes(Restart::restore(resumed,record).state),stateBytes(split.state));
} REPORT_FE_EXCEPTIONS
void TestStructuralTrussDynamicsFE::canonicalMovingAnchor_data()
{QTest::addColumn<uint>("steps");QTest::newRow("1000")<<1000u;QTest::newRow("10000")<<10000u;}
void TestStructuralTrussDynamicsFE::canonicalMovingAnchor() try
{
    QFETCH(uint,steps);const double dt=.00012345678901234567,finish=double(steps)*dt;
    auto m=axial(3,dt,steps);m.getProblemSetup().getStructuralInteractionSetup().getPrescribedMotions().push_back(
        {"moving anchor",{0,RDofType::UX},{{0,0,.001,0},{finish,.001*finish,.001,0}}});bind(m);Adapter a(m);auto e=fresh(a);
    while(e.completed<steps){const auto t=a.trial(e.state,e.completed,e.cursor);e.state=t.state;e.cursor=t.nextCommandIndex;++e.completed;
        QCOMPARE(e.state.time,double(e.completed)*dt);}
    QCOMPARE(e.state.time,finish);QVERIFY(near(e.state.displacements[0],.001*finish,1e-12));
    QVERIFY(near(e.state.velocities[0],.001,1e-12));QVERIFY(near(e.state.accelerations[0],0,1e-11));
    QVERIFY(near(e.state.displacements[3],.035+.08*finish,2e-10));
    auto candidate=a.candidateModel(m,e.state,steps,e.cursor);QCOMPARE(candidate.getTimeSolver().getComputedTime(),finish);
    const auto restored=Restart::restore(a,Restart::build(a,e.state,steps,e.cursor));QCOMPARE(stateBytes(restored.state),stateBytes(e.state));
} REPORT_FE_EXCEPTIONS
void TestStructuralTrussDynamicsFE::obliquePhysicalModel() try
{
    RModel m;m.addNode(RNode(-1,0,0));m.addNode(RNode(1,0,0));m.addNode(RNode(0,.5,0));const double rest=std::sqrt(1.25);
    for(uint i=0;i<2;++i){RElement e(R_ELEMENT_TRUSS1);e.setNodeId(0,i);e.setNodeId(1,2);m.addElement(e,true,i);m.setElementPhysicalId(i,701+i);
        material(m.getLine(i),40,2/rest);m.getLine(i).addBoundaryCondition(translations(true,false,true));}
    for(uint i=0;i<2;++i){point(m,i,i);m.getPoint(i).addBoundaryCondition(translations(true,true,true));}
    point(m,2,2);initial(m,2,.05,.07,1);grid(m,.0005,80);
    auto &s=m.getProblemSetup().getStructuralTrussDynamicsSetup();s.enabled=true;s.passiveMembers={{701,rest},{702,rest}};s.absoluteForceTolerance=1e-11;s.relativeForceTolerance=1e-11;bind(m);
    Adapter a(m);const double mass=a.system().nodalLumpedMass[2];QVERIFY(near(mass,2,1e-14));
    const auto restoring=[&](double y){const double length=std::hypot(1,y);return 80/rest*(length-rest)*y/length;};
    // Independent scalar average-acceleration Newmark solve for the free y coordinate.
    double y=.55,v=.07,acc=-restoring(y)/mass;
    for(uint step=0;step<80;++step){const double h=double(step+1)*.0005-double(step)*.0005;
        const double predictor=y+h*v+.25*h*h*acc;double next=predictor;
        for(int i=0;i<12;++i){const double length=std::hypot(1,next),f=restoring(next);
            const double tangent=80/rest*(1-rest/(length*length*length));
            const double residual=4*mass/(h*h)*(next-predictor)+f;
            next-=residual/(4*mass/(h*h)+tangent);}
        const double nextAcc=4*(next-predictor)/(h*h);v+=.5*h*(acc+nextAcc);y=next;acc=nextAcc;}
    const auto end=run(a,80,fresh(a));QVERIFY(near(.5+end.state.displacements[7],y,2e-9));QVERIFY(near(end.state.velocities[7],v,2e-8));
} REPORT_FE_EXCEPTIONS
void TestStructuralTrussDynamicsFE::staleSource_data()
{
    QTest::addColumn<int>("mutation");const char *names[]={"density","elastic-modulus","area","geometry","commands","motion","parallel","dt"};
    for(int i=0;i<8;++i)QTest::newRow(names[i])<<i;
}
void TestStructuralTrussDynamicsFE::staleSource() try
{
    QFETCH(int,mutation);auto original=axial();Adapter a(original);auto changed=original;
    switch(mutation){
    case 0:changed.getLine(0).getMaterial().get(changed.getLine(0).getMaterial().findPosition(RMaterialProperty::Density)).add(293.15,5);break;
    case 1:changed.getLine(0).getMaterial().get(changed.getLine(0).getMaterial().findPosition(RMaterialProperty::ModulusOfElasiticity)).add(293.15,41);break;
    case 2:changed.getLine(0).setCrossArea(1.1);break;case 3:changed.getNode(1).setX(1.01);break;
    case 4:changed.getProblemSetup().getStructuralTrussDynamicsSetup().heldCommands={{2,{},{},{{202,{{1,0,0}}}}}};break;
    case 5:changed.getProblemSetup().getStructuralInteractionSetup().getPrescribedMotions().push_back({"new anchor",{0,RDofType::UX},{{0,0,0,0},{.4,.01,0,0}}});break;
    case 6:changed.getProblemSetup().getParallelExecutionConfig().setMode(RParallelExecutionMode::SerialDeterministic);break;
    case 7:changed.getProblemSetup().getStructuralDynamicsSetup().setTimeStep(.005);break;}
    QVERIFY(rejects([&]{a.validateCurrentSource(changed);}));
    // A newly rebound but changed scene still cannot use the previously frozen adapter.
    bind(changed);QVERIFY(rejects([&]{a.validateCurrentSource(changed);}));
    QCOMPARE(original.getProblemSetup().getStructuralTrussDynamicsSetup().sourceFingerprint,a.sourceFingerprint());
} REPORT_FE_EXCEPTIONS
void TestStructuralTrussDynamicsFE::corruptRestore_data()
{
    QTest::addColumn<int>("mutation");const char *names[]={"system","node-map","member-map","port-map","port-kind","cursor","nominal-dt","run-target","clock","actual-electrical","actual-cable","actual-force","derived-voltage","derived-charge","derived-enthalpy","initial-energy"};
    for(int i=0;i<16;++i)QTest::newRow(names[i])<<i;
}
void TestStructuralTrussDynamicsFE::corruptRestore() try
{
    QFETCH(int,mutation);auto m=mixed();Adapter a(m);const auto e=run(a,6,fresh(a));const auto valid=Restart::build(a,e.state,6,e.cursor);auto bad=valid;
    switch(mutation){case 0:bad.systemIdentity[0]=char(bad.systemIdentity[0]^1);break;
    case 1:std::swap(bad.orderedPhysicalNodeIds[0],bad.orderedPhysicalNodeIds[1]);break;
    case 2:std::swap(bad.orderedMemberPhysicalIds[0],bad.orderedMemberPhysicalIds[1]);break;
    case 3:std::swap(bad.orderedPortIds[0],bad.orderedPortIds[1]);break;
    case 4:bad.memberInputKinds[1]=2;break;case 5:bad.nextCommandIndex=0;break;
    case 6:bad.fixedTimeStep=std::nextafter(bad.fixedTimeStep,1.0);break;case 7:++bad.runTargetSteps;break;
    case 8:bad.time=std::nextafter(bad.time,1.0);break;case 9:bad.memberInputValues[1]+=.01;break;
    case 10:bad.memberInputValues[3]+=.01;break;case 11:bad.externalForce[3]+=.01;break;
    case 12:bad.portVoltages[0]+=.01;break;case 13:bad.portCharges[1]+=.01;break;case 14:bad.portEnthalpies[0]+=.01;break;
    case 15:bad.initialTotalEnergy+=.01;break;}
    QVERIFY(rejects([&]{Restart::restore(a,bad);}));QCOMPARE(stateBytes(Restart::restore(a,valid).state),stateBytes(e.state));
} REPORT_FE_EXCEPTIONS
void TestStructuralTrussDynamicsFE::legacyAndUnsupportedGates() try
{
    auto m=axial();const auto original=m.getProblemSetup().getStructuralTrussDynamicsSetup().sourceFingerprint;
    QVERIFY(rejects([&]{RStructuralStaticTrussModel::build(m,true);}));
    // Disabled legacy validators are intentionally no-ops. Their actual solve
    // entries must refuse this enabled native model rather than run it.
    QVERIFY(rejects([&]{RStructuralCableModelAdapter::solve(m);}));QVERIFY(rejects([&]{RStructuralPiezoModelAdapter::solve(m);}));
    auto oldCable=m;oldCable.getProblemSetup().getStructuralCableSetup().enabled=true;bind(oldCable);
    QVERIFY(rejects([&]{Adapter rejected(oldCable);}));
    auto oldPiezo=m;oldPiezo.getProblemSetup().getStructuralPiezoSetup().enabled=true;bind(oldPiezo);
    QVERIFY(rejects([&]{Adapter rejected(oldPiezo);}));
    auto consistent=m;consistent.getProblemSetup().setMassMatrixType(RMassMatrixType::Consistent);bind(consistent);
    QVERIFY(rejects([&]{Adapter rejected(consistent);}));
    QCOMPARE(consistent.getProblemSetup().getMassMatrixType(),RMassMatrixType::Consistent);
    QCOMPARE(m.getProblemSetup().getStructuralTrussDynamicsSetup().sourceFingerprint,original);
    QVERIFY(m.getProblemSetup().getStructuralTrussDynamicsSetup().enabled);
} REPORT_FE_EXCEPTIONS
void TestStructuralTrussDynamicsFE::tinySourceContradictions_data()
{
    QTest::addColumn<int>("mutation");QTest::newRow("duplicate-displacement-1e-14-m")<<0;
    QTest::newRow("duplicate-velocity-1e-14-m-s")<<1;QTest::newRow("support-IC-1e-14-m")<<2;
    QTest::newRow("fixed-BC-conflict-1e-14-m")<<3;
}
void TestStructuralTrussDynamicsFE::tinySourceContradictions() try
{
    QFETCH(int,mutation);auto m=axial(1,.000001,20,true);initial(m,1,0,0);force(m,1,0);
    auto &active=m.getProblemSetup().getStructuralTrussDynamicsSetup().activeMembers[0];active.axialStiffness=1e9;active.initialValue=0;
    if(mutation<2){point(m,1,2);initial(m,2,mutation==0?1e-14:0,mutation==1?1e-14:0);}
    else if(mutation==2)initial(m,0,1e-14,0);
    else {point(m,0,2);auto bc=translations(true,true,true);
        component(bc,R_VARIABLE_DISPLACEMENT_X,1e-14);m.getPoint(2).addBoundaryCondition(bc);}
    bind(m);const auto fingerprint=m.getProblemSetup().getStructuralTrussDynamicsSetup().sourceFingerprint;
    QVERIFY(rejects([&]{Adapter invalid(m);}));
    QCOMPARE(m.getProblemSetup().getStructuralTrussDynamicsSetup().sourceFingerprint,fingerprint);
} REPORT_FE_EXCEPTIONS
void TestStructuralTrussDynamicsFE::tinySlackSupportRestartCorruption() try
{
    auto model=axial(3,.002,4);Adapter adapter(model);
    const auto accepted=run(adapter,2,fresh(adapter));
    const auto valid=Restart::build(adapter,accepted.state,2,accepted.cursor);
    QCOMPARE(valid.storedEnergy,0.0);QCOMPARE(valid.reactions[0],0.0);
    for(bool velocity:{false,true})
    {
        auto corrupted=valid;
        if(velocity) corrupted.v[0]=1e-14;else corrupted.q[0]=1e-14;
        QVERIFY(rejects([&]{Restart::restore(adapter,corrupted);}));
    }
    QCOMPARE(stateBytes(Restart::restore(adapter,valid).state),stateBytes(accepted.state));
} REPORT_FE_EXCEPTIONS
#undef REPORT_FE_EXCEPTIONS
QTEST_APPLESS_MAIN(TestStructuralTrussDynamicsFE)
#include "tst_structural_truss_dynamics_fe.moc"

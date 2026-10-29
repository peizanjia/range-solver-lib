#include <QtTest>
#include <rml_model.h>
#include <rml_boundary_condition.h>
#include <rml_initial_condition.h>
#include <rml_environment_condition.h>
#include "rstructural_piezo_layered_model_adapter.h"
#include <cmath>
#include <functional>
#include <algorithm>
#include <stdexcept>

namespace {
using Setup=RStructuralPiezoLayeredSetup;
using Adapter=RStructuralPiezoLayeredModelAdapter;
using Dof=RStructuralPiezoPhysicalDof;
bool fails(const std::function<void()> &f) {try {f();} catch(...) {return true;} return false;}
bool close(double a,double b) {return std::abs(a-b)<=1e-9*std::max({1.,std::abs(a),std::abs(b)});}
RModel scene(bool plate=false) {
    RModel m;
    const double xy[4][2]={{0,0},{1,0},{1,1},{0,1}};
    uint count=plate?4:3;
    for(uint i=0;i<count;++i) {m.addNode(plate?RNode(xy[i][0],xy[i][1],0):RNode(double(i),0,0)); m.setNodePhysicalId(i,101+11*i);}
    for(uint i=0;i<(plate?1u:2u);++i) {RElement e(plate?R_ELEMENT_QUAD1:R_ELEMENT_BEAM1);
        for(uint j=0;j<(plate?4u:2u);++j) e.setNodeId(j,plate?j:i+j);
        m.addElement(e,true,0); m.setElementPhysicalId(i,501+13*i);}
    // This fixture declares its complete loading through the layered setup.
    // New groups carry legacy gravity/temperature defaults; remove those
    // explicitly since removeAllEnvironmentConditions() reinstates them.
    for(uint i=0;i<(plate?m.getNSurfaces():m.getNLines());++i) {
        RElementGroup &group=plate?static_cast<RElementGroup &>(m.getSurface(i)):static_cast<RElementGroup &>(m.getLine(i));
        group.removeEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
        group.removeEnvironmentCondition(R_ENVIRONMENT_CONDITION_TEMPERATURE);
    }
    m.setProblemTaskTree(RProblemTaskItem::generateDefaultTree(R_PROBLEM_STRESS));
    m.getProblemSetup().setMassMatrixType(RMassMatrixType::Consistent);
    auto &d=m.getProblemSetup().getStructuralDynamicsSetup(); d.setEnabled(true); d.setTimeStep(.01); d.setNSteps(10);
    m.getTimeSolver().setInputStartTime(.01); return m;
}
Setup config(bool plate=false) {
    Setup s; s.enabled=true; if(plate) s.family=RStructuralPiezoLayeredModelFamily::FlatPlate;
    auto native=RStructuralPiezoContinuumMaterial::isotropic(30,.2,7);
    RStructuralPiezoLayeredMaterial host; host.materialId="host"; host.elasticityOrCompliance=native.elasticity; host.density=7;
    host.sourceId="independent-example"; host.sourceSha256=QString(64,'a');
    auto active=host; active.materialId="piezo"; active.piezoelectric=true; active.coupling[12]=2; active.coupling[13]=2;
    active.permittivity={{4,0,0,0,4,0,0,0,4}}; s.materials={host,active};
    RStructuralPiezoLayeredPort drive; drive.readout.portId="drive"; drive.initialValue=.12;
    RStructuralPiezoLayeredPort sense; sense.readout.portId="sense"; sense.readout.role=RStructuralPiezoPortRole::Sensor;
    sense.readout.readout=RStructuralPiezoReadoutMode::OpenCircuitVoltage; sense.initialValue=.03;
    s.ports={drive}; if(!plate) s.ports.push_back(sense);
    for(uint i=0;i<(plate?1u:2u);++i) {
        RStructuralPiezoLayeredModelLayer h; h.materialId="host"; h.lowerZ=.1; h.upperZ=.3;
        RStructuralPiezoLayeredModelLayer p; p.materialId="piezo"; p.lowerZ=.3; p.upperZ=.4; p.portId=i==0?"drive":"sense";
        s.members.push_back({501+13*i,{h,p}});
    }
    for(uint node:plate?std::vector<uint>{101,134}:std::vector<uint>{101}) {
        for(Dof d:plate?std::vector<Dof>{Dof::UX,Dof::UY,Dof::UZ,Dof::RX,Dof::RY}:std::vector<Dof>{Dof::UX,Dof::UZ,Dof::RY}) s.supports.push_back({node,d,0});
    }
    s.backgroundForces.push_back({plate?123u:123u,{{.2,0,-.1,0,0,0}}}); return s;
}
void explicitReference(Setup &s,const Adapter &a) {
    s.workpoint.available=true; auto q=a.mesh().reconstructPhysical(a.data().q0),f=a.mesh().reconstructPhysical(a.data().referenceForce);
    for(size_t i=0;i<a.mesh().nodeIds.size();++i) {
        RStructuralPiezoLayeredNodeWorkpoint node; node.physicalNodeId=a.mesh().nodeIds[i];
        std::copy(q.begin()+6*i,q.begin()+6*i+6,node.q0.begin()); std::copy(f.begin()+6*i,f.begin()+6*i+6,node.f0.begin()); s.workpoint.nodes.push_back(node);
    }
    for(const auto &p:a.electrical().referenceResponse().ports) s.workpoint.ports.push_back({p.portId,p.voltage,p.charge});
}
}
class TestStructuralPiezoLayeredModelAdapter:public QObject {
    Q_OBJECT
private slots:
    void multiBeamActualSourceOperators();
    void actualPlateAndPartialCoverage();
    void rotatedAffineGlobalConstraints();
    void certifiedExplicitReference();
    void unsupportedInitialFieldsAndMechanisms();
    void dampingAuthoritiesAndFullMassAcceleration();
    void sourceProvenanceAndLegacyGuards();
    void unmigratedLegacyConditions_data();
    void unmigratedLegacyConditions();
    void continuousElectrodeShortAndAmplifier();
    void strainChargeMaterialConversion();
};
void TestStructuralPiezoLayeredModelAdapter::multiBeamActualSourceOperators() {
    auto m=scene(); auto s=config(); Adapter a(m,s); const auto &d=a.data();
    QCOMPARE(a.mesh().nodeIds,(std::vector<uint>{101,112,123})); QCOMPARE(a.mesh().memberIds,(std::vector<uint>{501,514}));
    QCOMPARE(d.fullDofs,9u); QCOMPARE(d.freeDofs,6u); QCOMPARE(d.constraintRank,3u);
    QVERIFY(close(a.mesh().operators.totalMass,4.2)); QVERIFY(std::abs(a.mesh().operators.mass[2])>0);
    QVERIFY(d.forceResidual<1e-10); QVERIFY(d.initialAcceleration.size()==6);
    for(double x:d.initialAcceleration) QVERIFY(std::abs(x)<1e-9);
    for(double x:d.backgroundOffset) QVERIFY(std::abs(x)<1e-9);
    const auto &ports=a.electrical().referenceResponse().ports;
    QCOMPARE(ports[0].portId,QString("drive")); QVERIFY(close(ports[0].voltage,.12));
    QCOMPARE(ports[1].portId,QString("sense")); QVERIFY(close(ports[1].charge,.03));
    // Independently assemble the actual source connectivity with production mesh.
    std::vector<RStructuralPiezoLayeredMember> members;
    for(uint i=0;i<2;++i) {std::vector<RStructuralPiezoLayer> layers;
        for(const auto &source:s.members[i].layers) {RStructuralPiezoLayer l; l.material=RStructuralPiezoContinuumMaterial::isotropic(30,.2,7);
            if(!source.portId.isEmpty()) {l.material.piezoelectric=true; l.material.coupling[12]=2; l.material.coupling[13]=2; l.material.permittivity={{4,0,0,0,4,0,0,0,4}};}
            l.lowerZ=source.lowerZ; l.upperZ=source.upperZ; l.portId=source.portId; layers.push_back(l);}
        members.push_back({501+13*i,{101+11*i,112+11*i},layers});}
    auto native=RStructuralPiezoLayeredMesh::assemble(RStructuralPiezoLayeredFamily::Beam,{},{{101,{{0,0,0}}},{112,{{1,0,0}}},{123,{{2,0,0}}}},members);
    QCOMPARE(a.mesh().operators.mass,native.operators.mass); QCOMPARE(a.mesh().operators.stiffness,native.operators.stiffness);
    QCOMPARE(a.mesh().operators.coupling,native.operators.coupling); QCOMPARE(a.mesh().operators.capacitance,native.operators.capacitance);
}
void TestStructuralPiezoLayeredModelAdapter::actualPlateAndPartialCoverage() {
    auto m=scene(true); auto s=config(true); Adapter full(m,s); QCOMPARE(full.data().fullDofs,20u); QCOMPARE(full.data().freeDofs,10u);
    s.members[0].layers[1].coverage={{-1,0,-1,1}}; Adapter half(m,s);
    QVERIFY(close(half.mesh().operators.totalMass,1.75)); QVERIFY(close(half.mesh().operators.capacitance[0],full.mesh().operators.capacitance[0]/2));
    QCOMPARE(half.mesh().memberNodeIndices,(std::vector<std::vector<uint>>{{0,1,2,3}}));
}
void TestStructuralPiezoLayeredModelAdapter::rotatedAffineGlobalConstraints() {
    auto m=scene(); auto s=config(); const double c=std::sqrt(.5);
    s.frame.ex={{c,0,c}}; s.frame.ey={{0,1,0}}; s.frame.ez={{-c,0,c}};
    for(uint i=0;i<3;++i) m.getNode(i).set(c*i,0,c*i);
    s.supports[0].value=.01; s.supports[1].value=.02;
    s.supports.push_back({101,Dof::UY,0}); s.supports.push_back({101,Dof::RX,0});
    Adapter a(m,s); QCOMPARE(a.data().constraintRank,3u); QCOMPARE(a.data().redundantConstraints,2u);
    auto q=a.mesh().reconstructPhysical(a.data().q0); QVERIFY(close(q[0],.01)); QVERIFY(close(q[2],.02));
    auto z=a.reduceDisplacement(q),back=a.reconstruct(z); for(size_t i=0;i<q.size();++i) QVERIFY(close(q[i],back[i]));
    for(size_t i=0;i<a.data().freeDofs;++i) for(size_t j=0;j<a.data().freeDofs;++j) {
        double dot=0; for(size_t k=0;k<a.data().fullDofs;++k) dot+=a.data().T[k*a.data().freeDofs+i]*a.data().T[k*a.data().freeDofs+j]; QVERIFY(close(dot,i==j?1:0));}
    s.supports.back().value=1; QVERIFY(fails([&]{Adapter bad(m,s);}));
    // Different global rows can constrain the same supported rotation/shift.
    // Their equality is checked by affine elimination, rather than treating
    // each global axis as an independent generalized DOF.
    m=scene(); s=config(); s.frame.ex={{c,-c,0}}; s.frame.ey={{c,c,0}}; s.frame.ez={{0,0,1}};
    for(uint i=0;i<3;++i) m.getNode(i).set(c*i,-c*i,0);
    s.backgroundForces.clear(); s.supports.push_back({101,Dof::UY,0}); s.supports.push_back({101,Dof::RX,0});
    Adapter redundant(m,s); QCOMPARE(redundant.data().constraintRank,3u); QCOMPARE(redundant.data().redundantConstraints,2u);
    s.supports.back().value=1; QVERIFY(fails([&]{Adapter inconsistent(m,s);}));
}
void TestStructuralPiezoLayeredModelAdapter::certifiedExplicitReference() {
    auto m=scene(); auto s=config(); Adapter a(m,s); explicitReference(s,a); Adapter balanced(m,s);
    QCOMPARE(balanced.data().q0,a.data().q0); QVERIFY(balanced.data().electricResidual<1e-10);
    auto wrong=s; wrong.workpoint.nodes.back().f0[0]+=.1; QVERIFY(fails([&]{Adapter bad(m,wrong);}));
    wrong=s; wrong.workpoint.ports.back().Q0+=.1; QVERIFY(fails([&]{Adapter bad(m,wrong);}));
    wrong=s; wrong.workpoint.nodes.front().q0[0]=.1; QVERIFY(fails([&]{Adapter bad(m,wrong);}));
}
void TestStructuralPiezoLayeredModelAdapter::unsupportedInitialFieldsAndMechanisms() {
    auto m=scene(); auto s=config();
    auto wrong=s; wrong.backgroundForces.back().values[1]=.01; QVERIFY(fails([&]{Adapter bad(m,wrong);}));
    wrong=s; wrong.initialDisplacements={{123,{{0,0,0,1,0,0}}}}; QVERIFY(fails([&]{Adapter bad(m,wrong);}));
    wrong=s; wrong.initialVelocities={{101,{{1,0,0,0,0,0}}}}; QVERIFY(fails([&]{Adapter bad(m,wrong);}));
    wrong=s; wrong.supports.clear(); QVERIFY(fails([&]{Adapter bad(m,wrong);}));
}
void TestStructuralPiezoLayeredModelAdapter::dampingAuthoritiesAndFullMassAcceleration() {
    auto m=scene(); auto s=config(); s.damping.kind=RStructuralPiezoDampingKind::Rayleigh; s.damping.rayleighAlpha=.1; s.damping.rayleighBeta=.02;
    s.initialDisplacements={{123,{{.01,0,-.02,0,.03,0}}}}; Adapter a(m,s); const auto &d=a.data();
    std::vector<double> q=a.mesh().projectPhysical(a.reconstruct(d.initialZ),RStructuralPiezoPhysicalDomain::Displacement);
    auto force=a.electrical().evaluate(q,d.actuatorVoltages0).internalForce; std::vector<double> residual(d.freeDofs);
    for(size_t i=0;i<d.freeDofs;++i) {for(size_t j=0;j<d.freeDofs;++j) residual[i]+=d.mass[i*d.freeDofs+j]*d.initialAcceleration[j];
        for(size_t k=0;k<d.fullDofs;++k) residual[i]+=d.T[k*d.freeDofs+i]*(force[k]-d.backgroundForce[k]); QVERIFY(std::abs(residual[i])<1e-9);}
    auto imported=s; imported.damping={}; imported.damping.kind=RStructuralPiezoDampingKind::Imported;
    imported.damping.importedDimension=d.fullDofs; imported.damping.importedMatrix=d.fullDamping; imported.damping.sourceId="test-C"; imported.damping.sourceSha256=QString(64,'b');
    Adapter same(m,imported); QCOMPARE(same.data().fullDamping,d.fullDamping); QCOMPARE(same.data().damping,d.damping);
    imported.damping.importedMatrix[0]=-1; QVERIFY(fails([&]{Adapter bad(m,imported);}));
    auto modal=s; modal.damping={}; modal.damping.kind=RStructuralPiezoDampingKind::UniformModal; modal.damping.uniformRatio=.02;
    QVERIFY(fails([&]{Adapter missing(m,modal);})); Adapter supplied(m,modal,d.fullDamping); QCOMPARE(supplied.data().damping,d.damping);
}
void TestStructuralPiezoLayeredModelAdapter::sourceProvenanceAndLegacyGuards() {
    auto m=scene(); auto s=config(); s.bindToModel(m); Adapter a(m,s); QCOMPARE(a.data().materialProvenance[0].sourceId,QString("independent-example"));
    m.getLine(0).setCrossArea(100); QVERIFY(fails([&]{Adapter stale(m,s);})); s.sourceFingerprint.clear(); Adapter changed(m,s);
    QCOMPARE(changed.mesh().operators.mass,a.mesh().operators.mass); QVERIFY(changed.data().sourceFingerprint!=a.data().sourceFingerprint);
    // Existing host section/material changes invalidate identity but never add
    // another old host contribution to the actual layered production matrices.
    auto wrong=s; wrong.family=RStructuralPiezoLayeredModelFamily::FlatPlate; QVERIFY(fails([&]{Adapter bad(m,wrong);}));
    m.getLine(0).addBoundaryCondition(RBoundaryCondition(R_BOUNDARY_CONDITION_DISPLACEMENT)); QVERIFY(fails([&]{Adapter ignored(m,s);}));
}
void TestStructuralPiezoLayeredModelAdapter::unmigratedLegacyConditions_data() {
    QTest::addColumn<bool>("plate"); QTest::addColumn<int>("condition");
    const char *names[]={"defaults","boundary","initial","gravity","temperature"};
    for(bool plate:{false,true}) for(int condition=0;condition<5;++condition) {
        const QByteArray name=QByteArray(plate?"plate":"beam")+"-"+names[condition];
        QTest::newRow(name.constData())<<plate<<condition;
    }
}
void TestStructuralPiezoLayeredModelAdapter::unmigratedLegacyConditions() {
    QFETCH(bool,plate); QFETCH(int,condition);
    auto m=scene(plate); auto s=config(plate);
    Adapter admitted(m,s);
    RElementGroup &group=plate?static_cast<RElementGroup &>(m.getSurface(0)):static_cast<RElementGroup &>(m.getLine(0));
    QCOMPARE(group.getNBoundaryConditions(),0u); QCOMPARE(group.getNInitialConditions(),0u); QCOMPARE(group.getNEnvironmentConditions(),0u);
    switch(condition) {
    case 0: group.removeAllEnvironmentConditions(); QCOMPARE(group.getNEnvironmentConditions(),2u); break;
    case 1: group.addBoundaryCondition(RBoundaryCondition(R_BOUNDARY_CONDITION_DISPLACEMENT)); break;
    case 2: group.addInitialCondition(RInitialCondition(R_INITIAL_CONDITION_VELOCITY)); break;
    case 3: group.addEnvironmentCondition(REnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION)); break;
    case 4: group.addEnvironmentCondition(REnvironmentCondition(R_ENVIRONMENT_CONDITION_TEMPERATURE)); break;
    }
    // Rejection must come from the authoritative-condition guard, rather than
    // a stale fingerprint or some later numerical validation.
    QVERIFY(s.sourceFingerprint.isEmpty());
    try {Adapter rejected(m,s); QFAIL("Unmigrated legacy condition was accepted");}
    catch(const std::invalid_argument &error) {
        QCOMPARE(QString::fromUtf8(error.what()),QString("Migrate legacy group conditions to the authoritative layered setup before compilation"));
    }
}
void TestStructuralPiezoLayeredModelAdapter::continuousElectrodeShortAndAmplifier() {
    auto m=scene(); auto s=config(); s.members[1].layers[1].portId="drive"; s.ports.resize(1);
    Adapter continuous(m,s); QCOMPARE(continuous.mesh().operators.portIds,(std::vector<QString>{"drive"}));
    auto separate=config(); Adapter two(m,separate);
    QVERIFY(close(continuous.mesh().operators.capacitance[0],two.mesh().operators.capacitance[0]+two.mesh().operators.capacitance[3]));
    for(auto mode:{RStructuralPiezoReadoutMode::ShortCircuitCharge,RStructuralPiezoReadoutMode::IdealChargeAmplifier}) {
        auto source=config(); source.ports[1].readout.readout=mode; source.ports[1].initialValue=0;
        if(mode==RStructuralPiezoReadoutMode::IdealChargeAmplifier) {
            source.ports[1].readout.Cf=2; source.ports[1].readout.Qref=.2; source.ports[1].readout.amplifierPolarity=1; source.ports[1].readout.VoutBias=.7;
        }
        Adapter a(m,source); const auto &response=a.electrical().referenceResponse().ports[1]; QVERIFY(close(response.voltage,0));
        QCOMPARE(a.electrical().condensed().stiffness,a.mesh().operators.stiffness);
        if(mode==RStructuralPiezoReadoutMode::IdealChargeAmplifier) {
            QCOMPARE(a.data().portMetadata[1].readout.amplifierPolarity,1); QVERIFY(close(a.data().portMetadata[1].readout.VoutBias,.7));
            QVERIFY(close(response.amplifierVoltage,-(response.charge-.2)/2));
        }
    }
}
void TestStructuralPiezoLayeredModelAdapter::strainChargeMaterialConversion() {
    auto m=scene(); auto s=config(); auto &active=s.materials[1];
    active.form=RStructuralPiezoMaterialForm::StrainCharge; active.elasticityOrCompliance={}; active.coupling={};
    for(size_t i=0;i<6;++i) active.elasticityOrCompliance[6*i+i]=.02;
    active.coupling[12]=.01; active.permittivity={{4,0,0,0,4,0,0,0,4}};
    Adapter converted(m,s);
    auto native=RStructuralPiezoContinuumMaterial::fromStrainCharge(active.elasticityOrCompliance,active.coupling,active.permittivity,active.density);
    auto equivalent=s; equivalent.materials[1].form=RStructuralPiezoMaterialForm::StressCharge;
    equivalent.materials[1].elasticityOrCompliance=native.elasticity; equivalent.materials[1].coupling=native.coupling; equivalent.materials[1].permittivity=native.permittivity;
    Adapter direct(m,equivalent); QCOMPARE(converted.mesh().operators.stiffness,direct.mesh().operators.stiffness);
    QCOMPARE(converted.mesh().operators.coupling,direct.mesh().operators.coupling); QCOMPARE(converted.mesh().operators.capacitance,direct.mesh().operators.capacitance);
}
QTEST_APPLESS_MAIN(TestStructuralPiezoLayeredModelAdapter)
#include "tst_structural_piezo_layered_model_adapter.moc"

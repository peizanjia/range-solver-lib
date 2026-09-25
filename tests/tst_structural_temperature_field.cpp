#include <QtTest>

#include <limits>

#include <rml_model.h>

#include "rstructural_temperature_field.h"
#include "rstructural_coupled_field_snapshot.h"

namespace
{
RModel mesh(const std::vector<QString> &nodeIds,bool reverseElementOrder = false)
{
    RModel model;
    model.setNNodes(uint(nodeIds.size()));
    for (uint index=0;index<nodeIds.size();index++)
    {
        const QString &id = nodeIds[index];
        const double x = id == QStringLiteral("A") ? 0.0
                       : id == QStringLiteral("B") ? 1.0 : 2.0;
        model.getNode(index).set(x,0.0,0.0);
    }
    const auto indexOf = [&](const QString &id) {
        for (uint index=0;index<nodeIds.size();index++)
            if (nodeIds[index] == id) return index;
        return uint(nodeIds.size());
    };
    auto addRod = [&](const QString &left,const QString &right) {
        RElement element(R_ELEMENT_TRUSS1);
        element.setNodeId(0,indexOf(left));
        element.setNodeId(1,indexOf(right));
        model.addElement(element);
    };
    if (reverseElementOrder)
    {
        addRod(QStringLiteral("B"),QStringLiteral("C"));
        addRod(QStringLiteral("A"),QStringLiteral("B"));
    }
    else
    {
        addRod(QStringLiteral("A"),QStringLiteral("B"));
        addRod(QStringLiteral("B"),QStringLiteral("C"));
    }
    return model;
}
}

class TestStructuralTemperatureField : public QObject
{
    Q_OBJECT
private slots:
    void reordersByPhysicalIdWithoutChangingTheField();
    void rejectsGeometricChangeBelowLegacyFiveDigitPrecision();
    void rejectsInvalidIdsAndTemperatures();
    void rejectsOrientedConnectivityChange();
    void coupledSnapshotKeepsPhysicalIdsAndMissingFieldsExplicit();
    void coupledSnapshotUsesPersistedModelIds();
    void coupledSnapshotReadsOptionalMotionAndThermalFields();
};

void TestStructuralTemperatureField::coupledSnapshotReadsOptionalMotionAndThermalFields()
{
    RModel model = mesh({"A","B","C"});
    RVariable temperature(R_VARIABLE_TEMPERATURE,R_VARIABLE_APPLY_NODE);
    temperature.resize(1,3);
    for (uint i=0;i<3;i++) temperature.setValue(0,i,293.15+double(i));
    model.addVariable(temperature);
    RVariable velocity(R_VARIABLE_VELOCITY,R_VARIABLE_APPLY_NODE);
    velocity.resize(3,3);
    RVariable acceleration(R_VARIABLE_ACCELERATION,R_VARIABLE_APPLY_NODE);
    acceleration.resize(3,3);
    for (uint i=0;i<3;i++)
    {
        for (uint component=0;component<3;component++)
        {
            velocity.setValue(component,i,double(10*i+component));
            acceleration.setValue(component,i,double(100*i+component));
        }
    }
    model.addVariable(velocity);
    model.addVariable(acceleration);
    RVariable heatFlux(R_VARIABLE_HEAT_FLUX,R_VARIABLE_APPLY_ELEMENT);
    heatFlux.resize(3,2);
    RVariable heatSource(R_VARIABLE_HEAT,R_VARIABLE_APPLY_ELEMENT);
    heatSource.resize(1,2);
    for (uint i=0;i<2;i++)
    {
        for (uint component=0;component<3;component++)
            heatFlux.setValue(component,i,double(10*i+component+1));
        heatSource.setValue(0,i,double(i+5));
    }
    model.addVariable(heatFlux);
    model.addVariable(heatSource);

    const auto snapshot = RStructuralCoupledFieldSnapshot::capture(model);
    const QString nodeId = QString::number(model.getNodePhysicalId(1));
    const auto node = snapshot.selectNodes({nodeId}).front();
    QCOMPARE(node.velocityMetresPerSecond.value(),(std::array<double,3>{10.0,11.0,12.0}));
    QCOMPARE(node.accelerationMetresPerSecondSquared.value(),
             (std::array<double,3>{100.0,101.0,102.0}));
    QVERIFY(!node.displacementMetres.has_value());
    const QString elementId = QString::number(model.getElementPhysicalId(1));
    const auto element = snapshot.selectElements({elementId}).front();
    QCOMPARE(element.heatFluxWattsPerSquareMetre,std::vector<double>({11.0,12.0,13.0}));
    QCOMPARE(element.heatSourceWatts,std::vector<double>({6.0}));
    QVERIFY(element.stressPascals.empty());
}

void TestStructuralTemperatureField::coupledSnapshotUsesPersistedModelIds()
{
    RModel model = mesh({"A","B","C"});
    model.setNodePhysicalId(0,10);
    model.setNodePhysicalId(1,20);
    model.setNodePhysicalId(2,30);
    model.setElementPhysicalId(0,40);
    model.setElementPhysicalId(1,50);
    RVariable temperature(R_VARIABLE_TEMPERATURE,R_VARIABLE_APPLY_NODE);
    temperature.resize(1,3);
    for (uint i=0;i<3;i++) temperature.setValue(0,i,300.0+double(i));
    model.addVariable(temperature);
    const auto snapshot = RStructuralCoupledFieldSnapshot::capture(model);
    QCOMPARE(snapshot.selectNodes({"20"})[0].temperatureKelvin,301.0);
    QCOMPARE(snapshot.selectElements({"50"}).size(),std::size_t(1));
    snapshot.verifyMesh(model);
    RModel changed(model);
    changed.setNodePhysicalId(1,21);
    QVERIFY_EXCEPTION_THROWN(snapshot.verifyMesh(changed),std::invalid_argument);
}

void TestStructuralTemperatureField::coupledSnapshotKeepsPhysicalIdsAndMissingFieldsExplicit()
{
    const std::vector<QString> ids = {"A","B","C"};
    RModel first = mesh(ids);
    RVariable temperature(R_VARIABLE_TEMPERATURE,R_VARIABLE_APPLY_NODE);
    temperature.resize(1,3);
    for (uint i=0;i<3;i++) temperature.setValue(0,i,300.0+10.0*i);
    first.addVariable(temperature);
    RVariable displacement(R_VARIABLE_DISPLACEMENT,R_VARIABLE_APPLY_NODE);
    displacement.resize(3,3);
    for (uint i=0;i<3;i++)
    {
        displacement.setValue(0,i,0.001*i);
        displacement.setValue(1,i,0.0);
        displacement.setValue(2,i,0.0);
    }
    first.addVariable(displacement);
    const RVariableType stressTypes[6] = {
        R_VARIABLE_STRESS_X,R_VARIABLE_STRESS_Y,R_VARIABLE_STRESS_Z,
        R_VARIABLE_STRESS_YZ,R_VARIABLE_STRESS_XZ,R_VARIABLE_STRESS_XY
    };
    for (uint component=0;component<6;component++)
    {
        RVariable stress(stressTypes[component],R_VARIABLE_APPLY_ELEMENT);
        stress.resize(1,2);
        stress.setValue(0,0,11.0+component);
        stress.setValue(0,1,22.0+component);
        first.addVariable(stress);
    }
    first.getTimeSolver().setComputedTime(1.0);
    const auto before = RStructuralCoupledFieldSnapshot::capture(
        first,ids,{"AB","BC"});
    QCOMPARE(before.selectNodes({"B"})[0].temperatureKelvin,310.0);
    QVERIFY(!before.selectNodes({"B"})[0].velocityMetresPerSecond.has_value());
    QVERIFY(!before.selectNodes({"B"})[0].temperatureRateKelvinPerSecond.has_value());
    QCOMPARE(before.selectElements({"BC"})[0].stressPascals,
             std::vector<double>({22.0,23.0,24.0,25.0,26.0,27.0}));

    const std::vector<QString> reorderedIds = {"C","A","B"};
    RModel reordered = mesh(reorderedIds,true);
    RVariable nextTemperature(R_VARIABLE_TEMPERATURE,R_VARIABLE_APPLY_NODE);
    nextTemperature.resize(1,3);
    nextTemperature.setValue(0,0,322.0);
    nextTemperature.setValue(0,1,302.0);
    nextTemperature.setValue(0,2,312.0);
    reordered.addVariable(nextTemperature);
    reordered.getTimeSolver().setComputedTime(2.0);
    const auto after = RStructuralCoupledFieldSnapshot::capture(
        reordered,reorderedIds,{"BC","AB"},&before);
    before.verifyMesh(reordered,reorderedIds,{"BC","AB"});
    QCOMPARE(after.meshDigest(),before.meshDigest());
    QCOMPARE(after.selectNodes({"B"})[0].temperatureRateKelvinPerSecond.value(),2.0);
    QVERIFY(!after.selectNodes({"B"})[0].displacementMetres.has_value());
    QVERIFY_EXCEPTION_THROWN(after.selectNodes({"B","B"}),std::invalid_argument);
    reordered.getNode(2).set(1.0+1.0e-9,0.0,0.0);
    QVERIFY_EXCEPTION_THROWN(after.verifyMesh(reordered,reorderedIds,{"BC","AB"}),
                             std::invalid_argument);
}

void TestStructuralTemperatureField::reordersByPhysicalIdWithoutChangingTheField()
{
    const std::vector<QString> sourceIds = {"A","B","C"};
    const auto source = mesh(sourceIds);
    const auto field = RStructuralTemperatureField::capture(source,sourceIds,{300.0,350.0,400.0});
    const std::vector<QString> reorderedIds = {"C","A","B"};
    const auto reordered = mesh(reorderedIds,true);
    QCOMPARE(field.mapTo(reordered,reorderedIds),std::vector<double>({400.0,300.0,350.0}));
    QCOMPARE(field.meshDigest().size(),64);
}

void TestStructuralTemperatureField::rejectsGeometricChangeBelowLegacyFiveDigitPrecision()
{
    const std::vector<QString> ids = {"A","B","C"};
    const auto field = RStructuralTemperatureField::capture(mesh(ids),ids,{300.0,350.0,400.0});
    auto changed = mesh(ids);
    changed.getNode(1).set(1.0+1.0e-9,0.0,0.0);
    QVERIFY_EXCEPTION_THROWN(field.mapTo(changed,ids),std::invalid_argument);
}

void TestStructuralTemperatureField::rejectsInvalidIdsAndTemperatures()
{
    const std::vector<QString> ids = {"A","B","C"};
    const auto source = mesh(ids);
    QVERIFY_EXCEPTION_THROWN(RStructuralTemperatureField::capture(
        source,{"A","B","B"},{300.0,350.0,400.0}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralTemperatureField::capture(
        source,ids,{300.0,std::numeric_limits<double>::quiet_NaN(),400.0}),
        std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralTemperatureField::capture(
        source,ids,{300.0,-1.0,400.0}),std::invalid_argument);
}

void TestStructuralTemperatureField::rejectsOrientedConnectivityChange()
{
    const std::vector<QString> ids = {"A","B","C"};
    const auto field = RStructuralTemperatureField::capture(mesh(ids),ids,{300.0,350.0,400.0});
    RModel reversed;
    reversed.setNNodes(3);
    for (uint i=0;i<3;i++) reversed.getNode(i).set(double(i),0.0,0.0);
    RElement first(R_ELEMENT_TRUSS1);
    first.setNodeId(0,1);
    first.setNodeId(1,0);
    reversed.addElement(first);
    RElement second(R_ELEMENT_TRUSS1);
    second.setNodeId(0,1);
    second.setNodeId(1,2);
    reversed.addElement(second);
    QVERIFY_EXCEPTION_THROWN(field.mapTo(reversed,ids),std::invalid_argument);
}

QTEST_APPLESS_MAIN(TestStructuralTemperatureField)
#include "tst_structural_temperature_field.moc"

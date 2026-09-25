#include <QtTest>

#include <limits>

#include <rml_model.h>

#include "rstructural_temperature_field.h"

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
};

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

#include <QtTest>

#include <cmath>
#include <limits>
#include <stdexcept>

#include "rstructural_control_force_map.h"

class TestStructuralControlForceMap : public QObject
{
    Q_OBJECT

private slots:
    void sparseMappingAndPhysicalPower();
    void rejectInvalidCommandsAndMap();
    void heldCommandLedgerTrialRejectAndRestore();
};

void TestStructuralControlForceMap::sparseMappingAndPhysicalPower()
{
    const RStructuralControlForceMap map(3,2.0,4.0,{
        {"tip-x",{{0,1.0},{2,-0.5}}},
        {"tip-y",{{1,1.0}}}
    });
    const std::vector<RStructuralControlCommand> commands = {{"tip-x",3.0},{"tip-y",-2.0}};
    RRVector reduced(3,0.0);
    map.addForce(commands,reduced);
    QCOMPARE(reduced[0],6.0);
    QCOMPARE(reduced[1],-4.0);
    QCOMPARE(reduced[2],-3.0);
    RRVector velocity(3,0.0);
    velocity[0] = 4.0; velocity[1] = 8.0; velocity[2] = 12.0;
    const auto powers = map.powers(commands,velocity);
    QCOMPARE(powers.size(),size_t(2));
    QCOMPARE(powers[0],-1.5); // 3*(4-0.5*12)/4 physical watts.
    QCOMPARE(powers[1],-4.0);
    // Sparse weighted force and power satisfy physical duality under scaling.
    const double scaledRate = RRVector::dot(reduced,velocity);
    QVERIFY(std::abs(scaledRate/(2.0*4.0)-(powers[0]+powers[1])) < 1.0e-12);
}

void TestStructuralControlForceMap::rejectInvalidCommandsAndMap()
{
    QVERIFY_EXCEPTION_THROWN(RStructuralControlForceMap(2,1.0,1.0,
        {{"same",{{0,1.0}}},{"same",{{1,1.0}}}}),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralControlForceMap(2,1.0,1.0,
        {{"bad",{{2,1.0}}}}),std::invalid_argument);
    const RStructuralControlForceMap map(2,1.0,1.0,{{"a",{{0,1.0}}}});
    RRVector force(2,0.0);
    QVERIFY_EXCEPTION_THROWN(map.addForce({{"b",1.0}},force),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(map.addForce({{"a",std::numeric_limits<double>::infinity()}},force),
                             std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(map.powers({{"a",1.0}},RRVector(1,0.0)),std::invalid_argument);
    QCOMPARE(force[0],0.0);
}

void TestStructuralControlForceMap::heldCommandLedgerTrialRejectAndRestore()
{
    const RStructuralControlForceMap map(1,2.0,4.0,{{"force-x",{{0,1.0}}}});
    RStructuralControlLedger ledger(map);
    QCOMPARE(ledger.accepted()[0].cumulativeWork,0.0);
    const auto candidate = ledger.trial(map,{{"force-x",3.0}},RRVector(1,4.0),RRVector(1,8.0),0.5);
    QCOMPARE(candidate[0].endpointPower,6.0);
    QCOMPARE(candidate[0].cumulativeWork,2.25); // 0.5*(3+6)*0.5 J.
    QCOMPARE(ledger.accepted()[0].cumulativeWork,0.0);
    ledger.reject();
    QVERIFY_EXCEPTION_THROWN(ledger.candidate(),std::logic_error);
    QCOMPARE(ledger.accepted()[0].cumulativeWork,0.0);
    ledger.trial(map,{{"force-x",3.0}},RRVector(1,4.0),RRVector(1,8.0),0.5);
    ledger.accept();
    QCOMPARE(ledger.accepted()[0].cumulativeWork,2.25);
    const auto checkpoint = ledger.accepted();
    ledger.trial(map,{{"force-x",-1.0}},RRVector(1,8.0),RRVector(1,0.0),0.5);
    ledger.accept();
    QVERIFY(ledger.accepted()[0].cumulativeWork < checkpoint[0].cumulativeWork);
    ledger.restore(checkpoint);
    QCOMPARE(ledger.accepted()[0].command,3.0);
    QCOMPARE(ledger.accepted()[0].cumulativeWork,2.25);
}

QTEST_APPLESS_MAIN(TestStructuralControlForceMap)

#include "tst_structural_control_force_map.moc"

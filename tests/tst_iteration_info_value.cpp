#include <QtTest>

#include "riterationinfo.h"
#include "riterationinfovalue.h"

class TestIterationInfoValue : public QObject
{
    Q_OBJECT

private slots:

    void valueConstructor();
    void copyAndAssignment();
    void stagnationIsNotConvergence();
};

void TestIterationInfoValue::valueConstructor()
{
    RIterationInfoValue info("residual", 1.25);
    QCOMPARE(info.getName(), QString("residual"));
    QCOMPARE(info.getValue(), 1.25);
}

void TestIterationInfoValue::copyAndAssignment()
{
    RIterationInfoValue info("error", 0.5);

    RIterationInfoValue copy(info);
    QCOMPARE(copy.getName(), QString("error"));
    QCOMPARE(copy.getValue(), 0.5);

    RIterationInfoValue assigned;
    assigned = info;
    QCOMPARE(assigned.getName(), QString("error"));
    QCOMPARE(assigned.getValue(), 0.5);
}

void TestIterationInfoValue::stagnationIsNotConvergence()
{
    RIterationInfo info;
    info.setConvergenceValue(1.0e-8);
    info.setError(1.0);
    info.setError(1.0);
    QVERIFY(!info.hasConverged());

    info.setError(1.0e-9);
    QVERIFY(info.hasConverged());
}

QTEST_APPLESS_MAIN(TestIterationInfoValue)

#include "tst_iteration_info_value.moc"

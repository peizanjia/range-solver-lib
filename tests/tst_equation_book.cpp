#include <QtTest>

#include <rbl_utils.h>

#include "requationbook.h"

class TestEquationBook : public QObject
{
    Q_OBJECT

    private slots:

        void numbersMixedNodalDofsDeterministically();
        void disablingDofConsolidatesEquations();
        void resizeClearsPreviousState();
};

void TestEquationBook::numbersMixedNodalDofsDeterministically()
{
    REquationBook book(3);
    for (RDofType dof : {RDofType::UX,RDofType::UY,RDofType::UZ})
    {
        book.enable(0,dof);
        book.enable(2,dof);
    }
    for (RDofType dof : {RDofType::UX,RDofType::UY,RDofType::UZ,
                         RDofType::RX,RDofType::RY,RDofType::RZ})
    {
        book.enable(1,dof);
    }

    QCOMPARE(book.getNEnabled(),12u);
    uint equation = 0;
    QVERIFY(book.getValue(0,RDofType::UX,equation));
    QCOMPARE(equation,0u);
    QVERIFY(book.getValue(0,RDofType::UZ,equation));
    QCOMPARE(equation,2u);
    QVERIFY(book.getValue(1,RDofType::UX,equation));
    QCOMPARE(equation,3u);
    QVERIFY(book.getValue(1,RDofType::RZ,equation));
    QCOMPARE(equation,8u);
    QVERIFY(book.getValue(2,RDofType::UX,equation));
    QCOMPARE(equation,9u);
    QVERIFY(!book.getValue(2,RDofType::RX,equation));
    QCOMPARE(equation,RConstants::eod);

    QCOMPARE(book.getEnabledDofs(1),
             QList<RDofType>({RDofType::UX,RDofType::UY,RDofType::UZ,
                              RDofType::RX,RDofType::RY,RDofType::RZ}));
}

void TestEquationBook::disablingDofConsolidatesEquations()
{
    REquationBook book(2);
    for (uint nodeId=0;nodeId<2;nodeId++)
    {
        book.enable(nodeId,RDofType::UX);
        book.enable(nodeId,RDofType::UY);
    }
    book.disable(0,RDofType::UY);

    QCOMPARE(book.getNEnabled(),3u);
    uint equation = 0;
    QVERIFY(book.getValue(1,RDofType::UX,equation));
    QCOMPARE(equation,1u);
    QVERIFY(book.getValue(1,RDofType::UY,equation));
    QCOMPARE(equation,2u);
}

void TestEquationBook::resizeClearsPreviousState()
{
    REquationBook book(1);
    book.enable(0,RDofType::UX);
    QCOMPARE(book.getNEnabled(),1u);

    book.resize(2);
    QCOMPARE(book.getNNodes(),2u);
    QCOMPARE(book.getNEnabled(),0u);
    QVERIFY(!book.isEnabled(0,RDofType::UX));
}

QTEST_APPLESS_MAIN(TestEquationBook)

#include "tst_equation_book.moc"

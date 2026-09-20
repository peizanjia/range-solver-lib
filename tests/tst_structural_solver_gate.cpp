#include <QtTest>

#include <rbl_error.h>
#include <rml_model.h>

#include "rsolver.h"

class TestStructuralSolverGate : public QObject
{
    Q_OBJECT

private slots:

    void rejectsEmptyProblemTask();
    void defaultPolicyPreservesRegisteredNonStructuralSolver();
    void rejectsNonStructuralProblemBeforeMutatingResults();
    void rejectsMixedStructuralAndHeatProblem();
    void rejectsElementWithoutStructuralFormulation();
};

void TestStructuralSolverGate::rejectsEmptyProblemTask()
{
    RModel model;

    try
    {
        RSolver solver(model,QString(),QString());
        Q_UNUSED(solver);
        QFAIL("An empty problem task must be rejected.");
    }
    catch (const RError &error)
    {
        QCOMPARE(error.getType(), RError::Type::InvalidInput);
        QVERIFY(error.getMessage().contains("No solver problem task"));
    }
}

void TestStructuralSolverGate::defaultPolicyPreservesRegisteredNonStructuralSolver()
{
    RModel model;
    model.setProblemTaskTree(RProblemTaskItem::generateDefaultTree(R_PROBLEM_HEAT));

    RSolver solver(model,QString(),QString());
    Q_UNUSED(solver);
}

void TestStructuralSolverGate::rejectsNonStructuralProblemBeforeMutatingResults()
{
    RModel model;
    model.setProblemTaskTree(RProblemTaskItem::generateDefaultTree(R_PROBLEM_HEAT));
    model.addVariable(RVariable(R_VARIABLE_TEMPERATURE));
    QCOMPARE(model.getNVariables(), 1u);

    try
    {
        RSolver solver(model,QString(),QString(),RSolver::ExecutionPolicy::StructuralProduct);
        Q_UNUSED(solver);
        QFAIL("A heat-transfer problem must be rejected by the structural solver.");
    }
    catch (const RError &error)
    {
        QCOMPARE(error.getType(), RError::Type::InvalidInput);
        QVERIFY(error.getMessage().contains(RProblem::getName(R_PROBLEM_HEAT)));
    }

    QCOMPARE(model.getNVariables(), 1u);
}

void TestStructuralSolverGate::rejectsMixedStructuralAndHeatProblem()
{
    RModel model;
    model.setProblemTaskTree(
        RProblemTaskItem::generateDefaultTree(R_PROBLEM_STRESS | R_PROBLEM_HEAT));

    try
    {
        RSolver solver(model,QString(),QString(),RSolver::ExecutionPolicy::StructuralProduct);
        Q_UNUSED(solver);
        QFAIL("A mixed structural and heat-transfer task must be rejected as a whole.");
    }
    catch (const RError &error)
    {
        QCOMPARE(error.getType(), RError::Type::InvalidInput);
        QVERIFY(error.getMessage().contains(RProblem::getName(R_PROBLEM_HEAT)));
    }
}

void TestStructuralSolverGate::rejectsElementWithoutStructuralFormulation()
{
    RModel model;
    model.setProblemTaskTree(RProblemTaskItem::generateDefaultTree(R_PROBLEM_STRESS));
    model.addElement(RElement(R_ELEMENT_BEAM1));

    try
    {
        RSolver solver(model,QString(),QString(),RSolver::ExecutionPolicy::StructuralProduct);
        Q_UNUSED(solver);
        QFAIL("An element without a structural formulation must be rejected.");
    }
    catch (const RError &error)
    {
        QCOMPARE(error.getType(), RError::Type::InvalidInput);
        QVERIFY(error.getMessage().contains(RElement::getId(R_ELEMENT_BEAM1)));
        QVERIFY(error.getMessage().contains("no enabled formulation"));
    }
}

QTEST_APPLESS_MAIN(TestStructuralSolverGate)

#include "tst_structural_solver_gate.moc"

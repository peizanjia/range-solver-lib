#include <QtTest>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QTemporaryDir>

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <limits>

#include <omp.h>

#include "rml_model.h"
#include "rml_file_utils.h"
#include "rsolver.h"

namespace
{
constexpr double length = 2.0;
constexpr double initialDisplacement = 1.0e-4;

RModel makeTruss()
{
    RModel model;
    model.setNNodes(2);
    model.getNode(0).set(0.0,0.0,0.0);
    model.getNode(1).set(length,0.0,0.0);
    RElement bar(R_ELEMENT_TRUSS1);
    bar.setNodeId(0,0);
    bar.setNodeId(1,1);
    model.addElement(bar,true,0);

    RMaterial material(RMaterial::Solid);
    RMaterialProperty modulus(RMaterialProperty::ModulusOfElasiticity);
    modulus.add(293.15,2.0e11);
    material.add(modulus);
    RMaterialProperty poisson(RMaterialProperty::PoissonRatio);
    poisson.add(293.15,0.0);
    material.add(poisson);
    RMaterialProperty density(RMaterialProperty::Density);
    density.add(293.15,8000.0);
    material.add(density);
    RMaterialProperty expansion(RMaterialProperty::ThermalExpansionCoefficient);
    expansion.add(293.15,0.0);
    material.add(expansion);
    RLine &line = model.getLine(0);
    line.setCrossArea(0.01);
    line.setMaterial(material);

    RBoundaryCondition transverse(R_BOUNDARY_CONDITION_DISPLACEMENT);
    transverse.getComponent(transverse.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).setEnabled(false);
    transverse.getComponent(transverse.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).add(0.0,0.0);
    transverse.getComponent(transverse.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).add(0.0,0.0);
    line.addBoundaryCondition(transverse);

    RElement root(R_ELEMENT_POINT);
    root.setNodeId(0,0);
    model.addElement(root,true,0);
    RBoundaryCondition support(R_BOUNDARY_CONDITION_DISPLACEMENT);
    support.getComponent(support.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).add(0.0,0.0);
    support.getComponent(support.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).add(0.0,0.0);
    support.getComponent(support.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).add(0.0,0.0);
    model.getPoint(0).addBoundaryCondition(support);

    RElement tip(R_ELEMENT_POINT);
    tip.setNodeId(0,1);
    model.addElement(tip,true,1);
    RInitialCondition initial(R_INITIAL_CONDITION_DISPLACEMENT);
    initial.getComponent(initial.findComponentPosition(R_VARIABLE_DISPLACEMENT_X))
           .add(0.0,initialDisplacement);
    model.getPoint(1).addInitialCondition(initial);

    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    auto &linear = model.getMatrixSolverConf(RMatrixSolverConf::CG);
    linear.setNInnerIterations(100);
    linear.setNOuterIterations(1000);
    linear.setSolverCvgValue(1.0e-13);
    linear.setOutputFrequency(0);
    auto &dynamics = model.getProblemSetup().getStructuralDynamicsSetup();
    dynamics.setEnabled(true);
    dynamics.setMethod(RStructuralTimeIntegrationMethod::Newmark);
    dynamics.setTimeStep(0.025/std::sqrt(3.0*2.0e11/(8000.0*length*length)));
    dynamics.setNSteps(1);
    dynamics.setOutputFrequency(1);
    return model;
}

bool runSolver(RModel &model, QString &error, const QString &target = QString())
{
    try
    {
        RSolver solver(model,target,QString());
        solver.run();
        return true;
    }
    catch (const RError &failure)
    {
        error = failure.getMessage();
    }
    catch (const std::exception &failure)
    {
        error = QString::fromUtf8(failure.what());
    }
    return false;
}

bool runSolver(RSolver &solver, QString &error)
{
    error.clear();
    try
    {
        solver.run();
        return true;
    }
    catch (const RError &failure)
    {
        error = failure.getMessage();
    }
    catch (const std::exception &failure)
    {
        error = QString::fromUtf8(failure.what());
    }
    return false;
}

double tipDisplacement(const RModel &model)
{
    const uint position = model.findVariable(R_VARIABLE_DISPLACEMENT);
    return position == RConstants::eod ? std::numeric_limits<double>::quiet_NaN() :
                                         model.getVariable(position).getValue(0,1);
}

std::array<double,3> tipKinematics(const RModel &model)
{
    std::array<double,3> values{};
    uint index = 0;
    for (RVariableType type : {R_VARIABLE_DISPLACEMENT,R_VARIABLE_VELOCITY,R_VARIABLE_ACCELERATION})
    {
        const uint position = model.findVariable(type);
        values[index++] = position == RConstants::eod ? std::numeric_limits<double>::quiet_NaN() :
                           model.getVariable(position).getValue(0,1);
    }
    return values;
}

void makeSecondStepFail(RModel &model)
{
    // The first accepted step and its persisted checkpoint are valid. On a
    // restart, t0 and t1 still have a fixed support but t2 moves it. The
    // current baseline rejects moving supports during second-order dynamics.
    const double dt = model.getProblemSetup().getStructuralDynamicsSetup().getTimeStep();
    auto &support = model.getPoint(0).getBoundaryCondition(R_BOUNDARY_CONDITION_DISPLACEMENT);
    auto &x = support.getComponent(support.findComponentPosition(R_VARIABLE_DISPLACEMENT_X));
    x.add(dt,0.0);
    x.add(2.0*dt,0.01);
    model.getProblemSetup().setRestart(true);
    model.getProblemSetup().getStructuralDynamicsSetup().setNSteps(1);
}
}

class TestStructuralDynamicsFailureAtomicity : public QObject
{
    Q_OBJECT

private slots:
    void initTestCase();
    void failedRestartDoesNotAdvanceCheckpoint();
    void failedRestartRetainsPhysicalScaleAndPreviousResult();
    void failedResultWriteRollsBackAcceptedState();
    void sameSolverCanRetryAfterMidRunFailure();
    void resumePreservesEarlierRecord();
};

void TestStructuralDynamicsFailureAtomicity::initTestCase()
{
    omp_set_num_threads(1);
}

void TestStructuralDynamicsFailureAtomicity::failedRestartDoesNotAdvanceCheckpoint()
{
    RModel model = makeTruss();
    QString error;
    QVERIFY2(runSolver(model,error),qPrintable(error));
    const RStructuralDynamicsRestartData accepted =
        model.getProblemSetup().getStructuralDynamicsRestartData();
    QVERIFY(accepted.available);
    QCOMPARE(accepted.completedSteps,uint(1));
    makeSecondStepFail(model);

    error.clear();
    QVERIFY2(!runSolver(model,error),"A moving support should reject the trial step.");
    QVERIFY2(error.contains("Moving prescribed",Qt::CaseInsensitive),qPrintable(error));
    const auto &after = model.getProblemSetup().getStructuralDynamicsRestartData();
    QCOMPARE(after.completedSteps,accepted.completedSteps);
    QCOMPARE(after.time,accepted.time);
    QCOMPARE(after.fixedTimeStep,accepted.fixedTimeStep);
    QCOMPARE(after.systemFingerprint,accepted.systemFingerprint);
    QCOMPARE(after.q.size(),accepted.q.size());
    QCOMPARE(after.v.size(),accepted.v.size());
    QCOMPARE(after.a.size(),accepted.a.size());
    for (std::size_t i=0;i<accepted.q.size();++i)
    {
        QCOMPARE(after.q[i],accepted.q[i]);
        QCOMPARE(after.v[i],accepted.v[i]);
        QCOMPARE(after.a[i],accepted.a[i]);
    }
}

void TestStructuralDynamicsFailureAtomicity::failedRestartRetainsPhysicalScaleAndPreviousResult()
{
    RModel model = makeTruss();
    QString error;
    QVERIFY2(runSolver(model,error),qPrintable(error));
    const double acceptedTip = tipDisplacement(model);
    QVERIFY(std::isfinite(acceptedTip));
    const double physicalX = model.getNode(1).getX();
    const double physicalDt = model.getTimeSolver().getInputTimeStepSize();
    QCOMPARE(physicalX,length);
    makeSecondStepFail(model);

    error.clear();
    QVERIFY2(!runSolver(model,error),"A moving support should reject the trial step.");
    QVERIFY2(error.contains("Moving prescribed",Qt::CaseInsensitive),qPrintable(error));
    // A failed step must not leave the caller's mesh, stored result or time
    // input in internal solver units. This is intentionally a red regression
    // test for exception-safe scale restoration.
    QCOMPARE(model.getNode(1).getX(),physicalX);
    QCOMPARE(model.getTimeSolver().getInputTimeStepSize(),physicalDt);
    QVERIFY(std::isfinite(tipDisplacement(model)));
    QVERIFY(std::abs(tipDisplacement(model)-acceptedTip) < 1.0e-12*initialDisplacement);
    QCOMPARE(model.getProblemSetup().getStructuralDynamicsRestartData().completedSteps,uint(1));
}

void TestStructuralDynamicsFailureAtomicity::failedResultWriteRollsBackAcceptedState()
{
    QTemporaryDir temporary;
    QVERIFY(temporary.isValid());
    const QString latest = temporary.filePath("accepted.rbm");
    RModel model = makeTruss();
    QString error;
    QVERIFY2(runSolver(model,error,latest),qPrintable(error));
    QVERIFY(QFileInfo::exists(latest));
    const QString firstRecord = RFileUtils::getFileNameWithTimeStep(latest,1);
    auto readBytes = [](const QString &path)
    {
        QFile file(path);
        return file.open(QIODevice::ReadOnly) ? file.readAll() : QByteArray();
    };
    const QByteArray latestBytes = readBytes(latest);
    const QByteArray firstRecordBytes = readBytes(firstRecord);
    QVERIFY(!latestBytes.isEmpty());
    QVERIFY(!firstRecordBytes.isEmpty());

    const auto before = model.getProblemSetup().getStructuralDynamicsRestartData();
    const auto kinematics = tipKinematics(model);
    const uint acceptedStep = model.getTimeSolver().getCurrentTimeStep();
    const double acceptedTime = model.getTimeSolver().getCurrentTime();
    const double acceptedComputedTime = model.getTimeSolver().getComputedTime();
    QVERIFY(before.available);
    QCOMPARE(before.completedSteps,uint(1));
    for (double value : kinematics) QVERIFY(std::isfinite(value));

    // A regular file used as the parent of the second output path guarantees
    // ENOTDIR on all supported systems. The failure comes from writeResults,
    // after solve/process/store, and must not alter the first accepted files.
    model.getProblemSetup().setRestart(true);
    model.getProblemSetup().getStructuralDynamicsSetup().setNSteps(1);
    const QString blockingParent = temporary.filePath("blocked-parent");
    QFile blocker(blockingParent);
    QVERIFY(blocker.open(QIODevice::WriteOnly));
    QCOMPARE(blocker.write("block"),qint64(5));
    blocker.close();
    const QString invalidTarget = blockingParent+"/second.rbm";
    RSolver resumingSolver(model,invalidTarget,QString());
    error.clear();
    QVERIFY2(!runSolver(resumingSolver,error),"The blocked output parent unexpectedly succeeded.");
    QVERIFY2(error.contains("Failed to open the file",Qt::CaseInsensitive),qPrintable(error));

    // Publishing another record must never corrupt the latest accepted file.
    QVERIFY(QFileInfo::exists(latest));
    QCOMPARE(readBytes(latest),latestBytes);
    QCOMPARE(readBytes(firstRecord),firstRecordBytes);
    RModel persisted;
    persisted.read(latest);
    const auto persistedKinematics = tipKinematics(persisted);
    for (std::size_t i=0;i<kinematics.size();++i)
        QCOMPARE(persistedKinematics[i],kinematics[i]);

    // A failed write is still a failed *step*: in-memory results, restart
    // checkpoint and time must all describe the same previously accepted t1.
    const auto &after = model.getProblemSetup().getStructuralDynamicsRestartData();
    QCOMPARE(after.completedSteps,before.completedSteps);
    QCOMPARE(after.time,before.time);
    QCOMPARE(after.q.size(),before.q.size());
    for (std::size_t i=0;i<before.q.size();++i)
    {
        QCOMPARE(after.q[i],before.q[i]);
        QCOMPARE(after.v[i],before.v[i]);
        QCOMPARE(after.a[i],before.a[i]);
    }
    const auto afterKinematics = tipKinematics(model);
    for (std::size_t i=0;i<kinematics.size();++i)
        QCOMPARE(afterKinematics[i],kinematics[i]);
    QCOMPARE(model.getTimeSolver().getCurrentTimeStep(),acceptedStep);
    QCOMPARE(model.getTimeSolver().getCurrentTime(),acceptedTime);
    QCOMPARE(model.getTimeSolver().getComputedTime(),acceptedComputedTime);

    // The same solver object must be reusable after publication is repaired;
    // this also exercises restoration of its private integrator/force history.
    QVERIFY(blocker.rename(temporary.filePath("formerly-blocking-file")));
    QVERIFY(QDir().mkpath(blockingParent));
    error.clear();
    QVERIFY2(runSolver(resumingSolver,error),qPrintable(error));
    QCOMPARE(model.getProblemSetup().getStructuralDynamicsRestartData().completedSteps,uint(2));
}

void TestStructuralDynamicsFailureAtomicity::sameSolverCanRetryAfterMidRunFailure()
{
    RModel reference = makeTruss();
    reference.getProblemSetup().getStructuralDynamicsSetup().setNSteps(2);
    QString error;
    QVERIFY2(runSolver(reference,error),qPrintable(error));
    const auto expected = tipKinematics(reference);

    RModel model = makeTruss();
    auto &dynamics = model.getProblemSetup().getStructuralDynamicsSetup();
    dynamics.setNSteps(2);
    const double dt = dynamics.getTimeStep();
    auto &support = model.getPoint(0).getBoundaryCondition(R_BOUNDARY_CONDITION_DISPLACEMENT);
    auto &x = support.getComponent(support.findComponentPosition(R_VARIABLE_DISPLACEMENT_X));
    x.add(dt,0.0);
    x.add(2.0*dt,0.01);
    RSolver solver(model,QString(),QString());
    QVERIFY2(!runSolver(solver,error),"The second-step moving support unexpectedly succeeded.");
    QVERIFY2(error.contains("Moving prescribed",Qt::CaseInsensitive),qPrintable(error));
    QCOMPARE(model.getProblemSetup().getStructuralDynamicsRestartData().completedSteps,uint(1));
    QCOMPARE(model.getTimeSolver().getCurrentTimeStep(),uint(0));

    // Repair the input and resume with the *same* RSolver. Its integrator
    // already existed before the failed step and must retain accepted t1.
    // Model rollback replaces entity-group storage, so reacquire references
    // instead of using the pre-failure RConditionComponent reference `x`.
    auto &repairedSupport = model.getPoint(0).getBoundaryCondition(R_BOUNDARY_CONDITION_DISPLACEMENT);
    auto &repairedX = repairedSupport.getComponent(
        repairedSupport.findComponentPosition(R_VARIABLE_DISPLACEMENT_X));
    repairedX.clear();
    repairedX.add(0.0,0.0);
    model.getProblemSetup().setRestart(true);
    model.getProblemSetup().getStructuralDynamicsSetup().setNSteps(1);
    error.clear();
    QVERIFY2(runSolver(solver,error),qPrintable(error));
    QCOMPARE(model.getProblemSetup().getStructuralDynamicsRestartData().completedSteps,uint(2));
    const auto actual = tipKinematics(model);
    for (std::size_t i=0;i<expected.size();++i)
    {
        const double scale = std::max({1.0e-12,std::abs(actual[i]),std::abs(expected[i])});
        QVERIFY(std::abs(actual[i]-expected[i]) < 1.0e-9*scale+1.0e-12);
    }
}

void TestStructuralDynamicsFailureAtomicity::resumePreservesEarlierRecord()
{
    QTemporaryDir temporary;
    QVERIFY(temporary.isValid());
    const QString latest = temporary.filePath("history.rbm");
    RModel model = makeTruss();
    QString error;
    QVERIFY2(runSolver(model,error,latest),qPrintable(error));
    const QString firstRecord = RFileUtils::getFileNameWithTimeStep(latest,1);
    QFile first(firstRecord);
    QVERIFY(first.open(QIODevice::ReadOnly));
    const QByteArray acceptedBytes = first.readAll();
    QVERIFY(!acceptedBytes.isEmpty());
    first.close();

    model.getProblemSetup().setRestart(true);
    model.getProblemSetup().getStructuralDynamicsSetup().setNSteps(1);
    error.clear();
    QVERIFY2(runSolver(model,error,latest),qPrintable(error));
    // Only the latest link advances. A previously committed record is an
    // immutable restart point, not a file to rewrite/delete before solving.
    QVERIFY(QFileInfo::exists(firstRecord));
    QVERIFY(first.open(QIODevice::ReadOnly));
    QCOMPARE(first.readAll(),acceptedBytes);
    RModel oldState;
    oldState.read(firstRecord);
    QCOMPARE(oldState.getProblemSetup().getStructuralDynamicsRestartData().completedSteps,uint(1));
    RModel newest;
    newest.read(latest);
    QCOMPARE(newest.getProblemSetup().getStructuralDynamicsRestartData().completedSteps,uint(2));
}

QTEST_APPLESS_MAIN(TestStructuralDynamicsFailureAtomicity)

#include "tst_structural_dynamics_failure_atomicity.moc"

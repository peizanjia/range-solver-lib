#include <QtTest>

#include <QElapsedTimer>
#include <QFile>
#include <QTemporaryDir>

#include <omp.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <exception>
#include <limits>
#include <stdexcept>
#include <vector>

#include "rml_model.h"
#include "rparallel_execution_context.h"
#include "rsolver.h"

namespace
{
constexpr uint nElements = 32;
constexpr uint nSteps = 10;
constexpr uint nRepetitions = 3;
constexpr double length = 3.2;
constexpr double area = 0.005;
constexpr double youngModulus = 2.0e11;
constexpr double density = 8000.0;
constexpr double initialTipDisplacement = 1.0e-5;

struct RunResult
{
    QString error;
    double elapsedMs = 0.0;
    std::vector<double> displacement;
    std::vector<double> velocity;
    std::vector<double> acceleration;
    RStructuralDynamicsRestartData checkpoint;
    double mechanicalEnergy = 0.0;
    RParallelExecutionReport parallelReport;
    uint reportCalls = 0;
};

struct Policy
{
    const char *name;
    RParallelExecutionMode mode;
    uint threads;
};

RModel makeChain(const Policy &policy)
{
    RModel model;
    model.setNNodes(nElements+1);
    for (uint node=0;node<=nElements;node++)
    {
        model.getNode(node).set(length*double(node)/double(nElements),0.0,0.0);
    }
    for (uint elementId=0;elementId<nElements;elementId++)
    {
        RElement element(R_ELEMENT_TRUSS1);
        element.setNodeId(0,elementId);
        element.setNodeId(1,elementId+1);
        model.addElement(element,true,0);
    }

    RMaterial material(RMaterial::Solid);
    RMaterialProperty modulus(RMaterialProperty::ModulusOfElasiticity);
    modulus.add(293.15,youngModulus);
    material.add(modulus);
    RMaterialProperty poisson(RMaterialProperty::PoissonRatio);
    poisson.add(293.15,0.0);
    material.add(poisson);
    RMaterialProperty massDensity(RMaterialProperty::Density);
    massDensity.add(293.15,density);
    material.add(massDensity);
    RMaterialProperty expansion(RMaterialProperty::ThermalExpansionCoefficient);
    expansion.add(293.15,0.0);
    material.add(expansion);
    RLine &line = model.getLine(0);
    line.setCrossArea(area);
    line.setMaterial(material);

    // Truss transverse DOFs have no stiffness; the root also fixes axial motion.
    RBoundaryCondition transverse(R_BOUNDARY_CONDITION_DISPLACEMENT);
    transverse.getComponent(transverse.findComponentPosition(R_VARIABLE_DISPLACEMENT_X))
              .setEnabled(false);
    transverse.getComponent(transverse.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y))
              .add(0.0,0.0);
    transverse.getComponent(transverse.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z))
              .add(0.0,0.0);
    line.addBoundaryCondition(transverse);

    RElement rootPoint(R_ELEMENT_POINT);
    rootPoint.setNodeId(0,0);
    model.addElement(rootPoint,true,0);
    RBoundaryCondition fixed(R_BOUNDARY_CONDITION_DISPLACEMENT);
    for (RVariableType component : {R_VARIABLE_DISPLACEMENT_X,
                                    R_VARIABLE_DISPLACEMENT_Y,
                                    R_VARIABLE_DISPLACEMENT_Z})
    {
        fixed.getComponent(fixed.findComponentPosition(component)).add(0.0,0.0);
    }
    model.getPoint(0).addBoundaryCondition(fixed);

    RElement tipPoint(R_ELEMENT_POINT);
    tipPoint.setNodeId(0,nElements);
    model.addElement(tipPoint,true,1);
    RInitialCondition initial(R_INITIAL_CONDITION_DISPLACEMENT);
    initial.getComponent(initial.findComponentPosition(R_VARIABLE_DISPLACEMENT_X))
           .add(0.0,initialTipDisplacement);
    model.getPoint(1).addInitialCondition(initial);

    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    auto &linearSolver = model.getMatrixSolverConf(RMatrixSolverConf::CG);
    linearSolver.setOutputFrequency(0);
    linearSolver.setSolverCvgValue(1.0e-12);
    linearSolver.setNInnerIterations(250);
    linearSolver.setNOuterIterations(3000);
    auto &dynamics = model.getProblemSetup().getStructuralDynamicsSetup();
    dynamics.setEnabled(true);
    dynamics.setMethod(RStructuralTimeIntegrationMethod::Newmark);
    dynamics.setTimeStep(2.0e-6);
    dynamics.setNSteps(nSteps);
    dynamics.setOutputFrequency(nSteps);
    dynamics.setRayleighAlpha(25.0);

    auto &parallel = model.getProblemSetup().getParallelExecutionConfig();
    parallel.setMode(policy.mode);
    if (policy.mode == RParallelExecutionMode::FixedThreads)
    {
        parallel.setRequestedThreads(policy.threads);
    }
    parallel.setDeterministicReduction(true);
    return model;
}

bool extractNodalValues(const RModel &model, RVariableType type,
                        std::vector<double> &values)
{
    const uint index = model.findVariable(type);
    if (index == RConstants::eod)
    {
        return false;
    }
    values.reserve(nElements+1);
    for (uint node=0;node<=nElements;node++)
    {
        const double value = model.getVariable(index).getValue(0,node);
        if (!std::isfinite(value))
        {
            return false;
        }
        values.push_back(value);
    }
    return true;
}

double computeMechanicalEnergy(const RunResult &result)
{
    const double elementLength = length/double(nElements);
    const double elementMass = density*area*elementLength;
    const double elementStiffness = youngModulus*area/elementLength;
    double energy = 0.0;
    for (uint element=0;element<nElements;element++)
    {
        const double vi = result.velocity[element];
        const double vj = result.velocity[element+1];
        const double extension = result.displacement[element+1]-result.displacement[element];
        // v^T (rho*A*Le/6 [2 1; 1 2]) v / 2 + EA/(2*Le) * extension^2.
        energy += elementMass*(vi*vi+vi*vj+vj*vj)/6.0
                +0.5*elementStiffness*extension*extension;
    }
    return energy;
}

RunResult runCase(const Policy &policy)
{
    RunResult result;
    try
    {
        RModel model = makeChain(policy);
        QElapsedTimer timer;
        timer.start();
        RSolver solver(model,QString(),QString());
        solver.setParallelExecutionReportObserver(
            [&result](const RParallelExecutionReport &report)
            {
                result.parallelReport = report;
                ++result.reportCalls;
            });
        solver.run();
        result.elapsedMs = double(timer.nsecsElapsed())/1.0e6;

        if (result.reportCalls != 1 || result.parallelReport.mode != policy.mode ||
            result.parallelReport.backend != QStringLiteral("OpenMP") ||
            result.parallelReport.actualThreads < 1 ||
            result.parallelReport.actualThreads > result.parallelReport.appliedThreads)
        {
            result.error = "Solver did not report its scoped parallel context exactly once";
            return result;
        }

        if (!extractNodalValues(model,R_VARIABLE_DISPLACEMENT,result.displacement) ||
            !extractNodalValues(model,R_VARIABLE_VELOCITY,result.velocity) ||
            !extractNodalValues(model,R_VARIABLE_ACCELERATION,result.acceleration))
        {
            result.error = "Missing or non-finite nodal q/v/a output";
            return result;
        }
        result.checkpoint = model.getProblemSetup().getStructuralDynamicsRestartData();
        if (!result.checkpoint.available || result.checkpoint.completedSteps != nSteps)
        {
            result.error = "Dynamics checkpoint is missing or incomplete";
            return result;
        }
        result.mechanicalEnergy = computeMechanicalEnergy(result);
        if (!std::isfinite(result.mechanicalEnergy))
        {
            result.error = "Non-finite mechanical energy";
        }
    }
    catch (const RError &failure)
    {
        result.error = failure.getMessage();
    }
    catch (const std::exception &failure)
    {
        result.error = QString::fromUtf8(failure.what());
    }
    return result;
}

// All policies must agree numerically, but floating-point reduction order can
// differ between thread counts. The returned fraction is relative to the
// explicit mixed absolute/relative tolerance, not a claimed speedup.
QString compareScalar(const QString &label, double actual, double expected,
                      double absTolerance, double relTolerance, double &worstFraction)
{
    const double scale = std::max(std::abs(actual),std::abs(expected));
    const double tolerance = absTolerance+relTolerance*scale;
    const double fraction = std::abs(actual-expected)/tolerance;
    if (!std::isfinite(fraction))
    {
        return QString("%1 is non-finite (%2 vs %3)").arg(label).arg(actual).arg(expected);
    }
    worstFraction = std::max(worstFraction,fraction);
    if (fraction > 1.0)
    {
        return QString("%1 mismatch: %2 vs %3 (limit %4)")
            .arg(label).arg(actual,0,'g',17).arg(expected,0,'g',17)
            .arg(tolerance,0,'g',6);
    }
    return {};
}

QString compareVector(const QString &label, const std::vector<double> &actual,
                      const std::vector<double> &expected, double absTolerance,
                      double relTolerance, double &worstFraction)
{
    if (actual.size() != expected.size())
    {
        return QString("%1 size mismatch: %2 vs %3")
            .arg(label).arg(qulonglong(actual.size())).arg(qulonglong(expected.size()));
    }
    for (std::size_t i=0;i<actual.size();i++)
    {
        const QString error = compareScalar(QString("%1[%2]").arg(label).arg(qulonglong(i)),
                                            actual[i],expected[i],absTolerance,
                                            relTolerance,worstFraction);
        if (!error.isEmpty())
        {
            return error;
        }
    }
    return {};
}

QString compareRuns(const RunResult &actual, const RunResult &expected,
                    double &worstFraction)
{
    worstFraction = 0.0;
    struct VectorCheck
    {
        const char *name;
        const std::vector<double> *actual;
        const std::vector<double> *expected;
        double absTolerance;
    };
    const std::array<VectorCheck,7> vectors = {{
        {"nodal displacement",&actual.displacement,&expected.displacement,1.0e-12},
        {"nodal velocity",&actual.velocity,&expected.velocity,1.0e-9},
        {"nodal acceleration",&actual.acceleration,&expected.acceleration,1.0e-5},
        {"reduced displacement",&actual.checkpoint.q,&expected.checkpoint.q,1.0e-12},
        {"reduced velocity",&actual.checkpoint.v,&expected.checkpoint.v,1.0e-9},
        {"reduced acceleration",&actual.checkpoint.a,&expected.checkpoint.a,1.0e-5},
        {"reduced external force",&actual.checkpoint.externalForce,
         &expected.checkpoint.externalForce,1.0e-8}
    }};
    for (const VectorCheck &check : vectors)
    {
        const QString error = compareVector(check.name,*check.actual,*check.expected,
                                            check.absTolerance,1.0e-8,worstFraction);
        if (!error.isEmpty())
        {
            return error;
        }
    }
    struct ScalarCheck
    {
        const char *name;
        double actual;
        double expected;
        double absTolerance;
    };
    const std::array<ScalarCheck,5> scalars = {{
        {"physical mechanical energy",actual.mechanicalEnergy,expected.mechanicalEnergy,1.0e-8},
        {"initial mechanical energy",actual.checkpoint.initialMechanicalEnergy,
         expected.checkpoint.initialMechanicalEnergy,1.0e-8},
        {"external work",actual.checkpoint.cumulativeExternalWork,
         expected.checkpoint.cumulativeExternalWork,1.0e-8},
        {"damping dissipation",actual.checkpoint.cumulativeDampingDissipation,
         expected.checkpoint.cumulativeDampingDissipation,1.0e-8},
        {"time",actual.checkpoint.time,expected.checkpoint.time,1.0e-14}
    }};
    for (const ScalarCheck &check : scalars)
    {
        const QString error = compareScalar(check.name,check.actual,check.expected,
                                            check.absTolerance,1.0e-8,worstFraction);
        if (!error.isEmpty())
        {
            return error;
        }
    }
    return {};
}

double median(std::array<double,nRepetitions> values)
{
    std::sort(values.begin(),values.end());
    return values[nRepetitions/2];
}
}

class TestStructuralDynamicsParallelRepro : public QObject
{
    Q_OBJECT

private slots:
    void policiesAgreeAcrossRepeatedRuns();
    void observerFailureStopsBeforeTimeStepOrDiskWrite();
};

void TestStructuralDynamicsParallelRepro::policiesAgreeAcrossRepeatedRuns()
{
    const std::array<Policy,3> policies = {{{"serial",RParallelExecutionMode::SerialDeterministic,0},
                                            {"threads-2",RParallelExecutionMode::FixedThreads,2},
                                            {"threads-4",RParallelExecutionMode::FixedThreads,4}}};
    std::array<double,3> medians = {};
    RunResult serialReference;
    for (std::size_t policyIndex=0;policyIndex<policies.size();policyIndex++)
    {
        const Policy &policy = policies[policyIndex];
        std::array<double,nRepetitions> times = {};
        RunResult repeatReference;
        double maxCrossFraction = 0.0;
        double maxRepeatFraction = 0.0;
        for (uint repeat=0;repeat<nRepetitions;repeat++)
        {
            RunResult result = runCase(policy);
            QVERIFY2(result.error.isEmpty(),
                     qPrintable(QString("%1 run %2: %3")
                                .arg(QString::fromLatin1(policy.name)).arg(repeat+1)
                                .arg(result.error)));
            QVERIFY(result.mechanicalEnergy > 0.0);
            QVERIFY(result.checkpoint.initialMechanicalEnergy > 0.0);
            QVERIFY(result.checkpoint.cumulativeDampingDissipation > 0.0);
            times[repeat] = result.elapsedMs;
            if (policyIndex == 0 && repeat == 0)
            {
                serialReference = result;
            }
            if (repeat == 0)
            {
                repeatReference = result;
            }
            double fraction = 0.0;
            QString difference = compareRuns(result,serialReference,fraction);
            QVERIFY2(difference.isEmpty(),
                     qPrintable(QString("%1 run %2 vs serial: %3")
                                .arg(QString::fromLatin1(policy.name)).arg(repeat+1)
                                .arg(difference)));
            maxCrossFraction = std::max(maxCrossFraction,fraction);
            difference = compareRuns(result,repeatReference,fraction);
            QVERIFY2(difference.isEmpty(),
                     qPrintable(QString("%1 run %2 drift: %3")
                                .arg(QString::fromLatin1(policy.name)).arg(repeat+1)
                                .arg(difference)));
            maxRepeatFraction = std::max(maxRepeatFraction,fraction);
        }
        medians[policyIndex] = median(times);
        // Use stdout rather than Qt's message handler: the application logger
        // may filter QtInfoMsg, while CI should retain this diagnostic baseline.
        std::fprintf(stdout,
                     "BENCH %s: median %.3f ms; serial/this %.3fx; "
                     "max cross/repeat tolerance fractions %.4g/%.4g\n",
                     policy.name,medians[policyIndex],
                     medians[policyIndex] > 0.0
                         ? medians[0]/medians[policyIndex] : 0.0,
                     maxCrossFraction,maxRepeatFraction);
        std::fflush(stdout);
    }
    // Wall-clock numbers are diagnostic only: CI load and thread availability
    // vary, so this test intentionally has no minimum speedup assertion.
}

void TestStructuralDynamicsParallelRepro::observerFailureStopsBeforeTimeStepOrDiskWrite()
{
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString output = directory.filePath(QStringLiteral("observer-abort.ram"));
    const Policy policy{"serial",RParallelExecutionMode::SerialDeterministic,0};
    RModel model = makeChain(policy);
    RSolver solver(model,output,QString());
    const uint stepBefore = model.getTimeSolver().getCurrentTimeStep();
    const double computedTimeBefore = model.getTimeSolver().getComputedTime();
    const bool timeEnabledBefore = model.getTimeSolver().getEnabled();
    const int ompMaxBefore = omp_get_max_threads();
    uint reportCalls = 0;
    RParallelExecutionReport observed;
    solver.setParallelExecutionReportObserver(
        [&](const RParallelExecutionReport &report)
        {
            observed = report;
            ++reportCalls;
            throw std::runtime_error("observer refused run");
        });
    bool rejected = false;
    try
    {
        solver.run();
    }
    catch (const std::runtime_error &error)
    {
        rejected = QString::fromUtf8(error.what()) == QStringLiteral("observer refused run");
    }
    QVERIFY(rejected);
    QCOMPARE(reportCalls,1u);
    QCOMPARE(observed.mode,RParallelExecutionMode::SerialDeterministic);
    QCOMPARE(observed.actualThreads,1u);
    QCOMPARE(omp_get_max_threads(),ompMaxBefore);
    QCOMPARE(model.getTimeSolver().getCurrentTimeStep(),stepBefore);
    QCOMPARE(model.getTimeSolver().getComputedTime(),computedTimeBefore);
    QCOMPARE(model.getTimeSolver().getEnabled(),timeEnabledBefore);
    QVERIFY(!model.getProblemSetup().getStructuralDynamicsRestartData().available);
    QVERIFY(!QFile::exists(output));
}

QTEST_APPLESS_MAIN(TestStructuralDynamicsParallelRepro)
#include "tst_structural_dynamics_parallel_repro.moc"

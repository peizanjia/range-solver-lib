#include <QtTest>

#include <QFile>
#include <QFileInfo>
#include <QTemporaryDir>

#include <condition_variable>
#include <chrono>
#include <exception>
#include <mutex>
#include <thread>

#include <rbl_error.h>
#include <rml_monitoring_point.h>
#include <rml_model.h>

#include "rsolver.h"

namespace
{
void setComponent(RBoundaryCondition &condition,RVariableType type,double value)
{
    const uint position = condition.findComponentPosition(type);
    Q_ASSERT(position != RConstants::eod);
    RConditionComponent &component = condition.getComponent(position);
    component.setEnabled(true);
    component.clear();
    component.add(0.0,value);
}

void setMaterialProperty(RMaterial &material,RMaterialProperty::Type type,double value)
{
    RMaterialProperty property(type);
    property.add(293.15,value);
    material.add(property);
}

RModel transientRod(const QString &monitorFileName)
{
    RModel model;
    model.setNNodes(2);
    model.getNode(0).set(0.0,0.0,0.0);
    model.getNode(1).set(1.0,0.0,0.0);

    RElement rod(R_ELEMENT_TRUSS1);
    rod.setNodeId(0,0);
    rod.setNodeId(1,1);
    model.addElement(rod,true,0);

    RMaterial material(RMaterial::Solid);
    setMaterialProperty(material,RMaterialProperty::ModulusOfElasiticity,2.0e11);
    setMaterialProperty(material,RMaterialProperty::PoissonRatio,0.0);
    setMaterialProperty(material,RMaterialProperty::Density,8000.0);
    setMaterialProperty(material,RMaterialProperty::ThermalExpansionCoefficient,0.0);
    setMaterialProperty(material,RMaterialProperty::ThermalConductivity,2.0);
    setMaterialProperty(material,RMaterialProperty::HeatCapacity,500.0);
    setMaterialProperty(material,RMaterialProperty::Emissivity,0.0);
    RLine &line = model.getLine(0);
    line.setCrossArea(0.01);
    line.setMaterial(material);

    RBoundaryCondition lateral(R_BOUNDARY_CONDITION_DISPLACEMENT);
    lateral.getComponent(lateral.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).setEnabled(false);
    setComponent(lateral,R_VARIABLE_DISPLACEMENT_Y,0.0);
    setComponent(lateral,R_VARIABLE_DISPLACEMENT_Z,0.0);
    line.addBoundaryCondition(lateral);

    for (uint node=0;node<2;node++)
    {
        RElement point(R_ELEMENT_POINT);
        point.setNodeId(0,node);
        model.addElement(point,true,node);
        model.getPoint(node).setVolume(0.0);

        RBoundaryCondition temperature(R_BOUNDARY_CONDITION_TEMPERATURE);
        setComponent(temperature,R_VARIABLE_TEMPERATURE,293.15);
        model.getPoint(node).addBoundaryCondition(temperature);

        if (node == 0)
        {
            RBoundaryCondition fixed(R_BOUNDARY_CONDITION_DISPLACEMENT);
            setComponent(fixed,R_VARIABLE_DISPLACEMENT_X,0.0);
            setComponent(fixed,R_VARIABLE_DISPLACEMENT_Y,0.0);
            setComponent(fixed,R_VARIABLE_DISPLACEMENT_Z,0.0);
            model.getPoint(node).addBoundaryCondition(fixed);
        }
    }

    RProblemTaskItem sequence;
    sequence.addChild(RProblemTaskItem(R_PROBLEM_HEAT));
    sequence.addChild(RProblemTaskItem(R_PROBLEM_STRESS));
    model.setProblemTaskTree(sequence);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setOutputFrequency(0);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setSolverCvgValue(1.0e-12);

    RMonitoringPoint monitor(R_VARIABLE_TEMPERATURE,RR3Vector(0.5,0.0,0.0));
    monitor.setOutputFileName(monitorFileName);
    model.getMonitoringPointManager().push_back(monitor);

    RTimeSolver &timeSolver = model.getTimeSolver();
    timeSolver.setEnabled(true);
    timeSolver.setInputNTimeSteps(1);
    timeSolver.setInputStartTime(0.0);
    timeSolver.setInputTimeStepSize(1.0);
    timeSolver.setOutputFrequency(1);
    return model;
}

QByteArray readBytes(const QString &fileName)
{
    QFile file(fileName);
    return file.open(QIODevice::ReadOnly) ? file.readAll() : QByteArray();
}
}

class TestStructuralMonitorOutputConcurrency : public QObject
{
    Q_OBJECT
private slots:
    void transientRunsHoldSharedMonitorLockForWholeRun();
};

void TestStructuralMonitorOutputConcurrency::transientRunsHoldSharedMonitorLockForWholeRun()
{
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString monitorFileName = directory.filePath(QStringLiteral("shared-monitor.txt"));
    const QString firstModelFileName = directory.filePath(QStringLiteral("first-model.ram"));
    const QString secondModelFileName = directory.filePath(QStringLiteral("second-model.ram"));

    RModel firstModel = transientRod(monitorFileName);
    RModel secondModel = transientRod(monitorFileName);
    RSolver firstSolver(firstModel,firstModelFileName,QString());
    RSolver secondSolver(secondModel,secondModelFileName,QString());

    std::mutex mutex;
    std::condition_variable condition;
    bool observerEntered = false;
    bool resumeFirstRun = false;
    std::exception_ptr firstFailure;
    firstSolver.setParallelExecutionReportObserver(
        [&mutex,&condition,&observerEntered,&resumeFirstRun](const RParallelExecutionReport &)
        {
            std::unique_lock<std::mutex> lock(mutex);
            observerEntered = true;
            condition.notify_one();
            condition.wait(lock,[&resumeFirstRun]() { return resumeFirstRun; });
        });

    std::thread firstRun([&firstSolver,&firstFailure]()
    {
        try
        {
            firstSolver.run();
        }
        catch (...)
        {
            firstFailure = std::current_exception();
        }
    });

    bool didReachObserver = false;
    {
        std::unique_lock<std::mutex> lock(mutex);
        didReachObserver = condition.wait_for(lock,std::chrono::seconds(10),
                                              [&observerEntered]() { return observerEntered; });
    }
    if (!didReachObserver)
    {
        {
            std::lock_guard<std::mutex> lock(mutex);
            resumeFirstRun = true;
        }
        condition.notify_one();
        firstRun.join();
        QFAIL("First transient run did not reach its synchronized observer.");
    }

    bool secondRunWasRejectedForTheSharedMonitor = false;
    QString secondRunFailure;
    try
    {
        secondSolver.run();
    }
    catch (const RError &error)
    {
        secondRunFailure = error.getMessage();
        secondRunWasRejectedForTheSharedMonitor =
            secondRunFailure.contains(QStringLiteral("Transient monitoring output is locked"));
    }
    catch (...)
    {
        secondRunFailure = QStringLiteral("Second transient run failed with a non-RError exception.");
    }
    const bool secondModelWasPublished = QFileInfo::exists(secondModelFileName);

    {
        std::lock_guard<std::mutex> lock(mutex);
        resumeFirstRun = true;
    }
    condition.notify_one();
    firstRun.join();

    QVERIFY2(secondRunWasRejectedForTheSharedMonitor,
             qPrintable(QStringLiteral("Second run did not fail on the shared monitor lock: %1")
                            .arg(secondRunFailure)));
    QVERIFY(!secondModelWasPublished);

    if (firstFailure)
    {
        try
        {
            std::rethrow_exception(firstFailure);
        }
        catch (const std::exception &error)
        {
            QFAIL(qPrintable(QStringLiteral("First transient run failed: %1").arg(error.what())));
        }
        catch (...)
        {
            QFAIL("First transient run failed with a non-standard exception.");
        }
    }
    QVERIFY(QFileInfo::exists(firstModelFileName));
    QCOMPARE(readBytes(monitorFileName).count('\n'),2);

    // The run-level lock is released when the first transient run exits.
    secondSolver.run();
    QVERIFY(QFileInfo::exists(secondModelFileName));
    QCOMPARE(readBytes(monitorFileName).count('\n'),4);
}

QTEST_GUILESS_MAIN(TestStructuralMonitorOutputConcurrency)

#include "tst_structural_monitor_output_concurrency.moc"

#include <QtTest>

#include <QDir>
#include <QFileInfo>
#include <QTemporaryDir>

#include <cmath>
#include <thread>
#include <utility>
#include <vector>

#include "rml_file_utils.h"
#include "rml_monitoring_point.h"
#include "rml_model.h"
#include "rsolver.h"

namespace
{
constexpr double length = 2.0;
constexpr double initialDisplacement = 1.0e-4;

RModel makeTruss(uint steps, uint outputFrequency)
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

    RElement fixed(R_ELEMENT_POINT);
    fixed.setNodeId(0,0);
    model.addElement(fixed,true,0);
    RBoundaryCondition support(R_BOUNDARY_CONDITION_DISPLACEMENT);
    for (RVariableType type : {R_VARIABLE_DISPLACEMENT_X,
                               R_VARIABLE_DISPLACEMENT_Y,
                               R_VARIABLE_DISPLACEMENT_Z})
    {
        support.getComponent(support.findComponentPosition(type)).add(0.0,0.0);
    }
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
    dynamics.setNSteps(steps);
    dynamics.setOutputFrequency(outputFrequency);
    return model;
}

void setComponent(RBoundaryCondition &condition, RVariableType type, double value)
{
    const uint position = condition.findComponentPosition(type);
    condition.getComponent(position).setEnabled(true);
    condition.getComponent(position).clear();
    condition.getComponent(position).add(0.0,value);
}

RModel makeTransientHeatStressRod()
{
    RModel model;
    model.setNNodes(3);
    for (uint node=0;node<3;++node)
        model.getNode(node).set(double(node),0.0,0.0);
    for (uint segment=0;segment<2;++segment)
    {
        RElement bar(R_ELEMENT_TRUSS1);
        bar.setNodeId(0,segment);
        bar.setNodeId(1,segment+1);
        model.addElement(bar,true,0);
    }

    RMaterial material(RMaterial::Solid);
    for (const auto &property : {
             std::pair{RMaterialProperty::ModulusOfElasiticity,2.0e11},
             std::pair{RMaterialProperty::PoissonRatio,0.0},
             std::pair{RMaterialProperty::Density,8000.0},
             std::pair{RMaterialProperty::ThermalExpansionCoefficient,0.0},
             std::pair{RMaterialProperty::ThermalConductivity,2.0},
             std::pair{RMaterialProperty::HeatCapacity,500.0},
             std::pair{RMaterialProperty::Emissivity,0.0}})
    {
        RMaterialProperty entry(property.first);
        entry.add(293.15,property.second);
        material.add(entry);
    }
    RLine &line = model.getLine(0);
    line.setCrossArea(0.01);
    line.setMaterial(material);

    RBoundaryCondition lateral(R_BOUNDARY_CONDITION_DISPLACEMENT);
    lateral.getComponent(lateral.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).setEnabled(false);
    setComponent(lateral,R_VARIABLE_DISPLACEMENT_Y,0.0);
    setComponent(lateral,R_VARIABLE_DISPLACEMENT_Z,0.0);
    line.addBoundaryCondition(lateral);

    for (uint node=0;node<3;++node)
    {
        RElement point(R_ELEMENT_POINT);
        point.setNodeId(0,node);
        model.addElement(point,true,node);
        model.getPoint(node).setVolume(0.0);
        if (node == 0 || node == 2)
        {
            RBoundaryCondition temperature(R_BOUNDARY_CONDITION_TEMPERATURE);
            setComponent(temperature,R_VARIABLE_TEMPERATURE,node == 0 ? 300.0 : 400.0);
            model.getPoint(node).addBoundaryCondition(temperature);
        }
        if (node == 0)
        {
            RBoundaryCondition support(R_BOUNDARY_CONDITION_DISPLACEMENT);
            setComponent(support,R_VARIABLE_DISPLACEMENT_X,0.0);
            setComponent(support,R_VARIABLE_DISPLACEMENT_Y,0.0);
            setComponent(support,R_VARIABLE_DISPLACEMENT_Z,0.0);
            model.getPoint(node).addBoundaryCondition(support);
        }
    }
    RProblemTaskItem sequence;
    sequence.addChild(RProblemTaskItem(R_PROBLEM_HEAT));
    sequence.addChild(RProblemTaskItem(R_PROBLEM_STRESS));
    model.setProblemTaskTree(sequence);
    auto &linear = model.getMatrixSolverConf(RMatrixSolverConf::CG);
    linear.setOutputFrequency(0);
    linear.setSolverCvgValue(1.0e-12);

    RTimeSolver &time = model.getTimeSolver();
    time.setEnabled(true);
    time.setInputNTimeSteps(1);
    time.setInputStartTime(0.0);
    time.setInputTimeStepSize(500000.0);
    time.setOutputFrequency(1);
    return model;
}

void configureZeroControl(RModel &model)
{
    auto &setup = model.getProblemSetup().getStructuralControlSetup();
    setup.enabled = true;
    setup.controllerId = QStringLiteral("result-frame-test-controller");
    setup.deterministicSeed = 17;

    RStructuralActuator actuator;
    actuator.id = QStringLiteral("tip-x");
    actuator.target.nodeId = 1;
    actuator.variable = R_VARIABLE_FORCE_X;
    actuator.unit = RVariable::getUnits(R_VARIABLE_FORCE_X);
    actuator.minimum = -5.0e4;
    actuator.maximum = 5.0e4;
    model.getProblemSetup().getStructuralInteractionSetup().getActuators().push_back(actuator);

    RStructuralSensor sensor;
    sensor.id = QStringLiteral("tip-displacement");
    sensor.target.nodeId = 1;
    sensor.quantity = RStructuralSensorQuantity::Displacement;
    sensor.direction = RR3Vector(1.0,0.0,0.0);
    sensor.unit = RStructuralInteractionSetup::getSensorUnit(sensor.quantity);
    model.getProblemSetup().getStructuralInteractionSetup().getSensors().push_back(sensor);
}

QString runFailure(RSolver &solver)
{
    try
    {
        solver.run();
    }
    catch (const RError &error)
    {
        return error.getMessage();
    }
    catch (const std::exception &error)
    {
        return QString::fromUtf8(error.what());
    }
    return QString();
}
}

class TestNativeResultFrameObserver : public QObject
{
    Q_OBJECT

private slots:
    void reportsOnlyPersistedSparseFramesOnOwnerThread();
    void activeRunUsesTheSameFrameObserver();
    void rejectedStepDoesNotPublishItsFrame();
    void enclosingRollbackDiscardsAQueuedPersistedFrame();
};

void TestNativeResultFrameObserver::reportsOnlyPersistedSparseFramesOnOwnerThread()
{
    QTemporaryDir temporary;
    QVERIFY(temporary.isValid());
    const QString logicalPath = temporary.filePath(QStringLiteral("history.rbm"));
    RModel model = makeTruss(4,2);
    RSolver solver(model,logicalPath,QString());

    const std::thread::id ownerThread = std::this_thread::get_id();
    std::vector<RSolverResultFrame> frames;
    QString callbackError;
    bool sameThread = true;
    solver.setResultFrameObserver([&](const RSolverResultFrame &frame)
    {
        sameThread = sameThread && std::this_thread::get_id() == ownerThread;
        frames.push_back(frame);
        if (!QFileInfo::exists(frame.actualRecordPath) ||
            !QFileInfo::exists(frame.logicalPath))
        {
            callbackError = QStringLiteral("Observer ran before the result record or logical link existed");
            return;
        }
        RModel persisted;
        try
        {
            persisted.read(frame.actualRecordPath);
        }
        catch (const RError &error)
        {
            callbackError = error.getMessage();
            return;
        }
        const auto &checkpoint = persisted.getProblemSetup().getStructuralDynamicsRestartData();
        if (!checkpoint.available || checkpoint.completedSteps != frame.recordNo ||
            checkpoint.time != frame.checkpointTime)
        {
            callbackError = QStringLiteral("Published record does not contain its accepted checkpoint");
        }
    });

    const QString failure = runFailure(solver);
    QVERIFY2(failure.isEmpty(),qPrintable(failure));
    QVERIFY2(callbackError.isEmpty(),qPrintable(callbackError));
    QVERIFY(sameThread);
    QCOMPARE(frames.size(),std::size_t(2));
    QCOMPARE(frames[0].recordNo,uint(2));
    QCOMPARE(frames[1].recordNo,uint(4));
    for (const RSolverResultFrame &frame : frames)
    {
        QCOMPARE(frame.nativeAnalysis,QStringLiteral("transient"));
        QCOMPARE(frame.mode,-1);
        QCOMPARE(frame.logicalPath,QFileInfo(logicalPath).absoluteFilePath());
        QCOMPARE(frame.actualRecordPath,
                 QFileInfo(RFileUtils::getFileNameWithTimeStep(logicalPath,frame.recordNo))
                     .absoluteFilePath());
        QVERIFY(std::isfinite(frame.checkpointTime));
    }
}

void TestNativeResultFrameObserver::rejectedStepDoesNotPublishItsFrame()
{
    QTemporaryDir temporary;
    QVERIFY(temporary.isValid());
    const QString logicalPath = temporary.filePath(QStringLiteral("history.rbm"));
    const QString blockedRecord = RFileUtils::getFileNameWithTimeStep(logicalPath,2);
    QVERIFY(QDir().mkpath(blockedRecord));

    RModel model = makeTruss(2,1);
    RSolver solver(model,logicalPath,QString());
    std::vector<RSolverResultFrame> frames;
    solver.setResultFrameObserver([&frames](const RSolverResultFrame &frame)
    {
        frames.push_back(frame);
    });

    const QString failure = runFailure(solver);
    QVERIFY2(!failure.isEmpty(),"The blocked second output record unexpectedly succeeded");
    QCOMPARE(frames.size(),std::size_t(1));
    QCOMPARE(frames.front().recordNo,uint(1));
    QCOMPARE(model.getProblemSetup().getStructuralDynamicsRestartData().completedSteps,uint(1));
    QCOMPARE(model.getTimeSolver().getCurrentTimeStep(),uint(0));
    QVERIFY(QFileInfo::exists(frames.front().actualRecordPath));
    QVERIFY(QFileInfo(blockedRecord).isDir());
}

void TestNativeResultFrameObserver::activeRunUsesTheSameFrameObserver()
{
    QTemporaryDir temporary;
    QVERIFY(temporary.isValid());
    RModel model = makeTruss(2,1);
    configureZeroControl(model);
    RSolver solver(model,temporary.filePath(QStringLiteral("active.rbm")),QString());
    solver.setStructuralControlCallback(
        QStringLiteral("result-frame-test-controller"),QStringLiteral("test-revision"),
        [](const RStructuralControlObservation &)
        {
            RStructuralControlDecision decision;
            decision.requestedCommands = {{QStringLiteral("tip-x"),0.0}};
            return decision;
        });

    std::vector<RSolverResultFrame> frames;
    solver.setResultFrameObserver([&frames](const RSolverResultFrame &frame)
    {
        frames.push_back(frame);
    });
    const QString failure = runFailure(solver);
    QVERIFY2(failure.isEmpty(),qPrintable(failure));
    QCOMPARE(frames.size(),std::size_t(2));
    QCOMPARE(frames[0].recordNo,uint(1));
    QCOMPARE(frames[1].recordNo,uint(2));
}

void TestNativeResultFrameObserver::enclosingRollbackDiscardsAQueuedPersistedFrame()
{
    QTemporaryDir temporary;
    QVERIFY(temporary.isValid());
    const QString logicalPath = temporary.filePath(QStringLiteral("coupled.ram"));
    const QString firstMonitorPath = temporary.filePath(QStringLiteral("monitor-first.txt"));
    const QString blockedMonitorPath = temporary.filePath(QStringLiteral("monitor-blocked.txt"));
    QFile previousMonitor(firstMonitorPath);
    QVERIFY(previousMonitor.open(QIODevice::WriteOnly));
    QCOMPARE(previousMonitor.write(QByteArrayLiteral("previous accepted monitor data\n")),
             qint64(QByteArrayLiteral("previous accepted monitor data\n").size()));
    previousMonitor.close();
    QVERIFY(QDir().mkpath(blockedMonitorPath));

    RModel model = makeTransientHeatStressRod();
    RMonitoringPoint firstPoint(R_VARIABLE_TEMPERATURE,RR3Vector(1.0,0.0,0.0));
    firstPoint.setOutputFileName(firstMonitorPath);
    RMonitoringPoint blockedPoint(R_VARIABLE_TEMPERATURE,RR3Vector(1.0,0.0,0.0));
    blockedPoint.setOutputFileName(blockedMonitorPath);
    model.getMonitoringPointManager().push_back(firstPoint);
    model.getMonitoringPointManager().push_back(blockedPoint);

    RSolver solver(model,logicalPath,QString());
    std::vector<RSolverResultFrame> frames;
    solver.setResultFrameObserver([&frames](const RSolverResultFrame &frame)
    {
        frames.push_back(frame);
    });

    QVERIFY_EXCEPTION_THROWN(solver.run(),RError);
    QVERIFY(frames.empty());
    QVERIFY(!QFileInfo::exists(logicalPath));
    QVERIFY(!QFileInfo::exists(RFileUtils::getFileNameWithTimeStep(logicalPath,1)));
    QVERIFY(QFileInfo(blockedMonitorPath).isDir());

    QVERIFY(QDir(blockedMonitorPath).removeRecursively());
    solver.run();
    QCOMPARE(frames.size(),std::size_t(1));
    QCOMPARE(frames.front().recordNo,uint(1));
    QCOMPARE(frames.front().nativeAnalysis,QStringLiteral("transient"));
    QVERIFY(QFileInfo::exists(frames.front().actualRecordPath));
    RModel persisted;
    persisted.read(frames.front().actualRecordPath);
    QCOMPARE(persisted.getTimeSolver().getComputedTime(),frames.front().checkpointTime);
}

QTEST_APPLESS_MAIN(TestNativeResultFrameObserver)

#include "tst_native_result_frame_observer.moc"

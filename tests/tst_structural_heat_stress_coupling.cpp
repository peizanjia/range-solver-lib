#include <QtTest>

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QLockFile>
#include <QTemporaryDir>

#include <cmath>
#include <limits>

#include <rbl_error.h>
#include <rml_file_utils.h>
#include <rml_monitoring_point.h>
#include <rml_model.h>

#include "rsolver.h"

namespace
{
constexpr double referenceKelvin = 293.15;
constexpr double modulus = 2.0e11;
constexpr double expansion = 1.2e-5;

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
    RMaterialProperty entry(type);
    entry.add(referenceKelvin,value);
    material.add(entry);
}

RModel heatedRod(bool fixRight,
                 bool temperatureDependentExpansion = false,
                 bool elasticModulusDefined = true,
                 double thermalExpansion = expansion,
                 bool densityDefined = true)
{
    RModel model;
    model.setNNodes(3);
    for (uint node=0;node<3;node++)
        model.getNode(node).set(double(node),0.0,0.0);
    for (uint segment=0;segment<2;segment++)
    {
        RElement rod(R_ELEMENT_TRUSS1);
        rod.setNodeId(0,segment);
        rod.setNodeId(1,segment+1);
        model.addElement(rod,true,0);
    }

    RMaterial material(RMaterial::Solid);
    if (elasticModulusDefined)
        setMaterialProperty(material,RMaterialProperty::ModulusOfElasiticity,modulus);
    setMaterialProperty(material,RMaterialProperty::PoissonRatio,0.0);
    if (densityDefined)
        setMaterialProperty(material,RMaterialProperty::Density,8000.0);
    if (temperatureDependentExpansion)
    {
        RMaterialProperty variableAlpha(RMaterialProperty::ThermalExpansionCoefficient);
        variableAlpha.add(referenceKelvin,1.0e-5);
        variableAlpha.add(referenceKelvin+100.0,3.0e-5);
        material.add(variableAlpha);
    }
    else
    {
        setMaterialProperty(material,RMaterialProperty::ThermalExpansionCoefficient,thermalExpansion);
    }
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

    for (uint node=0;node<3;node++)
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
        if (node == 0 || (fixRight && node == 2))
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
    return model;
}

double nodal(const RModel &model,RVariableType type,uint component,uint node)
{
    const uint variable = model.findVariable(type);
    if (variable == RConstants::eod) return std::numeric_limits<double>::quiet_NaN();
    return model.getVariable(variable).getValue(component,node);
}

double elemental(const RModel &model,RVariableType type,uint component,uint element)
{
    const uint variable = model.findVariable(type);
    if (variable == RConstants::eod) return std::numeric_limits<double>::quiet_NaN();
    return model.getVariable(variable).getValue(component,element);
}

QByteArray readBytes(const QString &path)
{
    QFile file(path);
    return file.open(QIODevice::ReadOnly) ? file.readAll() : QByteArray();
}

void seedAcceptedArtifact(const QString &latestPath)
{
    RModel accepted;
    accepted.setName(QStringLiteral("accepted-before-heat-stress"));
    accepted.getTimeSolver().setEnabled(true);
    accepted.getTimeSolver().setTimes({0.01});
    accepted.write(latestPath);
}
}

class TestStructuralHeatStressCoupling : public QObject
{
    Q_OBJECT
private slots:
    void heatThenFreeExpansion();
    void heatThenConstrainedCompression();
    void temperatureDependentExpansionUsesIntegratedAlpha();
    void persistedTemperatureChangesModalFrequency();
    void stressFailureDoesNotPublishAndCanRetry();
    void resultWriteFailureRollsBackAndCanRetry();
    void transientSequenceAdvancesHeatThenStressAtEachTime();
    void transientStressFailureRollsBackAndPreservesRecordsForRetry();
    void transientPublicationFailureRetriesSameTimeIndex();
    void transientRestartRejectsUnpublishedFutureRecords();
    void transientRestartFindsSixDigitOutputRecord();
    void transientOutputLockRejectsConcurrentWriter();
    void transientRejectsMonitoringOutputBeforeStep();
    void unsupportedTransientSequenceIsRejectedBeforeOutputCleanup();
};

void TestStructuralHeatStressCoupling::heatThenFreeExpansion()
{
    RModel model = heatedRod(false);
    RSolver solver(model,QString(),QString());
    solver.run();
    QVERIFY(std::fabs(nodal(model,R_VARIABLE_TEMPERATURE,0,0)-300.0) < 1.0e-8);
    QVERIFY(std::fabs(nodal(model,R_VARIABLE_TEMPERATURE,0,1)-350.0) < 1.0e-8);
    QVERIFY(std::fabs(nodal(model,R_VARIABLE_TEMPERATURE,0,2)-400.0) < 1.0e-8);
    const double expectedTip = expansion*((325.0-referenceKelvin)+(375.0-referenceKelvin));
    QVERIFY(std::fabs(nodal(model,R_VARIABLE_DISPLACEMENT,0,2)-expectedTip) < 1.0e-7);
    QVERIFY(std::fabs(elemental(model,R_VARIABLE_STRESS_X,0,0)) < 1.0e5);
    QVERIFY(std::fabs(elemental(model,R_VARIABLE_STRESS_X,0,1)) < 1.0e5);
}

void TestStructuralHeatStressCoupling::heatThenConstrainedCompression()
{
    RModel model = heatedRod(true);
    RSolver solver(model,QString(),QString());
    solver.run();
    const double expectedStress = -modulus*expansion*(350.0-referenceKelvin);
    QVERIFY(std::fabs(nodal(model,R_VARIABLE_DISPLACEMENT,0,2)) < 1.0e-9);
    QVERIFY(std::fabs(elemental(model,R_VARIABLE_STRESS_X,0,0)-expectedStress) < 1.0e5);
    QVERIFY(std::fabs(elemental(model,R_VARIABLE_STRESS_X,0,1)-expectedStress) < 1.0e5);
}

void TestStructuralHeatStressCoupling::temperatureDependentExpansionUsesIntegratedAlpha()
{
    RModel model = heatedRod(false,true);
    RSolver solver(model,QString(),QString());
    solver.run();
    const double delta0 = 325.0-referenceKelvin;
    const double delta1 = 375.0-referenceKelvin;
    // alpha(T)=1e-5+2e-7*(T-Tref); its integral is not alpha(Tend)*DeltaT.
    const double expectedTip = 1.0e-5*(delta0+delta1)
                             + 1.0e-7*(delta0*delta0+delta1*delta1);
    QVERIFY(std::fabs(nodal(model,R_VARIABLE_DISPLACEMENT,0,2)-expectedTip) < 1.0e-7);
    QVERIFY(std::fabs(elemental(model,R_VARIABLE_STRESS_X,0,0)) < 1.0e5);
    QVERIFY(std::fabs(elemental(model,R_VARIABLE_STRESS_X,0,1)) < 1.0e5);
}

void TestStructuralHeatStressCoupling::persistedTemperatureChangesModalFrequency()
{
    auto frequencyAt = [](double kelvin) {
        RModel model = heatedRod(false);
        RMaterial material = model.getLine(0).getMaterial();
        RMaterialProperty youngsModulus(RMaterialProperty::ModulusOfElasiticity);
        youngsModulus.add(300.0,modulus);
        youngsModulus.add(400.0,0.5*modulus);
        material.add(youngsModulus);
        model.getLine(0).setMaterial(material);

        const uint position = model.addVariable(R_VARIABLE_TEMPERATURE);
        RVariable &temperature = model.getVariable(position);
        temperature.setApplyType(R_VARIABLE_APPLY_NODE);
        temperature.resize(1,model.getNNodes());
        for (uint node=0;node<model.getNNodes();node++)
            temperature.setValue(0,node,kelvin);

        model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS_MODAL));
        RModalSetup &modal = model.getProblemSetup().getModalSetup();
        modal.setMethod(R_MODAL_MULTIPLE_MODES);
        modal.setNModesToExtract(1);
        modal.setNIterations(300);
        modal.setConvergenceValue(1.0e-12);
        RSolver solver(model,QString(),QString());
        solver.run();
        return modal.getFrequency();
    };

    const double cold = frequencyAt(300.0);
    const double hot = frequencyAt(400.0);
    QVERIFY(std::isfinite(cold) && cold > 0.0);
    QVERIFY(std::isfinite(hot) && hot > 0.0);
    const double expectedRatio = std::sqrt(0.5);
    QVERIFY2(std::fabs(hot/cold-expectedRatio) < 0.01,
             qPrintable(QString("Thermal modal ratio = %1, expected %2")
                        .arg(hot/cold).arg(expectedRatio)));
}

void TestStructuralHeatStressCoupling::stressFailureDoesNotPublishAndCanRetry()
{
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString latestPath = directory.filePath("coupled.ram");
    seedAcceptedArtifact(latestPath);
    const QByteArray acceptedLink = readBytes(latestPath);
    QVERIFY(!acceptedLink.isEmpty());
    const QString trialRecord = RFileUtils::getFileNameWithTimeStep(latestPath,0);
    QVERIFY(!QFileInfo::exists(trialRecord));

    RModel model = heatedRod(false,false,false);
    RBoundaryCondition load(R_BOUNDARY_CONDITION_FORCE);
    setComponent(load,R_VARIABLE_FORCE_X,1000.0);
    model.getPoint(2).addBoundaryCondition(load);
    RSolver solver(model,latestPath,QString());

    // Heat succeeds, then the loaded zero-stiffness stress problem fails.
    // Neither the in-memory temperature nor a new file record may escape.
    QVERIFY_EXCEPTION_THROWN(solver.run(),RError);
    QCOMPARE(readBytes(latestPath),acceptedLink);
    QVERIFY(!QFileInfo::exists(trialRecord));
    QVERIFY(model.findVariable(R_VARIABLE_TEMPERATURE) == RConstants::eod);
    QVERIFY(model.findVariable(R_VARIABLE_DISPLACEMENT) == RConstants::eod);
    RModel persistedBeforeRetry;
    persistedBeforeRetry.read(latestPath);
    QCOMPARE(persistedBeforeRetry.getName(),QStringLiteral("accepted-before-heat-stress"));
    QVERIFY(persistedBeforeRetry.findVariable(R_VARIABLE_TEMPERATURE) == RConstants::eod);

    // Repair the model in place and retry on the same RSolver. Restored task
    // counters must make both sub-solvers run as a fresh coupled execution.
    RMaterial material = model.getLine(0).getMaterial();
    setMaterialProperty(material,RMaterialProperty::ModulusOfElasiticity,modulus);
    model.getLine(0).setMaterial(material);
    solver.run();

    QVERIFY(QFileInfo::exists(trialRecord));
    RModel persistedAfterRetry;
    persistedAfterRetry.read(latestPath);
    QVERIFY(persistedAfterRetry.findVariable(R_VARIABLE_TEMPERATURE) != RConstants::eod);
    QVERIFY(persistedAfterRetry.findVariable(R_VARIABLE_DISPLACEMENT) != RConstants::eod);
    QVERIFY(persistedAfterRetry.findVariable(R_VARIABLE_STRESS_X) != RConstants::eod);
}

void TestStructuralHeatStressCoupling::resultWriteFailureRollsBackAndCanRetry()
{
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString latestPath = directory.filePath("coupled.ram");
    seedAcceptedArtifact(latestPath);
    const QByteArray acceptedLink = readBytes(latestPath);
    QVERIFY(!acceptedLink.isEmpty());

    const QString trialRecord = RFileUtils::getFileNameWithTimeStep(latestPath,0);
    QVERIFY(QDir().mkpath(trialRecord));
    RModel model = heatedRod(false);
    RSolver solver(model,latestPath,QString());

    // Both solves succeed, but an existing directory blocks the final atomic
    // model-record write. The accepted link and pre-run model state survive.
    QVERIFY_EXCEPTION_THROWN(solver.run(),RError);
    QCOMPARE(readBytes(latestPath),acceptedLink);
    QVERIFY(QFileInfo(trialRecord).isDir());
    QVERIFY(model.findVariable(R_VARIABLE_TEMPERATURE) == RConstants::eod);
    QVERIFY(model.findVariable(R_VARIABLE_DISPLACEMENT) == RConstants::eod);
    RModel persistedBeforeRetry;
    persistedBeforeRetry.read(latestPath);
    QCOMPARE(persistedBeforeRetry.getName(),QStringLiteral("accepted-before-heat-stress"));

    QVERIFY(QDir(trialRecord).removeRecursively());
    solver.run();
    QVERIFY(QFileInfo::exists(trialRecord));
    RModel persistedAfterRetry;
    persistedAfterRetry.read(latestPath);
    QVERIFY(persistedAfterRetry.findVariable(R_VARIABLE_TEMPERATURE) != RConstants::eod);
    QVERIFY(persistedAfterRetry.findVariable(R_VARIABLE_STRESS_X) != RConstants::eod);
}

void TestStructuralHeatStressCoupling::transientSequenceAdvancesHeatThenStressAtEachTime()
{
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString latestPath = directory.filePath("transient-coupled.ram");
    // Keep the geometry fixed here so this time-marching check isolates the
    // coupled transaction and the exact transient Heat FE recurrence.
    RModel model = heatedRod(false,false,true,0.0);
    RTimeSolver &time = model.getTimeSolver();
    constexpr double dt = 500000.0;
    time.setEnabled(true);
    time.setInputNTimeSteps(2);
    time.setInputStartTime(0.0);
    time.setInputTimeStepSize(dt);
    time.setOutputFrequency(1);
    time.setTimeMarchApproximation(R_TIME_MARCH_BACKWARD);
    RSolver solver(model,latestPath,QString());
    solver.run();

    const double massDiagonal = 40000.0/3.0*2.0;
    const double centerConductance = 2.0*2.0*0.01;
    const double boundaryConductance = 2.0*0.01;
    auto nextCenterTemperature = [&](double previous) {
        return (massDiagonal*previous + dt*boundaryConductance*(300.0+400.0)) /
               (massDiagonal + dt*centerConductance);
    };
    const double firstCenter = nextCenterTemperature(referenceKelvin);
    const double secondCenter = nextCenterTemperature(firstCenter);
    QVERIFY(std::fabs(nodal(model,R_VARIABLE_TEMPERATURE,0,1)-secondCenter) < 1.0e-7);
    QVERIFY(std::fabs(nodal(model,R_VARIABLE_DISPLACEMENT,0,2)) < 1.0e-12);
    QVERIFY(std::fabs(elemental(model,R_VARIABLE_STRESS_X,0,0)) < 1.0e-8);
    QVERIFY(std::fabs(elemental(model,R_VARIABLE_STRESS_X,0,1)) < 1.0e-8);
    QCOMPARE(time.getCurrentTimeStep(),uint(1));
    QCOMPARE(time.getComputedTime(),2.0*dt);

    const QString firstRecord = RFileUtils::getFileNameWithTimeStep(latestPath,1);
    const QString secondRecord = RFileUtils::getFileNameWithTimeStep(latestPath,2);
    QVERIFY(QFileInfo::exists(firstRecord));
    QVERIFY(QFileInfo::exists(secondRecord));
    RModel persistedFirst;
    persistedFirst.read(firstRecord);
    const double persistedFirstCenter = nodal(persistedFirst,R_VARIABLE_TEMPERATURE,0,1);
    QVERIFY2(std::fabs(persistedFirstCenter-firstCenter) < 1.0e-3,
             qPrintable(QString("First stored coupled temperature = %1 K, expected %2 K")
                        .arg(persistedFirstCenter,0,'g',16).arg(firstCenter,0,'g',16)));
    QVERIFY(std::fabs(nodal(persistedFirst,R_VARIABLE_DISPLACEMENT,0,2)) < 1.0e-12);
    QVERIFY(std::fabs(elemental(persistedFirst,R_VARIABLE_STRESS_X,0,0)) < 1.0e-8);
    RModel persistedLast;
    persistedLast.read(latestPath);
    QVERIFY(std::fabs(nodal(persistedLast,R_VARIABLE_TEMPERATURE,0,1)-secondCenter) < 1.0e-3);
    QVERIFY(std::fabs(nodal(persistedLast,R_VARIABLE_DISPLACEMENT,0,2)) < 1.0e-12);
    QVERIFY(std::fabs(elemental(persistedLast,R_VARIABLE_STRESS_X,0,0)) < 1.0e-8);
}

void TestStructuralHeatStressCoupling::transientStressFailureRollsBackAndPreservesRecordsForRetry()
{
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString latestPath = directory.filePath("transient-coupled.ram");
    seedAcceptedArtifact(latestPath);
    const QByteArray acceptedLink = readBytes(latestPath);
    const QString acceptedRecord = RFileUtils::getFileNameWithTimeStep(latestPath,1);
    const QByteArray acceptedRecordBytes = readBytes(acceptedRecord);
    const QString futureRecord = RFileUtils::getFileNameWithTimeStep(latestPath,2);
    QVERIFY(!acceptedLink.isEmpty());
    QVERIFY(!acceptedRecordBytes.isEmpty());

    RModel model = heatedRod(false,false,false,0.0,false);
    RBoundaryCondition load(R_BOUNDARY_CONDITION_FORCE);
    setComponent(load,R_VARIABLE_FORCE_X,1000.0);
    model.getPoint(2).addBoundaryCondition(load);
    RTimeSolver &time = model.getTimeSolver();
    constexpr double dt = 500000.0;
    time.setEnabled(true);
    time.setInputNTimeSteps(2);
    time.setInputStartTime(0.0);
    time.setInputTimeStepSize(dt);
    time.setOutputFrequency(1);
    RSolver solver(model,latestPath,QString());

    QVERIFY_EXCEPTION_THROWN(solver.run(),RError);
    QCOMPARE(readBytes(latestPath),acceptedLink);
    QCOMPARE(readBytes(acceptedRecord),acceptedRecordBytes);
    QCOMPARE(time.getCurrentTimeStep(),uint(0));
    QCOMPARE(time.getComputedTime(),0.0);
    QVERIFY(model.findVariable(R_VARIABLE_TEMPERATURE) == RConstants::eod);
    QVERIFY(model.findVariable(R_VARIABLE_DISPLACEMENT) == RConstants::eod);

    RMaterial material = model.getLine(0).getMaterial();
    setMaterialProperty(material,RMaterialProperty::ModulusOfElasiticity,modulus);
    model.getLine(0).setMaterial(material);
    solver.run();
    QVERIFY(QFileInfo::exists(acceptedRecord));
    QVERIFY(QFileInfo::exists(futureRecord));
    RModel persisted;
    persisted.read(latestPath);
    QVERIFY(persisted.findVariable(R_VARIABLE_TEMPERATURE) != RConstants::eod);
    QVERIFY(persisted.findVariable(R_VARIABLE_STRESS_X) != RConstants::eod);
    QCOMPARE(persisted.getTimeSolver().getComputedTime(),2.0*dt);
}

void TestStructuralHeatStressCoupling::transientPublicationFailureRetriesSameTimeIndex()
{
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString latestPath = directory.filePath("transient-coupled.ram");
    RModel model = heatedRod(false,false,true,0.0);
    RTimeSolver &time = model.getTimeSolver();
    constexpr double dt = 500000.0;
    time.setEnabled(true);
    time.setInputNTimeSteps(1);
    time.setInputStartTime(0.0);
    time.setInputTimeStepSize(dt);
    time.setOutputFrequency(1);
    time.setTimeMarchApproximation(R_TIME_MARCH_BACKWARD);
    const QString futureRecord = RFileUtils::getFileNameWithTimeStep(latestPath,2);
    QFile futureFile(futureRecord);
    QVERIFY(futureFile.open(QIODevice::WriteOnly));
    const QByteArray futureSentinel("future-step-sentinel");
    QCOMPARE(futureFile.write(futureSentinel),qint64(futureSentinel.size()));
    futureFile.close();
    RSolver solver(model,latestPath,QString());
    solver.run();

    const double firstCenter =
        (40000.0/3.0*2.0*referenceKelvin + dt*2.0*0.01*(300.0+400.0)) /
        (40000.0/3.0*2.0 + dt*2.0*2.0*0.01);
    const QString firstRecord = RFileUtils::getFileNameWithTimeStep(latestPath,1);
    const QByteArray acceptedRecord = readBytes(firstRecord);
    const QByteArray acceptedLink = readBytes(latestPath);
    QVERIFY(!acceptedRecord.isEmpty());
    QVERIFY(!acceptedLink.isEmpty());
    QVERIFY(std::fabs(nodal(model,R_VARIABLE_TEMPERATURE,0,1)-firstCenter) < 1.0e-7);

    // Force RModel::write to commit the second-step record and then fail when
    // publishing the latest link. The transaction must restore the prior
    // future record and keep the accepted first-step model/time in memory.
    QVERIFY(QFile::remove(latestPath));
    QVERIFY(QDir().mkpath(latestPath));
    model.getProblemSetup().setRestart(true);
    QVERIFY_EXCEPTION_THROWN(solver.run(),RError);
    QCOMPARE(readBytes(futureRecord),futureSentinel);
    QVERIFY(QDir(latestPath).exists());
    QCOMPARE(readBytes(firstRecord),acceptedRecord);
    QVERIFY(std::fabs(nodal(model,R_VARIABLE_TEMPERATURE,0,1)-firstCenter) < 1.0e-7);
    QCOMPARE(time.getCurrentTimeStep(),uint(1));
    QCOMPARE(time.getComputedTime(),dt);

    QVERIFY(QDir(latestPath).removeRecursively());
    QFile restoredLink(latestPath);
    QVERIFY(restoredLink.open(QIODevice::WriteOnly));
    QCOMPARE(restoredLink.write(acceptedLink),qint64(acceptedLink.size()));
    restoredLink.close();
    solver.run();
    QCOMPARE(time.getCurrentTimeStep(),uint(1));
    QCOMPARE(time.getComputedTime(),2.0*dt);
    QCOMPARE(readBytes(firstRecord),acceptedRecord);
    QVERIFY(QFileInfo(futureRecord).isFile());
    QVERIFY(readBytes(futureRecord) != futureSentinel);
    RModel persisted;
    persisted.read(latestPath);
    QVERIFY(std::fabs(nodal(persisted,R_VARIABLE_TEMPERATURE,0,1)-
                      (40000.0/3.0*2.0*firstCenter + dt*2.0*0.01*(300.0+400.0)) /
                      (40000.0/3.0*2.0 + dt*2.0*2.0*0.01)) < 1.0e-3);
}

void TestStructuralHeatStressCoupling::transientRestartRejectsUnpublishedFutureRecords()
{
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString latestPath = directory.filePath("transient-restart.ram");
    constexpr double dt = 500000.0;

    RModel accepted;
    accepted.setName(QStringLiteral("accepted-restart-state"));
    RTimeSolver &acceptedTime = accepted.getTimeSolver();
    acceptedTime.setEnabled(true);
    acceptedTime.setTimes({dt,2.0*dt});
    acceptedTime.setCurrentTimeStep(1);
    acceptedTime.setComputedTime(2.0*dt);
    accepted.write(latestPath);
    const QByteArray acceptedLink = readBytes(latestPath);
    const QString acceptedRecord = RFileUtils::getFileNameWithTimeStep(latestPath,2);
    const QByteArray acceptedRecordBytes = readBytes(acceptedRecord);
    QVERIFY(!acceptedLink.isEmpty());
    QVERIFY(!acceptedRecordBytes.isEmpty());

    // Step three is the immediate next publication and is transactionally
    // replaceable. The later in-schedule records cannot be safely cleaned
    // until the run finishes, so restart must fail before starting a trial.
    const QString immediateRecord = RFileUtils::getFileNameWithTimeStep(latestPath,3);
    const QString laterRecord = RFileUtils::getFileNameWithTimeStep(latestPath,4);
    const QString finalRecord = RFileUtils::getFileNameWithTimeStep(latestPath,5);
    const QByteArray immediateSentinel("immediate-trial-sentinel");
    const QByteArray laterSentinel("stale-future-four");
    const QByteArray finalSentinel("stale-future-five");
    for (const auto &entry : {qMakePair(immediateRecord,immediateSentinel),
                              qMakePair(laterRecord,laterSentinel),
                              qMakePair(finalRecord,finalSentinel)})
    {
        QFile file(entry.first);
        QVERIFY(file.open(QIODevice::WriteOnly));
        QCOMPARE(file.write(entry.second),qint64(entry.second.size()));
    }

    RModel model = heatedRod(false,false,true,0.0);
    RTimeSolver &time = model.getTimeSolver();
    time.setEnabled(true);
    time.setTimes({dt,2.0*dt});
    time.setCurrentTimeStep(1);
    time.setComputedTime(2.0*dt);
    time.setInputNTimeSteps(3);
    time.setInputStartTime(2.0*dt);
    time.setInputTimeStepSize(dt);
    time.setOutputFrequency(1);
    model.getProblemSetup().setRestart(true);
    RSolver solver(model,latestPath,QString());

    QVERIFY_EXCEPTION_THROWN(solver.run(),RError);
    QCOMPARE(readBytes(latestPath),acceptedLink);
    QCOMPARE(readBytes(acceptedRecord),acceptedRecordBytes);
    QCOMPARE(readBytes(immediateRecord),immediateSentinel);
    QCOMPARE(readBytes(laterRecord),laterSentinel);
    QCOMPARE(readBytes(finalRecord),finalSentinel);
    QCOMPARE(time.getCurrentTimeStep(),uint(1));
    QCOMPARE(time.getComputedTime(),2.0*dt);
}

void TestStructuralHeatStressCoupling::transientRestartFindsSixDigitOutputRecord()
{
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString latestPath = directory.filePath("transient-six-digit.ram");
    constexpr double dt = 500000.0;

    RModel accepted;
    accepted.getTimeSolver().setEnabled(true);
    accepted.getTimeSolver().setTimes({dt});
    accepted.write(latestPath);
    const QByteArray acceptedLink = readBytes(latestPath);
    const QString acceptedRecord = RFileUtils::getFileNameWithTimeStep(latestPath,1);
    const QByteArray acceptedRecordBytes = readBytes(acceptedRecord);
    QVERIFY(!acceptedLink.isEmpty());
    QVERIFY(!acceptedRecordBytes.isEmpty());

    const QString sixDigitRecord = RFileUtils::getFileNameWithTimeStep(latestPath,100000);
    QFile staleFile(sixDigitRecord);
    QVERIFY(staleFile.open(QIODevice::WriteOnly));
    const QByteArray staleSentinel("stale-step-100000");
    QCOMPARE(staleFile.write(staleSentinel),qint64(staleSentinel.size()));
    staleFile.close();

    RModel model = heatedRod(false,false,true,0.0);
    RTimeSolver &time = model.getTimeSolver();
    time.setEnabled(true);
    time.setTimes({dt});
    time.setCurrentTimeStep(0);
    time.setComputedTime(dt);
    time.setInputNTimeSteps(99999);
    time.setInputStartTime(dt);
    time.setInputTimeStepSize(dt);
    time.setOutputFrequency(1);
    model.getProblemSetup().setRestart(true);
    RSolver solver(model,latestPath,QString());

    QVERIFY_EXCEPTION_THROWN(solver.run(),RError);
    QCOMPARE(readBytes(latestPath),acceptedLink);
    QCOMPARE(readBytes(acceptedRecord),acceptedRecordBytes);
    QCOMPARE(readBytes(sixDigitRecord),staleSentinel);
    QCOMPARE(time.getCurrentTimeStep(),uint(0));
    QCOMPARE(time.getComputedTime(),dt);
}

void TestStructuralHeatStressCoupling::transientOutputLockRejectsConcurrentWriter()
{
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString latestPath = directory.filePath("transient-locked.ram");
    RModel model = heatedRod(false,false,true,0.0);
    RTimeSolver &time = model.getTimeSolver();
    time.setEnabled(true);
    time.setInputNTimeSteps(2);
    time.setInputStartTime(0.0);
    time.setInputTimeStepSize(500000.0);
    time.setOutputFrequency(1);
    const QString lockPath = latestPath+QStringLiteral(".transient.lock");
    QLockFile competingWriter(lockPath);
    competingWriter.setStaleLockTime(0);
    QVERIFY(competingWriter.tryLock(0));
    RSolver solver(model,latestPath,QString());

    QVERIFY_EXCEPTION_THROWN(solver.run(),RError);
    QCOMPARE(time.getCurrentTimeStep(),uint(0));
    QCOMPARE(time.getComputedTime(),0.0);
    QVERIFY(!QFileInfo::exists(RFileUtils::getFileNameWithTimeStep(latestPath,1)));
    QVERIFY(!QFileInfo::exists(latestPath));
    competingWriter.unlock();
    QVERIFY(!QFileInfo::exists(lockPath));
}

void TestStructuralHeatStressCoupling::transientRejectsMonitoringOutputBeforeStep()
{
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString latestPath = directory.filePath("transient-monitor.ram");
    const QString monitorPath = directory.filePath("monitor.txt");
    const QByteArray monitorSentinel("accepted-monitor-data");
    QFile monitorFile(monitorPath);
    QVERIFY(monitorFile.open(QIODevice::WriteOnly));
    QCOMPARE(monitorFile.write(monitorSentinel),qint64(monitorSentinel.size()));
    monitorFile.close();

    RModel model = heatedRod(false,false,true,0.0);
    RMonitoringPoint point(R_VARIABLE_TEMPERATURE,RR3Vector(1.0,0.0,0.0));
    model.getMonitoringPointManager().push_back(point);
    model.getMonitoringPointManager().setOutputFileName(monitorPath);
    RTimeSolver &time = model.getTimeSolver();
    time.setEnabled(true);
    time.setInputNTimeSteps(2);
    time.setInputStartTime(0.0);
    time.setInputTimeStepSize(500000.0);
    time.setOutputFrequency(1);
    RSolver solver(model,latestPath,QString());

    QVERIFY_EXCEPTION_THROWN(solver.run(),RError);
    QCOMPARE(readBytes(monitorPath),monitorSentinel);
    QVERIFY(!QFileInfo::exists(latestPath));
    QVERIFY(!QFileInfo::exists(RFileUtils::getFileNameWithTimeStep(latestPath,1)));
    QCOMPARE(time.getCurrentTimeStep(),uint(0));
    QCOMPARE(time.getComputedTime(),0.0);
    QVERIFY(model.findVariable(R_VARIABLE_TEMPERATURE) == RConstants::eod);
}

void TestStructuralHeatStressCoupling::unsupportedTransientSequenceIsRejectedBeforeOutputCleanup()
{
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString latestPath = directory.filePath("coupled.ram");
    seedAcceptedArtifact(latestPath);
    const QByteArray acceptedLink = readBytes(latestPath);
    const QString acceptedRecord = RFileUtils::getFileNameWithTimeStep(latestPath,1);
    const QByteArray acceptedRecordBytes = readBytes(acceptedRecord);
    QVERIFY(!acceptedLink.isEmpty());
    QVERIFY(!acceptedRecordBytes.isEmpty());

    RModel model = heatedRod(false);
    model.getTimeSolver().setEnabled(true);
    RProblemTaskItem unsupported = model.getProblemTaskTree();
    unsupported.setNIterations(2);
    model.setProblemTaskTree(unsupported);
    RSolver solver(model,latestPath,QString());
    QVERIFY_EXCEPTION_THROWN(solver.run(),RError);

    // The unsupported transient shape is rejected before the legacy old-record
    // cleanup can rewrite or remove an accepted artifact.
    QCOMPARE(readBytes(latestPath),acceptedLink);
    QCOMPARE(readBytes(acceptedRecord),acceptedRecordBytes);
}

QTEST_APPLESS_MAIN(TestStructuralHeatStressCoupling)
#include "tst_structural_heat_stress_coupling.moc"

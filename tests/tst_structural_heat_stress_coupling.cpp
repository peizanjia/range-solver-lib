#include <QtTest>

#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QTemporaryDir>

#include <cmath>
#include <limits>

#include <rbl_error.h>
#include <rml_file_utils.h>
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
                 bool elasticModulusDefined = true)
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
        setMaterialProperty(material,RMaterialProperty::ThermalExpansionCoefficient,expansion);
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
    void transientSequenceIsRejectedBeforeOutputCleanup();
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

void TestStructuralHeatStressCoupling::transientSequenceIsRejectedBeforeOutputCleanup()
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
    model.getTimeSolver().setTimes({0.01,0.02});
    RSolver solver(model,latestPath,QString());
    QVERIFY_EXCEPTION_THROWN(solver.run(),RError);

    // The unsupported transient shape is rejected before the legacy old-record
    // cleanup can rewrite or remove an accepted artifact.
    QCOMPARE(readBytes(latestPath),acceptedLink);
    QCOMPARE(readBytes(acceptedRecord),acceptedRecordBytes);
}

QTEST_APPLESS_MAIN(TestStructuralHeatStressCoupling)
#include "tst_structural_heat_stress_coupling.moc"

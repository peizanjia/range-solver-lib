#include <QtTest>

#include <QTemporaryDir>

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <utility>
#include <vector>

#include "rml_model.h"
#include "rsolver.h"

namespace
{
constexpr double length = 2.0;
constexpr double area = 0.01;
constexpr double youngModulus = 2.0e11;
constexpr double density = 8000.0;
constexpr double initialDisplacement = 1.0e-4;
constexpr double stepSize = 1.0e-5;
constexpr uint firstSegmentSteps = 8;
constexpr uint secondSegmentSteps = 12;

RModel makeTruss(uint steps)
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

    // The axial truss has no transverse stiffness; remove those DOFs.
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
    tipPoint.setNodeId(0,1);
    model.addElement(tipPoint,true,1);
    RInitialCondition initial(R_INITIAL_CONDITION_DISPLACEMENT);
    initial.getComponent(initial.findComponentPosition(R_VARIABLE_DISPLACEMENT_X))
           .add(0.0,initialDisplacement);
    model.getPoint(1).addInitialCondition(initial);

    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    auto &solver = model.getMatrixSolverConf(RMatrixSolverConf::CG);
    solver.setOutputFrequency(0);
    solver.setSolverCvgValue(1.0e-13);
    solver.setNInnerIterations(100);
    solver.setNOuterIterations(1000);
    auto &dynamics = model.getProblemSetup().getStructuralDynamicsSetup();
    dynamics.setEnabled(true);
    dynamics.setMethod(RStructuralTimeIntegrationMethod::Newmark);
    dynamics.setTimeStep(stepSize);
    dynamics.setNSteps(steps);
    dynamics.setOutputFrequency(steps);
    dynamics.setRayleighAlpha(100.0); // Nonzero dissipated energy checks bookkeeping.
    return model;
}

void addMovingSupport(RModel &model)
{
    RStructuralPrescribedMotion motion;
    motion.id = "restart-base-motion";
    motion.target = {0,RDofType::UX};
    const double endTime = (firstSegmentSteps+secondSegmentSteps)*stepSize;
    // Nontrivial 17-digit data detects lossy ASCII serialization and a
    // fingerprint that rounds control trajectories like material operators.
    motion.knots = {{0.0,0.0,0.031415926535897932,1012.345678901234},
                    {endTime,0.000023456789012345678,0.11728182845904524,
                     -345.6789012345678}};
    model.getProblemSetup().getStructuralInteractionSetup()
         .getPrescribedMotions().push_back(motion);
}

bool solve(RModel &model, QString &error)
{
    error.clear();
    try
    {
        RSolver solver(model,QString(),QString());
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

bool load(RModel &model, const QString &path, QString &error)
{
    error.clear();
    try
    {
        model.read(path);
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

double value(const RModel &model, RVariableType type)
{
    const uint index = model.findVariable(type);
    if (index == RConstants::eod)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return model.getVariable(index).getValue(0,1);
}

double mechanicalEnergy(const RModel &model)
{
    const double mass = density*area*length/3.0;
    const double stiffness = youngModulus*area/length;
    const double q = value(model,R_VARIABLE_DISPLACEMENT);
    const double v = value(model,R_VARIABLE_VELOCITY);
    return 0.5*mass*v*v + 0.5*stiffness*q*q;
}

bool close(double actual, double expected, double absoluteTolerance = 1.0e-12)
{
    const double scale = std::max({1.0,std::abs(actual),std::abs(expected)});
    return std::isfinite(actual) && std::isfinite(expected) &&
           std::abs(actual-expected) <= absoluteTolerance*scale;
}

bool closeAtPhysicalScale(double actual, double expected, double characteristicScale,
                          double relativeTolerance = 1.0e-8)
{
    const double scale = std::max({characteristicScale,std::abs(actual),std::abs(expected)});
    return std::isfinite(actual) && std::isfinite(expected) &&
           std::isfinite(scale) && scale > 0.0 &&
           std::abs(actual-expected) <= relativeTolerance*scale;
}
}

class TestStructuralDynamicsRestartIO : public QObject
{
    Q_OBJECT

private slots:
    void diskRestartMatchesUninterrupted_data();
    void diskRestartMatchesUninterrupted();
    void modalDampingDiskRestartMatchesUninterrupted_data();
    void modalDampingDiskRestartMatchesUninterrupted();
    void modalDampingRejectsChangedRatio();
    void movingSupportDiskRestartMatchesUninterrupted_data();
    void movingSupportDiskRestartMatchesUninterrupted();
    void movingSupportRejectsChangedKnot();
    void rejectsTamperedCheckpoint();
    void rejectsChangedMesh();
};

void TestStructuralDynamicsRestartIO::diskRestartMatchesUninterrupted_data()
{
    QTest::addColumn<bool>("binary");
    QTest::newRow("ascii-ram") << false;
    QTest::newRow("binary-rbm") << true;
}

void TestStructuralDynamicsRestartIO::diskRestartMatchesUninterrupted()
{
    QFETCH(bool,binary);
    QTemporaryDir directory;
    QVERIFY(directory.isValid());

    RModel uninterrupted = makeTruss(firstSegmentSteps+secondSegmentSteps);
    RModel firstSegment = makeTruss(firstSegmentSteps);
    QString error;
    QVERIFY2(solve(uninterrupted,error),qPrintable(error));
    QVERIFY2(solve(firstSegment,error),qPrintable(error));
    const auto saved = firstSegment.getProblemSetup().getStructuralDynamicsRestartData();
    QVERIFY(saved.available);
    QCOMPARE(saved.completedSteps,firstSegmentSteps);
    QVERIFY(saved.cumulativeDampingDissipation > 0.0);

    // The native model formats are .ram (ASCII) and .rbm (binary), not .rml.
    // No link file is used: this exercises writing and reading the actual data.
    const QString path = directory.filePath(binary ? "accepted.rbm" : "accepted.ram");
    QCOMPARE(firstSegment.write(path,false),path);
    RModel resumed;
    QVERIFY2(load(resumed,path,error),qPrintable(error));
    QCOMPARE(resumed.getNNodes(),firstSegment.getNNodes());
    QCOMPARE(resumed.getNElements(),firstSegment.getNElements());
    QCOMPARE(resumed.getNPoints(),firstSegment.getNPoints());
    QCOMPARE(resumed.getNLines(),firstSegment.getNLines());
    for (uint point=0;point<resumed.getNPoints();point++)
    {
        QCOMPARE(resumed.getPoint(point).getNBoundaryConditions(),
                 firstSegment.getPoint(point).getNBoundaryConditions());
        QCOMPARE(resumed.getPoint(point).getNInitialConditions(),
                 firstSegment.getPoint(point).getNInitialConditions());
    }
    QCOMPARE(resumed.getLine(0).getNBoundaryConditions(),
             firstSegment.getLine(0).getNBoundaryConditions());
    QCOMPARE(resumed.getLine(0).getCrossArea(),firstSegment.getLine(0).getCrossArea());
    QVERIFY(resumed.getLine(0).getMaterial().propertiesEqual(
        firstSegment.getLine(0).getMaterial()));
    const auto &readBack = resumed.getProblemSetup().getStructuralDynamicsRestartData();
    QVERIFY(readBack.available);
    QCOMPARE(readBack.completedSteps,firstSegmentSteps);
    QCOMPARE(readBack.meshFingerprint,saved.meshFingerprint);
    QCOMPARE(readBack.dofFingerprint,saved.dofFingerprint);
    QCOMPARE(readBack.systemFingerprint,saved.systemFingerprint);
    QCOMPARE(readBack.q.size(),saved.q.size());
    for (std::size_t i=0;i<saved.q.size();++i)
    {
        QVERIFY(close(readBack.q[i],saved.q[i]));
        QVERIFY(close(readBack.v[i],saved.v[i]));
        QVERIFY(close(readBack.a[i],saved.a[i]));
    }

    resumed.getProblemSetup().setRestart(true);
    resumed.getProblemSetup().getStructuralDynamicsSetup().setNSteps(secondSegmentSteps);
    resumed.getProblemSetup().getStructuralDynamicsSetup().setOutputFrequency(secondSegmentSteps);
    QVERIFY2(solve(resumed,error),qPrintable(error));
    const auto &expected = uninterrupted.getProblemSetup().getStructuralDynamicsRestartData();
    const auto &actual = resumed.getProblemSetup().getStructuralDynamicsRestartData();
    QCOMPARE(actual.completedSteps,firstSegmentSteps+secondSegmentSteps);
    QVERIFY(close(actual.time,expected.time));
    QCOMPARE(actual.q.size(),expected.q.size());
    for (std::size_t i=0;i<expected.q.size();++i)
    {
        QVERIFY(close(actual.q[i],expected.q[i],1.0e-10));
        QVERIFY(close(actual.v[i],expected.v[i],1.0e-10));
        QVERIFY(close(actual.a[i],expected.a[i],1.0e-10));
        QVERIFY(close(actual.externalForce[i],expected.externalForce[i],1.0e-10));
    }
    QVERIFY(close(actual.cumulativeExternalWork,expected.cumulativeExternalWork,1.0e-10));
    QVERIFY(close(actual.cumulativeDampingDissipation,
                  expected.cumulativeDampingDissipation,1.0e-10));
    QVERIFY(close(actual.initialMechanicalEnergy,expected.initialMechanicalEnergy,1.0e-10));
    for (RVariableType type : {R_VARIABLE_DISPLACEMENT,
                               R_VARIABLE_VELOCITY,
                               R_VARIABLE_ACCELERATION})
    {
        QVERIFY2(close(value(resumed,type),value(uninterrupted,type),1.0e-10),
                 qPrintable(QString("Disk restart mismatch for variable %1").arg(int(type))));
    }
    QVERIFY(close(mechanicalEnergy(resumed),mechanicalEnergy(uninterrupted),1.0e-10));
}

void TestStructuralDynamicsRestartIO::modalDampingDiskRestartMatchesUninterrupted_data()
{
    QTest::addColumn<bool>("binary");
    QTest::newRow("ascii-ram") << false;
    QTest::newRow("binary-rbm") << true;
}

void TestStructuralDynamicsRestartIO::modalDampingDiskRestartMatchesUninterrupted()
{
    QFETCH(bool,binary);
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    constexpr double dampingRatio = 0.05000000000000001;
    auto enableModal = [dampingRatio](RModel &model)
    {
        auto &setup = model.getProblemSetup().getStructuralDynamicsSetup();
        setup.setRayleighAlpha(0.0);
        setup.setModalDampingRatios({dampingRatio});
    };
    RModel uninterrupted = makeTruss(firstSegmentSteps+secondSegmentSteps);
    RModel firstSegment = makeTruss(firstSegmentSteps);
    enableModal(uninterrupted);
    enableModal(firstSegment);
    QString error;
    QVERIFY2(solve(uninterrupted,error),qPrintable(error));
    QVERIFY2(solve(firstSegment,error),qPrintable(error));
    const QString path = directory.filePath(binary ? "modal-accepted.rbm" : "modal-accepted.ram");
    QCOMPARE(firstSegment.write(path,false),path);
    RModel resumed;
    QVERIFY2(load(resumed,path,error),qPrintable(error));
    const auto &loadedRatios = resumed.getProblemSetup()
        .getStructuralDynamicsSetup().getModalDampingRatios();
    QCOMPARE(loadedRatios.size(),std::size_t(1));
    QCOMPARE(loadedRatios.front(),dampingRatio);
    resumed.getProblemSetup().setRestart(true);
    auto &resumedSetup = resumed.getProblemSetup().getStructuralDynamicsSetup();
    resumedSetup.setNSteps(secondSegmentSteps);
    resumedSetup.setOutputFrequency(secondSegmentSteps);
    QVERIFY2(solve(resumed,error),qPrintable(error));

    const auto &expected = uninterrupted.getProblemSetup().getStructuralDynamicsRestartData();
    const auto &actual = resumed.getProblemSetup().getStructuralDynamicsRestartData();
    QCOMPARE(actual.completedSteps,firstSegmentSteps+secondSegmentSteps);
    QCOMPARE(actual.systemFingerprint,expected.systemFingerprint);
    const double omega = std::sqrt(3.0*youngModulus/(density*length*length));
    const double velocityScale = initialDisplacement*omega;
    const double accelerationScale = initialDisplacement*omega*omega;
    const double forceScale = youngModulus*area/length*initialDisplacement;
    const double energyScale = 0.5*forceScale*initialDisplacement;
    QVERIFY(closeAtPhysicalScale(actual.time,expected.time,stepSize,1.0e-10));
    QVERIFY(closeAtPhysicalScale(actual.fixedTimeStep,expected.fixedTimeStep,stepSize,1.0e-10));
    QCOMPARE(actual.q.size(),expected.q.size());
    QCOMPARE(actual.v.size(),expected.v.size());
    QCOMPARE(actual.a.size(),expected.a.size());
    QCOMPARE(actual.externalForce.size(),expected.externalForce.size());
    for (std::size_t i=0;i<expected.q.size();++i)
    {
        QVERIFY(closeAtPhysicalScale(actual.q[i],expected.q[i],initialDisplacement));
        QVERIFY(closeAtPhysicalScale(actual.v[i],expected.v[i],velocityScale));
        QVERIFY(closeAtPhysicalScale(actual.a[i],expected.a[i],accelerationScale));
        QVERIFY(closeAtPhysicalScale(actual.externalForce[i],expected.externalForce[i],
                                     forceScale,1.0e-12));
    }
    QVERIFY(closeAtPhysicalScale(actual.cumulativeExternalWork,
                                 expected.cumulativeExternalWork,energyScale));
    QVERIFY(closeAtPhysicalScale(actual.cumulativeDampingDissipation,
                                 expected.cumulativeDampingDissipation,energyScale));
    QVERIFY(closeAtPhysicalScale(actual.initialMechanicalEnergy,
                                 expected.initialMechanicalEnergy,energyScale));
    for (const auto &[type,scale] : {std::pair<RVariableType,double>{R_VARIABLE_DISPLACEMENT,initialDisplacement},
                                     {R_VARIABLE_VELOCITY,velocityScale},
                                     {R_VARIABLE_ACCELERATION,accelerationScale}})
    {
        QVERIFY2(closeAtPhysicalScale(value(resumed,type),value(uninterrupted,type),scale),
                 qPrintable(QString("Modal disk restart mismatch for variable %1").arg(int(type))));
    }
}

void TestStructuralDynamicsRestartIO::modalDampingRejectsChangedRatio()
{
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    RModel firstSegment = makeTruss(firstSegmentSteps);
    auto &setup = firstSegment.getProblemSetup().getStructuralDynamicsSetup();
    setup.setRayleighAlpha(0.0);
    setup.setModalDampingRatios({0.05000000000000001});
    QString error;
    QVERIFY2(solve(firstSegment,error),qPrintable(error));
    const QString path = directory.filePath("modal-checkpoint.ram");
    QCOMPARE(firstSegment.write(path,false),path);
    RModel altered;
    QVERIFY2(load(altered,path,error),qPrintable(error));
    altered.getProblemSetup().getStructuralDynamicsSetup()
           .setModalDampingRatios({0.05000000000000002});
    altered.getProblemSetup().setRestart(true);
    auto &alteredSetup = altered.getProblemSetup().getStructuralDynamicsSetup();
    alteredSetup.setNSteps(secondSegmentSteps);
    QVERIFY2(!solve(altered,error),"Changed modal damping ratio was accepted by restart.");
    QVERIFY2(error.contains("checkpoint",Qt::CaseInsensitive),qPrintable(error));
}

void TestStructuralDynamicsRestartIO::movingSupportDiskRestartMatchesUninterrupted_data()
{
    QTest::addColumn<bool>("binary");
    QTest::newRow("ascii-ram") << false;
    QTest::newRow("binary-rbm") << true;
}

void TestStructuralDynamicsRestartIO::movingSupportDiskRestartMatchesUninterrupted()
{
    QFETCH(bool,binary);
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    RModel uninterrupted = makeTruss(firstSegmentSteps+secondSegmentSteps);
    RModel firstSegment = makeTruss(firstSegmentSteps);
    addMovingSupport(uninterrupted);
    addMovingSupport(firstSegment);
    QString error;
    QVERIFY2(solve(uninterrupted,error),qPrintable(error));
    QVERIFY2(solve(firstSegment,error),qPrintable(error));
    const auto &savedMotion = firstSegment.getProblemSetup()
        .getStructuralInteractionSetup().getPrescribedMotions().front();
    const QString path = directory.filePath(binary ? "motion-accepted.rbm" : "motion-accepted.ram");
    QCOMPARE(firstSegment.write(path,false),path);
    RModel resumed;
    QVERIFY2(load(resumed,path,error),qPrintable(error));
    const auto &loadedMotions = resumed.getProblemSetup()
        .getStructuralInteractionSetup().getPrescribedMotions();
    QCOMPARE(loadedMotions.size(),std::size_t(1));
    QCOMPARE(loadedMotions.front().knots.size(),savedMotion.knots.size());
    QCOMPARE(loadedMotions.front().knots.back().acceleration,
             savedMotion.knots.back().acceleration);
    QCOMPARE(loadedMotions.front().knots.front().velocity,
             savedMotion.knots.front().velocity);

    resumed.getProblemSetup().setRestart(true);
    resumed.getProblemSetup().getStructuralDynamicsSetup().setNSteps(secondSegmentSteps);
    resumed.getProblemSetup().getStructuralDynamicsSetup().setOutputFrequency(secondSegmentSteps);
    QVERIFY2(solve(resumed,error),qPrintable(error));
    const auto &expected = uninterrupted.getProblemSetup().getStructuralDynamicsRestartData();
    const auto &actual = resumed.getProblemSetup().getStructuralDynamicsRestartData();
    QCOMPARE(actual.completedSteps,firstSegmentSteps+secondSegmentSteps);
    QCOMPARE(actual.systemFingerprint,expected.systemFingerprint);
    for (std::size_t i=0;i<expected.q.size();++i)
    {
        QVERIFY(close(actual.q[i],expected.q[i],1.0e-11));
        QVERIFY(close(actual.v[i],expected.v[i],1.0e-11));
        QVERIFY(close(actual.a[i],expected.a[i],1.0e-8));
        QVERIFY(close(actual.externalForce[i],expected.externalForce[i],1.0e-8));
    }
    for (RVariableType type : {R_VARIABLE_DISPLACEMENT,R_VARIABLE_VELOCITY,
                               R_VARIABLE_ACCELERATION})
    {
        QVERIFY2(close(value(resumed,type),value(uninterrupted,type),1.0e-9),
                 qPrintable(QString("Moving support restart mismatch in %1").arg(int(type))));
        const uint index = resumed.findVariable(type);
        QVERIFY(index != RConstants::eod);
        QVERIFY(close(resumed.getVariable(index).getValue(0,0),
                      uninterrupted.getVariable(index).getValue(0,0),1.0e-9));
    }
}

void TestStructuralDynamicsRestartIO::movingSupportRejectsChangedKnot()
{
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    RModel firstSegment = makeTruss(firstSegmentSteps);
    addMovingSupport(firstSegment);
    QString error;
    QVERIFY2(solve(firstSegment,error),qPrintable(error));
    const QString path = directory.filePath("motion-checkpoint.ram");
    QCOMPARE(firstSegment.write(path,false),path);
    RModel altered;
    QVERIFY2(load(altered,path,error),qPrintable(error));
    altered.getProblemSetup().getStructuralInteractionSetup()
           .getPrescribedMotions().front().knots.back().acceleration += 1.0e-6;
    altered.getProblemSetup().setRestart(true);
    altered.getProblemSetup().getStructuralDynamicsSetup().setNSteps(secondSegmentSteps);
    QVERIFY2(!solve(altered,error),"A modified control knot was accepted by restart.");
    QVERIFY2(error.contains("checkpoint",Qt::CaseInsensitive),qPrintable(error));
}

void TestStructuralDynamicsRestartIO::rejectsTamperedCheckpoint()
{
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    RModel firstSegment = makeTruss(firstSegmentSteps);
    QString error;
    QVERIFY2(solve(firstSegment,error),qPrintable(error));
    const QString path = directory.filePath("accepted.rbm");
    QCOMPARE(firstSegment.write(path,false),path);
    RModel altered;
    QVERIFY2(load(altered,path,error),qPrintable(error));
    auto &checkpoint = altered.getProblemSetup().getStructuralDynamicsRestartData();
    QVERIFY(checkpoint.available);
    checkpoint.systemFingerprint[0] = checkpoint.systemFingerprint[0] == QLatin1Char('a')
                                    ? QLatin1Char('b') : QLatin1Char('a');
    altered.getProblemSetup().setRestart(true);
    altered.getProblemSetup().getStructuralDynamicsSetup().setNSteps(1);
    QVERIFY2(!solve(altered,error),"A changed checkpoint system fingerprint was accepted.");
    QVERIFY2(error.contains("checkpoint",Qt::CaseInsensitive),qPrintable(error));
}

void TestStructuralDynamicsRestartIO::rejectsChangedMesh()
{
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    RModel firstSegment = makeTruss(firstSegmentSteps);
    QString error;
    QVERIFY2(solve(firstSegment,error),qPrintable(error));
    const QString path = directory.filePath("accepted.rbm");
    QCOMPARE(firstSegment.write(path,false),path);
    RModel altered;
    QVERIFY2(load(altered,path,error),qPrintable(error));
    altered.getNode(1).set(length+0.1,0.0,0.0);
    altered.getProblemSetup().setRestart(true);
    altered.getProblemSetup().getStructuralDynamicsSetup().setNSteps(1);
    QVERIFY2(!solve(altered,error),"A changed mesh was accepted for checkpoint restart.");
    QVERIFY2(error.contains("checkpoint",Qt::CaseInsensitive),qPrintable(error));
}

QTEST_APPLESS_MAIN(TestStructuralDynamicsRestartIO)
#include "tst_structural_dynamics_restart_io.moc"

#include <QtTest>

#include <QTemporaryDir>

#include <algorithm>
#include <array>
#include <cmath>
#include <exception>
#include <limits>
#include <vector>

#include <rbl_error.h>

#include <rml_boundary_condition.h>
#include <rml_element_group.h>
#include <rml_file_utils.h>
#include <rml_line.h>
#include <rml_material.h>
#include <rml_model.h>
#include <rml_problem_task_item.h>
#include <rml_structural_cable_setup.h>
#include <rml_structural_dynamics_restart_data.h>
#include <rml_structural_dynamics_setup.h>
#include <rml_variable.h>

#include "rsolver.h"
#include "rstructural_cable_model_adapter.h"

namespace
{
constexpr double cableLength = 2.0;
constexpr double cableArea = 0.01;
constexpr double youngModulus = 20000.0;
constexpr double density = 10.0;
constexpr double prestress = 1.0;
constexpr double initialVelocity = 0.01;
constexpr double timeStep = 0.005;
constexpr uint totalSteps = 6;
constexpr uint firstSegmentSteps = 3;

void setComponent(RBoundaryCondition &condition,
                  RVariableType type,
                  bool enabled,
                  double value = 0.0)
{
    const uint position = condition.findComponentPosition(type);
    if (position == RConstants::eod)
    {
        throw std::logic_error("Boundary condition is missing a cable component.");
    }
    RConditionComponent &component = condition.getComponent(position);
    component.setEnabled(enabled);
    component.clear();
    if (enabled) component.add(0.0,value);
}

RBoundaryCondition displacementCondition(const std::array<bool,3> &enabled)
{
    RBoundaryCondition condition(R_BOUNDARY_CONDITION_DISPLACEMENT);
    const RVariableType types[3] = {
        R_VARIABLE_DISPLACEMENT_X,R_VARIABLE_DISPLACEMENT_Y,R_VARIABLE_DISPLACEMENT_Z
    };
    for (uint axis=0;axis<3;axis++) setComponent(condition,types[axis],enabled[axis]);
    return condition;
}

RBoundaryCondition forceCondition(double x)
{
    RBoundaryCondition condition(R_BOUNDARY_CONDITION_FORCE);
    setComponent(condition,R_VARIABLE_FORCE_X,true,x);
    setComponent(condition,R_VARIABLE_FORCE_Y,true,0.0);
    setComponent(condition,R_VARIABLE_FORCE_Z,true,0.0);
    return condition;
}

void addPoint(RModel &model, uint node)
{
    RElement point(R_ELEMENT_POINT);
    point.setNodeId(0,node);
    model.addElement(point,true,node);
}

void setZeroGravity(RModel &model)
{
    const RVariableType types[3] = {
        R_VARIABLE_G_ACCELERATION_X,R_VARIABLE_G_ACCELERATION_Y,R_VARIABLE_G_ACCELERATION_Z
    };
    for (uint groupId=0;groupId<model.getNElementGroups();groupId++)
    {
        RElementGroup *group = model.getElementGroupPtr(groupId);
        if (!group) continue;
        for (RVariableType type : types)
        {
            REnvironmentCondition &condition =
                group->getEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
            const uint position = condition.findComponentPosition(type);
            if (position == RConstants::eod)
            {
                throw std::logic_error("Environment condition is missing a gravity component.");
            }
            RConditionComponent &component = condition.getComponent(position);
            component.setEnabled(true);
            component.clear();
            component.add(0.0,0.0);
        }
    }
}

void setInitialVelocity(RModel &model, double value)
{
    const uint position = model.addVariable(R_VARIABLE_VELOCITY);
    RVariable &velocity = model.getVariable(position);
    velocity.setApplyType(R_VARIABLE_APPLY_NODE);
    velocity.resize(3,model.getNNodes());
    velocity.setValue(0,1,value);
}

RModel makeVibratingCable(uint steps,
                          double memberPrestress = prestress,
                          double tipForce = prestress,
                          double startingVelocity = initialVelocity,
                          double dt = timeStep)
{
    RModel model;
    model.setNNodes(2);
    model.getNode(0).set(0.0,0.0,0.0);
    model.getNode(1).set(cableLength,0.0,0.0);

    RElement cable(R_ELEMENT_TRUSS1);
    cable.setNodeId(0,0);
    cable.setNodeId(1,1);
    model.addElement(cable,true,0);

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
    RLine &line = model.getLine(0);
    line.setCrossArea(cableArea);
    line.setMaterial(material);

    RBoundaryCondition transverse(R_BOUNDARY_CONDITION_DISPLACEMENT);
    setComponent(transverse,R_VARIABLE_DISPLACEMENT_X,false);
    setComponent(transverse,R_VARIABLE_DISPLACEMENT_Y,true,0.0);
    setComponent(transverse,R_VARIABLE_DISPLACEMENT_Z,true,0.0);
    line.addBoundaryCondition(transverse);

    addPoint(model,0);
    addPoint(model,1);
    model.getPoint(0).addBoundaryCondition(displacementCondition({{true,true,true}}));
    model.getPoint(1).addBoundaryCondition(displacementCondition({{false,true,true}}));
    model.getPoint(1).addBoundaryCondition(forceCondition(tipForce));
    setZeroGravity(model);
    setInitialVelocity(model,startingVelocity);

    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    RStructuralCableSetup &cableSetup = model.getProblemSetup().getStructuralCableSetup();
    cableSetup.enabled = true;
    cableSetup.absoluteForceTolerance = 1.0e-12;
    cableSetup.relativeForceTolerance = 1.0e-12;
    RStructuralCableMember member;
    member.elementId = 0;
    member.referenceLength = cableLength;
    member.initialTension = memberPrestress;
    cableSetup.members.push_back(member);
    cableSetup.bindToModel(model);

    RStructuralDynamicsSetup &dynamics =
        model.getProblemSetup().getStructuralDynamicsSetup();
    dynamics.setEnabled(true);
    dynamics.setMethod(RStructuralTimeIntegrationMethod::Newmark);
    dynamics.setNewmarkGamma(0.5);
    dynamics.setNewmarkBeta(0.25);
    dynamics.setTimeStep(dt);
    dynamics.setNSteps(steps);
    dynamics.setOutputFrequency(1);
    return model;
}

bool solve(RModel &model, QString &error, const QString &outputPath = QString())
{
    error.clear();
    try
    {
        RSolver solver(model,outputPath,QString());
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

double nodalValue(const RModel &model, RVariableType type, uint component, uint node)
{
    const uint position = model.findVariable(type);
    if (position == RConstants::eod)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return model.getVariable(position).getValue(component,node);
}

bool closeScaled(double actual, double expected, double scale, double relativeTolerance)
{
    const double bound = relativeTolerance*std::max({scale,std::fabs(actual),std::fabs(expected)});
    return std::isfinite(actual) && std::isfinite(expected) &&
           std::fabs(actual-expected) <= bound;
}

double vibrationEnergy(double q,double v)
{
    const double stiffness = youngModulus*cableArea/cableLength;
    const double endMass = density*cableArea*cableLength/2.0;
    const double tension = prestress+stiffness*q;
    const double strainEnergy = tension*tension/(2.0*stiffness);
    // The constant dead force balances prestress, so subtract its potential.
    return strainEnergy+0.5*endMass*v*v-prestress*q;
}

void compareCheckpoints(const RStructuralDynamicsRestartData &a,
                        const RStructuralDynamicsRestartData &b)
{
    QCOMPARE(a.available,true);
    QCOMPARE(b.available,true);
    QCOMPARE(a.completedSteps,b.completedSteps);
    QCOMPARE(a.meshFingerprint,b.meshFingerprint);
    QCOMPARE(a.dofFingerprint,b.dofFingerprint);
    QCOMPARE(a.systemFingerprint,b.systemFingerprint);
    QCOMPARE(a.q.size(),b.q.size());
    for (std::size_t i=0;i<a.q.size();i++)
    {
        QVERIFY(closeScaled(a.q[i],b.q[i],1.0,1.0e-11));
        QVERIFY(closeScaled(a.v[i],b.v[i],1.0,1.0e-11));
        QVERIFY(closeScaled(a.a[i],b.a[i],1.0,1.0e-11));
        QVERIFY(closeScaled(a.externalForce[i],b.externalForce[i],1.0,1.0e-11));
    }
    QVERIFY(closeScaled(a.time,b.time,timeStep*totalSteps,1.0e-11));
    QVERIFY(closeScaled(a.cumulativeExternalWork,b.cumulativeExternalWork,1.0,1.0e-11));
    QVERIFY(closeScaled(a.initialMechanicalEnergy,b.initialMechanicalEnergy,1.0,1.0e-11));
    QVERIFY(closeScaled(a.cumulativeDampingDissipation,b.cumulativeDampingDissipation,1.0,1.0e-11));
}
}

class TestStructuralCableRestart : public QObject
{
    Q_OBJECT

private slots:
    void persistedSplitRunMatchesFullVibration_data();
    void persistedSplitRunMatchesFullVibration();
    void changedCableSetupRejectsRestart();
    void changedAcceptedStepOrTimeRejectsRestart();
    void rejectedCandidateDoesNotAdvanceAcceptedCableCheckpoint();
    void enabledStructuralControlIsRejectedByCableAdapter();
};

void TestStructuralCableRestart::persistedSplitRunMatchesFullVibration_data()
{
    QTest::addColumn<bool>("binary");
    QTest::newRow("ascii-ram") << false;
    QTest::newRow("binary-rbm") << true;
}

void TestStructuralCableRestart::persistedSplitRunMatchesFullVibration()
{
    QFETCH(bool,binary);
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString extension = binary ? QStringLiteral("rbm") : QStringLiteral("ram");
    const QString fullPath = directory.filePath(QStringLiteral("full.")+extension);
    const QString splitPath = directory.filePath(QStringLiteral("split.")+extension);

    RModel uninterrupted = makeVibratingCable(totalSteps);
    RModel firstSegment = makeVibratingCable(firstSegmentSteps);
    QString error;
    QVERIFY2(solve(uninterrupted,error,fullPath),qPrintable(error));
    QVERIFY2(solve(firstSegment,error,splitPath),qPrintable(error));
    QVERIFY(firstSegment.getProblemSetup().getStructuralDynamicsRestartData().available);
    QCOMPARE(firstSegment.getProblemSetup().getStructuralDynamicsRestartData().completedSteps,
             firstSegmentSteps);

    RModel resumed;
    QVERIFY2(load(resumed,splitPath,error),qPrintable(error));
    QVERIFY(resumed.getProblemSetup().getStructuralCableSetup().enabled);
    QCOMPARE(resumed.getProblemSetup().getStructuralCableSetup().members.size(),std::size_t(1));
    resumed.getProblemSetup().setRestart(true);
    resumed.getProblemSetup().getStructuralDynamicsSetup().setNSteps(
        totalSteps-firstSegmentSteps);
    QVERIFY2(solve(resumed,error,splitPath),qPrintable(error));

    const double stiffness = youngModulus*cableArea/cableLength;
    const double endMass = density*cableArea*cableLength/2.0;
    const double omega = std::sqrt(stiffness/endMass);
    const double discreteAngle = 2.0*std::atan(omega*timeStep/2.0);
    const double expectedEnergy = prestress*prestress/(2.0*stiffness)+
        0.5*endMass*initialVelocity*initialVelocity;

    for (uint step=0;step<totalSteps;step++)
    {
        const uint recordNumber = step+1;
        RModel fullState;
        RModel splitState;
        QVERIFY2(load(fullState,RFileUtils::getFileNameWithTimeStep(fullPath,recordNumber),error),
                 qPrintable(error));
        QVERIFY2(load(splitState,RFileUtils::getFileNameWithTimeStep(splitPath,recordNumber),error),
                 qPrintable(error));
        const auto &fullCheckpoint =
            fullState.getProblemSetup().getStructuralDynamicsRestartData();
        const auto &splitCheckpoint =
            splitState.getProblemSetup().getStructuralDynamicsRestartData();
        compareCheckpoints(fullCheckpoint,splitCheckpoint);
        QCOMPARE(fullCheckpoint.completedSteps,recordNumber);

        // Checkpoint vectors are the exact accepted state (17-digit ASCII
        // persistence); legacy result variables use QTextStream's default
        // precision and are checked separately as display values below.
        const double q = fullCheckpoint.q[3];
        const double v = fullCheckpoint.v[3];
        const double displayedQ = nodalValue(fullState,R_VARIABLE_DISPLACEMENT,0,1);
        const double displayedV = nodalValue(fullState,R_VARIABLE_VELOCITY,0,1);
        const double tension = prestress+stiffness*q;
        const double displayedTension = nodalValue(fullState,R_VARIABLE_FORCE,0,1);
        const double displayTolerance = binary ? 2.0e-7 : 1.0e-5;
        const double time = double(recordNumber)*timeStep;
        const double expectedQ = initialVelocity/omega*std::sin(double(recordNumber)*discreteAngle);
        const double expectedV = initialVelocity*std::cos(double(recordNumber)*discreteAngle);
        const double expectedTension = prestress+stiffness*expectedQ;
        QVERIFY2(closeScaled(q,expectedQ,initialVelocity/omega,2.0e-7),
                 qPrintable(QString("cable displacement differs from the Newmark vibration benchmark at t=%1")
                                .arg(time,0,'g',16)));
        QVERIFY2(closeScaled(v,expectedV,initialVelocity,2.0e-7),
                 qPrintable(QString("cable velocity differs from the Newmark vibration benchmark at t=%1")
                                .arg(time,0,'g',16)));
        QVERIFY(closeScaled(tension,expectedTension,prestress,2.0e-7));
        QVERIFY(closeScaled(vibrationEnergy(q,v),expectedEnergy,expectedEnergy,2.0e-6));
        QVERIFY(closeScaled(displayedQ,q,initialVelocity/omega,displayTolerance));
        QVERIFY(closeScaled(displayedV,v,initialVelocity,displayTolerance));
        QVERIFY(closeScaled(displayedTension,tension,prestress,displayTolerance));
        QVERIFY(closeScaled(nodalValue(splitState,R_VARIABLE_FORCE,0,1),
                            tension,prestress,displayTolerance));
    }

    const auto &fullFinal = uninterrupted.getProblemSetup().getStructuralDynamicsRestartData();
    const auto &splitFinal = resumed.getProblemSetup().getStructuralDynamicsRestartData();
    compareCheckpoints(fullFinal,splitFinal);
    QCOMPARE(splitFinal.completedSteps,totalSteps);
}

void TestStructuralCableRestart::changedCableSetupRejectsRestart()
{
    RModel first = makeVibratingCable(1);
    QString error;
    QVERIFY2(solve(first,error),qPrintable(error));
    QVERIFY(first.getProblemSetup().getStructuralDynamicsRestartData().available);

    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString path = directory.filePath(QStringLiteral("accepted.ram"));
    QCOMPARE(first.write(path,false),path);
    RModel resumed;
    QVERIFY2(load(resumed,path,error),qPrintable(error));
    resumed.getProblemSetup().setRestart(true);
    resumed.getProblemSetup().getStructuralDynamicsSetup().setNSteps(1);
    resumed.getProblemSetup().getStructuralCableSetup().members[0].initialTension += 0.25;

    bool rejected = false;
    try
    {
        RStructuralCableModelAdapter::validate(resumed);
    }
    catch (const RError &failure)
    {
        rejected = failure.getMessage().contains("system fingerprint",Qt::CaseInsensitive);
    }
    QVERIFY(rejected);
}

void TestStructuralCableRestart::changedAcceptedStepOrTimeRejectsRestart()
{
    RModel accepted = makeVibratingCable(1);
    QString error;
    QVERIFY2(solve(accepted,error),qPrintable(error));
    QVERIFY(accepted.getProblemSetup().getStructuralDynamicsRestartData().available);

    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    const QString path = directory.filePath(QStringLiteral("accepted.ram"));
    QCOMPARE(accepted.write(path,false),path);
    RModel resumedByStep;
    QVERIFY2(load(resumedByStep,path,error),qPrintable(error));
    resumedByStep.getProblemSetup().setRestart(true);
    resumedByStep.getProblemSetup().getStructuralDynamicsSetup().setNSteps(1);
    ++resumedByStep.getProblemSetup().getStructuralDynamicsRestartData().completedSteps;
    QVERIFY2(!solve(resumedByStep,error),
             "A checkpoint with a forged accepted-step index must be rejected.");
    QVERIFY2(error.contains("accepted time step",Qt::CaseInsensitive),qPrintable(error));
    QCOMPARE(resumedByStep.getProblemSetup().getStructuralDynamicsRestartData().completedSteps,uint(2));

    RModel resumedByTime;
    QVERIFY2(load(resumedByTime,path,error),qPrintable(error));
    resumedByTime.getProblemSetup().setRestart(true);
    resumedByTime.getProblemSetup().getStructuralDynamicsSetup().setNSteps(1);
    resumedByTime.getProblemSetup().getStructuralDynamicsRestartData().time += timeStep;
    QVERIFY2(!solve(resumedByTime,error),
             "A checkpoint with a forged accepted time must be rejected.");
    QVERIFY2(error.contains("accepted time step",Qt::CaseInsensitive),qPrintable(error));
}

void TestStructuralCableRestart::rejectedCandidateDoesNotAdvanceAcceptedCableCheckpoint()
{
    RModel acceptedModel = makeVibratingCable(1);
    QString error;
    QVERIFY2(solve(acceptedModel,error),qPrintable(error));
    const auto accepted = acceptedModel.getProblemSetup().getStructuralDynamicsRestartData();
    QVERIFY(accepted.available);
    QCOMPARE(accepted.completedSteps,uint(1));

    RStructuralCableDynamicsResult rejectedCandidate;
    rejectedCandidate.displacements = accepted.q;
    rejectedCandidate.velocities = accepted.v;
    rejectedCandidate.accelerations = accepted.a;
    rejectedCandidate.internalForces.assign(accepted.q.size(),0.0);
    rejectedCandidate.elementTensions.assign(1,-1.0);
    bool rejected = false;
    try
    {
        RStructuralCableModelAdapter::publishDynamics(acceptedModel,rejectedCandidate);
    }
    catch (const RError &failure)
    {
        rejected = failure.getMessage().contains("nonnegative",Qt::CaseInsensitive);
    }
    QVERIFY(rejected);

    const auto &afterFailure =
        acceptedModel.getProblemSetup().getStructuralDynamicsRestartData();
    compareCheckpoints(accepted,afterFailure);
}

void TestStructuralCableRestart::enabledStructuralControlIsRejectedByCableAdapter()
{
    RModel model = makeVibratingCable(1);
    model.getProblemSetup().getStructuralControlSetup().enabled = true;
    bool rejected = false;
    try
    {
        RStructuralCableModelAdapter::validate(model);
    }
    catch (const RError &failure)
    {
        rejected = failure.getMessage().contains("Stage-6 structural control",Qt::CaseInsensitive);
    }
    QVERIFY(rejected);
}

QTEST_MAIN(TestStructuralCableRestart)
#include "tst_structural_cable_restart.moc"

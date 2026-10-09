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
#include "rstructural_control_conditioner.h"

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

namespace Runtime = RStructuralControlRuntimeParameters;
const QString baseRevision(64,QLatin1Char('a'));
const QString proposalId = "12345678-1234-4234-8234-123456789abc";

void configureControl(RModel &model)
{
    auto &setup = model.getProblemSetup().getStructuralControlSetup();
    setup.enabled = true;
    setup.controllerId = "runtime-velocity";
    RStructuralActuator actuator;
    actuator.id = "tip-x";
    actuator.target.nodeId = 1;
    actuator.variable = R_VARIABLE_FORCE_X;
    actuator.unit = RVariable::getUnits(R_VARIABLE_FORCE_X);
    actuator.minimum = -5.0e4;
    actuator.maximum = 5.0e4;
    model.getProblemSetup().getStructuralInteractionSetup().getActuators().push_back(actuator);
    RStructuralSensor sensor;
    sensor.id = "tip-velocity";
    sensor.target.nodeId = 1;
    sensor.quantity = RStructuralSensorQuantity::Velocity;
    sensor.direction = RR3Vector(1.0,0.0,0.0);
    sensor.unit = RStructuralInteractionSetup::getSensorUnit(sensor.quantity);
    model.getProblemSetup().getStructuralInteractionSetup().getSensors().push_back(sensor);
}

Runtime::State freshState()
{
    return Runtime::fresh(baseRevision,1.0e4,{{"tip-x",-5.0e4,5.0e4,0.0}});
}

RStructuralControlDecision velocityControl(const RStructuralControlObservation &observation,
                                          bool switching = true)
{
    auto state = observation.priorRuntimeParametersState.isEmpty() ? freshState() :
        Runtime::deserialize(observation.priorRuntimeParametersState);
    RStructuralControlDecision result;
    if (observation.priorRuntimeParametersState.isEmpty())
        result.initialRuntimeParametersState = Runtime::serialize(state);
    auto parameters = state.parameters;
    if (switching && observation.acceptedSeq == 2)
    {
        parameters = Runtime::candidate(state,4.0e4,{{"tip-x",-1200.0,1200.0,5.0e7}});
        result.runtimeParameterCandidate = parameters;
        result.runtimeParameterProposalId = proposalId;
    }
    result.requestedCommands = {{"tip-x",-parameters.velocityGainNsm*observation.sensors[0].value}};
    result.nextControllerState = QByteArray::number(observation.acceptedSeq+1);
    return result;
}

bool run(RSolver &solver,QString &error)
{
    error.clear();
    try { solver.run(); return true; }
    catch (const RError &failure) { error = failure.getMessage(); }
    catch (const std::exception &failure) { error = QString::fromUtf8(failure.what()); }
    return false;
}

struct ScalarReference
{
    // Independently reduced consistent-mass bar: m=rho*A*L/3, k=E*A/L.
    double q = initialDisplacement;
    double v = 0.0;
    double a = -1.0e9*initialDisplacement/(8000.0*0.01*length/3.0);
    double force = 0.0;
    double work = 0.0;
    void step(double nextForce,double dt)
    {
        const double mass = 8000.0*0.01*length/3.0;
        const double stiffness = 2.0e11*0.01/length;
        // The held force changes at the accepted boundary. Re-equilibrate
        // a_n, then hand-evaluate beta=1/4, gamma=1/2 Newmark.
        a = (nextForce-stiffness*q)/mass;
        const double qPredict = q+dt*v+0.25*dt*dt*a;
        const double vPredict = v+0.5*dt*a;
        const double nextA = (nextForce-stiffness*qPredict)/(mass+0.25*dt*dt*stiffness);
        const double nextQ = qPredict+0.25*dt*dt*nextA;
        const double nextV = vPredict+0.5*dt*nextA;
        work += nextForce*(nextQ-q);
        q = nextQ; v = nextV; a = nextA; force = nextForce;
    }
};

bool close(double actual,double expected)
{
    return std::isfinite(actual) && std::isfinite(expected) &&
        std::abs(actual-expected) <= 1.0e-9*std::max({1.0e-12,std::abs(actual),std::abs(expected)})+1.0e-12;
}

double nodal(const RModel &model,RVariableType type)
{
    const uint position = model.findVariable(type);
    return position == RConstants::eod ? std::numeric_limits<double>::quiet_NaN() :
        model.getVariable(position).getValue(0,1);
}
}

class TestStructuralRuntimeParametersFE : public QObject
{
    Q_OBJECT
private slots:
    void initTestCase() { omp_set_num_threads(1); }
    void switchedForcesAndWorkMatchHandNewmark();
    void restartRestoresParametersBeforeObservation();
    void failuresRollback_data();
    void failuresRollback();
    void invalidUpdates_data();
    void invalidUpdates();
    void legacyEmptyEnvelopeRemainsV1();
};

void TestStructuralRuntimeParametersFE::switchedForcesAndWorkMatchHandNewmark()
{
    RModel model = makeTruss();
    configureControl(model);
    model.getProblemSetup().getStructuralDynamicsSetup().setNSteps(7);
    ScalarReference reference;
    const double dt = model.getProblemSetup().getStructuralDynamicsSetup().getTimeStep();
    RSolver solver(model,QString(),QString());
    std::vector<RStructuralControlStepRecord> records;
    std::vector<RStructuralDynamicsRestartData> checkpoints;
    solver.setStructuralControlCallback("runtime-velocity",baseRevision,
        [](const auto &observation){return velocityControl(observation);});
    solver.setStructuralControlStepObserver([&](const auto &record){records.push_back(record);});
    solver.setAcceptedStructuralStepObserver([&](const RModel &accepted)
    { checkpoints.push_back(accepted.getProblemSetup().getStructuralDynamicsRestartData()); });
    QString error;
    QVERIFY2(run(solver,error),qPrintable(error));
    QCOMPARE(records.size(),std::size_t(7));
    QCOMPARE(checkpoints.size(),records.size());
    for (std::size_t i=0;i<records.size();++i)
    {
        const auto &record = records[i];
        const double gain = i<2 ? 1.0e4 : 4.0e4;
        const double request = -gain*reference.v;
        const double magnitude = i<2 ? 5.0e4 : 1200.0;
        const double delta = i<2 ? std::numeric_limits<double>::infinity() : 5.0e7*dt;
        const double force = std::clamp(request,std::max(-magnitude,reference.force-delta),
                                       std::min(magnitude,reference.force+delta));
        QVERIFY(close(record.sensors[0].value,reference.v));
        QVERIFY(close(record.requestedCommands[0].value,request));
        QVERIFY(close(record.conditionedCommands[0].value,request));
        QVERIFY(close(record.appliedChannels[0].command,force));
        reference.step(force,dt);
        const auto &checkpoint = checkpoints[i];
        // Checkpoint vectors retain solver scaling; nodal observations and
        // the independent scalar recurrence use physical SI units.
        const double scale = model.findNodeScale();
        QVERIFY(close(checkpoint.q[0],reference.q*scale));
        QVERIFY(close(checkpoint.v[0],reference.v*scale));
        QVERIFY(close(checkpoint.a[0],reference.a*scale));
        QVERIFY(close(record.appliedChannels[0].cumulativeWork,reference.work));
        QVERIFY(close(checkpoint.cumulativeExternalWork,reference.work*scale*scale));
        QCOMPARE(record.appliedSeq,quint64(i+1));
        const auto runtime = Runtime::deserialize(record.runtimeParametersState);
        QCOMPARE(runtime.parameters.revision,quint64(i<2 ? 0 : 1));
        const auto envelope = RStructuralControlStateEnvelope::deserialize(checkpoint.controlControllerState);
        QCOMPARE(envelope.runtimeParametersState,record.runtimeParametersState);
        if (i>=2)
        {
            QCOMPARE(runtime.history.size(),std::size_t(1));
            QCOMPARE(runtime.history[0].appliedSeq,quint64(3));
            QVERIFY(close(runtime.history[0].acceptedTimeS,3.0*dt));
        }
    }
    QCOMPARE(model.getProblemSetup().getStructuralInteractionSetup().getActuators()[0].minimum,-5.0e4);
    QCOMPARE(model.getProblemSetup().getStructuralInteractionSetup().getActuators()[0].maximum,5.0e4);
}

void TestStructuralRuntimeParametersFE::restartRestoresParametersBeforeObservation()
{
    QTemporaryDir temporary;
    QVERIFY(temporary.isValid());
    const auto path = temporary.filePath("runtime.rbm");
    RModel model = makeTruss();
    configureControl(model);
    model.getProblemSetup().getStructuralDynamicsSetup().setNSteps(3);
    RSolver first(model,path,QString());
    first.setStructuralControlCallback("runtime-velocity",baseRevision,
        [](const auto &observation){return velocityControl(observation);});
    QString error;
    QVERIFY2(run(first,error),qPrintable(error));
    const auto before = model.getProblemSetup().getStructuralDynamicsRestartData();
    const auto expectedEnvelope = RStructuralControlStateEnvelope::deserialize(before.controlControllerState);
    RModel resumed;
    resumed.read(path);
    resumed.getProblemSetup().setRestart(true);
    resumed.getProblemSetup().getStructuralDynamicsSetup().setNSteps(1);
    RSolver second(resumed,QString(),QString());
    std::vector<RStructuralControlObservation> observations;
    second.setStructuralControlCallback("runtime-velocity",baseRevision,[&](const auto &observation)
    { observations.push_back(observation); return velocityControl(observation); });
    QVERIFY2(run(second,error),qPrintable(error));
    QCOMPARE(observations.size(),std::size_t(1));
    QCOMPARE(observations[0].acceptedSeq,quint64(3));
    QCOMPARE(observations[0].priorControllerState,QByteArray("3"));
    QCOMPARE(observations[0].priorRuntimeParametersState,expectedEnvelope.runtimeParametersState);
    const auto after = resumed.getProblemSetup().getStructuralDynamicsRestartData();
    QCOMPARE(RStructuralControlStateEnvelope::deserialize(after.controlControllerState).runtimeParametersState,
             expectedEnvelope.runtimeParametersState);
    ScalarReference reference;
    const double scale = resumed.findNodeScale();
    reference.q = before.q[0]/scale; reference.v = before.v[0]/scale; reference.a = before.a[0]/scale;
    reference.force = before.controlLastCommands[0]; reference.work = before.cumulativeExternalWork/(scale*scale);
    const double dt = resumed.getProblemSetup().getStructuralDynamicsSetup().getTimeStep();
    const double force = std::clamp(-4.0e4*reference.v,
        std::max(-1200.0,reference.force-5.0e7*dt),std::min(1200.0,reference.force+5.0e7*dt));
    reference.step(force,dt);
    QVERIFY(close(after.q[0],reference.q*scale));
    QVERIFY(close(after.v[0],reference.v*scale));
    QVERIFY(close(after.a[0],reference.a*scale));
    QVERIFY(close(after.cumulativeExternalWork,reference.work*scale*scale));
    QVERIFY(close(after.cumulativeControlWork[0],reference.work));
}

void TestStructuralRuntimeParametersFE::failuresRollback_data()
{
    QTest::addColumn<QString>("failureKind");
    QTest::newRow("compute") << QString("compute");
    QTest::newRow("solve") << QString("solve");
    QTest::newRow("late-result-write") << QString("write");
}

void TestStructuralRuntimeParametersFE::failuresRollback()
{
    QFETCH(QString,failureKind);
    QTemporaryDir temporary;
    QVERIFY(temporary.isValid());
    RModel model = makeTruss(); configureControl(model);
    if (failureKind=="solve")
    {
        auto &actuator = model.getProblemSetup().getStructuralInteractionSetup().getActuators()[0];
        actuator.minimum = -std::numeric_limits<double>::max();
        actuator.maximum = std::numeric_limits<double>::max();
    }
    model.getProblemSetup().getStructuralDynamicsSetup().setNSteps(2);
    QString error;
    {
        RSolver first(model,QString(),QString());
        first.setStructuralControlCallback("runtime-velocity",baseRevision,
            [](const auto &observation){return velocityControl(observation);});
        QVERIFY2(run(first,error),qPrintable(error));
    }
    const auto before = model.getProblemSetup().getStructuralDynamicsRestartData();
    const std::array<double,3> kinematics = {nodal(model,R_VARIABLE_DISPLACEMENT),
        nodal(model,R_VARIABLE_VELOCITY),nodal(model,R_VARIABLE_ACCELERATION)};
    model.getProblemSetup().setRestart(true);
    model.getProblemSetup().getStructuralDynamicsSetup().setNSteps(1);
    const auto blockingParent = temporary.filePath("blocked-parent");
    QFile blocker(blockingParent); QVERIFY(blocker.open(QIODevice::WriteOnly));
    QCOMPARE(blocker.write("block"),qint64(5)); blocker.close();
    RSolver retry(model,failureKind=="write" ? blockingParent+"/trial.rbm" : QString(),QString());
    bool failing = true;
    std::vector<RStructuralControlObservation> observations;
    std::vector<RStructuralControlStepRecord> records;
    retry.setStructuralControlCallback("runtime-velocity",baseRevision,[&](const auto &observation)
    {
        observations.push_back(observation);
        auto decision = velocityControl(observation);
        if (failing && failureKind=="compute") throw std::runtime_error("Injected compute failure");
        if (failing && failureKind=="solve")
        {
            const auto state = Runtime::deserialize(observation.priorRuntimeParametersState);
            const double largest = std::numeric_limits<double>::max();
            decision.runtimeParameterCandidate = Runtime::candidate(state,4.0e4,
                {{"tip-x",-largest,largest,0.0}});
            decision.requestedCommands[0].value = 1.0e308;
        }
        return decision;
    });
    retry.setStructuralControlStepObserver([&](const auto &record){records.push_back(record);});
    QVERIFY2(!run(retry,error),qPrintable(error));
    QVERIFY2(!error.isEmpty(),"Failed trial did not report its error");
    qInfo("Runtime-parameter failure fixture %s: %s", qPrintable(failureKind), qPrintable(error));
    if (failureKind=="solve")
        QVERIFY2(error.contains("Dynamics linear solve returned non-finite values"),qPrintable(error));
    QCOMPARE(records.size(),std::size_t(0));
    const auto after = model.getProblemSetup().getStructuralDynamicsRestartData();
    QCOMPARE(after.completedSteps,before.completedSteps);
    QCOMPARE(after.time,before.time);
    QCOMPARE(after.controlControllerState,before.controlControllerState);
    QCOMPARE(after.controlLastCommands,before.controlLastCommands);
    QCOMPARE(after.controlPowers,before.controlPowers);
    QCOMPARE(after.cumulativeControlWork,before.cumulativeControlWork);
    QCOMPARE(after.cumulativeExternalWork,before.cumulativeExternalWork);
    QCOMPARE(after.q,before.q); QCOMPARE(after.v,before.v); QCOMPARE(after.a,before.a);
    QCOMPARE(nodal(model,R_VARIABLE_DISPLACEMENT),kinematics[0]);
    QCOMPARE(nodal(model,R_VARIABLE_VELOCITY),kinematics[1]);
    QCOMPARE(nodal(model,R_VARIABLE_ACCELERATION),kinematics[2]);
    if (failureKind=="write")
    {
        QVERIFY(blocker.rename(temporary.filePath("old-blocker")));
        QVERIFY(QDir().mkpath(blockingParent));
    }
    failing = false;
    QVERIFY2(run(retry,error),qPrintable(error));
    QCOMPARE(records.size(),std::size_t(1));
    QCOMPARE(observations.size(),std::size_t(2));
    QCOMPARE(observations[1].priorRuntimeParametersState,observations[0].priorRuntimeParametersState);
    QCOMPARE(observations[1].priorControllerState,observations[0].priorControllerState);
    QCOMPARE(observations[1].acceptedSeq,observations[0].acceptedSeq);
    const auto accepted = Runtime::deserialize(records[0].runtimeParametersState);
    QCOMPARE(accepted.parameters.revision,quint64(1));
    QCOMPARE(accepted.history[0].appliedSeq,quint64(3));
}

void TestStructuralRuntimeParametersFE::invalidUpdates_data()
{
    QTest::addColumn<QString>("invalidKind");
    for (const auto *kind : {"wrong-base","wrong-id","wrong-revision","noncanonical-bootstrap",
                            "bootstrap-history","zero-model-rate","other-method","amplitude","empty-intersection",
                            "time-mutation","catalog-mutation","result-mutation","alpha-to-newmark-mutation"})
        QTest::newRow(kind) << QString(kind);
}

void TestStructuralRuntimeParametersFE::invalidUpdates()
{
    QFETCH(QString,invalidKind);
    RModel model = makeTruss(); configureControl(model);
    if (invalidKind=="zero-model-rate")
        model.getProblemSetup().getStructuralInteractionSetup().getActuators()[0].rateLimit = 0.0;
    if (invalidKind=="other-method" || invalidKind=="alpha-to-newmark-mutation")
        model.getProblemSetup().getStructuralDynamicsSetup().setMethod(RStructuralTimeIntegrationMethod::GeneralizedAlpha);
    RSolver solver(model,QString(),QString());
    bool notified = false;
    solver.setStructuralControlStepObserver([&](const auto &){notified=true;});
    solver.setStructuralControlCallback("runtime-velocity",baseRevision,[&](const auto &observation)
    {
        auto result = velocityControl(observation);
        auto state = freshState();
        if (invalidKind=="wrong-base")
            result.initialRuntimeParametersState = Runtime::serialize(
                Runtime::fresh(QString(64,'b'),1.0e4,state.parameters.actuatorLimits));
        if (invalidKind=="wrong-id")
            result.initialRuntimeParametersState = Runtime::serialize(
                Runtime::fresh(baseRevision,1.0e4,{{"other-id",-100.0,100.0,0.0}}));
        if (invalidKind=="amplitude")
            result.initialRuntimeParametersState = Runtime::serialize(
                Runtime::fresh(baseRevision,1.0e4,{{"tip-x",-1.0e5,1.0e5,0.0}}));
        if (invalidKind=="empty-intersection")
            result.initialRuntimeParametersState = Runtime::serialize(
                Runtime::fresh(baseRevision,1.0e4,{{"tip-x",100.0,200.0,1.0}}));
        if (invalidKind=="noncanonical-bootstrap") result.initialRuntimeParametersState += ' ';
        if (invalidKind=="bootstrap-history")
            result.initialRuntimeParametersState = Runtime::serialize(Runtime::acceptedChange(
                state,proposalId,Runtime::candidate(state,2.0e4,state.parameters.actuatorLimits),1,1.0));
        if (invalidKind=="wrong-revision")
        {
            const auto future = Runtime::acceptedChange(state,proposalId,
                Runtime::candidate(state,2.0e4,state.parameters.actuatorLimits),1,1.0);
            result.runtimeParameterCandidate = Runtime::candidate(future,3.0e4,state.parameters.actuatorLimits);
            result.runtimeParameterProposalId = proposalId;
        }
        if (invalidKind=="time-mutation") model.getTimeSolver().setComputedTime(999.0);
        if (invalidKind=="alpha-to-newmark-mutation")
            model.getProblemSetup().getStructuralDynamicsSetup().setMethod(RStructuralTimeIntegrationMethod::Newmark);
        if (invalidKind=="catalog-mutation")
            model.getProblemSetup().getStructuralInteractionSetup().getActuators()[0].maximum = 1.0e6;
        if (invalidKind=="result-mutation")
        {
            RVariable pressure(R_VARIABLE_PRESSURE);
            pressure.setApplyType(R_VARIABLE_APPLY_NODE);
            pressure.resize(1,model.getNNodes());
            pressure.setValue(0,1,321.0);
            model.addVariable(pressure);
        }
        return result;
    });
    QString error;
    QVERIFY2(!run(solver,error),qPrintable(error));
    QVERIFY(!notified);
    QVERIFY(!model.getProblemSetup().getStructuralDynamicsRestartData().available);
}

void TestStructuralRuntimeParametersFE::legacyEmptyEnvelopeRemainsV1()
{
    RModel model = makeTruss(); configureControl(model);
    RSolver solver(model,QString(),QString());
    solver.setStructuralControlCallback("runtime-velocity","legacy-v1",
        [](const auto &){return RStructuralControlDecision{{{"tip-x",0.0}},QByteArray("legacy")};});
    QString error;
    QVERIFY2(run(solver,error),qPrintable(error));
    const auto &checkpoint = model.getProblemSetup().getStructuralDynamicsRestartData();
    const auto envelope = RStructuralControlStateEnvelope::deserialize(checkpoint.controlControllerState);
    QVERIFY(envelope.runtimeParametersState.isEmpty());
    QCOMPARE(checkpoint.controlControllerState,
             (RStructuralControlStateEnvelope{envelope.controllerState,envelope.conditionerState}).serialize());
    QCOMPARE(checkpoint.controlControllerState.mid(4,4),QByteArray::fromHex("01000000"));
}

QTEST_APPLESS_MAIN(TestStructuralRuntimeParametersFE)
#include "tst_structural_runtime_parameters_fe.moc"

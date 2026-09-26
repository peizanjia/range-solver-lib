#include <QtTest>

#include <algorithm>
#include <cmath>
#include <exception>
#include <limits>
#include <stdexcept>
#include <vector>

#include "rml_model.h"
#include "rbl_logger.h"
#include "rsolver.h"
#include "rstructural_control_conditioner.h"

namespace
{
constexpr double youngModulus = 2.0e11;
constexpr double massDensity = 8000.0;
constexpr double length = 2.0;
constexpr double crossArea = 0.01;
constexpr double initialDisplacement = 1.0e-4;

// The fixed-free, one-element truss has K=EA/L and a free-node consistent
// mass M=rho*AL/3. This gives an independent analytic check on the FE result.
double singleElementFrequency()
{
    return std::sqrt(3.0*youngModulus/(massDensity*length*length));
}

RBoundaryCondition fixedDisplacement()
{
    RBoundaryCondition bc(R_BOUNDARY_CONDITION_DISPLACEMENT);
    bc.getComponent(bc.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).add(0.0,0.0);
    bc.getComponent(bc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).add(0.0,0.0);
    bc.getComponent(bc.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).add(0.0,0.0);
    return bc;
}

RModel makeTruss(uint nElements, bool tieLastTwo = false)
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
    RMaterialProperty density(RMaterialProperty::Density);
    density.add(293.15,massDensity);
    material.add(density);
    RMaterialProperty expansion(RMaterialProperty::ThermalExpansionCoefficient);
    expansion.add(293.15,0.0);
    material.add(expansion);
    RLine &line = model.getLine(0);
    line.setCrossArea(crossArea);
    line.setMaterial(material);

    // A truss has only axial stiffness. Constrain the unused transverse DOFs
    // so that the full structural equation system remains nonsingular.
    RBoundaryCondition transverse(R_BOUNDARY_CONDITION_DISPLACEMENT);
    transverse.getComponent(transverse.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).setEnabled(false);
    transverse.getComponent(transverse.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).add(0.0,0.0);
    transverse.getComponent(transverse.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).add(0.0,0.0);
    line.addBoundaryCondition(transverse);

    RElement fixedPoint(R_ELEMENT_POINT);
    fixedPoint.setNodeId(0,0);
    model.addElement(fixedPoint,true,0);
    model.getPoint(0).addBoundaryCondition(fixedDisplacement());

    // Point entities are needed for node-specific initial conditions. The
    // tied pair starts on the constraint manifold, q_last=q_previous.
    const uint firstInitialized = tieLastTwo ? nElements-1 : nElements;
    for (uint node=firstInitialized;node<=nElements;node++)
    {
        RElement point(R_ELEMENT_POINT);
        point.setNodeId(0,node);
        model.addElement(point,true,1+node-firstInitialized);
        RInitialCondition initial(R_INITIAL_CONDITION_DISPLACEMENT);
        initial.getComponent(initial.findComponentPosition(R_VARIABLE_DISPLACEMENT_X))
               .add(0.0,initialDisplacement);
        model.getPoint(1+node-firstInitialized).addInitialCondition(initial);
    }

    if (tieLastTwo)
    {
        RLinearConstraint constraint;
        constraint.id = "dynamic-last-segment-tie";
        constraint.terms = {{{nElements,RDofType::UX},1.0},
                            {{nElements-1,RDofType::UX},-1.0}};
        model.getProblemSetup().getStructuralInteractionSetup()
             .getLinearConstraints().push_back(constraint);
    }

    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setOutputFrequency(0);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setSolverCvgValue(1.0e-13);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setNInnerIterations(100);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setNOuterIterations(1000);
    RStructuralDynamicsSetup &dynamics = model.getProblemSetup().getStructuralDynamicsSetup();
    dynamics.setEnabled(true);
    dynamics.setMethod(RStructuralTimeIntegrationMethod::Newmark);
    dynamics.setTimeStep(0.025/singleElementFrequency());
    dynamics.setNSteps(40);
    dynamics.setOutputFrequency(40);
    return model;
}

bool runSolver(RModel &model, QString &error,
               std::vector<std::vector<RStructuralSupportSample>> *supportHistory = nullptr,
               std::vector<RStructuralSupportSample> *lastSupport = nullptr)
{
    try
    {
        RSolver solver(model,QString(),QString());
        if (supportHistory || lastSupport)
        {
            solver.setStructuralSupportResponseObserver(
                [supportHistory](const std::vector<RStructuralSupportSample> &samples)
                {
                    if (supportHistory) supportHistory->push_back(samples);
                });
        }
        solver.run();
        if (lastSupport) *lastSupport = solver.getStructuralSupportSamples();
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

double nodalValue(const RModel &model, RVariableType type, uint node)
{
    const uint index = model.findVariable(type);
    if (index == RConstants::eod)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return model.getVariable(index).getValue(0,node);
}

double dampedDisplacement(double omega, double dampingRatio, double time)
{
    const double dampedOmega = omega*std::sqrt(1.0-dampingRatio*dampingRatio);
    return initialDisplacement*std::exp(-dampingRatio*omega*time)
         *(std::cos(dampedOmega*time)
           +dampingRatio*omega/dampedOmega*std::sin(dampedOmega*time));
}

void configureControl(RModel &model)
{
    auto &setup = model.getProblemSetup().getStructuralControlSetup();
    setup.enabled = true;
    setup.controllerId = "constant-force";
    setup.deterministicSeed = 17;
    RStructuralActuator actuator;
    actuator.id = "tip-x";
    actuator.target.nodeId = model.getNNodes()-1;
    actuator.variable = R_VARIABLE_FORCE_X;
    actuator.unit = RVariable::getUnits(R_VARIABLE_FORCE_X);
    actuator.minimum = -5.0e4;
    actuator.maximum = 5.0e4;
    model.getProblemSetup().getStructuralInteractionSetup().getActuators().push_back(actuator);
    RStructuralSensor sensor;
    sensor.id = "tip-displacement";
    sensor.target.nodeId = model.getNNodes()-1;
    sensor.quantity = RStructuralSensorQuantity::Displacement;
    sensor.direction = RR3Vector(1.0,0.0,0.0);
    sensor.unit = RStructuralInteractionSetup::getSensorUnit(sensor.quantity);
    model.getProblemSetup().getStructuralInteractionSetup().getSensors().push_back(sensor);
}

RStructuralControlDecision constantControl(const RStructuralControlObservation &observation)
{
    RStructuralControlDecision decision;
    decision.requestedCommands = {{"tip-x",1.0e5}};
    decision.nextControllerState = QByteArray::number(observation.priorControllerState.toUInt()+1);
    return decision;
}

const RStructuralSupportSample *findSupport(
    const std::vector<RStructuralSupportSample> &samples, uint node, RDofType dof)
{
    for (const RStructuralSupportSample &sample : samples)
    {
        if (sample.target.nodeId == node && sample.target.dof == dof) return &sample;
    }
    return nullptr;
}
}

class TestStructuralDynamicsFE : public QObject
{
    Q_OBJECT

private slots:
    void freeVibrationMatchesAnalyticAfterManySteps();
    void controlClampsAndMatchesAnalyticHeldForce();
    void controlRestartMatchesAndObserverFailureIsReported();
    void zeroControlMatchesPassiveAcrossMethods();
    void velocityFeedbackReducesMechanicalEnergy();
    void controlRateLimitAndRestartFailClosed();
    void controlCallbackCannotMutateCatalog();
    void controlCallbackCannotMutateClockOrConditioning();
    void controlCallbackCannotMutateUnrelatedResults();
    void conditionedDelayHoldNoiseRestartsAfterRejectedTrial();
    void rawLegacyControlStateFailsClosed();
    void twoChannelMappingIsLinearAndRestartsDeterministically();
    void rayleighDampingMatchesAnalyticAfterManySteps();
    void modalDampingMatchesAnalyticAfterManySteps();
    void twoModeDampingMatchesAnalyticProjection();
    void modalDampingWorksWithHomogeneousMpc();
    void unsupportedModalDampingCombinationsFailBeforeFirstStep();
    void stationarySupportReactionIncludesInertia();
    void movingSupportReactionPowerAndAppliedSupportForce();
    void supportGetterRequiresObserver();
    void supportObserverExceptionsDoNotRevokeAcceptedStep();
    void supportObserverRejectsModalAndMpc();
    void movingSupportMatchesConsistentMassAnalytic();
    void movingSupportRayleighMassDampingMatchesAnalytic();
    void unsupportedMovingSupportInputsFailBeforeFirstStep();
    void tiedNodesStayTiedDuringDynamics();
    void fixedThreadCountsProduceSameState();
    void restartMatchesUninterruptedRun();
    void restartRejectsChangedMeshAndTimeStep();
    void restartRejectsChangedOperatorAndMpc();
};

void TestStructuralDynamicsFE::freeVibrationMatchesAnalyticAfterManySteps()
{
    RModel model = makeTruss(1);
    QString error;
    QVERIFY2(runSolver(model,error),qPrintable(error));
    const double omega = singleElementFrequency();
    const double time = 40.0*model.getProblemSetup().getStructuralDynamicsSetup().getTimeStep();
    const double expectedQ = initialDisplacement*std::cos(omega*time);
    const double expectedV = -initialDisplacement*omega*std::sin(omega*time);
    const double q = nodalValue(model,R_VARIABLE_DISPLACEMENT,1);
    const double v = nodalValue(model,R_VARIABLE_VELOCITY,1);
    const double a = nodalValue(model,R_VARIABLE_ACCELERATION,1);
    QVERIFY2(std::isfinite(q) && std::isfinite(v) && std::isfinite(a),
             "The FE result must contain displacement, velocity and acceleration.");
    QVERIFY(std::abs(q-expectedQ) < 2.0e-4*initialDisplacement);
    QVERIFY(std::abs(v-expectedV) < 2.0e-4*initialDisplacement*omega);
    QVERIFY(std::abs(a+omega*omega*q) < 1.0e-7*omega*omega*initialDisplacement);
}

void TestStructuralDynamicsFE::controlClampsAndMatchesAnalyticHeldForce()
{
    RModel model = makeTruss(1);
    configureControl(model);
    std::vector<RStructuralControlStepRecord> records;
    RSolver solver(model,QString(),QString());
    solver.setStructuralControlCallback("constant-force","constant-force-test-v1",constantControl);
    solver.setStructuralControlStepObserver(
        [&records](const RStructuralControlStepRecord &record) { records.push_back(record); });
    try { solver.run(); }
    catch (const RError &error) { RLogger::warning("control FE test error: %s\n",qPrintable(error.getMessage())); QFAIL(qPrintable(error.getMessage())); }
    catch (const std::exception &error) { RLogger::warning("control FE test error: %s\n",error.what()); QFAIL(error.what()); }
    QCOMPARE(records.size(),size_t(40));
    QCOMPARE(records.front().requestedCommands[0].value,1.0e5);
    QCOMPARE(records.front().appliedChannels[0].command,5.0e4);
    QCOMPARE(records.front().sensors[0].sensorId,QString("tip-displacement"));
    QVERIFY(std::abs(records.front().sensors[0].value-initialDisplacement) < 1.0e-12);
    const double omega = singleElementFrequency();
    const double time = 40.0*model.getProblemSetup().getStructuralDynamicsSetup().getTimeStep();
    const double equilibrium = 5.0e4/(youngModulus*crossArea/length);
    const double expected = equilibrium+(initialDisplacement-equilibrium)*std::cos(omega*time);
    const double actual = nodalValue(model,R_VARIABLE_DISPLACEMENT,1);
    QVERIFY(std::abs(actual-expected) < 2.0e-4*initialDisplacement);
    const auto &checkpoint = model.getProblemSetup().getStructuralDynamicsRestartData();
    QVERIFY(checkpoint.controlAvailable);
    QCOMPARE(checkpoint.controlLastCommands[0],5.0e4);
    QCOMPARE(RStructuralControlStateEnvelope::deserialize(
                 checkpoint.controlControllerState).controllerState,QByteArray("40"));
    const double expectedWork = 5.0e4*(actual-initialDisplacement);
    QVERIFY(std::abs(checkpoint.cumulativeControlWork[0]-expectedWork) < 1.0e-6);
    QCOMPARE(solver.getStructuralControlSamples()[0].cumulativeWork,
             checkpoint.cumulativeControlWork[0]);
}

void TestStructuralDynamicsFE::controlRestartMatchesAndObserverFailureIsReported()
{
    RModel uninterrupted = makeTruss(1);
    configureControl(uninterrupted);
    RSolver reference(uninterrupted,QString(),QString());
    reference.setStructuralControlCallback("constant-force","constant-force-test-v1",constantControl);
    reference.run();

    RModel resumed = makeTruss(1);
    configureControl(resumed);
    resumed.getProblemSetup().getStructuralDynamicsSetup().setNSteps(17);
    {
        RSolver first(resumed,QString(),QString());
        first.setStructuralControlCallback("constant-force","constant-force-test-v1",constantControl);
        first.run();
    }
    QCOMPARE(RStructuralControlStateEnvelope::deserialize(resumed.getProblemSetup()
                 .getStructuralDynamicsRestartData().controlControllerState).controllerState,
             QByteArray("17"));
    resumed.getProblemSetup().setRestart(true);
    resumed.getProblemSetup().getStructuralDynamicsSetup().setNSteps(23);
    {
        RSolver second(resumed,QString(),QString());
        second.setStructuralControlCallback("constant-force","constant-force-test-v1",constantControl);
        second.run();
    }
    QVERIFY(std::abs(nodalValue(resumed,R_VARIABLE_DISPLACEMENT,1)-
                     nodalValue(uninterrupted,R_VARIABLE_DISPLACEMENT,1)) < 1.0e-12);
    QCOMPARE(RStructuralControlStateEnvelope::deserialize(resumed.getProblemSetup()
                 .getStructuralDynamicsRestartData().controlControllerState).controllerState,
             QByteArray("40"));

    RModel audited = makeTruss(1);
    configureControl(audited);
    RSolver auditSolver(audited,QString(),QString());
    auditSolver.setStructuralControlCallback("constant-force","constant-force-test-v1",constantControl);
    auditSolver.setStructuralControlStepObserver([](const RStructuralControlStepRecord &)
        { throw std::runtime_error("audit unavailable"); });
    bool failed = false;
    try { auditSolver.run(); }
    catch (const RError &error)
    {
        failed = true;
        QVERIFY(error.getMessage().contains("observer",Qt::CaseInsensitive));
    }
    QVERIFY(failed);
    QCOMPARE(audited.getProblemSetup().getStructuralDynamicsRestartData().completedSteps,uint(1));
}

void TestStructuralDynamicsFE::zeroControlMatchesPassiveAcrossMethods()
{
    for (RStructuralTimeIntegrationMethod method :
         {RStructuralTimeIntegrationMethod::Newmark,
          RStructuralTimeIntegrationMethod::GeneralizedAlpha,
          RStructuralTimeIntegrationMethod::CentralDifference})
    {
        RModel passive = makeTruss(1);
        RModel controlled = makeTruss(1);
        passive.getProblemSetup().getStructuralDynamicsSetup().setMethod(method);
        controlled.getProblemSetup().getStructuralDynamicsSetup().setMethod(method);
        if (method == RStructuralTimeIntegrationMethod::GeneralizedAlpha)
        {
            passive.getProblemSetup().getStructuralDynamicsSetup().setRhoInfinity(0.6);
            controlled.getProblemSetup().getStructuralDynamicsSetup().setRhoInfinity(0.6);
        }
        configureControl(controlled);
        QString error;
        QVERIFY2(runSolver(passive,error),qPrintable(error));
        RSolver solver(controlled,QString(),QString());
        solver.setStructuralControlCallback("constant-force","zero-command-test-v1",
            [](const RStructuralControlObservation &observation)
            {
                RStructuralControlDecision decision;
                decision.requestedCommands = {{"tip-x",0.0}};
                decision.nextControllerState = observation.priorControllerState;
                return decision;
            });
        try { solver.run(); }
        catch (const RError &failure) { QFAIL(qPrintable(failure.getMessage())); }
        const auto &passiveState = passive.getProblemSetup().getStructuralDynamicsRestartData();
        const auto &controlledState = controlled.getProblemSetup().getStructuralDynamicsRestartData();
        QVERIFY(passiveState.available && controlledState.available);
        QCOMPARE(passiveState.q.size(),controlledState.q.size());
        for (size_t i=0;i<passiveState.q.size();i++)
        {
            QVERIFY(std::abs(passiveState.q[i]-controlledState.q[i]) < 1.0e-12);
            QVERIFY(std::abs(passiveState.v[i]-controlledState.v[i]) < 1.0e-12);
            QVERIFY(std::abs(passiveState.a[i]-controlledState.a[i]) < 1.0e-12);
        }
        QCOMPARE(controlledState.controlLastCommands[0],0.0);
        QCOMPARE(controlledState.cumulativeControlWork[0],0.0);
    }
}

void TestStructuralDynamicsFE::velocityFeedbackReducesMechanicalEnergy()
{
    RModel passive = makeTruss(1);
    RModel controlled = makeTruss(1);
    configureControl(controlled);
    auto &catalog = controlled.getProblemSetup().getStructuralInteractionSetup();
    catalog.getActuators()[0].minimum = -1.0e5;
    catalog.getActuators()[0].maximum = 1.0e5;
    RStructuralSensor velocity;
    velocity.id = "tip-velocity";
    velocity.target.nodeId = 1;
    velocity.quantity = RStructuralSensorQuantity::Velocity;
    velocity.direction = RR3Vector(1.0,0.0,0.0);
    velocity.unit = RStructuralInteractionSetup::getSensorUnit(velocity.quantity);
    catalog.getSensors().push_back(velocity);
    QString error;
    QVERIFY2(runSolver(passive,error),qPrintable(error));
    const double modalMass = massDensity*crossArea*length/3.0;
    const double stiffness = youngModulus*crossArea/length;
    const double damper = 2.0*0.3*modalMass*singleElementFrequency();
    RSolver solver(controlled,QString(),QString());
    solver.setStructuralControlCallback("constant-force","velocity-feedback-test-v1",
        [damper](const RStructuralControlObservation &observation)
        {
            RStructuralControlDecision decision;
            double speed = 0.0;
            for (const auto &sensor : observation.sensors)
                if (sensor.sensorId == "tip-velocity") speed = sensor.value;
            decision.requestedCommands = {{"tip-x",-damper*speed}};
            decision.nextControllerState = observation.priorControllerState;
            return decision;
        });
    try { solver.run(); }
    catch (const RError &failure) { QFAIL(qPrintable(failure.getMessage())); }
    auto energy = [modalMass,stiffness](const RModel &model)
    {
        const double q = nodalValue(model,R_VARIABLE_DISPLACEMENT,1);
        const double v = nodalValue(model,R_VARIABLE_VELOCITY,1);
        return 0.5*modalMass*v*v+0.5*stiffness*q*q;
    };
    QVERIFY(energy(controlled) < 0.8*energy(passive));
    QVERIFY(controlled.getProblemSetup().getStructuralDynamicsRestartData()
                .cumulativeControlWork[0] < 0.0);
}

void TestStructuralDynamicsFE::controlRateLimitAndRestartFailClosed()
{
    RModel model = makeTruss(1);
    configureControl(model);
    auto &dynamics = model.getProblemSetup().getStructuralDynamicsSetup();
    dynamics.setNSteps(4);
    model.getProblemSetup().getStructuralInteractionSetup().getActuators()[0].rateLimit =
        1.0e4/dynamics.getTimeStep();
    std::vector<double> applied;
    RSolver solver(model,QString(),QString());
    solver.setStructuralControlCallback("constant-force","rate-limit-test-v1",constantControl);
    solver.setStructuralControlStepObserver([&applied](const RStructuralControlStepRecord &record)
        { applied.push_back(record.appliedChannels[0].command); });
    try { solver.run(); }
    catch (const RError &failure) { QFAIL(qPrintable(failure.getMessage())); }
    QCOMPARE(applied.size(),size_t(4));
    for (size_t i=0;i<applied.size();i++)
        QVERIFY(std::abs(applied[i]-1.0e4*double(i+1)) < 1.0e-7);
    const auto saved = model.getProblemSetup().getStructuralDynamicsRestartData();
    model.getProblemSetup().setRestart(true);
    dynamics.setNSteps(1);
    {
        RSolver changed(model,QString(),QString());
        changed.setStructuralControlCallback("constant-force","changed-provider-revision",constantControl);
        bool rejected = false;
        try { changed.run(); }
        catch (const RError &failure)
        {
            rejected = failure.getMessage().contains("control",Qt::CaseInsensitive);
        }
        QVERIFY(rejected);
    }
    QCOMPARE(model.getProblemSetup().getStructuralDynamicsRestartData().completedSteps,
             saved.completedSteps);
    model.getProblemSetup().getStructuralControlSetup().enabled = false;
    model.getProblemSetup().getStructuralControlSetup().controllerId.clear();
    model.getProblemSetup().getStructuralControlSetup().deterministicSeed = 0;
    {
        RSolver disabled(model,QString(),QString());
        bool rejected = false;
        try { disabled.run(); }
        catch (const RError &failure)
        {
            rejected = failure.getMessage().contains("control",Qt::CaseInsensitive);
        }
        QVERIFY(rejected);
    }
}

void TestStructuralDynamicsFE::controlCallbackCannotMutateCatalog()
{
    RModel model = makeTruss(1);
    configureControl(model);
    const double originalLimit = model.getProblemSetup().getStructuralInteractionSetup()
        .getActuators()[0].maximum;
    RSolver solver(model,QString(),QString());
    solver.setStructuralControlCallback("constant-force","mutating-provider-test-v1",
        [&model](const RStructuralControlObservation &)
        {
            model.getProblemSetup().getStructuralInteractionSetup()
                .getActuators()[0].maximum = 1.0e8;
            return RStructuralControlDecision{{{"tip-x",1.0e5}},QByteArray()};
        });
    bool rejected = false;
    try { solver.run(); }
    catch (const RError &failure)
    {
        rejected = failure.getMessage().contains("catalog",Qt::CaseInsensitive);
    }
    QVERIFY(rejected);
    QCOMPARE(model.getProblemSetup().getStructuralInteractionSetup().getActuators()[0].maximum,
             originalLimit);
    QVERIFY(!model.getProblemSetup().getStructuralDynamicsRestartData().available);
}

void TestStructuralDynamicsFE::controlCallbackCannotMutateClockOrConditioning()
{
    for (const bool changeClock : {true,false})
    {
        RModel model = makeTruss(1);
        configureControl(model);
        auto &conditioning = model.getProblemSetup().getStructuralControlSetup()
            .channelConditioning;
        conditioning.push_back({QStringLiteral("tip-x"),0,1,0.0});
        RSolver solver(model,QString(),QString());
        bool callbackCalled = false;
        solver.setStructuralControlCallback(
            "constant-force","mutating-time-or-setup-v1",
            [&model,&callbackCalled,changeClock](const RStructuralControlObservation &)
            {
                callbackCalled = true;
                if (changeClock)
                    model.getTimeSolver().setCurrentTimeStep(1);
                else
                    ++model.getProblemSetup().getStructuralControlSetup()
                        .channelConditioning.front().delaySteps;
                return RStructuralControlDecision{{{"tip-x",1.0e4}},QByteArray()};
            });
        QString error;
        try { solver.run(); }
        catch (const RError &failure) { error = failure.getMessage(); }
        QVERIFY(callbackCalled);
        QVERIFY2(error.contains(QStringLiteral("frozen time solver or control setup"),
                                Qt::CaseInsensitive),qPrintable(error));
        QVERIFY(!model.getProblemSetup().getStructuralDynamicsRestartData().available);
        QCOMPARE(model.getProblemSetup().getStructuralControlSetup()
                     .channelConditioning.front().delaySteps,0u);
    }
}

void TestStructuralDynamicsFE::controlCallbackCannotMutateUnrelatedResults()
{
    RModel model = makeTruss(1);
    configureControl(model);
    RSolver solver(model,QString(),QString());
    bool first = true;
    std::vector<QByteArray> priorStates;
    solver.setStructuralControlCallback("constant-force","mutating-results-v1",
        [&model,&first,&priorStates](const RStructuralControlObservation &observation)
        {
            priorStates.push_back(observation.priorControllerState);
            if (first)
            {
                first = false;
                RVariable pressure(R_VARIABLE_PRESSURE);
                pressure.setApplyType(R_VARIABLE_APPLY_NODE);
                pressure.resize(1,model.getNNodes());
                pressure.setValue(0,1,321.0);
                model.addVariable(pressure);
            }
            return RStructuralControlDecision{{{"tip-x",1.0e4}},QByteArrayLiteral("accepted-guard")};
        });
    QString error;
    try { solver.run(); }
    catch (const RError &failure) { error = failure.getMessage(); }
    QVERIFY2(error.contains("frozen result fields",Qt::CaseInsensitive),qPrintable(error));
    QVERIFY(!model.getProblemSetup().getStructuralDynamicsRestartData().available);
    QCOMPARE(model.findVariable(R_VARIABLE_PRESSURE),RConstants::eod);
    try { solver.run(); }
    catch (const RError &failure) { QFAIL(qPrintable(failure.getMessage())); }
    QVERIFY(priorStates.size()>=2);
    QCOMPARE(priorStates[0],priorStates[1]);
    QVERIFY(model.getProblemSetup().getStructuralDynamicsRestartData().controlAvailable);
    QCOMPARE(model.findVariable(R_VARIABLE_PRESSURE),RConstants::eod);
}

void TestStructuralDynamicsFE::conditionedDelayHoldNoiseRestartsAfterRejectedTrial()
{
    auto configured = []()
    {
        RModel model = makeTruss(1);
        configureControl(model);
        model.getProblemSetup().getStructuralDynamicsSetup().setNSteps(6);
        model.getProblemSetup().getStructuralControlSetup().channelConditioning =
            {{"tip-x",2,2,500.0}};
        return model;
    };
    auto provider = [](const RStructuralControlObservation &observation)
    {
        RStructuralControlDecision decision;
        decision.requestedCommands = {{"tip-x",1.0e4}};
        decision.nextControllerState = QByteArray::number(observation.priorControllerState.toUInt()+1);
        return decision;
    };
    RModel baseline = configured();
    std::vector<RStructuralControlStepRecord> expected;
    {
        RSolver solver(baseline,QString(),QString());
        solver.setStructuralControlCallback("constant-force","conditioned-test-v1",provider);
        solver.setStructuralControlStepObserver([&expected](const RStructuralControlStepRecord &record)
            { expected.push_back(record); });
        try { solver.run(); }
        catch (const RError &failure) { QFAIL(qPrintable(failure.getMessage())); }
    }
    QCOMPARE(expected.size(),size_t(6));
    QCOMPARE(expected[0].conditionedCommands[0].value,0.0);
    QCOMPARE(expected[1].conditionedCommands[0].value,0.0);
    QVERIFY(expected[2].conditionedCommands[0].value != 0.0);
    QCOMPARE(expected[2].conditionedCommands[0].value,
             expected[3].conditionedCommands[0].value);

    RModel resumed = configured();
    std::vector<RStructuralControlStepRecord> actual;
    {
        RSolver first(resumed,QString(),QString());
        first.setStructuralControlCallback("constant-force","conditioned-test-v1",
            [provider](const RStructuralControlObservation &observation)
            {
                if (observation.priorControllerState == "2")
                    throw std::runtime_error("injected controller failure");
                return provider(observation);
            });
        first.setStructuralControlStepObserver([&actual](const RStructuralControlStepRecord &record)
            { actual.push_back(record); });
        bool failed = false;
        try { first.run(); }
        catch (const RError &failure)
        {
            failed = failure.getMessage().contains("injected controller failure");
        }
        QVERIFY(failed);
    }
    const auto &interrupted = resumed.getProblemSetup().getStructuralDynamicsRestartData();
    QCOMPARE(interrupted.completedSteps,uint(2));
    QCOMPARE(RStructuralControlStateEnvelope::deserialize(
                 interrupted.controlControllerState).controllerState,QByteArray("2"));
    resumed.getProblemSetup().setRestart(true);
    resumed.getProblemSetup().getStructuralDynamicsSetup().setNSteps(4);
    {
        RSolver second(resumed,QString(),QString());
        second.setStructuralControlCallback("constant-force","conditioned-test-v1",provider);
        second.setStructuralControlStepObserver([&actual](const RStructuralControlStepRecord &record)
            { actual.push_back(record); });
        try { second.run(); }
        catch (const RError &failure) { QFAIL(qPrintable(failure.getMessage())); }
    }
    QCOMPARE(actual.size(),expected.size());
    for (size_t i=0;i<expected.size();i++)
    {
        QCOMPARE(actual[i].conditionedCommands[0].value,expected[i].conditionedCommands[0].value);
        QCOMPARE(actual[i].appliedChannels[0].command,expected[i].appliedChannels[0].command);
        QCOMPARE(actual[i].appliedChannels[0].cumulativeWork,
                 expected[i].appliedChannels[0].cumulativeWork);
    }
    QVERIFY(std::abs(nodalValue(resumed,R_VARIABLE_DISPLACEMENT,1)-
                     nodalValue(baseline,R_VARIABLE_DISPLACEMENT,1)) < 1.0e-12);
    QCOMPARE(RStructuralControlStateEnvelope::deserialize(resumed.getProblemSetup()
                 .getStructuralDynamicsRestartData().controlControllerState).controllerState,
             QByteArray("6"));
}

void TestStructuralDynamicsFE::rawLegacyControlStateFailsClosed()
{
    RModel model = makeTruss(1);
    configureControl(model);
    model.getProblemSetup().getStructuralDynamicsSetup().setNSteps(2);
    {
        RSolver first(model,QString(),QString());
        first.setStructuralControlCallback("constant-force","legacy-reject-test-v1",constantControl);
        first.run();
    }
    RStructuralDynamicsRestartData &checkpoint = model.getProblemSetup()
        .getStructuralDynamicsRestartData();
    checkpoint.controlControllerState = QByteArray("2"); // Pre-RCCE opaque layout.
    model.getProblemSetup().setRestart(true);
    model.getProblemSetup().getStructuralDynamicsSetup().setNSteps(1);
    RSolver resumed(model,QString(),QString());
    resumed.setStructuralControlCallback("constant-force","legacy-reject-test-v1",constantControl);
    bool rejected = false;
    try { resumed.run(); }
    catch (const RError &failure)
    {
        rejected = failure.getMessage().contains("composite",Qt::CaseInsensitive) ||
                   failure.getMessage().contains("control",Qt::CaseInsensitive);
    }
    QVERIFY(rejected);
    QCOMPARE(checkpoint.completedSteps,uint(2));
}

void TestStructuralDynamicsFE::twoChannelMappingIsLinearAndRestartsDeterministically()
{
    auto runWith = [](RModel &model, const std::vector<double> &forces,
                      std::vector<RStructuralControlStepRecord> *records = nullptr)
    {
        std::vector<QString> ids;
        for (const auto &actuator : model.getProblemSetup().getStructuralInteractionSetup().getActuators())
            ids.push_back(actuator.id);
        QCOMPARE(ids.size(),forces.size());
        RSolver solver(model,QString(),QString());
        solver.setStructuralControlCallback("constant-force","two-channel-test-v1",
            [ids,forces](const RStructuralControlObservation &)
            {
                RStructuralControlDecision decision;
                for (size_t i=0;i<ids.size();i++)
                    decision.requestedCommands.push_back({ids[i],forces[i]});
                return decision;
            });
        if (records)
            solver.setStructuralControlStepObserver([records](const RStructuralControlStepRecord &record)
                { records->push_back(record); });
        solver.run();
    };
    RModel passive = makeTruss(2);
    QString error;
    QVERIFY2(runSolver(passive,error),qPrintable(error));

    RModel tip = makeTruss(2);
    configureControl(tip);
    runWith(tip,{2.0e4});
    RModel middle = makeTruss(2);
    configureControl(middle);
    middle.getProblemSetup().getStructuralInteractionSetup().getActuators()[0].id = "mid-x";
    middle.getProblemSetup().getStructuralInteractionSetup().getActuators()[0].target.nodeId = 1;
    runWith(middle,{-1.0e4});

    auto makeBoth = []()
    {
        RModel model = makeTruss(2);
        configureControl(model);
        RStructuralActuator middleActuator = model.getProblemSetup()
            .getStructuralInteractionSetup().getActuators()[0];
        middleActuator.id = "mid-x";
        middleActuator.target.nodeId = 1;
        model.getProblemSetup().getStructuralInteractionSetup().getActuators().push_back(middleActuator);
        RStructuralSensor middleSensor = model.getProblemSetup()
            .getStructuralInteractionSetup().getSensors()[0];
        middleSensor.id = "mid-displacement";
        middleSensor.target.nodeId = 1;
        model.getProblemSetup().getStructuralInteractionSetup().getSensors().push_back(middleSensor);
        return model;
    };
    RModel both = makeBoth();
    std::vector<RStructuralControlStepRecord> records;
    runWith(both,{2.0e4,-1.0e4},&records);
    QCOMPARE(records.size(),size_t(40));
    QCOMPARE(records.front().sensors.size(),size_t(2));
    QCOMPARE(records.front().appliedChannels[0].actuatorId,QString("tip-x"));
    QCOMPARE(records.front().appliedChannels[1].actuatorId,QString("mid-x"));
    for (uint node : {1u,2u})
    {
        const double expected = nodalValue(tip,R_VARIABLE_DISPLACEMENT,node)+
                                nodalValue(middle,R_VARIABLE_DISPLACEMENT,node)-
                                nodalValue(passive,R_VARIABLE_DISPLACEMENT,node);
        QVERIFY(std::abs(nodalValue(both,R_VARIABLE_DISPLACEMENT,node)-expected) < 1.0e-10);
    }
    RModel resumed = makeBoth();
    resumed.getProblemSetup().getStructuralDynamicsSetup().setNSteps(13);
    runWith(resumed,{2.0e4,-1.0e4});
    resumed.getProblemSetup().setRestart(true);
    resumed.getProblemSetup().getStructuralDynamicsSetup().setNSteps(27);
    runWith(resumed,{2.0e4,-1.0e4});
    for (uint node : {1u,2u})
    {
        QVERIFY(std::abs(nodalValue(both,R_VARIABLE_DISPLACEMENT,node)-
                         nodalValue(resumed,R_VARIABLE_DISPLACEMENT,node)) < 1.0e-12);
        QVERIFY(std::abs(nodalValue(both,R_VARIABLE_VELOCITY,node)-
                         nodalValue(resumed,R_VARIABLE_VELOCITY,node)) < 1.0e-10);
    }
    const auto &a = both.getProblemSetup().getStructuralDynamicsRestartData();
    const auto &b = resumed.getProblemSetup().getStructuralDynamicsRestartData();
    QCOMPARE(a.cumulativeControlWork.size(),size_t(2));
    for (size_t i=0;i<2;i++)
        QVERIFY(std::abs(a.cumulativeControlWork[i]-b.cumulativeControlWork[i]) < 1.0e-9);
}

void TestStructuralDynamicsFE::rayleighDampingMatchesAnalyticAfterManySteps()
{
    constexpr double dampingRatio = 0.08;
    RModel model = makeTruss(1);
    const double omega = singleElementFrequency();
    // C=alpha*M, alpha=2*zeta*omega for this single-mode model.
    model.getProblemSetup().getStructuralDynamicsSetup()
         .setRayleighAlpha(2.0*dampingRatio*omega);
    QString error;
    QVERIFY2(runSolver(model,error),qPrintable(error));
    const double time = 40.0*model.getProblemSetup().getStructuralDynamicsSetup().getTimeStep();
    const double expected = dampedDisplacement(omega,dampingRatio,time);
    const double actual = nodalValue(model,R_VARIABLE_DISPLACEMENT,1);
    QVERIFY(std::isfinite(actual));
    QVERIFY2(std::abs(actual-expected) < 3.0e-4*initialDisplacement,
             qPrintable(QString("Damped FE q=%1, analytic q=%2").arg(actual).arg(expected)));
}

void TestStructuralDynamicsFE::modalDampingMatchesAnalyticAfterManySteps()
{
    constexpr double dampingRatio = 0.08;
    RModel model = makeTruss(1);
    model.getProblemSetup().getStructuralDynamicsSetup()
         .setModalDampingRatios({dampingRatio});
    QString error;
    QVERIFY2(runSolver(model,error),qPrintable(error));
    const double omega = singleElementFrequency();
    const double time = 40.0*model.getProblemSetup().getStructuralDynamicsSetup().getTimeStep();
    const double expected = dampedDisplacement(omega,dampingRatio,time);
    const double actual = nodalValue(model,R_VARIABLE_DISPLACEMENT,1);
    QVERIFY(std::isfinite(actual));
    QVERIFY2(std::abs(actual-expected) < 3.0e-4*initialDisplacement,
             qPrintable(QString("Modal-damped FE q=%1, analytic q=%2")
                        .arg(actual).arg(expected)));
    const auto &checkpoint = model.getProblemSetup().getStructuralDynamicsRestartData();
    QVERIFY(checkpoint.available);
    QVERIFY(checkpoint.cumulativeDampingDissipation > 0.0);
}

void TestStructuralDynamicsFE::twoModeDampingMatchesAnalyticProjection()
{
    // Two equal-length consistent-mass elements, with node 0 held, give
    // K=k[[2,-1],[-1,1]] and M=m[[4,1],[1,2]], where m=rho*A*Le/6.
    // The exact generalized eigenvalues are (5 +/- 3*sqrt(2))*k/(7*m).
    constexpr double lowDamping = 0.08;
    RModel model = makeTruss(2);
    const double elementLength = length/2.0;
    const double k = youngModulus*crossArea/elementLength;
    const double m = massDensity*crossArea*elementLength/6.0;
    const double root = std::sqrt(2.0);
    const double eigenvalues[2] = {(5.0-3.0*root)*k/(7.0*m),
                                    (5.0+3.0*root)*k/(7.0*m)};
    auto &setup = model.getProblemSetup().getStructuralDynamicsSetup();
    setup.setModalDampingRatios({lowDamping,0.0});
    setup.setTimeStep(0.01/std::sqrt(eigenvalues[1]));
    setup.setNSteps(120);
    setup.setOutputFrequency(120);

    QString error;
    QVERIFY2(runSolver(model,error),qPrintable(error));
    const double time = 120.0*setup.getTimeStep();
    double expected[2] = {0.0,0.0};
    for (uint mode=0;mode<2;mode++)
    {
        const double lambda = eigenvalues[mode];
        const double shape2 = (2.0*k-4.0*m*lambda)/(k+m*lambda);
        const double massNorm = std::sqrt(m*(4.0+2.0*shape2+2.0*shape2*shape2));
        const double phi[2] = {1.0/massNorm,shape2/massNorm};
        // q(0)=[0,u0], v(0)=0; modal coordinate is phi^T*M*q(0).
        const double amplitude = m*initialDisplacement*(phi[0]+2.0*phi[1]);
        const double zeta = mode == 0 ? lowDamping : 0.0;
        const double omega = std::sqrt(lambda);
        const double dampedOmega = omega*std::sqrt(1.0-zeta*zeta);
        const double response = std::exp(-zeta*omega*time)*
            (std::cos(dampedOmega*time)+zeta*omega/dampedOmega*std::sin(dampedOmega*time));
        for (uint node=0;node<2;node++) expected[node] += phi[node]*amplitude*response;
    }
    for (uint node=1;node<=2;node++)
    {
        const double actual = nodalValue(model,R_VARIABLE_DISPLACEMENT,node);
        QVERIFY2(std::isfinite(actual) &&
                 std::abs(actual-expected[node-1]) < 5.0e-4*initialDisplacement,
                 qPrintable(QString("Two-mode FE node %1: q=%2, analytic=%3")
                            .arg(node).arg(actual).arg(expected[node-1])));
    }
}

void TestStructuralDynamicsFE::modalDampingWorksWithHomogeneousMpc()
{
    RModel model = makeTruss(2,true);
    constexpr double dampingRatio = 0.08;
    model.getProblemSetup().getStructuralDynamicsSetup()
         .setModalDampingRatios({dampingRatio});
    QString error;
    QVERIFY2(runSolver(model,error),qPrintable(error));
    // T=[1,1]^T maps the tied physical DOFs onto one reduced DOF. This
    // yields K_red=k and M_red=8m, hence the one-element frequency below.
    const double time = 40.0*model.getProblemSetup().getStructuralDynamicsSetup().getTimeStep();
    const double expected = dampedDisplacement(singleElementFrequency(),dampingRatio,time);
    for (uint node=1;node<=2;node++)
    {
        const double actual = nodalValue(model,R_VARIABLE_DISPLACEMENT,node);
        QVERIFY2(std::isfinite(actual) && std::abs(actual-expected) < 3.0e-4*initialDisplacement,
                 qPrintable(QString("Modal MPC node %1: q=%2, analytic=%3")
                            .arg(node).arg(actual).arg(expected)));
    }
}

void TestStructuralDynamicsFE::unsupportedModalDampingCombinationsFailBeforeFirstStep()
{
    RModel mixed = makeTruss(1);
    auto &mixedSetup = mixed.getProblemSetup().getStructuralDynamicsSetup();
    mixedSetup.setModalDampingRatios({0.05});
    mixedSetup.setRayleighAlpha(1.0);
    QString error;
    QVERIFY2(!runSolver(mixed,error),"Modal and Rayleigh damping were silently combined.");
    QVERIFY2(error.contains("Modal and Rayleigh",Qt::CaseInsensitive),qPrintable(error));
    QVERIFY(!mixed.getProblemSetup().getStructuralDynamicsRestartData().available);

    RModel moving = makeTruss(1);
    auto &movingSetup = moving.getProblemSetup().getStructuralDynamicsSetup();
    movingSetup.setModalDampingRatios({0.05});
    const double duration = movingSetup.getNSteps()*movingSetup.getTimeStep();
    RStructuralPrescribedMotion motion;
    motion.id = "modal-motion-unsupported";
    motion.target = {0,RDofType::UX};
    motion.knots = {{0.0,0.0,0.0,0.0},{duration,0.0,0.0,0.0}};
    moving.getProblemSetup().getStructuralInteractionSetup()
          .getPrescribedMotions().push_back(motion);
    QVERIFY2(!runSolver(moving,error),"Modal damping with moving support was accepted.");
    QVERIFY2(error.contains("moving support",Qt::CaseInsensitive),qPrintable(error));
    QVERIFY(!moving.getProblemSetup().getStructuralDynamicsRestartData().available);
}

void TestStructuralDynamicsFE::stationarySupportReactionIncludesInertia()
{
    RModel model = makeTruss(1);
    std::vector<std::vector<RStructuralSupportSample>> history;
    std::vector<RStructuralSupportSample> last;
    QString error;
    QVERIFY2(runSolver(model,error,&history,&last),qPrintable(error));
    QCOMPARE(history.size(),std::size_t(40));
    const RStructuralSupportSample *support = findSupport(last,0,RDofType::UX);
    QVERIFY(support != nullptr);
    QVERIFY(support->motionId.isEmpty());
    const double k = youngModulus*crossArea/length;
    const double q = nodalValue(model,R_VARIABLE_DISPLACEMENT,1);
    // M_cf=M_ff/2, a_f=-k*q/M_ff; K_cf=-k.
    const double expected = -1.5*k*q;
    QVERIFY2(std::abs(support->reaction-expected) < 1.0e-7*k*initialDisplacement,
             qPrintable(QString("Support-on-structure R=%1, M*a+K*u=%2")
                        .arg(support->reaction).arg(expected)));
    QCOMPARE(support->velocity,0.0);
    QCOMPARE(support->power,0.0);
    QCOMPARE(support->cumulativeWork,0.0);
}

void TestStructuralDynamicsFE::movingSupportReactionPowerAndAppliedSupportForce()
{
    constexpr double speed = 0.1;
    constexpr double appliedSupportForce = 1234.567;
    auto makeMoving = [speed](double appliedForce)
    {
        RModel model = makeTruss(1);
        model.getPoint(1).removeAllInitialConditions();
        const auto &setup = model.getProblemSetup().getStructuralDynamicsSetup();
        const double duration = setup.getNSteps()*setup.getTimeStep();
        RStructuralPrescribedMotion motion;
        motion.id = "linear-base-x";
        motion.target = {0,RDofType::UX};
        motion.knots = {{0.0,0.0,speed,0.0},{duration,speed*duration,speed,0.0}};
        model.getProblemSetup().getStructuralInteractionSetup()
             .getPrescribedMotions().push_back(motion);
        if (appliedForce != 0.0)
        {
            RBoundaryCondition force(R_BOUNDARY_CONDITION_FORCE);
            force.getComponent(force.findComponentPosition(R_VARIABLE_FORCE_X))
                 .add(0.0,appliedForce);
            model.getPoint(0).addBoundaryCondition(force);
        }
        return model;
    };
    RModel model = makeMoving(0.0);
    RModel loaded = makeMoving(appliedSupportForce);
    std::vector<RStructuralSupportSample> last;
    std::vector<RStructuralSupportSample> loadedLast;
    QString error;
    QVERIFY2(runSolver(model,error,nullptr,&last),qPrintable(error));
    QVERIFY2(runSolver(loaded,error,nullptr,&loadedLast),qPrintable(error));
    const RStructuralSupportSample *support = findSupport(last,0,RDofType::UX);
    const RStructuralSupportSample *loadedSupport = findSupport(loadedLast,0,RDofType::UX);
    QVERIFY(support != nullptr && loadedSupport != nullptr);
    QCOMPARE(support->motionId,QString("linear-base-x"));
    const double omega = singleElementFrequency();
    const double time = 40.0*model.getProblemSetup().getStructuralDynamicsSetup().getTimeStep();
    const double totalMass = massDensity*crossArea*length;
    const double expectedReaction = 0.5*totalMass*speed*omega*std::sin(omega*time);
    const double expectedWork = 0.5*totalMass*speed*speed*(1.0-std::cos(omega*time));
    QVERIFY2(std::abs(support->reaction-expectedReaction) < 5.0e-4*expectedReaction,
             qPrintable(QString("Moving support R=%1, analytic=%2")
                        .arg(support->reaction).arg(expectedReaction)));
    QVERIFY(std::abs(support->power-support->reaction*speed) < 1.0e-10*expectedReaction*speed);
    QVERIFY(std::abs(support->cumulativeWork-expectedWork) < 1.0e-3*expectedWork);
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_DISPLACEMENT,1)-
                     nodalValue(loaded,R_VARIABLE_DISPLACEMENT,1)) < 1.0e-12);
    QVERIFY(std::abs((loadedSupport->reaction-support->reaction)+appliedSupportForce) < 1.0e-7);
    QVERIFY(std::abs((loadedSupport->power-support->power)+appliedSupportForce*speed) < 1.0e-8);
    QVERIFY(std::abs((loadedSupport->cumulativeWork-support->cumulativeWork)+
                     appliedSupportForce*speed*time) < 1.0e-6);
}

void TestStructuralDynamicsFE::supportGetterRequiresObserver()
{
    RModel model = makeTruss(1);
    RSolver solver(model,QString(),QString());
    QVERIFY_EXCEPTION_THROWN(solver.getStructuralSupportSamples(),RError);
}

void TestStructuralDynamicsFE::supportObserverExceptionsDoNotRevokeAcceptedStep()
{
    RModel model = makeTruss(1);
    model.getProblemSetup().getStructuralDynamicsSetup().setNSteps(1);
    RSolver solver(model,QString(),QString());
    uint notifications = 0;
    solver.setStructuralSupportResponseObserver(
        [&notifications](const std::vector<RStructuralSupportSample> &)
        {
            notifications++;
            throw std::runtime_error("consumer failure after support publication");
        });
    solver.run();
    QCOMPARE(notifications,uint(1));
    QCOMPARE(model.getProblemSetup().getStructuralDynamicsRestartData().completedSteps,uint(1));
    QVERIFY(model.getProblemSetup().getStructuralDynamicsRestartData().supportResponseAvailable);
    const auto samples = solver.getStructuralSupportSamples();
    QVERIFY(findSupport(samples,0,RDofType::UX) != nullptr);
}

void TestStructuralDynamicsFE::supportObserverRejectsModalAndMpc()
{
    RModel modal = makeTruss(1);
    modal.getProblemSetup().getStructuralDynamicsSetup().setModalDampingRatios({0.05});
    QString error;
    std::vector<std::vector<RStructuralSupportSample>> history;
    QVERIFY2(!runSolver(modal,error,&history),"Modal support reaction sampling was accepted.");
    QVERIFY2(error.contains("modal damping",Qt::CaseInsensitive),qPrintable(error));
    QVERIFY(history.empty());

    RModel mpc = makeTruss(2,true);
    QVERIFY2(!runSolver(mpc,error,&history),"MPC support reaction sampling was accepted.");
    QVERIFY2(error.contains("MPC",Qt::CaseInsensitive),qPrintable(error));
    QVERIFY(history.empty());
}

void TestStructuralDynamicsFE::movingSupportMatchesConsistentMassAnalytic()
{
    constexpr double baseAcceleration = 1000.0;
    RModel model = makeTruss(1);
    model.getPoint(1).removeAllInitialConditions();
    const RStructuralDynamicsSetup &setup = model.getProblemSetup().getStructuralDynamicsSetup();
    const double duration = setup.getNSteps()*setup.getTimeStep();
    RStructuralPrescribedMotion motion;
    motion.id = "quadratic-base";
    motion.target = {0,RDofType::UX};
    motion.knots = {{0.0,0.0,0.0,baseAcceleration},
                    {duration,0.5*baseAcceleration*duration*duration,
                     baseAcceleration*duration,baseAcceleration}};
    model.getProblemSetup().getStructuralInteractionSetup()
         .getPrescribedMotions().push_back(motion);

    QString error;
    QVERIFY2(runSolver(model,error),qPrintable(error));
    const double omega = singleElementFrequency();
    const double stiffness = youngModulus*crossArea/length;
    const double freeMass = massDensity*crossArea*length/3.0;
    const double couplingMass = massDensity*crossArea*length/6.0;
    const double phase = omega*duration;
    const double expectedQ = 0.5*baseAcceleration*duration*duration
        -(freeMass+couplingMass)*baseAcceleration/stiffness*(1.0-std::cos(phase));
    const double expectedV = baseAcceleration*duration
        -(freeMass+couplingMass)*baseAcceleration/stiffness*omega*std::sin(phase);
    const double expectedA = baseAcceleration
        -(freeMass+couplingMass)*baseAcceleration/freeMass*std::cos(phase);
    const double q = nodalValue(model,R_VARIABLE_DISPLACEMENT,1);
    const double v = nodalValue(model,R_VARIABLE_VELOCITY,1);
    const double a = nodalValue(model,R_VARIABLE_ACCELERATION,1);
    QVERIFY2(std::abs(q-expectedQ) < 4.0e-8,
             qPrintable(QString("moving support q=%1 exact=%2").arg(q,0,'g',16).arg(expectedQ,0,'g',16)));
    QVERIFY2(std::abs(v-expectedV) < 2.0e-4,
             qPrintable(QString("moving support v=%1 exact=%2").arg(v,0,'g',16).arg(expectedV,0,'g',16)));
    QVERIFY2(std::abs(a-expectedA) < 0.5,
             qPrintable(QString("moving support a=%1 exact=%2").arg(a,0,'g',16).arg(expectedA,0,'g',16)));
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_DISPLACEMENT,0)
                     -0.5*baseAcceleration*duration*duration) < 1.0e-11);
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_VELOCITY,0)
                     -baseAcceleration*duration) < 1.0e-9);
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_ACCELERATION,0)
                     -baseAcceleration) < 1.0e-6);
}

void TestStructuralDynamicsFE::movingSupportRayleighMassDampingMatchesAnalytic()
{
    constexpr double baseVelocity = 1.0;
    RModel model = makeTruss(1);
    model.getPoint(1).removeAllInitialConditions();
    RStructuralDynamicsSetup &setup = model.getProblemSetup().getStructuralDynamicsSetup();
    const double omega = singleElementFrequency();
    const double alpha = 0.1*omega;
    setup.setRayleighAlpha(alpha);
    const double duration = setup.getNSteps()*setup.getTimeStep();
    RStructuralPrescribedMotion motion;
    motion.id = "linear-base-with-mass-damping";
    motion.target = {0,RDofType::UX};
    motion.knots = {{0.0,0.0,baseVelocity,0.0},
                    {duration,baseVelocity*duration,baseVelocity,0.0}};
    model.getProblemSetup().getStructuralInteractionSetup()
         .getPrescribedMotions().push_back(motion);

    QString error;
    QVERIFY2(runSolver(model,error),qPrintable(error));
    const double freeMass = massDensity*crossArea*length/3.0;
    const double couplingMass = massDensity*crossArea*length/6.0;
    const double stiffness = youngModulus*crossArea/length;
    const double zeta = alpha/(2.0*omega);
    const double dampedOmega = omega*std::sqrt(1.0-zeta*zeta);
    const double relativeSteady = -alpha*baseVelocity*(freeMass+couplingMass)/stiffness;
    const double cosine = -relativeSteady;
    const double sine = (-baseVelocity+zeta*omega*cosine)/dampedOmega;
    const double expected = baseVelocity*duration+relativeSteady
        +std::exp(-zeta*omega*duration)
         *(cosine*std::cos(dampedOmega*duration)+sine*std::sin(dampedOmega*duration));
    const double actual = nodalValue(model,R_VARIABLE_DISPLACEMENT,1);
    QVERIFY2(std::abs(actual-expected) < 4.0e-8,
             qPrintable(QString("damped moving support q=%1 exact=%2")
                        .arg(actual,0,'g',16).arg(expected,0,'g',16)));
}

void TestStructuralDynamicsFE::unsupportedMovingSupportInputsFailBeforeFirstStep()
{
    auto addZeroMotion = [](RModel &model)
    {
        const auto &dynamics = model.getProblemSetup().getStructuralDynamicsSetup();
        const double duration = dynamics.getNSteps()*dynamics.getTimeStep();
        RStructuralPrescribedMotion motion;
        motion.id = "unsupported-base";
        motion.target = {0,RDofType::UX};
        motion.knots = {{0.0,0.0,0.0,0.0},{duration,0.0,0.0,0.0}};
        model.getProblemSetup().getStructuralInteractionSetup()
             .getPrescribedMotions().push_back(motion);
    };
    QString error;
    RModel withMpc = makeTruss(2,true);
    addZeroMotion(withMpc);
    QVERIFY2(!runSolver(withMpc,error),"MPC and moving support were silently combined.");
    QVERIFY2(error.contains("MPC/coupling",Qt::CaseInsensitive),qPrintable(error));

    RModel delayedBoundary = makeTruss(1);
    addZeroMotion(delayedBoundary);
    const double duration = delayedBoundary.getProblemSetup()
        .getStructuralDynamicsSetup().getNSteps()*delayedBoundary.getProblemSetup()
        .getStructuralDynamicsSetup().getTimeStep();
    auto &component = delayedBoundary.getPoint(0).getBoundaryCondition(0).getComponent(
        delayedBoundary.getPoint(0).getBoundaryCondition(0)
            .findComponentPosition(R_VARIABLE_DISPLACEMENT_X));
    component.add(duration/2.0,0.0);
    component.add(duration,1.0e-5);
    QVERIFY2(!runSolver(delayedBoundary,error),"A delayed nonzero baseline BC was accepted.");
    QVERIFY2(error.contains("baseline BC knot",Qt::CaseInsensitive),qPrintable(error));

    RModel nearZeroBoundary = makeTruss(1);
    addZeroMotion(nearZeroBoundary);
    auto &nearZeroComponent = nearZeroBoundary.getPoint(0).getBoundaryCondition(0).getComponent(
        nearZeroBoundary.getPoint(0).getBoundaryCondition(0)
            .findComponentPosition(R_VARIABLE_DISPLACEMENT_X));
    nearZeroComponent.add(0.0,5.0e-13);
    QVERIFY2(!runSolver(nearZeroBoundary,error),
             "A small nonzero baseline BC was silently added to the moving support force.");
    QVERIFY2(error.contains("baseline",Qt::CaseInsensitive),qPrintable(error));
    QVERIFY2(!nearZeroBoundary.getProblemSetup().getStructuralDynamicsRestartData().available,
             "A nonzero baseline BC must fail before accepting the first step.");

    RModel shortDomain = makeTruss(1);
    addZeroMotion(shortDomain);
    shortDomain.getProblemSetup().getStructuralInteractionSetup()
        .getPrescribedMotions().front().knots.back().time = duration/2.0;
    QVERIFY2(!runSolver(shortDomain,error),"A motion knot domain shorter than the solve was accepted.");
    QVERIFY2(error.contains("outside its time domain",Qt::CaseInsensitive),qPrintable(error));

    RModel nearEndDomain = makeTruss(1);
    addZeroMotion(nearEndDomain);
    // 1e-13 s is inside the old solver's 1e-12 preflight window but far
    // outside the model evaluator's 16-ulp endpoint rounding tolerance.
    nearEndDomain.getProblemSetup().getStructuralInteractionSetup()
        .getPrescribedMotions().front().knots.back().time = duration-1.0e-13;
    QVERIFY2(!runSolver(nearEndDomain,error),
             "A slightly short motion domain advanced until the final step.");
    QVERIFY2(error.contains("outside its time domain",Qt::CaseInsensitive),qPrintable(error));
    QVERIFY2(!nearEndDomain.getProblemSetup().getStructuralDynamicsRestartData().available,
             "A motion domain failure must happen before accepting the first step.");

    RModel localFrame = makeTruss(1);
    addZeroMotion(localFrame);
    localFrame.getPoint(0).removeAllBoundaryConditions();
    RBoundaryCondition normal(R_BOUNDARY_CONDITION_DISPLACEMENT_NORMAL);
    normal.setLocalDirection(RR3Vector(1.0,0.0,0.0));
    normal.getComponent(normal.findComponentPosition(R_VARIABLE_DISPLACEMENT)).add(0.0,0.0);
    localFrame.getPoint(0).addBoundaryCondition(normal);
    QVERIFY2(!runSolver(localFrame,error),"A local-frame motion was accepted.");
    QVERIFY2(error.contains("local",Qt::CaseInsensitive),qPrintable(error));
}

void TestStructuralDynamicsFE::tiedNodesStayTiedDuringDynamics()
{
    RModel model = makeTruss(2,true);
    model.getProblemSetup().getStructuralDynamicsSetup().setNSteps(12);
    QString error;
    QVERIFY2(runSolver(model,error),qPrintable(error));
    const double middle = nodalValue(model,R_VARIABLE_DISPLACEMENT,1);
    const double tip = nodalValue(model,R_VARIABLE_DISPLACEMENT,2);
    const double middleVelocity = nodalValue(model,R_VARIABLE_VELOCITY,1);
    const double tipVelocity = nodalValue(model,R_VARIABLE_VELOCITY,2);
    QVERIFY(std::isfinite(middle) && std::isfinite(tip));
    QVERIFY(std::isfinite(middleVelocity) && std::isfinite(tipVelocity));
    QVERIFY(std::abs(tip-middle) < 1.0e-9*initialDisplacement);
    QVERIFY(std::abs(tipVelocity-middleVelocity)
            < 1.0e-9*initialDisplacement*singleElementFrequency());
    QVERIFY(std::abs(tip) > 0.5*initialDisplacement);
}

void TestStructuralDynamicsFE::fixedThreadCountsProduceSameState()
{
    RModel serial = makeTruss(6);
    RModel parallel = makeTruss(6);
    for (RModel *model : {&serial,&parallel})
    {
        model->getProblemSetup().getStructuralDynamicsSetup().setNSteps(8);
    }
    auto &serialPolicy = serial.getProblemSetup().getParallelExecutionConfig();
    serialPolicy.setMode(RParallelExecutionMode::FixedThreads);
    serialPolicy.setRequestedThreads(1);
    auto &parallelPolicy = parallel.getProblemSetup().getParallelExecutionConfig();
    parallelPolicy.setMode(RParallelExecutionMode::FixedThreads);
    parallelPolicy.setRequestedThreads(2);

    QString error;
    QVERIFY2(runSolver(serial,error),qPrintable(error));
    error.clear();
    QVERIFY2(runSolver(parallel,error),qPrintable(error));
    for (uint node=0;node<serial.getNNodes();node++)
    {
        for (RVariableType variable : {R_VARIABLE_DISPLACEMENT,
                                       R_VARIABLE_VELOCITY,
                                       R_VARIABLE_ACCELERATION})
        {
            const double oneThread = nodalValue(serial,variable,node);
            const double twoThreads = nodalValue(parallel,variable,node);
            QVERIFY(std::isfinite(oneThread) && std::isfinite(twoThreads));
            const double scale = std::max({1.0e-10,std::abs(oneThread),std::abs(twoThreads)});
            QVERIFY2(std::abs(oneThread-twoThreads) < 1.0e-9*scale+1.0e-13,
                     qPrintable(QString("Parallel result mismatch at node %1, variable %2: %3 vs %4")
                                .arg(node).arg(int(variable)).arg(oneThread).arg(twoThreads)));
        }
    }
}

void TestStructuralDynamicsFE::restartMatchesUninterruptedRun()
{
    RModel uninterrupted = makeTruss(2,true);
    RModel resumed = makeTruss(2,true);
    QString error;
    QVERIFY2(runSolver(uninterrupted,error),qPrintable(error));
    resumed.getProblemSetup().getStructuralDynamicsSetup().setNSteps(17);
    QVERIFY2(runSolver(resumed,error),qPrintable(error));
    const auto &saved = resumed.getProblemSetup().getStructuralDynamicsRestartData();
    QVERIFY(saved.available);
    QCOMPARE(saved.completedSteps,uint(17));
    QCOMPARE(saved.q.size(),saved.v.size());
    resumed.getProblemSetup().setRestart(true);
    resumed.getProblemSetup().getStructuralDynamicsSetup().setNSteps(23);
    error.clear();
    QVERIFY2(runSolver(resumed,error),qPrintable(error));
    QCOMPARE(resumed.getProblemSetup().getStructuralDynamicsRestartData().completedSteps,uint(40));
    for (uint node=0;node<resumed.getNNodes();node++)
    {
        for (RVariableType variable : {R_VARIABLE_DISPLACEMENT,R_VARIABLE_VELOCITY,R_VARIABLE_ACCELERATION})
        {
            const double expected = nodalValue(uninterrupted,variable,node);
            const double actual = nodalValue(resumed,variable,node);
            const double scale = std::max({1.0e-10,std::fabs(expected),std::fabs(actual)});
            QVERIFY2(std::fabs(actual-expected) < 1.0e-8*scale+1.0e-13,
                     qPrintable(QString("Restart mismatch at node %1 variable %2: %3 vs %4")
                                .arg(node).arg(int(variable)).arg(actual).arg(expected)));
        }
    }
}

void TestStructuralDynamicsFE::restartRejectsChangedMeshAndTimeStep()
{
    RModel model = makeTruss(1);
    model.getProblemSetup().getStructuralDynamicsSetup().setNSteps(2);
    QString error;
    QVERIFY2(runSolver(model,error),qPrintable(error));
    model.getProblemSetup().setRestart(true);
    model.getProblemSetup().getStructuralDynamicsSetup().setNSteps(1);
    const double oldStep = model.getProblemSetup().getStructuralDynamicsSetup().getTimeStep();
    model.getProblemSetup().getStructuralDynamicsSetup().setTimeStep(oldStep*1.1);
    QVERIFY2(!runSolver(model,error),"Restart with changed time step unexpectedly succeeded");
    QVERIFY(error.contains("checkpoint",Qt::CaseInsensitive));
    model.getProblemSetup().getStructuralDynamicsSetup().setTimeStep(oldStep);
    model.getNode(1).set(length+0.01,0.0,0.0);
    error.clear();
    QVERIFY2(!runSolver(model,error),"Restart with changed mesh unexpectedly succeeded");
    QVERIFY(error.contains("checkpoint",Qt::CaseInsensitive));
}

void TestStructuralDynamicsFE::restartRejectsChangedOperatorAndMpc()
{
    RModel saved = makeTruss(2,true);
    saved.getProblemSetup().getStructuralDynamicsSetup().setNSteps(2);
    QString error;
    QVERIFY2(runSolver(saved,error),qPrintable(error));
    RModel changedDamping(saved);
    changedDamping.getProblemSetup().setRestart(true);
    changedDamping.getProblemSetup().getStructuralDynamicsSetup().setNSteps(1);
    changedDamping.getProblemSetup().getStructuralDynamicsSetup().setRayleighAlpha(1.0);
    error.clear();
    QVERIFY2(!runSolver(changedDamping,error),"Restart with changed damping unexpectedly succeeded");
    QVERIFY(error.contains("operator",Qt::CaseInsensitive));

    RModel changedMpc(saved);
    changedMpc.getProblemSetup().setRestart(true);
    changedMpc.getProblemSetup().getStructuralDynamicsSetup().setNSteps(1);
    changedMpc.getProblemSetup().getStructuralInteractionSetup()
              .getLinearConstraints()[0].terms[1].coefficient = -0.5;
    error.clear();
    QVERIFY2(!runSolver(changedMpc,error),"Restart with changed MPC unexpectedly succeeded");
    QVERIFY(error.contains("operator",Qt::CaseInsensitive));
}

QTEST_APPLESS_MAIN(TestStructuralDynamicsFE)

#include "tst_structural_dynamics_fe.moc"

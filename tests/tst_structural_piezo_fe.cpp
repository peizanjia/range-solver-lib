#include <QtTest>
#include <QDir>
#include <QElapsedTimer>
#include <QFileInfo>
#include <QProcess>
#include <QTemporaryDir>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <stdexcept>
#include <vector>

#include <rbl_error.h>

#include <rml_boundary_condition.h>
#include <rml_line.h>
#include <rml_model.h>
#include <rml_problem_setup.h>
#include <rml_problem_task_item.h>
#include <rml_structural_cable_setup.h>
#include <rml_structural_dynamics_restart_data.h>
#include <rml_structural_piezo_setup.h>
#include <rml_variable.h>

#include "rsolver.h"
#include "rstructural_piezo_model_adapter.h"

namespace
{
void setComponent(RBoundaryCondition &condition,
                  RVariableType type,
                  bool enabled,
                  double value = 0.0)
{
    const uint position = condition.findComponentPosition(type);
    if (position == RConstants::eod)
        throw std::logic_error("The boundary condition is missing a required component.");
    RConditionComponent &component = condition.getComponent(position);
    component.setEnabled(enabled);
    component.clear();
    if (enabled) component.add(0.0,value);
}

RBoundaryCondition displacementCondition(const std::array<bool,3> &enabled,
                                         const std::array<double,3> &values)
{
    RBoundaryCondition condition(R_BOUNDARY_CONDITION_DISPLACEMENT);
    const RVariableType types[3] = {
        R_VARIABLE_DISPLACEMENT_X,R_VARIABLE_DISPLACEMENT_Y,R_VARIABLE_DISPLACEMENT_Z
    };
    for (uint axis=0;axis<3;axis++) setComponent(condition,types[axis],enabled[axis],values[axis]);
    return condition;
}

RBoundaryCondition forceCondition(const std::array<double,3> &values)
{
    RBoundaryCondition condition(R_BOUNDARY_CONDITION_FORCE);
    const RVariableType types[3] = {R_VARIABLE_FORCE_X,R_VARIABLE_FORCE_Y,R_VARIABLE_FORCE_Z};
    for (uint axis=0;axis<3;axis++) setComponent(condition,types[axis],true,values[axis]);
    return condition;
}

void setEnvironmentValue(RElementGroup &group,
                         REnvironmentConditionType conditionType,
                         RVariableType variableType,
                         double value)
{
    REnvironmentCondition &condition = group.getEnvironmentCondition(conditionType);
    const uint position = condition.findComponentPosition(variableType);
    if (position == RConstants::eod)
        throw std::logic_error("The environment condition is missing a required component.");
    RConditionComponent &component = condition.getComponent(position);
    component.setEnabled(true);
    component.clear();
    component.add(0.0,value);
}

void setZeroGravity(RModel &model)
{
    const RVariableType gravityTypes[3] = {
        R_VARIABLE_G_ACCELERATION_X,R_VARIABLE_G_ACCELERATION_Y,R_VARIABLE_G_ACCELERATION_Z
    };
    for (uint groupId=0;groupId<model.getNElementGroups();groupId++)
    {
        RElementGroup *group = model.getElementGroupPtr(groupId);
        if (!group) continue;
        for (RVariableType type : gravityTypes)
            setEnvironmentValue(*group,R_ENVIRONMENT_CONDITION_G_ACCELERATION,type,0.0);
    }
}

void addPointElement(RModel &model,uint node,uint pointGroup)
{
    RElement point(R_ELEMENT_POINT);
    point.setNodeId(0,node);
    model.addElement(point,true,pointGroup);
}

void assignRodMaterial(RLine &line,double area,double youngModulus)
{
    line.setCrossArea(area);
    RMaterial material(RMaterial::Solid);
    RMaterialProperty modulus(RMaterialProperty::ModulusOfElasiticity);
    modulus.add(293.15,youngModulus);
    material.add(modulus);
    line.setMaterial(material);
}

RModel makeRod(double voltage,
               double coupling,
               double permittivity,
               bool reversedEndpoints = false,
               bool prescribeRightAxial = false,
               double rightAxialDisplacement = 0.0,
               double rightAxialForce = 0.0,
               uint segmentCount = 1,
               bool onlyFirstPort = false)
{
    constexpr double segmentLength = 2.0;
    RModel model;
    model.setNNodes(segmentCount+1);
    for (uint node=0;node<=segmentCount;node++)
        model.getNode(node).set(segmentLength*node,0.0,0.0);

    for (uint segment=0;segment<segmentCount;segment++)
    {
        RElement rod(R_ELEMENT_TRUSS1);
        const uint left = reversedEndpoints ? segment+1 : segment;
        const uint right = reversedEndpoints ? segment : segment+1;
        rod.setNodeId(0,left);
        rod.setNodeId(1,right);
        model.addElement(rod,true,0);
    }
    assignRodMaterial(model.getLine(0),0.5,10.0);
    for (uint node=0;node<=segmentCount;node++) addPointElement(model,node,node);

    model.getPoint(0).addBoundaryCondition(
        displacementCondition({{true,true,true}},{{0.0,0.0,0.0}}));
    for (uint node=1;node<segmentCount;node++)
        model.getPoint(node).addBoundaryCondition(
            displacementCondition({{false,true,true}},{{0.0,0.0,0.0}}));
    model.getPoint(segmentCount).addBoundaryCondition(displacementCondition(
        {{prescribeRightAxial,true,true}},
        {{rightAxialDisplacement,0.0,0.0}}));
    if (rightAxialForce != 0.0)
        model.getPoint(segmentCount).addBoundaryCondition(
            forceCondition({{rightAxialForce,0.0,0.0}}));

    setZeroGravity(model);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    RStructuralPiezoSetup &setup = model.getProblemSetup().getStructuralPiezoSetup();
    setup.enabled = true;
    const uint portCount = onlyFirstPort ? std::min(1u,segmentCount) : segmentCount;
    for (uint segment=0;segment<portCount;segment++)
    {
        RStructuralPiezoPort port;
        port.id = QString("port-%1").arg(segment);
        port.elementId = segment;
        port.coupling = reversedEndpoints ? -coupling : coupling;
        port.permittivity = permittivity;
        port.voltage = reversedEndpoints ? -voltage : voltage;
        setup.ports.push_back(port);
    }
    setup.bindToModel(model);
    return model;
}

// Independent static fixtures: one mechanical material group per bar. Passive
// materials contain only E and A, without fictitious electrical properties.
RModel makeMixedNetwork(bool triangle)
{
    RModel model;
    model.setNNodes(3);
    model.getNode(0).set(0.0,0.0,0.0);
    model.getNode(1).set(2.0,0.0,0.0);
    model.getNode(2).set(triangle ? 0.0 : 4.0,triangle ? 2.0 : 0.0,0.0);
    const uint count = triangle ? 3 : 2;
    const std::array<double,3> areas = {{0.5,0.8,0.3}};
    const std::array<double,3> moduli = {{10.0,20.0,30.0}};
    for (uint member=0;member<count;member++)
    {
        RElement rod(R_ELEMENT_TRUSS1);
        rod.setNodeId(0,member);
        rod.setNodeId(1,(member+1)%3);
        model.addElement(rod,true,member);
        assignRodMaterial(model.getLine(member),areas[member],moduli[member]);
    }
    for (uint node=0;node<3;node++) addPointElement(model,node,node);
    model.getPoint(0).addBoundaryCondition(
        displacementCondition({{true,true,true}},{{0.0,0.0,0.0}}));
    model.getPoint(1).addBoundaryCondition(
        displacementCondition({{false,true,true}},{{0.0,0.0,0.0}}));
    model.getPoint(2).addBoundaryCondition(
        displacementCondition({{false,!triangle,true}},{{0.0,0.0,0.0}}));
    model.getPoint(2).addBoundaryCondition(
        forceCondition(triangle ? std::array<double,3>{{0.04,-0.02,0.0}}
                                : std::array<double,3>{{0.025,0.0,0.0}}));
    setZeroGravity(model);
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    auto &setup = model.getProblemSetup().getStructuralPiezoSetup();
    setup.enabled = true;
    RStructuralPiezoPort port;
    port.id = "active-base";
    port.elementId = 0;
    port.coupling = 2.0;
    port.permittivity = 4.0;
    port.voltage = 0.03;
    setup.ports.push_back(port);
    for (uint member=1;member<count;member++) setup.passiveElementIds.push_back(member);
    setup.bindToModel(model);
    return model;
}

double nodalValue(const RModel &model,RVariableType type,uint component,uint node)
{
    const uint position = model.findVariable(type);
    if (position == RConstants::eod) return std::numeric_limits<double>::quiet_NaN();
    return model.getVariable(position).getValue(component,node);
}

double elementValue(const RModel &model,RVariableType type,uint element)
{
    const uint position = model.findVariable(type);
    if (position == RConstants::eod) return std::numeric_limits<double>::quiet_NaN();
    return model.getVariable(position).getValue(0,element);
}

bool sameResultVariables(const RModel &first,const RModel &second)
{
    if (first.getNVariables() != second.getNVariables()) return false;
    for (uint variable=0;variable<first.getNVariables();variable++)
    {
        const auto &a = first.getVariable(variable);
        const auto &b = second.getVariable(variable);
        if (a.getType() != b.getType() || a.getApplyType() != b.getApplyType() ||
            a.getNVectors() != b.getNVectors() || a.getNValues() != b.getNValues()) return false;
        for (uint component=0;component<a.getNVectors();component++)
            for (uint value=0;value<a.getNValues();value++)
                if (a.getValue(component,value) != b.getValue(component,value)) return false;
    }
    return true;
}

struct NativeSolverRun
{
    bool started = false;
    bool finished = false;
    QProcess::ExitStatus status = QProcess::CrashExit;
    int exitCode = -1;
    QByteArray diagnostics;
};

NativeSolverRun runNativeSolver(const QString &program,const QString &input,
                               const QString &directory)
{
    QProcess process;
    process.setWorkingDirectory(directory);
    process.setProcessChannelMode(QProcess::SeparateChannels);
    QElapsedTimer timer;
    timer.start();
    // Independent argv: no shell command interpolation or GUI/LLM dependency.
    process.start(program,QStringList{QStringLiteral("--file=")+input});
    NativeSolverRun run;
    run.started = process.waitForStarted(5000);
    if (run.started)
        run.finished = process.waitForFinished(std::max(1,30000-int(timer.elapsed())));
    if (!run.finished && process.state() != QProcess::NotRunning)
    {
        process.kill();
        process.waitForFinished(1000);
    }
    run.status = process.exitStatus();
    run.exitCode = process.exitCode();
    run.diagnostics = "input="+input.toUtf8()+"\nprogram="+program.toUtf8()+
        "\nprocess="+process.errorString().toUtf8()+"\nstdout:\n"+process.readAllStandardOutput()+
        "\nstderr:\n"+process.readAllStandardError();
    return run;
}

bool runSolver(RModel &model,QString &error)
{
    try
    {
        RSolver solver(model,QString(),QString());
        solver.run();
        return true;
    }
    catch (const RError &exception)
    {
        error = exception.getMessage();
    }
    catch (const std::exception &exception)
    {
        error = QString::fromUtf8(exception.what());
    }
    return false;
}

bool adapterRejects(const std::function<void()> &call,QString &message)
{
    try
    {
        call();
    }
    catch (const RError &error)
    {
        message = error.getMessage();
        return true;
    }
    catch (const std::exception &error)
    {
        message = QString::fromUtf8(error.what());
        return true;
    }
    return false;
}
}

class TestStructuralPiezoFE : public QObject
{
    Q_OBJECT

private slots:
    void freeVoltageActuationMatchesRodSolutionAndRightCharge();
    void clampedVoltagePublishesReactionAndCharge();
    void shortCircuitForceAndChargeMatchAnalyticRod();
    void polingReversalChangesTheExpectedSigns();
    void endpointReversalPreservesGlobalResponse();
    void rejectStaleMissingAndMixedSetups();
    void failedLinearSolveDoesNotPublishPiezoResults();
    void publishRejectsNonFiniteChargeBeforeMutation();
    void obliqueThreeDimensionalRodMatchesFreeAndBlockedAnalyticSolutions();
    void mixedSeriesMatchesDistinctMaterialForceChargeAndEnergy();
    void mixedClosedTriangleMatchesIndependentEquilibriumSolution();
    void rejectInvalidActivePassiveClassification();
    void publishRejectsElectricalAvailabilityTamperingBeforeMutation();
    void publishRejectsChangedPhysicalInputsAndForeignResults();
    void publishRejectsFiniteResultTamperingAndHandBuiltResults();
    void nodalThermalProvenanceTracksActualTemperatureAndResultFields();
    void nativeCliMixedTriangleRoundTrip_data();
    void nativeCliMixedTriangleRoundTrip();
    void nativeCliRejectsUnclassifiedMemberWithoutSuccessResult();
};

void TestStructuralPiezoFE::freeVoltageActuationMatchesRodSolutionAndRightCharge()
{
    constexpr double area = 0.5;
    constexpr double length = 2.0;
    constexpr double cE = 10.0;
    constexpr double e = 2.0;
    constexpr double epsS = 4.0;
    constexpr double voltage = 0.03;
    RModel model = makeRod(voltage,e,epsS);
    QString error;
    QVERIFY2(runSolver(model,error),error.toUtf8().constData());

    const double expectedExtension = -e*voltage/cE;
    const double expectedCharge = -(area/length)*(e*expectedExtension-epsS*voltage);
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_DISPLACEMENT,0,1)-expectedExtension) < 1.0e-10);
    QVERIFY(std::abs(elementValue(model,R_VARIABLE_ELECTRIC_CHARGE,0)-expectedCharge) < 1.0e-11);
    QVERIFY(std::abs(elementValue(model,R_VARIABLE_STRESS_NORMAL,0)) < 1.0e-8);
    QCOMPARE(model.getVariable(model.findVariable(R_VARIABLE_ELECTRIC_CHARGE)).getApplyType(),
             R_VARIABLE_APPLY_ELEMENT);
    QCOMPARE(model.getVariable(model.findVariable(R_VARIABLE_ELECTRIC_CHARGE)).getNValues(),
             model.getNElements());
    // The model has two point elements after its rod; charge is only assigned
    // to the explicit port element.
    QCOMPARE(elementValue(model,R_VARIABLE_ELECTRIC_CHARGE,1),0.0);
    QCOMPARE(elementValue(model,R_VARIABLE_ELECTRIC_CHARGE,2),0.0);
}

void TestStructuralPiezoFE::clampedVoltagePublishesReactionAndCharge()
{
    constexpr double area = 0.5;
    constexpr double length = 2.0;
    constexpr double e = 2.0;
    constexpr double epsS = 4.0;
    constexpr double voltage = 0.03;
    RModel model = makeRod(voltage,e,epsS,false,true,0.0);
    QString error;
    QVERIFY2(runSolver(model,error),error.toUtf8().constData());

    const double expectedForce = e*area/length*voltage;
    const double expectedRightCharge = epsS*area/length*voltage;
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_DISPLACEMENT,0,0)) < 1.0e-12);
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_DISPLACEMENT,0,1)) < 1.0e-12);
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_FORCE,0,0)+expectedForce) < 1.0e-10);
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_FORCE,0,1)-expectedForce) < 1.0e-10);
    QVERIFY(std::abs(elementValue(model,R_VARIABLE_ELECTRIC_CHARGE,0)-expectedRightCharge) < 1.0e-11);
    QVERIFY(std::abs(elementValue(model,R_VARIABLE_STRESS_NORMAL,0)-expectedForce/area) < 1.0e-10);
}

void TestStructuralPiezoFE::shortCircuitForceAndChargeMatchAnalyticRod()
{
    constexpr double area = 0.5;
    constexpr double length = 2.0;
    constexpr double cE = 10.0;
    constexpr double e = 2.0;
    constexpr double force = 0.025;
    RModel model = makeRod(0.0,e,4.0,false,false,0.0,force);
    QString error;
    QVERIFY2(runSolver(model,error),error.toUtf8().constData());

    const double expectedExtension = force*length/(cE*area);
    const double expectedRightCharge = -(area/length)*e*expectedExtension;
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_DISPLACEMENT,0,1)-expectedExtension) < 1.0e-10);
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_FORCE,0,0)+force) < 1.0e-10);
    QVERIFY(std::abs(nodalValue(model,R_VARIABLE_FORCE,0,1)-force) < 1.0e-10);
    QVERIFY(std::abs(elementValue(model,R_VARIABLE_ELECTRIC_CHARGE,0)-expectedRightCharge) < 1.0e-11);
    QVERIFY(std::abs(elementValue(model,R_VARIABLE_STRESS_NORMAL,0)-force/area) < 1.0e-10);
}

void TestStructuralPiezoFE::polingReversalChangesTheExpectedSigns()
{
    constexpr double imposedExtension = 0.01;
    const auto positiveFixed = RStructuralPiezoModelAdapter::solve(
        makeRod(0.0,2.0,4.0,false,true,imposedExtension));
    const auto negativeFixed = RStructuralPiezoModelAdapter::solve(
        makeRod(0.0,-2.0,4.0,false,true,imposedExtension));
    QVERIFY(std::abs(positiveFixed.memberRightCharges[0]+0.005) < 1.0e-12);
    QVERIFY(std::abs(negativeFixed.memberRightCharges[0]-0.005) < 1.0e-12);
    QVERIFY(std::abs(positiveFixed.memberAxialForces[0]-negativeFixed.memberAxialForces[0]) < 1.0e-12);

    constexpr double voltage = 0.03;
    const auto positiveFree = RStructuralPiezoModelAdapter::solve(
        makeRod(voltage,2.0,4.0));
    const auto negativeFree = RStructuralPiezoModelAdapter::solve(
        makeRod(voltage,-2.0,4.0));
    QVERIFY(std::abs(positiveFree.displacements[3]+0.006) < 1.0e-10);
    QVERIFY(std::abs(negativeFree.displacements[3]-0.006) < 1.0e-10);
    // At voltage-fixed mechanical free extension, the capacitive term and e*u
    // both retain sign, so the electrode charge does not simply reverse.
    QVERIFY(std::abs(positiveFree.memberRightCharges[0]-negativeFree.memberRightCharges[0]) < 1.0e-11);
}

void TestStructuralPiezoFE::endpointReversalPreservesGlobalResponse()
{
    constexpr double voltage = 0.03;
    const auto forward = RStructuralPiezoModelAdapter::solve(
        makeRod(voltage,2.0,4.0,false));
    const auto reversed = RStructuralPiezoModelAdapter::solve(
        makeRod(voltage,2.0,4.0,true));
    for (std::size_t dof=0;dof<forward.displacements.size();dof++)
        QVERIFY(std::abs(forward.displacements[dof]-reversed.displacements[dof]) < 1.0e-12);
    for (std::size_t dof=0;dof<forward.internalForces.size();dof++)
        QVERIFY(std::abs(forward.internalForces[dof]-reversed.internalForces[dof]) < 1.0e-12);
    QVERIFY(std::abs(forward.memberRightCharges[0]+reversed.memberRightCharges[0]) < 1.0e-11);
    QVERIFY(std::abs(forward.enthalpy-reversed.enthalpy) < 1.0e-12);
}

void TestStructuralPiezoFE::rejectStaleMissingAndMixedSetups()
{
    RModel stale = makeRod(0.03,2.0,4.0);
    stale.getNode(1).setX(2.1);
    QString message;
    QVERIFY(adapterRejects([&stale] { RStructuralPiezoModelAdapter::validate(stale); },message));
    QVERIFY(message.contains("different mesh",Qt::CaseInsensitive));

    RModel missing = makeRod(0.01,2.0,4.0,false,false,0.0,0.0,2,true);
    QVERIFY(adapterRejects([&missing] { RStructuralPiezoModelAdapter::validate(missing); },message));
    QVERIFY(message.contains("every TRUSS1 member",Qt::CaseInsensitive));

    RModel mixed = makeRod(0.01,2.0,4.0);
    mixed.getProblemSetup().getStructuralCableSetup().enabled = true;
    QVERIFY(adapterRejects([&mixed] { RStructuralPiezoModelAdapter::validate(mixed); },message));
    QVERIFY(message.contains("cannot be enabled together",Qt::CaseInsensitive));
}

void TestStructuralPiezoFE::failedLinearSolveDoesNotPublishPiezoResults()
{
    RModel model = makeRod(0.03,2.0,4.0);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setNOuterIterations(0);
    QString message;
    QVERIFY(adapterRejects([&model] { (void)RStructuralPiezoModelAdapter::solve(model); },message));
    QVERIFY(message.contains("residual certification",Qt::CaseInsensitive));
    QVERIFY(model.findVariable(R_VARIABLE_ELECTRIC_CHARGE) == RConstants::eod);
    QVERIFY(model.findVariable(R_VARIABLE_DISPLACEMENT) == RConstants::eod);

    QString runError;
    QVERIFY(!runSolver(model,runError));
    QVERIFY(model.findVariable(R_VARIABLE_ELECTRIC_CHARGE) == RConstants::eod);
    QVERIFY(model.findVariable(R_VARIABLE_DISPLACEMENT) == RConstants::eod);
}

void TestStructuralPiezoFE::publishRejectsNonFiniteChargeBeforeMutation()
{
    RModel model = makeRod(0.03,2.0,4.0);
    auto result = RStructuralPiezoModelAdapter::solve(model);
    result.memberRightCharges[0] = std::numeric_limits<double>::infinity();
    QString message;
    QVERIFY(adapterRejects([&model,&result] {
        RStructuralPiezoModelAdapter::publish(model,result);
    },message));
    QVERIFY(message.contains("non-finite",Qt::CaseInsensitive));
    QVERIFY(model.findVariable(R_VARIABLE_ELECTRIC_CHARGE) == RConstants::eod);
    QVERIFY(model.findVariable(R_VARIABLE_DISPLACEMENT) == RConstants::eod);
}

void TestStructuralPiezoFE::obliqueThreeDimensionalRodMatchesFreeAndBlockedAnalyticSolutions()
{
    // Length 2, direction (1,2,2)/3. Right y,z restraints remove transverse
    // zero-energy motion; the remaining x DOF still permits axial extension.
    for (bool blocked : {false,true})
    {
        for (bool reversed : {false,true})
        {
            RModel model = makeRod(0.03,2.0,4.0,reversed,blocked);
            model.getNode(1).set(2.0/3.0,4.0/3.0,4.0/3.0);
            model.getProblemSetup().getStructuralPiezoSetup().bindToModel(model);
            const auto result = RStructuralPiezoModelAdapter::solve(model);
            const double extension = blocked ? 0.0 : -2.0*0.03/10.0;
            const double force = blocked ? 2.0*0.5/2.0*0.03 : 0.0;
            const double charge = (reversed ? -1.0 : 1.0)*0.5/2.0*(4.0*0.03-2.0*extension);
            QVERIFY(std::abs(result.displacements[3]-3.0*extension) < 1.0e-10);
            for (uint axis=0;axis<3;axis++)
            {
                const double direction = axis == 0 ? 1.0/3.0 : 2.0/3.0;
                QVERIFY(std::abs(result.internalForces[axis]+direction*force) < 1.0e-10);
                QVERIFY(std::abs(result.internalForces[3+axis]-direction*force) < 1.0e-10);
            }
            QVERIFY(std::abs(result.memberAxialForces[0]-force) < 1.0e-10);
            QVERIFY(std::abs(result.memberRightCharges[0]-charge) < 1.0e-11);
            QVERIFY(result.memberElectricalPortAvailable[0]);
        }
    }
}

void TestStructuralPiezoFE::mixedSeriesMatchesDistinctMaterialForceChargeAndEnergy()
{
    RModel model = makeMixedNetwork(false);
    const auto result = RStructuralPiezoModelAdapter::solve(model);
    // N = applied end load through both rods. Active N = EA/L*d + eA/L*V;
    // passive N = EA/L*d, with independent, unequal EA values 5 and 16.
    const double force = 0.025;
    const double activeExtension = (force-2.0*0.5/2.0*0.03)/(10.0*0.5/2.0);
    const double passiveExtension = force*2.0/(20.0*0.8);
    const double charge = 0.5/2.0*(4.0*0.03-2.0*activeExtension);
    const double enthalpy = 0.5*0.5/2.0*(10.0*activeExtension*activeExtension+
        2.0*2.0*activeExtension*0.03-4.0*0.03*0.03)+0.5*force*passiveExtension;
    QVERIFY(std::abs(result.displacements[3]-activeExtension) < 1.0e-10);
    QVERIFY(std::abs(result.displacements[6]-activeExtension-passiveExtension) < 1.0e-10);
    for (double memberForce : result.memberAxialForces)
        QVERIFY(std::abs(memberForce-force) < 1.0e-10);
    QVERIFY(std::abs(result.internalForces[0]+force) < 1.0e-10);
    QVERIFY(std::abs(result.internalForces[3]) < 1.0e-10);
    QVERIFY(std::abs(result.internalForces[6]-force) < 1.0e-10);
    QVERIFY(std::abs(result.memberRightCharges[0]-charge) < 1.0e-11);
    QCOMPARE(result.memberRightCharges[1],0.0);
    QVERIFY(result.memberElectricalPortAvailable[0]);
    QVERIFY(!result.memberElectricalPortAvailable[1]);
    QVERIFY(std::abs(result.enthalpy-enthalpy) < 1.0e-12);
    RStructuralPiezoModelAdapter::publish(model,result);
    QVERIFY(std::abs(elementValue(model,R_VARIABLE_STRESS_NORMAL,1)-force/0.8) < 1.0e-10);
    QCOMPARE(elementValue(model,R_VARIABLE_ELECTRIC_CHARGE,1),0.0); // placeholder only
    QVERIFY(std::abs(elementValue(model,R_VARIABLE_ELECTRIC_CHARGE,0)-charge) < 1.0e-11);
}

void TestStructuralPiezoFE::mixedClosedTriangleMatchesIndependentEquilibriumSolution()
{
    RModel model = makeMixedNetwork(true);
    const auto result = RStructuralPiezoModelAdapter::solve(model);
    const double s = 1.0/std::sqrt(2.0);
    const double kBase = 10.0*0.5/2.0;
    const double kDiagonal = 20.0*0.8/std::sqrt(8.0);
    const double kVertical = 30.0*0.3/2.0;
    // Independently reduced stiffness in [u1x,u2x,u2y]:
    // [[kBase+h,-h,h],[-h,h,-h],[h,-h,h+kVertical]], h=kDiagonal/2.
    // All leading principal minors are positive, verifying these supports
    // remove planar rigid motions/mechanisms; z is restrained at every node.
    const double h = kDiagonal/2.0;
    QVERIFY(kBase+h > 0.0);
    QVERIFY(kBase*h > 0.0);
    QVERIFY(kBase*h*kVertical > 0.0);
    // Joint equilibrium, without reading adapter stiffness/recovery data:
    // -s*Ndiag=Fx, s*Ndiag+Nvertical=Fy, Nbase+s*Ndiag=0.
    const double fx = 0.04, fy = -0.02;
    const std::array<double,3> forces = {{fx,-fx/s,fy+fx}};
    const double u1x = (forces[0]-2.0*0.5/2.0*0.03)/kBase;
    const double u2y = forces[2]/kVertical;
    const double u2x = u1x+u2y-forces[1]/(s*kDiagonal);
    const std::array<double,9> displacements = {{0.0,0.0,0.0,u1x,0.0,0.0,u2x,u2y,0.0}};
    const std::array<double,9> internal = {{-fx,-forces[2],0.0,0.0,fx,0.0,fx,fy,0.0}};
    for (uint dof=0;dof<9;dof++)
    {
        QVERIFY(std::abs(result.displacements[dof]-displacements[dof]) < 1.0e-9);
        QVERIFY(std::abs(result.internalForces[dof]-internal[dof]) < 1.0e-9);
    }
    for (uint member=0;member<3;member++)
    {
        QVERIFY(std::abs(result.memberAxialForces[member]-forces[member]) < 1.0e-9);
        QVERIFY(result.memberElectricalPortAvailable[member] == (member == 0));
    }
    const double charge = 0.5/2.0*(4.0*0.03-2.0*u1x);
    const double energy = 0.5*0.5/2.0*(10.0*u1x*u1x+4.0*u1x*0.03-4.0*0.03*0.03)
        +0.5*forces[1]*forces[1]/kDiagonal+0.5*forces[2]*forces[2]/kVertical;
    QVERIFY(std::abs(result.memberRightCharges[0]-charge) < 1.0e-10);
    QCOMPARE(result.memberRightCharges[1],0.0);
    QCOMPARE(result.memberRightCharges[2],0.0);
    QVERIFY(std::abs(result.enthalpy-energy) < 1.0e-10);
    // Static reference only; a residual alone does not certify network stability.
}

void TestStructuralPiezoFE::rejectInvalidActivePassiveClassification()
{
    QString message;
    for (uint invalid=0;invalid<5;invalid++)
    {
        RModel model = makeMixedNetwork(false);
        auto &setup = model.getProblemSetup().getStructuralPiezoSetup();
        if (invalid == 0) setup.passiveElementIds.clear(); // unclassified rod
        if (invalid == 1) setup.passiveElementIds.push_back(1); // duplicate
        if (invalid == 2) setup.passiveElementIds[0] = 0; // overlaps active port
        if (invalid == 3) setup.passiveElementIds[0] = 2; // point, not TRUSS1
        if (invalid == 4) setup.electricalMode = RStructuralPiezoElectricalMode::NodalPotential;
        QVERIFY(adapterRejects([&model] { RStructuralPiezoModelAdapter::validate(model); },message));
        QVERIFY(model.findVariable(R_VARIABLE_DISPLACEMENT) == RConstants::eod);
    }
}

void TestStructuralPiezoFE::publishRejectsElectricalAvailabilityTamperingBeforeMutation()
{
    RModel model = makeMixedNetwork(false);
    const auto valid = RStructuralPiezoModelAdapter::solve(model);
    RStructuralPiezoModelAdapter::publish(model,valid);
    // Publish a real prior result, then attempt to corrupt both a mechanical
    // value and electrical identity. Rejection must preserve every old output.
    for (uint invalid=0;invalid<4;invalid++)
    {
        auto result = valid;
        result.displacements[3] = 999.0;
        if (invalid == 0) result.memberElectricalPortAvailable[0] = false;
        if (invalid == 1) result.memberElectricalPortAvailable[1] = true;
        if (invalid == 2) result.memberElectricalPortAvailable.pop_back();
        if (invalid == 3) result.memberRightCharges[1] = 1.0;
        QString message;
        QVERIFY(adapterRejects([&model,&result] {
            RStructuralPiezoModelAdapter::publish(model,result);
        },message));
        for (uint node=0;node<3;node++)
            for (uint component=0;component<3;component++)
            {
                QCOMPARE(nodalValue(model,R_VARIABLE_DISPLACEMENT,component,node),valid.displacements[3*node+component]);
                QCOMPARE(nodalValue(model,R_VARIABLE_FORCE,component,node),valid.internalForces[3*node+component]);
            }
        QVERIFY(std::abs(elementValue(model,R_VARIABLE_STRESS_NORMAL,0)-valid.memberAxialForces[0]/0.5) < 1.0e-12);
        QVERIFY(std::abs(elementValue(model,R_VARIABLE_STRESS_NORMAL,1)-valid.memberAxialForces[1]/0.8) < 1.0e-12);
        QCOMPARE(elementValue(model,R_VARIABLE_ELECTRIC_CHARGE,0),valid.memberRightCharges[0]);
        QCOMPARE(elementValue(model,R_VARIABLE_ELECTRIC_CHARGE,1),0.0);
    }
}

void TestStructuralPiezoFE::publishRejectsChangedPhysicalInputsAndForeignResults()
{
    RModel original = makeMixedNetwork(false);
    const auto result = RStructuralPiezoModelAdapter::solve(original);
    RStructuralPiezoModelAdapter::publish(original,result);
    const RModel baseline = original;
    RStructuralPiezoModelAdapter::publish(original,result);
    QVERIFY(sameResultVariables(original,baseline)); // output variables are not physical inputs
    for (uint changed=0;changed<8;changed++)
    {
        RModel model = original;
        auto &setup = model.getProblemSetup().getStructuralPiezoSetup();
        if (changed == 0) setup.ports[0].voltage = 0.04;
        if (changed == 1) assignRodMaterial(model.getLine(0),0.5,11.0);
        if (changed == 2)
            setComponent(model.getPoint(2).getBoundaryCondition(R_BOUNDARY_CONDITION_FORCE),
                         R_VARIABLE_FORCE_X,true,0.03);
        if (changed == 3)
        {
            model.getNode(2).setX(4.1);
            setup.bindToModel(model); // valid fresh binding, not merely stale mesh rejection
        }
        if (changed == 4)
            setComponent(model.getPoint(2).getBoundaryCondition(R_BOUNDARY_CONDITION_DISPLACEMENT),
                         R_VARIABLE_DISPLACEMENT_X,true,0.001);
        if (changed == 5) model.getLine(1).setCrossArea(0.9);
        if (changed == 6) setup.ports[0].coupling = 2.1;
        if (changed == 7) setup.ports[0].permittivity = 4.1;
        QString message;
        QVERIFY(adapterRejects([&model,&result] {
            RStructuralPiezoModelAdapter::publish(model,result);
        },message));
        QVERIFY2(message.contains("provenance",Qt::CaseInsensitive),message.toUtf8().constData());
        QVERIFY(sameResultVariables(model,baseline));
    }
    RModel foreign = makeMixedNetwork(false); // same dimensions and availability
    foreign.getProblemSetup().getStructuralPiezoSetup().ports[0].voltage = 0.05;
    const auto foreignResult = RStructuralPiezoModelAdapter::solve(foreign);
    RStructuralPiezoModelAdapter::publish(foreign,foreignResult);
    const RModel foreignBaseline = foreign;
    QString message;
    QVERIFY(adapterRejects([&foreign,&result] {
        RStructuralPiezoModelAdapter::publish(foreign,result);
    },message));
    QVERIFY(message.contains("provenance",Qt::CaseInsensitive));
    QVERIFY(sameResultVariables(foreign,foreignBaseline));
}

void TestStructuralPiezoFE::publishRejectsFiniteResultTamperingAndHandBuiltResults()
{
    RModel model = makeMixedNetwork(false);
    const auto valid = RStructuralPiezoModelAdapter::solve(model);
    RStructuralPiezoModelAdapter::publish(model,valid);
    const RModel baseline = model;
    for (uint changed=0;changed<9;changed++)
    {
        auto result = valid;
        if (changed == 0) result.displacements[3] += 0.01;
        if (changed == 1) result.internalForces[0] += 0.01;
        if (changed == 2) result.memberAxialForces[0] += 0.01;
        if (changed == 3) result.memberRightCharges[0] += 0.01;
        if (changed == 4) result.residualNorm += 0.01;
        if (changed == 5) result.enthalpy += 0.01;
        if (changed == 6) result.electricResidualNorm += 0.01;
        if (changed == 7) result.nodalPotentials.push_back(0.01);
        if (changed == 8) result.nodalConductorCharges.push_back(0.01);
        QString message;
        QVERIFY(adapterRejects([&model,&result] {
            RStructuralPiezoModelAdapter::publish(model,result);
        },message));
        QVERIFY2(message.contains("seal",Qt::CaseInsensitive),message.toUtf8().constData());
        QVERIFY(sameResultVariables(model,baseline));
    }
    RStructuralPiezoModelResult handmade;
    handmade.displacements = valid.displacements;
    handmade.internalForces = valid.internalForces;
    handmade.memberAxialForces = valid.memberAxialForces;
    handmade.memberRightCharges = valid.memberRightCharges;
    handmade.memberElectricalPortAvailable = valid.memberElectricalPortAvailable;
    handmade.residualNorm = valid.residualNorm;
    handmade.enthalpy = valid.enthalpy;
    handmade.electricalMode = valid.electricalMode;
    handmade.nodalPotentials = valid.nodalPotentials;
    handmade.nodalConductorCharges = valid.nodalConductorCharges;
    handmade.electricResidualNorm = valid.electricResidualNorm;
    QString message;
    QVERIFY(adapterRejects([&model,&handmade] {
        RStructuralPiezoModelAdapter::publish(model,handmade);
    },message));
    QVERIFY(message.contains("provenance",Qt::CaseInsensitive));
    QVERIFY(sameResultVariables(model,baseline));
}

void TestStructuralPiezoFE::nodalThermalProvenanceTracksActualTemperatureAndResultFields()
{
    // Direct adapter fixture consumes a pre-existing nodal heat field. No
    // thermal/dynamic solve is claimed by this static provenance regression.
    RModel model = makeRod(0.0,2.0,4.0);
    auto &setup = model.getProblemSetup().getStructuralPiezoSetup();
    setup.electricalMode = RStructuralPiezoElectricalMode::NodalPotential;
    setup.nodalBoundaryConditions.push_back({0,RStructuralPiezoBoundaryType::Potential,0.0});
    setup.nodalBoundaryConditions.push_back({1,RStructuralPiezoBoundaryType::Potential,0.03});
    RMaterialProperty alpha(RMaterialProperty::ThermalExpansionCoefficient);
    alpha.add(293.15,1.0e-3);
    model.getLine(0).getMaterial().add(alpha);
    RProblemTaskItem sequence;
    sequence.addChild(RProblemTaskItem(R_PROBLEM_HEAT));
    sequence.addChild(RProblemTaskItem(R_PROBLEM_STRESS));
    model.setProblemTaskTree(sequence);
    const uint temperaturePosition = model.addVariable(R_VARIABLE_TEMPERATURE);
    auto &temperature = model.getVariable(temperaturePosition);
    temperature.setApplyType(R_VARIABLE_APPLY_NODE);
    temperature.resize(1,2);
    temperature.setValue(0,0,393.15);
    temperature.setValue(0,1,393.15);
    setup.bindToModel(model);
    const auto valid = RStructuralPiezoModelAdapter::solve(model);
    RStructuralPiezoModelAdapter::publish(model,valid);
    const RModel baseline = model;
    RStructuralPiezoModelAdapter::publish(model,valid);
    QVERIFY(sameResultVariables(model,baseline));
    for (uint changed=0;changed<5;changed++)
    {
        RModel candidate = model;
        if (changed < 3)
        {
            auto &field = candidate.getVariable(candidate.findVariable(R_VARIABLE_TEMPERATURE));
            field.setValue(0,0,403.15);
            if (changed == 1) field.setValue(0,1,383.15); // identical mean, distinct field
            if (changed == 2)
            {
                field.setValue(0,0,393.15);
                auto &property = candidate.getLine(0).getMaterial().get(
                    candidate.getLine(0).getMaterial().findPosition(RMaterialProperty::ThermalExpansionCoefficient));
                property.clear();
                property.add(293.15,2.0e-3); // distinct integrated thermal strain
            }
        }
        else
        {
            if (changed == 3)
                candidate.getProblemSetup().getStructuralPiezoSetup().nodalBoundaryConditions[1].value = 0.04;
            else
            {
                auto &material = candidate.getLine(0).getMaterial();
                auto &modulus = material.get(material.findPosition(RMaterialProperty::ModulusOfElasiticity));
                modulus.clear();
                modulus.add(293.15,11.0);
            }
        }
        const RModel beforePublication = candidate;
        QString message;
        QVERIFY(adapterRejects([&candidate,&valid] {
            RStructuralPiezoModelAdapter::publish(candidate,valid);
        },message));
        QVERIFY(sameResultVariables(candidate,beforePublication));
    }
    for (uint changed=0;changed<2;changed++)
    {
        auto tampered = valid;
        if (changed == 0) tampered.nodalPotentials[1] += 0.01;
        else tampered.nodalConductorCharges[1] += 0.01;
        QString message;
        QVERIFY(adapterRejects([&model,&tampered] {
            RStructuralPiezoModelAdapter::publish(model,tampered);
        },message));
        QVERIFY(message.contains("seal",Qt::CaseInsensitive));
        QVERIFY(sameResultVariables(model,baseline));
    }
}

void TestStructuralPiezoFE::nativeCliMixedTriangleRoundTrip_data()
{
    QTest::addColumn<QString>("extension");
    QTest::newRow("ram") << QStringLiteral("ram");
    QTest::newRow("rbm") << QStringLiteral("rbm");
}

void TestStructuralPiezoFE::nativeCliMixedTriangleRoundTrip()
{
    QFETCH(QString,extension);
    const QString program = QString::fromUtf8(RANGE_PIEZO_FEA_SOLVER);
    QVERIFY2(QFileInfo(program).isFile(),qPrintable("Missing required native solver: "+program));
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    RModel source = makeMixedNetwork(true);
    const QString input = directory.filePath("mixed-triangle."+extension);
    QString resultPath;
    try
    {
        source.write(input,false); // actual standalone model, not a link wrapper
        resultPath = source.nativeWritePaths(input,true).first();
    }
    catch (const RError &error) { QFAIL(qPrintable(error.getMessage())); }
    QVERIFY(QFileInfo(input).isFile());
    QVERIFY2(resultPath != input,"Static CLI acceptance requires the actual numbered result target.");
    QVERIFY(!QFileInfo::exists(resultPath));
    const NativeSolverRun run = runNativeSolver(program,input,directory.path());
    QVERIFY2(run.started && run.finished,run.diagnostics.constData());
    QVERIFY2(run.status == QProcess::NormalExit && run.exitCode == 0,run.diagnostics.constData());
    QVERIFY2(QFileInfo(resultPath).isFile(),run.diagnostics.constData());
    RModel actual;
    try { actual.read(resultPath); } // load computed target, not the input/link
    catch (const RError &error) { QFAIL(qPrintable(error.getMessage())); }
    const auto &setup = actual.getProblemSetup().getStructuralPiezoSetup();
    QVERIFY(setup.enabled);
    QVERIFY(setup.electricalMode == RStructuralPiezoElectricalMode::LegacyTerminalVoltage);
    QCOMPARE(setup.ports.size(),std::size_t(1));
    QCOMPARE(setup.ports[0].elementId,0u);
    QCOMPARE(setup.ports[0].id,QString("active-base"));
    QCOMPARE(setup.ports[0].coupling,2.0);
    QCOMPARE(setup.ports[0].permittivity,4.0);
    QCOMPARE(setup.ports[0].voltage,0.03);
    QVERIFY(setup.passiveElementIds == std::vector<uint>({1u,2u}));
    QCOMPARE(setup.meshFingerprint,source.getProblemSetup().getStructuralPiezoSetup().meshFingerprint);
    // The persisted input identity supplies availability: only element 0 is an
    // electrical port. Zero charge slots on passive bars are not measurements.
    // Use the same independent joint-equilibrium derivation as the static test,
    // rather than using an adapter result to generate reference values.
    const double s = 1.0/std::sqrt(2.0);
    const double kBase = 10.0*0.5/2.0;
    const double kDiagonal = 20.0*0.8/std::sqrt(8.0);
    const double kVertical = 30.0*0.3/2.0;
    const double fx = 0.04, fy = -0.02;
    const std::array<double,3> forces = {{fx,-fx/s,fy+fx}};
    const std::array<double,3> areas = {{0.5,0.8,0.3}};
    const double u1x = (fx-2.0*0.5/2.0*0.03)/kBase;
    const double u2y = (fy+fx)/kVertical;
    const double u2x = u1x+u2y-forces[1]/(s*kDiagonal);
    const std::array<double,9> displacements = {{0.0,0.0,0.0,u1x,0.0,0.0,u2x,u2y,0.0}};
    const std::array<double,9> internal = {{-fx,-forces[2],0.0,0.0,fx,0.0,fx,fy,0.0}};
    for (uint node=0;node<3;node++)
        for (uint component=0;component<3;component++)
        {
            QVERIFY(std::abs(nodalValue(actual,R_VARIABLE_DISPLACEMENT,component,node)-
                             displacements[3*node+component]) < 1.0e-9);
            QVERIFY(std::abs(nodalValue(actual,R_VARIABLE_FORCE,component,node)-
                             internal[3*node+component]) < 1.0e-9);
        }
    for (uint member=0;member<3;member++)
    {
        const double stress = elementValue(actual,R_VARIABLE_STRESS_NORMAL,member);
        QVERIFY(std::abs(stress-forces[member]/areas[member]) < 1.0e-9);
        QVERIFY(std::abs(stress*areas[member]-forces[member]) < 1.0e-9);
    }
    const double charge = 0.5/2.0*(4.0*0.03-2.0*u1x);
    QVERIFY(std::abs(elementValue(actual,R_VARIABLE_ELECTRIC_CHARGE,0)-charge) < 1.0e-10);
    QCOMPARE(elementValue(actual,R_VARIABLE_ELECTRIC_CHARGE,1),0.0);
    QCOMPARE(elementValue(actual,R_VARIABLE_ELECTRIC_CHARGE,2),0.0);
    const auto &charges = actual.getVariable(actual.findVariable(R_VARIABLE_ELECTRIC_CHARGE));
    QCOMPARE(charges.getApplyType(),R_VARIABLE_APPLY_ELEMENT);
    QCOMPARE(charges.getNValues(),actual.getNElements());
}

void TestStructuralPiezoFE::nativeCliRejectsUnclassifiedMemberWithoutSuccessResult()
{
    const QString program = QString::fromUtf8(RANGE_PIEZO_FEA_SOLVER);
    QVERIFY2(QFileInfo(program).isFile(),qPrintable("Missing required native solver: "+program));
    QTemporaryDir directory;
    QVERIFY(directory.isValid());
    RModel source = makeMixedNetwork(true);
    source.getProblemSetup().getStructuralPiezoSetup().passiveElementIds.pop_back();
    const QString input = directory.filePath("unclassified-triangle.ram");
    const QString resultPath = source.nativeWritePaths(input,true).first();
    QVERIFY(resultPath != input);
    QString message;
    const bool writeRejected = adapterRejects([&source,&input] {
        source.write(input,false);
    },message);
    if (!writeRejected)
    {
        const NativeSolverRun run = runNativeSolver(program,input,directory.path());
        QVERIFY2(run.started && run.finished,run.diagnostics.constData());
        QVERIFY2(run.status == QProcess::NormalExit && run.exitCode != 0,run.diagnostics.constData());
        QVERIFY2(!run.diagnostics.contains("Model has been successfully saved"),run.diagnostics.constData());
    }
    else
    {
        QVERIFY2(message.contains("classif",Qt::CaseInsensitive) ||
                 message.contains("TRUSS1",Qt::CaseInsensitive),qPrintable(message));
    }
    QVERIFY(!QFileInfo::exists(resultPath));
    // Failure may leave an unchanged source, but never a successful computed
    // model. This is static CLI validation, not mixed dynamics or GUI coverage.
}

QTEST_MAIN(TestStructuralPiezoFE)
#include "tst_structural_piezo_fe.moc"

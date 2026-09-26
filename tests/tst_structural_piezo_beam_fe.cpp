#include <QtTest>
#include <QFile>
#include <QTemporaryDir>

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <vector>

#include <rbl_error.h>

#include <rml_boundary_condition.h>
#include <rml_environment_condition.h>
#include <rml_line.h>
#include <rml_matrix_solver_conf.h>
#include <rml_model.h>
#include <rml_point.h>
#include <rml_problem_task_item.h>
#include <rml_variable.h>

#include "rstructural_piezo_beam_fe_adapter.h"
#include "rstructural_typed_fe_control.h"

namespace
{
constexpr double beamLength = 1.0;
constexpr double youngModulus = 70.0e9;
constexpr double hostArea = 1.0e-4;
constexpr double hostIz = 1.0e-8;
constexpr double patchArea = 1.0e-5;
constexpr double patchThickness = 1.0e-3;
constexpr double patchOffsetY = 1.0e-2;
constexpr double patchYoungModulus = 60.0e9;
constexpr double piezoCoupling = -5.0;
constexpr double permittivity = 1.5e-8;
constexpr double drivenVoltage = 100.0;
constexpr double appliedTipForce = 3.0;

void setComponent(RBoundaryCondition &condition,RVariableType type,bool enabled,double value=0.0)
{
    const uint position = condition.findComponentPosition(type);
    if (position==RConstants::eod)
        throw std::logic_error("The point condition is missing a required component.");
    RConditionComponent &component = condition.getComponent(position);
    component.setEnabled(enabled);
    component.clear();
    if (enabled) component.add(0.0,value);
}

RBoundaryCondition translationCondition(const std::array<bool,3> &enabled)
{
    RBoundaryCondition condition(R_BOUNDARY_CONDITION_DISPLACEMENT);
    const RVariableType types[3] = {R_VARIABLE_DISPLACEMENT_X,
                                    R_VARIABLE_DISPLACEMENT_Y,
                                    R_VARIABLE_DISPLACEMENT_Z};
    for (uint axis=0;axis<3;axis++) setComponent(condition,types[axis],enabled[axis]);
    return condition;
}

RBoundaryCondition rotationCondition(const std::array<bool,3> &enabled)
{
    RBoundaryCondition condition(R_BOUNDARY_CONDITION_ROTATION);
    const RVariableType types[3] = {R_VARIABLE_ROTATION_X,
                                    R_VARIABLE_ROTATION_Y,
                                    R_VARIABLE_ROTATION_Z};
    for (uint axis=0;axis<3;axis++) setComponent(condition,types[axis],enabled[axis]);
    return condition;
}

RBoundaryCondition forceCondition(const std::array<double,3> &values)
{
    RBoundaryCondition condition(R_BOUNDARY_CONDITION_FORCE);
    const RVariableType types[3] = {R_VARIABLE_FORCE_X,R_VARIABLE_FORCE_Y,R_VARIABLE_FORCE_Z};
    for (uint axis=0;axis<3;axis++) setComponent(condition,types[axis],true,values[axis]);
    return condition;
}

RModel makeCantilever(uint segmentCount,double tipForce=0.0)
{
    RModel model;
    model.setNNodes(segmentCount+1);
    for (uint node=0;node<=segmentCount;node++)
        model.getNode(node).set(beamLength*double(node)/segmentCount,0.0,0.0);

    for (uint segment=0;segment<segmentCount;segment++)
    {
        RElement beam(R_ELEMENT_BEAM1);
        beam.setNodeId(0,segment);
        beam.setNodeId(1,segment+1);
        model.addElement(beam,true,0);
    }

    RLine &line = model.getLine(0);
    line.setCrossArea(hostArea);
    line.setSecondMomentY(hostIz);
    line.setSecondMomentZ(hostIz);
    line.setTorsionalConstant(1.0e-8);
    line.setLocalYDirection(RR3Vector(0.0,1.0,0.0));
    REnvironmentCondition &gravity =
        line.getEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
    const RVariableType gravityTypes[3] = {R_VARIABLE_G_ACCELERATION_X,
                                            R_VARIABLE_G_ACCELERATION_Y,
                                            R_VARIABLE_G_ACCELERATION_Z};
    for (RVariableType type : gravityTypes)
    {
        const uint position = gravity.findComponentPosition(type);
        RConditionComponent &component = gravity.getComponent(position);
        component.setEnabled(true);
        component.clear();
        component.add(0.0,0.0);
    }
    RMaterial material(RMaterial::Solid);
    RMaterialProperty modulus(RMaterialProperty::ModulusOfElasiticity);
    modulus.add(293.15,youngModulus);
    material.add(modulus);
    RMaterialProperty density(RMaterialProperty::Density);
    density.add(293.15,2700.0);
    material.add(density);
    line.setMaterial(material);

    for (uint node=0;node<=segmentCount;node++)
    {
        RElement point(R_ELEMENT_POINT);
        point.setNodeId(0,node);
        model.addElement(point,true,node);
        if (node==0)
        {
            model.getPoint(node).addBoundaryCondition(
                translationCondition({{true,true,true}}));
            model.getPoint(node).addBoundaryCondition(
                rotationCondition({{true,true,true}}));
        }
        else
        {
            model.getPoint(node).addBoundaryCondition(
                translationCondition({{false,false,true}}));
            model.getPoint(node).addBoundaryCondition(
                rotationCondition({{true,true,false}}));
        }
        if (node==segmentCount && tipForce!=0.0)
            model.getPoint(node).addBoundaryCondition(
                forceCondition({{0.0,tipForce,0.0}}));
    }

    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    RMatrixSolverConf &gmres = model.getMatrixSolverConf(RMatrixSolverConf::GMRES);
    gmres.setNInnerIterations(80);
    gmres.setNOuterIterations(20);
    gmres.setSolverCvgValue(1.0e-12);
    return model;
}

RStructuralPiezoBeamFEInput makeInput(const RModel &model,double voltage=drivenVoltage)
{
    RStructuralPiezoBeamFEInput input;
    const RLine &line = model.getLine(0);
    for (uint position=0;position<line.size();position++)
    {
        RStructuralPiezoBeamPatch patch;
        patch.elementId = line.get(position);
        patch.area = patchArea;
        patch.thickness = patchThickness;
        patch.centroidOffsetY = patchOffsetY;
        patch.youngModulus = patchYoungModulus;
        patch.coupling = piezoCoupling;
        patch.permittivity = permittivity;
        patch.prescribedValue = voltage;
        input.patches.push_back(patch);
    }
    input.translationalAbsoluteTolerance = 1.0e-7;
    input.translationalRelativeTolerance = 1.0e-10;
    input.rotationalAbsoluteTolerance = 1.0e-7;
    input.rotationalRelativeTolerance = 1.0e-10;
    input.electricAbsoluteTolerance = 1.0e-14;
    input.electricRelativeTolerance = 1.0e-10;
    return input;
}

bool rejects(const std::function<void()> &call)
{
    try
    {
        call();
    }
    catch (const RError &)
    {
        return true;
    }
    catch (const std::exception &)
    {
        return true;
    }
    return false;
}

double expectedCurvature(double voltage)
{
    const double totalAxialStiffness = youngModulus*hostArea+
        patchYoungModulus*patchArea;
    const double effectiveBendingStiffness = youngModulus*hostIz+
        patchYoungModulus*patchArea*patchOffsetY*patchOffsetY*
            youngModulus*hostArea/totalAxialStiffness;
    const double effectivePiezoMoment = piezoCoupling*patchArea*patchOffsetY*
        voltage/patchThickness*youngModulus*hostArea/totalAxialStiffness;
    return effectivePiezoMoment/effectiveBendingStiffness;
}

double expectedEnthalpy(const RStructuralPiezoBeamFEResult &result,double voltage)
{
    const std::size_t lastNode = result.displacements.size()/3-1;
    const double axialStrain = result.displacements[3*lastNode]/beamLength;
    const double curvature = result.rotations[3*lastNode+2]/beamLength;
    const double patchStrain = axialStrain-patchOffsetY*curvature;
    const double hostEnergy = 0.5*youngModulus*hostArea*axialStrain*axialStrain*beamLength+
        0.5*youngModulus*hostIz*curvature*curvature*beamLength;
    const double patchEnergy = 0.5*patchYoungModulus*patchArea*patchStrain*patchStrain*beamLength+
        piezoCoupling*patchArea/patchThickness*voltage*patchStrain*beamLength-
        0.5*permittivity*patchArea*beamLength/(patchThickness*patchThickness)*
            voltage*voltage;
    return hostEnergy+patchEnergy;
}
}

class TestStructuralPiezoBeamFE : public QObject
{
    Q_OBJECT

private slots:
    void beamMesh_data();
    void beamMesh();
    void pointForceUsesGlobalPointLoadAndBeamConventions();
    void candidatePublicationAndRollback();
    void rejectsUnsupportedGeometry();
    void rejectsUnsupportedPatchGravity();
    void typedVoltageAndChargeDriveActualBeamAndWork();
    void rejectsIgnoredTemperatureAndCompetingSetup();
    void typedDriverRejectsTimeVaryingNonControlLoad();
};

void TestStructuralPiezoBeamFE::beamMesh_data()
{
    QTest::addColumn<uint>("segments");
    QTest::newRow("2-elements") << 2u;
    QTest::newRow("4-elements") << 4u;
    QTest::newRow("8-elements") << 8u;
}

void TestStructuralPiezoBeamFE::beamMesh()
{
    QFETCH(uint,segments);
    RModel model = makeCantilever(segments);
    const RStructuralPiezoBeamFEInput voltageInput = makeInput(model);
    const uint variableCountBefore = model.getNVariables();
    const RStructuralPiezoBeamFEResult voltage =
        RStructuralPiezoBeamFEAdapter::solve(model,voltageInput);
    QCOMPARE(model.getNVariables(),variableCountBefore);

    const double curvature = expectedCurvature(drivenVoltage);
    const uint tip = segments;
    QVERIFY(std::fabs(voltage.rotations[3*tip+2]-curvature*beamLength)<
            2.0e-8*std::max(1.0,std::fabs(curvature*beamLength)));
    QVERIFY(std::fabs(voltage.displacements[3*tip+1]-0.5*curvature*beamLength*beamLength)<
            2.0e-8*std::max(1.0,std::fabs(curvature*beamLength*beamLength)));
    QVERIFY(voltage.translationalResidualNorm<1.0e-6);
    QVERIFY(voltage.rotationalResidualNorm<1.0e-6);
    QVERIFY(std::fabs(voltage.enthalpy-expectedEnthalpy(voltage,drivenVoltage))<
            1.0e-8*std::max(1.0,std::fabs(voltage.enthalpy)));

    const double step = 1.0;
    const RStructuralPiezoBeamFEResult plus = RStructuralPiezoBeamFEAdapter::solve(
        model,makeInput(model,drivenVoltage+step));
    const RStructuralPiezoBeamFEResult minus = RStructuralPiezoBeamFEAdapter::solve(
        model,makeInput(model,drivenVoltage-step));
    const double enthalpyDerivative = (plus.enthalpy-minus.enthalpy)/(2.0*step);
    const double totalCharge = std::accumulate(voltage.patchCharges.begin(),
                                               voltage.patchCharges.end(),0.0);
    QVERIFY(std::fabs(enthalpyDerivative+totalCharge)<
            2.0e-7*std::max(1.0e-12,std::fabs(totalCharge)));

    RStructuralPiezoBeamFEInput chargeInput = voltageInput;
    chargeInput.electricalMode = RStructuralPiezoBeamElectricalMode::PrescribedCharge;
    for (std::size_t patch=0;patch<chargeInput.patches.size();patch++)
        chargeInput.patches[patch].prescribedValue = voltage.patchCharges[patch];
    const RStructuralPiezoBeamFEResult charge =
        RStructuralPiezoBeamFEAdapter::solve(model,chargeInput);
    QCOMPARE(charge.patchVoltages.size(),voltage.patchVoltages.size());
    for (std::size_t node=0;node<voltage.displacements.size();node++)
        QVERIFY(std::fabs(charge.displacements[node]-voltage.displacements[node])<
                2.0e-8*std::max(1.0,std::fabs(voltage.displacements[node])));
    for (std::size_t patch=0;patch<voltage.patchVoltages.size();patch++)
        QVERIFY(std::fabs(charge.patchVoltages[patch]-voltage.patchVoltages[patch])<
                2.0e-8*std::max(1.0,std::fabs(voltage.patchVoltages[patch])));
    QVERIFY(charge.translationalResidualNorm<1.0e-6);
    QVERIFY(charge.rotationalResidualNorm<1.0e-6);
    QVERIFY(charge.electricResidualNorm<1.0e-12);
}

void TestStructuralPiezoBeamFE::pointForceUsesGlobalPointLoadAndBeamConventions()
{
    const double totalAxialStiffness = youngModulus*hostArea+
        patchYoungModulus*patchArea;
    const double effectiveBendingStiffness = youngModulus*hostIz+
        patchYoungModulus*patchArea*patchOffsetY*patchOffsetY*
            youngModulus*hostArea/totalAxialStiffness;
    const double actuatorMoment = piezoCoupling*patchArea*patchOffsetY*
        drivenVoltage/patchThickness*youngModulus*hostArea/totalAxialStiffness;
    const double expectedTipDeflection = appliedTipForce*beamLength*beamLength*beamLength/
        (3.0*effectiveBendingStiffness)+actuatorMoment*beamLength*beamLength/
        (2.0*effectiveBendingStiffness);
    const double expectedTipRotation = appliedTipForce*beamLength*beamLength/
        (2.0*effectiveBendingStiffness)+actuatorMoment*beamLength/effectiveBendingStiffness;
    const std::array<uint,3> segmentCounts = {{2u,4u,8u}};
    std::array<double,3> deflectionErrors{};
    for (std::size_t mesh=0;mesh<segmentCounts.size();mesh++)
    {
        const uint segments = segmentCounts[mesh];
        RModel model = makeCantilever(segments,appliedTipForce);
        const RStructuralPiezoBeamFEResult result =
            RStructuralPiezoBeamFEAdapter::solve(model,makeInput(model));
        const uint tip = segments;
        deflectionErrors[mesh] = std::fabs(
            result.displacements[3*tip+1]-expectedTipDeflection);
        QVERIFY(std::fabs(result.rotations[3*tip+2]-expectedTipRotation)<
                1.0e-8*std::max(1.0,std::fabs(expectedTipRotation)));
        QVERIFY(std::fabs(result.internalForces[6*tip+1]-appliedTipForce)<1.0e-6);
        QVERIFY(result.translationalResidualNorm<1.0e-6);
        QVERIFY(result.rotationalResidualNorm<1.0e-6);
    }
    // The offset patch needs a quadratic axial field as EB curvature varies.
    // BEAM1's linear axial interpolation converges monotonically to this
    // independent continuum tip-deflection reference under mesh refinement.
    QVERIFY(deflectionErrors[1]<0.26*deflectionErrors[0]);
    QVERIFY(deflectionErrors[2]<0.26*deflectionErrors[1]);
    QVERIFY(deflectionErrors[2]<5.0e-5*std::fabs(expectedTipDeflection));
}

void TestStructuralPiezoBeamFE::candidatePublicationAndRollback()
{
    RModel model = makeCantilever(2);
    RVariable sentinel(R_VARIABLE_DISPLACEMENT,R_VARIABLE_APPLY_NODE);
    sentinel.resize(3,model.getNNodes());
    sentinel.setValue(0,0,42.0);
    model.addVariable(sentinel);
    const RStructuralPiezoBeamFEInput input = makeInput(model);
    const RStructuralPiezoBeamFEResult result =
        RStructuralPiezoBeamFEAdapter::solve(model,input);
    QCOMPARE(model.getNVariables(),1u);

    const uint priorVariables = model.getNVariables();
    const auto rejectedWithoutMutation = [&](const RStructuralPiezoBeamFEResult &candidate)
    {
        QVERIFY(rejects([&]() { RStructuralPiezoBeamFEAdapter::publish(model,input,candidate); }));
        QCOMPARE(model.getNVariables(),priorVariables);
        QCOMPARE(model.getVariable(model.findVariable(R_VARIABLE_DISPLACEMENT)).getValue(0,0),42.0);
        QCOMPARE(model.findVariable(R_VARIABLE_ROTATION),uint(RConstants::eod));
        QCOMPARE(model.findVariable(R_VARIABLE_ELECTRIC_POTENTIAL),uint(RConstants::eod));
        QCOMPARE(model.findVariable(R_VARIABLE_ELECTRIC_CHARGE),uint(RConstants::eod));
    };
    const std::vector<std::function<void(RStructuralPiezoBeamFEResult &)>> finiteTampering = {
        [](RStructuralPiezoBeamFEResult &candidate) { candidate.displacements[0]+=1.0; },
        [](RStructuralPiezoBeamFEResult &candidate) { candidate.rotations[0]+=1.0; },
        [](RStructuralPiezoBeamFEResult &candidate) { candidate.internalForces[0]+=1.0; },
        [](RStructuralPiezoBeamFEResult &candidate) { candidate.patchVoltages[0]+=1.0; },
        [](RStructuralPiezoBeamFEResult &candidate) { candidate.patchCharges[0]+=1.0; },
        [](RStructuralPiezoBeamFEResult &candidate) { candidate.enthalpy+=1.0; },
        [](RStructuralPiezoBeamFEResult &candidate) { candidate.translationalResidualNorm+=1.0; },
        [](RStructuralPiezoBeamFEResult &candidate) { candidate.rotationalResidualNorm+=1.0; },
        [](RStructuralPiezoBeamFEResult &candidate) { candidate.translationalRoundoffBound+=1.0; },
        [](RStructuralPiezoBeamFEResult &candidate) { candidate.rotationalRoundoffBound+=1.0; },
        [](RStructuralPiezoBeamFEResult &candidate) { candidate.electricResidualNorm+=1.0; },
        [](RStructuralPiezoBeamFEResult &candidate) { candidate.electricRoundoffBound+=1.0; },
        [](RStructuralPiezoBeamFEResult &candidate)
        {
            candidate.modelFingerprint.append(QStringLiteral("-tampered"));
        }
    };
    for (const auto &tamper : finiteTampering)
    {
        RStructuralPiezoBeamFEResult candidate = result;
        tamper(candidate);
        rejectedWithoutMutation(candidate);
    }
    RStructuralPiezoBeamFEInput changedToleranceInput = input;
    changedToleranceInput.rotationalAbsoluteTolerance *= 2.0;
    QVERIFY(rejects([&]()
    {
        RStructuralPiezoBeamFEAdapter::publish(model,changedToleranceInput,result);
    }));
    QCOMPARE(model.getNVariables(),priorVariables);
    QCOMPARE(model.getVariable(model.findVariable(R_VARIABLE_DISPLACEMENT)).getValue(0,0),42.0);
    QCOMPARE(model.findVariable(R_VARIABLE_ROTATION),uint(RConstants::eod));

    // A later invalid batch entry cannot publish an earlier valid output.
    RVariable stagedRotation(R_VARIABLE_ROTATION,R_VARIABLE_APPLY_NODE);
    stagedRotation.resize(3,model.getNNodes());
    RVariable invalidType;
    QVERIFY(rejects([&]() {
        model.upsertVariablesAtomic({stagedRotation,invalidType});
    }));
    QCOMPARE(model.getNVariables(),priorVariables);
    QCOMPARE(model.findVariable(R_VARIABLE_ROTATION),uint(RConstants::eod));
    QCOMPARE(model.getVariable(model.findVariable(R_VARIABLE_DISPLACEMENT)).getValue(0,0),42.0);

    RStructuralPiezoBeamFEAdapter::publish(model,input,result);
    QCOMPARE(model.findVariable(R_VARIABLE_DISPLACEMENT),uint(0));
    const RVariable &displacement = model.getVariable(model.findVariable(R_VARIABLE_DISPLACEMENT));
    QCOMPARE(displacement.getApplyType(),R_VARIABLE_APPLY_NODE);
    QCOMPARE(displacement.getNVectors(),3u);
    QVERIFY(std::fabs(displacement.getValue(1,2)-result.displacements[7])<1.0e-14);
    const RVariable &potential = model.getVariable(model.findVariable(R_VARIABLE_ELECTRIC_POTENTIAL));
    QCOMPARE(potential.getApplyType(),R_VARIABLE_APPLY_ELEMENT);
    QCOMPARE(potential.getNValues(),model.getNElements());
    for (std::size_t patch=0;patch<input.patches.size();patch++)
        QVERIFY(std::fabs(potential.getValue(0,input.patches[patch].elementId)-
                          result.patchVoltages[patch])<1.0e-12);
}

void TestStructuralPiezoBeamFE::rejectsUnsupportedGeometry()
{
    RModel model = makeCantilever(2);
    model.getNode(1).setY(0.01);
    const RStructuralPiezoBeamFEInput input = makeInput(model);
    QVERIFY(rejects([&]() { RStructuralPiezoBeamFEAdapter::validate(model,input); }));
    RModel wrongOrder = makeCantilever(2);
    RStructuralPiezoBeamFEInput wrongInput = makeInput(wrongOrder);
    wrongInput.patches[0].thickness = 0.0;
    QVERIFY(rejects([&]() { RStructuralPiezoBeamFEAdapter::solve(wrongOrder,wrongInput); }));
    QCOMPARE(wrongOrder.findVariable(R_VARIABLE_DISPLACEMENT),uint(RConstants::eod));
}

void TestStructuralPiezoBeamFE::rejectsUnsupportedPatchGravity()
{
    RModel model = makeCantilever(2);
    REnvironmentCondition &gravity = model.getLine(0).getEnvironmentCondition(
        R_ENVIRONMENT_CONDITION_G_ACCELERATION);
    const uint position = gravity.findComponentPosition(R_VARIABLE_G_ACCELERATION_Y);
    QVERIFY(position!=RConstants::eod);
    RConditionComponent &component = gravity.getComponent(position);
    component.setEnabled(true);
    component.clear();
    component.add(0.0,1.0);
    const RStructuralPiezoBeamFEInput input = makeInput(model);
    QVERIFY(rejects([&]() { RStructuralPiezoBeamFEAdapter::validate(model,input); }));
}

void TestStructuralPiezoBeamFE::rejectsIgnoredTemperatureAndCompetingSetup()
{
    for (uint variant=0;variant<6;variant++)
    {
        RModel model=makeCantilever(2);
        const auto input=makeInput(model);
        if (variant==0)
        {
            auto &temperature=model.getLine(0).getEnvironmentCondition(R_ENVIRONMENT_CONDITION_TEMPERATURE);
            auto &component=temperature.getComponent(temperature.findComponentPosition(R_VARIABLE_TEMPERATURE));
            component.setEnabled(true);
            component.clear();
            component.add(0.0,350.0);
        }
        else if (variant==1)
        {
            RVariable temperature(R_VARIABLE_TEMPERATURE);
            temperature.setApplyType(R_VARIABLE_APPLY_NODE);
            temperature.resize(1,model.getNNodes());
            temperature.setValue(0,1,350.0);
            model.addVariable(temperature);
        }
        else if (variant==2) model.getProblemSetup().getStructuralCableSetup().enabled=true;
        else if (variant==3) model.getProblemSetup().getStructuralPiezoSetup().enabled=true;
        else if (variant==4) model.getProblemSetup().getStructuralInteractionSetup().getLinearConstraints().push_back(RLinearConstraint());
        else model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS_MODAL));
        QVERIFY(rejects([&]{RStructuralPiezoBeamFEAdapter::solve(model,input);}));
        QCOMPARE(model.findVariable(R_VARIABLE_DISPLACEMENT),RConstants::eod);
    }
}

void TestStructuralPiezoBeamFE::typedDriverRejectsTimeVaryingNonControlLoad()
{
    RModel model=makeCantilever(1,1.0);
    auto input=makeInput(model,0.0);
    RBoundaryCondition *load=nullptr;
    for (uint i=0;i<model.getPoint(1).getNBoundaryConditions();i++)
        if (model.getPoint(1).getBoundaryCondition(i).getType()==R_BOUNDARY_CONDITION_FORCE)
            load=&model.getPoint(1).getBoundaryCondition(i);
    QVERIFY(load!=nullptr);
    const uint position=load->findComponentPosition(R_VARIABLE_FORCE_Y);
    QVERIFY(position!=RConstants::eod);
    load->getComponent(position).add(1.0,2.0);
    RStructuralTypedControlChannel channel;
    channel.channelId="electrode";
    channel.quantity=RStructuralTypedControlQuantity::Voltage;
    RStructuralTypedControlSession session(1,0.1,{channel});
    QVERIFY(rejects([&]{RStructuralTypedFEControlDriver::piezoBeam(model,input,session,{{"electrode",0,0,0}});}));
}

void TestStructuralPiezoBeamFE::typedVoltageAndChargeDriveActualBeamAndWork()
{
    using Quantity=RStructuralTypedControlQuantity;
    const double curvature=expectedCurvature(drivenVoltage);
    const double strain=(patchYoungModulus*patchArea*patchOffsetY*curvature-
        piezoCoupling*patchArea*drivenVoltage/patchThickness)/
        (youngModulus*hostArea+patchYoungModulus*patchArea);
    const double patchStrain=strain-patchOffsetY*curvature;
    const double charge=(permittivity*patchArea*drivenVoltage/(patchThickness*patchThickness)-
                         piezoCoupling*patchArea*patchStrain/patchThickness)*beamLength/2.0;
    for (const bool chargeMode : {false,true})
    {
        RModel model=makeCantilever(2);
        auto input=makeInput(model,0.0);
        if (chargeMode) input.electricalMode=RStructuralPiezoBeamElectricalMode::PrescribedCharge;
        const Quantity quantity=chargeMode ? Quantity::Charge : Quantity::Voltage;
        std::vector<RStructuralTypedControlChannel> catalog;
        std::vector<RStructuralTypedFEControlBinding> bindings;
        for (const auto &patch : input.patches)
        {
            const QString id=QString("electrode-%1").arg(patch.elementId);
            RStructuralTypedControlChannel channel;
            channel.channelId=id;
            channel.quantity=quantity;
            catalog.push_back(channel);
            bindings.push_back({id,patch.elementId,0,0});
        }
        RStructuralTypedControlSession session(123,0.1,catalog);
        auto driver=RStructuralTypedFEControlDriver::piezoBeam(model,input,session,bindings);
        const auto before=session.serialize(driver->state());
        QVERIFY(rejects([&]
        {
            driver->advance([&](const RStructuralControlObservation &)
            {
                RStructuralTypedControlDecision decision;
                for (const auto &channel : catalog)
                    decision.requestedCommands.push_back({channel.channelId,Quantity::Force,drivenVoltage});
                decision.nextControllerState="rejected";
                return decision;
            });
        }));
        QCOMPARE(session.serialize(driver->state()),before);
        RStructuralControlObservation seen;
        driver->advance([&](const RStructuralControlObservation &observation)
        {
            seen=observation;
            RStructuralTypedControlDecision decision;
            for (const auto &channel : catalog)
                decision.requestedCommands.push_back({channel.channelId,quantity,chargeMode ? charge : drivenVoltage});
            decision.nextControllerState="accepted";
            return decision;
        });
        QCOMPARE(seen.acceptedTime,0.0);
        QVERIFY(seen.priorControllerState.isEmpty());
        QCOMPARE(driver->state().completedSteps,quint64(1));
        QCOMPARE(driver->state().controllerState,QByteArray("accepted"));
        const RModel &accepted=driver->model();
        const RVariable &u=accepted.getVariable(accepted.findVariable(R_VARIABLE_DISPLACEMENT));
        QVERIFY(std::fabs(u.getValue(1,2)-0.5*curvature*beamLength*beamLength)<1.0e-10);
        double work=0.0;
        for (const auto &sample : driver->lastChannels())
        {
            QCOMPARE(sample.quantity,quantity);
            QVERIFY(std::fabs(sample.cumulativeWorkJ-0.5*drivenVoltage*charge)<1.0e-12);
            QVERIFY(std::fabs(sample.averagePowerW-sample.cumulativeWorkJ/0.1)<1.0e-12);
            work+=sample.cumulativeWorkJ;
        }
        QVERIFY(std::fabs(work-drivenVoltage*charge)<1.0e-12);
        QCOMPARE(model.findVariable(R_VARIABLE_DISPLACEMENT),RConstants::eod);
        QTemporaryDir directory;
        QVERIFY(directory.isValid());
        const QString checkpoint=directory.filePath("piezo.sfec");
        driver->writeCheckpoint(checkpoint);
        auto resumed=RStructuralTypedFEControlDriver::readCheckpoint(checkpoint,"cpp-typed-v1");
        QCOMPARE(session.serialize(resumed->state()),session.serialize(driver->state()));
        QVERIFY(rejects([&]{RStructuralTypedFEControlDriver::readCheckpoint(checkpoint,"other-controller");}));
        const auto second=[&](const RStructuralControlObservation &observation)
        {
            if (observation.priorControllerState!=QByteArray("accepted"))
                throw std::runtime_error("Missing accepted opaque state after restart.");
            RStructuralTypedControlDecision decision;
            for (const auto &channel : catalog)
                decision.requestedCommands.push_back({channel.channelId,quantity,
                    0.5*(chargeMode ? charge : drivenVoltage)});
            decision.nextControllerState="second";
            return decision;
        };
        driver->advance(second);
        resumed->advance(second);
        QCOMPARE(resumed->state().completedSteps,driver->state().completedSteps);
        QCOMPARE(resumed->state().controllerState,driver->state().controllerState);
        for (uint i=0;i<catalog.size();i++)
            QVERIFY(std::fabs(resumed->state().cumulativeWorkJ[i]-driver->state().cumulativeWorkJ[i])<1.0e-12);
        const auto preserved=session.serialize(resumed->state());
        QFile damaged(checkpoint);
        QVERIFY(damaged.open(QIODevice::ReadWrite));
        QVERIFY(damaged.seek(50));
        const char corrupt='X';
        QCOMPARE(damaged.write(&corrupt,1),qint64(1));
        damaged.close();
        QVERIFY(rejects([&]{RStructuralTypedFEControlDriver::readCheckpoint(checkpoint,"cpp-typed-v1");}));
        QCOMPARE(session.serialize(resumed->state()),preserved);
    }
}

QTEST_APPLESS_MAIN(TestStructuralPiezoBeamFE)
#include "tst_structural_piezo_beam_fe.moc"

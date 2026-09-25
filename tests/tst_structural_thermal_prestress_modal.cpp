#include <QtTest>

#include <cmath>
#include <limits>

#include <rbl_error.h>
#include <rml_model.h>
#include <rml_problem_task_item.h>

#include "rsolver.h"

namespace
{
constexpr double referenceTemperature = 293.15;
constexpr double coldTemperature = 280.0;
constexpr double youngsModulus = 2.0e11;
constexpr double density = 7800.0;
constexpr double thermalExpansion = 1.2e-5;
constexpr double beamLength = 1.0;
constexpr double area = 1.0e-4;
constexpr double secondMoment = 8.333333333333333e-10;

void setValue(RCondition &condition,RVariableType type,double value)
{
    const uint position = condition.findComponentPosition(type);
    Q_ASSERT(position != RConstants::eod);
    RConditionComponent &component = condition.getComponent(position);
    component.setEnabled(true);
    component.clear();
    component.add(0.0,value);
}

void disableValue(RCondition &condition,RVariableType type)
{
    const uint position = condition.findComponentPosition(type);
    Q_ASSERT(position != RConstants::eod);
    condition.getComponent(position).setEnabled(false);
}

void setProperty(RMaterial &material,RMaterialProperty::Type type,double value)
{
    RMaterialProperty property(type);
    property.add(referenceTemperature,value);
    material.add(property);
}

void setZeroGravity(RModel &model)
{
    const RVariableType components[] = {
        R_VARIABLE_G_ACCELERATION_X,
        R_VARIABLE_G_ACCELERATION_Y,
        R_VARIABLE_G_ACCELERATION_Z
    };
    for (uint groupId=0;groupId<model.getNElementGroups();groupId++)
    {
        RElementGroup *group = model.getElementGroupPtr(groupId);
        if (!group)
            continue;
        REnvironmentCondition &gravity =
            group->getEnvironmentCondition(R_ENVIRONMENT_CONDITION_G_ACCELERATION);
        for (RVariableType componentType : components)
        {
            setValue(gravity,componentType,0.0);
        }
    }
}

void setModalConfiguration(RModel &model)
{
    model.getProblemSetup().setMassMatrixType(RMassMatrixType::Consistent);
    RModalSetup &modal = model.getProblemSetup().getModalSetup();
    modal.setMethod(R_MODAL_MULTIPLE_MODES);
    modal.setNModesToExtract(2);
    modal.setNIterations(500);
    modal.setConvergenceValue(1.0e-12);
    modal.setMode(0);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setOutputFrequency(0);
    model.getMatrixSolverConf(RMatrixSolverConf::CG).setSolverCvgValue(1.0e-12);
}

RModel makeBeam(uint segments,double temperature,bool includeStress)
{
    RModel model;
    model.setNNodes(segments+1);
    for (uint node=0;node<=segments;node++)
    {
        model.getNode(node).set(beamLength*double(node)/double(segments),0.0,0.0);
    }
    for (uint segment=0;segment<segments;segment++)
    {
        RElement beam(R_ELEMENT_BEAM1);
        beam.setNodeId(0,segment);
        beam.setNodeId(1,segment+1);
        model.addElement(beam,true,0);
    }

    RMaterial material(RMaterial::Solid);
    setProperty(material,RMaterialProperty::ModulusOfElasiticity,youngsModulus);
    setProperty(material,RMaterialProperty::PoissonRatio,0.30);
    setProperty(material,RMaterialProperty::Density,density);
    setProperty(material,RMaterialProperty::ThermalExpansionCoefficient,thermalExpansion);
    setProperty(material,RMaterialProperty::ThermalConductivity,16.0);
    setProperty(material,RMaterialProperty::HeatCapacity,500.0);
    setProperty(material,RMaterialProperty::Emissivity,0.0);

    RLine &line = model.getLine(0);
    line.setName(QStringLiteral("Thermal prestress beam"));
    line.setCrossArea(area);
    // Split the two bending frequencies so the first transverse mode is
    // unique; the analytical reference uses the weaker local-z plane.
    line.setSecondMomentY(1.25*secondMoment);
    line.setSecondMomentZ(secondMoment);
    line.setTorsionalConstant(2.0*secondMoment);
    line.setSectionModulusY(2.0*secondMoment/0.01);
    line.setSectionModulusZ(2.0*secondMoment/0.01);
    line.setTorsionalSectionModulus(2.0*secondMoment/0.01);
    line.setMaterial(material);

    for (uint node=0;node<=segments;node++)
    {
        RElement point(R_ELEMENT_POINT);
        point.setNodeId(0,node);
        model.addElement(point,true,node);
        model.getPoint(node).setVolume(0.0);

        RBoundaryCondition thermal(R_BOUNDARY_CONDITION_TEMPERATURE);
        setValue(thermal,R_VARIABLE_TEMPERATURE,temperature);
        model.getPoint(node).addBoundaryCondition(thermal);

        if (node == 0 || node == segments)
        {
            RBoundaryCondition translation(R_BOUNDARY_CONDITION_DISPLACEMENT);
            setValue(translation,R_VARIABLE_DISPLACEMENT_X,0.0);
            setValue(translation,R_VARIABLE_DISPLACEMENT_Y,0.0);
            setValue(translation,R_VARIABLE_DISPLACEMENT_Z,0.0);
            model.getPoint(node).addBoundaryCondition(translation);
        }

        // Suppress the unrelated torsional rigid mode. RY and RZ stay free,
        // giving simply supported transverse bending at both beam ends.
        RBoundaryCondition axialSpin(R_BOUNDARY_CONDITION_ROTATION);
        setValue(axialSpin,R_VARIABLE_ROTATION_X,0.0);
        disableValue(axialSpin,R_VARIABLE_ROTATION_Y);
        disableValue(axialSpin,R_VARIABLE_ROTATION_Z);
        model.getPoint(node).addBoundaryCondition(axialSpin);
    }

    setZeroGravity(model);
    RProblemTaskItem sequence;
    sequence.addChild(RProblemTaskItem(R_PROBLEM_HEAT));
    if (includeStress)
        sequence.addChild(RProblemTaskItem(R_PROBLEM_STRESS));
    sequence.addChild(RProblemTaskItem(R_PROBLEM_STRESS_MODAL));
    model.setProblemTaskTree(sequence);
    setModalConfiguration(model);
    return model;
}

double solveFrequency(RModel &model)
{
    RSolver solver(model,QString(),QString());
    solver.run();
    return model.getProblemSetup().getModalSetup().getFrequency();
}

double beamColumnFrequency(double axialForce)
{
    const double waveNumber = RConstants::pi/beamLength;
    const double omegaSquared = youngsModulus*secondMoment/(density*area)
                                    * std::pow(waveNumber,4)
                              + axialForce/(density*area)*waveNumber*waveNumber;
    return std::sqrt(omegaSquared)/(2.0*RConstants::pi);
}

RModel makeUnsupportedThermalTruss()
{
    RModel model;
    model.setNNodes(3);
    model.getNode(0).set(0.0,0.0,0.0);
    model.getNode(1).set(0.5*beamLength,0.0,0.0);
    model.getNode(2).set(beamLength,0.0,0.0);
    RElement truss(R_ELEMENT_TRUSS2);
    truss.setNodeId(0,0);
    truss.setNodeId(1,2);
    truss.setNodeId(2,1);
    model.addElement(truss,true,0);

    RMaterial material(RMaterial::Solid);
    setProperty(material,RMaterialProperty::ModulusOfElasiticity,youngsModulus);
    setProperty(material,RMaterialProperty::PoissonRatio,0.30);
    setProperty(material,RMaterialProperty::Density,density);
    setProperty(material,RMaterialProperty::ThermalExpansionCoefficient,thermalExpansion);
    setProperty(material,RMaterialProperty::ThermalConductivity,16.0);
    setProperty(material,RMaterialProperty::HeatCapacity,500.0);
    setProperty(material,RMaterialProperty::Emissivity,0.0);
    RLine &line = model.getLine(0);
    line.setCrossArea(area);
    line.setMaterial(material);

    for (uint node=0;node<3;node++)
    {
        RElement point(R_ELEMENT_POINT);
        point.setNodeId(0,node);
        model.addElement(point,true,node);
        model.getPoint(node).setVolume(0.0);
        RBoundaryCondition thermal(R_BOUNDARY_CONDITION_TEMPERATURE);
        setValue(thermal,R_VARIABLE_TEMPERATURE,coldTemperature);
        model.getPoint(node).addBoundaryCondition(thermal);
        if (node == 0 || node == 2)
        {
            RBoundaryCondition translation(R_BOUNDARY_CONDITION_DISPLACEMENT);
            setValue(translation,R_VARIABLE_DISPLACEMENT_X,0.0);
            setValue(translation,R_VARIABLE_DISPLACEMENT_Y,0.0);
            setValue(translation,R_VARIABLE_DISPLACEMENT_Z,0.0);
            model.getPoint(node).addBoundaryCondition(translation);
        }
        else
        {
            RBoundaryCondition lateral(R_BOUNDARY_CONDITION_DISPLACEMENT);
            setValue(lateral,R_VARIABLE_DISPLACEMENT_Y,0.0);
            setValue(lateral,R_VARIABLE_DISPLACEMENT_Z,0.0);
            disableValue(lateral,R_VARIABLE_DISPLACEMENT_X);
            model.getPoint(node).addBoundaryCondition(lateral);
        }
    }

    setZeroGravity(model);
    RProblemTaskItem sequence;
    sequence.addChild(RProblemTaskItem(R_PROBLEM_HEAT));
    sequence.addChild(RProblemTaskItem(R_PROBLEM_STRESS));
    sequence.addChild(RProblemTaskItem(R_PROBLEM_STRESS_MODAL));
    model.setProblemTaskTree(sequence);
    setModalConfiguration(model);
    return model;
}
}

class TestStructuralThermalPrestressModal : public QObject
{
    Q_OBJECT

private slots:
    void uniformCoolingPrestressConvergesToBeamColumnFormula();
    void absentOrZeroAcceptedForceLeavesBaselineModalStiffness();
    void unsupportedAcceptedStressFailsClosed();
};

void TestStructuralThermalPrestressModal::uniformCoolingPrestressConvergesToBeamColumnFormula()
{
    const double expectedForce = youngsModulus*area*thermalExpansion
                               *(referenceTemperature-coldTemperature);
    const double expectedFrequency = beamColumnFrequency(expectedForce);
    double previousError = std::numeric_limits<double>::infinity();
    for (uint segments : {2u,4u,8u})
    {
        RModel model = makeBeam(segments,coldTemperature,true);
        const double frequency = solveFrequency(model);
        const double error = std::fabs(frequency-expectedFrequency);
        QVERIFY2(std::isfinite(frequency) && frequency > 0.0,
                 qPrintable(QString("%1-element beam returned frequency %2")
                            .arg(segments).arg(frequency)));
        QVERIFY2(error < previousError,
                 qPrintable(QString("Beam-column errors did not decrease: previous=%1, current=%2")
                            .arg(previousError).arg(error)));
        previousError = error;
        if (segments == 8)
        {
            QVERIFY2(error < 0.005*expectedFrequency,
                     qPrintable(QString("Fine-mesh frequency %1 differs from beam-column value %2 by %3%%")
                                .arg(frequency).arg(expectedFrequency)
                                .arg(100.0*error/expectedFrequency)));
        }
    }
}

void TestStructuralThermalPrestressModal::absentOrZeroAcceptedForceLeavesBaselineModalStiffness()
{
    RModel noStress = makeBeam(8,coldTemperature,false);
    const double noStressFrequency = solveFrequency(noStress);
    const double zeroForceFrequency = beamColumnFrequency(0.0);
    QVERIFY2(std::fabs(noStressFrequency-zeroForceFrequency) < 0.005*zeroForceFrequency,
             qPrintable(QString("Heat->Modal without a Stress field changed the baseline frequency: %1 vs %2")
                        .arg(noStressFrequency).arg(zeroForceFrequency)));

    RModel zeroStress = makeBeam(8,referenceTemperature,true);
    const double zeroStressFrequency = solveFrequency(zeroStress);
    QVERIFY2(std::fabs(zeroStressFrequency-zeroForceFrequency) < 0.005*zeroForceFrequency,
             qPrintable(QString("Zero thermal force changed the baseline frequency: %1 vs %2")
                        .arg(zeroStressFrequency).arg(zeroForceFrequency)));
}

void TestStructuralThermalPrestressModal::unsupportedAcceptedStressFailsClosed()
{
    RModel stressOnly = makeUnsupportedThermalTruss();
    RProblemTaskItem heatThenStress;
    heatThenStress.addChild(RProblemTaskItem(R_PROBLEM_HEAT));
    heatThenStress.addChild(RProblemTaskItem(R_PROBLEM_STRESS));
    stressOnly.setProblemTaskTree(heatThenStress);
    RSolver stressSolver(stressOnly,QString(),QString());
    stressSolver.run();

    RModel model = makeUnsupportedThermalTruss();
    RSolver solver(model,QString(),QString());
    bool rejectedUnsupportedElement = false;
    try
    {
        solver.run();
    }
    catch (const RError &error)
    {
        rejectedUnsupportedElement = error.getMessage().contains(
            "Thermal-prestress modal geometric stiffness supports only TRUSS1 and BEAM1");
    }
    QVERIFY2(rejectedUnsupportedElement,
             "TRUSS2 was not rejected during the modal preparation stage.");
}

QTEST_APPLESS_MAIN(TestStructuralThermalPrestressModal)
#include "tst_structural_thermal_prestress_modal.moc"

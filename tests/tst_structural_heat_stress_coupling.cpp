#include <QtTest>

#include <cmath>
#include <limits>

#include <rbl_error.h>
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

void property(RMaterial &material,RMaterialProperty::Type type,double value)
{
    RMaterialProperty entry(type);
    entry.add(referenceKelvin,value);
    material.add(entry);
}

RModel heatedRod(bool fixRight,bool temperatureDependentExpansion = false)
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
    property(material,RMaterialProperty::ModulusOfElasiticity,modulus);
    property(material,RMaterialProperty::PoissonRatio,0.0);
    property(material,RMaterialProperty::Density,8000.0);
    if (temperatureDependentExpansion)
    {
        RMaterialProperty variableAlpha(RMaterialProperty::ThermalExpansionCoefficient);
        variableAlpha.add(referenceKelvin,1.0e-5);
        variableAlpha.add(referenceKelvin+100.0,3.0e-5);
        material.add(variableAlpha);
    }
    else
    {
        property(material,RMaterialProperty::ThermalExpansionCoefficient,expansion);
    }
    property(material,RMaterialProperty::ThermalConductivity,2.0);
    property(material,RMaterialProperty::HeatCapacity,500.0);
    property(material,RMaterialProperty::Emissivity,0.0);
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
}

class TestStructuralHeatStressCoupling : public QObject
{
    Q_OBJECT
private slots:
    void heatThenFreeExpansion();
    void heatThenConstrainedCompression();
    void temperatureDependentExpansionUsesIntegratedAlpha();
    void persistedTemperatureChangesModalFrequency();
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

QTEST_APPLESS_MAIN(TestStructuralHeatStressCoupling)
#include "tst_structural_heat_stress_coupling.moc"

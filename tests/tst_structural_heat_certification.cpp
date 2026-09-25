#include <QtTest>

#include <cmath>
#include <limits>
#include <stdexcept>

#include <rbl_error.h>

#include <rml_boundary_condition.h>
#include <rml_line.h>
#include <rml_model.h>
#include <rml_point.h>
#include <rml_problem_task_item.h>
#include <rml_variable.h>

#include "rsolverheat.h"
#include "rsolvershareddata.h"

namespace
{

void setBoundaryValue(RBoundaryCondition &condition,RVariableType variableType,double value)
{
    const uint position = condition.findComponentPosition(variableType);
    if (position == RConstants::eod)
    {
        throw std::logic_error("Heat test boundary condition is missing a required variable.");
    }
    RConditionComponent &component = condition.getComponent(position);
    component.setEnabled(true);
    component.clear();
    component.add(0.0,value);
}

RModel makeHeatRod(uint nNodes,bool addInteriorHeat)
{
    RModel model;
    model.setNNodes(nNodes);
    for (uint node=0;node<nNodes;node++)
    {
        model.getNode(node).set(double(node),0.0,0.0);
    }

    for (uint segment=0;segment+1<nNodes;segment++)
    {
        RElement element(R_ELEMENT_TRUSS1);
        element.setNodeId(0,segment);
        element.setNodeId(1,segment+1);
        model.addElement(element,true,0);
    }

    RMaterial material(RMaterial::Solid);
    material.setName("Heat certification rod");
    RMaterialProperty conductivity(RMaterialProperty::ThermalConductivity);
    conductivity.add(293.15,2.0);
    material.add(conductivity);
    RMaterialProperty density(RMaterialProperty::Density);
    density.add(293.15,1.0);
    material.add(density);
    RMaterialProperty heatCapacity(RMaterialProperty::HeatCapacity);
    heatCapacity.add(293.15,1.0);
    material.add(heatCapacity);
    RMaterialProperty emissivity(RMaterialProperty::Emissivity);
    emissivity.add(293.15,0.0);
    material.add(emissivity);
    RLine &line = model.getLine(0);
    line.setCrossArea(1.0);
    line.setMaterial(material);

    for (uint node=0;node<nNodes;node++)
    {
        RElement pointElement(R_ELEMENT_POINT);
        pointElement.setNodeId(0,node);
        model.addElement(pointElement,true,node);
        RPoint &point = model.getPoint(node);
        point.setVolume(0.0);

        if (node == 0 || node+1 == nNodes)
        {
            RBoundaryCondition temperature(R_BOUNDARY_CONDITION_TEMPERATURE);
            setBoundaryValue(temperature,R_VARIABLE_TEMPERATURE,node == 0 ? 300.0 : 400.0);
            point.addBoundaryCondition(temperature);
        }
        if (addInteriorHeat && node == nNodes/2)
        {
            RBoundaryCondition heat(R_BOUNDARY_CONDITION_HEAT);
            setBoundaryValue(heat,R_VARIABLE_HEAT,100.0);
            point.addBoundaryCondition(heat);
        }
    }

    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_HEAT));
    model.getTimeSolver().setEnabled(false);
    RMatrixSolverConf &cg = model.getMatrixSolverConf(RMatrixSolverConf::CG);
    cg.setNOuterIterations(1000);
    cg.setSolverCvgValue(1.0e-12);
    cg.setOutputFrequency(0);
    return model;
}

double nodeTemperature(const RModel &model,uint node)
{
    const uint position = model.findVariable(R_VARIABLE_TEMPERATURE);
    if (position == RConstants::eod)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return model.getVariable(position).getValue(0,node);
}

}

class TestStructuralHeatCertification : public QObject
{
    Q_OBJECT

private slots:

    void steadyRodProducesCertifiedAnalyticTemperatures();
    void allFixedRodCertifiesZeroSizedSystem();
    void underiteratedCgDoesNotPublishTemperature();
    void invalidToleranceResetsPriorConvergence();
};

void TestStructuralHeatCertification::steadyRodProducesCertifiedAnalyticTemperatures()
{
    RModel model = makeHeatRod(5,false);
    RSolverSharedData sharedData;
    RSolverHeat heat(&model,QString(),QString(),sharedData);

    QVERIFY(!heat.hasConverged());
    heat.run(true,0);
    QVERIFY(heat.hasConverged());

    for (uint node=0;node<model.getNNodes();node++)
    {
        const double expected = 300.0 + 25.0*double(node);
        QVERIFY2(std::fabs(nodeTemperature(model,node)-expected) < 1.0e-7,
                 qPrintable(QString("Node %1 temperature %2, expected %3")
                            .arg(node).arg(nodeTemperature(model,node),0,'g',16).arg(expected)));
    }
}

void TestStructuralHeatCertification::allFixedRodCertifiesZeroSizedSystem()
{
    RModel model = makeHeatRod(2,false);
    RSolverSharedData sharedData;
    RSolverHeat heat(&model,QString(),QString(),sharedData);

    heat.run(true,0);
    QVERIFY(heat.hasConverged());
    QCOMPARE(nodeTemperature(model,0),300.0);
    QCOMPARE(nodeTemperature(model,1),400.0);
}

void TestStructuralHeatCertification::underiteratedCgDoesNotPublishTemperature()
{
    RModel model = makeHeatRod(7,true);
    const double sentinel = 321.25;
    const uint temperaturePosition = model.addVariable(R_VARIABLE_TEMPERATURE);
    RVariable &temperature = model.getVariable(temperaturePosition);
    temperature.setApplyType(R_VARIABLE_APPLY_NODE);
    temperature.resize(1,model.getNNodes());
    for (uint node=0;node<model.getNNodes();node++)
    {
        temperature.setValue(0,node,sentinel);
    }

    model.getMatrixSolverConf(RMatrixSolverConf::CG).setNOuterIterations(1);
    RSolverSharedData sharedData;
    RSolverHeat heat(&model,QString(),QString(),sharedData);

    QString errorMessage;
    try
    {
        heat.run(true,0);
    }
    catch (const RError &error)
    {
        errorMessage = error.getMessage();
    }

    QVERIFY2(!errorMessage.isEmpty(),"Underiterated heat CG result was accepted.");
    QVERIFY2(errorMessage.contains("did not converge",Qt::CaseInsensitive),qPrintable(errorMessage));
    QVERIFY(!heat.hasConverged());
    QCOMPARE(model.findVariable(R_VARIABLE_TEMPERATURE),temperaturePosition);
    for (uint node=0;node<model.getNNodes();node++)
    {
        QVERIFY2(std::fabs(nodeTemperature(model,node)-sentinel) < 1.0e-12,
                 qPrintable(QString("Node %1 published %2 instead of preserving %3")
                            .arg(node).arg(nodeTemperature(model,node),0,'g',16).arg(sentinel)));
    }
}

void TestStructuralHeatCertification::invalidToleranceResetsPriorConvergence()
{
    RModel model = makeHeatRod(3,false);
    RSolverSharedData sharedData;
    RSolverHeat heat(&model,QString(),QString(),sharedData);
    heat.run(true,0);
    QVERIFY(heat.hasConverged());

    model.getMatrixSolverConf(RMatrixSolverConf::CG).setSolverCvgValue(
        std::numeric_limits<double>::quiet_NaN());
    QString errorMessage;
    try
    {
        heat.run(false,1);
    }
    catch (const RError &error)
    {
        errorMessage = error.getMessage();
    }
    QVERIFY2(errorMessage.contains("tolerance",Qt::CaseInsensitive),qPrintable(errorMessage));
    QVERIFY(!heat.hasConverged());
}

QTEST_APPLESS_MAIN(TestStructuralHeatCertification)
#include "tst_structural_heat_certification.moc"

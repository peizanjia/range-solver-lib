#include <QtTest>

#include <cmath>
#include <stdexcept>

#include <rbl_error.h>

#include <rml_model.h>
#include <rml_point.h>
#include <rml_problem_task_item.h>
#include <rml_surface.h>
#include <rml_variable.h>
#include <rml_volume.h>

#include "rsolverheat.h"
#include "rsolvershareddata.h"

namespace
{

enum class CellType
{
    Quad,
    Quad2,
    Tetra
};

struct PulseResponse
{
    double nodalIncrementSum = 0.0;
    double referenceWeightedEnergy = 0.0;
};

void addConstantProperty(RMaterial &material,RMaterialProperty::Type type,double value)
{
    RMaterialProperty property(type);
    property.add(RVariable::getInitValue(R_VARIABLE_TEMPERATURE),value);
    material.add(property);
}

RMaterial thermalMaterial()
{
    RMaterial material(RMaterial::Solid);
    material.setName("Reference mass pulse test material");
    addConstantProperty(material,RMaterialProperty::ThermalConductivity,0.0);
    addConstantProperty(material,RMaterialProperty::Density,1.0);
    addConstantProperty(material,RMaterialProperty::HeatCapacity,1.0);
    addConstantProperty(material,RMaterialProperty::Emissivity,0.0);
    return material;
}

RModel makePulseModel(CellType type,bool deformed,bool degenerateReference = false)
{
    RModel model;
    const uint nNodes = type == CellType::Quad2 ? 9U : 4U;
    const bool isSurface = type != CellType::Tetra;
    model.setNNodes(nNodes);
    if (isSurface)
    {
        if (degenerateReference)
        {
            model.getNode(0).set(0.0,0.0,0.0);
            model.getNode(1).set(1.0,0.0,0.0);
            model.getNode(2).set(1.0,0.0,0.0);
            model.getNode(3).set(0.0,0.0,0.0);
        }
        else
        {
            model.getNode(0).set(0.0,0.0,0.0);
            model.getNode(1).set(1.0,0.0,0.0);
            model.getNode(2).set(1.0,1.0,0.0);
            model.getNode(3).set(0.0,1.0,0.0);
        }
        if (type == CellType::Quad2)
        {
            model.getNode(4).set(0.5,0.0,0.0);
            model.getNode(5).set(1.0,0.5,0.0);
            model.getNode(6).set(0.5,1.0,0.0);
            model.getNode(7).set(0.0,0.5,0.0);
            model.getNode(8).set(0.5,0.5,0.0);
        }
        RElement quad(type == CellType::Quad2 ? R_ELEMENT_QUAD2 : R_ELEMENT_QUAD1);
        for (uint node=0;node<quad.size();node++)
        {
            quad.setNodeId(node,node);
        }
        model.addElement(quad,true,0);
        RSurface &surface = model.getSurface(0);
        surface.setThickness(1.0);
        surface.setMaterial(thermalMaterial());
    }
    else
    {
        model.getNode(0).set(0.0,0.0,0.0);
        model.getNode(1).set(1.0,0.0,0.0);
        model.getNode(2).set(0.0,1.0,0.0);
        model.getNode(3).set(0.0,0.0,1.0);
        RElement tetra(R_ELEMENT_TETRA1);
        for (uint node=0;node<4;node++)
        {
            tetra.setNodeId(node,node);
        }
        model.addElement(tetra,true,0);
        model.getVolume(0).setMaterial(thermalMaterial());
    }

    for (uint node=0;node<nNodes;node++)
    {
        RElement point(R_ELEMENT_POINT);
        point.setNodeId(0,node);
        model.addElement(point,true,node);
        model.getPoint(node).setVolume(0.0);
    }

    const uint temperaturePosition = model.addVariable(R_VARIABLE_TEMPERATURE);
    RVariable &temperature = model.getVariable(temperaturePosition);
    temperature.setApplyType(R_VARIABLE_APPLY_NODE);
    temperature.resize(1,model.getNNodes());
    for (uint node=0;node<model.getNNodes();node++)
    {
        temperature.setValue(0,node,300.0);
    }

    if (deformed)
    {
        const uint displacementPosition = model.addVariable(R_VARIABLE_DISPLACEMENT);
        RVariable &displacement = model.getVariable(displacementPosition);
        displacement.setApplyType(R_VARIABLE_APPLY_NODE);
        displacement.resize(3,model.getNNodes());
        for (uint node=0;node<model.getNNodes();node++)
        {
            const RNode &position = model.getNode(node);
            double ux = position.getX();
            double uy = position.getY();
            double uz = position.getZ();
            if (degenerateReference)
            {
                // Move the initially collapsed QUAD1 into a valid current
                // square while preserving its degenerate reference geometry.
                if (node >= 2)
                {
                    uy = 1.0;
                }
            }
            displacement.setValue(0,node,ux);
            displacement.setValue(1,node,uy);
            displacement.setValue(2,node,uz);
        }
    }

    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_HEAT));
    RTimeSolver &time = model.getTimeSolver();
    time.setEnabled(true);
    time.setTimes({1.0});
    time.setInputTimeStepSize(1.0);
    time.setCurrentTimeStep(0);
    time.setComputedTime(0.0);
    time.setTimeMarchApproximation(R_TIME_MARCH_BACKWARD);

    RMatrixSolverConf &cg = model.getMatrixSolverConf(RMatrixSolverConf::CG);
    cg.setNOuterIterations(1000);
    cg.setSolverCvgValue(1.0e-12);
    cg.setOutputFrequency(0);

    return model;
}

PulseResponse solvePulse(CellType type,bool referenceMass)
{
    RModel model = makePulseModel(type,true);
    const double initialTemperature = 300.0;
    RSolverSharedData sharedData;
    RSolverHeat heat(&model,QString(),QString(),sharedData);
    heat.setControlledPointHeat({{model.getPoint(0).get(0),1.0}});
    heat.setReferenceSurfaceVolumeThermalMass(referenceMass);
    heat.run(false,0);
    if (!heat.hasConverged())
    {
        throw std::runtime_error("Reference thermal mass pulse did not converge.");
    }

    const uint temperaturePosition = model.findVariable(R_VARIABLE_TEMPERATURE);
    if (temperaturePosition == RConstants::eod)
    {
        throw std::runtime_error("Reference thermal mass pulse did not publish nodal temperature.");
    }
    const RVariable &temperature = model.getVariable(temperaturePosition);
    PulseResponse response;
    for (uint node=0;node<model.getNNodes();node++)
    {
        response.nodalIncrementSum += temperature.getValue(0,node)-initialTemperature;
    }

    const RElement &element = model.getElement(0);
    const double referenceDetJ = type == CellType::Tetra ? 1.0 : 0.25;
    for (uint point=0;point<RElement::getNIntegrationPoints(element.getType());point++)
    {
        const RElementShapeFunction &shape = RElement::getShapeFunction(element.getType(),point);
        const RRVector &N = shape.getN();
        for (uint node=0;node<element.size();node++)
        {
            response.referenceWeightedEnergy += N[node]*referenceDetJ*shape.getW()
                *(temperature.getValue(0,element.getNodeId(node))-initialTemperature);
        }
    }
    return response;
}

}

class TestStructuralReferenceThermalMass : public QObject
{
    Q_OBJECT

private slots:

    void quadPulseUsesReferenceCapacityUnderDilation();
    void quad2PulseUsesReferenceCapacityUnderDilation();
    void tetraPulseUsesReferenceCapacityUnderDilation();
    void degenerateReferenceGeometryDoesNotPublishCandidate();
};

void TestStructuralReferenceThermalMass::quadPulseUsesReferenceCapacityUnderDilation()
{
    // The unit QUAD1 has total capacity 1 J/K. A one joule point pulse gives
    // sum(delta T)=4 K for its four equal row-sum weights. A 2x in-plane
    // dilation would incorrectly reduce that sum to 1 K if current area were
    // used for capacity.
    const PulseResponse referenceResult = solvePulse(CellType::Quad,true);
    const PulseResponse currentResult = solvePulse(CellType::Quad,false);
    QVERIFY(std::fabs(referenceResult.nodalIncrementSum-4.0) < 1.0e-7);
    QVERIFY(std::fabs(currentResult.nodalIncrementSum-1.0) < 1.0e-7);
    QVERIFY(std::fabs(referenceResult.referenceWeightedEnergy-1.0) < 1.0e-8);
    QVERIFY(std::fabs(currentResult.referenceWeightedEnergy-0.25) < 1.0e-8);
}

void TestStructuralReferenceThermalMass::quad2PulseUsesReferenceCapacityUnderDilation()
{
    const PulseResponse referenceResult = solvePulse(CellType::Quad2,true);
    const PulseResponse currentResult = solvePulse(CellType::Quad2,false);
    QVERIFY(std::fabs(referenceResult.referenceWeightedEnergy-1.0) < 1.0e-8);
    QVERIFY(std::fabs(currentResult.referenceWeightedEnergy-0.25) < 1.0e-8);
}

void TestStructuralReferenceThermalMass::tetraPulseUsesReferenceCapacityUnderDilation()
{
    // The unit-coordinate TETRA1 has volume 1/6 m^3, hence total capacity
    // 1/6 J/K and sum(delta T)=24 K for one joule. A 2x isotropic dilation
    // would incorrectly reduce this sum to 3 K when current volume is used.
    const PulseResponse referenceResult = solvePulse(CellType::Tetra,true);
    const PulseResponse currentResult = solvePulse(CellType::Tetra,false);
    QVERIFY(std::fabs(referenceResult.nodalIncrementSum-24.0) < 1.0e-7);
    QVERIFY(std::fabs(currentResult.nodalIncrementSum-3.0) < 1.0e-7);
    QVERIFY(std::fabs(referenceResult.referenceWeightedEnergy-1.0) < 1.0e-8);
    QVERIFY(std::fabs(currentResult.referenceWeightedEnergy-0.125) < 1.0e-8);
}

void TestStructuralReferenceThermalMass::degenerateReferenceGeometryDoesNotPublishCandidate()
{
    RModel model = makePulseModel(CellType::Quad,true,true);
    const uint temperaturePosition = model.findVariable(R_VARIABLE_TEMPERATURE);
    QVERIFY(temperaturePosition != RConstants::eod);
    RSolverSharedData sharedData;
    RSolverHeat heat(&model,QString(),QString(),sharedData);
    const uint sourceElement = model.getPoint(0).get(0);
    heat.setControlledPointHeat({{sourceElement,1.0}});
    heat.setReferenceSurfaceVolumeThermalMass(true);

    QString errorMessage;
    try
    {
        heat.run(false,0);
    }
    catch (const RError &error)
    {
        errorMessage = error.getMessage();
    }

    QVERIFY2(errorMessage.contains("degenerate",Qt::CaseInsensitive),
             qPrintable(errorMessage));
    QVERIFY(!heat.hasConverged());
    const RVariable &temperature = model.getVariable(temperaturePosition);
    for (uint node=0;node<model.getNNodes();node++)
    {
        QCOMPARE(temperature.getValue(0,node),300.0);
    }
}

QTEST_APPLESS_MAIN(TestStructuralReferenceThermalMass)
#include "tst_structural_reference_thermal_mass.moc"

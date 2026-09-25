#include <QtTest>

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include <rbl_error.h>
#include <rml_boundary_condition.h>
#include <rml_environment_condition.h>
#include <rml_model.h>
#include <rml_problem_task_item.h>
#include <rml_variable.h>

#include "rsolver.h"

namespace
{
constexpr double beamLength = 0.30;
constexpr double beamHeight = 0.01;
constexpr double bottomTemperature = 300.0;
constexpr double topTemperature = 400.0;
constexpr double environmentTemperature = 0.5*(bottomTemperature+topTemperature);
constexpr double youngsModulus = 70.0e9;
constexpr double poissonRatio = 0.30;
constexpr double thermalExpansion = 1.2e-5;

template <typename Condition>
void setBoundaryValue(Condition &condition,RVariableType type,double value)
{
    const uint position = condition.findComponentPosition(type);
    if (position == RConstants::eod)
    {
        throw std::logic_error("Thermal bending test boundary condition is missing a component.");
    }
    RConditionComponent &component = condition.getComponent(position);
    component.setEnabled(true);
    component.clear();
    component.add(0.0,value);
}

void setDisabled(RBoundaryCondition &condition,RVariableType type)
{
    const uint position = condition.findComponentPosition(type);
    if (position == RConstants::eod)
    {
        throw std::logic_error("Thermal bending test boundary condition is missing a component.");
    }
    condition.getComponent(position).setEnabled(false);
}

RModel makeCantilever(uint axialSlices,uint thicknessLayers,double elementPoissonRatio)
{
    if (axialSlices < 2 || thicknessLayers < 2)
    {
        throw std::invalid_argument("Thermal bending mesh requires at least two cells in each direction.");
    }

    // A QUAD2 cell uses a shared 3x3 grid of corner, edge-midpoint, and
    // center nodes. The quadratic x interpolation represents beam curvature
    // without the severe shear locking observed with the first-order QUAD1.
    const uint nodesPerRow = 2*axialSlices+1;
    const uint nodeRows = 2*thicknessLayers+1;
    const uint nodeCount = nodesPerRow*nodeRows;
    RModel model;
    model.setNNodes(nodeCount);
    for (uint j=0;j<nodeRows;j++)
    {
        const double y = -0.5*beamHeight + beamHeight*double(j)/double(nodeRows-1);
        for (uint i=0;i<nodesPerRow;i++)
        {
            const double x = beamLength*double(i)/double(nodesPerRow-1);
            model.getNode(j*nodesPerRow+i).set(x,y,0.0);
        }
    }

    for (uint j=0;j<thicknessLayers;j++)
    {
        for (uint i=0;i<axialSlices;i++)
        {
            const uint bottomLeft = (2*j)*nodesPerRow+2*i;
            const uint bottomRight = bottomLeft+2;
            const uint topLeft = bottomLeft+2*nodesPerRow;
            const uint topRight = topLeft+2;
            RElement element(R_ELEMENT_QUAD2);
            element.setNodeId(0,bottomLeft);
            element.setNodeId(1,bottomRight);
            element.setNodeId(2,topRight);
            element.setNodeId(3,topLeft);
            element.setNodeId(4,bottomLeft+1);
            element.setNodeId(5,bottomRight+nodesPerRow);
            element.setNodeId(6,topLeft+1);
            element.setNodeId(7,bottomLeft+nodesPerRow);
            element.setNodeId(8,bottomLeft+nodesPerRow+1);
            model.addElement(element,true,0);
        }
    }

    RMaterial material(RMaterial::Solid);
    material.setName("Thermal bending plane-stress material");
    auto addProperty = [&material](RMaterialProperty::Type type,double value) {
        RMaterialProperty property(type);
        // RMaterialProperty starts with a default sample at the temperature
        // variable's initial value (293.15 K). Replace that sample so a
        // constant property remains constant over the 300–400 K field.
        property.add(RVariable::getInitValue(R_VARIABLE_TEMPERATURE),value);
        material.add(property);
    };
    addProperty(RMaterialProperty::ModulusOfElasiticity,youngsModulus);
    addProperty(RMaterialProperty::PoissonRatio,elementPoissonRatio);
    addProperty(RMaterialProperty::Density,2700.0);
    addProperty(RMaterialProperty::ThermalExpansionCoefficient,thermalExpansion);
    addProperty(RMaterialProperty::ThermalConductivity,1.0);
    addProperty(RMaterialProperty::HeatCapacity,900.0);
    addProperty(RMaterialProperty::Emissivity,0.0);

    RSurface &surface = model.getSurface(0);
    surface.setName("Thermal bending cantilever");
    surface.setThickness(1.0);
    surface.setStructuralFormulation(RElementFormulation::PlaneStress);
    surface.setMaterial(material);

    REnvironmentCondition environment(R_ENVIRONMENT_CONDITION_TEMPERATURE);
    setBoundaryValue(environment,R_VARIABLE_TEMPERATURE,environmentTemperature);
    surface.addEnvironmentCondition(environment);

    // A 2-D surface has no out-of-plane bending degree of freedom in this
    // formulation; restrain only that component and leave both in-plane
    // directions free away from the clamped face.
    RBoundaryCondition planar(R_BOUNDARY_CONDITION_DISPLACEMENT);
    setDisabled(planar,R_VARIABLE_DISPLACEMENT_X);
    setDisabled(planar,R_VARIABLE_DISPLACEMENT_Y);
    setBoundaryValue(planar,R_VARIABLE_DISPLACEMENT_Z,0.0);
    surface.addBoundaryCondition(planar);

    // Match the model convention of one POINT group per node. A corner owns
    // both its thermal and mechanical conditions on that single group; do not
    // create duplicate point elements for the same physical node.
    auto addPointToGroup = [&model](uint node,uint group) {
        RElement point(R_ELEMENT_POINT);
        point.setNodeId(0,node);
        model.addElement(point,true,group);
        model.getPoint(group).setVolume(0.0);
    };
    for (uint j=0;j<nodeRows;j++)
    {
        for (uint i=0;i<nodesPerRow;i++)
        {
            const uint node = j*nodesPerRow+i;
            addPointToGroup(node,node);
            if (j == 0)
            {
                RBoundaryCondition bottom(R_BOUNDARY_CONDITION_TEMPERATURE);
                setBoundaryValue(bottom,R_VARIABLE_TEMPERATURE,bottomTemperature);
                model.getPoint(node).addBoundaryCondition(bottom);
            }
            if (j+1 == nodeRows)
            {
                RBoundaryCondition top(R_BOUNDARY_CONDITION_TEMPERATURE);
                setBoundaryValue(top,R_VARIABLE_TEMPERATURE,topTemperature);
                model.getPoint(node).addBoundaryCondition(top);
            }
            if (i == 0)
            {
                RBoundaryCondition clamp(R_BOUNDARY_CONDITION_DISPLACEMENT);
                setBoundaryValue(clamp,R_VARIABLE_DISPLACEMENT_X,0.0);
                setBoundaryValue(clamp,R_VARIABLE_DISPLACEMENT_Y,0.0);
                setBoundaryValue(clamp,R_VARIABLE_DISPLACEMENT_Z,0.0);
                model.getPoint(node).addBoundaryCondition(clamp);
            }
        }
    }

    RProblemTaskItem sequence;
    sequence.addChild(RProblemTaskItem(R_PROBLEM_HEAT));
    sequence.addChild(RProblemTaskItem(R_PROBLEM_STRESS));
    model.setProblemTaskTree(sequence);
    model.getTimeSolver().setEnabled(false);
    RMatrixSolverConf &cg = model.getMatrixSolverConf(RMatrixSolverConf::CG);
    cg.setNInnerIterations(500);
    cg.setNOuterIterations(5000);
    cg.setSolverCvgValue(1.0e-12);
    cg.setOutputFrequency(0);
    return model;
}

double value(const RModel &model,RVariableType type,uint component,uint entity)
{
    const uint position = model.findVariable(type);
    if (position == RConstants::eod)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return model.getVariable(position).getValue(component,entity);
}

struct ThermalBendingResponse
{
    double tipDeflection;
    double interiorCurvature;
    double minimumElementTemperature;
    double maximumElementTemperature;
    double minimumElementTemperatureDelta;
    double maximumElementTemperatureDelta;
    double minimumMeanAlpha;
    double maximumMeanAlpha;
};

ThermalBendingResponse solveThermalBending(uint axialSlices,uint thicknessLayers,
                                           double elementPoissonRatio = poissonRatio)
{
    RModel model = makeCantilever(axialSlices,thicknessLayers,elementPoissonRatio);
    RSolver solver(model,QString(),QString());
    solver.run();

    const uint nodesPerRow = 2*axialSlices+1;
    const uint nodeRows = 2*thicknessLayers+1;
    for (uint j=0;j<nodeRows;j++)
    {
        const double expectedTemperature = bottomTemperature
            + (topTemperature-bottomTemperature)*double(j)/double(nodeRows-1);
        for (uint i=0;i<nodesPerRow;i++)
        {
            const uint node = j*nodesPerRow+i;
            const double temperature = value(model,R_VARIABLE_TEMPERATURE,0,node);
            if (!std::isfinite(temperature) || std::fabs(temperature-expectedTemperature) > 1.0e-4)
            {
                throw std::runtime_error(QString("Surface heat profile mismatch at node %1 (i=%2, j=%3): got %4 K, expected %5 K.")
                                         .arg(node).arg(i).arg(j).arg(temperature,0,'g',16)
                                         .arg(expectedTemperature,0,'g',16).toStdString());
            }
        }
    }

    const uint centerlineRow = thicknessLayers;
    const uint halfSpanNode = centerlineRow*nodesPerRow+axialSlices;
    const uint threeQuarterSpanNode = centerlineRow*nodesPerRow+(3*axialSlices)/2;
    const uint tipCenterNode = centerlineRow*nodesPerRow+2*axialSlices;
    const double halfSpanDeflection = value(model,R_VARIABLE_DISPLACEMENT,1,halfSpanNode);
    const double threeQuarterSpanDeflection = value(model,R_VARIABLE_DISPLACEMENT,1,threeQuarterSpanNode);
    const double tipDeflection = value(model,R_VARIABLE_DISPLACEMENT,1,tipCenterNode);
    const double spanStep = 0.25*beamLength;
    const double interiorCurvature = -(tipDeflection-2.0*threeQuarterSpanDeflection+halfSpanDeflection)
                                   / (spanStep*spanStep);

    // Reconstruct the element data sent to the structural solve from the
    // checked nodal field. QUAD2's 3x3 nodal mean equals its center value for
    // this linear profile, so these bounds expose any temperature-transfer
    // or thermal-load scaling discrepancy at the element level.
    double minimumElementTemperature = std::numeric_limits<double>::infinity();
    double maximumElementTemperature = -std::numeric_limits<double>::infinity();
    const RSurface &surface = model.getSurface(0);
    for (uint position=0;position<surface.size();position++)
    {
        const RElement &element = model.getElement(surface.get(position));
        double elementTemperature = 0.0;
        for (uint localNode=0;localNode<element.size();localNode++)
        {
            elementTemperature += value(model,R_VARIABLE_TEMPERATURE,0,element.getNodeId(localNode));
        }
        elementTemperature /= double(element.size());
        minimumElementTemperature = std::min(minimumElementTemperature,elementTemperature);
        maximumElementTemperature = std::max(maximumElementTemperature,elementTemperature);
    }
    const RMaterialProperty &alphaProperty = surface.getMaterial().get(
        surface.getMaterial().findPosition(RMaterialProperty::ThermalExpansionCoefficient));
    const double minimumDelta = minimumElementTemperature-environmentTemperature;
    const double maximumDelta = maximumElementTemperature-environmentTemperature;
    const double minimumAlpha = alphaProperty.get(minimumElementTemperature);
    const double maximumAlpha = alphaProperty.get(maximumElementTemperature);
    return {tipDeflection,interiorCurvature,minimumElementTemperature,maximumElementTemperature,
            minimumDelta,maximumDelta,minimumAlpha,maximumAlpha};
}
}

class TestStructuralThermalBending : public QObject
{
    Q_OBJECT
private slots:
    void planeStressCantileverMatchesThermalCurvatureAcrossMeshes();
};

void TestStructuralThermalBending::planeStressCantileverMatchesThermalCurvatureAcrossMeshes()
{
    constexpr double gradient = (topTemperature-bottomTemperature)/beamHeight;
    constexpr double curvature = thermalExpansion*gradient;
    const double expectedTip = -0.5*curvature*beamLength*beamLength;

    QString error;
    ThermalBendingResponse coarse{};
    ThermalBendingResponse fine{};
    try
    {
        coarse = solveThermalBending(24,6);
        fine = solveThermalBending(48,12);
    }
    catch (const RError &solverError)
    {
        error = solverError.getMessage();
    }
    catch (const std::exception &exception)
    {
        error = QString::fromUtf8(exception.what());
    }
    QVERIFY2(error.isEmpty(),qPrintable(error));

    ThermalBendingResponse fineZeroPoisson{};
    try
    {
        fineZeroPoisson = solveThermalBending(48,12,0.0);
    }
    catch (const RError &solverError)
    {
        error = solverError.getMessage();
    }
    catch (const std::exception &exception)
    {
        error = QString::fromUtf8(exception.what());
    }
    QVERIFY2(error.isEmpty(),qPrintable(error));
    qInfo().noquote() << QString("Thermal FE diagnostic: E=%1 Pa, alpha=1.2e-5 1/K, env=350 K, "
                                 "fine element T=[%2,%3] K, deltaT=[%4,%5] K, mean alpha=[%6,%7] 1/K; "
                                 "coarse 24x6: kappa=%8 1/m (%9%%), tip=%10 m (%11%%); "
                                 "fine 48x12: kappa=%12 1/m (%13%%), tip=%14 m (%15%%); "
                                 "mesh tip change=%16%%; nu=0: kappa=%17 1/m, tip=%18 m")
                         .arg(youngsModulus,0,'g',8)
                         .arg(fine.minimumElementTemperature,0,'g',12)
                         .arg(fine.maximumElementTemperature,0,'g',12)
                         .arg(fine.minimumElementTemperatureDelta,0,'g',12)
                         .arg(fine.maximumElementTemperatureDelta,0,'g',12)
                         .arg(fine.minimumMeanAlpha,0,'g',12)
                         .arg(fine.maximumMeanAlpha,0,'g',12)
                         .arg(coarse.interiorCurvature,0,'g',12)
                         .arg(100.0*std::fabs((coarse.interiorCurvature-curvature)/curvature),0,'g',8)
                         .arg(coarse.tipDeflection,0,'g',12)
                         .arg(100.0*std::fabs((coarse.tipDeflection-expectedTip)/expectedTip),0,'g',8)
                         .arg(fine.interiorCurvature,0,'g',12)
                         .arg(100.0*std::fabs((fine.interiorCurvature-curvature)/curvature),0,'g',8)
                         .arg(fine.tipDeflection,0,'g',12)
                         .arg(100.0*std::fabs((fine.tipDeflection-expectedTip)/expectedTip),0,'g',8)
                         .arg(100.0*std::fabs((fine.tipDeflection-coarse.tipDeflection)/expectedTip),0,'g',8)
                         .arg(fineZeroPoisson.interiorCurvature,0,'g',12)
                         .arg(fineZeroPoisson.tipDeflection,0,'g',12);

    QVERIFY2(std::isfinite(coarse.tipDeflection) && std::isfinite(fine.tipDeflection) &&
             std::isfinite(coarse.interiorCurvature) && std::isfinite(fine.interiorCurvature),
             "The heat-stress FE solution returned a non-finite tip displacement.");
    QVERIFY2(coarse.tipDeflection*expectedTip > 0.0 && fine.tipDeflection*expectedTip > 0.0 &&
             coarse.interiorCurvature*curvature > 0.0 && fine.interiorCurvature*curvature > 0.0,
             qPrintable(QString("Tip deflection signs do not follow thermal curvature: coarse=%1, fine=%2, expected=%3")
                        .arg(coarse.tipDeflection,0,'g',16).arg(fine.tipDeflection,0,'g',16).arg(expectedTip,0,'g',16)));

    const double coarseTipError = std::fabs((coarse.tipDeflection-expectedTip)/expectedTip);
    const double fineTipError = std::fabs((fine.tipDeflection-expectedTip)/expectedTip);
    const double coarseCurvatureError = std::fabs((coarse.interiorCurvature-curvature)/curvature);
    const double fineCurvatureError = std::fabs((fine.interiorCurvature-curvature)/curvature);
    const double refinementChange = std::fabs((fine.tipDeflection-coarse.tipDeflection)/expectedTip);
    QVERIFY2(fineCurvatureError < 0.05,
             qPrintable(QString("Fine interior curvature %1 differs from alpha*gradient %2 by %3%%")
                        .arg(fine.interiorCurvature,0,'g',16).arg(curvature,0,'g',16)
                        .arg(100.0*fineCurvatureError)));
    QVERIFY2(coarseCurvatureError < 0.10,
             qPrintable(QString("Coarse interior curvature %1 differs from alpha*gradient %2 by %3%%")
                        .arg(coarse.interiorCurvature,0,'g',16).arg(curvature,0,'g',16)
                        .arg(100.0*coarseCurvatureError)));
    QVERIFY2(coarseTipError < 0.15,
             qPrintable(QString("Coarse thermal bending tip %1 differs from Euler-Bernoulli value %2 by %3%%")
                        .arg(coarse.tipDeflection,0,'g',16).arg(expectedTip,0,'g',16).arg(100.0*coarseTipError)));
    QVERIFY2(fineTipError < 0.10,
             qPrintable(QString("Fine thermal bending tip %1 differs from Euler-Bernoulli value %2 by %3%%")
                        .arg(fine.tipDeflection,0,'g',16).arg(expectedTip,0,'g',16).arg(100.0*fineTipError)));
    QVERIFY2(refinementChange < 0.08,
             qPrintable(QString("Tip displacement changed by %1%% between the two mesh resolutions")
                        .arg(100.0*refinementChange)));
}

QTEST_APPLESS_MAIN(TestStructuralThermalBending)
#include "tst_structural_thermal_bending.moc"

#include <QtTest>

#include <cmath>
#include <exception>
#include <limits>

#include "rml_model.h"
#include "rsolver.h"

namespace
{
constexpr double youngModulus = 2.0e11;
constexpr double density = 8000.0;
constexpr double beamLength = 2.0;
constexpr double area = 0.01;
constexpr double iy = 2.0e-5;
constexpr double iz = 5.0e-5;
constexpr uint slices = 3;

RModel makeCantilever()
{
    RModel model;
    model.setNNodes(slices+1);
    for (uint node=0;node<=slices;node++)
    {
        model.getNode(node).set(beamLength*double(node)/double(slices),0.0,0.0);
    }
    for (uint elementId=0;elementId<slices;elementId++)
    {
        RElement element(R_ELEMENT_BEAM1);
        element.setNodeId(0,elementId);
        element.setNodeId(1,elementId+1);
        model.addElement(element,true,0);
    }

    RMaterial material(RMaterial::Solid);
    RMaterialProperty modulus(RMaterialProperty::ModulusOfElasiticity);
    modulus.add(293.15,youngModulus);
    material.add(modulus);
    RMaterialProperty poisson(RMaterialProperty::PoissonRatio);
    poisson.add(293.15,0.3);
    material.add(poisson);
    RMaterialProperty massDensity(RMaterialProperty::Density);
    massDensity.add(293.15,density);
    material.add(massDensity);
    RMaterialProperty expansion(RMaterialProperty::ThermalExpansionCoefficient);
    expansion.add(293.15,0.0);
    material.add(expansion);

    RLine &line = model.getLine(0);
    line.setCrossArea(area);
    line.setSecondMomentY(iy);
    line.setSecondMomentZ(iz);
    line.setTorsionalConstant(iy+iz);
    line.setSectionModulusY(4.0e-4);
    line.setSectionModulusZ(5.0e-4);
    line.setTorsionalSectionModulus(6.0e-4);
    line.setLocalYDirection(RR3Vector(0.0,1.0,0.0));
    line.setMaterial(material);

    RElement root(R_ELEMENT_POINT);
    root.setNodeId(0,0);
    model.addElement(root,true,0);
    RBoundaryCondition displacement(R_BOUNDARY_CONDITION_DISPLACEMENT);
    displacement.getComponent(displacement.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).add(0.0,0.0);
    displacement.getComponent(displacement.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).add(0.0,0.0);
    displacement.getComponent(displacement.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).add(0.0,0.0);
    model.getPoint(0).addBoundaryCondition(displacement);
    RBoundaryCondition rotation(R_BOUNDARY_CONDITION_ROTATION);
    rotation.getComponent(rotation.findComponentPosition(R_VARIABLE_ROTATION_X)).add(0.0,0.0);
    rotation.getComponent(rotation.findComponentPosition(R_VARIABLE_ROTATION_Y)).add(0.0,0.0);
    rotation.getComponent(rotation.findComponentPosition(R_VARIABLE_ROTATION_Z)).add(0.0,0.0);
    model.getPoint(0).addBoundaryCondition(rotation);

    model.getProblemSetup().setMassMatrixType(RMassMatrixType::Consistent);
    RMatrixSolverConf &linear = model.getMatrixSolverConf(RMatrixSolverConf::CG);
    linear.setNInnerIterations(200);
    linear.setNOuterIterations(2000);
    linear.setSolverCvgValue(1.0e-12);
    linear.setOutputFrequency(0);
    return model;
}

bool runSolver(RModel &model, QString &error)
{
    try
    {
        RSolver solver(model,QString(),QString());
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

double value(const RModel &model, RVariableType type, uint component, uint node)
{
    const uint index = model.findVariable(type);
    if (index == RConstants::eod)
    {
        return std::numeric_limits<double>::quiet_NaN();
    }
    return model.getVariable(index).getValue(component,node);
}
}

class TestStructuralDynamicsBeam : public QObject
{
    Q_OBJECT

private slots:
    void firstBendingModeMatchesHalfPeriodTransient();
};

void TestStructuralDynamicsBeam::firstBendingModeMatchesHalfPeriodTransient()
{
    RModel model = makeCantilever();
    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS_MODAL));
    RModalSetup &modal = model.getProblemSetup().getModalSetup();
    modal.setMethod(R_MODAL_MULTIPLE_MODES);
    modal.setNModesToExtract(1);
    modal.setNIterations(400);
    modal.setConvergenceValue(1.0e-11);
    QString error;
    QVERIFY2(runSolver(model,error),qPrintable(error));

    constexpr double firstRoot = 1.875104068711961;
    const double analyticalFrequency = firstRoot*firstRoot/(2.0*std::acos(-1.0))
        *std::sqrt(youngModulus*iy/(density*area*std::pow(beamLength,4)));
    const double frequency = modal.getFrequency();
    QVERIFY2(std::isfinite(frequency) && frequency > 0.0,"Modal frequency is not positive and finite.");
    QVERIFY2(std::abs(frequency-analyticalFrequency) < 0.05*analyticalFrequency,
             qPrintable(QString("Beam mode %1 Hz; continuum cantilever %2 Hz")
                        .arg(frequency).arg(analyticalFrequency)));

    // The modal displacement and rotation fields are retained as the dynamic
    // initial state. Starting on one eigenvector isolates its frequency.
    const double tipDisplacement0 = value(model,R_VARIABLE_DISPLACEMENT,2,slices);
    const double tipRotation0 = value(model,R_VARIABLE_ROTATION,1,slices);
    QVERIFY(std::isfinite(tipDisplacement0) && std::abs(tipDisplacement0) > 1.0e-12);
    QVERIFY(std::isfinite(tipRotation0) && std::abs(tipRotation0) > 1.0e-12);

    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    RStructuralDynamicsSetup &dynamics = model.getProblemSetup().getStructuralDynamicsSetup();
    dynamics.setEnabled(true);
    dynamics.setMethod(RStructuralTimeIntegrationMethod::Newmark);
    dynamics.setTimeStep(0.5/(40.0*frequency));
    dynamics.setNSteps(40);
    dynamics.setOutputFrequency(40);
    error.clear();
    QVERIFY2(runSolver(model,error),qPrintable(error));

    const double tipDisplacement = value(model,R_VARIABLE_DISPLACEMENT,2,slices);
    const double tipRotation = value(model,R_VARIABLE_ROTATION,1,slices);
    const double tipAngularVelocity = value(model,R_VARIABLE_ANGULAR_VELOCITY,1,slices);
    const double tipAngularAcceleration = value(model,R_VARIABLE_ANGULAR_ACCELERATION,1,slices);
    const double omega = 2.0*std::acos(-1.0)*frequency;
    QVERIFY(std::isfinite(tipDisplacement) && std::isfinite(tipRotation));
    QVERIFY(std::isfinite(tipAngularVelocity) && std::isfinite(tipAngularAcceleration));
    QVERIFY2(std::abs(tipDisplacement+tipDisplacement0) < 0.03*std::abs(tipDisplacement0),
             qPrintable(QString("Half-period tip displacement %1, initial %2")
                        .arg(tipDisplacement).arg(tipDisplacement0)));
    QVERIFY2(std::abs(tipRotation+tipRotation0) < 0.03*std::abs(tipRotation0),
             qPrintable(QString("Half-period tip rotation %1, initial %2")
                        .arg(tipRotation).arg(tipRotation0)));
    QVERIFY(std::abs(tipAngularVelocity) < 0.05*omega*std::abs(tipRotation0));
    QVERIFY(std::abs(tipAngularAcceleration+omega*omega*tipRotation)
            < 0.05*omega*omega*std::abs(tipRotation0));

    // The clamped end remains fixed in every stored kinematic field.
    for (RVariableType type : {R_VARIABLE_DISPLACEMENT,R_VARIABLE_ROTATION,
                               R_VARIABLE_VELOCITY,R_VARIABLE_ACCELERATION,
                               R_VARIABLE_ANGULAR_VELOCITY,R_VARIABLE_ANGULAR_ACCELERATION})
    {
        for (uint component=0;component<3;component++)
        {
            const double rootValue = value(model,type,component,0);
            QVERIFY(std::isfinite(rootValue));
            QVERIFY2(std::abs(rootValue) < 1.0e-8,
                     qPrintable(QString("Root %1 component %2 = %3")
                                .arg(int(type)).arg(component).arg(rootValue)));
        }
    }
}

QTEST_APPLESS_MAIN(TestStructuralDynamicsBeam)

#include "tst_structural_dynamics_beam.moc"

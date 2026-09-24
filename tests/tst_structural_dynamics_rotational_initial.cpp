#include <QtTest>

#include <cmath>
#include <exception>
#include <limits>

#include "rml_model.h"
#include "rsolver.h"

namespace
{
constexpr double youngModulus = 2.0e11;
constexpr double poissonRatio = 0.3;
constexpr double density = 8000.0;
constexpr double length = 2.0;
constexpr double area = 0.01;
constexpr double secondMomentY = 2.0e-5;
constexpr double secondMomentZ = 5.0e-5;
constexpr double polarConstant = secondMomentY+secondMomentZ;

// The clamped BEAM1 torsion branch is independent of bending and translation.
// Its free-end K=GJ/L and consistent M=rho*(Iy+Iz)*L/3.  Taking J=Iy+Iz
// therefore gives a one-DOF oscillator with an independently known omega.
double torsionFrequency()
{
    const double shearModulus = youngModulus/(2.0*(1.0+poissonRatio));
    const double stiffness = shearModulus*polarConstant/length;
    const double mass = density*(secondMomentY+secondMomentZ)*length/3.0;
    return std::sqrt(stiffness/mass);
}

RModel makeTorsionalCantilever(RInitialConditionType initialType, double initialValue)
{
    RModel model;
    model.addNode(RNode(0.0,0.0,0.0));
    model.addNode(RNode(length,0.0,0.0));
    RElement beam(R_ELEMENT_BEAM1);
    beam.setNodeId(0,0);
    beam.setNodeId(1,1);
    model.addElement(beam,true,0);

    RMaterial material(RMaterial::Solid);
    RMaterialProperty modulus(RMaterialProperty::ModulusOfElasiticity);
    modulus.add(293.15,youngModulus);
    material.add(modulus);
    RMaterialProperty poisson(RMaterialProperty::PoissonRatio);
    poisson.add(293.15,poissonRatio);
    material.add(poisson);
    RMaterialProperty massDensity(RMaterialProperty::Density);
    massDensity.add(293.15,density);
    material.add(massDensity);
    RMaterialProperty expansion(RMaterialProperty::ThermalExpansionCoefficient);
    expansion.add(293.15,0.0);
    material.add(expansion);

    RLine &line = model.getLine(0);
    line.setCrossArea(area);
    line.setSecondMomentY(secondMomentY);
    line.setSecondMomentZ(secondMomentZ);
    line.setTorsionalConstant(polarConstant);
    line.setSectionModulusY(4.0e-4);
    line.setSectionModulusZ(5.0e-4);
    line.setTorsionalSectionModulus(6.0e-4);
    line.setLocalYDirection(RR3Vector(0.0,1.0,0.0));
    line.setMaterial(material);

    RElement root(R_ELEMENT_POINT);
    root.setNodeId(0,0);
    model.addElement(root,true,0);
    RBoundaryCondition heldDisplacement(R_BOUNDARY_CONDITION_DISPLACEMENT);
    heldDisplacement.getComponent(heldDisplacement.findComponentPosition(R_VARIABLE_DISPLACEMENT_X)).add(0.0,0.0);
    heldDisplacement.getComponent(heldDisplacement.findComponentPosition(R_VARIABLE_DISPLACEMENT_Y)).add(0.0,0.0);
    heldDisplacement.getComponent(heldDisplacement.findComponentPosition(R_VARIABLE_DISPLACEMENT_Z)).add(0.0,0.0);
    model.getPoint(0).addBoundaryCondition(heldDisplacement);
    RBoundaryCondition heldRotation(R_BOUNDARY_CONDITION_ROTATION);
    heldRotation.getComponent(heldRotation.findComponentPosition(R_VARIABLE_ROTATION_X)).add(0.0,0.0);
    heldRotation.getComponent(heldRotation.findComponentPosition(R_VARIABLE_ROTATION_Y)).add(0.0,0.0);
    heldRotation.getComponent(heldRotation.findComponentPosition(R_VARIABLE_ROTATION_Z)).add(0.0,0.0);
    model.getPoint(0).addBoundaryCondition(heldRotation);

    RElement tip(R_ELEMENT_POINT);
    tip.setNodeId(0,1);
    model.addElement(tip,true,1);
    RInitialCondition initial(initialType);
    const RVariableType component = initialType == R_INITIAL_CONDITION_ROTATION
        ? R_VARIABLE_ROTATION_X : R_VARIABLE_ANGULAR_VELOCITY_X;
    initial.getComponent(initial.findComponentPosition(component)).add(0.0,initialValue);
    model.getPoint(1).addInitialCondition(initial);

    model.setProblemTaskTree(RProblemTaskItem(R_PROBLEM_STRESS));
    model.getProblemSetup().setMassMatrixType(RMassMatrixType::Consistent);
    RMatrixSolverConf &linear = model.getMatrixSolverConf(RMatrixSolverConf::CG);
    linear.setNInnerIterations(200);
    linear.setNOuterIterations(2000);
    linear.setSolverCvgValue(1.0e-13);
    linear.setOutputFrequency(0);
    RStructuralDynamicsSetup &dynamics = model.getProblemSetup().getStructuralDynamicsSetup();
    dynamics.setEnabled(true);
    dynamics.setMethod(RStructuralTimeIntegrationMethod::Newmark);
    dynamics.setTimeStep(0.1/torsionFrequency());
    dynamics.setNSteps(1);
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
    const uint position = model.findVariable(type);
    return position == RConstants::eod ? std::numeric_limits<double>::quiet_NaN()
                                     : model.getVariable(position).getValue(component,node);
}

void checkClose(double actual, double expected, double scale, const char *name)
{
    QVERIFY2(std::isfinite(actual),name);
    QVERIFY2(std::abs(actual-expected) <= 2.0e-7*scale,
             qPrintable(QString("%1: actual=%2 expected=%3")
                        .arg(QString::fromLatin1(name)).arg(actual,0,'g',17).arg(expected,0,'g',17)));
}
}

class TestStructuralDynamicsRotationalInitial : public QObject
{
    Q_OBJECT

private slots:
    void initialRotationEntersTorsionalDisplacement();
    void initialAngularVelocityEntersTorsionalVelocity();
};

void TestStructuralDynamicsRotationalInitial::initialRotationEntersTorsionalDisplacement()
{
    constexpr double theta0 = 1.0e-3;
    RModel model = makeTorsionalCantilever(R_INITIAL_CONDITION_ROTATION,theta0);
    const double omega = torsionFrequency();
    const double dt = model.getProblemSetup().getStructuralDynamicsSetup().getTimeStep();
    const double x2 = omega*omega*dt*dt;
    const double denominator = 1.0+0.25*x2;
    const double expectedTheta = theta0*(1.0-0.25*x2)/denominator;
    const double expectedVelocity = -omega*omega*dt*theta0/denominator;

    QString error;
    QVERIFY2(runSolver(model,error),qPrintable(error));
    checkClose(value(model,R_VARIABLE_ROTATION,0,1),expectedTheta,theta0,"tip rotation");
    checkClose(value(model,R_VARIABLE_ANGULAR_VELOCITY,0,1),expectedVelocity,
               omega*theta0,"tip angular velocity");
    checkClose(value(model,R_VARIABLE_ANGULAR_ACCELERATION,0,1),
               -omega*omega*expectedTheta,omega*omega*theta0,"tip angular acceleration");
    checkClose(value(model,R_VARIABLE_ROTATION,0,0),0.0,theta0,"fixed-end rotation");
}

void TestStructuralDynamicsRotationalInitial::initialAngularVelocityEntersTorsionalVelocity()
{
    constexpr double angularVelocity0 = 0.25;
    RModel model = makeTorsionalCantilever(R_INITIAL_CONDITION_ANGULAR_VELOCITY,
                                           angularVelocity0);
    const double omega = torsionFrequency();
    const double dt = model.getProblemSetup().getStructuralDynamicsSetup().getTimeStep();
    const double x2 = omega*omega*dt*dt;
    const double denominator = 1.0+0.25*x2;
    const double expectedTheta = dt*angularVelocity0/denominator;
    const double expectedVelocity = angularVelocity0*(1.0-0.25*x2)/denominator;

    QString error;
    QVERIFY2(runSolver(model,error),qPrintable(error));
    checkClose(value(model,R_VARIABLE_ROTATION,0,1),expectedTheta,
               dt*angularVelocity0,"tip rotation from initial velocity");
    checkClose(value(model,R_VARIABLE_ANGULAR_VELOCITY,0,1),expectedVelocity,
               angularVelocity0,"tip angular velocity");
    checkClose(value(model,R_VARIABLE_ANGULAR_ACCELERATION,0,1),
               -omega*omega*expectedTheta,omega*angularVelocity0,"tip angular acceleration");
    checkClose(value(model,R_VARIABLE_ANGULAR_VELOCITY,0,0),0.0,
               angularVelocity0,"fixed-end angular velocity");
}

QTEST_APPLESS_MAIN(TestStructuralDynamicsRotationalInitial)

#include "tst_structural_dynamics_rotational_initial.moc"

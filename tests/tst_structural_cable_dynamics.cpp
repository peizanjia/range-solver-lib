#include <QtTest>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "rstructural_cable_dynamics.h"

namespace
{
using Node = RStructuralCableVector3;
using Member = RStructuralCableEquilibriumMember;
using Result = RStructuralCableDynamicsResult;

std::vector<bool> fixedMask(std::size_t nodeCount)
{
    return std::vector<bool>(3*nodeCount,true);
}

std::vector<double> internalForces(const std::vector<Node> &nodes,
                                   const std::vector<Member> &members,
                                   const std::vector<double> &displacements)
{
    std::vector<double> forces(displacements.size(),0.0);
    for (const auto &member : members)
    {
        Node node1{},node2{};
        for (std::size_t axis=0;axis<3;axis++)
        {
            node1[axis] = nodes[member.node1][axis]+displacements[member.node1*3+axis];
            node2[axis] = nodes[member.node2][axis]+displacements[member.node2*3+axis];
        }
        const auto response = RStructuralCableKernel::evaluate(
            node1,node2,member.referenceLength,member.initialTension,member.axialRigidity);
        const std::size_t dofs[6] = {
            member.node1*3,member.node1*3+1,member.node1*3+2,
            member.node2*3,member.node2*3+1,member.node2*3+2
        };
        for (std::size_t local=0;local<6;local++)
            forces[dofs[local]] += response.internalForce[local];
    }
    return forces;
}

double freeResidualNorm(const std::vector<double> &forces,
                        const std::vector<double> &accelerations,
                        const std::vector<double> &masses,
                        const std::vector<double> &loads,
                        const std::vector<bool> &prescribed)
{
    double norm = 0.0;
    for (std::size_t dof=0;dof<forces.size();dof++)
    {
        if (prescribed[dof]) continue;
        const double residual = forces[dof]+masses[dof/3]*accelerations[dof]-loads[dof];
        norm = std::hypot(norm,residual);
    }
    return norm;
}

RStructuralCableDynamicsOptions optionsFor(double dt)
{
    RStructuralCableDynamicsOptions options;
    options.dt = dt;
    options.absForceTol = 1.0e-12;
    options.relForceTol = 1.0e-12;
    return options;
}

double integrateLinearOscillator(double dt)
{
    const std::vector<Node> nodes = {{{0.0,0.0,0.0}},{{1.0,0.0,0.0}}};
    const std::vector<Member> members = {{0,1,1.0,10.0,100.0}};
    const std::vector<double> masses = {0.0,1.0};
    std::vector<double> displacement(6,0.0),velocity(6,0.0),acceleration(6,0.0);
    displacement[3] = 1.0e-3;
    acceleration[3] = -0.1; // Independent linearized acceleration: -EA/L * u / m.
    std::vector<bool> prescribed = fixedMask(nodes.size());
    prescribed[3] = false;
    const std::vector<double> prescribedValues(6,0.0);
    std::vector<double> loads(6,0.0);
    loads[3] = 10.0; // Balances the cable's N0 at the zero-displacement equilibrium.
    constexpr double finalTime = 0.2;
    const std::size_t steps = std::size_t(std::llround(finalTime/dt));
    for (std::size_t i=0;i<steps;i++)
    {
        const Result result = RStructuralCableDynamicsSolver::step(
            nodes,members,masses,displacement,velocity,acceleration,prescribed,
            prescribedValues,loads,optionsFor(dt));
        displacement = result.displacements;
        velocity = result.velocities;
        acceleration = result.accelerations;
    }
    return displacement[3];
}
}

class TestStructuralCableDynamics : public QObject
{
    Q_OBJECT

private slots:
    void smallOscillationAboutBalancedPretensionMatchesLinearizedResponse();
    void prescribedMotionKeepsExactKinematicsAndSlackCanRetension();
    void failedStepLeavesAcceptedStateUnchanged();
    void timeStepRefinementConvergesForAxialOscillator();
    void unforcedSymmetricCableOscillationBalancesEnergy();
};

void TestStructuralCableDynamics::smallOscillationAboutBalancedPretensionMatchesLinearizedResponse()
{
    constexpr double halfSpan = 1.0;
    constexpr double sag = 0.4;
    constexpr double referenceLength = 1.0;
    constexpr double initialTension = 10.0;
    constexpr double axialRigidity = 100.0;
    constexpr double mass = 1.0;
    constexpr double perturbation = 1.0e-5;
    const std::vector<Node> nodes = {
        {{-halfSpan,0.0,0.0}},{{0.0,0.0,0.0}},{{halfSpan,0.0,0.0}}
    };
    const std::vector<Member> members = {
        {0,1,referenceLength,initialTension,axialRigidity},
        {1,2,referenceLength,initialTension,axialRigidity}
    };
    std::vector<double> displacement(9,0.0),velocity(9,0.0),acceleration(9,0.0);
    displacement[4] = -sag+perturbation;
    std::vector<bool> prescribed = fixedMask(nodes.size());
    prescribed[4] = false;
    const std::vector<double> prescribedValues(9,0.0);
    const double currentLength = std::hypot(halfSpan,sag);
    const double balancedTension = initialTension+
        axialRigidity/referenceLength*(currentLength-referenceLength);
    const double balancedLoad = -2.0*balancedTension*sag/currentLength;
    std::vector<double> loads(9,0.0);
    loads[4] = balancedLoad;

    const std::vector<double> previousInternal = internalForces(nodes,members,displacement);
    acceleration[4] = -(previousInternal[4]-loads[4])/mass;
    const double directionY = sag/currentLength;
    const double oneCableTangent = axialRigidity/referenceLength*directionY*directionY+
        (balancedTension/currentLength)*(1.0-directionY*directionY);
    const double tangent = 2.0*oneCableTangent;
    const double dt = 0.03;
    const auto options = optionsFor(dt);
    const double predictorPerturbation = perturbation+
        dt*dt*(0.5-options.beta)*acceleration[4];
    const double inertialTangent = mass/(options.beta*dt*dt);
    const double expectedPerturbation = inertialTangent*predictorPerturbation/
        (tangent+inertialTangent);

    const Result result = RStructuralCableDynamicsSolver::step(
        nodes,members,{0.0,mass,0.0},displacement,velocity,acceleration,
        prescribed,prescribedValues,loads,options);
    QVERIFY(std::abs(result.displacements[4]-(-sag+expectedPerturbation)) < 2.0e-9);
    QVERIFY(result.iterations > 0);
    QVERIFY(result.residualNorm < 1.0e-10);
    const std::vector<double> freshForces = internalForces(nodes,members,result.displacements);
    const double freshResidual = freeResidualNorm(
        freshForces,result.accelerations,{0.0,mass,0.0},loads,prescribed);
    QVERIFY(std::abs(freshResidual-result.residualNorm) < 1.0e-12);
    QVERIFY(std::abs(result.reactions[4]) <= 1.0e-10);
    QVERIFY(std::abs(result.accelerations[4]-
        (result.displacements[4]-(-sag+perturbation)-dt*velocity[4]-
         dt*dt*(0.5-options.beta)*acceleration[4])/(options.beta*dt*dt)) < 1.0e-10);
}

void TestStructuralCableDynamics::prescribedMotionKeepsExactKinematicsAndSlackCanRetension()
{
    const std::vector<Node> nodes = {{{0.0,0.0,0.0}},{{1.0,0.0,0.0}}};
    const std::vector<Member> members = {{0,1,1.0,0.0,100.0}};
    const std::vector<double> masses = {0.0,1.0};
    const std::vector<bool> prescribed = fixedMask(nodes.size());
    std::vector<double> displacement(6,0.0),velocity(6,0.0),acceleration(6,0.0);
    std::vector<double> prescribedValues(6,0.0);
    prescribedValues[3] = -0.2;
    const std::vector<double> loads(6,0.0);
    const auto options = optionsFor(0.1);
    const Result slack = RStructuralCableDynamicsSolver::step(
        nodes,members,masses,displacement,velocity,acceleration,prescribed,
        prescribedValues,loads,options);
    QCOMPARE(slack.displacements[3],-0.2);
    QCOMPARE(slack.elementTensions[0],0.0);
    QCOMPARE(slack.residualNorm,0.0);
    QVERIFY(std::abs(slack.velocities[3]-options.gamma*options.dt*slack.accelerations[3]) < 1.0e-12);
    QVERIFY(std::abs(slack.supportWorkIncrement-8.0) < 1.0e-12);
    QVERIFY(std::abs(slack.energyBalanceError) < 1.0e-12);

    prescribedValues[3] = 0.2;
    const Result retensioned = RStructuralCableDynamicsSolver::step(
        nodes,members,masses,slack.displacements,slack.velocities,slack.accelerations,
        prescribed,prescribedValues,loads,options);
    QCOMPARE(retensioned.displacements[3],0.2);
    QVERIFY(std::abs(retensioned.elementTensions[0]-20.0) < 1.0e-12);
    QCOMPARE(retensioned.residualNorm,0.0);
    for (std::size_t dof=0;dof<6;dof++)
        QCOMPARE(retensioned.displacements[dof],prescribedValues[dof]);
}

void TestStructuralCableDynamics::failedStepLeavesAcceptedStateUnchanged()
{
    const std::vector<Node> nodes = {{{0.0,0.0,0.0}},{{1.0,0.0,0.0}}};
    const std::vector<Member> members = {{0,1,1.0,0.0,100.0}};
    const std::vector<double> masses = {0.0,0.0};
    const std::vector<double> displacement(6,0.0),velocity(6,0.0),acceleration(6,0.0);
    const std::vector<bool> prescribed = {true,true,true,false,true,true};
    const std::vector<double> prescribedValues(6,0.0);
    std::vector<double> loads(6,0.0);
    loads[3] = 1.0;
    const auto displacementBefore = displacement;
    const auto velocityBefore = velocity;
    const auto accelerationBefore = acceleration;
    QVERIFY_EXCEPTION_THROWN(RStructuralCableDynamicsSolver::step(
        nodes,members,masses,displacement,velocity,acceleration,prescribed,
        prescribedValues,loads,optionsFor(0.1)),std::runtime_error);
    QVERIFY(displacement == displacementBefore);
    QVERIFY(velocity == velocityBefore);
    QVERIFY(acceleration == accelerationBefore);

    auto unstable = optionsFor(0.1);
    unstable.beta = 0.24;
    QVERIFY_EXCEPTION_THROWN(RStructuralCableDynamicsSolver::step(
        nodes,members,masses,displacement,velocity,acceleration,prescribed,
        prescribedValues,loads,unstable),std::invalid_argument);
}

void TestStructuralCableDynamics::timeStepRefinementConvergesForAxialOscillator()
{
    const double exact = 1.0e-3*std::cos(10.0*0.2);
    const double coarseError = std::abs(integrateLinearOscillator(0.02)-exact);
    const double mediumError = std::abs(integrateLinearOscillator(0.01)-exact);
    const double fineError = std::abs(integrateLinearOscillator(0.005)-exact);
    QVERIFY(coarseError > mediumError);
    QVERIFY(mediumError > fineError);
    QVERIFY(fineError < 0.35*mediumError);
}

void TestStructuralCableDynamics::unforcedSymmetricCableOscillationBalancesEnergy()
{
    const std::vector<Node> nodes = {
        {{-1.0,0.0,0.0}},{{0.0,0.0,0.0}},{{1.0,0.0,0.0}}
    };
    const std::vector<Member> members = {
        {0,1,0.8,10.0,100.0},{1,2,0.8,10.0,100.0}
    };
    const std::vector<double> masses = {0.0,1.0,0.0};
    std::vector<double> displacement(9,0.0),velocity(9,0.0),acceleration(9,0.0);
    velocity[3] = 0.01;
    std::vector<bool> prescribed = fixedMask(nodes.size());
    prescribed[3] = false;
    const std::vector<double> prescribedValues(9,0.0);
    const std::vector<double> loads(9,0.0);
    const auto options = optionsFor(0.001);
    double maximumEnergyBalanceError = 0.0;
    double initialMechanicalEnergy = 0.0;
    for (std::size_t step=0;step<20;step++)
    {
        const Result result = RStructuralCableDynamicsSolver::step(
            nodes,members,masses,displacement,velocity,acceleration,prescribed,
            prescribedValues,loads,options);
        maximumEnergyBalanceError = std::max(
            maximumEnergyBalanceError,std::abs(result.energyBalanceError));
        if (step == 0)
            initialMechanicalEnergy = result.previousStrainEnergy+result.previousKineticEnergy;
        displacement = result.displacements;
        velocity = result.velocities;
        acceleration = result.accelerations;
    }
    QVERIFY(initialMechanicalEnergy > 0.0);
    QVERIFY(maximumEnergyBalanceError < 1.0e-9*initialMechanicalEnergy);
    const std::vector<double> finalForces = internalForces(nodes,members,displacement);
    QVERIFY(freeResidualNorm(finalForces,acceleration,masses,loads,prescribed) < 1.0e-10);
}

QTEST_APPLESS_MAIN(TestStructuralCableDynamics)

#include "tst_structural_cable_dynamics.moc"

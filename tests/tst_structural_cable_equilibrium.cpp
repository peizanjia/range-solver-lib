#include <QtTest>

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

#include "rstructural_cable_equilibrium.h"

namespace
{
using Node = RStructuralCableVector3;
using Member = RStructuralCableEquilibriumMember;

std::vector<bool> fixedMask(std::size_t nodeCount)
{
    return std::vector<bool>(3*nodeCount,true);
}

std::vector<double> recomputeInternalForces(
    const std::vector<Node> &nodes,
    const std::vector<Member> &members,
    const std::vector<double> &displacements)
{
    std::vector<double> forces(displacements.size(),0.0);
    for (const auto &member : members)
    {
        Node node1{}, node2{};
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

double freshFreeResidualNorm(const std::vector<double> &internalForces,
                             const std::vector<double> &loads,
                             const std::vector<bool> &prescribedMask)
{
    double norm = 0.0;
    for (std::size_t dof=0;dof<internalForces.size();dof++)
    {
        if (!prescribedMask[dof]) norm = std::hypot(norm,internalForces[dof]-loads[dof]);
    }
    return norm;
}
}

class TestStructuralCableEquilibrium : public QObject
{
    Q_OBJECT

private slots:
    void axialCableBalancesKnownLoad();
    void smallForceScaleStillSolvesTheLinearCorrection();
    void manufacturedSkewCableEquilibriumUsesAllFreeComponents();
    void symmetricTwoSegmentNetworkRecoversKnownSag();
    void highPrestressSagSurvivesPotentialRoundoff();
    void initiallyUnbalancedPrestressRelaxesToSlackCutoff();
    void allPrescribedSlackAndRetensionStates();
    void rejectUnderiterationInvalidInputAndLoadedSlackMechanism();
};

void TestStructuralCableEquilibrium::axialCableBalancesKnownLoad()
{
    const std::vector<Node> nodes = {{{0.0,0.0,0.0}},{{1.0,0.0,0.0}}};
    const std::vector<Member> members = {{0,1,1.0,5.0,100.0}};
    const std::vector<bool> prescribed = {true,true,true,false,true,true};
    const std::vector<double> prescribedValues(6,0.0);
    std::vector<double> loads(6,0.0);
    loads[3] = 7.0;

    const auto result = RStructuralCableEquilibriumSolver::solve(
        nodes,members,prescribed,prescribedValues,loads);
    QVERIFY(std::abs(result.displacements[3]-0.02) < 1.0e-11);
    QVERIFY(std::abs(result.elementTensions[0]-7.0) < 1.0e-10);
    QVERIFY(result.residualNorm <= 1.0e-8);
    QVERIFY(std::abs(result.internalForces[3]-7.0) < 1.0e-10);
    QVERIFY(std::abs(result.reactions[0]+7.0) < 1.0e-10);
    QVERIFY(std::abs(result.reactions[3]) < 1.0e-10);
    QVERIFY(std::isfinite(result.energy));
    const auto freshInternal = recomputeInternalForces(nodes,members,result.displacements);
    const double freshResidual = freshFreeResidualNorm(freshInternal,loads,prescribed);
    QVERIFY(std::abs(result.residualNorm-freshResidual) < 1.0e-12);
    for (std::size_t dof=0;dof<freshInternal.size();dof++)
        QVERIFY(std::abs(result.internalForces[dof]-freshInternal[dof]) < 1.0e-12);
}

void TestStructuralCableEquilibrium::smallForceScaleStillSolvesTheLinearCorrection()
{
    const std::vector<Node> nodes = {{{0.0,0.0,0.0}},{{1.0,0.0,0.0}}};
    const std::vector<Member> members = {{0,1,1.0,1.0e-15,1.0e-6}};
    const std::vector<bool> prescribed = {true,true,true,false,true,true};
    const std::vector<double> prescribedValues(6,0.0);
    std::vector<double> loads(6,0.0);
    loads[3] = 2.0e-15;
    RStructuralCableEquilibriumOptions options;
    options.absForceTol = 1.0e-20;
    options.relForceTol = 0.0;

    const auto result = RStructuralCableEquilibriumSolver::solve(
        nodes,members,prescribed,prescribedValues,loads,options);
    QVERIFY(std::abs(result.displacements[3]-1.0e-9) < 2.0e-16);
    QVERIFY(std::abs(result.elementTensions[0]-2.0e-15) < 1.0e-22);
    QVERIFY(result.residualNorm <= options.absForceTol);
}

void TestStructuralCableEquilibrium::manufacturedSkewCableEquilibriumUsesAllFreeComponents()
{
    constexpr double referenceLength = 1.0;
    constexpr double initialTension = 10.0;
    constexpr double axialRigidity = 100.0;
    const Node targetDisplacement{{0.1,0.2,-0.15}};
    const Node targetRelativePosition{{1.0+targetDisplacement[0],
                                       targetDisplacement[1],
                                       targetDisplacement[2]}};
    const double targetLength = std::hypot(
        std::hypot(targetRelativePosition[0],targetRelativePosition[1]),
        targetRelativePosition[2]);
    const double targetTension = initialTension+
        axialRigidity/referenceLength*(targetLength-referenceLength);
    const std::vector<Node> nodes = {{{0.0,0.0,0.0}},{{1.0,0.0,0.0}}};
    const std::vector<Member> members = {
        {0,1,referenceLength,initialTension,axialRigidity}
    };
    const std::vector<bool> prescribed = {true,true,true,false,false,false};
    const std::vector<double> prescribedValues(6,0.0);
    std::vector<double> loads(6,0.0);
    for (std::size_t axis=0;axis<3;axis++)
        loads[3+axis] = targetTension*targetRelativePosition[axis]/targetLength;

    const auto result = RStructuralCableEquilibriumSolver::solve(
        nodes,members,prescribed,prescribedValues,loads);
    for (std::size_t axis=0;axis<3;axis++)
        QVERIFY(std::abs(result.displacements[3+axis]-targetDisplacement[axis]) < 2.0e-8);
    QVERIFY(result.iterations > 1);
    const auto freshInternal = recomputeInternalForces(nodes,members,result.displacements);
    const double freshResidual = freshFreeResidualNorm(freshInternal,loads,prescribed);
    QVERIFY(std::abs(result.residualNorm-freshResidual) < 1.0e-12);
}

void TestStructuralCableEquilibrium::symmetricTwoSegmentNetworkRecoversKnownSag()
{
    constexpr double halfSpan = 1.0;
    constexpr double sag = 0.4;
    constexpr double referenceLength = 1.0;
    constexpr double initialTension = 10.0;
    constexpr double axialRigidity = 100.0;
    const double targetLength = std::hypot(halfSpan,sag);
    const double targetTension = initialTension+
        (axialRigidity/referenceLength)*(targetLength-referenceLength);
    const double balancingDeadLoad = -2.0*targetTension*sag/targetLength;
    const std::vector<Node> nodes = {
        {{-halfSpan,0.0,0.0}},{{0.0,0.0,0.0}},{{halfSpan,0.0,0.0}}
    };
    const std::vector<Member> members = {
        {0,1,referenceLength,initialTension,axialRigidity},
        {1,2,referenceLength,initialTension,axialRigidity}
    };
    std::vector<bool> prescribed = fixedMask(nodes.size());
    prescribed[4] = false; // Only the center-node transverse coordinate is free.
    const std::vector<double> prescribedValues(9,0.0);
    std::vector<double> loads(9,0.0);
    loads[4] = balancingDeadLoad;

    const auto result = RStructuralCableEquilibriumSolver::solve(
        nodes,members,prescribed,prescribedValues,loads);
    QVERIFY(std::abs(result.displacements[4]+sag) < 2.0e-9);
    QVERIFY(std::abs(result.elementTensions[0]-targetTension) < 2.0e-7);
    QVERIFY(std::abs(result.elementTensions[1]-targetTension) < 2.0e-7);
    QVERIFY(result.residualNorm <= 1.0e-8+1.0e-8*std::abs(balancingDeadLoad));
    QVERIFY(result.iterations > 1);
    const auto freshInternal = recomputeInternalForces(nodes,members,result.displacements);
    const double freshResidual = freshFreeResidualNorm(freshInternal,loads,prescribed);
    QVERIFY(std::abs(result.residualNorm-freshResidual) < 1.0e-12);
}

void TestStructuralCableEquilibrium::highPrestressSagSurvivesPotentialRoundoff()
{
    // This fixture reaches a physical residual just above the outer tolerance
    // while the final Newton energy decrease (~1e-14 J) is below rounding in
    // the roughly 1.24e4 J total potential. The bounded energy-roundoff rule
    // may accept that step only because the independently computed force
    // residual drops materially; the outer certificate remains unchanged.
    const std::vector<Node> nodes = {
        {{-1.0,0.0,0.0}},{{0.0,0.0,0.0}},{{1.0,0.0,0.0}}
    };
    const std::vector<Member> members = {
        {0,1,1.25,100.0,1.0},
        {1,2,1.25,100.0,1.0}
    };
    std::vector<bool> prescribed = fixedMask(nodes.size());
    prescribed[4] = false; // Only the center-node transverse coordinate is free.
    const std::vector<double> prescribedValues(9,0.0);
    std::vector<double> loads(9,0.0);
    loads[4] = -120.0;

    const auto result = RStructuralCableEquilibriumSolver::solve(
        nodes,members,prescribed,prescribedValues,loads);
    QVERIFY(std::abs(result.displacements[4]+0.75) < 2.0e-8);
    QVERIFY(std::abs(result.elementTensions[0]-100.0) < 2.0e-6);
    QVERIFY(std::abs(result.elementTensions[1]-100.0) < 2.0e-6);
    QVERIFY(result.residualNorm <= 1.0e-8+1.0e-8*120.0);
    const auto freshInternal = recomputeInternalForces(nodes,members,result.displacements);
    const double freshResidual = freshFreeResidualNorm(freshInternal,loads,prescribed);
    QVERIFY(std::abs(result.residualNorm-freshResidual) < 1.0e-12);
}

void TestStructuralCableEquilibrium::initiallyUnbalancedPrestressRelaxesToSlackCutoff()
{
    const std::vector<Node> nodes = {{{0.0,0.0,0.0}},{{1.0,0.0,0.0}}};
    const std::vector<Member> members = {{0,1,1.0,10.0,100.0}};
    const std::vector<bool> prescribed = {true,true,true,false,false,false};
    const std::vector<double> prescribedValues(6,0.0);
    const std::vector<double> loads(6,0.0);
    const auto result = RStructuralCableEquilibriumSolver::solve(
        nodes,members,prescribed,prescribedValues,loads);
    QVERIFY(std::abs(result.displacements[3]+0.1) < 1.0e-11);
    QVERIFY(result.elementTensions[0] <= 1.0e-8);
    QVERIFY(result.residualNorm <= 1.0e-8);
    QVERIFY(result.energy >= 0.0);
    QVERIFY(result.iterations > 0);
}

void TestStructuralCableEquilibrium::allPrescribedSlackAndRetensionStates()
{
    const std::vector<Node> nodes = {{{0.0,0.0,0.0}},{{1.0,0.0,0.0}}};
    const auto prescribed = fixedMask(nodes.size());
    std::vector<double> values(6,0.0);
    values[3] = -0.5;
    const std::vector<double> noLoads(6,0.0);
    const auto slack = RStructuralCableEquilibriumSolver::solve(
        nodes,{{0,1,1.0,0.0,100.0}},prescribed,values,noLoads);
    QCOMPARE(slack.iterations,std::size_t(0));
    QCOMPARE(slack.residualNorm,0.0); // No free DOFs; this is not a stability claim.
    QCOMPARE(slack.elementTensions[0],0.0);
    for (double force : slack.internalForces) QCOMPARE(force,0.0);

    values[3] = 0.2;
    std::vector<double> loaded(6,0.0);
    loaded[3] = 3.0;
    const auto retensioned = RStructuralCableEquilibriumSolver::solve(
        nodes,{{0,1,1.0,5.0,100.0}},prescribed,values,loaded);
    QCOMPARE(retensioned.iterations,std::size_t(0));
    QCOMPARE(retensioned.elementTensions[0],25.0);
    QCOMPARE(retensioned.internalForces[0],-25.0);
    QCOMPARE(retensioned.internalForces[3],25.0);
    QCOMPARE(retensioned.reactions[3],22.0);
    QVERIFY(std::abs(retensioned.energy-2.525) < 1.0e-14);
}

void TestStructuralCableEquilibrium::rejectUnderiterationInvalidInputAndLoadedSlackMechanism()
{
    const std::vector<Node> sagNodes = {
        {{-1.0,0.0,0.0}},{{0.0,0.0,0.0}},{{1.0,0.0,0.0}}
    };
    const std::vector<Member> sagMembers = {
        {0,1,1.0,10.0,100.0},{1,2,1.0,10.0,100.0}
    };
    std::vector<bool> sagMask = fixedMask(sagNodes.size());
    sagMask[4] = false;
    std::vector<double> sagLoads(9,0.0);
    const double length = std::hypot(1.0,0.4);
    const double tension = 10.0+100.0*(length-1.0);
    sagLoads[4] = -2.0*tension*0.4/length;
    RStructuralCableEquilibriumOptions oneIteration;
    oneIteration.maxIterations = 1;
    QVERIFY_EXCEPTION_THROWN(RStructuralCableEquilibriumSolver::solve(
        sagNodes,sagMembers,sagMask,std::vector<double>(9,0.0),sagLoads,oneIteration),
        std::runtime_error);

    const std::vector<Node> twoNodes = {{{0.0,0.0,0.0}},{{1.0,0.0,0.0}}};
    const std::vector<Member> oneMember = {{0,1,1.0,0.0,100.0}};
    const auto fixed = fixedMask(twoNodes.size());
    std::vector<double> values(6,0.0);
    QVERIFY_EXCEPTION_THROWN(RStructuralCableEquilibriumSolver::solve(
        twoNodes,oneMember,fixed,std::vector<double>(5,0.0),std::vector<double>(6,0.0)),
        std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralCableEquilibriumSolver::solve(
        twoNodes,{{0,0,1.0,0.0,100.0}},fixed,values,std::vector<double>(6,0.0)),
        std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralCableEquilibriumSolver::solve(
        twoNodes,{{0,1,0.0,0.0,100.0}},fixed,values,std::vector<double>(6,0.0)),
        std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralCableEquilibriumSolver::solve(
        std::vector<Node>{{{0.0,0.0,0.0}},{{0.0,0.0,0.0}}},oneMember,fixed,values,
        std::vector<double>(6,0.0)),std::invalid_argument);
    std::vector<double> nonfiniteLoads(6,0.0);
    nonfiniteLoads[3] = std::numeric_limits<double>::infinity();
    QVERIFY_EXCEPTION_THROWN(RStructuralCableEquilibriumSolver::solve(
        twoNodes,oneMember,fixed,values,nonfiniteLoads),std::invalid_argument);

    // At the tension cutoff, the slack tangent is zero. A transverse or axial
    // load on the free endpoint therefore exposes a singular loaded mechanism.
    const std::vector<bool> oneFixed = {true,true,true,false,false,false};
    std::vector<double> loadedSlack(6,0.0);
    loadedSlack[3] = 1.0;
    QVERIFY_EXCEPTION_THROWN(RStructuralCableEquilibriumSolver::solve(
        twoNodes,oneMember,oneFixed,std::vector<double>(6,0.0),loadedSlack),
        std::runtime_error);
}

QTEST_APPLESS_MAIN(TestStructuralCableEquilibrium)

#include "tst_structural_cable_equilibrium.moc"

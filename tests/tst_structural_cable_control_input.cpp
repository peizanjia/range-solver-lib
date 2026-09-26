#include <QtTest>

#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

#include "rstructural_cable_control_input.h"

namespace
{
using Vector3 = RStructuralCableVector3;
using Vector6 = RStructuralCableVector6;

RStructuralCableKernelResult evaluate(const Vector6 &position,
                                      double referenceLength,
                                      double initialTension,
                                      double axialRigidity)
{
    return RStructuralCableKernel::evaluate(
        {{position[0],position[1],position[2]}},
        {{position[3],position[4],position[5]}},
        referenceLength,initialTension,axialRigidity);
}

Vector6 orientationOne()
{
    return {{0.2,-0.3,0.6,1.4,0.8,1.7}};
}

Vector6 orientationTwo()
{
    return {{-1.1,0.3,-0.8,-0.2,1.5,0.4}};
}

Vector6 orientationThree()
{
    return {{0.4,1.6,-0.2,1.7,0.1,1.1}};
}

bool close(double actual, double expected, double tolerance = 3.0e-6)
{
    return std::abs(actual-expected) <= tolerance*(1.0+std::abs(expected));
}
}

class TestStructuralCableControlInput : public QObject
{
    Q_OBJECT

private slots:
    void scalarInputsMatchForceAndEnergyFiniteDifferences();
    void anchorTranslationMatchesForceAndVirtualWorkDerivatives();
    void preloadedAnchorMapIncludesGeometricStiffness();
    void slackIsZeroAndKinkIsRejected();
    void rejectInvalidParametersAndEndpoint();
};

void TestStructuralCableControlInput::scalarInputsMatchForceAndEnergyFiniteDifferences()
{
    const std::array<Vector6,3> orientations{{orientationOne(),orientationTwo(),orientationThree()}};
    constexpr double referenceLength = 1.2;
    constexpr double initialTension = 1.5;
    constexpr double axialRigidity = 42.0;
    constexpr double hLength = 1.0e-6;
    constexpr double hTension = 1.0e-5;

    for (const Vector6 &position : orientations)
    {
        const auto reference = RStructuralCableControlInput::referenceLength(
            {{position[0],position[1],position[2]}},
            {{position[3],position[4],position[5]}},
            referenceLength,initialTension,axialRigidity);
        const auto prestress = RStructuralCableControlInput::initialTension(
            {{position[0],position[1],position[2]}},
            {{position[3],position[4],position[5]}},
            referenceLength,initialTension,axialRigidity);
        QVERIFY(evaluate(position,referenceLength,initialTension,axialRigidity).taut);

        const auto lengthPlus = evaluate(position,referenceLength+hLength,
                                         initialTension,axialRigidity);
        const auto lengthMinus = evaluate(position,referenceLength-hLength,
                                          initialTension,axialRigidity);
        const auto tensionPlus = evaluate(position,referenceLength,
                                          initialTension+hTension,axialRigidity);
        const auto tensionMinus = evaluate(position,referenceLength,
                                           initialTension-hTension,axialRigidity);

        for (std::size_t row=0;row<6;row++)
        {
            const double lengthForceMap = -
                (lengthPlus.internalForce[row]-lengthMinus.internalForce[row])/(2.0*hLength);
            const double tensionForceMap = -
                (tensionPlus.internalForce[row]-tensionMinus.internalForce[row])/(2.0*hTension);
            QVERIFY2(close(reference.forceMap[row],lengthForceMap),
                     "Reference-length force map differs from an independent finite difference");
            QVERIFY2(close(prestress.forceMap[row],tensionForceMap),
                     "Prestress force map differs from an independent finite difference");
        }

        const double lengthWork = (lengthPlus.energy-lengthMinus.energy)/(2.0*hLength);
        const double tensionWork = (tensionPlus.energy-tensionMinus.energy)/(2.0*hTension);
        QVERIFY2(close(reference.potentialDerivative,lengthWork),
                 "Reference-length potential derivative violates finite-difference virtual work");
        QVERIFY2(close(prestress.potentialDerivative,tensionWork),
                 "Prestress potential derivative violates finite-difference virtual work");
    }
}

void TestStructuralCableControlInput::anchorTranslationMatchesForceAndVirtualWorkDerivatives()
{
    const std::array<Vector6,3> orientations{{orientationOne(),orientationTwo(),orientationThree()}};
    constexpr double referenceLength = 1.2;
    constexpr double initialTension = 1.5;
    constexpr double axialRigidity = 42.0;
    constexpr double h = 1.0e-6;

    for (const Vector6 &base : orientations)
    {
        QVERIFY(evaluate(base,referenceLength,initialTension,axialRigidity).taut);
        for (const auto endpoint : {RStructuralCableControlInputEndpoint::Node1,
                                    RStructuralCableControlInputEndpoint::Node2})
        {
            const auto linearization = RStructuralCableControlInput::anchorTranslation(
                {{base[0],base[1],base[2]}},
                {{base[3],base[4],base[5]}},
                referenceLength,initialTension,axialRigidity,endpoint);
            const std::size_t offset = endpoint == RStructuralCableControlInputEndpoint::Node1 ? 0 : 3;
            for (std::size_t coordinate=0;coordinate<3;coordinate++)
            {
                Vector6 plus = base;
                Vector6 minus = base;
                plus[offset+coordinate] += h;
                minus[offset+coordinate] -= h;
                const auto forcePlus = evaluate(plus,referenceLength,initialTension,axialRigidity);
                const auto forceMinus = evaluate(minus,referenceLength,initialTension,axialRigidity);
                for (std::size_t row=0;row<6;row++)
                {
                    const double forceMap = -(forcePlus.internalForce[row]-
                                              forceMinus.internalForce[row])/(2.0*h);
                    QVERIFY2(close(linearization.forceMap[row][coordinate],forceMap),
                             "Anchor force map differs from an independent force finite difference");
                }

                const double potentialDerivative = (forcePlus.energy-forceMinus.energy)/(2.0*h);
                QVERIFY2(close(linearization.potentialDerivative[coordinate],potentialDerivative),
                         "Anchor potential derivative violates finite-difference virtual work");
            }
        }
    }
}

void TestStructuralCableControlInput::preloadedAnchorMapIncludesGeometricStiffness()
{
    const Vector3 node1{{-0.2,0.3,0.4}};
    const Vector3 node2{{1.4,1.5,1.6}};
    const Vector3 difference{{node2[0]-node1[0],node2[1]-node1[1],node2[2]-node1[2]}};
    const double length = std::hypot(std::hypot(difference[0],difference[1]),difference[2]);
    constexpr double preload = 6.0;
    constexpr double axialRigidity = 100.0;
    const auto response = RStructuralCableKernel::evaluate(
        node1,node2,length,preload,axialRigidity);
    QVERIFY(response.taut);
    QCOMPARE(response.tension,preload);

    const auto linearization = RStructuralCableControlInput::anchorTranslation(
        node1,node2,length,preload,axialRigidity,
        RStructuralCableControlInputEndpoint::Node2);
    const double directionY = difference[1]/length;
    const double directionZ = difference[2]/length;
    const double expected = -(axialRigidity/length-preload/length)*directionY*directionZ;
    QVERIFY(close(linearization.forceMap[4][2],expected,1.0e-12));
    QVERIFY(std::abs(linearization.forceMap[4][2]) > 1.0e-3);
}

void TestStructuralCableControlInput::slackIsZeroAndKinkIsRejected()
{
    const Vector3 node1{{0.0,0.0,0.0}};
    const Vector3 slackNode2{{0.4,-0.1,0.2}};
    const auto slackLength = RStructuralCableControlInput::referenceLength(
        node1,slackNode2,1.0,0.2,4.0);
    const auto slackTension = RStructuralCableControlInput::initialTension(
        node1,slackNode2,1.0,0.2,4.0);
    const auto slackAnchor = RStructuralCableControlInput::anchorTranslation(
        node1,slackNode2,1.0,0.2,4.0,RStructuralCableControlInputEndpoint::Node2);
    for (double value : slackLength.forceMap) QCOMPARE(value,0.0);
    for (double value : slackTension.forceMap) QCOMPARE(value,0.0);
    QCOMPARE(slackLength.potentialDerivative,0.0);
    QCOMPARE(slackTension.potentialDerivative,0.0);
    for (const auto &row : slackAnchor.forceMap)
        for (double value : row) QCOMPARE(value,0.0);
    for (double value : slackAnchor.potentialDerivative) QCOMPARE(value,0.0);

    const Vector3 kinkNode2{{0.5,0.0,0.0}};
    QVERIFY_EXCEPTION_THROWN(
        RStructuralCableControlInput::referenceLength(node1,kinkNode2,1.0,2.0,4.0),
        std::domain_error);
    QVERIFY_EXCEPTION_THROWN(
        RStructuralCableControlInput::initialTension(node1,kinkNode2,1.0,2.0,4.0),
        std::domain_error);
    QVERIFY_EXCEPTION_THROWN(
        RStructuralCableControlInput::anchorTranslation(
            node1,kinkNode2,1.0,2.0,4.0,RStructuralCableControlInputEndpoint::Node1),
        std::domain_error);
}

void TestStructuralCableControlInput::rejectInvalidParametersAndEndpoint()
{
    const Vector3 node1{{0.0,0.0,0.0}};
    const Vector3 node2{{1.5,0.4,0.2}};
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double infinity = std::numeric_limits<double>::infinity();

    QVERIFY_EXCEPTION_THROWN(
        RStructuralCableControlInput::referenceLength(node1,node2,nan,1.0,20.0),
        std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(
        RStructuralCableControlInput::initialTension(node1,node2,1.0,infinity,20.0),
        std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(
        RStructuralCableControlInput::anchorTranslation(
            node1,node2,1.0,1.0,0.0,RStructuralCableControlInputEndpoint::Node1),
        std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(
        RStructuralCableControlInput::anchorTranslation(
            node1,node2,1.0,1.0,20.0,static_cast<RStructuralCableControlInputEndpoint>(99)),
        std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(
        RStructuralCableControlInput::referenceLength(
            {{nan,0.0,0.0}},node2,1.0,1.0,20.0),
        std::invalid_argument);
}

QTEST_APPLESS_MAIN(TestStructuralCableControlInput)

#include "tst_structural_cable_control_input.moc"

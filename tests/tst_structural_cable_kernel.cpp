#include <QtTest>

#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>

#include "rstructural_cable_kernel.h"

namespace
{
using Vector3 = RStructuralCableVector3;
using Vector6 = RStructuralCableVector6;
using Matrix3 = std::array<std::array<double,3>,3>;

RStructuralCableKernelResult evaluate(const Vector6 &position,
                                      double referenceLength = 1.2,
                                      double initialTension = 3.0,
                                      double axialRigidity = 80.0)
{
    return RStructuralCableKernel::evaluate(
        {{position[0],position[1],position[2]}},
        {{position[3],position[4],position[5]}},
        referenceLength,initialTension,axialRigidity);
}

Vector6 skewStretchedPosition()
{
    return {{0.1,-0.2,0.4,1.25,0.73,1.05}};
}

Vector6 rotate(const Vector6 &position, const Matrix3 &rotation)
{
    Vector6 result{};
    for (std::size_t node=0;node<2;node++)
    {
        for (std::size_t i=0;i<3;i++)
        {
            for (std::size_t j=0;j<3;j++)
                result[node*3+i] += rotation[i][j]*position[node*3+j];
        }
    }
    return result;
}

Vector6 rotateVector(const Vector6 &vector, const Matrix3 &rotation)
{
    return rotate(vector,rotation);
}
}

class TestStructuralCableKernel : public QObject
{
    Q_OBJECT

private slots:
    void consistentTangentMatchesAllSixForceDerivatives();
    void internalForceIsEnergyGradientAndTangentIsSymmetric();
    void rigidTranslationAndRotationObjectivity();
    void slackKinkAndRetension();
    void balancedStraightSegmentsHaveAnalyticTransverseStiffness();
    void rejectInvalidAndDegenerateInputs();
};

void TestStructuralCableKernel::consistentTangentMatchesAllSixForceDerivatives()
{
    const Vector6 position = skewStretchedPosition();
    const auto base = evaluate(position);
    QVERIFY(base.taut);
    const double h = 1.0e-6;
    for (std::size_t column=0;column<6;column++)
    {
        Vector6 plus = position;
        Vector6 minus = position;
        plus[column] += h;
        minus[column] -= h;
        const auto forcePlus = evaluate(plus).internalForce;
        const auto forceMinus = evaluate(minus).internalForce;
        for (std::size_t row=0;row<6;row++)
        {
            const double finiteDifference = (forcePlus[row]-forceMinus[row])/(2.0*h);
            const double expected = base.tangent[row][column];
            QVERIFY2(std::abs(finiteDifference-expected) <= 2.0e-6*(1.0+std::abs(expected)),
                     "Consistent cable tangent differs from a force finite difference");
        }
    }
}

void TestStructuralCableKernel::internalForceIsEnergyGradientAndTangentIsSymmetric()
{
    const Vector6 position = skewStretchedPosition();
    const auto base = evaluate(position);
    const double h = 1.0e-6;
    for (std::size_t coordinate=0;coordinate<6;coordinate++)
    {
        Vector6 plus = position;
        Vector6 minus = position;
        plus[coordinate] += h;
        minus[coordinate] -= h;
        const double energyGradient = (evaluate(plus).energy-evaluate(minus).energy)/(2.0*h);
        QVERIFY2(std::abs(energyGradient-base.internalForce[coordinate]) <=
                     2.0e-7*(1.0+std::abs(base.internalForce[coordinate])),
                 "Cable internal force differs from the energy gradient");
    }
    for (std::size_t i=0;i<6;i++)
    {
        for (std::size_t j=0;j<6;j++)
            QVERIFY(std::abs(base.tangent[i][j]-base.tangent[j][i]) < 1.0e-12);
    }
}

void TestStructuralCableKernel::rigidTranslationAndRotationObjectivity()
{
    const Vector6 position = skewStretchedPosition();
    Vector6 translated = position;
    for (std::size_t node=0;node<2;node++)
    {
        translated[node*3] += 10.0;
        translated[node*3+1] -= 5.0;
        translated[node*3+2] += 7.0;
    }
    const auto original = evaluate(position);
    const auto shifted = evaluate(translated);
    QVERIFY(std::abs(original.length-shifted.length) < 1.0e-13);
    QVERIFY(std::abs(original.energy-shifted.energy) < 1.0e-12);
    for (std::size_t i=0;i<6;i++)
    {
        QVERIFY(std::abs(original.internalForce[i]-shifted.internalForce[i]) < 1.0e-12);
        for (std::size_t j=0;j<6;j++)
            QVERIFY(std::abs(original.tangent[i][j]-shifted.tangent[i][j]) < 1.0e-12);
    }

    // Proper 90-degree rotation about z.
    const Matrix3 rotation{{{{0.0,-1.0,0.0}},{{1.0,0.0,0.0}},{{0.0,0.0,1.0}}}};
    const auto rotated = evaluate(rotate(position,rotation));
    const Vector6 expectedForce = rotateVector(original.internalForce,rotation);
    for (std::size_t i=0;i<6;i++)
        QVERIFY(std::abs(rotated.internalForce[i]-expectedForce[i]) < 1.0e-12);
    for (std::size_t i=0;i<6;i++)
    {
        for (std::size_t j=0;j<6;j++)
        {
            double expected = 0.0;
            for (std::size_t a=0;a<6;a++)
            {
                for (std::size_t b=0;b<6;b++)
                {
                    const double left = rotation[i%3][a%3];
                    const double right = rotation[j%3][b%3];
                    if ((i < 3) != (a < 3) || (j < 3) != (b < 3)) continue;
                    expected += left*original.tangent[a][b]*right;
                }
            }
            QVERIFY(std::abs(rotated.tangent[i][j]-expected) < 1.0e-12);
        }
    }
}

void TestStructuralCableKernel::slackKinkAndRetension()
{
    const auto compressed = RStructuralCableKernel::evaluate(
        {{0.0,0.0,0.0}},{{0.4,0.0,0.0}},1.0,2.0,4.0);
    QVERIFY(!compressed.taut);
    QVERIFY(compressed.trialTension < 0.0);
    QCOMPARE(compressed.tension,0.0);
    QCOMPARE(compressed.energy,0.0);
    for (double value : compressed.internalForce) QCOMPARE(value,0.0);
    for (const auto &row : compressed.tangent)
        for (double value : row) QCOMPARE(value,0.0);

    // Ntrial = 2 + (4/1)*(0.5-1) = 0. The selected kink tangent is slack.
    const auto kink = RStructuralCableKernel::evaluate(
        {{0.0,0.0,0.0}},{{0.5,0.0,0.0}},1.0,2.0,4.0);
    QCOMPARE(kink.trialTension,0.0);
    QVERIFY(!kink.taut);
    QCOMPARE(kink.tension,0.0);
    QCOMPARE(kink.energy,0.0);
    for (const auto &row : kink.tangent)
        for (double value : row) QCOMPARE(value,0.0);

    const auto retensioned = RStructuralCableKernel::evaluate(
        {{0.0,0.0,0.0}},{{1.0,0.0,0.0}},1.0,2.0,4.0);
    QVERIFY(retensioned.taut);
    QCOMPARE(retensioned.tension,2.0);
    QCOMPARE(retensioned.energy,0.5); // N^2/(2k) with k=4.
    QCOMPARE(retensioned.internalForce[0],-2.0);
    QCOMPARE(retensioned.internalForce[3],2.0);
    QCOMPARE(retensioned.tangent[1][1],2.0); // N/L geometric stiffness.
}

void TestStructuralCableKernel::balancedStraightSegmentsHaveAnalyticTransverseStiffness()
{
    // Two equal taut segments meet at the origin and pull in opposite x
    // directions. Their nodal forces cancel; each contributes N/L in y.
    const Vector3 center{{0.0,0.0,0.0}};
    const auto left = RStructuralCableKernel::evaluate(
        {{-2.0,0.0,0.0}},center,2.0,5.0,100.0);
    const auto right = RStructuralCableKernel::evaluate(
        center,{{2.0,0.0,0.0}},2.0,5.0,100.0);
    QVERIFY(left.taut && right.taut);
    QVERIFY(std::abs(left.internalForce[3]+right.internalForce[0]) < 1.0e-14);
    QVERIFY(std::abs(left.internalForce[4]+right.internalForce[1]) < 1.0e-14);
    const double centralTransverseStiffness = left.tangent[4][4]+right.tangent[1][1];
    QVERIFY(std::abs(centralTransverseStiffness-2.0*5.0/2.0) < 1.0e-14);
}

void TestStructuralCableKernel::rejectInvalidAndDegenerateInputs()
{
    const Vector3 origin{{0.0,0.0,0.0}};
    const Vector3 end{{1.0,0.0,0.0}};
    QVERIFY_EXCEPTION_THROWN(RStructuralCableKernel::evaluate(origin,end,0.0,0.0,1.0),
                             std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralCableKernel::evaluate(origin,end,1.0,-1.0,1.0),
                             std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralCableKernel::evaluate(origin,end,1.0,0.0,0.0),
                             std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralCableKernel::evaluate(origin,origin,1.0,0.0,1.0),
                             std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralCableKernel::evaluate(
        {{std::numeric_limits<double>::quiet_NaN(),0.0,0.0}},end,1.0,0.0,1.0),
        std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralCableKernel::evaluate(
        origin,{{std::numeric_limits<double>::infinity(),0.0,0.0}},1.0,0.0,1.0),
        std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralCableKernel::evaluate(
        origin,end,std::numeric_limits<double>::min(),0.0,
        std::numeric_limits<double>::max()),std::invalid_argument);
    QVERIFY_EXCEPTION_THROWN(RStructuralCableKernel::evaluate(
        origin,end,1.0,std::numeric_limits<double>::max(),
        std::numeric_limits<double>::max()),std::invalid_argument);
}

QTEST_APPLESS_MAIN(TestStructuralCableKernel)

#include "tst_structural_cable_kernel.moc"

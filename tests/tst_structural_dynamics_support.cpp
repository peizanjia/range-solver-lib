#include <QtTest>

#include <cmath>
#include <stdexcept>
#include <vector>

#include "rstructural_dynamics_support.h"

namespace
{
RSparseMatrix sparse(uint rows)
{
    RSparseMatrix matrix;
    matrix.setNRows(rows);
    return matrix;
}

RSparseMatrix scalar(double value)
{
    RSparseMatrix matrix = sparse(1);
    if (value != 0.0)
    {
        matrix.addValue(0,0,value);
    }
    return matrix;
}

RStructuralDynamicsConstraintMotion zeroParticular(uint count)
{
    return {RRVector(count,0.0),RRVector(count,0.0),RRVector(count,0.0)};
}
}

class TestStructuralDynamicsSupport : public QObject
{
    Q_OBJECT

private slots:
    void consistentMassSupportAcceleration();
    void supportVelocityDamping();
    void movingMpcParticularAndReconstruction();
    void rejectIncorrectDimensions();
};

void TestStructuralDynamicsSupport::consistentMassSupportAcceleration()
{
    // Two-node linear 1D element: M = rho A L / 6 * [2 1; 1 2].
    // Fix node c and retain node f. With initial q=v=f=0 and gddot=A,
    // M_ff a_f = -M_fc A, hence a_f = -A/2 (not zero).
    constexpr double elementMass = 6.0;
    constexpr double supportAcceleration = 4.0;
    const RSparseMatrix Mff = scalar(elementMass/3.0);
    const RSparseMatrix Mfc = scalar(elementMass/6.0);
    RRMatrix T(1,1,1.0);
    const RStructuralDynamicsSupportMotion support{RRVector(1,0.0),
                                                   RRVector(1,0.0),
                                                   RRVector(1,supportAcceleration)};
    const RRVector force = RStructuralDynamicsSupportCoupling::effectiveForce(
        scalar(10.0),sparse(1),Mff,scalar(-10.0),sparse(1),Mfc,
        T,RRVector(1,0.0),support,zeroParticular(1));
    const double freeAcceleration = force[0]/Mff.findValue(0,0);
    QVERIFY(std::abs(freeAcceleration+supportAcceleration/2.0) < 1.0e-14);

    const auto physical = RStructuralDynamicsSupportCoupling::reconstruct(
        T,RRVector(1,0.0),RRVector(1,0.0),RRVector(1,freeAcceleration),
        support,zeroParticular(1));
    QCOMPARE(physical.freeAcceleration[0],freeAcceleration);
    QCOMPARE(physical.support.acceleration[0],supportAcceleration);
}

void TestStructuralDynamicsSupport::supportVelocityDamping()
{
    // Off-diagonal damping is an independent forcing channel: -C_fc gdot.
    RRMatrix T(1,1,1.0);
    const RStructuralDynamicsSupportMotion support{RRVector(1,0.0),
                                                   RRVector(1,2.5),
                                                   RRVector(1,0.0)};
    const RRVector force = RStructuralDynamicsSupportCoupling::effectiveForce(
        sparse(1),scalar(7.0),scalar(4.0),
        sparse(1),scalar(3.0),sparse(1),
        T,RRVector(1,0.0),support,zeroParticular(1));
    QCOMPARE(force.size(),std::size_t(1));
    QCOMPARE(force[0],-7.5);
}

void TestStructuralDynamicsSupport::movingMpcParticularAndReconstruction()
{
    // q_f = [1; 2] z + p(t): each p derivative contributes to effective force.
    RRMatrix T(2,1,0.0);
    T[0][0] = 1.0;
    T[1][0] = 2.0;
    RSparseMatrix Kff = sparse(2), Cff = sparse(2), Mff = sparse(2);
    Kff.addValue(0,0,2.0); Kff.addValue(1,1,3.0);
    Cff.addValue(0,0,4.0); Cff.addValue(1,1,5.0);
    Mff.addValue(0,0,6.0); Mff.addValue(1,1,7.0);
    const auto particular = RStructuralDynamicsConstraintMotion{
        RRVector(std::vector<double>{1.0,2.0}),
        RRVector(std::vector<double>{3.0,4.0}),
        RRVector(std::vector<double>{5.0,6.0})};
    const auto support = RStructuralDynamicsSupportMotion{
        RRVector(0),RRVector(0),RRVector(0)};
    const RRVector force = RStructuralDynamicsSupportCoupling::effectiveForce(
        Kff,Cff,Mff,sparse(2),sparse(2),sparse(2),T,
        RRVector(2,0.0),support,particular);
    // Physical residual = [-44, -68]; T^T residual = -44 - 2*68.
    QCOMPARE(force[0],-180.0);

    const auto physical = RStructuralDynamicsSupportCoupling::reconstruct(
        T,RRVector(1,8.0),RRVector(1,9.0),RRVector(1,10.0),
        support,particular);
    QCOMPARE(physical.freeDisplacement[0],9.0);
    QCOMPARE(physical.freeDisplacement[1],18.0);
    QCOMPARE(physical.freeVelocity[0],12.0);
    QCOMPARE(physical.freeVelocity[1],22.0);
    QCOMPARE(physical.freeAcceleration[0],15.0);
    QCOMPARE(physical.freeAcceleration[1],26.0);
}

void TestStructuralDynamicsSupport::rejectIncorrectDimensions()
{
    RRMatrix T(1,1,1.0);
    const auto support = RStructuralDynamicsSupportMotion{
        RRVector(1,0.0),RRVector(0),RRVector(1,0.0)};
    QVERIFY_EXCEPTION_THROWN(RStructuralDynamicsSupportCoupling::effectiveForce(
        sparse(1),sparse(1),scalar(1.0),
        sparse(1),sparse(1),sparse(1),T,
        RRVector(1,0.0),support,zeroParticular(1)),std::invalid_argument);
}

QTEST_GUILESS_MAIN(TestStructuralDynamicsSupport)
#include "tst_structural_dynamics_support.moc"

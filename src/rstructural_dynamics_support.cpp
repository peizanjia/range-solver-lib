#include "rstructural_dynamics_support.h"

#include <cmath>
#include <stdexcept>

namespace
{
void checkVector(const RRVector &vector, uint expected, const char *name)
{
    if (vector.size() != expected)
    {
        throw std::invalid_argument(name);
    }
    for (double value : vector)
    {
        if (!std::isfinite(value))
        {
            throw std::invalid_argument(name);
        }
    }
}

void checkMatrix(const RSparseMatrix &matrix, uint rows, uint columns,
                 const char *name)
{
    if (matrix.getNRows() != rows)
    {
        throw std::invalid_argument(name);
    }
    for (uint row=0; row<rows; ++row)
    {
        for (uint position=0; position<matrix.getNColumns(row); ++position)
        {
            if (matrix.getVector(row).getIndex(position) >= columns ||
                !std::isfinite(matrix.getValue(row,position)))
            {
                throw std::invalid_argument(name);
            }
        }
    }
}

uint checkNullspace(const RRMatrix &nullspace)
{
    const uint rows = nullspace.getNRows();
    const uint columns = nullspace.getNColumns();
    if (rows == 0 || columns == 0)
    {
        throw std::invalid_argument("MPC nullspace must have free and reduced DOFs");
    }
    for (uint row=0; row<rows; ++row)
    {
        if (nullspace[row].size() != columns)
        {
            throw std::invalid_argument("MPC nullspace has inconsistent row lengths");
        }
        for (uint column=0; column<columns; ++column)
        {
            if (!std::isfinite(nullspace[row][column]))
            {
                throw std::invalid_argument("MPC nullspace contains a non-finite value");
            }
        }
    }
    return columns;
}

void subtractProduct(RRVector &target, const RSparseMatrix &matrix,
                     const RRVector &vector)
{
    for (uint row=0; row<matrix.getNRows(); ++row)
    {
        for (uint position=0; position<matrix.getNColumns(row); ++position)
        {
            target[row] -= matrix.getValue(row,position)*
                           vector[matrix.getVector(row).getIndex(position)];
        }
    }
}

RRVector lift(const RRMatrix &nullspace, const RRVector &reduced,
              const RRVector &particular)
{
    RRVector physical(particular);
    for (uint row=0; row<nullspace.getNRows(); ++row)
    {
        for (uint column=0; column<nullspace.getNColumns(); ++column)
        {
            physical[row] += nullspace[row][column]*reduced[column];
        }
    }
    return physical;
}
}

RRVector RStructuralDynamicsSupportCoupling::effectiveForce(
        const RSparseMatrix &freeStiffness,
        const RSparseMatrix &freeDamping,
        const RSparseMatrix &freeMass,
        const RSparseMatrix &supportStiffness,
        const RSparseMatrix &supportDamping,
        const RSparseMatrix &supportMass,
        const RRMatrix &nullspace,
        const RRVector &freeForce,
        const RStructuralDynamicsSupportMotion &support,
        const RStructuralDynamicsConstraintMotion &particular)
{
    const uint reducedCount = checkNullspace(nullspace);
    const uint freeCount = nullspace.getNRows();
    const uint supportCount = uint(support.displacement.size());
    checkVector(freeForce,freeCount,"Free force size/value mismatch");
    checkVector(support.displacement,supportCount,"Support displacement size/value mismatch");
    checkVector(support.velocity,supportCount,"Support velocity size/value mismatch");
    checkVector(support.acceleration,supportCount,"Support acceleration size/value mismatch");
    checkVector(particular.displacement,freeCount,"MPC displacement size/value mismatch");
    checkVector(particular.velocity,freeCount,"MPC velocity size/value mismatch");
    checkVector(particular.acceleration,freeCount,"MPC acceleration size/value mismatch");
    checkMatrix(freeStiffness,freeCount,freeCount,"Free stiffness size/value mismatch");
    checkMatrix(freeDamping,freeCount,freeCount,"Free damping size/value mismatch");
    checkMatrix(freeMass,freeCount,freeCount,"Free mass size/value mismatch");
    checkMatrix(supportStiffness,freeCount,supportCount,"Support stiffness size/value mismatch");
    checkMatrix(supportDamping,freeCount,supportCount,"Support damping size/value mismatch");
    checkMatrix(supportMass,freeCount,supportCount,"Support mass size/value mismatch");

    RRVector physicalForce(freeForce);
    subtractProduct(physicalForce,supportStiffness,support.displacement);
    subtractProduct(physicalForce,supportDamping,support.velocity);
    subtractProduct(physicalForce,supportMass,support.acceleration);
    subtractProduct(physicalForce,freeStiffness,particular.displacement);
    subtractProduct(physicalForce,freeDamping,particular.velocity);
    subtractProduct(physicalForce,freeMass,particular.acceleration);

    RRVector reducedForce(reducedCount,0.0);
    for (uint row=0; row<freeCount; ++row)
    {
        for (uint column=0; column<reducedCount; ++column)
        {
            reducedForce[column] += nullspace[row][column]*physicalForce[row];
        }
    }
    return reducedForce;
}

RStructuralDynamicsPhysicalState RStructuralDynamicsSupportCoupling::reconstruct(
        const RRMatrix &nullspace,
        const RRVector &reducedDisplacement,
        const RRVector &reducedVelocity,
        const RRVector &reducedAcceleration,
        const RStructuralDynamicsSupportMotion &support,
        const RStructuralDynamicsConstraintMotion &particular)
{
    const uint reducedCount = checkNullspace(nullspace);
    const uint freeCount = nullspace.getNRows();
    const uint supportCount = uint(support.displacement.size());
    checkVector(reducedDisplacement,reducedCount,"Reduced displacement size/value mismatch");
    checkVector(reducedVelocity,reducedCount,"Reduced velocity size/value mismatch");
    checkVector(reducedAcceleration,reducedCount,"Reduced acceleration size/value mismatch");
    checkVector(support.displacement,supportCount,"Support displacement size/value mismatch");
    checkVector(support.velocity,supportCount,"Support velocity size/value mismatch");
    checkVector(support.acceleration,supportCount,"Support acceleration size/value mismatch");
    checkVector(particular.displacement,freeCount,"MPC displacement size/value mismatch");
    checkVector(particular.velocity,freeCount,"MPC velocity size/value mismatch");
    checkVector(particular.acceleration,freeCount,"MPC acceleration size/value mismatch");

    RStructuralDynamicsPhysicalState physical;
    physical.freeDisplacement = lift(nullspace,reducedDisplacement,particular.displacement);
    physical.freeVelocity = lift(nullspace,reducedVelocity,particular.velocity);
    physical.freeAcceleration = lift(nullspace,reducedAcceleration,particular.acceleration);
    physical.support = support;
    return physical;
}

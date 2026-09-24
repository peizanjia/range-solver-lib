#ifndef RSTRUCTURAL_DYNAMICS_SUPPORT_H
#define RSTRUCTURAL_DYNAMICS_SUPPORT_H

#include <rbl_rmatrix.h>
#include <rml_sparse_matrix.h>

//! Prescribed support displacement and its physical time derivatives.
struct RStructuralDynamicsSupportMotion
{
    RRVector displacement;
    RRVector velocity;
    RRVector acceleration;
};

//! Particular solution of the free-DOF constraints, q_f = T z + p(t).
struct RStructuralDynamicsConstraintMotion
{
    RRVector displacement;
    RRVector velocity;
    RRVector acceleration;
};

//! Physical displacement, velocity and acceleration after MPC reconstruction.
struct RStructuralDynamicsPhysicalState
{
    RRVector freeDisplacement;
    RRVector freeVelocity;
    RRVector freeAcceleration;
    RStructuralDynamicsSupportMotion support;
};

//! Algebraic moving-support coupling for M qddot + C qdot + K q = f.
//!
//! `free` denotes unknown DOFs after prescribed support DOFs have been
//! partitioned out. MPCs further reduce q_f = T z + p(t). All matrices retain
//! one consistent unit system (which may be solver-scaled). This class does
//! not own time histories or a linear solver, so a controller can evaluate it
//! repeatedly during a trial.
class RStructuralDynamicsSupportCoupling
{
    public:
        //! Return T^T [f_f - K_fc g - C_fc gdot - M_fc gddot
        //!              - K_ff p - C_ff pdot - M_ff pddot].
        static RRVector effectiveForce(const RSparseMatrix &freeStiffness,
                                       const RSparseMatrix &freeDamping,
                                       const RSparseMatrix &freeMass,
                                       const RSparseMatrix &supportStiffness,
                                       const RSparseMatrix &supportDamping,
                                       const RSparseMatrix &supportMass,
                                       const RRMatrix &nullspace,
                                       const RRVector &freeForce,
                                       const RStructuralDynamicsSupportMotion &support,
                                       const RStructuralDynamicsConstraintMotion &particular);

        //! Reconstruct the free and prescribed physical states from z/v/a.
        static RStructuralDynamicsPhysicalState reconstruct(const RRMatrix &nullspace,
                                                            const RRVector &reducedDisplacement,
                                                            const RRVector &reducedVelocity,
                                                            const RRVector &reducedAcceleration,
                                                            const RStructuralDynamicsSupportMotion &support,
                                                            const RStructuralDynamicsConstraintMotion &particular);
};

#endif // RSTRUCTURAL_DYNAMICS_SUPPORT_H

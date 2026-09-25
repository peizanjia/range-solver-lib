#ifndef RSTRUCTURAL_CABLE_EQUILIBRIUM_H
#define RSTRUCTURAL_CABLE_EQUILIBRIUM_H

#include <cstddef>
#include <vector>

#include "rstructural_cable_kernel.h"

//! One tension-only member in the pure sparse cable-network equilibrium model.
struct RStructuralCableEquilibriumMember
{
    std::size_t node1 = 0;
    std::size_t node2 = 0;
    double referenceLength = 0.0;
    double initialTension = 0.0;
    double axialRigidity = 0.0; //!< EA.
};

//! Nonlinear equilibrium controls. Convergence uses
//! ||f_internal-f_external||_free <= absForceTol + relForceTol*||f_external||_free.
//! The sparse tangent solve uses checked Jacobi-preconditioned CG.
struct RStructuralCableEquilibriumOptions
{
    std::size_t maxIterations = 50;
    std::size_t maxLineSearchIterations = 24;
    double absForceTol = 1.0e-8;
    double relForceTol = 1.0e-8;
};

//! Converged equilibrium state. Nodal arrays are flattened xyz in node order.
struct RStructuralCableEquilibriumResult
{
    std::vector<double> displacements;
    std::vector<double> internalForces;
    //! Full-DOF residual/reaction convention: internalForces - deadNodalLoads.
    std::vector<double> reactions;
    //! One nonnegative tension per input member, in input order.
    std::vector<double> elementTensions;
    //! Euclidean residual norm in newtons on unconstrained DOFs.
    double residualNorm = 0.0;
    //! Total potential: cable energies minus dead-load work.
    double energy = 0.0;
    //! Number of accepted Newton updates (zero for an initially balanced state).
    std::size_t iterations = 0;
};

//! Local, side-effect-free nonlinear force-balance solve for cable networks.
//!
//! Prescribed values and dead loads are flattened [node0 xyz, node1 xyz, ...].
//! Free DOFs start at zero displacement; prescribed values remain exact at all
//! trial points. The total potential uses the cable-kernel potential and
//! -f_ext dot u. Newton steps use Armijo backtracking; when that energy
//! comparison is unresolved at floating-point precision, a step may pass only
//! within 32 machine epsilons of the accumulated variable-potential scale and
//! with at least 10% physical residual reduction (or by reaching the same
//! outer force tolerance). This bounded tie-break does not change the final
//! residual certificate or make total energy evaluation exact. Slack
//! equilibria are force-balanced only: they carry no stability, uniqueness,
//! or mechanism certificate. This service is a static equilibrium kernel,
//! not a nonlinear dynamics implementation. Its local Jacobi-CG checks
//! curvature and recomputes the linear residual because the general matrix
//! solver has no result status for convergence certification.
class RStructuralCableEquilibriumSolver
{
    public:
        //! Throws std::invalid_argument for malformed/non-finite input and
        //! std::runtime_error for singular solves, exhausted line search, or
        //! failure to converge. No partial candidate state is returned.
        static RStructuralCableEquilibriumResult solve(
            const std::vector<RStructuralCableVector3> &referenceNodes,
            const std::vector<RStructuralCableEquilibriumMember> &members,
            const std::vector<bool> &prescribedMask,
            const std::vector<double> &prescribedDisplacements,
            const std::vector<double> &deadNodalLoads,
            const RStructuralCableEquilibriumOptions &options =
                RStructuralCableEquilibriumOptions{});
};

#endif // RSTRUCTURAL_CABLE_EQUILIBRIUM_H

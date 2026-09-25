#ifndef RSTRUCTURAL_CABLE_DYNAMICS_H
#define RSTRUCTURAL_CABLE_DYNAMICS_H

#include <cstddef>
#include <vector>

#include "rstructural_cable_equilibrium.h"

//! Newmark controls for one implicit, undamped cable-network step.
struct RStructuralCableDynamicsOptions
{
    double dt = 0.0;
    double beta = 0.25;
    double gamma = 0.5;
    std::size_t maxIterations = 50;
    std::size_t maxLineSearchIterations = 24;
    double absForceTol = 1.0e-8;
    double relForceTol = 1.0e-8;
};

//! Candidate state and mechanical work diagnostics for one accepted step.
struct RStructuralCableDynamicsResult
{
    std::vector<double> displacements;
    std::vector<double> velocities;
    std::vector<double> accelerations;
    std::vector<double> internalForces;
    //! Dynamic residual on each DOF. It equals the applied support force
    //! on prescribed DOFs under M*a + grad(U) = external + support;
    //! free entries are certified near zero.
    std::vector<double> reactions;
    //! One cable tension per input member, in input order.
    std::vector<double> elementTensions;
    double residualNorm = 0.0;
    double previousStrainEnergy = 0.0;
    double strainEnergy = 0.0;
    double previousKineticEnergy = 0.0;
    double kineticEnergy = 0.0;
    //! Constant dead/control force work over the step: f dot (u[n+1]-u[n]).
    double externalWorkIncrement = 0.0;
    //! Trapezoidal work by prescribed supports, using the endpoint dynamic
    //! residual as support force.
    double supportWorkIncrement = 0.0;
    //! Delta(strain+kinetic energy) minus external and support work. This is
    //! a diagnostic; nonlinear Newmark does not guarantee exact conservation.
    double energyBalanceError = 0.0;
    std::size_t iterations = 0;
};

//! Side-effect-free implicit Newmark step for a lumped-mass cable network.
//!
//! Nodal histories and loads use flattened xyz order. `nodalLumpedMass` has
//! one nonnegative mass per node, shared by its three translational DOFs.
//! Prescribed displacements are enforced exactly at the new time. The solver
//! uses the cable kernel's tension-only force and consistent material plus
//! geometric tangent, with no artificial stiffness in slack members. It
//! returns a candidate only after Newton convergence and a fresh free-residual
//! check. A thrown exception leaves all caller-owned accepted-step state
//! unchanged. Supported Newmark parameters satisfy gamma >= 0.5 and
//! beta >= (gamma+0.5)^2/4.
class RStructuralCableDynamicsSolver
{
    public:
        //! Throws std::invalid_argument for malformed/non-finite input and
        //! std::runtime_error for singular tangents, exhausted line search,
        //! or failure to converge. There is no partial-result output.
        static RStructuralCableDynamicsResult step(
            const std::vector<RStructuralCableVector3> &referenceNodes,
            const std::vector<RStructuralCableEquilibriumMember> &members,
            const std::vector<double> &nodalLumpedMass,
            const std::vector<double> &previousDisplacements,
            const std::vector<double> &previousVelocities,
            const std::vector<double> &previousAccelerations,
            const std::vector<bool> &prescribedMask,
            const std::vector<double> &prescribedDisplacements,
            const std::vector<double> &deadNodalLoads,
            const RStructuralCableDynamicsOptions &options);
};

#endif // RSTRUCTURAL_CABLE_DYNAMICS_H

#ifndef RSTRUCTURAL_PIEZO_COUPLED_1D_H
#define RSTRUCTURAL_PIEZO_COUPLED_1D_H

#include <QtGlobal>

#include <vector>

//! One small-strain axial piezoelectric rod in a 1-D nodal mesh.
//!
//! The directed connectivity node0 -> node1 defines the local piezo axis.
//! `coupling` is signed in that local direction. Reversing both endpoints
//! while preserving the same physical poling direction therefore also
//! requires reversing `coupling`.
struct RStructuralPiezoCoupled1DMember
{
    uint node0 = 0;
    uint node1 = 0;
    double area = 0.0;
    double youngModulus = 0.0;
    double coupling = 0.0;
    double permittivity = 0.0;
};

//! Complete input for a linear coupled 1-D piezoelectric FE solve.
//!
//! Mechanical loads use positive global-x force. Nodal conductor charges use
//! the physical free-charge sign of the terminal kernel (the charge at local
//! node1 is `-A*D`). At free-potential nodes, the equation is
//! `computedConductorCharge = prescribedConductorCharge`; at prescribed-
//! potential nodes, the returned conductor charge is the electrode reaction.
//! The supplied charge at a prescribed-potential node must be zero to avoid an
//! ambiguous mixture of voltage and charge boundary conditions.
struct RStructuralPiezoCoupled1DInput
{
    //! Global-x coordinate for each node. Members must have nonzero length.
    std::vector<double> coordinates;
    std::vector<RStructuralPiezoCoupled1DMember> members;

    //! All nodal vectors and masks have one entry per coordinate.
    std::vector<double> mechanicalLoads;
    std::vector<double> prescribedConductorCharges;
    std::vector<bool> displacementPrescribed;
    std::vector<double> prescribedDisplacements;
    std::vector<bool> potentialPrescribed;
    std::vector<double> prescribedPotentials;

    //! The mixed Hessian is indefinite, so the core solves it with GMRES.
    //! Each outer cycle performs up to gmresRestart Arnoldi iterations.
    unsigned int gmresRestart = 40;
    unsigned int gmresCycles = 20;
    double gmresRelativeTolerance = 1.0e-12;

    //! Independent physical residual tolerances. Units are N and C.
    double mechanicalAbsoluteTolerance = 1.0e-10;
    double mechanicalRelativeTolerance = 1.0e-10;
    double electricAbsoluteTolerance = 1.0e-14;
    double electricRelativeTolerance = 1.0e-10;
};

//! Coupled static result in the same node/member order as the input.
struct RStructuralPiezoCoupled1DResult
{
    std::vector<double> displacements;
    std::vector<double> potentials;
    //! Mechanical internal force in positive global-x coordinates.
    std::vector<double> internalForces;
    //! Physical free charge on each conductor node, including Dirichlet reactions.
    std::vector<double> conductorCharges;
    //! Signed force conjugate to elongation along each member's directed local axis.
    std::vector<double> memberAxialForces;
    //! Physical charge at each member's local node1 terminal.
    std::vector<double> memberRightCharges;
    //! Sum of element electric enthalpies.
    double enthalpy = 0.0;
    //! Norms of free mechanical-force and free electrical-charge residuals.
    double mechanicalResidualNorm = 0.0;
    double electricResidualNorm = 0.0;
};

//! Side-effect-free coupled small-strain 1-D piezoelectric finite-element solve.
//!
//! Both fields are solved together with sparse GMRES. The core does not use
//! conjugate gradients on the indefinite electric-enthalpy Hessian and does
//! not mutate a model or publish outputs. Every result is gated by freshly
//! recovered mechanical and electrical residuals in their own physical units.
class RStructuralPiezoCoupled1D
{
    public:
        static RStructuralPiezoCoupled1DResult solve(
            const RStructuralPiezoCoupled1DInput &input);
};

#endif // RSTRUCTURAL_PIEZO_COUPLED_1D_H

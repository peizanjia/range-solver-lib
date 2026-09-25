#ifndef RSTRUCTURAL_CABLE_KERNEL_H
#define RSTRUCTURAL_CABLE_KERNEL_H

#include <array>

using RStructuralCableVector3 = std::array<double,3>;
using RStructuralCableVector6 = std::array<double,6>;
using RStructuralCableMatrix6 = std::array<std::array<double,6>,6>;

//! Result for a two-node, 3-D tension-only truss/cable element.
struct RStructuralCableKernelResult
{
    double length = 0.0;
    double trialTension = 0.0;
    double tension = 0.0;
    double energy = 0.0;
    bool taut = false;
    RStructuralCableVector6 internalForce{};
    RStructuralCableMatrix6 tangent{};
};

//! Pure engineering-strain cable law and its consistent material/geometric tangent.
//!
//! The supplied initialTension is an element prestress parameter. It need not
//! describe a globally equilibrated structure; a prestressed cable network
//! must be force-balanced before it is treated as an unloaded equilibrium.
//! The conservative potential is max(0,Ntrial)^2/(2 EA/L0), with energy zero
//! at the slack/taut cutoff. At L0 it includes the prestress reference value
//! N0^2/(2 EA/L0); this offset does not affect force. Prestress forces still
//! need global nodal balance before the structure is an unloaded equilibrium.
class RStructuralCableKernel
{
    public:
        //! Positions and result ordering are node 1 xyz, followed by node 2 xyz.
        //! axialRigidity is EA. At the Ntrial == 0 kink, the slack (zero)
        //! tangent is selected explicitly.
        static RStructuralCableKernelResult evaluate(
            const RStructuralCableVector3 &node1Position,
            const RStructuralCableVector3 &node2Position,
            double referenceLength,
            double initialTension,
            double axialRigidity);
};

#endif // RSTRUCTURAL_CABLE_KERNEL_H

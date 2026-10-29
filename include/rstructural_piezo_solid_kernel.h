#ifndef RSTRUCTURAL_PIEZO_SOLID_KERNEL_H
#define RSTRUCTURAL_PIEZO_SOLID_KERNEL_H

#include "rstructural_piezo_layered_mesh.h"

//! Small-strain three-dimensional continuum with engineering Voigt strain.
//! Coordinates/material are in the same Cartesian frame. Nodes follow the
//! native TETRA1 / HEXA1 ordering; HEX8 admission is affine parallelepipeds.
//! G couples to NODAL potentials: force=K*u+G*phi, Q=C*phi-G^T*u.
//! No electrical inertia; the mesh owns equipotential electrodes/condensation.
struct RStructuralPiezoSolidNodalOperators
{
    RStructuralPiezoElementOperators mechanical;
    std::vector<double> coupling; //!< mechanical.dofs by nodes, row-major.
    std::vector<double> dielectric; //!< nodes by nodes, constant-potential nullspace.
};
class RStructuralPiezoSolidKernel
{
public:
    static RStructuralPiezoSolidNodalOperators assemble(
        const std::vector<std::array<double,3>> &nodes,
        const RStructuralPiezoContinuumMaterial &material);
};

//! Full 3D physical translations; rotations/moments are unsupported. Active
//! solids sharing a port form one electrically continuous domain. Different
//! ports are electrically insulated even if their mechanical nodes coincide.
//! Electrode faces lie on the port's common frame-z lowerZ/upperZ planes.
//! All other electrical boundary faces have the natural zero-free-charge BC.
RStructuralPiezoLayeredMeshOperators rStructuralPiezoAssembleSolidMesh(
    const RStructuralPiezoPlanarFrame &frame,
    const std::vector<RStructuralPiezoLayeredNode> &nodes,
    const std::vector<RStructuralPiezoLayeredMember> &members,
    double relativeTolerance=1e-10);

#endif

#ifndef RSTRUCTURAL_PIEZO_LAYERED_MESH_H
#define RSTRUCTURAL_PIEZO_LAYERED_MESH_H

#include "rstructural_piezo_element_operators.h"

//! Homogeneous beam/plate or solid mesh; physical coordinates/vectors are global SI.
enum class RStructuralPiezoLayeredFamily { Beam, Plate, Solid };
enum class RStructuralPiezoPhysicalDomain { Displacement, Velocity, ForceAndMoment };
struct RStructuralPiezoPlanarFrame
{
    std::array<double,3> origin{{0,0,0}},ex{{1,0,0}},ey{{0,1,0}},ez{{0,0,1}};
};
struct RStructuralPiezoLayeredNode
{
    uint physicalId = 0;
    std::array<double,3> position{};
};
struct RStructuralPiezoLayeredMember
{
    uint physicalId = 0;
    std::vector<uint> nodeIds; //!< Beam: 2; plate: 4 CCW; solid: native TET4 or HEX8 order.
    std::vector<RStructuralPiezoLayer> layers;
};
struct RStructuralPiezoLayeredMeshOperators
{
    RStructuralPiezoLayeredFamily family = RStructuralPiezoLayeredFamily::Beam;
    RStructuralPiezoElementOperators operators;
    std::vector<uint> nodeIds,memberIds; //!< Ascending stable physical IDs.
    std::vector<std::vector<uint>> memberNodeIndices; //!< Ordered element connectivity in nodeIds.
    //! Row-major (6*nodes) by dofs. qPhysical=P*q; supported M is never lumped.
    std::vector<double> physicalTransform;
    //! P^T*x; reject unsupported components separately per translation/rotation
    //! (or force/moment) block with |residual| <= relativeTolerance*|block|.
    //! No absolute tolerance mixes units. Zero blocks must be exactly zero.
    std::vector<double> projectPhysical(const std::vector<double> &physical,
        RStructuralPiezoPhysicalDomain domain,double relativeTolerance=1.0e-10) const;
    std::vector<double> reconstructPhysical(const std::vector<double> &generalized) const;
};
class RStructuralPiezoLayeredMesh
{
public:
    //! Atomic value-returning assembly using the production layer kernels. Frame
    //! axes are right handed and orthonormal; beam plane is x/z, plate plane x/y.
    //! Beam material x follows each chord, z=chord x common ey. Plate material
    //! axes follow the common frame. No torsion, other-plane bending or drilling.
    //! Hosts may meet along boundaries but must not overlap in their interiors.
    //! Solid delegates to the full 3D kernel with three global translations per
    //! node and frame-oriented material axes. It admits homogeneous TET4 or
    //! affine HEX8 meshes; one layer declares each whole solid volume/material.
    static RStructuralPiezoLayeredMeshOperators assemble(
        RStructuralPiezoLayeredFamily family,const RStructuralPiezoPlanarFrame &frame,
        const std::vector<RStructuralPiezoLayeredNode> &nodes,
        const std::vector<RStructuralPiezoLayeredMember> &members,
        double relativeTolerance=1.0e-10,double shearCorrection=5.0/6.0);
};
#endif

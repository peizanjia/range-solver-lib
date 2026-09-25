#ifndef RSTRUCTURAL_STATIC_TRUSS_MODEL_H
#define RSTRUCTURAL_STATIC_TRUSS_MODEL_H

#include <QtGlobal>

#include <array>
#include <vector>

class RModel;

//! One reference-configuration TRUSS1 member with its assigned line data.
struct RStructuralStaticTrussMember
{
    uint elementId = 0;
    uint node1 = 0;
    uint node2 = 0;
    double area = 0.0;
    double youngModulus = 0.0;
    //! Geometric length between reference mesh coordinates (not cable setup L0).
    double referenceLength = 0.0;
};

//! Validated common model inputs for the bounded static structural TRUSS1 paths.
struct RStructuralStaticTrussModelData
{
    std::vector<std::array<double,3>> nodes;
    //! Sorted by ascending elementId; point elements are intentionally omitted.
    std::vector<RStructuralStaticTrussMember> members;
    //! Prescribed translations and dead nodal loads are flattened xyz by node.
    std::vector<bool> prescribedMask;
    std::vector<double> prescribedDisplacements;
    std::vector<double> deadNodalLoads;
};

//! Shared static TRUSS1 model validation, boundary translation and mechanics output.
class RStructuralStaticTrussModel
{
    public:
        //! Validate common static restrictions and extract model data.
        //! Optional modes preserve all shared mesh/BC checks while allowing
        //! a caller-owned dynamic or steady Heat->Piezo solve.
        static RStructuralStaticTrussModelData build(
            const RModel &model,bool allowStructuralDynamics = false,
            bool allowThermalPiezo = false);

        //! Clear stale structural and electric-charge output fields.
        static void clearResults(RModel &model);

        //! Publish finite converged mechanics. `internalForces` are nodal internal
        //! forces rather than reactions. Signed member axial forces are ordered
        //! by ascending elementId; cable tensions are naturally nonnegative.
        static void publishMechanical(
            RModel &model,
            const std::vector<double> &displacements,
            const std::vector<double> &internalForces,
            const std::vector<double> &memberAxialForces,
            bool allowStructuralDynamics = false);
};

#endif // RSTRUCTURAL_STATIC_TRUSS_MODEL_H

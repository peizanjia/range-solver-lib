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
    //! Positive for legacy/passive/cable members. Zero for native effective
    //! active members, whose independent k/g/Ce law has no required material E.
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
    private:
        enum class ExtractionPolicy { Legacy, NativeMixedDynamics };
        static RStructuralStaticTrussModelData buildWithPolicy(
            const RModel &model,bool allowStructuralDynamics,
            bool allowThermalPiezo,ExtractionPolicy policy);
        static void publishWithPolicy(
            RModel &model,const std::vector<double> &displacements,
            const std::vector<double> &internalForces,
            const std::vector<double> &memberAxialForces,
            bool allowStructuralDynamics,ExtractionPolicy policy);

    public:
        //! Validate common static restrictions and extract model data.
        //! Optional modes preserve all shared mesh/BC checks while allowing
        //! a caller-owned dynamic or steady Heat->Piezo solve.
        static RStructuralStaticTrussModelData build(
            const RModel &model,bool allowStructuralDynamics = false,
            bool allowThermalPiezo = false);

        //! Explicit native extraction policy: dedicated restart, source q/v IC
        //! and analytic translational motion are owned by the mixed adapter.
        //! Every geometry/area/global BC/environment restriction still applies.
        static RStructuralStaticTrussModelData buildMixedDynamics(const RModel &model);

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

        //! Native output policy uses effective active laws and clears stale
        //! nodal electric potential; independent ports live in native records.
        //! The mixed adapter must certify the complete accepted state first.
        //! Analytic motion owns absolute q at the accepted physical time;
        //! this publisher never compares moving q with the static BC baseline.
        static void publishMechanicalMixed(
            RModel &model,const std::vector<double> &displacements,
            const std::vector<double> &internalForces,
            const std::vector<double> &memberAxialForces);
};

#endif // RSTRUCTURAL_STATIC_TRUSS_MODEL_H

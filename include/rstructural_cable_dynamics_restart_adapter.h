#ifndef RSTRUCTURAL_CABLE_DYNAMICS_RESTART_ADAPTER_H
#define RSTRUCTURAL_CABLE_DYNAMICS_RESTART_ADAPTER_H

#include <QString>

#include <vector>

#include "rstructural_cable_dynamics.h"

class RModel;
struct RStructuralDynamicsRestartData;

//! Persists the accepted full-nodal state of the bounded cable Newmark path.
//!
//! The existing structural restart record is reused; cable dynamics bind it
//! to their own physical data and keep the nodal result variables as display
//! output. The checkpoint vectors are authoritative for resumed integration.
class RStructuralCableDynamicsRestartAdapter
{
    public:
        //! Hash the exact data used by the cable transient FE path.
        static QString computeSystemFingerprint(
            const RModel &model,
            const std::vector<RStructuralCableVector3> &referenceNodes,
            const std::vector<RStructuralCableEquilibriumMember> &members,
            const std::vector<double> &elementAreas,
            const std::vector<double> &nodalLumpedMass,
            const std::vector<bool> &prescribedMask,
            const std::vector<double> &prescribedDisplacements,
            const std::vector<double> &externalForces,
            const RStructuralCableDynamicsOptions &options);

        //! Validate a persisted checkpoint against the current model setup.
        static void validateAcceptedState(
            const RModel &model,
            const QString &systemFingerprint,
            const std::vector<double> &externalForces);

        //! Use the exact persisted accepted state for the next trial step.
        static void restoreAcceptedState(
            const RModel &model,
            const QString &systemFingerprint,
            const std::vector<double> &externalForces,
            std::vector<double> &displacements,
            std::vector<double> &velocities,
            std::vector<double> &accelerations);

        //! Build a checkpoint candidate without changing the accepted model.
        static RStructuralDynamicsRestartData stageAcceptedState(
            const RModel &model,
            const RStructuralCableDynamicsResult &result,
            const QString &systemFingerprint,
            const std::vector<double> &externalForces);

        //! Commit a fully validated candidate after model result publication.
        static void commitAcceptedState(
            RModel &model,
            RStructuralDynamicsRestartData candidate);
};

#endif // RSTRUCTURAL_CABLE_DYNAMICS_RESTART_ADAPTER_H

#ifndef RSTRUCTURAL_CABLE_MODEL_ADAPTER_H
#define RSTRUCTURAL_CABLE_MODEL_ADAPTER_H

#include <QString>

#include "rstructural_cable_dynamics.h"

class RModel;

//! Validates, translates and publishes the bounded static/transient cable FE model.
class RStructuralCableModelAdapter
{
    public:
        //! Validate the current mutable model against the cable solver scope.
        static void validate(const RModel &model);

        //! Solve one converged static equilibrium without publishing results.
        static RStructuralCableEquilibriumResult solve(const RModel &model);

        //! Advance one Newmark-average-acceleration step without publishing it.
        static RStructuralCableDynamicsResult solveDynamics(
            const RModel &model,
            const std::vector<double> &previousDisplacements,
            const std::vector<double> &previousVelocities,
            const std::vector<double> &previousAccelerations,
            const std::vector<double> &additionalNodalForces = {},
            bool refreshAcceptedAcceleration = false);

        //! Initialize the fresh run's accepted acceleration from cable equilibrium.
        static void initializeFreshStartAcceleration(
            const RModel &model,
            const std::vector<double> &displacements,
            std::vector<double> &accelerations);

        //! Fingerprint the exact physical data used by the cable Newmark FE path.
        static QString computeDynamicsSystemFingerprint(const RModel &model);

        //! Return the flattened xyz mask of fixed cable translations.
        static std::vector<bool> prescribedDofMask(const RModel &model);

        //! Remove prior cable outputs so a failed repeat solve cannot leave a stale result.
        static void clearResults(RModel &model);

        //! Publish the converged state using the legacy FE result variables.
        static void publish(RModel &model,
                            const RStructuralCableEquilibriumResult &result);

        //! Publish one converged cable transient state and its nodal kinematics.
        static void publishDynamics(RModel &model,
                                    const RStructuralCableDynamicsResult &result);
};

#endif // RSTRUCTURAL_CABLE_MODEL_ADAPTER_H

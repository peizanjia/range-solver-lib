#ifndef RSTRUCTURAL_CABLE_MODEL_ADAPTER_H
#define RSTRUCTURAL_CABLE_MODEL_ADAPTER_H

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
            const std::vector<double> &previousAccelerations);

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

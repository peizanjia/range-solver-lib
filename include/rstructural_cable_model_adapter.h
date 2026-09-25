#ifndef RSTRUCTURAL_CABLE_MODEL_ADAPTER_H
#define RSTRUCTURAL_CABLE_MODEL_ADAPTER_H

#include "rstructural_cable_equilibrium.h"

class RModel;

//! Validates, translates and publishes the bounded static cable FE model.
class RStructuralCableModelAdapter
{
    public:
        //! Validate the current mutable model against the cable solver scope.
        static void validate(const RModel &model);

        //! Solve one converged static equilibrium without publishing results.
        static RStructuralCableEquilibriumResult solve(const RModel &model);

        //! Remove prior cable outputs so a failed repeat solve cannot leave a stale result.
        static void clearResults(RModel &model);

        //! Publish the converged state using the legacy FE result variables.
        static void publish(RModel &model,
                            const RStructuralCableEquilibriumResult &result);
};

#endif // RSTRUCTURAL_CABLE_MODEL_ADAPTER_H

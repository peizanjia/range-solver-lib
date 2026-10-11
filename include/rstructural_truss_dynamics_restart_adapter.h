#ifndef RSTRUCTURAL_TRUSS_DYNAMICS_RESTART_ADAPTER_H
#define RSTRUCTURAL_TRUSS_DYNAMICS_RESTART_ADAPTER_H

#include <rml_structural_truss_dynamics_restart_data.h>
#include "rstructural_truss_dynamics_model_adapter.h"

struct RStructuralTrussDynamicsRestoredState
{
    RStructuralTrussDynamicsState state;
    uint completedSteps = 0;
    uint nextCommandIndex = 0;
};

class RStructuralTrussDynamicsRestartAdapter
{
public:
    static RStructuralTrussDynamicsRestartData build(
        const RStructuralTrussDynamicsModelAdapter &adapter,
        const RStructuralTrussDynamicsState &state,
        uint completedSteps,uint nextCommandIndex);
    static RStructuralTrussDynamicsRestoredState restore(
        const RStructuralTrussDynamicsModelAdapter &adapter,
        const RStructuralTrussDynamicsRestartData &record);
    //! Exact binary64 complete-record SHA256 compatibility seal, not authentication.
    static QString stateFingerprint(const RStructuralTrussDynamicsRestartData &record);
};

#endif

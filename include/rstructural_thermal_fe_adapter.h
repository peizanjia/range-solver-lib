#ifndef RSTRUCTURAL_THERMAL_FE_ADAPTER_H
#define RSTRUCTURAL_THERMAL_FE_ADAPTER_H

#include "rstructural_thermal_control.h"

class RModel;

//! Same-mesh, in-process bridge from W heater commands to transient RModel FE.
//!
//! Each advance runs the transient heat FE once per scheduler thermal substep
//! and quasistatic stress FE on the final substep. It keeps candidate changes in a
//! private copy until the entire macrostep is certified. No result files are
//! written; the caller persists the accepted model and serialized scheduler
//! state together. IDs bind to ordered connectivity and reference coordinates
//! quantized to 1e-12 m, so mesh changes above that threshold fail closed.
class RStructuralThermalFEAdapter
{
    public:
        static QString nodeId(const RModel &model,unsigned int nodeIndex);
        //! Construct the scheduler from versioned RProblemSetup configuration.
        static RStructuralThermalControlScheduler schedulerFromModel(const RModel &model);
        static RStructuralThermalControlScheduler::Plant capture(const RModel &model);
        static RStructuralThermalControlScheduler::Result advance(
            RModel &model,
            const RStructuralThermalControlScheduler &scheduler,
            const RStructuralThermalControlState &accepted);
        //! Atomically publish one same-mesh model and scheduler state bundle.
        static void writeCheckpoint(
            const QString &fileName,const RModel &model,
            const RStructuralThermalControlScheduler &scheduler,
            const RStructuralThermalControlState &accepted);
        //! Validate the complete bundle before replacing either output.
        static void readCheckpoint(
            const QString &fileName,RModel &model,
            const RStructuralThermalControlScheduler &scheduler,
            RStructuralThermalControlState &accepted);
};

#endif // RSTRUCTURAL_THERMAL_FE_ADAPTER_H

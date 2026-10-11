#ifndef RSTRUCTURAL_NATIVE_CONTROL_H
#define RSTRUCTURAL_NATIVE_CONTROL_H

#include <rml_structural_native_control.h>
#include "rstructural_truss_dynamics_model_adapter.h"

struct RStructuralNativeControlProposal
{
    RStructuralNativeControlState state;
    RStructuralTrussDynamicsInputs completeInputs;
};

//! Pure built-in PI provider, N=1. No mutable history or device I/O.
//! Declared environment identity is limited to the canonical manifest below.
class RStructuralNativeControl
{
public:
    explicit RStructuralNativeControl(const RStructuralTrussDynamicsModelAdapter &adapter);
    static void stampBuiltinIdentity(RStructuralNativeControlSpec &spec);
    static QByteArray builtinEnvironmentManifest();
    RStructuralNativeControlState initialize() const;
    RStructuralNativeControlFrame sample(const RStructuralTrussDynamicsState &state,
                                         uint completedFE) const;
    RStructuralNativeControlProposal propose(const RStructuralNativeControlState &accepted,
        const RStructuralTrussDynamicsState &physical,uint completedFE) const;
    RStructuralTrussDynamicsInputs completeInputs(const std::vector<double> &values) const;
    //! Certify the latest pure proposal, clock, envelopes, limits and action phase.
    //! Historical left samples are checkpoint integrity data, not authenticated FE history.
    void validateAccepted(const RStructuralNativeControlState &state,uint completedFE) const;

private:
    const RStructuralTrussDynamicsModelAdapter &adapter;
    RStructuralNativeControlSpec spec;
    std::vector<std::size_t> sensorTargets,actionTargets,loopSensors,loopActions;
    std::vector<double> initialValues;
    struct PureProposal { std::vector<double> requested; QByteArray nextState; };
    PureProposal step(const RStructuralNativeControlFrame &frame,const QByteArray &state) const;
    RStructuralNativeControlActionRecord limit(std::size_t action,double previous,
                                               double requested,double physicalDt) const;
    void validateFrame(const RStructuralNativeControlFrame &frame,uint completedFE) const;
};

#endif // RSTRUCTURAL_NATIVE_CONTROL_H

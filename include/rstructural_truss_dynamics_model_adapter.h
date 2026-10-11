#ifndef RSTRUCTURAL_TRUSS_DYNAMICS_MODEL_ADAPTER_H
#define RSTRUCTURAL_TRUSS_DYNAMICS_MODEL_ADAPTER_H

#include <map>
#include <memory>
#include <QString>
#include <rml_structural_truss_dynamics_setup.h>
#include <rml_structural_native_control.h>
#include "rstructural_truss_dynamics.h"

class RModel;

struct RStructuralTrussDynamicsModelMaps
{
    std::vector<uint> orderedPhysicalNodeIds;
    std::vector<uint> orderedMemberPhysicalIds;
    std::vector<QString> orderedPortIds;
    std::map<uint,std::size_t> nodeByPhysicalId;
    std::map<uint,std::size_t> memberByPhysicalId;
    std::map<QString,std::size_t> memberByPortId;
};

struct RStructuralTrussDynamicsAction
{
    RStructuralTrussDynamicsInputs completeInputs;
    uint nextCommandIndex = 0;
};

struct RStructuralTrussDynamicsModelTrial
{
    RStructuralTrussDynamicsState state;
    RStructuralTrussDynamicsInputs completeInputs;
    uint nextCommandIndex = 0;
    RStructuralTrussDynamicsResult result;
};

//! Frozen native mixed model. All operations produce independent candidates.
//! Result variables never supply initial conditions or restore authority.
class RStructuralTrussDynamicsModelAdapter
{
public:
    explicit RStructuralTrussDynamicsModelAdapter(const RModel &model);
    const RStructuralTrussDynamics &driver() const { return *frozenDriver; }
    const RStructuralTrussDynamicsSystem &system() const { return frozenSystem; }
    const RStructuralTrussDynamicsModelMaps &maps() const { return frozenMaps; }
    const RStructuralTrussDynamicsModelMaps &getMaps() const { return frozenMaps; }
    const QString &sourceFingerprint() const { return frozenFingerprint; }
    double initialTime() const { return frozenInitialTime; }
    double nominalTimeStep() const { return frozenDt; }
    uint targetSteps() const { return frozenTargetSteps; }
    const RStructuralTrussDynamicsOptions &options() const { return frozenOptions; }
    const RStructuralNativeControlSpec &nativeControlSpec() const { return frozenControlSpec; }
    const RStructuralTrussDynamicsInputs &initialCompleteInputs() const { return initialInputs; }
    double targetTime(uint completedSteps) const;
    uint commandCursor(uint completedSteps) const;
    RStructuralTrussDynamicsState initialize() const;
    void validateCurrentSource(const RModel &model) const;
    RStructuralTrussDynamicsAction actionForStep(
        const RStructuralTrussDynamicsInputs &acceptedInputs,
        uint nextCommandIndex,uint stepIndex) const;
    RStructuralTrussDynamicsAction replayInputs(uint completedSteps) const;
    RStructuralTrussDynamicsModelTrial trial(
        const RStructuralTrussDynamicsState &state,
        uint completedSteps,uint nextCommandIndex) const;
    //! LocalV3 only: caller supplies every physical input; no table replay.
    RStructuralTrussDynamicsModelTrial trialWithCompleteInputs(
        const RStructuralTrussDynamicsState &state,uint completedSteps,
        const RStructuralTrussDynamicsInputs &inputs) const;
    //! Validate and publish into an independent copy of source.
    RModel candidateModel(const RModel &source,
        const RStructuralTrussDynamicsState &state,
        uint completedSteps,uint nextCommandIndex) const;
    RModel candidateModel(const RModel &source,
        const RStructuralTrussDynamicsState &state,
        uint completedSteps,uint nextCommandIndex,
        const RStructuralNativeControlState &controlState) const;

private:
    RStructuralTrussDynamicsSetup frozenSetup;
    RStructuralNativeControlSpec frozenControlSpec;
    RStructuralTrussDynamicsSystem frozenSystem;
    RStructuralTrussDynamicsModelMaps frozenMaps;
    RStructuralTrussDynamicsInputs initialInputs;
    std::vector<double> initialQ,initialV;
    QString frozenFingerprint;
    double frozenInitialTime = 0.0,frozenDt = 0.0;
    uint frozenTargetSteps = 0;
    RStructuralTrussDynamicsOptions frozenOptions;
    std::unique_ptr<RStructuralTrussDynamics> frozenDriver;
};

#endif

#ifndef RSTRUCTURAL_CONTROL_PROTOCOL_H
#define RSTRUCTURAL_CONTROL_PROTOCOL_H

#include <QByteArray>
#include <QString>

#include <vector>
#include <optional>

#include "rstructural_control_force_map.h"
#include "rstructural_control_runtime_parameters.h"

//! Physical-unit, stable-ID observations at an accepted solver boundary.
struct RStructuralControlSensorSample
{
    QString sensorId;
    double value = 0.0;
    QString unit;
};

struct RStructuralControlObservation
{
    double acceptedTime = 0.0;
    double nextTime = 0.0;
    quint64 deterministicSeed = 0;
    std::vector<RStructuralControlSensorSample> sensors;
    //! State after the previous accepted interval; callbacks must not mutate
    //! external state because the same observation may be replayed on retry.
    QByteArray priorControllerState;
    QByteArray priorRuntimeParametersState;
    quint64 acceptedSeq = 0;
};

struct RStructuralControlDecision
{
    //! Exactly one requested physical-newton value per catalog actuator ID.
    std::vector<RStructuralControlCommand> requestedCommands;
    //! Candidate state, published only after the result step commits.
    QByteArray nextControllerState;
    //! Canonical revision-zero bootstrap, only while no runtime state exists.
    QByteArray initialRuntimeParametersState;
    std::optional<RStructuralControlRuntimeParameters::Parameters> runtimeParameterCandidate;
    QString runtimeParameterProposalId;
};

struct RStructuralControlStepRecord
{
    double acceptedTime = 0.0;
    double nextTime = 0.0;
    std::vector<RStructuralControlSensorSample> sensors;
    std::vector<RStructuralControlCommand> requestedCommands;
    //! After deterministic delay/hold/noise, before amplitude/rate limiting.
    std::vector<RStructuralControlCommand> conditionedCommands;
    std::vector<RStructuralControlChannelSample> appliedChannels;
    QByteArray controllerStateSha256;
    //! Authoritative accepted snapshot, filled only after the step commits.
    QByteArray runtimeParametersState;
    quint64 appliedSeq = 0;
};

#endif // RSTRUCTURAL_CONTROL_PROTOCOL_H

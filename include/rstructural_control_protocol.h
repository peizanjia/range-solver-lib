#ifndef RSTRUCTURAL_CONTROL_PROTOCOL_H
#define RSTRUCTURAL_CONTROL_PROTOCOL_H

#include <QByteArray>
#include <QString>

#include <vector>

#include "rstructural_control_force_map.h"

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
};

struct RStructuralControlDecision
{
    //! Exactly one requested physical-newton value per catalog actuator ID.
    std::vector<RStructuralControlCommand> requestedCommands;
    //! Candidate state, published only after the result step commits.
    QByteArray nextControllerState;
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
};

#endif // RSTRUCTURAL_CONTROL_PROTOCOL_H

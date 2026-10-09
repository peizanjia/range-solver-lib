#ifndef RSTRUCTURAL_CONTROL_CONDITIONER_H
#define RSTRUCTURAL_CONTROL_CONDITIONER_H

#include <QByteArray>
#include <QString>

#include <vector>

#include "rstructural_control_force_map.h"

//! Stage-6 fixed-step conditioning policy for one ordered force channel.
//! noiseAmplitude is a uniform physical-newton half-width, not Gaussian sigma.
struct RStructuralControlConditioningPolicy
{
    QString actuatorId;
    uint delaySteps = 0;
    uint holdSteps = 1;
    double noiseAmplitude = 0.0;
};

//! Pure, deterministic command conditioning. It does not mutate a controller,
//! solver or model; callers publish the candidate state only on step accept.
//! The resulting command must still pass the solver's amplitude/rate limiter.
class RStructuralControlConditioner
{
    public:
        struct State
        {
            quint64 completedSteps = 0;
            QByteArray configurationSha256;
            std::vector<double> heldCommands;
            std::vector<std::vector<double>> delayedCommands;
        };

        struct Result
        {
            State candidate;
            std::vector<RStructuralControlCommand> appliedCommands;
        };

        RStructuralControlConditioner(quint64 deterministicSeed,
                                      std::vector<RStructuralControlConditioningPolicy> policies);
        State initialState() const;
        Result advance(const State &accepted,
                       const std::vector<RStructuralControlCommand> &requested) const;
        QByteArray serialize(const State &state) const;
        State deserialize(const QByteArray &bytes) const;
        const QByteArray &configurationSha256() const { return digest; }

    private:
        quint64 seed;
        std::vector<RStructuralControlConditioningPolicy> policies;
        QByteArray digest;
        void validate(const State &state) const;
};

//! Portable checkpoint envelope. The first component is byte-exact opaque
//! controller-private state; the second is conditioner state. JSON validation
//! belongs to a transport that explicitly negotiates JSON, not this kernel.
//! This is an opt-in building block, not an enabled FE capability by itself.
struct RStructuralControlStateEnvelope
{
    QByteArray controllerState;
    QByteArray conditionerState;
    //! Optional accepted runtime parameter state. Empty retains v1 exact bytes;
    //! a canonical nonempty component selects v2 without changing legacy state.
    QByteArray runtimeParametersState;

    QByteArray serialize() const;
    static RStructuralControlStateEnvelope deserialize(const QByteArray &bytes);
};

#endif // RSTRUCTURAL_CONTROL_CONDITIONER_H

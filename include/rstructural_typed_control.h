#ifndef RSTRUCTURAL_TYPED_CONTROL_H
#define RSTRUCTURAL_TYPED_CONTROL_H

#include <QByteArray>
#include <QString>

#include <functional>
#include <limits>
#include <vector>

#include "rstructural_control_protocol.h"

//! Physical quantity carried by a typed control channel.
enum class RStructuralTypedControlQuantity : quint8
{
    Force = 0,                    //!< N
    CableReferenceLength = 1,     //!< m
    CableInitialTension = 2,      //!< N
    CableAnchorDisplacement = 3,  //!< m
    Voltage = 4,                  //!< V
    Charge = 5                    //!< C
};

//! Return the canonical SI unit symbol for a typed control quantity.
QString rStructuralTypedControlUnit(RStructuralTypedControlQuantity quantity);

//! A stable-ID command whose value is expressed in its declared quantity.
struct RStructuralTypedControlCommand
{
    QString channelId;
    RStructuralTypedControlQuantity quantity = RStructuralTypedControlQuantity::Force;
    double value = 0.0;
};

//! Limits and fixed-step conditioning policy for one typed channel.
struct RStructuralTypedControlChannel
{
    QString channelId;
    RStructuralTypedControlQuantity quantity = RStructuralTypedControlQuantity::Force;
    double minimum = -std::numeric_limits<double>::infinity();
    double maximum = std::numeric_limits<double>::infinity();
    double rateLimit = std::numeric_limits<double>::infinity(); //!< Quantity units per second.
    double initialCommand = 0.0;
    uint delaySteps = 0;
    uint holdSteps = 1;
    double noiseAmplitude = 0.0; //!< Uniform half-width in this channel's quantity units.
};

//! Typed controller output; state is candidate data until its FE step accepts.
struct RStructuralTypedControlDecision
{
    std::vector<RStructuralTypedControlCommand> requestedCommands;
    QByteArray nextControllerState;
};

//! Per-channel typed values and generic work accounting for a candidate step.
struct RStructuralTypedControlChannelSample
{
    QString channelId;
    RStructuralTypedControlQuantity quantity = RStructuralTypedControlQuantity::Force;
    double requestedValue = 0.0;
    double conditionedValue = 0.0;
    double appliedValue = 0.0;
    double cumulativeWorkJ = 0.0;
    double averagePowerW = 0.0; //!< Work increment / fixedDt for this interval.
};

//! Pure, fixed-step typed control state for direct FE adapters.
//!
//! A session is bound to its ordered catalog, deterministic seed and fixed dt.
//! Its clock is initialTime + integerStep*dt with one std::fma rounding.
//! Poorly representable intervals are rejected rather than changing dt.
//! The optional controllerIdentity should identify the controller
//! implementation and revision because an opaque callback cannot be inspected.
//! Call makeObservation() at each trial so a retried callback sees the same
//! previously accepted opaque controller state. advance() returns an
//! unpublished candidate; the FE driver commits candidate.state only after
//! the physical step accepts. withWork() applies the FE driver's separate
//! per-channel work increments in joules. averagePowerW is the current
//! interval's work increment divided by fixedDt.
//!
//! Channel rateLimit is expressed in quantity units per second. Amplitude
//! limits are applied together with rate limits by their feasible interval.
//! Delay, hold, and uniform noise are applied before those physical limits.
class RStructuralTypedControlSession
{
    public:
        struct State
        {
            double initialTime = 0.0;
            double acceptedTime = 0.0;
            quint64 completedSteps = 0;
            QByteArray configurationSha256;
            QByteArray controllerState;
            std::vector<double> lastAppliedValues;
            std::vector<double> heldValues;
            std::vector<std::vector<double>> delayedValues;
            std::vector<double> cumulativeWorkJ;
        };

        //! Candidate data is local to the caller until its FE step accepts.
        struct Candidate
        {
            State state;
            std::vector<RStructuralTypedControlChannelSample> channels;
        };

        using ControllerCallback =
            std::function<RStructuralTypedControlDecision(const RStructuralControlObservation &)>;

        RStructuralTypedControlSession(quint64 deterministicSeed,
                                       double fixedDt,
                                       std::vector<RStructuralTypedControlChannel> catalog,
                                       QString controllerIdentity = QStringLiteral("cpp-typed-v1"));

        State initialState(double initialTime = 0.0,
                           QByteArray controllerState = {}) const;

        //! Build an existing Stage-6 observation using this accepted state.
        RStructuralControlObservation makeObservation(
            const State &accepted,
            std::vector<RStructuralControlSensorSample> sensors) const;

        //! Condition an already computed typed decision without mutating state.
        Candidate advance(const State &accepted,
                          const RStructuralControlObservation &observation,
                          const RStructuralTypedControlDecision &decision) const;

        //! Convenience callback trial; retries receive the same accepted state.
        Candidate trial(const State &accepted,
                        std::vector<RStructuralControlSensorSample> sensors,
                        const ControllerCallback &controller) const;

        //! Add FE-computed work increments (J) without changing accepted state.
        Candidate withWork(const Candidate &candidate,
                           const std::vector<double> &workIncrementsJ) const;

        //! Byte-exact little-endian state envelope with SHA-256 integrity check.
        QByteArray serialize(const State &state) const;
        State deserialize(const QByteArray &bytes) const;

        const QByteArray &configurationSha256() const { return this->digest; }
        const QString &controllerIdentity() const { return this->providerIdentity; }
        quint64 deterministicSeed() const { return this->seed; }
        double fixedDt() const { return this->timeStep; }
        const std::vector<RStructuralTypedControlChannel> &catalog() const
        {
            return this->channels;
        }

    private:
        quint64 seed;
        double timeStep;
        std::vector<RStructuralTypedControlChannel> channels;
        QString providerIdentity;
        QByteArray digest;

        void validateState(const State &state) const;
        void validateCandidate(const Candidate &candidate) const;
};

#endif // RSTRUCTURAL_TYPED_CONTROL_H

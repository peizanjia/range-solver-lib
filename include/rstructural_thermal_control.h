#ifndef RSTRUCTURAL_THERMAL_CONTROL_H
#define RSTRUCTURAL_THERMAL_CONTROL_H

#include <QByteArray>
#include <QString>

#include <functional>
#include <limits>
#include <vector>

#include "rstructural_thermo_coupled_1d.h"

//! Bang-bang thermostat and heater binding, all values in SI units.
struct RStructuralThermalControlChannelConfig
{
    QString actuatorId;
    QString heaterNodeId;
    QString temperatureSensorId;
    QString temperatureSensorNodeId;
    double turnOnBelowKelvin = 0.0;
    double turnOffAtOrAboveKelvin = 0.0;
    double onPowerWatts = 0.0;
    //! Lower bound for an on target; off remains 0 W and slew transitions may
    //! pass through [0, minimumPowerWatts].
    double minimumPowerWatts = 0.0;
    double maximumPowerWatts = 0.0;
    //! Maximum absolute command change per second; infinity disables rate limiting.
    double maximumPowerSlewWattsPerSecond = std::numeric_limits<double>::infinity();
};

struct RStructuralThermalControlConfig
{
    double thermalSubstepSeconds = 0.0;
    unsigned int thermalSubstepsPerMechanicalMacrostep = 1;
    unsigned int mechanicalMacrostepsPerControllerSample = 1;
    std::vector<RStructuralThermalControlChannelConfig> channels;
    //! Optional subset of nodes recorded as sparse temperature/displacement output.
    std::vector<QString> sparseObservationNodeIds;
};

//! A power command delivered to a stable actuator ID.
struct RStructuralThermalHeaterPower
{
    QString actuatorId;
    double watts = 0.0;
};

struct RStructuralThermalControlChannelState
{
    QString actuatorId;
    bool heaterOn = false;
    double requestedPowerWatts = 0.0;
    double conditionedPowerWatts = 0.0;
    double appliedPowerWatts = 0.0;
    double cumulativeRequestedEnergyJoules = 0.0;
    double cumulativeConditionedEnergyJoules = 0.0;
    double cumulativeAppliedEnergyJoules = 0.0;
};

//! Serializable accepted state. `plant` is keyed by stable physical node IDs.
struct RStructuralThermalControlState
{
    QByteArray configurationSha256;
    quint64 completedMechanicalMacrosteps = 0;
    quint64 completedControllerSamples = 0;
    std::vector<RStructuralThermoCoupled1DNodeState> plant;
    std::vector<RStructuralThermalControlChannelState> channels;
};

struct RStructuralThermalControlChannelSample
{
    QString actuatorId;
    double requestedPowerWatts = 0.0;
    double conditionedPowerWatts = 0.0;
    double appliedAveragePowerWatts = 0.0;
    double requestedEnergyJoules = 0.0;
    double conditionedEnergyJoules = 0.0;
    double appliedEnergyJoules = 0.0;
};

struct RStructuralThermalFieldObservation
{
    QString nodeId;
    double temperatureKelvin = 0.0;
    double axialDisplacementMetres = 0.0;
};

struct RStructuralThermalSensorObservation
{
    QString sensorId;
    QString nodeId;
    double temperatureKelvin = 0.0;
};

struct RStructuralThermalControlObservation
{
    //! Every accepted node, sorted by stable node ID.
    std::vector<RStructuralThermalFieldObservation> fullField;
    //! Configured sparse node subset, sorted by stable node ID.
    std::vector<RStructuralThermalFieldObservation> sparseField;
    //! Temperature-only control sensors, sorted by stable sensor ID.
    std::vector<RStructuralThermalSensorObservation> temperatureSensors;
};

struct RStructuralThermalControlSubstepResult
{
    std::vector<RStructuralThermoCoupled1DNodeState> candidatePlant;
    //! Actual power accepted by the backend, one entry per configured actuator.
    //! A slew transition may be below the on-target minimum, including 0 W.
    std::vector<RStructuralThermalHeaterPower> appliedHeaters;
};

//! In-process zero-order-hold heater control with transactional multirate steps.
//!
//! Each call advances one mechanical macrostep made of an integer number of
//! thermal substeps. The callback receives `solveMechanical == true` only on
//! the final thermal substep. Controller samples occur every configured
//! integer number of mechanical macrosteps. Callback results are candidates;
//! this scheduler never mutates or publishes the supplied accepted state.
//! The caller commits the returned candidate only after `advance` succeeds.
//!
//! The callback seam can be adapted to a model backend. The existing
//! RStructuralThermoCoupled1D API is included in the state contract, but that
//! coupled API itself has no independent thermal-only/mechanical-only solve.
class RStructuralThermalControlScheduler
{
    public:
        using Plant = std::vector<RStructuralThermoCoupled1DNodeState>;
        using SubstepFunction = std::function<RStructuralThermalControlSubstepResult(
            const Plant &acceptedPlant,
            const std::vector<RStructuralThermalHeaterPower> &conditionedHeaters,
            double dtSeconds,
            bool solveMechanical)>;

        struct Result
        {
            RStructuralThermalControlState candidate;
            std::vector<RStructuralThermalControlChannelSample> channels;
            RStructuralThermalControlObservation observation;
            bool controllerSampled = false;
        };

        explicit RStructuralThermalControlScheduler(RStructuralThermalControlConfig config);

        RStructuralThermalControlState initialState(const Plant &initialPlant) const;
        Result advance(const RStructuralThermalControlState &accepted,
                       const SubstepFunction &advanceSubstep) const;

        QByteArray serialize(const RStructuralThermalControlState &state) const;
        RStructuralThermalControlState deserialize(const QByteArray &bytes) const;
        const QByteArray &configurationSha256() const { return digest; }
        const RStructuralThermalControlConfig &getConfig() const { return config; }

    private:
        RStructuralThermalControlConfig config;
        QByteArray digest;

        void validateState(const RStructuralThermalControlState &state) const;
        RStructuralThermalControlObservation observe(const Plant &plant) const;
};

#endif // RSTRUCTURAL_THERMAL_CONTROL_H

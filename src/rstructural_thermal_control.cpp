#include "rstructural_thermal_control.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include <QCryptographicHash>
#include <QDataStream>
#include <QIODevice>
#include <QSet>

namespace
{
constexpr quint32 stateMagic = 0x52544843u; // RTHC
constexpr quint32 stateVersion = 1u;
constexpr qsizetype stateDigestBytes = 32;
constexpr qsizetype maxStateBytes = 64 * 1024 * 1024;
constexpr quint32 maxStateRecords = 1000000u;

void configureStream(QDataStream &stream)
{
    stream.setVersion(QDataStream::Qt_6_0);
    stream.setByteOrder(QDataStream::LittleEndian);
    stream.setFloatingPointPrecision(QDataStream::DoublePrecision);
}

bool validStableId(const QString &id)
{
    return !id.isEmpty() && id == id.trimmed();
}

void writeChannelConfig(QDataStream &stream,
                        const RStructuralThermalControlChannelConfig &channel)
{
    stream << channel.actuatorId << channel.heaterNodeId
           << channel.temperatureSensorId << channel.temperatureSensorNodeId
           << channel.turnOnBelowKelvin << channel.turnOffAtOrAboveKelvin
           << channel.onPowerWatts << channel.minimumPowerWatts
           << channel.maximumPowerWatts << channel.maximumPowerSlewWattsPerSecond;
}

void validatePlant(const std::vector<RStructuralThermoCoupled1DNodeState> &plant)
{
    if (plant.empty() || plant.size() > maxStateRecords)
        throw std::invalid_argument("Thermal control plant must contain a supported number of nodes.");
    QSet<QString> ids;
    for (const auto &node : plant)
    {
        if (!validStableId(node.nodeId) || ids.contains(node.nodeId))
            throw std::invalid_argument("Thermal control plant node IDs must be unique stable IDs.");
        if (!std::isfinite(node.temperatureKelvin) || node.temperatureKelvin <= 0.0 ||
            !std::isfinite(node.axialDisplacement))
            throw std::invalid_argument("Thermal control plant state must be finite with positive kelvin.");
        ids.insert(node.nodeId);
    }
}

void sortPlant(std::vector<RStructuralThermoCoupled1DNodeState> &plant)
{
    std::sort(plant.begin(),plant.end(),[](const auto &left,const auto &right) {
        return left.nodeId < right.nodeId;
    });
}

std::vector<RStructuralThermalHeaterPower> powersFromState(
    const std::vector<RStructuralThermalControlChannelState> &channels)
{
    std::vector<RStructuralThermalHeaterPower> powers;
    powers.reserve(channels.size());
    for (const auto &channel : channels)
        powers.push_back({channel.actuatorId,channel.conditionedPowerWatts});
    return powers;
}

void validateChannelState(const RStructuralThermalControlChannelState &state,
                          const RStructuralThermalControlChannelConfig &config)
{
    if (state.actuatorId != config.actuatorId ||
        !std::isfinite(state.requestedPowerWatts) ||
        !std::isfinite(state.conditionedPowerWatts) ||
        !std::isfinite(state.appliedPowerWatts) ||
        !std::isfinite(state.cumulativeRequestedEnergyJoules) ||
        !std::isfinite(state.cumulativeConditionedEnergyJoules) ||
        !std::isfinite(state.cumulativeAppliedEnergyJoules) ||
        state.cumulativeRequestedEnergyJoules < 0.0 ||
        state.cumulativeConditionedEnergyJoules < 0.0 ||
        state.cumulativeAppliedEnergyJoules < 0.0 ||
        state.requestedPowerWatts < 0.0 ||
        state.requestedPowerWatts !=
            (state.heaterOn ? std::clamp(config.onPowerWatts,
                                         config.minimumPowerWatts,
                                         config.maximumPowerWatts) : 0.0) ||
        state.conditionedPowerWatts < 0.0 ||
        state.conditionedPowerWatts > config.maximumPowerWatts ||
        state.appliedPowerWatts < 0.0 ||
        state.appliedPowerWatts > config.maximumPowerWatts)
        throw std::invalid_argument("Thermal control channel state is invalid or outside its power bounds.");
}
}

RStructuralThermalControlScheduler::RStructuralThermalControlScheduler(
    RStructuralThermalControlConfig configIn)
    : config(std::move(configIn))
{
    if (!std::isfinite(config.thermalSubstepSeconds) || config.thermalSubstepSeconds <= 0.0 ||
        config.thermalSubstepsPerMechanicalMacrostep == 0 ||
        config.mechanicalMacrostepsPerControllerSample == 0 || config.channels.empty() ||
        config.channels.size() > maxStateRecords ||
        config.sparseObservationNodeIds.size() > maxStateRecords)
        throw std::invalid_argument("Thermal control rates and channel count must be positive and bounded.");
    const double mechanicalPeriod = config.thermalSubstepSeconds *
                                   config.thermalSubstepsPerMechanicalMacrostep;
    if (!std::isfinite(mechanicalPeriod) ||
        !std::isfinite(mechanicalPeriod*config.mechanicalMacrostepsPerControllerSample))
        throw std::invalid_argument("Thermal and controller intervals must remain finite in seconds.");

    std::sort(config.channels.begin(),config.channels.end(),[](const auto &left,const auto &right) {
        return left.actuatorId < right.actuatorId;
    });
    QSet<QString> actuatorIds;
    QSet<QString> sensorIds;
    for (const auto &channel : config.channels)
    {
        if (!validStableId(channel.actuatorId) || !validStableId(channel.heaterNodeId) ||
            !validStableId(channel.temperatureSensorId) ||
            !validStableId(channel.temperatureSensorNodeId) ||
            actuatorIds.contains(channel.actuatorId) || sensorIds.contains(channel.temperatureSensorId))
            throw std::invalid_argument("Thermal channels require unique actuator and sensor IDs.");
        if (!std::isfinite(channel.turnOnBelowKelvin) || channel.turnOnBelowKelvin <= 0.0 ||
            !std::isfinite(channel.turnOffAtOrAboveKelvin) ||
            channel.turnOffAtOrAboveKelvin < channel.turnOnBelowKelvin ||
            !std::isfinite(channel.onPowerWatts) || channel.onPowerWatts < 0.0 ||
            !std::isfinite(channel.minimumPowerWatts) || channel.minimumPowerWatts < 0.0 ||
            !std::isfinite(channel.maximumPowerWatts) ||
            channel.maximumPowerWatts < channel.minimumPowerWatts ||
            channel.onPowerWatts < channel.minimumPowerWatts ||
            channel.onPowerWatts > channel.maximumPowerWatts ||
            std::isnan(channel.maximumPowerSlewWattsPerSecond) ||
            channel.maximumPowerSlewWattsPerSecond < 0.0)
            throw std::invalid_argument("Thermal channel thresholds, power and slew bounds are invalid.");
        actuatorIds.insert(channel.actuatorId);
        sensorIds.insert(channel.temperatureSensorId);
    }
    std::sort(config.sparseObservationNodeIds.begin(),config.sparseObservationNodeIds.end());
    for (std::size_t i=0;i<config.sparseObservationNodeIds.size();i++)
    {
        if (!validStableId(config.sparseObservationNodeIds[i]) ||
            (i > 0 && config.sparseObservationNodeIds[i] == config.sparseObservationNodeIds[i-1]))
            throw std::invalid_argument("Sparse observation node IDs must be unique stable IDs.");
    }

    QByteArray encoded;
    QDataStream stream(&encoded,QIODevice::WriteOnly);
    configureStream(stream);
    stream << quint32(1) << config.thermalSubstepSeconds
           << quint32(config.thermalSubstepsPerMechanicalMacrostep)
           << quint32(config.mechanicalMacrostepsPerControllerSample)
           << quint32(config.channels.size());
    for (const auto &channel : config.channels) writeChannelConfig(stream,channel);
    stream << quint32(config.sparseObservationNodeIds.size());
    for (const QString &nodeId : config.sparseObservationNodeIds) stream << nodeId;
    if (stream.status() != QDataStream::Ok)
        throw std::runtime_error("Could not encode thermal control configuration.");
    digest = QCryptographicHash::hash(encoded,QCryptographicHash::Sha256);
}

RStructuralThermalControlState RStructuralThermalControlScheduler::initialState(
    const Plant &initialPlant) const
{
    validatePlant(initialPlant);
    RStructuralThermalControlState state;
    state.configurationSha256 = digest;
    state.plant = initialPlant;
    sortPlant(state.plant);
    for (const auto &channel : config.channels)
    {
        state.channels.push_back({channel.actuatorId,false,0.0,0.0,0.0,0.0,0.0,0.0});
    }
    validateState(state);
    return state;
}

void RStructuralThermalControlScheduler::validateState(
    const RStructuralThermalControlState &state) const
{
    if (state.configurationSha256 != digest || state.channels.size() != config.channels.size())
        throw std::invalid_argument("Thermal control state is bound to a different configuration.");
    const quint64 expectedSamples =
        state.completedMechanicalMacrosteps/config.mechanicalMacrostepsPerControllerSample +
        (state.completedMechanicalMacrosteps%config.mechanicalMacrostepsPerControllerSample != 0 ? 1u : 0u);
    if (state.completedControllerSamples != expectedSamples)
        throw std::invalid_argument("Thermal control sample count does not match accepted macrosteps.");
    validatePlant(state.plant);
    if (!std::is_sorted(state.plant.begin(),state.plant.end(),[](const auto &left,const auto &right) {
            return left.nodeId < right.nodeId;
        }))
        throw std::invalid_argument("Thermal control state nodes must be sorted by stable ID.");
    QSet<QString> nodeIds;
    for (const auto &node : state.plant) nodeIds.insert(node.nodeId);
    for (std::size_t i=0;i<config.channels.size();i++)
    {
        const auto &channelConfig = config.channels[i];
        const auto &channelState = state.channels[i];
        validateChannelState(channelState,channelConfig);
        if (!nodeIds.contains(channelConfig.heaterNodeId) ||
            !nodeIds.contains(channelConfig.temperatureSensorNodeId))
            throw std::invalid_argument("Thermal control actuator or sensor node is absent from the plant.");
    }
    for (const QString &nodeId : config.sparseObservationNodeIds)
    {
        if (!nodeIds.contains(nodeId))
            throw std::invalid_argument("Sparse observation node is absent from the plant.");
    }
}

RStructuralThermalControlObservation RStructuralThermalControlScheduler::observe(
    const Plant &plant) const
{
    RStructuralThermalControlObservation result;
    result.fullField.reserve(plant.size());
    for (const auto &node : plant)
        result.fullField.push_back({node.nodeId,node.temperatureKelvin,node.axialDisplacement});

    result.sparseField.reserve(config.sparseObservationNodeIds.size());
    for (const QString &nodeId : config.sparseObservationNodeIds)
    {
        const auto it = std::lower_bound(plant.begin(),plant.end(),nodeId,
                                         [](const auto &node,const QString &id) {
                                             return node.nodeId < id;
                                         });
        result.sparseField.push_back({it->nodeId,it->temperatureKelvin,it->axialDisplacement});
    }

    result.temperatureSensors.reserve(config.channels.size());
    for (const auto &channel : config.channels)
    {
        const auto it = std::lower_bound(plant.begin(),plant.end(),channel.temperatureSensorNodeId,
                                         [](const auto &node,const QString &id) {
                                             return node.nodeId < id;
                                         });
        result.temperatureSensors.push_back({channel.temperatureSensorId,it->nodeId,
                                             it->temperatureKelvin});
    }
    std::sort(result.temperatureSensors.begin(),result.temperatureSensors.end(),
              [](const auto &left,const auto &right) { return left.sensorId < right.sensorId; });
    return result;
}

RStructuralThermalControlScheduler::Result RStructuralThermalControlScheduler::advance(
    const RStructuralThermalControlState &accepted,
    const SubstepFunction &advanceSubstep) const
{
    validateState(accepted);
    if (!advanceSubstep)
        throw std::invalid_argument("Thermal control scheduler requires a substep callback.");
    if (accepted.completedMechanicalMacrosteps == std::numeric_limits<quint64>::max())
        throw std::overflow_error("Thermal control macrostep counter is exhausted.");

    Result result;
    result.candidate = accepted;
    const double macrostepSeconds = config.thermalSubstepSeconds *
                                   config.thermalSubstepsPerMechanicalMacrostep;
    const bool sample = accepted.completedMechanicalMacrosteps %
                        config.mechanicalMacrostepsPerControllerSample == 0;
    result.controllerSampled = sample;
    if (sample)
    {
        if (result.candidate.completedControllerSamples == std::numeric_limits<quint64>::max())
            throw std::overflow_error("Thermal control sample counter is exhausted.");
        const double controllerPeriodSeconds = macrostepSeconds *
                                               config.mechanicalMacrostepsPerControllerSample;
        for (std::size_t i=0;i<config.channels.size();i++)
        {
            const auto &channelConfig = config.channels[i];
            auto &channelState = result.candidate.channels[i];
            const auto sensor = std::lower_bound(accepted.plant.begin(),accepted.plant.end(),
                                                 channelConfig.temperatureSensorNodeId,
                                                 [](const auto &node,const QString &id) {
                                                     return node.nodeId < id;
                                                 });
            if (sensor == accepted.plant.end() || sensor->nodeId != channelConfig.temperatureSensorNodeId)
                throw std::invalid_argument("Thermal control sensor node is absent from the accepted plant.");
            if (sensor->temperatureKelvin < channelConfig.turnOnBelowKelvin)
                channelState.heaterOn = true;
            else if (sensor->temperatureKelvin >= channelConfig.turnOffAtOrAboveKelvin)
                channelState.heaterOn = false;

            channelState.requestedPowerWatts = channelState.heaterOn ? channelConfig.onPowerWatts : 0.0;
            const double boundedRequest = channelState.heaterOn
                ? std::clamp(channelState.requestedPowerWatts,
                             channelConfig.minimumPowerWatts,
                             channelConfig.maximumPowerWatts)
                : 0.0;
            double conditioned = boundedRequest;
            if (std::isfinite(channelConfig.maximumPowerSlewWattsPerSecond))
            {
                const double maximumChange = channelConfig.maximumPowerSlewWattsPerSecond *
                                             controllerPeriodSeconds;
                conditioned = std::clamp(boundedRequest,
                                         channelState.conditionedPowerWatts-maximumChange,
                                         channelState.conditionedPowerWatts+maximumChange);
            }
            channelState.conditionedPowerWatts = std::clamp(conditioned,
                                                            0.0,
                                                            channelConfig.maximumPowerWatts);
        }
        ++result.candidate.completedControllerSamples;
    }

    result.channels.resize(config.channels.size());
    for (std::size_t i=0;i<config.channels.size();i++)
        result.channels[i].actuatorId = config.channels[i].actuatorId;

    Plant workingPlant = accepted.plant;
    const auto heldCommands = powersFromState(result.candidate.channels);
    for (unsigned int substep=0;substep<config.thermalSubstepsPerMechanicalMacrostep;substep++)
    {
        const bool solveMechanical = substep+1 == config.thermalSubstepsPerMechanicalMacrostep;
        // Pass an isolated copy: even an ill-behaved backend cannot mutate the
        // accepted or scheduler-owned candidate through this argument.
        const Plant callbackInput = workingPlant;
        auto substepResult = advanceSubstep(callbackInput,heldCommands,
                                            config.thermalSubstepSeconds,solveMechanical);
        validatePlant(substepResult.candidatePlant);
        sortPlant(substepResult.candidatePlant);
        if (substepResult.candidatePlant.size() != workingPlant.size())
            throw std::runtime_error("Thermal backend changed the plant node set during a step.");
        for (std::size_t nodeIndex=0;nodeIndex<workingPlant.size();nodeIndex++)
        {
            const auto &before = workingPlant[nodeIndex];
            const auto &after = substepResult.candidatePlant[nodeIndex];
            if (before.nodeId != after.nodeId)
                throw std::runtime_error("Thermal backend changed a stable plant node ID during a step.");
            if (!solveMechanical && before.axialDisplacement != after.axialDisplacement)
                throw std::runtime_error("Thermal backend changed displacement outside a mechanical macrostep.");
        }

        if (substepResult.appliedHeaters.size() != config.channels.size())
            throw std::runtime_error("Thermal backend returned an incomplete applied-heater ledger.");
        std::vector<double> applied(config.channels.size(),0.0);
        std::vector<bool> seen(config.channels.size(),false);
        for (const auto &actual : substepResult.appliedHeaters)
        {
            auto channelIt = std::lower_bound(config.channels.begin(),config.channels.end(),actual.actuatorId,
                                              [](const auto &channel,const QString &id) {
                                                  return channel.actuatorId < id;
                                              });
            if (channelIt == config.channels.end() || channelIt->actuatorId != actual.actuatorId ||
                !std::isfinite(actual.watts))
                throw std::runtime_error("Thermal backend returned an unknown or non-finite heater power.");
            const std::size_t channelIndex = std::size_t(channelIt-config.channels.begin());
            if (seen[channelIndex] || actual.watts < 0.0 ||
                actual.watts > channelIt->maximumPowerWatts)
                throw std::runtime_error("Thermal backend returned duplicate or out-of-bounds heater power.");
            seen[channelIndex] = true;
            applied[channelIndex] = actual.watts;
        }
        if (std::find(seen.begin(),seen.end(),false) != seen.end())
            throw std::runtime_error("Thermal backend omitted an applied heater channel.");

        for (std::size_t i=0;i<config.channels.size();i++)
        {
            auto &channelState = result.candidate.channels[i];
            auto &sampleState = result.channels[i];
            const double dt = config.thermalSubstepSeconds;
            const double requestedEnergy = channelState.requestedPowerWatts*dt;
            const double conditionedEnergy = channelState.conditionedPowerWatts*dt;
            const double appliedEnergy = applied[i]*dt;
            channelState.appliedPowerWatts = applied[i];
            channelState.cumulativeRequestedEnergyJoules += requestedEnergy;
            channelState.cumulativeConditionedEnergyJoules += conditionedEnergy;
            channelState.cumulativeAppliedEnergyJoules += appliedEnergy;
            sampleState.requestedPowerWatts = channelState.requestedPowerWatts;
            sampleState.conditionedPowerWatts = channelState.conditionedPowerWatts;
            sampleState.requestedEnergyJoules += requestedEnergy;
            sampleState.conditionedEnergyJoules += conditionedEnergy;
            sampleState.appliedEnergyJoules += appliedEnergy;
        }
        workingPlant = std::move(substepResult.candidatePlant);
    }

    for (std::size_t i=0;i<config.channels.size();i++)
        result.channels[i].appliedAveragePowerWatts =
            result.channels[i].appliedEnergyJoules/macrostepSeconds;
    result.candidate.plant = std::move(workingPlant);
    ++result.candidate.completedMechanicalMacrosteps;
    for (std::size_t i=0;i<config.channels.size();i++)
        validateChannelState(result.candidate.channels[i],config.channels[i]);
    validateState(result.candidate);
    result.observation = observe(result.candidate.plant);
    return result;
}

QByteArray RStructuralThermalControlScheduler::serialize(
    const RStructuralThermalControlState &state) const
{
    validateState(state);
    QByteArray payload;
    QDataStream stream(&payload,QIODevice::WriteOnly);
    configureStream(stream);
    stream << stateMagic << stateVersion << state.configurationSha256
           << state.completedMechanicalMacrosteps << state.completedControllerSamples
           << quint32(state.plant.size());
    for (const auto &node : state.plant)
        stream << node.nodeId << node.temperatureKelvin << node.axialDisplacement;
    stream << quint32(state.channels.size());
    for (const auto &channel : state.channels)
        stream << channel.actuatorId << channel.heaterOn << channel.requestedPowerWatts
               << channel.conditionedPowerWatts << channel.appliedPowerWatts
               << channel.cumulativeRequestedEnergyJoules
               << channel.cumulativeConditionedEnergyJoules
               << channel.cumulativeAppliedEnergyJoules;
    if (stream.status() != QDataStream::Ok || payload.size() > maxStateBytes-stateDigestBytes)
        throw std::runtime_error("Could not encode thermal control state within its size limit.");
    payload.append(QCryptographicHash::hash(payload,QCryptographicHash::Sha256));
    return payload;
}

RStructuralThermalControlState RStructuralThermalControlScheduler::deserialize(
    const QByteArray &bytes) const
{
    if (bytes.size() < stateDigestBytes || bytes.size() > maxStateBytes)
        throw std::invalid_argument("Thermal control state is truncated or exceeds its size limit.");
    const QByteArray payload = bytes.left(bytes.size()-stateDigestBytes);
    const QByteArray checksum = bytes.right(stateDigestBytes);
    if (QCryptographicHash::hash(payload,QCryptographicHash::Sha256) != checksum)
        throw std::invalid_argument("Thermal control state checksum is invalid.");

    QDataStream stream(payload);
    configureStream(stream);
    quint32 magic = 0;
    quint32 version = 0;
    RStructuralThermalControlState state;
    quint32 nodeCount = 0;
    stream >> magic >> version >> state.configurationSha256
           >> state.completedMechanicalMacrosteps >> state.completedControllerSamples >> nodeCount;
    if (stream.status() != QDataStream::Ok || magic != stateMagic || version != stateVersion ||
        nodeCount == 0 || nodeCount > maxStateRecords)
        throw std::invalid_argument("Thermal control state header is invalid or truncated.");
    state.plant.resize(nodeCount);
    for (auto &node : state.plant)
        stream >> node.nodeId >> node.temperatureKelvin >> node.axialDisplacement;
    quint32 channelCount = 0;
    stream >> channelCount;
    if (stream.status() != QDataStream::Ok || channelCount > maxStateRecords)
        throw std::invalid_argument("Thermal control channel record count is invalid.");
    state.channels.resize(channelCount);
    for (auto &channel : state.channels)
        stream >> channel.actuatorId >> channel.heaterOn >> channel.requestedPowerWatts
               >> channel.conditionedPowerWatts >> channel.appliedPowerWatts
               >> channel.cumulativeRequestedEnergyJoules
               >> channel.cumulativeConditionedEnergyJoules
               >> channel.cumulativeAppliedEnergyJoules;
    if (stream.status() != QDataStream::Ok || !stream.atEnd())
        throw std::invalid_argument("Thermal control state is truncated or has trailing records.");
    validateState(state);
    return state;
}

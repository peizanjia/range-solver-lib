#include "rstructural_control_force_map.h"

#include <QSet>

#include <cmath>
#include <stdexcept>

RStructuralControlForceMap::RStructuralControlForceMap(
    uint reducedSize, double forceScale, double velocityScale,
    std::vector<RStructuralControlChannelMap> channels)
    : reducedSize(reducedSize), forceScale(forceScale), velocityScale(velocityScale),
      channels(std::move(channels))
{
    if (reducedSize == 0 || !std::isfinite(forceScale) || forceScale <= 0.0 ||
        !std::isfinite(velocityScale) || velocityScale <= 0.0)
        throw std::invalid_argument("Invalid structural control map dimensions or scales");
    QSet<QString> ids;
    for (const auto &channel : this->channels)
    {
        if (channel.actuatorId.isEmpty() || ids.contains(channel.actuatorId) ||
            channel.reducedWeights.empty())
            throw std::invalid_argument("Structural control actuator IDs and maps must be nonempty and unique");
        ids.insert(channel.actuatorId);
        for (const auto &weight : channel.reducedWeights)
        {
            if (weight.first >= reducedSize || !std::isfinite(weight.second))
                throw std::invalid_argument("Structural control map has invalid reduced weight");
        }
    }
}

void RStructuralControlForceMap::validateCommands(
    const std::vector<RStructuralControlCommand> &commands) const
{
    if (commands.size() != this->channels.size())
        throw std::invalid_argument("Structural control requires one command per mapped actuator");
    for (size_t i=0;i<commands.size();i++)
    {
        if (commands[i].actuatorId != this->channels[i].actuatorId ||
            !std::isfinite(commands[i].value))
            throw std::invalid_argument("Structural control command ID/order/value is invalid");
    }
}

void RStructuralControlForceMap::addForce(
    const std::vector<RStructuralControlCommand> &commands, RRVector &reducedForce) const
{
    this->validateCommands(commands);
    if (reducedForce.size() != this->reducedSize)
        throw std::invalid_argument("Structural control force vector has wrong size");
    for (size_t i=0;i<commands.size();i++)
    {
        const double scaled = commands[i].value*this->forceScale;
        if (!std::isfinite(scaled))
            throw std::invalid_argument("Structural control command overflows scaled force");
        for (const auto &weight : this->channels[i].reducedWeights)
            reducedForce[weight.first] += scaled*weight.second;
    }
}

std::vector<double> RStructuralControlForceMap::powers(
    const std::vector<RStructuralControlCommand> &commands,
    const RRVector &reducedVelocity) const
{
    this->validateCommands(commands);
    if (reducedVelocity.size() != this->reducedSize)
        throw std::invalid_argument("Structural control velocity vector has wrong size");
    std::vector<double> result(commands.size(),0.0);
    for (size_t i=0;i<commands.size();i++)
    {
        double velocity = 0.0;
        for (const auto &weight : this->channels[i].reducedWeights)
            velocity += weight.second*reducedVelocity[weight.first];
        result[i] = commands[i].value*velocity/this->velocityScale;
        if (!std::isfinite(result[i]))
            throw std::invalid_argument("Structural control physical power is non-finite");
    }
    return result;
}

RStructuralControlLedger::RStructuralControlLedger(const RStructuralControlForceMap &map)
{
    for (const auto &channel : map.getChannels())
        this->committed.push_back({channel.actuatorId,0.0,0.0,0.0});
}

const std::vector<RStructuralControlChannelSample> &RStructuralControlLedger::trial(
    const RStructuralControlForceMap &map,
    const std::vector<RStructuralControlCommand> &commands,
    const RRVector &startingVelocity,
    const RRVector &endingVelocity,
    double dt)
{
    this->pendingReady = false;
    if (!std::isfinite(dt) || dt <= 0.0 || commands.size() != this->committed.size())
        throw std::invalid_argument("Invalid structural control ledger step");
    const auto startPower = map.powers(commands,startingVelocity);
    const auto endPower = map.powers(commands,endingVelocity);
    auto next = this->committed;
    for (size_t i=0;i<next.size();i++)
    {
        if (commands[i].actuatorId != next[i].actuatorId)
            throw std::invalid_argument("Structural control ledger channel order changed");
        next[i].command = commands[i].value;
        next[i].endpointPower = endPower[i];
        next[i].cumulativeWork += 0.5*(startPower[i]+endPower[i])*dt;
        if (!std::isfinite(next[i].cumulativeWork))
            throw std::invalid_argument("Structural control work is non-finite");
    }
    this->pending = std::move(next);
    this->pendingReady = true;
    return this->pending;
}

const std::vector<RStructuralControlChannelSample> &RStructuralControlLedger::candidate() const
{
    if (!this->pendingReady) throw std::logic_error("No candidate structural control ledger");
    return this->pending;
}

void RStructuralControlLedger::accept()
{
    if (!this->pendingReady) throw std::logic_error("No candidate structural control ledger");
    this->committed = std::move(this->pending);
    this->pendingReady = false;
}

void RStructuralControlLedger::reject()
{
    this->pending.clear();
    this->pendingReady = false;
}

void RStructuralControlLedger::restore(
    const std::vector<RStructuralControlChannelSample> &samples)
{
    if (samples.size() != this->committed.size())
        throw std::invalid_argument("Structural control checkpoint channel count changed");
    for (size_t i=0;i<samples.size();i++)
    {
        if (samples[i].actuatorId != this->committed[i].actuatorId ||
            !std::isfinite(samples[i].command) || !std::isfinite(samples[i].endpointPower) ||
            !std::isfinite(samples[i].cumulativeWork))
            throw std::invalid_argument("Invalid structural control checkpoint sample");
    }
    this->committed = samples;
    this->reject();
}

#ifndef RSTRUCTURAL_CONTROL_FORCE_MAP_H
#define RSTRUCTURAL_CONTROL_FORCE_MAP_H

#include <QString>

#include <utility>
#include <vector>

#include <rml_sparse_matrix.h>

//! Solver-internal sparse image of one stable Stage-4 nodal actuator ID.
//! Weights map a scalar physical-newton command into reduced solver DOFs,
//! after constraints and local equation transformations are resolved once.
struct RStructuralControlChannelMap
{
    QString actuatorId;
    std::vector<std::pair<uint,double>> reducedWeights;
};

struct RStructuralControlCommand
{
    QString actuatorId;
    double value = 0.0; //!< Physical newtons along the catalog direction.
};

//! O(sum of mapped nonzeros) RHS injection and power sampling, no dense matrix.
class RStructuralControlForceMap
{
    public:
        RStructuralControlForceMap(uint reducedSize,
                                   double forceScale,
                                   double velocityScale,
                                   std::vector<RStructuralControlChannelMap> channels);

        const std::vector<RStructuralControlChannelMap> &getChannels() const { return channels; }
        uint getReducedSize() const { return reducedSize; }
        void addForce(const std::vector<RStructuralControlCommand> &commands,
                      RRVector &reducedForce) const;
        //! Physical watts at a velocity endpoint for a held command.
        std::vector<double> powers(const std::vector<RStructuralControlCommand> &commands,
                                   const RRVector &reducedVelocity) const;

    private:
        uint reducedSize;
        double forceScale;
        double velocityScale;
        std::vector<RStructuralControlChannelMap> channels;
        void validateCommands(const std::vector<RStructuralControlCommand> &commands) const;
};

struct RStructuralControlChannelSample
{
    QString actuatorId;
    double command = 0.0;        //!< N, held over the just accepted interval.
    double endpointPower = 0.0;  //!< W at the accepted interval endpoint.
    double cumulativeWork = 0.0; //!< J, trapezoid of held-command endpoint powers.
};

//! Trial/accept/reject ledger, independent of a controller's opaque state.
class RStructuralControlLedger
{
    public:
        explicit RStructuralControlLedger(const RStructuralControlForceMap &map);
        const std::vector<RStructuralControlChannelSample> &accepted() const { return committed; }
        const std::vector<RStructuralControlChannelSample> &trial(
            const RStructuralControlForceMap &map,
            const std::vector<RStructuralControlCommand> &commands,
            const RRVector &startingVelocity,
            const RRVector &endingVelocity,
            double dt);
        const std::vector<RStructuralControlChannelSample> &candidate() const;
        void accept();
        void reject();
        void restore(const std::vector<RStructuralControlChannelSample> &samples);

    private:
        std::vector<RStructuralControlChannelSample> committed;
        std::vector<RStructuralControlChannelSample> pending;
        bool pendingReady = false;
};

#endif // RSTRUCTURAL_CONTROL_FORCE_MAP_H

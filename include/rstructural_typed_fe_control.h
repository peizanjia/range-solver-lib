#ifndef RSTRUCTURAL_TYPED_FE_CONTROL_H
#define RSTRUCTURAL_TYPED_FE_CONTROL_H

#include <functional>
#include <memory>
#include <vector>

#include "rstructural_piezo_beam_fe_adapter.h"
#include "rstructural_typed_control.h"

class RModel;

//! One catalog channel bound to a native member, patch, or prescribed anchor.
struct RStructuralTypedFEControlBinding
{
    QString channelId;
    uint elementId = 0;
    uint nodeId = 0;
    uint axis = 0;
};

//! Accepted-boundary controller bridge to actual quasistatic cable / piezo FE.
//!
//! The driver owns its model, controller state and output ledger. A callback or
//! solve failure leaves the whole accepted snapshot unchanged. It uses typed
//! Stage-6 observations and commands, without converting m, V or C into N.
//! This bridge is quasistatic: it supplies no inertia, electrical circuit,
//! current spike or arbitrary RSolver task dispatch. The separate cable
//! Newmark NodalForce path retains its dynamic contract.
//!
//! Cable commands change L0 with EA/N0 fixed, change N0 with L0/EA fixed, or
//! change an existing prescribed point-anchor displacement. Cable work is
//! the change in stationary total potential along catalog-ordered changes.
//! Piezo work is midpoint(V) * delta(Q) for the reversible linear path between
//! equilibria, with Q the terminal conductor charge. Power is interval-average
//! work/dt, not an instantaneous electrical current measurement.
class RStructuralTypedFEControlDriver
{
    public:
        using Callback = std::function<RStructuralTypedControlDecision(
            const RStructuralControlObservation &)>;

        static std::unique_ptr<RStructuralTypedFEControlDriver> cable(
            const RModel &model,const RStructuralTypedControlSession &session,
            const std::vector<RStructuralTypedFEControlBinding> &bindings);
        static std::unique_ptr<RStructuralTypedFEControlDriver> piezoBeam(
            const RModel &model,const RStructuralPiezoBeamFEInput &input,
            const RStructuralTypedControlSession &session,
            const std::vector<RStructuralTypedFEControlBinding> &bindings);

        ~RStructuralTypedFEControlDriver();
        RStructuralTypedFEControlDriver(const RStructuralTypedFEControlDriver &) = delete;
        RStructuralTypedFEControlDriver &operator=(const RStructuralTypedFEControlDriver &) = delete;

        const RModel &model() const;
        const RStructuralTypedControlSession::State &state() const;
        const RStructuralTypedControlSession &session() const;
        //! Last local invocation's record; empty after construction/restore.
        const std::vector<RStructuralTypedControlChannelSample> &lastChannels() const;
        RStructuralControlObservation observation() const;
        void advance(const Callback &callback);
        //! Atomic bundle of the native binary model, typed catalog, bindings,
        //! patch inputs and accepted controller/conditioning/work state.
        void writeCheckpoint(const QString &fileName) const;
        //! Check integrity, identity and current physical inputs, then freshly
        //! certify equilibrium before returning a separately owned driver.
        static std::unique_ptr<RStructuralTypedFEControlDriver> readCheckpoint(
            const QString &fileName,const QString &expectedControllerIdentity);

    private:
        struct Impl;
        std::unique_ptr<Impl> impl;
        explicit RStructuralTypedFEControlDriver(std::unique_ptr<Impl> implementation);
};

#endif

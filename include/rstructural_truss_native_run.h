#ifndef RSTRUCTURAL_TRUSS_NATIVE_RUN_H
#define RSTRUCTURAL_TRUSS_NATIVE_RUN_H

#include <memory>

#include "raccepted_step_publication.h"
#include "rstructural_truss_dynamics_model_adapter.h"

class RModel;
class RResults;
class RStructuralNativeControl;

//! Accepted physical state for the dedicated RSolver native mixed run.
//! File transactions and notifications belong to RSolver, not this facade.
//! Its original source, target and command schedule remain frozen across runs.
class RStructuralTrussNativeRun
{
private:
    struct AcceptedState;

public:
    class PreparedStep
    {
    public:
        ~PreparedStep();
        PreparedStep(const PreparedStep &) = delete;
        PreparedStep &operator=(const PreparedStep &) = delete;

        RModel &candidateModel() noexcept;
        const RModel &candidateModel() const noexcept;
        const RAcceptedStepIdentity &identity() const noexcept;
        const RStructuralTrussDynamicsResult &result() const noexcept;
        uint completedSteps() const noexcept;
        uint nextCommandIndex() const noexcept;

    private:
        friend class RStructuralTrussNativeRun;
        PreparedStep();
        const RStructuralTrussNativeRun *owner = nullptr;
        const AcceptedState *baseState = nullptr;
        std::unique_ptr<AcceptedState> nextState;
        std::unique_ptr<RModel> model;
        std::unique_ptr<RResults> sealedResults;
        RAcceptedStepIdentity stepIdentity;
        RStructuralTrussDynamicsResult stepResult;
        uint completed = 0;
        uint cursor = 0;
        bool committed = false;
    };

    //! getRestart() selects dedicated complete-record restore; otherwise only
    //! source initial conditions supply the initial physical state.
    explicit RStructuralTrussNativeRun(const RModel &source);
    ~RStructuralTrussNativeRun();
    RStructuralTrussNativeRun(const RStructuralTrussNativeRun &) = delete;
    RStructuralTrussNativeRun &operator=(const RStructuralTrussNativeRun &) = delete;

    const RStructuralTrussDynamicsModelAdapter &adapter() const noexcept;
    const RStructuralTrussDynamicsState &acceptedState() const noexcept;
    const RStructuralNativeControlState &acceptedControlState() const noexcept;
    uint completedSteps() const noexcept;
    uint nextCommandIndex() const noexcept;
    uint targetSteps() const noexcept;
    bool hasRemainingSteps() const noexcept;
    bool isRetryable() const noexcept;

    //! Incomplete external rollback permanently closes this instance to trials.
    //! Its last accepted physical state remains available for diagnostics.
    void markNonRetryable() noexcept;
    void validateCurrentSource(const RModel &source) const;

    //! Independent full accepted baseline, including original canonical grid.
    RModel acceptedModel(const RModel &source) const;

    //! All physical checks, state/model copies and identity construction occur
    //! before external publication. Failure never advances accepted state.
    //! A known invocation boundary requests a final checkpoint notification
    //! without changing the original physical run target or source identity.
    std::unique_ptr<PreparedStep> prepare(const RModel &source,
                                         bool invocationBoundary = false) const;

    //! Validate a prepared candidate before any irreversible publication.
    //! A candidate must belong to this run and to its current accepted state.
    void validatePrepared(const PreparedStep &step) const;
    //! Same certificate for the host's assigned/staged output model.
    void validatePrepared(const PreparedStep &step,const RModel &model) const;

    //! PRECONDITION: validatePrepared(step) succeeded and no other step was
    //! committed since; the final recovery-authority publish has succeeded.
    //! Only unique_ptr/scalar promotion occurs here. The old accepted state
    //! stays owned by step until destruction outside the promotion interval.
    void commit(PreparedStep &step) noexcept;

private:
    RStructuralTrussDynamicsModelAdapter frozenAdapter;
    std::unique_ptr<RStructuralNativeControl> controlProvider;
    std::unique_ptr<AcceptedState> accepted;
    bool retryable = true;
};

#endif // RSTRUCTURAL_TRUSS_NATIVE_RUN_H

#ifndef RACCEPTED_STEP_PUBLICATION_H
#define RACCEPTED_STEP_PUBLICATION_H

#include <QtGlobal>
#include <memory>

class RModel;

struct RAcceptedStepIdentity
{
    quint64 sequence = 0;
    double timeS = 0.0;
    double dtS = 0.0;
    //! Last accepted step of this invocation (including a bounded pause),
    //! or the original physical run's final step; force a recovery checkpoint.
    //! An invocation boundary never changes the model's original run target.
    bool finalStep = false;
};

//! Prepared external publication. A recovery-authority participant is the
//! LAST publication: all physical/controller/result candidates and other
//! publications must already be ready. After its successful publish(), the
//! host performs only proven noexcept promotion and commit(), with no I/O,
//! allocation, observer, hook or cancellation check between those operations.
//! A failed publish remains abortable; rollback reports whether the previous
//! authority was restored. Rollback after a successful final publication is
//! compensation outside that abort/reopen guarantee. commit performs no I/O.
//! Participant implementations and fault hooks must not mutate the host's
//! candidate model during publish()/commit(). Tokens own staged snapshots;
//! physical actions enter through the candidate preparation API, not hooks.
class RAcceptedStepPublication
{
public:
    virtual ~RAcceptedStepPublication() = default;
    virtual void publish() = 0;
    virtual bool rollback() noexcept = 0;
    virtual void commit() noexcept = 0;
};

class RAcceptedStepParticipant
{
public:
    virtual ~RAcceptedStepParticipant() = default;
    virtual std::unique_ptr<RAcceptedStepPublication> stage(
        const RModel &candidateModel, const RAcceptedStepIdentity &identity) = 0;
};

#endif // RACCEPTED_STEP_PUBLICATION_H

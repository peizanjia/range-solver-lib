#include "rstructural_truss_native_run.h"

#include <cassert>
#include <cstring>
#include <utility>
#include <rbl_error.h>
#include <rml_model.h>

#include "rstructural_truss_dynamics_restart_adapter.h"
#include "rstructural_native_control.h"

namespace
{
[[noreturn]] void fail(const char *message)
{
    throw RError(RError::Type::InvalidInput,R_ERROR_REF,"%s",message);
}

bool samePreparedResults(const RResults &left,const RResults &right)
{
    const auto sameBits = [](double a,double b)
    { return std::memcmp(&a,&b,sizeof(double)) == 0; };
    if (left.getNNodes() != right.getNNodes() || left.getNElements() != right.getNElements() ||
        left.getNVariables() != right.getNVariables()) return false;
    for (uint i=0;i<left.getNVariables();++i)
    {
        const auto &a = left.getVariable(i);
        const auto &b = right.getVariable(i);
        if (a.getType() != b.getType() || a.getApplyType() != b.getApplyType() ||
            a.getName() != b.getName() || a.getUnits() != b.getUnits() ||
            a.getNVectors() != b.getNVectors() || a.getNValues() != b.getNValues()) return false;
        const auto &ad = a.getVariableData();
        const auto &bd = b.getVariableData();
        if (!sameBits(ad.getMinDisplayValue(),bd.getMinDisplayValue()) ||
            !sameBits(ad.getMaxDisplayValue(),bd.getMaxDisplayValue()) ||
            ad.getDisplayValueRange() != bd.getDisplayValueRange() ||
            ad.getValueRangeName() != bd.getValueRangeName() ||
            !sameBits(ad.getScale(),bd.getScale())) return false;
        for (uint component=0;component<a.getNVectors();++component)
            for (uint value=0;value<a.getNValues();++value)
                if (!sameBits(a.getValue(component,value),b.getValue(component,value))) return false;
    }
    return true;
}
}

struct RStructuralTrussNativeRun::AcceptedState
{
    RStructuralTrussDynamicsState physical;
    RStructuralNativeControlState modelControlState;
    uint completed = 0;
    uint cursor = 0;
};

RStructuralTrussNativeRun::PreparedStep::PreparedStep() = default;
RStructuralTrussNativeRun::PreparedStep::~PreparedStep() = default;

RModel &RStructuralTrussNativeRun::PreparedStep::candidateModel() noexcept
{
    return *model;
}

const RModel &RStructuralTrussNativeRun::PreparedStep::candidateModel() const noexcept
{
    return *model;
}

const RAcceptedStepIdentity &RStructuralTrussNativeRun::PreparedStep::identity() const noexcept
{
    return stepIdentity;
}

const RStructuralTrussDynamicsResult &RStructuralTrussNativeRun::PreparedStep::result() const noexcept
{
    return stepResult;
}

uint RStructuralTrussNativeRun::PreparedStep::completedSteps() const noexcept
{
    return completed;
}

uint RStructuralTrussNativeRun::PreparedStep::nextCommandIndex() const noexcept
{
    return cursor;
}

RStructuralTrussNativeRun::RStructuralTrussNativeRun(const RModel &source)
    : frozenAdapter(source)
{
    if (frozenAdapter.nativeControlSpec().source == RStructuralNativeControlSource::LocalV3)
        controlProvider = std::make_unique<RStructuralNativeControl>(frozenAdapter);
    auto initial = std::make_unique<AcceptedState>();
    if (source.getProblemSetup().getRestart())
    {
        auto restored = RStructuralTrussDynamicsRestartAdapter::restore(
            frozenAdapter,source.getProblemSetup().getStructuralTrussDynamicsRestartData(),
            source.getProblemSetup().getStructuralNativeControlState());
        initial->physical = std::move(restored.state);
        initial->modelControlState = std::move(restored.controlState);
        initial->completed = restored.completedSteps;
        initial->cursor = restored.nextCommandIndex;
    }
    else
    {
        initial->physical = frozenAdapter.initialize();
        if (controlProvider) initial->modelControlState = controlProvider->initialize();
        initial->cursor = frozenAdapter.commandCursor(0);
    }
    accepted = std::move(initial);
}

RStructuralTrussNativeRun::~RStructuralTrussNativeRun() = default;

const RStructuralTrussDynamicsModelAdapter &RStructuralTrussNativeRun::adapter() const noexcept
{
    return frozenAdapter;
}

const RStructuralTrussDynamicsState &RStructuralTrussNativeRun::acceptedState() const noexcept
{
    return accepted->physical;
}

const RStructuralNativeControlState &RStructuralTrussNativeRun::acceptedControlState() const noexcept
{
    return accepted->modelControlState;
}

uint RStructuralTrussNativeRun::completedSteps() const noexcept
{
    return accepted->completed;
}

uint RStructuralTrussNativeRun::nextCommandIndex() const noexcept
{
    return accepted->cursor;
}

uint RStructuralTrussNativeRun::targetSteps() const noexcept
{
    return frozenAdapter.targetSteps();
}

bool RStructuralTrussNativeRun::hasRemainingSteps() const noexcept
{
    return accepted->completed < frozenAdapter.targetSteps();
}

bool RStructuralTrussNativeRun::isRetryable() const noexcept
{
    return retryable;
}

void RStructuralTrussNativeRun::markNonRetryable() noexcept
{
    retryable = false;
}

void RStructuralTrussNativeRun::validateCurrentSource(const RModel &source) const
{
    frozenAdapter.validateCurrentSource(source);
}

RModel RStructuralTrussNativeRun::acceptedModel(const RModel &source) const
{
    return frozenAdapter.candidateModel(source,accepted->physical,
                                        accepted->completed,accepted->cursor,accepted->modelControlState);
}

std::unique_ptr<RStructuralTrussNativeRun::PreparedStep>
RStructuralTrussNativeRun::prepare(const RModel &source,bool invocationBoundary) const
{
    if (!retryable)
        fail("Native mixed run cannot retry after an incomplete publication rollback.");
    // Validate even an exhausted run so a caller cannot silently reuse a
    // changed source by observing only the old accepted cursor.
    validateCurrentSource(source);
    if (!hasRemainingSteps())
        fail("Native mixed run has already reached its original frozen target.");

    RStructuralTrussDynamicsModelTrial trial;
    RStructuralNativeControlState nextControlState;
    if (controlProvider)
    {
        const auto proposal = controlProvider->propose(accepted->modelControlState,accepted->physical,accepted->completed);
        validateCurrentSource(source); // Provider cannot alter the frozen source across a callback boundary.
        nextControlState = proposal.state;
        trial = frozenAdapter.trialWithCompleteInputs(accepted->physical,accepted->completed,proposal.completeInputs);
    }
    else trial = frozenAdapter.trial(accepted->physical,accepted->completed,accepted->cursor);
    std::unique_ptr<PreparedStep> prepared(new PreparedStep);
    prepared->owner = this;
    prepared->baseState = accepted.get();
    prepared->completed = accepted->completed+1;
    prepared->cursor = trial.nextCommandIndex;
    prepared->nextState = std::make_unique<AcceptedState>();
    prepared->nextState->physical = std::move(trial.state);
    prepared->nextState->modelControlState = std::move(nextControlState);
    prepared->nextState->completed = prepared->completed;
    prepared->nextState->cursor = prepared->cursor;
    prepared->stepResult = std::move(trial.result);
    prepared->model = std::make_unique<RModel>(frozenAdapter.candidateModel(
        source,prepared->nextState->physical,prepared->completed,prepared->cursor,prepared->nextState->modelControlState));
    prepared->sealedResults = std::make_unique<RResults>(*prepared->model);
    prepared->stepIdentity.sequence = prepared->completed;
    prepared->stepIdentity.timeS = prepared->nextState->physical.time;
    // Publication/recovery identity carries the exact frozen contract dt.
    // The adapter's representable grid interval is local to the physical trial.
    prepared->stepIdentity.dtS = frozenAdapter.nominalTimeStep();
    prepared->stepIdentity.finalStep = invocationBoundary ||
        prepared->completed == frozenAdapter.targetSteps();
    return prepared;
}

void RStructuralTrussNativeRun::validatePrepared(const PreparedStep &step) const
{
    if (!step.model)
        fail("Native mixed prepared step has no complete candidate model.");
    validatePrepared(step,*step.model);
}

void RStructuralTrussNativeRun::validatePrepared(const PreparedStep &step,const RModel &model) const
{
    if (!retryable)
        fail("Native mixed run cannot publish after an incomplete publication rollback.");
    if (step.owner != this || step.baseState != accepted.get() || step.committed ||
        !step.nextState || !step.model || !step.sealedResults || step.completed != accepted->completed+1 ||
        step.completed > frozenAdapter.targetSteps())
        fail("Native mixed prepared step no longer belongs to the current accepted run state.");
    validateCurrentSource(model);
    if (!samePreparedResults(model,*step.sealedResults))
        fail("Native mixed staged model results or display metadata differ from the prepared physical candidate.");
    // candidateModel is mutable for host-owned monitoring/output preparation.
    // Certify its full physical recovery record again before publication so
    // such preparation cannot accidentally replace the candidate authority.
    const auto expected = RStructuralTrussDynamicsRestartAdapter::build(
        frozenAdapter,step.nextState->physical,step.completed,step.cursor,step.nextState->modelControlState);
    const auto &record = model.getProblemSetup().getStructuralTrussDynamicsRestartData();
    if (RStructuralTrussDynamicsRestartAdapter::stateFingerprint(record,model.getProblemSetup().getStructuralNativeControlState()) !=
        RStructuralTrussDynamicsRestartAdapter::stateFingerprint(expected,step.nextState->modelControlState))
        fail("Native mixed prepared model recovery record differs from its complete physical candidate.");
    const auto &time = model.getTimeSolver();
    const auto &times = time.getTimes();
    if (!time.getEnabled() || time.getComputedTime() != step.stepIdentity.timeS ||
        time.getCurrentTimeStep() != step.completed-1 || times.size() != frozenAdapter.targetSteps() ||
        time.getInputNTimeSteps() != frozenAdapter.targetSteps() ||
        time.getInputTimeStepSize() != frozenAdapter.nominalTimeStep())
        fail("Native mixed prepared model clock differs from its complete physical candidate.");
    for (uint i=0;i<frozenAdapter.targetSteps();++i)
        if (times[i] != frozenAdapter.targetTime(i+1))
            fail("Native mixed prepared model no longer contains its original canonical run grid.");
}

void RStructuralTrussNativeRun::commit(PreparedStep &step) noexcept
{
    // Contract checks only: all throwing physical/source/model checks were
    // completed before the final recovery-authority publication.
    assert(retryable && step.owner == this && step.baseState == accepted.get());
    assert(!step.committed && step.nextState && step.model);
    static_assert(noexcept(accepted.swap(step.nextState)),
                  "Native mixed accepted promotion must not throw.");
    accepted.swap(step.nextState);
    step.committed = true;
}

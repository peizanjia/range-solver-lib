#ifndef RSTRUCTURAL_PIEZO_LAYERED_NATIVE_RUN_H
#define RSTRUCTURAL_PIEZO_LAYERED_NATIVE_RUN_H
#include "rstructural_piezo_layered_run_state.h"
#include "raccepted_step_publication.h"
#include <memory>
class RResults;
struct RStructuralPiezoLayeredControllerProposal
{
    std::vector<double> requestedVoltages;
    QByteArray nextOpaqueState,lastFrameBytes;
};
//! Value-returning callbacks are deterministic for explicit state. bindForceState
//! only selects an immutable query ledger; it never promotes state. Identities bind compatibility,
//! not authenticity of arbitrary code. Configuration bytes must exactly match source.
struct RStructuralPiezoLayeredRuntimeCallbacks
{
    std::vector<double> authoritativeFullC;
    QByteArray excitationConfiguration,controllerConfiguration,designConfiguration;
    QString sourceFingerprint,forcingIdentity,controllerIdentity,backendIdentity;
    RStructuralPiezoLayeredDynamics::BackgroundForce perturbationForce;
    //! Freeze this interval into an explicit candidate; never promote accepted state.
    std::function<QByteArray(double,uint,const QByteArray &)> prepareFE;
    //! Select a frozen ledger for read-only force queries and certification.
    std::function<void(const QByteArray &)> bindForceState;
    //! Optional interval force, distinct from the preceding accepted boundary.
    std::function<std::vector<double>(double)> intervalPerturbationForce;
    uint controlEveryFESteps=1;
    std::vector<RStructuralPiezoLayeredVoltageEnvelope> envelopes;
    QByteArray initialOpaqueControllerState;
    std::function<RStructuralPiezoLayeredControllerProposal(const RStructuralPiezoLayeredDynamicsState &,uint,uint,const QByteArray &)> propose;
    std::function<QByteArray(const RStructuralPiezoLayeredDynamicsState &,uint,uint,const std::vector<double> &,const QByteArray &)> postFE;
    std::function<void(const QByteArray &,uint,uint,uint)> validateOpaque;
    //! Certify duplicated backend histories against the actual accepted plant
    //! boundary, including previous-applied commands. Never mutates history.
    std::function<void(const RStructuralPiezoLayeredRunState &)> validateAcceptedBoundary;
};
class RStructuralPiezoLayeredNativeRun
{
public:
    class PreparedStep
    {
    public:
        ~PreparedStep();
        PreparedStep(const PreparedStep &)=delete;
        PreparedStep &operator=(const PreparedStep &)=delete;
        RModel &candidateModel() noexcept;
        const RModel &candidateModel() const noexcept;
        const RAcceptedStepIdentity &identity() const noexcept;
        const RStructuralPiezoLayeredDynamicsState &result() const noexcept;
        const RStructuralPiezoLayeredRunState &runtime() const noexcept;
        uint completedSteps() const noexcept;
    private:
        friend class RStructuralPiezoLayeredNativeRun;
        PreparedStep();
        const RStructuralPiezoLayeredNativeRun *owner=nullptr;
        const RStructuralPiezoLayeredRunState *base=nullptr;
        std::unique_ptr<RStructuralPiezoLayeredRunState> next;
        std::unique_ptr<const RStructuralPiezoLayeredRunState> sealedRuntime;
        std::unique_ptr<RModel> model;
        std::unique_ptr<RResults> sealedResults;
        RAcceptedStepIdentity stepIdentity;
        bool committed=false;
    };
    explicit RStructuralPiezoLayeredNativeRun(const RModel &source,
        RStructuralPiezoLayeredRuntimeCallbacks callbacks={},RMatrixSolverConf solverConf=RMatrixSolverConf(RMatrixSolverConf::GMRES));
    ~RStructuralPiezoLayeredNativeRun();
    RStructuralPiezoLayeredNativeRun(const RStructuralPiezoLayeredNativeRun &)=delete;
    RStructuralPiezoLayeredNativeRun &operator=(const RStructuralPiezoLayeredNativeRun &)=delete;
    const RStructuralPiezoLayeredModelAdapter &adapter() const noexcept;
    const RStructuralPiezoLayeredDynamicsState &acceptedState() const noexcept;
    const RStructuralPiezoLayeredRunState &acceptedRuntime() const noexcept;
    uint completedSteps() const noexcept;
    uint targetSteps() const noexcept;
    bool hasRemainingSteps() const noexcept;
    bool isRetryable() const noexcept;
    void markNonRetryable() noexcept;
    void validateCurrentSource(const RModel &) const;
    //! Validate a complete accepted MODEL against this frozen source/driver,
    //! including its restart codec and original canonical clock. Returns an
    //! independent state; never advances accepted physics/controller history.
    RStructuralPiezoLayeredRunState certifyAcceptedRuntime(const RModel &) const;
    //! The caller must first certify/bind the record's source authority to this
    //! frozen context. Validates all runtime/physical/backend guards without
    //! an RModel clock or source declaration; never changes accepted state.
    void certifyAcceptedRuntime(const RStructuralPiezoLayeredRunState &) const;
    RModel acceptedModel(const RModel &) const;
    //! Reconstruct a complete MODEL from independently certified restart bytes.
    //! Checks the frozen source, physical/controller boundary and original grid;
    //! leaves this run's accepted state unchanged.
    RModel modelForAcceptedRuntime(const RModel &,const RStructuralPiezoLayeredRunState &) const;
    std::unique_ptr<PreparedStep> prepare(const RModel &,bool invocationBoundary=false) const;
    void validatePrepared(const PreparedStep &) const;
    void validatePrepared(const PreparedStep &,const RModel &) const;
    void commit(PreparedStep &) noexcept;
private:
    RStructuralPiezoLayeredRuntimeCallbacks callbacks;
    RStructuralPiezoLayeredModelAdapter frozenAdapter;
    std::unique_ptr<RStructuralPiezoLayeredDynamics> driver;
    std::unique_ptr<RStructuralPiezoLayeredRunState> accepted;
    QByteArray frozenConfiguration;
    QString frozenSource;
    double initialTime=0,dt=0;
    uint target=0;
    bool retryable=true;
    double targetTime(uint) const;
    void certifyRuntime(const RStructuralPiezoLayeredRunState &) const;
    RModel candidate(const RModel &,const RStructuralPiezoLayeredRunState &) const;
};
#endif

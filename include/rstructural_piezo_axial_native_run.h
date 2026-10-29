#ifndef RSTRUCTURAL_PIEZO_AXIAL_NATIVE_RUN_H
#define RSTRUCTURAL_PIEZO_AXIAL_NATIVE_RUN_H
#include "rstructural_piezo_axial_run_state.h"
#include "raccepted_step_publication.h"
#include <functional>
#include <memory>
class RResults;
struct RStructuralPiezoAxialControllerProposal
{
    std::vector<double> requestedVoltages;
    QByteArray nextOpaqueState,lastFrameBytes;
};
//! Value-returning callbacks are deterministic for explicit state. bindForceState
//! only selects an immutable query ledger; it never promotes state. Recipe and identities bind
//! compatibility only; the application must certify recipe-to-backend semantics.
struct RStructuralPiezoAxialRuntimeCallbacks
{
    //! Actuator vectors follow adapter.maps().orderedPortIds, filtered role != Sensor.
    std::vector<double> authoritativeFullC,initialActuatorVoltages;
    QByteArray recipeConfiguration;
    QString sourceFingerprint,forcingIdentity,controllerIdentity,backendIdentity;
    std::function<std::vector<double>(double)> perturbationForce;
    //! Freeze this interval into an explicit candidate; never promote accepted state.
    std::function<QByteArray(double,uint,const QByteArray &)> prepareFE;
    //! Select a frozen ledger for read-only force queries and certification.
    std::function<void(const QByteArray &)> bindForceState;
    //! Optional interval force, distinct from the preceding accepted boundary.
    std::function<std::vector<double>(double)> intervalPerturbationForce;
    uint controlEveryFESteps=1;
    std::vector<RStructuralPiezoAxialVoltageEnvelope> envelopes;
    QByteArray initialOpaqueControllerState;
    std::function<RStructuralPiezoAxialControllerProposal(const RStructuralTrussDynamicsState &,const std::vector<RStructuralPiezoNativeReadoutSample> &,uint,uint,const QByteArray &)> propose;
    std::function<QByteArray(const RStructuralTrussDynamicsState &,const std::vector<RStructuralPiezoNativeReadoutSample> &,uint,uint,const std::vector<double> &,const QByteArray &)> postFE;
    std::function<void(const QByteArray &,uint,uint,uint)> validateOpaque;
    //! Required for controllers: certify duplicate physical/controller histories,
    //! including previousApplied, actual applied actions, and all counters.
    std::function<void(const RStructuralPiezoAxialRunState &)> validateAcceptedBoundary;
};
class RStructuralPiezoAxialNativeRun
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
        const RStructuralTrussDynamicsState &result() const noexcept;
        const RStructuralPiezoAxialRunState &runtime() const noexcept;
        uint completedSteps() const noexcept;
    private:
        friend class RStructuralPiezoAxialNativeRun;
        PreparedStep();
        const RStructuralPiezoAxialNativeRun *owner=nullptr;
        const RStructuralPiezoAxialRunState *base=nullptr;
        std::unique_ptr<RStructuralPiezoAxialRunState> next;
        std::unique_ptr<const RStructuralPiezoAxialRunState> sealedRuntime;
        std::unique_ptr<RModel> model;
        std::unique_ptr<RResults> sealedResults;
        RAcceptedStepIdentity stepIdentity;
        bool committed=false;
    };
    explicit RStructuralPiezoAxialNativeRun(const RModel &source,RStructuralPiezoAxialRuntimeCallbacks callbacks={});
    ~RStructuralPiezoAxialNativeRun();
    RStructuralPiezoAxialNativeRun(const RStructuralPiezoAxialNativeRun &)=delete;
    RStructuralPiezoAxialNativeRun &operator=(const RStructuralPiezoAxialNativeRun &)=delete;
    const RStructuralTrussDynamicsModelAdapter &adapter() const noexcept;
    const RStructuralTrussDynamicsSystem &system() const noexcept;
    const RStructuralTrussDynamicsState &acceptedState() const noexcept;
    const RStructuralPiezoAxialRunState &acceptedRuntime() const noexcept;
    uint completedSteps() const noexcept;
    uint targetSteps() const noexcept;
    bool hasRemainingSteps() const noexcept;
    bool isRetryable() const noexcept;
    void markNonRetryable() noexcept;
    void validateCurrentSource(const RModel &) const;
    //! Validate a complete accepted MODEL against this frozen source/driver,
    //! including its restart codec and original canonical clock. Returns an
    //! independent state; never advances accepted physics/controller history.
    RStructuralPiezoAxialRunState certifyAcceptedRuntime(const RModel &) const;
    //! The caller must first certify/bind the record's source authority to this
    //! frozen context. Validates all runtime/physical/backend guards without
    //! an RModel clock or source declaration; never changes accepted state.
    void certifyAcceptedRuntime(const RStructuralPiezoAxialRunState &) const;
    RModel acceptedModel(const RModel &) const;
    //! Pure reconstruction; verifies the frozen source, complete physical and
    //! controller boundary and canonical original MODEL grid before returning.
    RModel modelForAcceptedRuntime(const RModel &,const RStructuralPiezoAxialRunState &) const;
    std::unique_ptr<PreparedStep> prepare(const RModel &,bool invocationBoundary=false) const;
    void validatePrepared(const PreparedStep &) const;
    void validatePrepared(const PreparedStep &,const RModel &) const;
    void commit(PreparedStep &) noexcept;
private:
    RStructuralPiezoAxialRuntimeCallbacks callbacks;
    RStructuralTrussDynamicsModelAdapter frozenAdapter;
    RStructuralTrussDynamicsSystem actualSystem;
    std::unique_ptr<RStructuralTrussDynamics> driver;
    std::unique_ptr<RStructuralPiezoNativeReadout> readout;
    std::unique_ptr<RStructuralPiezoAxialRunState> accepted;
    RStructuralTrussDynamicsState initialPhysical;
    std::vector<QString> actuatorIds;
    std::vector<std::size_t> actuatorMembers;
    QByteArray frozenConfiguration;
    QString frozenSource;
    double initialTime=0,dt=0;
    uint target=0;
    bool retryable=true;
    double targetTime(uint) const;
    std::vector<double> totalForce(double) const;
    RStructuralTrussDynamicsInputs inputs(double,const std::vector<double> &) const;
    void certifyRuntime(const RStructuralPiezoAxialRunState &) const;
    RModel candidate(const RModel &,const RStructuralPiezoAxialRunState &) const;
};
#endif

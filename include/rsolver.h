#ifndef RSOLVER_H
#define RSOLVER_H

#include <functional>
#include <memory>
#include <vector>

#include <rml_model.h>

#include "rsolvergeneric.h"
#include "rstructural_support_response.h"
#include "rstructural_control_protocol.h"
#include "rsolvershareddata.h"
#include "raccepted_step_publication.h"
#include "rstructural_piezo_layered_native_run.h"
#include "rstructural_piezo_axial_native_run.h"

struct RParallelExecutionReport;
class RStructuralTrussNativeRun;

class RSolver
{

    public:

        using ParallelExecutionReportObserver =
            std::function<void(const RParallelExecutionReport &)>;
        using StructuralSupportResponseObserver =
            std::function<void(const std::vector<RStructuralSupportSample> &)>;
        using StructuralControlCallback =
            std::function<RStructuralControlDecision(const RStructuralControlObservation &)>;
        using StructuralControlStepObserver =
            std::function<void(const RStructuralControlStepRecord &)>;
        using ResultFrameObserver = std::function<void(const RSolverResultFrame &)>;
        using AcceptedStructuralStepObserver = std::function<void(const RModel &)>;

        enum class ExecutionPolicy
        {
            StructuralProduct,
            AllRegisteredSolvers
        };

    protected:

        //! Pointer to model object.
        RModel *pModel;
        //! Model file name.
        QString modelFileName;
        //! Convergence file name.
        QString convergenceFileName;
        //! Controls whether product capability restrictions are enforced.
        ExecutionPolicy executionPolicy;
        //! Shared variables.
        RSolverSharedData sharedData;
        //! Map of solvers.
        QMap<RProblemTypeMask,RSolverGeneric*> solvers;
        //! Map of sover type and execution count.
        QMap<RProblemType,uint> solversExecutionCount;
        //! Optional synchronous observation of the run's scoped OpenMP report.
        ParallelExecutionReportObserver parallelExecutionReportObserver;
        StructuralSupportResponseObserver structuralSupportResponseObserver;
        QString structuralControllerId;
        QString structuralControllerRevision;
        StructuralControlCallback structuralControlCallback;
        StructuralControlStepObserver structuralControlStepObserver;
        AcceptedStructuralStepObserver acceptedStructuralStepObserver;

    private:

        //! Internal initialization function.
        void _init();
        //! Preserve a failed transient Heat->Stress trial index for same-solver retry.
        bool transientHeatStressRetryPending;
        //! Dedicated accepted state survives failed trials and repeated run().
        std::unique_ptr<RStructuralTrussNativeRun> nativeMixedRun;
        std::unique_ptr<RStructuralPiezoLayeredNativeRun> nativeLayeredRun;
        RStructuralPiezoLayeredRuntimeCallbacks nativeLayeredCallbacks;
        std::unique_ptr<RStructuralPiezoAxialNativeRun> nativeAxialRun;
        RStructuralPiezoAxialRuntimeCallbacks nativeAxialCallbacks;
        std::shared_ptr<RAcceptedStepParticipant> acceptedStepParticipant;
        uint nativeMixedStepBudget = 0;
        uint nativeLayeredStepBudget = 0;
        uint nativeAxialStepBudget = 0;
        bool nativeAcceptedTraceOutput = false;
        //! Keep the exact model backup if restoring caller-owned storage fails.
        std::unique_ptr<RModel> nativeMixedRollbackModel;
        std::unique_ptr<RModel> nativeLayeredRollbackModel;
        std::unique_ptr<RModel> nativeAxialRollbackModel;
        template<class NativeRun> void runNativePiezo(NativeRun &run,uint stepBudget,
                                                     std::unique_ptr<RModel> &rollbackModel);

    public:

        //! Constructor.
        RSolver(RModel &model,
                const QString &modelFileName,
                const QString &convergenceFileName,
                ExecutionPolicy executionPolicy = ExecutionPolicy::AllRegisteredSolvers,
                RStructuralPiezoLayeredRuntimeCallbacks layeredCallbacks = {},
                RStructuralPiezoAxialRuntimeCallbacks axialCallbacks = {});

        //! Copy constructor (deleted to prevent double-free of solvers).
        RSolver(const RSolver &solver) = delete;

        //! Assignment operator (deleted to prevent double-free of solvers).
        RSolver & operator =(const RSolver &solver) = delete;

        //! Destructor.
        ~RSolver();

        //! Run solver.
        void run();

        //! Called once after the parallel context is created and before any
        //! time step or result write. The report's actualThreads is a scope
        //! probe, not a count measured separately in every kernel. An
        //! exception aborts run() before solving.
        void setParallelExecutionReportObserver(ParallelExecutionReportObserver observer);

        //! Accepted-step support samples, emitted after a successful result write.
        //! Observer exceptions are logged and cannot revoke an already committed step.
        void setStructuralSupportResponseObserver(StructuralSupportResponseObserver observer);

        //! Last accepted physical-unit support samples (a detached copy).
        //! Requires setStructuralSupportResponseObserver() before run(); an
        //! unsubscribed solve does not assemble constrained dynamic rows.
        std::vector<RStructuralSupportSample> getStructuralSupportSamples() const;

        //! Register a pure, replayable in-process controller before run().
        //! It is sampled once from each accepted state, never during assembly
        //! or an iterative linear solve. The ID must match model control setup.
        void setStructuralControlCallback(const QString &controllerId,
                                          const QString &controllerRevision,
                                          StructuralControlCallback callback);
        //! Notification after successful result write. An observer failure
        //! stops subsequent steps and is reported by run(), without undoing
        //! the already committed model artifact.
        void setStructuralControlStepObserver(StructuralControlStepObserver observer);
        std::vector<RStructuralControlChannelSample> getStructuralControlSamples() const;

        //! Observe actual model-result frames after the native file write completes.
        //! Callback is synchronous on the solver owner thread; observer failures
        //! are diagnostic only and cannot change solver success.
        void setResultFrameObserver(ResultFrameObserver observer);

        //! Observe each fully accepted structural-dynamics step, independent
        //! of result display frequency. Observer failures stop subsequent steps.
        void setAcceptedStructuralStepObserver(AcceptedStructuralStepObserver observer);

        //! Dedicated native mixed, layered or axial mode. Stage before output publication; its
        //! recovery-authority publish is last, followed by noexcept commits.
        void setAcceptedStepParticipant(std::shared_ptr<RAcceptedStepParticipant> participant);
        //! Complete accepted trace replaces repeated display MODEL writes. The
        //! physical grid and saved display frequency remain unchanged.
        void setNativeAcceptedTraceOutput(bool enabled);

        //! Maximum accepted FE steps per run() invocation; zero is unlimited.
        //! The original physical target and source identity remain unchanged.
        //! The last budget step requests a final checkpoint publication.
        void setNativeMixedStepBudget(uint stepBudget);
        void setNativeLayeredStepBudget(uint stepBudget);
        void setNativeAxialStepBudget(uint stepBudget);

        //! False after incomplete external/model rollback; further native
        //! trials are rejected and the last accepted state is preserved.
        bool isNativeMixedRetryable() const noexcept;
        bool isNativeLayeredRetryable() const noexcept;
        bool isNativeAxialRetryable() const noexcept;

    protected:

        //! Dedicated native dispatch; never enters the generic stress driver.
        void runNativeMixed();
        void runNativeLayered();
        void runNativeAxial();

        //! Run single solver.
        void runSingle();

        //! Run single solver for given problem task.
        //! Return task convergence status (true = converged).
        bool runProblemTask(const RProblemTaskItem &problemTaskItem, uint taskIteration);

};

#endif // RSOLVER_H

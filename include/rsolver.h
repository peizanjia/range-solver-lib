#ifndef RSOLVER_H
#define RSOLVER_H

#include <functional>
#include <vector>

#include <rml_model.h>

#include "rsolvergeneric.h"
#include "rstructural_support_response.h"
#include "rsolvershareddata.h"

struct RParallelExecutionReport;

class RSolver
{

    public:

        using ParallelExecutionReportObserver =
            std::function<void(const RParallelExecutionReport &)>;
        using StructuralSupportResponseObserver =
            std::function<void(const std::vector<RStructuralSupportSample> &)>;

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

    private:

        //! Internal initialization function.
        void _init();

    public:

        //! Constructor.
        RSolver(RModel &model,
                const QString &modelFileName,
                const QString &convergenceFileName,
                ExecutionPolicy executionPolicy = ExecutionPolicy::AllRegisteredSolvers);

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

    protected:

        //! Run single solver.
        void runSingle();

        //! Run single solver for given problem task.
        //! Return task convergence status (true = converged).
        bool runProblemTask(const RProblemTaskItem &problemTaskItem, uint taskIteration);

};

#endif // RSOLVER_H

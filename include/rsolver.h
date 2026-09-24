#ifndef RSOLVER_H
#define RSOLVER_H

#include <functional>

#include <rml_model.h>

#include "rsolvergeneric.h"
#include "rsolvershareddata.h"

struct RParallelExecutionReport;

class RSolver
{

    public:

        using ParallelExecutionReportObserver =
            std::function<void(const RParallelExecutionReport &)>;

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

    protected:

        //! Run single solver.
        void runSingle();

        //! Run single solver for given problem task.
        //! Return task convergence status (true = converged).
        bool runProblemTask(const RProblemTaskItem &problemTaskItem, uint taskIteration);

};

#endif // RSOLVER_H

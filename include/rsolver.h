#ifndef RSOLVER_H
#define RSOLVER_H

#include <rml_model.h>

#include "rsolvergeneric.h"
#include "rsolvershareddata.h"

class RSolver
{

    public:

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

    protected:

        //! Run single solver.
        void runSingle();

        //! Run single solver for given problem task.
        //! Return task convergence status (true = converged).
        bool runProblemTask(const RProblemTaskItem &problemTaskItem, uint taskIteration);

};

#endif // RSOLVER_H

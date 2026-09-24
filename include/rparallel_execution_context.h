#ifndef RPARALLEL_EXECUTION_CONTEXT_H
#define RPARALLEL_EXECUTION_CONTEXT_H

#include <rml_parallel_execution_config.h>

#include <QString>

#include <thread>

//! Effective OpenMP policy for one solver run on the creating host thread.
struct RParallelExecutionReport
{
    RParallelExecutionMode mode = RParallelExecutionMode::Auto;
    uint requestedThreads = 0; //!< Stored model request; zero means runtime selected.
    uint appliedThreads = 1;   //!< Value passed to the OpenMP run-time ICV.
    uint actualThreads = 1;    //!< Team size measured when entering the scope.
    uint availableThreads = 1; //!< Processors visible to OpenMP.
    QString backend = QStringLiteral("OpenMP");
    bool deterministicReductionRequested = true;
    bool deterministicReduction = true; //!< True only for a single worker.
    bool fallback = false;
    QString diagnostic;
};

//! Scoped, host-thread-local OpenMP configuration for an entire solver run.
//!
//! The context changes the OpenMP ICV of its creating host thread. A run must
//! start and end on that thread. Workers spawned by that run inherit the ICV.
//! Parallel kernels remain responsible for reproducible reductions and for
//! using the same context during assembly, matrix solution and time stepping.
class RParallelExecutionContext
{
    public:
        explicit RParallelExecutionContext(const RParallelExecutionConfig &config);
        ~RParallelExecutionContext();

        RParallelExecutionContext(const RParallelExecutionContext &) = delete;
        RParallelExecutionContext &operator=(const RParallelExecutionContext &) = delete;
        RParallelExecutionContext(RParallelExecutionContext &&) = delete;
        RParallelExecutionContext &operator=(RParallelExecutionContext &&) = delete;

        const RParallelExecutionReport &getReport() const;

    private:
        int previousMaxThreads = 1;
        int previousDynamic = 0;
        std::thread::id ownerThread;
        RParallelExecutionReport report;
};

#endif // RPARALLEL_EXECUTION_CONTEXT_H

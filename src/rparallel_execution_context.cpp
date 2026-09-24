#include "rparallel_execution_context.h"

#include <QStringList>
#include <QtGlobal>

#include <omp.h>

#include <algorithm>
#include <limits>

RParallelExecutionContext::RParallelExecutionContext(const RParallelExecutionConfig &config)
    : ownerThread(std::this_thread::get_id())
{
    config.validate();

    this->previousMaxThreads = std::max(1,omp_get_max_threads());
    this->previousDynamic = omp_get_dynamic();
    this->report.mode = config.getMode();
    this->report.requestedThreads = config.getRequestedThreads();
    this->report.availableThreads = uint(std::max(1,omp_get_num_procs()));
    this->report.deterministicReductionRequested = config.getDeterministicReduction();

    const bool nested = omp_in_parallel() != 0;
    const uint threadLimit = uint(std::max(1,omp_get_thread_limit()));
    const uint integerLimit = uint(std::numeric_limits<int>::max());
    const uint maxRuntimeThreads = std::min({threadLimit,integerLimit,
                                             this->report.availableThreads});
    uint requestedScopeThreads = 1;
    switch (config.getMode())
    {
        case RParallelExecutionMode::Auto:
            requestedScopeThreads = uint(this->previousMaxThreads);
            break;
        case RParallelExecutionMode::SerialDeterministic:
            requestedScopeThreads = 1;
            break;
        case RParallelExecutionMode::FixedThreads:
            requestedScopeThreads = config.getRequestedThreads();
            break;
    }
    this->report.appliedThreads = std::min(requestedScopeThreads,maxRuntimeThreads);

    // Dynamic team sizing would silently violate FixedThreads and make the
    // reported run configuration unstable. Restore both ICVs on scope exit.
    omp_set_dynamic(0);
    omp_set_num_threads(int(this->report.appliedThreads));

    int measuredThreads = 1;
#pragma omp parallel reduction(max:measuredThreads)
    {
        measuredThreads = omp_get_num_threads();
    }
    this->report.actualThreads = uint(std::max(1,measuredThreads));
    this->report.deterministicReduction = this->report.actualThreads == 1;

    QStringList reasons;
    if (this->report.appliedThreads != requestedScopeThreads)
    {
        this->report.fallback = true;
        reasons.append(QStringLiteral("Requested thread count exceeds available processors or the OpenMP runtime limit"));
    }
    if (this->report.actualThreads != this->report.appliedThreads)
    {
        this->report.fallback = true;
        reasons.append(QStringLiteral("OpenMP created %1 of %2 requested workers")
                       .arg(this->report.actualThreads)
                       .arg(this->report.appliedThreads));
    }
    if (nested && this->report.appliedThreads > 1)
    {
        reasons.append(QStringLiteral("Scope entered inside an existing OpenMP parallel region"));
    }
    if (this->report.deterministicReductionRequested &&
        !this->report.deterministicReduction)
    {
        reasons.append(QStringLiteral("Parallel reduction order must be made deterministic by each solver kernel"));
    }
    this->report.diagnostic = reasons.join(QStringLiteral("; "));
}

RParallelExecutionContext::~RParallelExecutionContext()
{
    Q_ASSERT(std::this_thread::get_id() == this->ownerThread);
    omp_set_num_threads(this->previousMaxThreads);
    omp_set_dynamic(this->previousDynamic);
}

const RParallelExecutionReport &RParallelExecutionContext::getReport() const
{
    return this->report;
}

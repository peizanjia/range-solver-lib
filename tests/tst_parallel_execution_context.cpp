#include <QtTest>

#include <rbl_error.h>

#include <omp.h>

#include <algorithm>

#include "rparallel_execution_context.h"

class TestParallelExecutionContext : public QObject
{
    Q_OBJECT

    private slots:
        void serialScopeRestoresRuntimeSettings();
        void fixedThreadsReportsActualBackendAndFallback();
        void autoScopeUsesRuntimeSelection();
        void invalidConfigurationDoesNotChangeRuntime();
};

void TestParallelExecutionContext::serialScopeRestoresRuntimeSettings()
{
    const int originalMaxThreads = omp_get_max_threads();
    const int originalDynamic = omp_get_dynamic();

    RParallelExecutionConfig config;
    config.setMode(RParallelExecutionMode::SerialDeterministic);
    {
        RParallelExecutionContext context(config);
        const RParallelExecutionReport &report = context.getReport();
        QCOMPARE(report.mode,RParallelExecutionMode::SerialDeterministic);
        QCOMPARE(report.requestedThreads,0u);
        QCOMPARE(report.appliedThreads,1u);
        QCOMPARE(report.actualThreads,1u);
        QCOMPARE(report.backend,QStringLiteral("OpenMP"));
        QVERIFY(report.deterministicReductionRequested);
        QVERIFY(report.deterministicReduction);
        QVERIFY(!report.fallback);
        QCOMPARE(omp_get_max_threads(),1);
        QCOMPARE(omp_get_dynamic(),0);
    }
    QCOMPARE(omp_get_max_threads(),originalMaxThreads);
    QCOMPARE(omp_get_dynamic(),originalDynamic);
}

void TestParallelExecutionContext::fixedThreadsReportsActualBackendAndFallback()
{
    RParallelExecutionConfig config;
    config.setMode(RParallelExecutionMode::FixedThreads);
    config.setRequestedThreads(2);
    config.setDeterministicReduction(false);
    {
        RParallelExecutionContext context(config);
        const RParallelExecutionReport &report = context.getReport();
        QCOMPARE(report.requestedThreads,2u);
        QCOMPARE(report.backend,QStringLiteral("OpenMP"));
        QVERIFY(report.appliedThreads >= 1u);
        QVERIFY(report.appliedThreads <= 2u);
        QVERIFY(report.actualThreads >= 1u);
        QVERIFY(report.actualThreads <= report.appliedThreads);
        QCOMPARE(report.deterministicReduction,report.actualThreads == 1u);
        QCOMPARE(report.fallback,report.actualThreads != 2u);
        if (report.fallback)
        {
            QVERIFY(!report.diagnostic.isEmpty());
        }
    }

    // The requested count is intentionally above the visible processor count.
    // It must be reported as a fallback without creating excessive workers.
    const uint processors = uint(qMax(1,omp_get_num_procs()));
    config.setRequestedThreads(processors + 1u);
    RParallelExecutionContext limited(config);
    const RParallelExecutionReport &report = limited.getReport();
    QCOMPARE(report.requestedThreads,processors + 1u);
    QCOMPARE(report.appliedThreads,uint(std::min(omp_get_num_procs(),
                                                 omp_get_thread_limit())));
    QVERIFY(report.fallback);
    QVERIFY(!report.diagnostic.isEmpty());
}

void TestParallelExecutionContext::autoScopeUsesRuntimeSelection()
{
    const int originalMaxThreads = omp_get_max_threads();
    RParallelExecutionConfig config;
    RParallelExecutionContext context(config);
    const RParallelExecutionReport &report = context.getReport();
    QCOMPARE(report.mode,RParallelExecutionMode::Auto);
    QCOMPARE(report.requestedThreads,0u);
    QCOMPARE(report.appliedThreads,uint(std::min({originalMaxThreads,
                                                   omp_get_num_procs(),
                                                   omp_get_thread_limit()})));
    QVERIFY(report.actualThreads >= 1u);
    QVERIFY(report.actualThreads <= report.appliedThreads);
    QCOMPARE(report.deterministicReduction,report.actualThreads == 1u);
    if (report.actualThreads > 1)
    {
        QVERIFY(report.diagnostic.contains(QStringLiteral("reduction order")));
    }
}

void TestParallelExecutionContext::invalidConfigurationDoesNotChangeRuntime()
{
    const int originalMaxThreads = omp_get_max_threads();
    const int originalDynamic = omp_get_dynamic();
    RParallelExecutionConfig config;
    config.setMode(RParallelExecutionMode::FixedThreads);
    QVERIFY_EXCEPTION_THROWN(static_cast<void>(RParallelExecutionContext{config}),RError);
    QCOMPARE(omp_get_max_threads(),originalMaxThreads);
    QCOMPARE(omp_get_dynamic(),originalDynamic);
}

QTEST_APPLESS_MAIN(TestParallelExecutionContext)

#include "tst_parallel_execution_context.moc"

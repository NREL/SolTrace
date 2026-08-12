#pragma once

#include <QPromise>

#include "job_run_common.h"

/// Simulation backend used by the threaded runner.
enum class ThreadRunnerBackend {
    Native,
    Embree,
    Optix,
};

/// Configuration for one threaded simulation run.
struct ThreadRunnerConfig {
    uint32_t            thread_count;
    ThreadRunnerBackend backend = ThreadRunnerBackend::Native;
};

/// Execute a simulation through the selected runner and report through promise.
void execute_thread_runner(QPromise<SimResult>&      promise,
                           SimDataPtr                data,
                           ThreadRunnerConfig const& config);

#pragma once

#include <QPromise>

#include "job_run_common.h"

enum class ThreadRunnerBackend {
    Native,
    Embree,
};

struct ThreadRunnerConfig {
    uint32_t            thread_count;
    ThreadRunnerBackend backend = ThreadRunnerBackend::Native;
};

void execute_thread_runner(QPromise<SimResult>&      promise,
                           SimDataPtr                data,
                           ThreadRunnerConfig const& config);

#pragma once

#include <QPromise>

#include "job_run_common.h"

#include "job_run_thread.h"

void execute_process_runner(QPromise<SimResult>&      promise,
                            SimDataPtr                data,
                            ThreadRunnerConfig const& config);


/// Checks if this process is configured to be a worker only process.
/// If so, executes, and then EXITS the process!
void check_if_process_worker(int argc, char* argv[]);

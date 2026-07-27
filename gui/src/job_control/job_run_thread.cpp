#include "job_run_thread.h"

#include "native_runner/native_runner.hpp"
#include "simulation_result.hpp"
#include "simulation_runner.hpp"
#include "utilities/math_utility.h"

#ifdef SOLTRACE_HAS_EMBREE_RUNNER
#include "embree_runner/embree_runner.hpp"
#endif

#ifdef SOLTRACE_HAS_OPTIX_RUNNER
#include "optix_runner/optix_runner.hpp"
#endif

#include <QtConcurrentRun>
#include <QtGlobal>

#include <memory>

// Qconcurrent will auto call start and finish on the promise

#define SOLTRACE_SECTION(FUNC, VALUE, TEXT)                                    \
    qDebug() << Q_FUNC_INFO << TEXT;                                           \
    promise.setProgressValueAndText(VALUE, TEXT);                              \
    promise.suspendIfRequested();                                              \
    if (promise.isCanceled()) { return; }                                      \
    result = current_runner->FUNC;                                             \
    if (result == SolTrace::Runner::RunnerStatus::ERROR) {                     \
        qWarning() << Q_FUNC_INFO << "failed at" << TEXT << (int)result;       \
        promise.emplaceResult(QString(TEXT " failed"));                        \
        return;                                                                \
    }

#define SECTION(VALUE, TEXT)                                                   \
    qDebug() << Q_FUNC_INFO << TEXT;                                           \
    promise.setProgressValueAndText(VALUE, TEXT);                              \
    promise.suspendIfRequested();                                              \
    if (promise.isCanceled()) {                                                \
        promise.emplaceResult("Cancelled at " TEXT);                           \
        return;                                                                \
    }

/// Function to wrap up run into a thread. We will be checking progress
/// OUTSIDE this function
static SolTrace::Runner::RunnerStatus
execute_solve_with(SolTrace::Runner::SimulationRunner* ptr) {
    qDebug() << "execute runner";
    return ptr->run_simulation();
}

static std::unique_ptr<SolTrace::Runner::SimulationRunner>
make_runner(ThreadRunnerBackend backend) {
#ifdef SOLTRACE_HAS_EMBREE_RUNNER
    if (backend == ThreadRunnerBackend::Embree) {
        return std::make_unique<SolTrace::EmbreeRunner::EmbreeRunner>();
    }
#endif

#ifdef SOLTRACE_HAS_OPTIX_RUNNER
    if (backend == ThreadRunnerBackend::Optix) {
        return std::make_unique<OptixRunner>();
    }
#endif

    return std::make_unique<SolTrace::NativeRunner::NativeRunner>();
}


/// Add results from a simulation. We don't use a constructor here as
/// it does not play well with the progress update and cancel concept

void construct_result(QPromise<SimResult>&                promise,
                      SimDataPtr                          exported_source,
                      SolTrace::Result::SimulationResult& result) {

    SECTION(90, "Building lookup tables");

    db::SimulationResultConversion opts {
        .result = result,
        .data   = *(exported_source->data),
        .map    = exported_source->element_map,
    };

    auto destination = db::SimulationResult::convert(opts);

    destination->database = std::move(exported_source->source_database);

    SECTION(100, "Done");

    promise.emplaceResult(std::move(destination));
}


void execute_thread_runner(QPromise<SimResult>&      promise,
                           SimDataPtr                data,
                           ThreadRunnerConfig const& config) {
    try {
        auto start_instant = std::chrono::high_resolution_clock::now();

        promise.setProgressRange(0, 100);

        auto current_runner = make_runner(config.backend);


        size_t thread_count = config.thread_count;

        qDebug() << "Starting count" << thread_count;

        if (thread_count == 0) {
            thread_count =
                std::max<size_t>(std::thread::hardware_concurrency(), 1);
        }

        // appears to be something odd in the sim code if the number of threads
        // is at or above the ray count. clamp.
        if (thread_count >= data->data->get_number_of_rays()) {
            thread_count = 1;
        }

        if (auto ptr = dynamic_cast<SolTrace::NativeRunner::NativeRunner*>(
                current_runner.get());
            ptr) {
            ptr->set_number_of_threads(thread_count);
            // ptr->disable_stages();

            qDebug() << "Native parameters" << thread_count;
        }

        SolTrace::Runner::RunnerStatus result;

        qDebug() << "Starting simulation with"
                 << data->data->get_number_of_rays()
                 << data->data->get_max_number_rays_traced();

        SOLTRACE_SECTION(initialize(), 0, "Starting simulation");


        SOLTRACE_SECTION(
            setup_simulation(data->data.get()), 0, "Setup simulation");

        qDebug() << Q_FUNC_INFO << "setup complete";

#if defined(Q_OS_WASM) && !defined(__EMSCRIPTEN_PTHREADS__)
        promise.setProgressValueAndText(10, "Running...");
        auto run_result = current_runner->run_simulation();
        promise.setProgressValueAndText(90, "Simulation complete");
#else
        // run simulation will indeed run, but we need to stuff it into another
        // thread so we can poll status.

        auto run_future =
            QtConcurrent::run(execute_solve_with, current_runner.get());

        double last_progress = -1;

        while (true) {
            if (run_future.isFinished()) { break; }

            // Normally polling would be The Wrong Thing, but the progress stuff
            // requires active checking

            // Docs say sleep is in nsecs
            // we want milliseconds
            QThread::msleep(500);

            double progress = -1;

            // If cancelled, we set the flag, wait, and bail
            if (promise.isCanceled()) {
                promise.setProgressValueAndText(100, "Cancelling...");

                // User cancelled, set flag
                current_runner->cancel_simulation();

                // Wait for completion
                run_future.result();

                // bail
                promise.emplaceResult(QStringLiteral("Cancelled"));
                qDebug() << Q_FUNC_INFO << "cancelled";
                return;
            }

            // Not cancelled, check and see how things are going
            current_runner->status_simulation(&progress);

            qDebug() << "raw" << progress;

            // assuming progress is 0-1
            progress = std::clamp(progress, 0.0, 1.0);


            if (last_progress != progress) {
                // we have reserved progress points 10 to 90 for sim run

                promise.setProgressValueAndText(std::lerp(10, 90, progress),
                                                "Running...");
                last_progress = progress;
                qDebug() << "sim progress " << std::lerp(10, 90, progress);
            }
        }

        auto run_result = run_future.result();
#endif

        auto end_instant = std::chrono::high_resolution_clock::now();

        qDebug() << Q_FUNC_INFO << "Run complete in: "
                 << std::chrono::duration<double>(end_instant - start_instant)
                        .count();

        switch (run_result) {
        case SolTrace::Runner::RunnerStatus::CANCEL:
            promise.setProgressValueAndText(100, "Cancelled");
            promise.emplaceResult(QStringLiteral("Cancelled"));
            return;
        case SolTrace::Runner::RunnerStatus::ERROR:
            promise.setProgressValueAndText(100, "Simulation failed");
            promise.emplaceResult(QStringLiteral("The simulation failed."));
            return;
        case SolTrace::Runner::RunnerStatus::RUNNING:
            // we really shouldnt get here
            qWarning("Sim result declares running, but was finished??");
            break;
        case SolTrace::Runner::RunnerStatus::SUCCESS: break;
        case SolTrace::Runner::RunnerStatus::TIMEOUT:
            promise.setProgressValueAndText(100, "Simulation timed out");
            promise.emplaceResult(QStringLiteral("The simulation timed out."));
            return;
        default:
            promise.setProgressValueAndText(100, "Simulation failed");
            promise.emplaceResult(
                QStringLiteral("The simulation failed for an unknown reason."));
            return;
        }

        qDebug() << Q_FUNC_INFO << "Build result database";


        SolTrace::Result::SimulationResult soltrace_result;

        SOLTRACE_SECTION(
            report_simulation(&soltrace_result, 100), 90, "Report simulation");

        qDebug() << Q_FUNC_INFO
                 << "Launched:" << current_runner->get_number_rays_launched();
        qDebug() << Q_FUNC_INFO
                 << "Traced:" << current_runner->get_number_rays_traced();

        return construct_result(promise, data, soltrace_result);

    } catch (std::exception& e) {
        promise.emplaceResult(QString(e.what()));
        return;
    }
}

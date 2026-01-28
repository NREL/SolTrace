#include <gtest/gtest.h>

#include <thread_manager.hpp>
#include <trace_logger.hpp>

#include <chrono>
#include <future>
#include <thread>

using SolTrace::NativeRunner::make_thread_manager;
using SolTrace::NativeRunner::make_trace_logger;

using SolTrace::NativeRunner::thread_manager_ptr;
using SolTrace::NativeRunner::ThreadManager;
using SolTrace::NativeRunner::trace_logger_ptr;

// "Runs" nsteps step where each step takes at least tlen milliseconds
// checking the manager to see if it needs to terminate early after
// each step
ThreadManager::ThreadStatus timed_task(thread_manager_ptr manager,
                                       trace_logger_ptr logger,
                                       unsigned thread_id,
                                       uint_fast64_t nsteps,
                                       uint_fast64_t tlen)
{
    uint_fast64_t count = 0;
    double thread_progress = 0.0;

    while (count < nsteps)
    {
        std::this_thread::sleep_for(std::chrono::milliseconds(tlen));
        ++count;
        thread_progress = static_cast<double>(count) / nsteps;
        manager->progress_update(thread_id, count);

        if (manager->terminate(thread_id))
        {
            return ThreadManager::ThreadStatus::CANCEL;
        }
    }

    return ThreadManager::ThreadStatus::SUCCESS;
}

// Waits for tlen milliseconds and returns an error
ThreadManager::ThreadStatus error_task(thread_manager_ptr manager,
                                       trace_logger_ptr logger,
                                       unsigned thread_id,
                                       uint_fast64_t tlen)
{
    std::this_thread::sleep_for(std::chrono::milliseconds(tlen));
    return ThreadManager::ThreadStatus::ERROR;
}

TEST(ThreadManager, CleanExit)
{
    const unsigned NTHREADS = 2;
    const uint_fast64_t NSTEPS = 5;
    const uint_fast64_t T_MS = 100;

    auto logger = make_trace_logger();
    auto manager = make_thread_manager(logger);

    manager->initialize();
    for (auto thid = 0; thid < NTHREADS; ++thid)
    {
        auto f = std::async(std::launch::async,
                            timed_task,
                            manager,
                            logger,
                            thid,
                            NSTEPS,
                            T_MS);
        manager->manage(thid, std::move(f));
    }

    std::this_thread::sleep_for(std::chrono::milliseconds(2 * T_MS));
    double progress = 0.0;
    ThreadManager::ThreadStatus sts = manager->status(&progress);
    EXPECT_EQ(sts, ThreadManager::ThreadStatus::RUNNING);
    EXPECT_GT(progress, 0.0);

    sts = manager->monitor_until_completion();
    EXPECT_EQ(sts, ThreadManager::ThreadStatus::SUCCESS);
}

TEST(ThreadManager, ErrorExit)
{
    const unsigned NTHREADS = 2;
    const uint_fast64_t NSTEPS = 5;
    const uint_fast64_t T_MS = 100;

    auto logger = make_trace_logger();
    auto manager = make_thread_manager(logger);

    manager->initialize();
    for (auto thid = 0; thid < NTHREADS; ++thid)
    {
        ThreadManager::future f;
        if (thid < NTHREADS - 1)
        {
            f = std::async(std::launch::async,
                           timed_task,
                           manager,
                           logger,
                           thid,
                           NSTEPS,
                           T_MS);
            manager->manage(thid, std::move(f));
        }
        else
        {
            f = std::async(std::launch::async,
                           error_task,
                           manager,
                           logger,
                           thid,
                           T_MS);
            manager->manage(thid, std::move(f));
        }
    }

    ThreadManager::ThreadStatus sts = manager->monitor_until_completion();
    EXPECT_EQ(sts, ThreadManager::ThreadStatus::ERROR);
}

TEST(ThreadManager, CancelExit)
{
    const unsigned NTHREADS = 2;
    const uint_fast64_t NSTEPS = 5;
    const uint_fast64_t T_MS = 100;

    auto logger = make_trace_logger();
    auto manager = make_thread_manager(logger);

    manager->initialize();
    for (auto thid = 0; thid < NTHREADS; ++thid)
    {
        auto f = std::async(std::launch::async,
                            timed_task,
                            manager,
                            logger,
                            thid,
                            NSTEPS,
                            T_MS);
        manager->manage(thid, std::move(f));
    }

    std::this_thread::sleep_for(std::chrono::milliseconds(2 * T_MS));

    manager->cancel();
    ThreadManager::ThreadStatus sts = manager->monitor_until_completion();
    EXPECT_EQ(sts, ThreadManager::ThreadStatus::CANCEL);
}

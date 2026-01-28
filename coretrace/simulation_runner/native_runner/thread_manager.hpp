#ifndef SOLTRACE_THREAD_MANAGER_H
#define SOLTRACE_THREAD_MANAGER_H

#include <future>
#include <map>
#include <memory>
#include <mutex>
#include <ostream>
#include <string>

#include <simulation_runner.hpp>

#include "trace_logger.hpp"

namespace SolTrace::NativeRunner
{

    class ThreadManager
    {
    public:
        using ThreadStatus = SolTrace::Runner::RunnerStatus;
        using future = std::future<SolTrace::Runner::RunnerStatus>;
        
        ThreadManager(trace_logger_ptr logger);
        ~ThreadManager();

        // Functions for whoever "hired" the manager. Pattern for use is
        // 1. Call initialize()
        // 2. Make calls to std::async() for each task/thread handing
        // the returned future to manage along with a unique id
        // 3. Call monitor_until_completion() to wait for everything
        // to finish.
        // NOTE: Tasks should return type of ThreadManager::ThreadStatus
        
        // Must call prior to manage and not while running!!
        void initialize();
        // Returns a thread id
        unsigned int manage(unsigned int id, future f);
        ThreadStatus monitor_until_completion();

        // For the worker threads that are being managed
        void progress_update(unsigned int id, double progress);
        bool terminate(unsigned int id);

        // For general public use -- thread safe and can be called whenever
        ThreadStatus status(double *progress = nullptr) const;
        void cancel() const;

    private:
        // unsigned int next_id;

        mutable std::mutex state_mutex;
        mutable ThreadStatus state;

        mutable std::mutex progress_mutex;
        std::map<unsigned int, double> progress;

        std::map<unsigned int, future> threads;

        trace_logger_ptr logger;
    };

    using thread_manager_ptr = std::shared_ptr<ThreadManager>;
    template <typename... Args>
    inline auto make_thread_manager(Args &&...args)
    {
        return std::make_shared<ThreadManager>(std::forward<Args>(args)...);
    }

} // namespace SolTrace::NativeRunner

#endif

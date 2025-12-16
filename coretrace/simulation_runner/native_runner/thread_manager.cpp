#include "thread_manager.hpp"

#include <chrono>
#include <future>
#include <sstream>
#include <thread>

namespace SolTrace::NativeRunner
{

    ThreadManager::ThreadManager()
    {
        this->initialize();
        return;
    }

    ThreadManager::~ThreadManager()
    {
        this->threads.clear();
        this->progress.clear();
        return;
    }

    unsigned int ThreadManager::manage(unsigned int id, future f)
    {
        assert(this->threads.find(id) == this->threads.end());

        this->threads[id] = std::move(f);
        this->progress[id] = 0.0;
        return id;
    }

    ThreadManager::ThreadStatus ThreadManager::monitor_until_completion()
    {
        ThreadStatus sts = ThreadStatus::SUCCESS;
        bool canceled = false;

        while (!this->threads.empty())
        {
            auto iter = this->threads.begin();
            while (iter != this->threads.end())
            {
                if (iter->second.wait_for(std::chrono::seconds(0)) == std::future_status::ready)
                {
                    ThreadStatus thread_status = iter->second.get();
                    if (thread_status == ThreadStatus::SUCCESS)
                    {
                        ; // Intentional no-op
                    }
                    else if (thread_status == ThreadStatus::ERROR)
                    {
                        // Something went wrong -- shut everything down
                        this->cancel();
                        // Return ERROR status
                        sts = ThreadStatus::ERROR;
                        canceled = true;
                    }
                    else if (thread_status == ThreadStatus::CANCEL)
                    {
                        // Return CANCEL if threads were canceled
                        // by an external call to cancel()
                        sts = canceled ? sts : thread_status;
                    }
                    else
                    {
                        // Thread terminated with something other than
                        // SUCCESS, ERROR, or CANCEL. This is unexpected.
                        std::stringstream ss;
                        ss << "Thread " << iter->first
                           << " returned status "
                           << status_string(thread_status)
                           << ". This is an error.";
                        this->error_log(ss.str());

                        // Unexpected behavior so we terminate
                        this->cancel();
                        sts = ThreadStatus::ERROR;
                        canceled = true;
                    }
                    // Task is complete so we stop tracking it.
                    // Also increments to the next spot
                    iter = this->threads.erase(iter);
                }
                else
                {
                    ++iter;
                }
            }

            std::this_thread::sleep_for(std::chrono::milliseconds(10));
        }

        return sts;
    }

    void ThreadManager::error_log(const std::string &msg)
    {
        std::lock_guard<std::mutex> lk(this->message_mutex);
        this->messages.push_back(msg);
        return;
    }

    void ThreadManager::progress_update(unsigned int id, double prog)
    {
        // std::cout << "Update from thread " << id << " reporting "
        //           << prog << " complete." << std::endl;
        std::lock_guard<std::mutex> lk(this->progress_mutex);
        this->progress[id] = prog;
        return;
    }

    bool ThreadManager::terminate(unsigned int id)
    {
        std::lock_guard<std::mutex> lk(this->state_mutex);
        return this->state != ThreadStatus::RUNNING;
    }

    ThreadManager::ThreadStatus ThreadManager::status(double *progress) const
    {
        ThreadStatus sts = ThreadStatus::ERROR;
        // Create isolated scope for lock guard
        {
            std::lock_guard<std::mutex> lk(this->state_mutex);
            sts = this->state;
        }

        if (progress != nullptr)
        {
            int k = 0;
            double avg = 0.0;
            std::lock_guard<std::mutex> lk(this->progress_mutex);
            for (auto it = this->progress.cbegin();
                 it != this->progress.cend();
                 ++it)
            {
                // std::cout << "Thread " << k
                //           << " progress " << it->second
                //           << std::endl;
                ++k;
                avg += (it->second - avg) / k;
            }
            *progress = avg;
        }

        return sts;
    }

    void ThreadManager::cancel() const
    {
        std::lock_guard<std::mutex> lk(this->state_mutex);
        this->state = ThreadStatus::CANCEL;
        return;
    }

    void ThreadManager::print_log(std::ostream &os) const
    {
        std::lock_guard<std::mutex> lk(this->message_mutex);
        for (auto iter = this->messages.cbegin();
             iter != this->messages.cend();
             ++iter)
        {
            os << *iter;
        }
        return;
    }

    void ThreadManager::initialize()
    {
        // Must call prior to manage and not while running!!
        // this->next_id = 0;
        {
            std::lock_guard<std::mutex> lk(this->state_mutex);
            this->state = ThreadStatus::RUNNING;
        }
        {
            std::lock_guard<std::mutex> lk(this->progress_mutex);
            this->progress.clear();
            this->threads.clear();
        }
        {
            std::lock_guard<std::mutex> lk(this->message_mutex);
            this->messages.clear();
        }
        return;
    }

} // namespace SolTrace::NativeRunner

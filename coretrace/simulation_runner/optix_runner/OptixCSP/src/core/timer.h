#ifndef TIMER_H
#define TIMER_H

//#include <ratio>
#include <chrono>

namespace OptixCSP {

    /// Chrono wrappers for the high-resolution timer.
    class Timer {
    public:
        Timer() : m_running(false), m_total(0) {}

        /// Start the timer.
        void start() {
            m_running = true;
            m_start = std::chrono::high_resolution_clock::now();
        }

        /// Stop the timer.
        /// The current duration from the last time the timer was stopped is added to the accumulated time.
        /// Call reset() to zero out the accumulated time. 
        void stop() {
            if (m_running) {
                m_running = false;
                m_total += std::chrono::high_resolution_clock::now() - m_start;
            }
        }

        /// Reset the total accumulated time (when repeating multiple start/stop).
        void reset() { m_total = std::chrono::duration<double>(0); }

        /// Return the time in milliseconds.
        /// If the timer was started but not stopped, this function returns the intermediate value and the timer keeps running.
        /// Otherwise, it returns the total accumulated time when the timer was last stopped.
        unsigned long long get_time_milli_sec() const {
            if (m_running) {
                auto now = std::chrono::high_resolution_clock::now();
                return std::chrono::duration_cast<std::chrono::milliseconds>(now - m_start).count();
            }
            return std::chrono::duration_cast<std::chrono::milliseconds>(m_total).count();
        }

        /// Return the time in seconds.
        /// If the timer was started but not stopped, this function returns the intermediate value and the timer keeps
        /// running. Otherwise, it returns the timer value when it was stopped.
        double get_time_sec() const {
            if (m_running) {
                std::chrono::duration<double> int_time = std::chrono::high_resolution_clock::now() - m_start;
                return int_time.count();
            }
            return m_total.count();
        }

        /// Get the last timer value, in seconds.
        double operator()() const { return get_time_sec(); }

    private:
        bool m_running;
        std::chrono::high_resolution_clock::time_point m_start;
        std::chrono::duration<double> m_total;
        // TODO: add counter for number of times started/stopped
    };

}
#endif

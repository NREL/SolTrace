#ifndef SOLTRACE_SIMULATION_RUNNER_H
#define SOLTRACE_SIMULATION_RUNNER_H

#include <map>
#include <memory>

#include "simulation_data.hpp"
#include "result_spec.hpp"

namespace SolTrace::Runner {

enum class RunnerStatus {
    CANCEL,
    ERROR,
    RUNNING,
    SUCCESS,
    TIMEOUT,
    UNKNOWN,
};

const std::map<RunnerStatus, std::string> STATUS_TO_STR {
    { RunnerStatus::CANCEL, "CANCEL" },   { RunnerStatus::ERROR, "ERROR" },
    { RunnerStatus::RUNNING, "RUNNING" }, { RunnerStatus::SUCCESS, "SUCCESS" },
    { RunnerStatus::TIMEOUT, "TIMEOUT" }, { RunnerStatus::UNKNOWN, "UNKNOWN" }
};

inline const std::string& status_string(const RunnerStatus sts) {
    auto item = STATUS_TO_STR.find(sts);
    if (item != STATUS_TO_STR.cend()) {
        return item->second;
    } else {
        return STATUS_TO_STR.find(RunnerStatus::UNKNOWN)->second;
    }
}

class SimulationRunner {
public:
    SimulationRunner() { };
    virtual ~SimulationRunner() { };

    // Disable copy constructor
    SimulationRunner(const SimulationRunner&) = delete;
    // Disable move constructor
    SimulationRunner(SimulationRunner&&) = delete;
    // Disable assignment operators
    SimulationRunner& operator=(const SimulationRunner&) = delete;
    SimulationRunner& operator=(SimulationRunner&&)      = delete;

    virtual RunnerStatus initialize() = 0;
    virtual RunnerStatus
    setup_simulation(const SolTrace::Data::SimulationData* data,
                     const SolTrace::Result::ResultSpec&   spec
                         = SolTrace::Result::RayHistorySpec{}) = 0;
    // TODO: Determine what can be "updated", that is changed
    virtual RunnerStatus
    update_simulation(const SolTrace::Data::SimulationData* data)      = 0;
    virtual RunnerStatus run_simulation()                              = 0;
    virtual RunnerStatus status_simulation(double* progress = nullptr) = 0;
    virtual RunnerStatus cancel_simulation()                           = 0;

    /// Creates and returns a result object of the type configured in
    /// setup_simulation(). Must be called after setup_simulation().
    /// Returns nullptr if the configured ResultType is not supported.
    virtual std::unique_ptr<SolTrace::Result::SimulationResult> create_result() = 0;

    virtual RunnerStatus
    report_simulation(SolTrace::Result::SimulationResult* result,
                      int                                 level_spec = 0) = 0;

    virtual uint_fast64_t get_number_rays_launched() const = 0;
    virtual uint_fast64_t get_number_rays_traced() const   = 0;

private:
};

} // namespace SolTrace::Runner

#endif

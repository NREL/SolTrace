#ifndef SOLTRACE_SIMULATION_RUNNER_H
#define SOLTRACE_SIMULATION_RUNNER_H

#include "simulation_data.hpp"
#include "simulation_result.hpp"

namespace SolTrace::Runner
{

    enum class RunnerStatus
    {
        CANCEL,
        ERROR,
        RUNNING,
        SUCCESS,
        TIMEOUT,
    };

    class SimulationRunner
    {
    public:
        SimulationRunner() {};
        virtual ~SimulationRunner() {};

        // Disable copy constructor
        SimulationRunner(const SimulationRunner &) = delete;
        // Disable move constructor
        SimulationRunner(SimulationRunner &&) = delete;
        // Disable assignment operators
        SimulationRunner &operator=(const SimulationRunner &) = delete;
        SimulationRunner &operator=(SimulationRunner &&) = delete;

        virtual RunnerStatus initialize() = 0;
        virtual RunnerStatus setup_simulation(const SolTrace::Data::SimulationData *data) = 0;
        // TODO: Determine what can be "updated", that is changed
        virtual RunnerStatus update_simulation(const SolTrace::Data::SimulationData *data) = 0;
        virtual RunnerStatus run_simulation() = 0;
        virtual RunnerStatus status_simulation(double *progress = nullptr) = 0;
        virtual RunnerStatus cancel_simulation() = 0;
        virtual RunnerStatus report_simulation(SolTrace::Result::SimulationResult *result,
                                               int level_spec) = 0;

    private:
    };

} // namespace SolTrace::Runner

#endif

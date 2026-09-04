#ifndef SOLTRACE_EMBREE_RUNNER_H
#define SOLTRACE_EMBREE_RUNNER_H

#include <embree4/rtcore.h>

#include <simulation_data.hpp>
#include <simulation_runner.hpp>

#include <native_runner.hpp>
#include <native_runner_types.hpp>

class grouped_results_EmbreeRunner_helper;

namespace SolTrace::EmbreeRunner
{
    using SolTrace::Runner::RunnerStatus;
    using SolTrace::Runner::SimulationRunner;

    // NativeRunner types that are used here that we want to make visible
    // through the EmbreeRunner namespace
    using SolTrace::NativeRunner::TRayData;
    using SolTrace::NativeRunner::TSystem;
    using SolTrace::NativeRunner::tstage_ptr;
    using SolTrace::NativeRunner::telement_ptr;

    class EmbreeRunner : public SolTrace::NativeRunner::NativeRunner
    {
    public:
        EmbreeRunner();
        virtual ~EmbreeRunner();

        EmbreeRunner(const EmbreeRunner &) = delete;
        EmbreeRunner(EmbreeRunner &&) = delete;

        virtual RunnerStatus setup_simulation(const SolTrace::Data::SimulationData *data) override;
        virtual RunnerStatus run_simulation() override;
        virtual RunnerStatus update_simulation(const SolTrace::Data::SimulationData *data) override;

        // TODO: Do we want loud errors when a user calls these?
        void disable_power_tower() = delete;
        void enable_power_tower() = delete;
        void disable_point_focus() = delete;
        void enable_point_focus() = delete;

    private:

        RTCDevice embree_device;
        RTCScene embree_scene;

        void clean_embree();

        // could use FRIEND_TEST macro, however to avoid linking gtest to prod, forward declare test class and make it a friend
        friend class ::grouped_results_EmbreeRunner_helper;
    };

} // namespace SolTrace::EmbreeRunner

#endif


#include "embree_runner.hpp"

#include "trace_embree.hpp"

#include <native_runner.hpp>
#include <native_runner_types.hpp>
#include <simulation_data_export.hpp>
#include <simulation_runner.hpp>

namespace SolTrace::EmbreeRunner
{
    using SolTrace::Runner::RunnerStatus;

    using SolTrace::NativeRunner::telement_ptr;
    using SolTrace::NativeRunner::TRayData;
    using SolTrace::NativeRunner::tstage_ptr;
    using SolTrace::NativeRunner::TSystem;

    using SolTrace::Result::SimulationResult;

    EmbreeRunner::EmbreeRunner() : NativeRunner(),
                                   embree_device(nullptr),
                                   embree_scene(nullptr)
    {
        return;
    }

    EmbreeRunner::~EmbreeRunner()
    {
        this->clean_embree();
        return;
    }

    RunnerStatus EmbreeRunner::setup_simulation(const SimulationData *data)
    {

        RunnerStatus sts = NativeRunner::setup_simulation(data);

        make_embree_scene(this->my_logger,
                          &this->tsys,
                          this->embree_device,
                          this->embree_scene,
                          this->number_of_threads);

        return sts;
    }

    RunnerStatus EmbreeRunner::update_simulation(const SimulationData *data)
    {
        // TODO: Do a more efficient implementation of this?
        this->clean_embree();
        NativeRunner::update_simulation(data);
        return RunnerStatus::SUCCESS;
    }

    RunnerStatus EmbreeRunner::run_simulation()
    {
        this->set_seeds();

        RunnerStatus sts = trace_embree(
            this->my_manager,
            this->my_logger,
            &this->tsys,
            // this->tsys.seed,
            this->seeds,
            this->number_of_threads,
            this->tsys.sim_raycount,
            this->tsys.sim_raymax,
            this->tsys.sim_errors_sunshape,
            this->tsys.sim_errors_optical,
            this->embree_scene);

        return sts;
    }

    void EmbreeRunner::clean_embree()
    {
        // Clean embree
        if (this->embree_scene != nullptr)
        {
            rtcReleaseScene(embree_scene);
            embree_scene = nullptr;
        }

        if (this->embree_device != nullptr)
        {
            rtcReleaseDevice(embree_device);
            embree_device = nullptr;
        }
        return;
    }

} // namespace SolTrace::EmbreeRunner

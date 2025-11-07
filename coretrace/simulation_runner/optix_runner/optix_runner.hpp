#ifndef SOLTRACE_OPTIX_RUNNER_H
#define SOLTRACE_OPTIX_RUNNER_H

#include "simulation_data.hpp"
#include "simulation_result.hpp"
#include "simulation_runner.hpp"
#include "core/soltrace_system.h"

// using SolTrace::Runner::RunnerStatus;

class OptixRunner : public SolTrace::Runner::SimulationRunner
{
public:
    OptixRunner();
    ~OptixRunner();

    virtual SolTrace::Runner::RunnerStatus initialize() override;
    virtual SolTrace::Runner::RunnerStatus setup_simulation(
        const SolTrace::Data::SimulationData *data) override;
    virtual SolTrace::Runner::RunnerStatus update_simulation(
        const SolTrace::Data::SimulationData *data) override;
    virtual SolTrace::Runner::RunnerStatus run_simulation() override;
    virtual SolTrace::Runner::RunnerStatus status_simulation(double *progress = nullptr) override;
    virtual SolTrace::Runner::RunnerStatus cancel_simulation() override;
    virtual SolTrace::Runner::RunnerStatus report_simulation(
        SolTrace::Result::SimulationResult *result,
        int level_spec) override;

    // Runner options
    // void disable_sun_shape_errors() { this->include_sun_shape_errors = false; }
    // void enable_sun_shape_errors() { this->include_sun_shape_errors = true; }
    // void disable_errors() { this->include_errors = false; }
    // void enable_errors() { this->include_errors = true; }

    // Runner accessors
    // const TSystem *get_system() const { return &this->tsys; }

private:
    OptixCSP::SolTraceSystem m_sys;

    const SolTrace::Data::SimulationData *m_simdata;
    SolTrace::Runner::RunnerStatus setup_parameters(
        const SolTrace::Data::SimulationData *data);
    SolTrace::Runner::RunnerStatus setup_sun(
        const SolTrace::Data::SimulationData *data);
    SolTrace::Runner::RunnerStatus setup_elements(
        const SolTrace::Data::SimulationData *data);

    // helper function, convert Vector3d to Optix::Vec3d
    OptixCSP::Vec3d ToVec3d(SolTrace::Data::Vector3d v);
};

#endif

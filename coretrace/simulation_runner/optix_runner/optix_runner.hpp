#ifndef SOLTRACE_OPTIX_RUNNER_H
#define SOLTRACE_OPTIX_RUNNER_H

#include <set>
#include <vector>

#include "simulation_data.hpp"
#include "simulation_result.hpp"
#include "simulation_runner.hpp"
#include "core/soltrace_system.h"
#include "core/timer.h"

// using SolTrace::Runner::RunnerStatus;

class OptixRunner : public SolTrace::Runner::SimulationRunner
{
public:
    OptixRunner();
    ~OptixRunner() override = default;

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

    SolTrace::Runner::RunnerStatus run_simulation_core();
    SolTrace::Runner::RunnerStatus get_hp_output(std::vector<float4>& hp_vec,
        std::vector<uint_fast64_t>& raynumber_vec, std::vector<int32_t>& element_id_vec);

    double get_sun_plane_area() const { return m_sys.get_sun_plane_area(); }

    uint_fast64_t get_N_sun_rays() const { return m_sys.get_N_sun_rays(); }
    inline uint_fast64_t get_number_rays_launched() const override {return get_N_sun_rays(); }
    inline uint_fast64_t get_number_rays_traced() const override {return m_sys.get_N_hit_rays(); }

    uint64_t get_N_run_iterations() const;

    void print_timing() const;

    void set_verbose(bool verbose);

    // Set the number of rays to launch for a trace in each optixLaunch call.
    // WARNING: The runner is forced to use this batch size regardless of available GPU memory!!!!
    // Setting a large batch size can cause device out of memory errors or degraded GPU performance.
    // Setting a small batch size can cause long run times. Care is required when using this function.
    void set_batch_size(uint_fast64_t batch_size);
    uint_fast64_t get_batch_size() const;

    /// Set the maximum ray interaction depth. Must be called before initialize().
    /// Depth is clamped to [2, 255] with a warning if either bound is exceeded. Defaults to DEFAULT_MAX_TRACE_DEPTH.
    void set_max_ray_depth(uint_fast64_t depth);
    uint8_t get_max_ray_depth() const { return m_sys.get_max_ray_depth(); }

    /// Returns the number of rays terminated by reaching max_depth without being absorbed.
    /// Valid after run_simulation() completes; resets to 0 at the start of each run.
    uint_fast64_t get_N_depth_exceeded_rays() const { return m_sys.get_N_depth_exceeded_rays(); }

    /// Enable or disable trimming of excess rays at the end of run() so that
    /// exactly the requested number of hit rays is returned.  Enabled by default.
    void set_trim_excess_rays(bool trim);
    bool get_trim_excess_rays() const;

    // Runner options
    // void disable_sun_shape_errors() { this->include_sun_shape_errors = false; }
    // void enable_sun_shape_errors() { this->include_sun_shape_errors = true; }
    // void disable_errors() { this->include_errors = false; }
    // void enable_errors() { this->include_errors = true; }

    // Runner accessors
    OptixCSP::SolTraceSystem *get_optix_system() { return &this->m_sys; }

    // group functions
    std::vector<std::set<int32_t>> &get_groups() { return m_groups; }
    void set_groups(const std::vector<std::set<int32_t>>& groups) { m_groups = groups; }
    void set_groups(const std::vector<std::set<uint_fast64_t>>& groups) 
    { 
        m_groups.clear();
        m_groups.reserve(groups.size());
        for (size_t i = 0; i < groups.size(); ++i)
        {
            std::set<int32_t> group_i;
            for (auto id : groups[i])
            {
                group_i.insert(static_cast<int32_t>(id));
            }
            m_groups.push_back(std::move(group_i));
        }
    }
    size_t get_num_groups() const { return m_groups.size(); }
    int32_t get_group(int32_t element_id);

private:
    OptixCSP::SolTraceSystem m_sys;

    std::vector<std::set<int32_t>> m_groups;

    const SolTrace::Data::SimulationData *m_simdata;
    SolTrace::Runner::RunnerStatus setup_parameters(
        const SolTrace::Data::SimulationData *data);
    SolTrace::Runner::RunnerStatus setup_sun(
        const SolTrace::Data::SimulationData *data);
    SolTrace::Runner::RunnerStatus setup_elements(
        const SolTrace::Data::SimulationData *data);

    // helper function, convert Vector3d to Optix::Vec3d
    OptixCSP::Vec3d ToVec3d(glm::dvec3 v);
    OptixCSP::Matrix33d ToMatrix33d(const glm::dmat3& mat);
    // helper function, convert SolTrace::Data::DistributionType to Optix::OpticalDistribution
    OptixCSP::OpticalDistribution to_optical_distribution(SolTrace::Data::DistributionType dt);
    
    // timers for report simulation
    OptixCSP::Timer m_timer_report;
    OptixCSP::Timer m_timer_get_output;
    OptixCSP::Timer m_timer_report_loop;
};

#endif

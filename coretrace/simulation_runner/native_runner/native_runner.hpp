#ifndef SOLTRACE_NATIVE_RUNNER_H
#define SOLTRACE_NATIVE_RUNNER_H

#include <future>
#include <map>
#include <vector>
#include <set>

#include "native_runner_types.hpp"
#include "simulation_runner.hpp"
#include "simulation_result.hpp"
#include "thread_manager.hpp"

class grouped_results_NativeRunner_helper;

namespace SolTrace::NativeRunner
{

    using SolTrace::Runner::RunnerStatus;
    using SolTrace::Runner::SimulationRunner;

    class NativeRunner : public SimulationRunner
    {
    public:
        NativeRunner();
        ~NativeRunner();

        // Disable copy constructor
        NativeRunner(const NativeRunner &) = delete;
        // Disable move constructor
        NativeRunner(NativeRunner &&) = delete;
        // Disable assignment operators
        NativeRunner &operator=(const NativeRunner &) = delete;
        NativeRunner &operator=(NativeRunner &&) = delete;

        virtual RunnerStatus initialize() override;
        virtual RunnerStatus setup_simulation(const SolTrace::Data::SimulationData *data) override;
        virtual RunnerStatus update_simulation(const SolTrace::Data::SimulationData *data) override;
        virtual RunnerStatus run_simulation() override;
        virtual RunnerStatus status_simulation(double *progress = nullptr) override;
        virtual RunnerStatus cancel_simulation() override;
        virtual RunnerStatus report_simulation(SolTrace::Result::SimulationResult *result,
                                               int level_spec) override;

        virtual uint_fast64_t get_number_rays_launched() const override
        {
            return tsys.SunRayCount;
        }
        virtual uint_fast64_t get_number_rays_traced() const override
        {
            // TODO: This could be wrong if we hit max number of rays before getting this many hits.
            // At the moment max number of rays is ignored though...
            // return tsys.sim_raycount;
            return tsys.SunRayCount > 0 ? tsys.sim_raycount : 0;
        }

        // Runner options
        void disable_power_tower() { this->as_power_tower = false; }
        void enable_power_tower() { this->as_power_tower = true; }
        void disable_point_focus() { this->tsys.sim_dynamic_group = false; }
        void enable_point_focus() { this->tsys.sim_dynamic_group = true; }

        // Only applies BEFORE setup_simulation()
        void disable_stages() { this->use_stages = false; }
        void enable_stages() { this->use_stages = true; }

        void set_newton_tolerance(double tol)
        {
            this->eparams.newton_tolerance = tol;
            return;
        }

        void set_newton_max_iters(uint_fast64_t max_iters)
        {
            this->eparams.newton_max_iters = max_iters;
            return;
        }

        void set_number_of_threads(uint_fast64_t nthr)
        {
            this->number_of_threads = nthr;
            this->seeds.clear();
            return;
        }

        void set_number_of_threads(uint_fast64_t nthr,
                                   const std::vector<unsigned int> &seeds)
        {
            if (nthr == seeds.size())
            {
                this->number_of_threads = nthr;
                this->seeds = seeds;
            }
            else
            {
                throw std::invalid_argument("Number of seeds must equal number of threads.");
            }
            return;
        }

        void print_log(std::ostream &os)
        {
            this->my_logger->print_log(os);
            return;
        }

        // Accessors
        int_fast64_t get_number_stages() const
        {
            return this->get_system()->StageList.size();
        }
        int_fast64_t get_number_elements() const
        {
            int_fast64_t nelems = 0;
            for (auto stage : this->get_system()->StageList)
            {
                nelems += stage->ElementList.size();
            }
            return nelems;
        }

        const TSystem *get_system() const
        {
            return &this->tsys;
        }

        // Helper functions
        RunnerStatus setup_parameters(const SolTrace::Data::SimulationData *data);
        RunnerStatus setup_sun(const SolTrace::Data::SimulationData *data);
        RunnerStatus setup_elements(const SolTrace::Data::SimulationData *data);

        // group functions
        std::vector<std::set<uint_fast64_t>> &get_groups() { return m_groups; }
        void set_groups(const std::vector<std::set<uint_fast64_t>>& groups) { m_groups = groups; }
        size_t get_num_groups() const { return m_groups.size(); }
        int32_t get_group(uint_fast64_t element_id);

    protected:
        // Use power tower speed ups
        bool as_power_tower;

        // Group elements in stages
        bool use_stages = true;

        // Number of threads to use when tracing
        uint_fast64_t number_of_threads;
        std::vector<unsigned int> seeds;
        // std::map<unsigned int, std::future<RunnerStatus> > procs;

        ElementParameters eparams;

        trace_logger_ptr my_logger;
        thread_manager_ptr my_manager;
        TSystem tsys;

        std::vector<std::set<uint_fast64_t>> m_groups;

        bool set_aperture_planes(TSystem *tsys);
        bool set_aperture_planes(tstage_ptr stage);
        bool aperture_plane(telement_ptr Element);

        void set_seeds();

        void check_supported_optical_distribution(
            SolTrace::Data::DistributionType dt);
        void check_supported_options(telement_ptr telem);

    private:
            // could use FRIEND_TEST macro, however to avoid linking gtest to prod, forward declare test class and make it a friend
            friend class ::grouped_results_NativeRunner_helper;
    };

} // namespace SolTrace::NativeRunner

#endif

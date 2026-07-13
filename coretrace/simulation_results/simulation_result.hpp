#ifndef SOLTRACE_SIMULATION_RESULT_H
#define SOLTRACE_SIMULATION_RESULT_H

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include <glm/vec3.hpp>

#include "element.hpp"

// RayHistoryResult headers
#include "records.hpp"

namespace SolTrace::Result
{
    enum class ResultType
    {
        RAY_HISTORY,   ///< Full per-ray interaction history (current default behavior)
        ELEMENT_STATS, ///< Aggregated per-element statistics only
        RAW_HITS,      ///< Raw compacted hit buffer, no aggregation
    };

    /**
     * Abstract base for all result types produced by a SimulationRunner.
     * Runners inspect get_result_type() at setup_simulation() time to configure
     * internal data structures, and validate the result at report_simulation() time.
     * Sun-plane metadata is stored here so all result types carry it without
     * requiring a down-cast.
     */
    class SimulationResult
    {
    public:
        virtual ~SimulationResult() = default;

        SimulationResult(const SimulationResult &) = delete;
        SimulationResult &operator=(const SimulationResult &) = delete;

        virtual ResultType get_result_type() const = 0;

        // Sun metadata — common to all result types
        void set_sun_ray_count(uint_fast64_t ray_count) { sun_ray_count = ray_count; }
        uint_fast64_t get_sun_ray_count() const { return sun_ray_count; }
        void set_sun_dimensions(double width, double height)
        {
            sun_width = width;
            sun_height = height;
        }
        void get_sun_dimensions(double &width, double &height) const
        {
            width = sun_width;
            height = sun_height;
        }
        void set_sun_A_box(double A) { A_sun_box = A; }
        double get_sun_A_box() const { return A_sun_box; }

    protected:
        SimulationResult() = default;

    private:
        uint_fast64_t sun_ray_count = 0;
        double sun_width = 0;
        double sun_height = 0;
        double A_sun_box = 0;
    };

    using RayRecordContainer = typename std::vector<ray_record_ptr>;
    using ElementRecordContainer = typename std::map<SolTrace::Data::element_id,
                                                     element_record_ptr>;

    class RayHistoryResult : public SimulationResult
    {
    public:
        RayHistoryResult();
        ~RayHistoryResult() override;

        ResultType get_result_type() const override { return ResultType::RAY_HISTORY; }

        // Functions for getting and analyzing results
        uint_fast64_t get_number_of_records() const
        {
            return this->ray_history.size();
        }

        const element_record_ptr get_element_record(SolTrace::Data::element_id elid) const;

        RayRecordContainer::const_iterator get_ray_record_iterator() const
        {
            return ray_history.cbegin();
        }
        bool is_at_end(RayRecordContainer::const_iterator citer) const
        {
            return citer == this->ray_history.cend();
        }
        ElementRecordContainer::const_iterator get_element_record_iterator() const
        {
            return element_view.cbegin();
        }
        bool is_at_end(ElementRecordContainer::const_iterator citer) const
        {
            return citer == this->element_view.cend();
        }

        // Functions for building up results (used by Runners)
        void add_ray_record(ray_record_ptr);

        // Functions for file IO
        void write_csv_file(std::string csv_name, int precision = 12) const;
        void write_csv_file(const char *csv_name, int precision = 12) const;

        // Legacy stuff -- TODO:
        // void results_to_legacy_csv(std::string csv_name,
        //                            SimulationData *data);

        // Operator overloads
        const ray_record_ptr &operator[](int_fast64_t idx) const;
        friend std::ostream &operator<<(std::ostream &os,
                                        const RayHistoryResult &simres);

    private:
        RayRecordContainer ray_history;
        ElementRecordContainer element_view;

        void add_element_view(const ray_record_ptr rp);
    };

} // namespace SolTrace::Result

#endif

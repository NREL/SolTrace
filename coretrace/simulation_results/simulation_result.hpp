#ifndef SOLTRACE_SIMULATION_RESULT_H
#define SOLTRACE_SIMULATION_RESULT_H

#include <cstdint>
#include <map>
#include <string>
#include <vector>

#include <glm/vec3.hpp>

#include "element.hpp"

// SimulationResult headers
#include "records.hpp"

namespace SolTrace::Result
{
using RayRecordContainer = typename std::vector<ray_record_ptr>;
using ElementRecordContainer =
    typename std::map<SolTrace::Data::element_id, element_record_ptr>;

class SimulationResult
{
public:
    SimulationResult();
    virtual ~SimulationResult();

    // Functions for getting and analyzing results
    uint_fast64_t get_number_of_records() const
    { return this->ray_history.size(); }

    const element_record_ptr
    get_element_record(SolTrace::Data::element_id elid) const;

    RayRecordContainer::const_iterator get_ray_record_iterator() const
    { return ray_history.cbegin(); }
    bool is_at_end(RayRecordContainer::const_iterator citer) const
    { return citer == this->ray_history.cend(); }
    ElementRecordContainer::const_iterator get_element_record_iterator() const
    { return element_view.cbegin(); }
    bool is_at_end(ElementRecordContainer::const_iterator citer) const
    { return citer == this->element_view.cend(); }

    // Functions for building up results (used by Runners)
    void add_ray_record(ray_record_ptr);

    // Functions for file IO
    void write_csv_file(std::string csv_name, int precision = 12) const;
    void write_csv_file(const char* csv_name, int precision = 12) const;

    // Legacy stuff -- TODO:
    // void results_to_legacy_csv(std::string csv_name,
    //                            SimulationData *data);

    // Operator overloads
    const ray_record_ptr& operator[](int_fast64_t idx) const;
    friend std::ostream&  operator<<(std::ostream&           os,
                                     const SimulationResult& simres);

    // Sun results
    // void set_sun_ray_count(uint_fast64_t ray_count) { this->sun_ray_count =
    // ray_count; }
    uint_fast64_t get_sun_ray_count() const { return this->sun_ray_count; }
    // void set_sun_dimensions(double width, double height) { this->sun_width =
    // width; this->sun_height = height; }
    void get_sun_dimensions(double& width, double& height) const
    {
        width  = this->sun_width;
        height = this->sun_height;
    }
    // void set_sun_A_box(double A) { this->A_sun_box = A; }
    double get_sun_A_box() const { return this->A_sun_box; }
    double get_ray_area_weight() const { return this->ray_area_weight; }

    void set_sun_sampling_stats(double        width,
                                double        height,
                                uint_fast64_t sun_ray_count);
    void set_sun_sampling_stats(double sun_A_box, uint_fast64_t sun_ray_count);

private:
    RayRecordContainer     ray_history;
    ElementRecordContainer element_view;

    void add_element_view(const ray_record_ptr rp);

    // Sun results
    uint_fast64_t sun_ray_count   = 0;
    double        sun_width       = -1.0;
    double        sun_height      = -1.0;
    double        A_sun_box       = -1.0;
    double        ray_area_weight = -1.0;
};

} // namespace SolTrace::Result

#endif

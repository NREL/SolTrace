/**
 * @file simulation_data.hpp
 * @brief Main simulation data container and management
 *
 * Defines the SimulationData class which serves as the main container
 * for all optical elements, stages, and simulation parameters.
 * Manages the complete optical system definition and provides
 * interfaces for ray tracing execution.
 *
 * @defgroup simulation Simulation Management
 * @{
 */

#ifndef SOLTRACE_SIMULATION_DATA_H
#define SOLTRACE_SIMULATION_DATA_H

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <set>
#include <vector>

#include "container.hpp"
#include "element.hpp"
#include "ray_source.hpp"
#include "simulation_parameters.hpp"

namespace SolTrace::Data {

class SimulationData
{
    friend void load_json_file(SimulationData& sd, std::string filename, std::string* upgrade_log);
    friend void load_json_cstr(SimulationData& sd, const char* json_str, std::string* upgrade_log);
    friend void json_to_simulation_data(SimulationData& sd, nlohmann::ordered_json& root, std::string* upgrade_log);

public:
    SimulationData();
    virtual ~SimulationData();

    /// @brief enforce that all user fields are set and coordinate
    ///        retations are calculated
    // for added robustness and user flexibility consilidate enforcing
    // that elements are ready for the runner to one function.
    // adds robustness because calling from Native/OptixRunner::setup_elements
    // is actually when the fields need to be set by.
    // the current implementation: 
    // - misses enforcing in SimulationData::replace_element
    // - misses enforcing data->get_number_of_elements() <= 0 in
    //   OptixRunner::setup_elements
    // - requires addtional checks that optical property sets exists when 
    //   implying they are unnecessary
    // adds flexibility to users by making SimulationData more freeform
    // before being sent to the runner
    void enforce_elements_ready();

    /// @brief Add the given RaySource to the simulation data
    /// @param src RaySource to add
    /// @return unique identifier for the RaySource
    ray_source_id add_ray_source(ray_source_ptr src)
    {
        return this->my_sources.add_item(src);
    }

    /// @brief Remove the RaySource corresponding to the unique identifier `id`
    /// @param id unique identifier of the RaySource to remove
    /// @return true if the element was removed, false otherwise
    auto remove_ray_source(ray_source_id id)
    {
        return this->my_sources.remove_item(id);
    }

    /// @brief Get the RaySource corresponding to the unique identifier `id`
    /// @param id unique identifier of the RaySource to get
    /// @return requested RaySource if found, nullptr otherwise
    ray_source_ptr get_ray_source(ray_source_id id) const
    {
        return this->my_sources.get_item(id);
    }

    ray_source_ptr get_ray_source() const
    {
        return this->my_sources.get_const_iterator()->second;
    }

    /// @brief Replace the RaySource with id `id` with the RaySource `src`
    /// @param id unique identifier of the RaySource to replace
    /// @param src RaySource to use in replacement
    /// @return true if the RaySource was replaced, false otherwise
    bool replace_ray_source(ray_source_id id, ray_source_ptr src)
    {
        return this->my_sources.replace_item(id, src);
    }

    uint_fast64_t get_number_of_ray_sources() const
    {
        return this->my_sources.get_number_of_items();
    }

    /// @brief Get an iterator that can be used to access all ray sources owned
    ///        by this SimulationData object.
    /// @return iterator
    RaySourceContainer::iterator get_ray_source_iterator()
    {
        return this->my_sources.get_iterator();
    }

    /// @brief Tests whether the given iterator is at the end
    /// @param iter iterator to test
    /// @return true if at end, false otherwise
    bool is_ray_source_at_end(RaySourceContainer::iterator it)
    {
        return this->my_sources.is_at_end(it);
    }

    /// @brief Add the given Element to the SimulationData. When adding a
    ///        CompositeElement, all subelements must already be present.
    ///        Any subelements added to the CompositeElement after calling
    ///        `add_element` will not be used in any ray tracing.
    /// @param el Element to add
    /// @return unique identifier for the given object
    element_id add_element(element_ptr el);
    inline element_id add_stage(element_ptr el)
    {
        return add_element(el);
    }

    /// @brief Remove the Element corresponding to the unique identifier `id`
    /// @param id unique identifier of element to remove
    /// @return number of elements removed
    uint_fast64_t remove_element(element_id id);

    /// @brief Get the Element corresponding to the unique identifier `id`
    /// @param id unique identifier of element to get
    /// @return requested Element (pointer) if found, nullptr otherwise
    element_ptr get_element(element_id id) const;

    /// @brief Replace the Element with id `id` with the Element `el`
    /// @param id unique identifier of element to replace
    /// @param el Element to use in replacement
    /// @return true if Element was replaced, false otherwise
    bool replace_element(element_id id, element_ptr el);

    void clear_elements()
    {
        this->my_elements.clear();
        this->my_groups.clear();
        this->number_of_elements = 0;
    }

    /// @brief Gives the number of elements owned by the SimulationData.
    ///        CompositeElements do not count toward this number.
    /// @return Number of elements owned by the SimulationData object
    uint_fast64_t get_number_of_elements() const
    {
        // return this->my_elements.size();
        return this->number_of_elements;
    }

    /// @brief Get an iterator that can be used to access all elements owned
    ///        by this SimulationData object.
    /// @return iterator
    ElementContainer::iterator get_iterator()
    {
        return this->my_elements.get_iterator();
        // return this->my_elements.begin();
    }

    /// @brief Get an const_iterator that can be used to access all elements
    ///        owned by this SimulationData object.
    /// @return iterator
    ElementContainer::const_iterator get_const_iterator() const
    {
        return this->my_elements.get_const_iterator();
        // return this->my_elements.cbegin();
    }

    /// @brief Tests whether the given iterator is at the end
    /// @param iter iterator to test
    /// @return true if at end, false otherwise
    bool is_at_end(ElementContainer::iterator iter)
    {
        return this->my_elements.is_at_end(iter);
        // return iter == this->my_elements.end();
    }
    bool is_at_end(ElementContainer::const_iterator citer) const
    {
        return this->my_elements.is_at_end(citer);
        // return citer == this->my_elements.end();
    }

    OpticalPropertySetReference add_optical_property_set(const OpticalPropertySet& opt_set);
    OpticalPropertySetReference find_or_add_optical_property_set(const OpticalPropertySet& opt_set);

    const OpticalPropertySet* get_optical_property_set(const Element& el) const;
    OpticalPropertySet* get_mutable_optical_property_set(const Element& el);

    /// @brief Remove the OpticalPropertySet corresponding to the unique identifier `id`
    /// @param id unique identifier of the OpticalPropertySet to remove
    /// @return true if the element was removed, false otherwise
    auto remove_optical_property_set(ray_source_id id)
    {
        return this->my_optical_property_sets.remove_item(id);
    }

    void clear_optical_property_sets()
    {
        this->my_optical_property_sets.clear();
    }

    /// @brief Get an iterator that can be used to access all 
    ///  optical property sets owned by this SimulationData object.
    /// @return iterator
    OpticalPropertySetContainer::iterator get_optics_iterator()
    {
        return this->my_optical_property_sets.get_iterator();
    }

    /// @brief Tests whether the given iterator is at the end
    /// @param iter iterator to test
    /// @return true if at end, false otherwise
    bool is_optics_at_end(OpticalPropertySetContainer::iterator it)
    {
        return this->my_optical_property_sets.is_at_end(it);
    }

    uint_fast64_t get_number_of_optocal_property_sets()
    {
        return this->my_optical_property_sets.get_number_of_items();
    }

    /// @brief Set the number of rays to trace
    /// @param nrays number of rays to trace
    void set_number_of_rays(uint_fast64_t nrays)
    {
        this->my_parameters.number_of_rays = nrays;
        return;
    }

    /// @brief Get the number of rays to trace
    /// @return number of rays to trace
    uint_fast64_t get_number_of_rays() const
    {
        return this->my_parameters.number_of_rays;
    }

    void set_max_rays_traced(uint_fast64_t nrays)
    {
        this->my_parameters.max_number_of_rays = nrays;
        return;
    }

    uint_fast64_t get_max_number_rays_traced() const
    {
        return this->my_parameters.max_number_of_rays;
    }

    void set_tolerance(double tolerance)
    {
        this->my_parameters.tolerance = tolerance;
        return;
    }
    double get_tolerance() const
    {
        return this->my_parameters.tolerance;
    }

    /// @brief Set the seed used for the random number generation
    /// @param seed seed to set
    void set_seed(uint_fast64_t seed)
    {
        this->my_parameters.seed = seed;
        return;
    }

    /// @brief Get the seed used for random number generation
    /// @return current seed
    uint_fast64_t get_seed() const
    {
        return this->my_parameters.seed;
    }

    void set_latitude(double latitude)
    {
        this->my_parameters.latitude = latitude;
        return;
    }
    double get_latitude() const
    {
        return this->my_parameters.latitude;
    }

    void set_longitude(double longitude)
    {
        this->my_parameters.longitude = longitude;
        return;
    }
    double get_longitude() const
    {
        return this->my_parameters.longitude;
    }

    void set_simulation_date(const Date &d)
    {
        this->my_parameters.sim_dt.my_date = d;
        return;
    }
    const Date &get_simulation_date() const
    {
        return this->my_parameters.sim_dt.my_date;
    }

    void set_simulation_datetime(const DateTime &dt)
    {
        this->my_parameters.sim_dt = dt;
        return;
    }
    const DateTime &get_simulation_datetime() const
    {
        return this->my_parameters.sim_dt;
    }

    void set_simulation_time(const Time &t)
    {
        this->my_parameters.sim_dt.my_time = t;
        return;
    }
    const Time &get_simulation_time() const
    {
        return this->my_parameters.sim_dt.my_time;
    }

    const bool get_include_sun_shape_errors() const
    {
        return this->my_parameters.include_sun_shape_errors;
    }
    void set_include_sun_shape_errors(bool e)
    {
        this->my_parameters.include_sun_shape_errors = e;
    }

    const bool get_include_optical_errors() const
    {
        return this->my_parameters.include_optical_errors;
    }
    void set_include_optical_errors(bool e)
    {
        this->my_parameters.include_optical_errors = e;
    }
    
    const bool get_as_power_tower() const
    {
        return this->my_parameters.as_power_tower;
    }
    void set_as_power_tower(bool pt)
    {
        this->my_parameters.as_power_tower = pt;
    }

    SimulationParameters &get_simulation_parameters()
    {
        return this->my_parameters;
    }
    const SimulationParameters &get_simulation_parameters() const
    {
        return this->my_parameters;
    }

    const std::vector<std::set<uint_fast64_t>>& get_groups() const
    {
        return this->my_groups;
    }

    void set_element_group(uint_fast64_t element_id, int32_t group);

    int update_simulation_positions();
    int update_simulation_positions(const Time &);
    int update_simulation_positions(const Date &);
    int update_simulation_positions(const DateTime &);

    bool import_from_file(const char *file_name);
    bool import_from_file(const std::string file_name);

    /// @brief Import simulation data from a JSON file.
    /// @param file_name Path to the JSON file to import.
    /// @param upgrade_log Optional pointer to a string that will be populated
    ///        with a human-readable description of any schema upgrades applied
    ///        during import. Left untouched if no upgrade was needed. Pass
    ///        nullptr (default) if this information is not needed.
    /// @throws std::runtime_error if the file cannot be read or parsed.
    /// Loads simulation data from the specified JSON file and updates the current simulation state.
    void import_json_file(const std::string file_name, std::string* upgrade_log = nullptr);

    /// @brief Import simulation data from a JSON string.
    /// @param json_str The JSON string to import.
    /// @param upgrade_log Optional pointer to a string that will be populated
    ///        with a human-readable description of any schema upgrades applied
    ///        during import. Left untouched if no upgrade was needed. Pass
    ///        nullptr (default) if this information is not needed.
    /// @throws std::runtime_error if the file cannot be read or parsed.
    /// Loads simulation data from the specified JSON string and updates the current simulation state.
    void import_json_string(const std::string json_str, std::string* upgrade_log = nullptr);
    void import_json_string(const char* json_str, std::string* upgrade_log = nullptr);

    /// @brief Export simulation data to a JSON file.
    /// @param file_name Path to the JSON file to write.
    /// @throws std::runtime_error if the file cannot be written.
    /// Serializes the current simulation data and writes it to the specified JSON file.
    void export_json_file(const std::string file_name);

    /// @brief Clear elements, ray sources, and optical property sets. 
    /// Designed for use by the file stinput and json file readers when there is an error.
    /// @param reset_parameters Optional bool to reset simulation parameters
    void clear(bool reset_parameters = false);

private:
    // mutable element_id next_element_id;

    uint_fast64_t number_of_elements;
    // index is the group number and the value is a set containing the element ids of the group
    std::vector<std::set<uint_fast64_t>> my_groups;

    ElementContainer my_elements;
    RaySourceContainer my_sources;
    SimulationParameters my_parameters;
    OpticalPropertySetContainer my_optical_property_sets;

    // void add_single_element(element_id key, element_ptr el);
    // void add_composite_element(element_id key, element_ptr el);
    // uint_fast64_t remove_single_element(ElementContainer::iterator iter);
    // uint_fast64_t remove_composite_element(element_ptr el);
    uint_fast64_t add_subelements(element_ptr el);
    uint_fast64_t remove_subelements(element_ptr el);

    const OpticalPropertySet* get_optical_property_set(optics_id id) const;
    OpticalPropertySet* get_optical_property_set(optics_id id);

    void insert_element_into_group(element_ptr el, int32_t group);
    void remove_element_from_group(element_ptr el);
};

} // namespace SolTrace::Data

/**
 * @}
 */

#endif

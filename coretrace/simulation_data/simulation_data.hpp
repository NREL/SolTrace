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

#include "container.hpp"
#include "element.hpp"
#include "ray_source.hpp"
#include "simulation_parameters.hpp"

namespace SolTrace::Data {

class SimulationData
{
    friend void load_json_file(SimulationData& sd, std::string filename, std::string* upgrade_log);

public:
    SimulationData();
    virtual ~SimulationData();

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

    SimulationParameters &get_simulation_parameters()
    {
        return this->my_parameters;
    }
    const SimulationParameters &get_simulation_parameters() const
    {
        return this->my_parameters;
    }

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
};

} // namespace SolTrace::Data

/**
 * @}
 */

#endif

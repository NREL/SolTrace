#include "simulation_data.hpp"

#include <cassert>
#include <sstream>

#include "composite_element.hpp"
#include "simdata_io.hpp"
#include "json_schema.hpp"

namespace SolTrace::Data {

SimulationData::SimulationData() : number_of_elements(0),
                                   my_elements(1),
                                   my_sources(0),
                                   my_optical_property_sets(0)
{
    return;
}

SimulationData::~SimulationData()
{
    return;
}

element_id SimulationData::add_element(element_ptr el)
{
    element_id id = ELEMENT_ERROR;

    if (el->get_id() != ELEMENT_ID_UNASSIGNED)
    {
        id = ELEMENT_ALREADY_REGISTERED;
    }
    else
    {
        id = this->my_elements.add_item(el);
        if (Element::is_success(id))
        {
            // Check that all fields required from the user have been specified
            el->enforce_user_fields_set();
            // Make sure coordinate stuff has been computed
            el->compute_coordinate_rotations();
            el->set_id(id);
            if (el->is_composite())
            {
                // The CompositeElement itself does not count toward the
                // number of elements since it does not impact the ray
                // tracing computation. So we do not increment the
                // number of elements here.
                uint_fast64_t n = this->add_subelements(el);
                assert(n == el->get_number_of_elements());
            }
            else
            {
                this->number_of_elements++;
            }
        }
        else
        {
            // TODO: Throw an error here
        }
    }

    return id;
}

// void SimulationData::add_single_element(element_id key,
//                                         element_ptr el)
// {
//     SingleElementMap::value_type to_insert(key, el);
//     auto result = this->my_elements.insert(to_insert);
//     assert(result.second);
//     return;
// }

uint_fast64_t SimulationData::add_subelements(element_ptr el)
{
    // CompositeElementMap::value_type to_insert(key, el);
    // auto result = this->composite_elements.insert(to_insert);
    // assert(result.second);

    composite_element_ptr cptr =
        std::dynamic_pointer_cast<CompositeElement>(el);
    assert(cptr != nullptr);

    uint_fast64_t nadded = cptr->get_number_of_elements();

    auto iter = cptr->get_const_iterator();
    while (!cptr->is_at_end(iter))
    {
        auto id = this->add_element(iter->second);
        assert(Element::is_success(id));
        // // Make sure coordinate stuff has been computed
        // iter->second->compute_coordinate_rotations();
        ++iter;
    }

    return nadded;
}

bool SimulationData::replace_element(element_id id, element_ptr el)
{
    bool success = false;
    assert(el->get_id() == ELEMENT_ID_UNASSIGNED);
    element_ptr old_el = this->get_element(id);
    if (old_el != nullptr)
    {
        success = this->my_elements.replace_item(id, el);
        if (success)
        {
            old_el->set_id(ELEMENT_ID_UNASSIGNED);
            if (old_el->is_composite())
            {
                this->remove_subelements(old_el);
            }
            else
            {
                this->number_of_elements--;
            }
            // Make sure coordinate stuff has been computed
            el->compute_coordinate_rotations();
            el->set_id(id);
            if (el->is_composite())
            {
                this->add_subelements(el);
            }
            else
            {
                this->number_of_elements++;
            }
            // this->number_of_elements -= old_el->get_number_of_elements();
            // this->number_of_elements += el->get_number_of_elements();
        }
    }

    return success;
}

element_ptr SimulationData::get_element(element_id id) const
{
    return this->my_elements.get_item(id);
    // element_ptr retval = nullptr;
    // auto iter1 = this->my_elements.find(id);
    // auto iter2 = this->composite_elements.find(id);
    // if (iter1 != this->my_elements.end())
    // {
    //     retval = iter1->second;
    // }
    // else if (iter2 != this->composite_elements.end())
    // {
    //     retval = iter2->second;
    // }
    // else
    // {
    //     // Intentional no-op
    // }
    // return retval;
}

uint_fast64_t SimulationData::remove_element(element_id id)
{

    uint_fast64_t nremoved = 0;
    element_ptr el = this->my_elements.get_item(id);

    if (el != nullptr)
    {
        this->my_elements.remove_item(id);
        el->set_id(ELEMENT_ID_UNASSIGNED);
        if (el->is_composite())
        {
            // The CompositeElement itself does not count toward the number
            // of elements since it does not impact the ray tracing computation.
            // So we do not decrement the number of elements here.
            nremoved = this->remove_subelements(el);
            assert(nremoved == el->get_number_of_elements());
        }
        else
        {
            this->number_of_elements--;
            nremoved = 1;
        }
    }

    return nremoved;

    // auto iter2 = this->composite_elements.find(id);
    // if (iter1 != this->my_elements.end())
    // {
    //     retval = remove_single_element(iter1);
    // }
    // else if (iter2 != this->composite_elements.end())
    // {
    //     retval = remove_composite_element(iter2);
    // }
    // else
    // {
    //     // Intentional no-op
    // }
    // return retval;
}

// uint_fast64_t SimulationData::remove_single_element(
//     SingleElementMap::iterator iter)
// {
//     iter->second->set_id(ELEMENT_ID_UNASSIGNED);
//     this->my_elements.erase(iter);
//     return 1;
// }

// uint_fast64_t SimulationData::remove_composite_element(element_ptr el)
// {
//     composite_element_ptr cptr = iter->second;
//     iter->second->set_id(ELEMENT_ID_UNASSIGNED);
//     this->composite_elements.erase(iter);
//     return this->remove_subelements(cptr);
// }

uint_fast64_t SimulationData::remove_subelements(element_ptr el)
{
    uint_fast64_t retval = 0;

    composite_element_ptr cptr =
        std::dynamic_pointer_cast<CompositeElement>(el);
    assert(cptr != nullptr);

    element_ptr ptr = nullptr;
    for (auto iter = cptr->get_iterator(); !cptr->is_at_end(iter); ++iter)
    {
        ptr = iter->second;
        retval += this->remove_element(ptr->get_id());
    }

    assert(retval == cptr->get_number_of_elements());

    return retval;
}

OpticalPropertySetReference SimulationData::add_optical_property_set(const OpticalPropertySet& opt_set)
{
    std::shared_ptr<OpticalPropertySet> ptr = std::make_shared<OpticalPropertySet>(opt_set);
    const optics_id id = this->my_optical_property_sets.add_item(ptr);
    return { id, ptr };
}

OpticalPropertySetReference SimulationData::find_or_add_optical_property_set(const OpticalPropertySet& opt_set)
{
    // Check if set already exists
    for (auto it = this->get_optics_iterator(); !this->is_optics_at_end(it); ++it)
    {
        const OpticalPropertySet& existing = *it->second;
        if (existing == opt_set)
        {
            const optics_id id = it->first;
            auto ptr = it->second;
            return { id, ptr };

        }
    }

    // Add if it doesn't exist
    return add_optical_property_set(opt_set);
}

const OpticalPropertySet* SimulationData::get_optical_property_set(const Element& el) const
{
    return get_optical_property_set(el.get_optical_property_set_id());
}

OpticalPropertySet* SimulationData::get_mutable_optical_property_set(const Element& el)
{
    return this->get_optical_property_set(el.get_optical_property_set_id());
}

const OpticalPropertySet* SimulationData::get_optical_property_set(optics_id id) const
{
    auto ptr = this->my_optical_property_sets.get_item(id);
    return ptr == nullptr ? nullptr : ptr.get();
}

OpticalPropertySet* SimulationData::get_optical_property_set(optics_id id)
{
    auto ptr = this->my_optical_property_sets.get_item(id);
    return ptr == nullptr ? nullptr : ptr.get();
}

int SimulationData::update_simulation_positions()
{
    int sts = 0;
    // TODO: Implement this
    return sts;
}

int SimulationData::update_simulation_positions(const Time &t)
{
    this->my_parameters.sim_dt.my_time = t;
    return this->update_simulation_positions();
}

int SimulationData::update_simulation_positions(const Date &d)
{
    this->my_parameters.sim_dt.my_date = d;
    return this->update_simulation_positions();
}

int SimulationData::update_simulation_positions(const DateTime &dt)
{
    this->my_parameters.sim_dt.set(dt);
    return this->update_simulation_positions();
}

bool SimulationData::import_from_file(const char *file_name)
{
    return load_stinput_file(*this, file_name);
}

bool SimulationData::import_from_file(const std::string file_name)
{
    return this->import_from_file(file_name.c_str());
}

void SimulationData::import_json_file(const std::string file_name, std::string* upgrade_log)
{
    load_json_file(*this, file_name, upgrade_log);
}

void SimulationData::export_json_file(const std::string file_name)
{
    write_json_file(*this, file_name);
}

void SimulationData::clear(bool reset_parameters)
{
    this->my_elements.clear();
    this->my_sources.clear();
    this->number_of_elements = 0;

    this->my_optical_property_sets.reset(0);

    if (reset_parameters)
        this->my_parameters = SimulationParameters();   // Reset
}

} // namespace SolTrace::Data

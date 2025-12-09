
#include "composite_element.hpp"

#include <sstream>

#include "single_element.hpp"



namespace SolTrace::Data
{

    CompositeElement::CompositeElement() : ElementBase(),
                                           number_of_elements(0),
                                           my_elements()
    {
        return;
    }

    CompositeElement::CompositeElement(const nlohmann::ordered_json& jnode) : ElementBase(jnode),
        number_of_elements(0),
        my_elements()
    {
        using json = nlohmann::ordered_json;

        // Common parameters are set in ElementBase constructor before here

        json jelements = jnode.at("elements");
        for (auto& [key, jelement] : jelements.items())
        {
            bool is_single = jelement.at("is_single");
            if (is_single)
            {
                element_ptr el = make_element<SingleElement>(jelement);
                this->add_element(el);
            }
            else
            {
                composite_element_ptr comp = make_element<CompositeElement>(jelement);
                this->add_element(comp);
            }
        }
        return;
    }

    CompositeElement::~CompositeElement()
    {
        this->clear();
        return;
    }

    void CompositeElement::disable() const
    {
        this->active = false;
        for (auto iter = this->get_const_iterator();
             !this->is_at_end(iter);
             ++iter)
        {
            iter->second->disable();
        }
        return;
    }

    void CompositeElement::enable() const
    {
        this->active = true;
        for (auto iter = this->get_const_iterator();
             !this->is_at_end(iter);
             ++iter)
        {
            iter->second->enable();
        }
        return;
    }

    void CompositeElement::mark_virtual() const
    {
        this->ElementBase::mark_virtual();
        for (auto iter = this->get_const_iterator();
             !this->is_at_end(iter);
             ++iter)
        {
            iter->second->mark_virtual();
        }
        return;
    }

    void CompositeElement::unmark_virtual() const
    {
        this->ElementBase::unmark_virtual();
        for (auto iter = this->get_const_iterator();
             !this->is_at_end(iter);
             ++iter)
        {
            iter->second->unmark_virtual();
        }
        return;
    }

    void CompositeElement::set_stage(int_fast64_t stage)
    {
        this->stage = stage;
        for (auto iter = this->get_iterator(); !this->is_at_end(iter); ++iter)
        {
            iter->second->set_stage(stage);
        }
        return;
    }

    element_id CompositeElement::add_element(element_ptr el)
    {
        // Cannot nest a StageElement in another CompositeElement
        // assert(!el->is_stage());
        // Disable adding CompositeElements unless `this` is a StageElement
        // assert(this->is_stage() || el->is_single());

        if (el->is_stage() ||
            (el->is_composite() && !this->is_stage()))
        {
            return ELEMENT_INVALID_SETUP;
        }

        el->enforce_user_fields_set();

        element_id id = this->my_elements.add_item(el);
        if (is_success(id))
        {
            this->number_of_elements += el->get_number_of_elements();
            el->set_reference_element(this);
        }
        return id;
    }

    uint_fast64_t CompositeElement::remove_element(element_id id)
    {
        element_ptr el = this->my_elements.get_item(id);
        uint_fast64_t nremoved = this->my_elements.remove_item(id);

        if (nremoved > 0)
        {
            nremoved = el->get_number_of_elements();
            this->number_of_elements -= nremoved;
        }

        return nremoved;
    }

    element_ptr CompositeElement::get_element(element_id id)
    {
        return this->my_elements.get_item(id);
    }

    bool CompositeElement::replace_element(element_id id, element_ptr el)
    {
        element_ptr old_el = this->my_elements.get_item(id);
        bool replaced = this->my_elements.replace_item(id, el);

        if (replaced)
        {
            this->number_of_elements -= old_el->get_number_of_elements();
            this->number_of_elements += el->get_number_of_elements();
        }

        return replaced;
    }

    void CompositeElement::clear()
    {
        this->number_of_elements = 0;
        this->my_elements.clear();
        return;
    }

    void CompositeElement::enforce_user_fields_set() const
    {
        ElementBase::enforce_user_fields_set();

        if (this->number_of_elements == 0)
        {
            std::stringstream ss;
            ss << "CompositeElement (Name: " << this->get_name()
               << ", UUID: " << this->get_id()
               << ") has no subelements.";
            throw std::invalid_argument(ss.str());
        }
        return;
    }

    // Stage and Composite should have the same function
    void CompositeElement::write_json(nlohmann::ordered_json& jnode) const
    {
        using json = nlohmann::ordered_json;

        // Composite/stage specific info
        jnode["is_composite"] = true;
        jnode["is_single"] = false;
        jnode["number_of_elements"] = this->get_number_of_elements();

        // Common properties of all elements
        this->write_common_json(jnode);

        // Loop through each element, writing json node
        json jelements;
        for (auto it = this->get_const_iterator(); !this->is_at_end(it); ++it)
        {
            json jelement;

            SolTrace::Data::element_id i = it->first;
            auto element = it->second;

            element->write_json(jelement);

            jelements[std::to_string(i)] = jelement;
        }
        jnode["elements"] = jelements;
    }

} // namespace SolTrace::Data

/**
 * @file composite_element.hpp
 * @brief Composite element class for complex optical systems
 *
 * Defines the CompositeElement class which can contain multiple
 * sub-elements, allowing for hierarchical optical system definitions.
 * Enables grouping of related optical elements for easier management
 * and coordinate transformations.
 *
 * @defgroup elements Optical Elements
 * @{
 */

#ifndef SOLTRACE_COMPOSITE_ELEMENT_H
#define SOLTRACE_COMPOSITE_ELEMENT_H

#include <memory>

#include "container.hpp"
#include "element.hpp"

namespace SolTrace::Data
{

    class CompositeElement : public ElementBase
    {
    public:
        /**
         * @brief Default constructor for composite element
         */
        CompositeElement();
        CompositeElement(const nlohmann::ordered_json& jnode);
        virtual ~CompositeElement();

        /**
         * @brief Disable this composite element and all sub-elements
         */
        virtual void disable() const override;

        /**
         * @brief Enable this composite element and all sub-elements
         */
        virtual void enable() const override;

        /**
         * @brief Check if this element is composite
         * @return Always returns true for composite elements
         */
        virtual bool is_composite() const override
        {
            return true;
        }

        virtual void mark_virtual() const override;
        virtual void unmark_virtual() const override;

        /**
         * @brief Set the stage number for this composite element
         * @param stage Stage number to assign
         */
        virtual void set_stage(int_fast64_t stage) override;

        /**
         * @brief Get the number of sub-elements in this composite
         * @return Number of elements contained in this composite
         */
        virtual uint_fast64_t get_number_of_elements() const override
        {
            // return this->my_elements.get_number_of_items();
            return this->number_of_elements;
        }

        // Element interface functions
        /**
         * @brief Get aperture pointer (always null for composite elements)
         * @return nullptr (composite elements don't have apertures)
         */
        virtual const aperture_ptr get_aperture() const override { return nullptr; }

        /**
         * @brief Get aperture pointer (always null for composite elements)
         * @return nullptr (composite elements don't have apertures)
         */
        virtual aperture_ptr get_aperture() override { return nullptr; }

        /**
         * @brief Set aperture (no-op for composite elements)
         * @param aperture_ptr Ignored for composite elements
         */
        virtual void set_aperture(aperture_ptr) override {}

        /**
         * @brief Get surface pointer (always null for composite elements)
         * @return nullptr (composite elements don't have surfaces)
         */
        virtual const surface_ptr get_surface() const override { return nullptr; }

        /**
         * @brief Get surface pointer (always null for composite elements)
         * @return nullptr (composite elements don't have surfaces)
         */
        virtual surface_ptr get_surface() override { return nullptr; }

        /**
         * @brief Set surface (no-op for composite elements)
         * @param surface_ptr Ignored for composite elements
         */
        virtual void set_surface(surface_ptr) override {}

        /**
         * @brief Get front optical properties (always null for composite elements)
         * @return nullptr (composite elements don't have optical properties)
         */
        virtual const OpticalProperties *get_front_optical_properties() const override
        {
            return nullptr;
        }

        /**
         * @brief Get front optical properties (always null for composite elements)
         * @return nullptr (composite elements don't have optical properties)
         */
        virtual OpticalProperties *get_front_optical_properties() override
        {
            return nullptr;
        }

        /**
         * @brief Set front optical properties (no-op for composite elements)
         * @param op Ignored for composite elements
         */
        virtual void set_front_optical_properties(const OpticalProperties &) override {}

        /**
         * @brief Get back optical properties (always null for composite elements)
         * @return nullptr (composite elements don't have optical properties)
         */
        virtual const OpticalProperties *get_back_optical_properties() const override
        {
            return nullptr;
        }

        /**
         * @brief Get back optical properties (always null for composite elements)
         * @return nullptr (composite elements don't have optical properties)
         */
        virtual OpticalProperties *get_back_optical_properties() override
        {
            return nullptr;
        }

        /**
         * @brief Set back optical properties (no-op for composite elements)
         * @param op Ignored for composite elements
         */
        virtual void set_back_optical_properties(const OpticalProperties &) override {};

        // CompositeElement accessors
        /**
         * @brief Add an element to this composite
         * @param el Shared pointer to element to add
         * @return Element ID of the added element
         */
        element_id add_element(element_ptr el);

        /**
         * @brief Remove an element from this composite
         * @param id Element ID to remove
         * @return Number of elements removed (0 or 1)
         */
        uint_fast64_t remove_element(element_id id);

        /**
         * @brief Get an element by ID
         * @param id Element ID to retrieve
         * @return Shared pointer to element, or null if not found
         */
        element_ptr get_element(element_id id);

        /**
         * @brief Replace an element with a new one
         * @param id Element ID to replace
         * @param el New element to insert
         * @return True if replacement was successful
         */
        bool replace_element(element_id id, element_ptr el);

        /**
         * @brief Remove all elements from this composite
         */
        void clear();

        // uint64_t get_total_number_of_elements() const
        // {
        //     return this->my_elements.get_total_number_of_items();
        // }

        /**
         * @brief Get iterator to beginning of element container
         * @return Iterator to first element
         */
        virtual ElementContainer::iterator get_iterator()
        {
            return this->my_elements.get_iterator();
        }
        virtual ElementContainer::const_iterator get_const_iterator() const
        {
            return this->my_elements.get_const_iterator();
        }
        virtual bool is_at_end(ElementContainer::iterator iter)
        {
            return this->my_elements.is_at_end(iter);
        }
        virtual bool is_at_end(ElementContainer::const_iterator citer) const
        {
            return this->my_elements.is_at_end(citer);
        }

        virtual void enforce_user_fields_set() const override;

        virtual void write_json(nlohmann::ordered_json& jnode) const override;

    private:
        uint_fast64_t number_of_elements;
        ElementContainer my_elements;
    };

    using composite_element_ptr = std::shared_ptr<CompositeElement>;

} // namespace SolTrace::Data

#endif

/**
 * @}
 */

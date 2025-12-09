/**
 * @file single_element.hpp
 * @brief Single optical element implementation
 *
 * Defines the SingleElement class for individual optical surfaces
 * with associated apertures, surfaces, and optical properties.
 * Represents the fundamental building block of optical systems
 * in SolTrace ray tracing simulations.
 *
 * @defgroup elements Optical Elements
 * @{
 */

#ifndef SOLTRACE_PLANE_H
#define SOLTRACE_PLANE_H

#include "element.hpp"

#include <memory>

#include "aperture.hpp"
#include "surface.hpp"

namespace SolTrace::Data {

class SingleElement : public ElementBase
{
public:
    SingleElement();
    SingleElement(const nlohmann::ordered_json& jnode);
    virtual ~SingleElement();

    virtual bool is_single() const override { return true; }

    virtual const aperture_ptr get_aperture() const override
    {
        return this->aperture;
    }
    virtual aperture_ptr get_aperture() override
    {
        return this->aperture;
    }
    virtual void set_aperture(aperture_ptr ap) override
    {
        this->aperture = ap;
        return;
    }

    virtual const surface_ptr get_surface() const override
    {
        return this->surface;
    }
    virtual surface_ptr get_surface() override
    {
        return this->surface;
    }
    virtual void set_surface(surface_ptr sp) override
    {
        this->surface = sp;
        return;
    }

    const OpticalProperties *get_front_optical_properties() const override
    {
        return &(this->optics_front);
    }
    OpticalProperties *get_front_optical_properties() override
    {
        return &(this->optics_front);
    }
    void set_front_optical_properties(const OpticalProperties &op) override
    {
        this->optics_front = op;
        return;
    }

    const OpticalProperties *get_back_optical_properties() const override
    {
        return &(this->optics_back);
    }
    OpticalProperties *get_back_optical_properties() override
    {
        return &(this->optics_back);
    }
    void set_back_optical_properties(const OpticalProperties &op) override
    {
        this->optics_back = op;
    }

    virtual void enforce_user_fields_set() const override;

    virtual void write_json(nlohmann::ordered_json& jnode) const override;

protected:
    aperture_ptr aperture;
    surface_ptr surface;

    OpticalProperties optics_front;
    OpticalProperties optics_back;
};

using single_element_ptr = typename std::shared_ptr<SingleElement>;

} // namespace SolTrace::Data

/**
 * @}
 */

#endif

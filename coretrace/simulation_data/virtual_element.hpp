/**
 * @file virtual_element.hpp
 * @brief Virtual element for non-interacting ray tracking
 *
 * Defines virtual elements that can track rays without optical
 * interaction, useful for coordinate system definitions and ray tracking.
 * Virtual elements allow rays to pass through without reflection,
 * refraction, or absorption, while still recording ray intersections.
 */

#ifndef SOLTRACE_VIRTUAL_ELEMENT_H
#define SOLTRACE_VIRTUAL_ELEMENT_H

#include "single_element.hpp"

namespace SolTrace::Data {

class VirtualElement : public SingleElement
{
public:
    VirtualElement();
    VirtualElement(const nlohmann::ordered_json& jnode,
        const OpticalPropertySetResolver& resolve_optics);
    virtual ~VirtualElement();

    virtual bool is_virtual() const override { return true; }

    //virtual void set_optical_property_set_id(optics_id) override {}
    //virtual optics_id get_optical_property_set_id() const override 
    //{
    //    return opt_id;  // TODO: return identifier saying it's permanent/virtual
    //}
private:
    optical_set_ptr owned_optical_property_set;
};

class VirtualPlane : public VirtualElement
{
public:
    VirtualPlane(double x_len, double y_len);
    virtual ~VirtualPlane();

    void set_aperture(aperture_ptr ap) override {}
    void set_surface(surface_ptr sp) override {}

private:
};

} // namespace SolTrace::Data

#endif

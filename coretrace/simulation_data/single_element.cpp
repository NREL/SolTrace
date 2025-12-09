
#include "single_element.hpp"

#include <exception>
#include <memory>
#include <sstream>

#include "aperture.hpp"
#include "element.hpp"
#include "vector3d.hpp"

namespace SolTrace::Data {

SingleElement::SingleElement() : ElementBase(),
                                 aperture(nullptr),
                                 surface(nullptr),
                                 optics_front(),
                                 optics_back()
{
    this->optics_front.set_ideal_absorption();
    this->optics_back.set_ideal_absorption();
    return;
}

SingleElement::SingleElement(const nlohmann::ordered_json& jnode) : ElementBase(jnode),
                                                                    aperture(nullptr),
                                                                    surface(nullptr),
                                                                    optics_front(),
                                                                    optics_back()
{
    this->set_aperture(Aperture::make_aperture_from_json(jnode.at("aperture")));
    this->set_surface(make_surface_from_json(jnode.at("surface")));
    this->set_front_optical_properties(OpticalProperties(jnode.at("optics_front")));
    this->set_back_optical_properties(OpticalProperties(jnode.at("optics_back")));
}

SingleElement::~SingleElement()
{
    this->aperture = nullptr;
    this->surface = nullptr;
    return;
}

void SingleElement::enforce_user_fields_set() const
{
    ElementBase::enforce_user_fields_set();

    if (this->aperture == nullptr)
    {
        std::stringstream ss;
        ss << "Element (Name: " << this->get_name()
           << ", UUID: " << this->get_id()
           << ") has no aperture.";
        throw std::invalid_argument(ss.str());
    }

    if (this->surface == nullptr)
    {
        std::stringstream ss;
        ss << "Element (Name: " << this->get_name()
           << ", UUID: " << this->get_id()
           << ") has no surface.";
        throw std::invalid_argument(ss.str());
    }

    return;
}

void SingleElement::write_json(nlohmann::ordered_json& jnode) const
{
    using json = nlohmann::ordered_json;
    
    // Write shared properties
    this->write_common_json(jnode);
    
    // Optical Properties
    json joptics_front, joptics_back;
    this->optics_front.write_json(joptics_front);
    this->optics_back.write_json(joptics_back);
    jnode["optics_front"] = joptics_front;
    jnode["optics_back"] = joptics_back;

    // Aperture
    json japerture;
    this->aperture->write_json(japerture);
    jnode["aperture"] = japerture;

    // Surface
    json jsurface;
    this->surface->write_json(jsurface);
    jnode["surface"] = jsurface;

    // Write element type
    jnode["is_single"] = true;
}

} // namespace SolTrace::Data

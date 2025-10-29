
#include "virtual_element.hpp"

namespace SolTrace::Data {

void VirtualElement::set_virtual_optics(OpticalProperties &op)
{
    op.my_type = InteractionType::REFRACTION;
    op.reflectivity = 0.0;
    op.slope_error = 0.0;
    op.specularity_error = 0.0;
    op.transmitivity = 1.0;
    op.refraction_index_front = 1.0;
    op.refraction_index_back = 1.0;
    return;
}

VirtualElement::VirtualElement() : SingleElement()
{
    set_virtual_optics(this->optics_front);
    set_virtual_optics(this->optics_back);
}

VirtualElement::~VirtualElement()
{
}

VirtualPlane::VirtualPlane(double x_len, double y_len) : VirtualElement()
{
    this->aperture = make_aperture<Rectangle>(x_len, y_len);
    this->surface = make_surface<Flat>();
}

VirtualPlane::~VirtualPlane()
{
}

} // namespace SolTrace::Data

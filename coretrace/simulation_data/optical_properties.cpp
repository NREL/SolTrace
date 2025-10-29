
// #include <exception>
// #include <sstream>

#include "optical_properties.hpp"

namespace SolTrace::Data
{

const std::string interaction_string(InteractionType it)
{
    if (it == InteractionType::REFLECTION)
    {
        return "Reflection";
    }
    else if(it == InteractionType::REFRACTION)
    {
        return "Refraction";
    }
    else
    {
        return "Unknown";
    }
}

std::ostream &operator<<(std::ostream &os,
    const OpticalProperties &op)
{
    os << "Type: " << interaction_string(op.my_type)
    << "\nTransmitivity: " << op.transmitivity
    << "\nReflectivity: " << op.reflectivity
    << "\nSlope Error: " << op.slope_error
    << "\nSpecularity Error: " << op.specularity_error
    << "\nRefraction Index Front: " << op.refraction_index_front
    << "\nRefraction Index Back: " << op.refraction_index_back;
    return os;
}

} // namespace SolTrace::Data

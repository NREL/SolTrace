
// #include <exception>
// #include <sstream>

#include "optical_properties.hpp"
#include "simdata_io.hpp"

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

OpticalProperties::OpticalProperties(const nlohmann::ordered_json& jnode)
{
    std::string interaction_string = jnode.at("my_type");
    this->my_type = get_enum_from_string(interaction_string, InteractionTypeMap, InteractionType::UNKNOWN);

    std::string distribution_string = jnode.at("error_distribution_type");
    this->error_distribution_type = get_enum_from_string(distribution_string, DistributionTypeMap, DistributionType::UNKNOWN);

    this->transmitivity = jnode.at("transmissivity");
    this->reflectivity = jnode.at("reflectivity");
    this->slope_error = jnode.at("slope_error");
    this->specularity_error = jnode.at("specularity_error");
    this->refraction_index_front = jnode.at("refraction_index_front");
    this->refraction_index_back = jnode.at("refraction_index_back");
}

void OpticalProperties::write_json(nlohmann::ordered_json& jnode) const
{
    jnode["my_type"] = InteractionTypeMap.at(this->my_type);
    jnode["error_distribution_type"] = DistributionTypeMap.at(this->error_distribution_type);
    jnode["transmissivity"] = this->transmitivity;
    jnode["reflectivity"] = this->reflectivity;
    jnode["slope_error"] = this->slope_error;
    jnode["specularity_error"] = this->specularity_error;
    jnode["refraction_index_front"] = this->refraction_index_front;
    jnode["refraction_index_back"] = this->refraction_index_back;
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

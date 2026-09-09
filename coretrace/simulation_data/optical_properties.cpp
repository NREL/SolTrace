
#include <cmath>

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

OpticalPropertySet::OpticalPropertiesFace::OpticalPropertiesFace(const nlohmann::ordered_json& jnode)
{
    std::string distribution_string = jnode.at("error_distribution_type");
    this->error_distribution_type = get_enum_from_string(distribution_string, DistributionTypeMap, DistributionType::UNKNOWN);

    this->transmissivity = jnode.at("transmissivity");
    this->reflectivity = jnode.at("reflectivity");
    this->slope_error = jnode.at("slope_error");
    this->specularity_error = jnode.at("specularity_error");

    validate();
}

void OpticalPropertySet::OpticalPropertiesFace::validate() const
{
    if (this->error_distribution_type == DistributionType::UNKNOWN)
        throw std::invalid_argument("Optical properties: unknown error distribution type");

    if (!std::isfinite(this->slope_error) || this->slope_error < 0.0)
        throw std::invalid_argument("Optical properties: slope error must be finite and non-negative");

    if (!std::isfinite(this->specularity_error) || this->specularity_error < 0.0)
        throw std::invalid_argument("Optical properties: specularity error must be finite and non-negative");
}

void OpticalPropertySet::OpticalPropertiesFace::write_json(nlohmann::ordered_json& jnode) const
{
    jnode["error_distribution_type"] = DistributionTypeMap.at(this->error_distribution_type);
    jnode["transmissivity"] = this->transmissivity;
    jnode["reflectivity"] = this->reflectivity;
    jnode["slope_error"] = this->slope_error;
    jnode["specularity_error"] = this->specularity_error;
}

bool OpticalPropertySet::OpticalPropertiesFace::operator==(const OpticalPropertiesFace& other) const
{
    return this->error_distribution_type == other.error_distribution_type &&
        this->transmissivity == other.transmissivity &&
        this->reflectivity == other.reflectivity &&
        this->slope_error == other.slope_error &&
        this->specularity_error == other.specularity_error;
}

bool OpticalPropertySet::OpticalPropertiesFace::operator!=(const OpticalPropertiesFace& other) const
{
    return !(*this == other);
}

std::ostream &operator<<(std::ostream &os,
    const OpticalPropertySet::OpticalPropertiesFace& op)
{
    os << "Transmissivity: " << op.transmissivity
       << "\nReflectivity: " << op.reflectivity
       << "\nSlope Error: " << op.slope_error
       << "\nSpecularity Error: " << op.specularity_error;
    return os;
}

OpticalPropertySet::OpticalPropertySet(const nlohmann::ordered_json& jnode)
    : front(jnode.at("front")), back(jnode.at("back"))
{
    std::string interaction_type_string = jnode.at("my_type");
    this->my_type = get_enum_from_string(interaction_type_string, InteractionTypeMap, InteractionType::UNKNOWN);

    this->refraction_index_front = jnode.at("refraction_index_front");
    this->refraction_index_back = jnode.at("refraction_index_back");

    this->my_name = jnode.at("my_name");
}

void OpticalPropertySet::write_json(nlohmann::ordered_json& jnode) const
{
    jnode["my_type"] = InteractionTypeMap.at(this->my_type);
    jnode["refraction_index_front"] = this->refraction_index_front;
    jnode["refraction_index_back"] = this->refraction_index_back;
    jnode["my_name"] = this->my_name;

    front.write_json(jnode["front"]);
    back.write_json(jnode["back"]);
}

bool OpticalPropertySet::operator==(const OpticalPropertySet& other) const
{
    return this->front == other.front &&
        this->back == other.back &&
        this->my_type == other.my_type &&
        this->refraction_index_front == other.refraction_index_front &&
        this->refraction_index_back == other.refraction_index_back &&
        this->my_name == other.my_name;
}

bool OpticalPropertySet::operator!=(const OpticalPropertySet& other) const
{
    return !(*this == other);
}

std::ostream& operator<<(std::ostream& os,
    const OpticalPropertySet& op)
{
    os << "Type: " << interaction_string(op.my_type)
       << "\nRefraction Index Front: " << op.refraction_index_front
       << "\nRefraction Index Back: " << op.refraction_index_back
       << "\nFront: " << op.front
       << "\nBack: " << op.back;
    return os;
}

} // namespace SolTrace::Data

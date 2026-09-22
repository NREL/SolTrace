
#include <cmath>

#include "optical_properties.hpp"
#include "simdata_io.hpp"

namespace SolTrace::Data
{

namespace
{
constexpr double kMradToRad = 1.0e-3;
}

const std::string interaction_string(InteractionType it)
{
    if (it == InteractionType::REFLECTION) { return "Reflection"; }
    else if (it == InteractionType::REFRACTION) { return "Refraction"; }
    else
    {
        return "Unknown";
    }
}

namespace
{
std::vector<AngularTablePoint>
read_angular_table(const nlohmann::ordered_json& jnode, const char* key)
{
    std::vector<AngularTablePoint> table;

    if (!jnode.contains(key)) return table;

    for (const auto& point_node : jnode.at(key))
    {
        table.push_back(AngularTablePoint { point_node.at("angle"),
                                            point_node.at("value") });
    }

    return table;
}

void write_angular_table(nlohmann::ordered_json&               jnode,
                         const char*                           key,
                         const std::vector<AngularTablePoint>& table)
{
    nlohmann::ordered_json table_node = nlohmann::ordered_json::array();

    for (const auto& point : table)
    {
        nlohmann::ordered_json point_node;
        point_node["angle"] = point.angle;
        point_node["value"] = point.value;
        table_node.push_back(point_node);
    }

    jnode[key] = table_node;
}
} // namespace

OpticalPropertySet::OpticalPropertiesFace::OpticalPropertiesFace(
    const nlohmann::ordered_json& jnode)
{
    std::string distribution_string = jnode.at("error_distribution_type");
    this->error_distribution_type   = get_enum_from_string(
        distribution_string, DistributionTypeMap, DistributionType::UNKNOWN);

    this->transmissivity    = jnode.at("transmissivity");
    this->reflectivity      = jnode.at("reflectivity");
    this->slope_error       = jnode.at("slope_error");
    this->specularity_error = jnode.at("specularity_error");

    this->use_reflectivity_table = jnode.value("use_reflectivity_table", false);
    this->use_transmissivity_table =
        jnode.value("use_transmissivity_table", false);
    this->reflectivity_table = read_angular_table(jnode, "reflectivity_table");
    this->transmissivity_table =
        read_angular_table(jnode, "transmissivity_table");

    validate();
}

namespace
{
void validate_angular_table(const std::vector<AngularTablePoint>& table,
                            bool                                  use_table,
                            const char*                           name)
{
    if (!use_table) return;

    if (table.empty())
        throw std::invalid_argument(std::string("Optical properties: ") + name +
                                    " table is enabled but empty");

    for (std::size_t i = 0; i < table.size(); ++i)
    {
        if (!std::isfinite(table[i].angle) || table[i].angle < 0.0)
            throw std::invalid_argument(
                std::string("Optical properties: ") + name +
                " table angle must be finite and non-negative");

        if (!std::isfinite(table[i].value))
            throw std::invalid_argument(std::string("Optical properties: ") +
                                        name + " table value must be finite");

        if (table[i].value < 0.0 || table[i].value > 1.0)
            throw std::invalid_argument(std::string("Optical properties: ") +
                                        name + " table value must be in [0,1]");

        if (i > 0 && table[i].angle <= table[i - 1].angle)
            throw std::invalid_argument(
                std::string("Optical properties: ") + name +
                " table angles must be strictly ascending");
    }
}
} // namespace

void OpticalPropertySet::OpticalPropertiesFace::validate() const
{
    if (this->error_distribution_type == DistributionType::UNKNOWN)
        throw std::invalid_argument(
            "Optical properties: unknown error distribution type");

    if (!std::isfinite(this->slope_error) || this->slope_error < 0.0)
        throw std::invalid_argument(
            "Optical properties: slope error must be finite and non-negative");

    if (!std::isfinite(this->specularity_error) ||
        this->specularity_error < 0.0)
        throw std::invalid_argument("Optical properties: specularity error "
                                    "must be finite and non-negative");

    validate_angular_table(
        this->reflectivity_table, this->use_reflectivity_table, "reflectivity");
    validate_angular_table(this->transmissivity_table,
                           this->use_transmissivity_table,
                           "transmissivity");
}

double OpticalPropertySet::OpticalPropertiesFace::get_reflectivity(
    double incident_angle_mrad) const
{
    return this->use_reflectivity_table
               ? interpolate_table(this->reflectivity_table,
                                   incident_angle_mrad)
               : this->reflectivity;
}

double OpticalPropertySet::OpticalPropertiesFace::get_transmissivity(
    double incident_angle_mrad) const
{
    return this->use_transmissivity_table
               ? interpolate_table(this->transmissivity_table,
                                   incident_angle_mrad)
               : this->transmissivity;
}

double OpticalPropertySet::OpticalPropertiesFace::interpolate_table(
    const std::vector<AngularTablePoint>& table,
    double                                incident_angle_mrad)
{
    // The table is sorted ascending by angle, i.e. descending by cos(angle),
    // so interpolation is performed linearly in cos(angle) rather than angle.
    if (table.size() == 1) return table.front().value;

    const double c              = std::cos(incident_angle_mrad * kMradToRad);
    const double c_at_min_angle = std::cos(table.front().angle * kMradToRad);
    const double c_at_max_angle = std::cos(table.back().angle * kMradToRad);

    if (c >= c_at_min_angle) return table.front().value;

    if (c <= c_at_max_angle) return table.back().value;

    std::size_t k = 1;
    while (std::cos(table[k].angle * kMradToRad) > c)
        ++k;

    const double c0 = std::cos(table[k - 1].angle * kMradToRad);
    const double c1 = std::cos(table[k].angle * kMradToRad);
    const double v0 = table[k - 1].value;
    const double v1 = table[k].value;

    return v0 + (c - c0) / (c1 - c0) * (v1 - v0);
}

void OpticalPropertySet::OpticalPropertiesFace::write_json(
    nlohmann::ordered_json& jnode) const
{
    jnode["error_distribution_type"] =
        DistributionTypeMap.at(this->error_distribution_type);
    jnode["transmissivity"]    = this->transmissivity;
    jnode["reflectivity"]      = this->reflectivity;
    jnode["slope_error"]       = this->slope_error;
    jnode["specularity_error"] = this->specularity_error;

    jnode["use_reflectivity_table"]   = this->use_reflectivity_table;
    jnode["use_transmissivity_table"] = this->use_transmissivity_table;
    write_angular_table(jnode, "reflectivity_table", this->reflectivity_table);
    write_angular_table(
        jnode, "transmissivity_table", this->transmissivity_table);
}

bool OpticalPropertySet::OpticalPropertiesFace::operator==(
    const OpticalPropertiesFace& other) const
{
    return this->error_distribution_type == other.error_distribution_type &&
           this->transmissivity == other.transmissivity &&
           this->reflectivity == other.reflectivity &&
           this->slope_error == other.slope_error &&
           this->specularity_error == other.specularity_error &&
           this->use_reflectivity_table == other.use_reflectivity_table &&
           this->use_transmissivity_table == other.use_transmissivity_table &&
           this->reflectivity_table == other.reflectivity_table &&
           this->transmissivity_table == other.transmissivity_table;
}

bool OpticalPropertySet::OpticalPropertiesFace::operator!=(
    const OpticalPropertiesFace& other) const
{ return !(*this == other); }

std::ostream& operator<<(std::ostream&                                    os,
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
    this->my_type                       = get_enum_from_string(
        interaction_type_string, InteractionTypeMap, InteractionType::UNKNOWN);

    this->refraction_index_front = jnode.at("refraction_index_front");
    this->refraction_index_back  = jnode.at("refraction_index_back");

    this->my_name = jnode.at("my_name");
}

void OpticalPropertySet::write_json(nlohmann::ordered_json& jnode) const
{
    jnode["my_type"]                = InteractionTypeMap.at(this->my_type);
    jnode["refraction_index_front"] = this->refraction_index_front;
    jnode["refraction_index_back"]  = this->refraction_index_back;
    jnode["my_name"]                = this->my_name;

    front.write_json(jnode["front"]);
    back.write_json(jnode["back"]);
}

bool OpticalPropertySet::operator==(const OpticalPropertySet& other) const
{
    return this->front == other.front && this->back == other.back &&
           this->my_type == other.my_type &&
           this->refraction_index_front == other.refraction_index_front &&
           this->refraction_index_back == other.refraction_index_back &&
           this->my_name == other.my_name;
}

bool OpticalPropertySet::operator!=(const OpticalPropertySet& other) const
{ return !(*this == other); }

std::ostream& operator<<(std::ostream& os, const OpticalPropertySet& op)
{
    os << "Type: " << interaction_string(op.my_type)
       << "\nRefraction Index Front: " << op.refraction_index_front
       << "\nRefraction Index Back: " << op.refraction_index_back
       << "\nFront: " << op.front << "\nBack: " << op.back;
    return os;
}

} // namespace SolTrace::Data

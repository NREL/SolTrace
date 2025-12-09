
#include "surface.hpp"
#include "simdata_io.hpp"

#include <vector>

namespace SolTrace::Data {

surface_ptr make_surface_from_type(SurfaceType type, const std::vector<double> &args)
{
    surface_ptr retval = nullptr;
    unsigned nargs = args.size();

    switch (type)
    {
    case SurfaceType::CONE:
        retval = nargs < 1 ? nullptr : make_surface<Cone>(args[0]);
        break;
    case SurfaceType::CYLINDER:
        retval = nargs < 1 ? nullptr : make_surface<Cylinder>(1.0 / args[0]);
        break;
    case SurfaceType::FLAT:
        retval = make_surface<Flat>();
        break;
    case SurfaceType::PARABOLA:
    {
        if (nargs < 2)
            retval = nullptr;
        else
        {
            double cx = args[0];
            double cy = args[1];
            double fx = 1.0 / (2.0 * cx);
            double fy = 1.0 / (2.0 * cy);
            retval = make_surface<Parabola>(fx, fy);
        }
        break;
    } 
    case SurfaceType::SPHERE:
        retval = nargs < 1 ? nullptr : make_surface<Sphere>(args[0]);
        break;
    case SurfaceType::HYPER:
    case SurfaceType::GENERAL_SPENCER_MURTY:
    case SurfaceType::TORUS:
    default:
        retval = nullptr; // Not implemented yet
        break;
    }

    return retval;
}

surface_ptr make_surface_from_json(const nlohmann::ordered_json& jnode)
{
    if (!jnode.contains("surface_type"))
        throw std::invalid_argument("Missing surface_type");
    std::string type_str = jnode.at("surface_type");
    SurfaceType surface_type = get_enum_from_string(type_str, SurfaceTypeMap, SurfaceType::SURFACE_UNKNOWN);
    if (surface_type == SurfaceType::SURFACE_UNKNOWN)
        throw std::invalid_argument("Unknown surface");
    switch (surface_type)
    {
        case SurfaceType::CONE:                 return make_surface<Cone>(jnode);
        case SurfaceType::CYLINDER:             return make_surface<Cylinder>(jnode);
        case SurfaceType::FLAT:                 return make_surface<Flat>(jnode);
        case SurfaceType::PARABOLA:             return make_surface<Parabola>(jnode);
        case SurfaceType::SPHERE:               return make_surface<Sphere>(jnode);
        default:
            throw std::invalid_argument("Unsupported surface_type: " + type_str);
    }
}

Cone::Cone(const nlohmann::ordered_json& jnode)
    : Surface(SurfaceType::CONE)
{
    this->half_angle = jnode.at("half_angle");
}

void Cone::write_json(nlohmann::ordered_json& jnode) const
{
    SurfaceType type = SurfaceType::CONE;
    jnode["surface_type"] = SurfaceTypeMap.at(type);
    jnode["half_angle"] = this->half_angle;
}

Cylinder::Cylinder(const nlohmann::ordered_json& jnode)
    : Surface(SurfaceType::CYLINDER)
{
    this->radius = jnode.at("radius");
}

void Cylinder::write_json(nlohmann::ordered_json& jnode) const
{
    SurfaceType type = SurfaceType::CYLINDER;
    jnode["surface_type"] = SurfaceTypeMap.at(type);
    jnode["radius"] = this->radius;
}

void Flat::write_json(nlohmann::ordered_json& jnode) const
{
    SurfaceType type = SurfaceType::FLAT;
    jnode["surface_type"] = SurfaceTypeMap.at(type);
}

Parabola::Parabola(const nlohmann::ordered_json& jnode)
    : Surface(SurfaceType::PARABOLA)
{
    this->focal_length_x = jnode.at("focal_length_x");
    this->focal_length_y = jnode.at("focal_length_y");
}

void Parabola::write_json(nlohmann::ordered_json& jnode) const
{
    SurfaceType type = SurfaceType::PARABOLA;
    jnode["surface_type"] = SurfaceTypeMap.at(type);
    jnode["focal_length_x"] = this->focal_length_x;
    jnode["focal_length_y"] = this->focal_length_y;
}

Sphere::Sphere(const nlohmann::ordered_json& jnode)
    : Surface(SurfaceType::SPHERE)
{
    this->vertex_curv = jnode.at("vertex_curv");
}

void Sphere::write_json(nlohmann::ordered_json& jnode) const
{
    SurfaceType type = SurfaceType::SPHERE;
    jnode["surface_type"] = SurfaceTypeMap.at(type);
    jnode["vertex_curv"] = this->vertex_curv;
}

} // namespace SolTrace::Data


#include "surface.hpp"
#include "simdata_io.hpp"

#include <cmath>
#include <stdexcept>
#include <vector>

#include <utilities.hpp>

namespace SolTrace::Data
{

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

    surface_ptr make_surface_from_json(const nlohmann::ordered_json &jnode)
    {
        if (!jnode.contains("surface_type"))
            throw std::invalid_argument("Missing surface_type");
        std::string type_str = jnode.at("surface_type");
        SurfaceType surface_type = get_enum_from_string(type_str, SurfaceTypeMap, SurfaceType::SURFACE_UNKNOWN);
        if (surface_type == SurfaceType::SURFACE_UNKNOWN)
            throw std::invalid_argument("Unknown surface");
        switch (surface_type)
        {
        case SurfaceType::CONE:
            return make_surface<Cone>(jnode);
        case SurfaceType::CYLINDER:
            return make_surface<Cylinder>(jnode);
        case SurfaceType::FLAT:
            return make_surface<Flat>(jnode);
        case SurfaceType::PARABOLA:
            return make_surface<Parabola>(jnode);
        case SurfaceType::SPHERE:
            return make_surface<Sphere>(jnode);
        default:
            throw std::invalid_argument("Unsupported surface_type: " + type_str);
        }
    }

    Cone::Cone(const nlohmann::ordered_json &jnode)
        : Surface(SurfaceType::CONE)
    {
        this->half_angle = jnode.at("half_angle");
        validate();
    }

    void Cone::validate() const
    {
        if (!std::isfinite(half_angle) || half_angle <= 0.0 || half_angle >= 90.0)
            throw std::invalid_argument("Cone: half_angle must be in (0, 90) degrees");
    }

    void Cone::bounding_box(const double x_minmax[2],
                            const double y_minmax[2],
                            double &z_min,
                            double &z_max) const
    {
        double theta = this->half_angle * D2R;
        double x_abs = abs_max(x_minmax, 2);
        double y_abs = abs_max(y_minmax, 2);
        z_min = 0.0;
        z_max = sqrt(x_abs * x_abs + y_abs * y_abs) / tan(theta);
        return;
    }

    surface_ptr Cone::make_copy() const
    {
        return make_surface<Cone>(*this);
    }

    void Cone::write_json(nlohmann::ordered_json &jnode) const
    {
        SurfaceType type = SurfaceType::CONE;
        jnode["surface_type"] = SurfaceTypeMap.at(type);
        jnode["half_angle"] = this->half_angle;
    }

    Cylinder::Cylinder(const nlohmann::ordered_json &jnode)
        : Surface(SurfaceType::CYLINDER)
    {
        this->radius = jnode.at("radius");
        validate();
    }

    void Cylinder::validate() const
    {
        if (!std::isfinite(radius) || radius <= 0.0)
            throw std::invalid_argument("Cylinder: radius must be positive");
    }

    void Cylinder::bounding_box(const double x_minmax[2],
                                const double y_minmax[2],
                                double &z_min,
                                double &z_max) const
    {
        double r = this->radius;
        // Debug check only. This should be caught upstream before
        // any bounding box calculations are done by a SimulationRunner.
        // See, e.g., native_runner/cylinder_calculator.cpp.
        assert(is_approx(x_minmax[0], -r, 1e-6));
        assert(is_approx(x_minmax[1], r, 1e-6));
        z_min = -r;
        z_max = r;
        return;
    }

    surface_ptr Cylinder::make_copy() const
    {
        return make_surface<Cylinder>(*this);
    }

    void Cylinder::write_json(nlohmann::ordered_json &jnode) const
    {
        SurfaceType type = SurfaceType::CYLINDER;
        jnode["surface_type"] = SurfaceTypeMap.at(type);
        jnode["radius"] = this->radius;
    }

    void Flat::bounding_box(const double x_minmax[2],
                            const double y_minmax[2],
                            double &z_min,
                            double &z_max) const
    {
        z_min = -1e-4;
        z_max = 1e-4;
        return;
    }

    surface_ptr Flat::make_copy() const
    {
        return make_surface<Flat>(*this);
    }

    void Flat::write_json(nlohmann::ordered_json &jnode) const
    {
        SurfaceType type = SurfaceType::FLAT;
        jnode["surface_type"] = SurfaceTypeMap.at(type);
    }

    Parabola::Parabola(const nlohmann::ordered_json &jnode)
        : Surface(SurfaceType::PARABOLA)
    {
        this->focal_length_x = jnode.at("focal_length_x");
        this->focal_length_y = jnode.at("focal_length_y");
        validate();
    }

    void Parabola::validate() const
    {
        if (std::isnan(focal_length_x) || std::isnan(focal_length_y))
            throw std::invalid_argument("Parabola: focal lengths cannot be NaN");
        if (std::isinf(focal_length_x) && std::isinf(focal_length_y))
            throw std::invalid_argument("Parabola: both focal lengths cannot be infinite");
    }

    void Parabola::bounding_box(const double x_minmax[2],
                                const double y_minmax[2],
                                double &z_min,
                                double &z_max) const
    {
        double cx = 0.5 / this->focal_length_x;
        double cy = 0.5 / this->focal_length_y;
        double x_max = abs_max(x_minmax, 2);
        double y_max = abs_max(y_minmax, 2);
        z_max = 0.5 * (cx * x_max * x_max + cy * y_max * y_max);

        if (x_minmax[0] <= 0.0 && 0.0 <= x_minmax[1] &&
            y_minmax[0] <= 0.0 && 0.0 <= y_minmax[1])
        {
            z_min = 0.0;
        }
        else
        {
            double x_min = abs_min(x_minmax, 2);
            double y_min = abs_min(y_minmax, 2);
            z_min = 0.5 * (cx * x_min * x_min + cy * y_min * y_min);
        }

        return;
    }

    surface_ptr Parabola::make_copy() const
    {
        return make_surface<Parabola>(*this);
    }

    void Parabola::write_json(nlohmann::ordered_json &jnode) const
    {
        SurfaceType type = SurfaceType::PARABOLA;
        jnode["surface_type"] = SurfaceTypeMap.at(type);
        jnode["focal_length_x"] = this->focal_length_x;
        jnode["focal_length_y"] = this->focal_length_y;
    }

    Sphere::Sphere(const nlohmann::ordered_json &jnode)
        : Surface(SurfaceType::SPHERE)
    {
        this->vertex_curv = jnode.at("vertex_curv");
        validate();
    }

    void Sphere::validate() const
    {
        if (!std::isfinite(vertex_curv) || vertex_curv <= 0.0)
            throw std::invalid_argument("Sphere: vertex_curv must be positive");
    }

    void Sphere::bounding_box(const double x_minmax[2],
                              const double y_minmax[2],
                              double &z_min,
                              double &z_max) const
    {
        z_min = 0.0;
        double R = 1.0 / this->vertex_curv;
        double x_max = abs_max(x_minmax, 2);
        double y_max = abs_max(y_minmax, 2);
        double rsq = x_max * x_max + y_max * y_max;
        z_max = R > sqrt(rsq) ? R - sqrt(R * R - rsq) : R;
        return;
    }

    surface_ptr Sphere::make_copy() const
    {
        return make_surface<Sphere>(*this);
    }

    void Sphere::write_json(nlohmann::ordered_json &jnode) const
    {
        SurfaceType type = SurfaceType::SPHERE;
        jnode["surface_type"] = SurfaceTypeMap.at(type);
        jnode["vertex_curv"] = this->vertex_curv;
    }

    double Cone::z(double x, double y) const
    {
        return sqrt(x * x + y * y) / tan(half_angle);
    }

    double Cylinder::z(double x, double) const
    {
        // TODO: Fix ? This is really only the top half of the cylinder.
        //       Cylinder breaks the model since it is a multi-valued function: each
        //       x value produces two z values. Returning only the positive root.
        return radius + sqrt(x * x + radius * radius);
    }

    double Parabola::z(double x, double y) const
    {
        // z(x,y) = (cx * x^2 + cy * y^2) / 2
        return x * x / focal_length_x + y * y / focal_length_y;
    }

    double Sphere::z(double x, double y) const
    {
        return vertex_curv * (x * x + y * y) /
               (1 + sqrt(1 - vertex_curv * vertex_curv * (x * x + y * y)));
    }

} // namespace SolTrace::Data

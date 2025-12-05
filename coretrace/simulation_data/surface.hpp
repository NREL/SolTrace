/**
 * @file surface.hpp
 * @brief Optical surface geometry definitions
 *
 * Defines various optical surface types (flat, parabolic, spherical, etc.)
 * and their mathematical representations for ray-surface intersection
 * calculations. Provides the foundation for ray tracing on different
 * geometric surface types used in concentrated solar power systems.
 *
 * @defgroup surfaces Optical Surfaces
 * @{
 */

#ifndef SOLTRACE_SURFACE_H
#define SOLTRACE_SURFACE_H

#include <memory>
#include <vector>
#include <nlohmann/json.hpp>

namespace SolTrace::Data {

enum SurfaceType
{
    CONE,
    CYLINDER,
    FLAT,
    PARABOLA,
    SPHERE,

    HYPER,
    GENERAL_SPENCER_MURTY,
    TORUS,

    SURFACE_UNKNOWN
};

inline const std::map<SurfaceType, std::string> SurfaceTypeMap =
{
    {SurfaceType::CONE, "CONE"},
    {SurfaceType::CYLINDER, "CYLINDER"},
    {SurfaceType::FLAT, "FLAT"},
    {SurfaceType::PARABOLA, "PARABOLA"},
    {SurfaceType::SPHERE, "SPHERE"},
    {SurfaceType::HYPER, "HYPER"},
    {SurfaceType::GENERAL_SPENCER_MURTY, "GENERAL_SPENCER_MURTY"},
    {SurfaceType::TORUS, "TORUS"},
    {SurfaceType::SURFACE_UNKNOWN, "SURFACE_UNKNOWN"}
};

struct Surface
{
public:
    SurfaceType my_type;

    Surface(SurfaceType st) : my_type(st) {}
    virtual ~Surface() {}

    SurfaceType get_type() { return my_type; }

    virtual void write_json(nlohmann::ordered_json& jnode) const = 0;

    virtual double z(double x, double y) const { return 0; }

    inline std::string get_type_string() const {
        switch (my_type) {
        case CONE: return "Cone";
        case CYLINDER: return "Cylinder";
        case FLAT: return "Flat";
        case PARABOLA: return "Parabola";
        case SPHERE: return "Sphere";
        case HYPER: return "Hyper";
        case GENERAL_SPENCER_MURTY: return "General Spencer Murty";
        case TORUS: return "Torus";
        case SURFACE_UNKNOWN: return "Unknown";
        }
        return "Unknown";
    }
};

struct Cone : public Surface
{
    // z(x,y) = sqrt(x^2 + y^2) / tan(theta)
    // where theta = half_angle
    double half_angle;
    Cone(double ha) : Surface(SurfaceType::CONE), half_angle(ha) {}
    Cone(const nlohmann::ordered_json& jnode);
    virtual ~Cone() {}
    virtual void write_json(nlohmann::ordered_json& jnode) const override;

    virtual double z(double x, double y) const override;
};

struct Cylinder : public Surface
{
    // x^2 + (z - r)^2 = r^2
    // where r = radius
    double radius;
    Cylinder(double r) : Surface(SurfaceType::CYLINDER), radius(r)
    {
    }
    Cylinder(const nlohmann::ordered_json& jnode);
    virtual ~Cylinder() {}
    virtual void write_json(nlohmann::ordered_json& jnode) const override;

    virtual double z(double x, double y) const override;
};

struct Flat : public Surface
{
    Flat() : Surface(SurfaceType::FLAT) {}
    Flat(const nlohmann::ordered_json& jnode) : Surface(SurfaceType::FLAT) {};
    virtual ~Flat() {}
    virtual void write_json(nlohmann::ordered_json& jnode) const override;
};

struct Parabola : public Surface
{
    // z(x,y) = (cx * x^2 + cy * y^2) / 2
    // TODO: Assuming that vertex_x_curv gives cx and
    // that vertex_y_curv gives cy
    double focal_length_x;
    double focal_length_y;

    Parabola(double focal_x, double focal_y) : Surface(SurfaceType::PARABOLA),
                                               focal_length_x(focal_x),
                                               focal_length_y(focal_y)
    {
    }
    Parabola(const nlohmann::ordered_json& jnode);
    virtual ~Parabola() {}
    virtual void write_json(nlohmann::ordered_json& jnode) const override;

    virtual double z(double x, double y) const override;
};

struct Sphere : public Surface
{
    // z(x,y) = c(x^2 + y^2) / [1 + sqrt(1 - c^2{x^2 + y^2})]
    // where c = 1/R.
    // TODO: This form seems to be unnecessarily complicated.
    // Could easily just use one of the equations
    // z(x,y) = (1 - sqrt(1 - c^2 (x^2 + y^2))) / c
    //        = R - sqrt(R^2 - (x^2 + y^2))
    // Need to check on this.
    double vertex_curv;

    Sphere(double curv) : Surface(SurfaceType::SPHERE),
                          vertex_curv(curv)
    {
    }
    Sphere(const nlohmann::ordered_json& jnode);
    virtual ~Sphere() {}
    virtual void write_json(nlohmann::ordered_json& jnode) const override;

    virtual double z(double x, double y) const override;
};

// TODO: Add other surface types. Documentation has the following:
// 1. Hyperboloid/Ellipsoid
// 2. Zernike Series
// 3. VSHOT data
// 4. Finite Element data
// 5. General Spencer & Murty Equation
// 6. Polynomial Series (rotationally symmetric)
// 7. Cubic Spline Interpolation (rotationally symmetric)

using surface_ptr = std::shared_ptr<Surface>;

template <typename S, typename... Args>
inline auto make_surface(Args &&...args)
{
    return std::make_shared<S>(std::forward<Args>(args)...);
}

surface_ptr make_surface_from_type(SurfaceType type,
                                   const std::vector<double> &args);

surface_ptr make_surface_from_json(const nlohmann::ordered_json& jnode);

} // namespace SolTrace::Data

/**
 * @}
 */

#endif

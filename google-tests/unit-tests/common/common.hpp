#ifndef SOLTRACE_UNITTESTS_COMMON_H
#define SOLTRACE_UNITTESTS_COMMON_H

#include <element.hpp>
#include <simulation_data_export.hpp>
#include <vector3d.hpp>

// Vectors exactly match component-wise
bool is_identical(const SolTrace::Data::Vector3d &x,
                  const SolTrace::Data::Vector3d &y);
// Each vector component are within `tol` of each other so ||x - y||_\infty <= tol
bool is_identical(const SolTrace::Data::Vector3d &x,
                  const SolTrace::Data::Vector3d &y,
                  double tol);
bool is_identical(const SolTrace::Data::Matrix3d &A,
                  const SolTrace::Data::Matrix3d &B);

// Convenience function for making element with all
// required fields are set. Used when the test does not
// care about the specifics of the element.
SolTrace::Data::element_ptr make_configured_element();

// Helper functions to create surfaces
inline std::shared_ptr<Cylinder> create_cylinder_surface(double radius = 1.0)
{
    auto cylinder = std::make_shared<Cylinder>(radius);
    return cylinder;
}

inline std::shared_ptr<Flat> create_flat_surface()
{
    auto flat = std::make_shared<Flat>();
    return flat;
}

inline std::shared_ptr<Parabola> create_parabola_surface(double focal_length = 1.0)
{
    auto parabola = std::make_shared<Parabola>(focal_length, focal_length);
    return parabola;
}

inline std::shared_ptr<Sphere> create_sphere_surface(double vertex_curvature = 0.1)
{
    auto sphere = std::make_shared<Sphere>(vertex_curvature);
    return sphere;
}

// Helper functions to create apertures
inline std::shared_ptr<Circle> create_circle_aperture(double radius = 1.0)
{
    auto circle = std::make_shared<Circle>(2.0 * radius);
    return circle;
}

// inline std::shared_ptr<Rectangle> create_rect_aperture(double a = 1.0, double b = 1.0)
// {
//     auto rect = std::make_shared<Rectangle>(a, b);
//     return rect;
// }

inline std::shared_ptr<Rectangle> create_rectangle_aperture(double x_length = 2.0,
                                                            double y_length = 2.0)
{
    auto rect = std::make_shared<Rectangle>(x_length, y_length);
    return rect;
}

#endif

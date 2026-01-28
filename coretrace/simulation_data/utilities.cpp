
#include "utilities.hpp"

#include <cassert>
#include <cmath>

#include "constants.hpp"
#include "vector3d.hpp"

namespace SolTrace::Data {

void project_onto_plane(const Vector3d &n,
                        const Vector3d &u,
                        Vector3d &uproj)
{
    double nmag = n.norm();
    double coef = -dot_product(n, u) / (nmag * nmag);
    vector_add(1.0, u, coef, n, uproj);
    return;
}

void project_onto_plane(const Vector3d &n,
                        Vector3d &u)
{
    double nmag = n.norm();
    double coef = -dot_product(n, u) / (nmag * nmag);
    // vector_add(1.0, u, coef, n, uproj);
    vector_add(coef, n, 1.0, u);
    return;
}

void project_onto_vector(const Vector3d &u,
                         const Vector3d &v,
                         Vector3d &vproj)
{
    double umag = u.norm();
    double coef = dot_product(u, v) / (umag * umag);
    vector_add(coef, u, 0.0, vproj);
    return;
}

void project_onto_vector(const Vector3d &u,
                         Vector3d &v)
{
    double umag = u.norm();
    double coef = dot_product(u, v) / (umag * umag);
    vector_add(coef, u, 0.0, v);
    return;
}

void rotate_vector_degrees(const Vector3d &k,
                           const Vector3d &v,
                           double theta,
                           Vector3d &vrot)
{
    theta *= D2R;
    return rotate_vector_radians(k, v, theta, vrot);
}

void rotate_vector_radians(const Vector3d &k,
                           const Vector3d &v,
                           double theta,
                           Vector3d &vrot)
{
    assert(fabs(k.norm() - 1.0) < 1e-12);

    Vector3d scratch;

    cross_product(k, v, scratch);
    vector_add(sin(theta), scratch, cos(theta), v, vrot);

    double coef = (1.0 - cos(theta)) * dot_product(k, v);
    vector_add(coef, k, 1.0, vrot);

    return;
}

void sun_position_vector_degrees(Vector3d &sun_pos,
                                 double azimuth,
                                 double elevation)
{
    return sun_position_vector_radians(sun_pos,
                                       D2R * azimuth,
                                       D2R * elevation);
}

void sun_position_vector_radians(Vector3d &sun_pos,
                                 double azimuth,
                                 double elevation)
{
    double x = sin(azimuth) * cos(elevation);
    double y = cos(azimuth) * cos(elevation);
    double z = sin(elevation);
    sun_pos.set_values(x, y, z);
    return;
}

} // namespace SolTrace::Data

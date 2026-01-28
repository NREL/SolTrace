
#ifndef SOLTRACE_UTILITIES_H
#define SOLTRACE_UTILITIES_H

#include <cstdint>
#include <limits>

#include "vector3d.hpp"

namespace SolTrace::Data {

template <typename R>
R abs_min(const R values[], uint_fast64_t size)
{
    R amin = std::numeric_limits<R>::max();
    for (uint_fast64_t k=0; k < size; ++k)
    {
        amin = std::min(amin, std::abs(values[k]));
    }
    return amin;
}

template <typename R>
R abs_max(const R values[], uint_fast64_t size)
{
    R amax = 0.0;
    for (uint_fast64_t k=0; k < size; ++k)
    {
        amax = std::max(amax, std::abs(values[k]));
    }
    return amax;
}

template <typename R>
bool is_approx(const R &x, const R &y, const R &atol=1e-6)
{
    return (fabs(x - y) < atol);
}

void project_onto_plane(const Vector3d &n,
                        const Vector3d &u,
                        Vector3d &uproj);

void project_onto_plane(const Vector3d &n,
                        Vector3d &u);

void project_onto_vector(const Vector3d &u,
                         const Vector3d &v,
                         Vector3d &vproj);

void project_onto_vector(const Vector3d &u,
                         Vector3d &v);

void rotate_vector_degrees(const Vector3d &k,
                           const Vector3d &v,
                           double theta,
                           Vector3d &vrot);

void rotate_vector_radians(const Vector3d &k,
                           const Vector3d &v,
                           double theta,
                           Vector3d &vrot);

void sun_position_vector_degrees(Vector3d &sun_pos,
                                 double azimuth,
                                 double elevation);

void sun_position_vector_radians(Vector3d &sun_pos,
                                 double azimuth,
                                 double elevation);

} // namespace SolTrace::Data

#endif

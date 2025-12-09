/**
 * @file ray_source.hpp
 * @brief Ray source definitions and solar models
 *
 * Defines ray source properties including sun shape models,
 * solar disk properties, and ray generation parameters. Provides
 * base classes for different types of ray sources used in
 * concentrated solar power simulations.
 */

#ifndef SOLTRACE_RAY_SOURCE_H
#define SOLTRACE_RAY_SOURCE_H

#include <limits>
#include <map>

#include "container.hpp"
#include "datetime.hpp"
#include "vector3d.hpp"

namespace SolTrace::Data {

// Add new enums to SunShapeMap too
enum class SunShape
{
    GAUSSIAN,
    PILLBOX,
    LIMBDARKENED,
    BUIE_CSR,
    USER_DEFINED,
    UNKNOWN
};

inline const std::map<SunShape, std::string> SunShapeMap =
{
    {SunShape::GAUSSIAN, "GAUSSIAN"},
    {SunShape::PILLBOX, "PILLBOX"},
    {SunShape::LIMBDARKENED, "LIMBDARKENED"},
    {SunShape::BUIE_CSR, "BUIE_CSR"},
    {SunShape::USER_DEFINED, "USER_DEFINED"},
    {SunShape::UNKNOWN, "UNKNOWN"}
};

class RaySource
{
public:
    RaySource() {}
    virtual ~RaySource() {}

    virtual const Vector3d &get_position() const = 0;
    virtual Vector3d &get_position() = 0;
    virtual void set_position(const Vector3d &) = 0;
    virtual void set_position(double, double, double) = 0;
    virtual void set_position(const DateTime &, double lat, double long) = 0;
    virtual SunShape get_shape() const = 0;
    virtual void set_shape(SunShape shape, double _sigma, double _half_width, double _csr,
        std::vector<double> _user_angle = {}, std::vector<double> _user_intensity = {}) = 0;

    double get_sigma()
    {
        return this->sigma;
    }
    double get_half_width()
    {
        return this->half_width;
    }
    double get_circumsolar_ratio()
    {
        return this->circumsolar_ratio;
    }
    void get_user_data(std::vector<double> &angle, std::vector<double> &intensity)
    {
        angle = this->user_angle;
        intensity = this->user_intensity;
        return;
    }
    virtual void calculate_buie_parameters(double& kappa, double& gamma) = 0;

protected:
    double sigma = std::numeric_limits<double>::quiet_NaN();
    double half_width = std::numeric_limits<double>::quiet_NaN();
    double circumsolar_ratio = std::numeric_limits<double>::quiet_NaN();
    std::vector<double> user_angle;
    std::vector<double> user_intensity;
};

using ray_source_id = std::int_fast64_t;
// using ray_source_ptr = typename std::shared_ptr<RaySource>;
// using RaySourceContainer = typename std::map<ray_source_id, ray_source_ptr>;
using RaySourceContainer = Container<ray_source_id, RaySource>;
using ray_source_ptr = RaySourceContainer::value_pointer;

template<typename C, typename... Args>
inline auto make_ray_source(Args&&... args)
{
    return RaySourceContainer::make_pointer<C>(std::forward<Args>(args)...);
}

} // namespace SolTrace::Data

#endif

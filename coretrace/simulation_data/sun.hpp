/**
 * @file sun.hpp
 * @brief Solar source modeling and sun shape definitions
 *
 * Defines solar source properties including sun shape models,
 * solar position calculations, and ray generation from solar disk.
 * Includes models for different sun shape distributions and
 * solar tracking calculations for CST systems.
 */

#ifndef SOLTRACE_SUN_H
#define SOLTRACE_SUN_H

#include "ray_source.hpp"
#include "datetime.hpp"
#include "vector3d.hpp"

namespace SolTrace::Data {

class Sun : public RaySource
{
public:
    Sun() { this->my_position.zero(); }
    virtual ~Sun() {}

    virtual const Vector3d &get_position() const
    {
        return this->my_position;
    }
    virtual Vector3d &get_position()
    {
        return this->my_position;
    }
    virtual void set_position(const Vector3d &pos)
    {
        this->my_position = pos;
        return;
    }
    virtual void set_position(double x, double y, double z)
    {
        this->my_position.set_values(x, y, z);
        return;
    }
    virtual void set_position(const DateTime &, double lat, double long) {}
    virtual SunShape get_shape() const
    {
        return this->my_shape;
    }
    virtual void set_shape(SunShape shape,
                           double _sigma,
                           double _half_width,
                           double _csr,        
                           std::vector<double> _user_angle = {},
                           std::vector<double> _user_intensity = {});
    virtual void calculate_buie_parameters(double& kappa, double& gamma);

private:
    void set_gaussian_distribution(double _sigma);
    void set_pillbox_distribution(double _half_width);
    void set_buie_csr_distribution(double _csr);
    void set_user_defined_distribution(std::vector<double> _user_angle,
                                       std::vector<double> _user_intensity);

    SunShape my_shape;
    Vector3d my_position;
};

} // namespace SolTrace::Data

#endif

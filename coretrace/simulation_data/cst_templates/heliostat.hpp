/**
 * @file heliostat.hpp
 * @brief Heliostat CST template for power tower systems
 *
 * Defines the Heliostat class which generates heliostat field elements
 * with configurable facet layouts, canting methods, and tracking capabilities
 * for concentrated solar power tower applications. Supports various canting
 * strategies and geometric configurations.
 *
 * @defgroup cst_templates CST Templates
 * @{
 */

#ifndef SOLTRACE_HELIOSTAT_H
#define SOLTRACE_HELIOSTAT_H

#include <vector>

#include "composite_element.hpp"
#include "single_element.hpp"

namespace SolTrace::Data {

class Heliostat : public CompositeElement
{
public:
    enum CantingType
    {
        UNSET = 0,
        NONE = 1,
        OFF_AXIS = 2,
        ON_AXIS = 3
    };

    Heliostat();
    ~Heliostat();

    void create_geometry();
    void update_geometry(double azimuth, double elevation);

    void set_aperture_size(double size_x, double size_y);
    void set_focal_length(double flen);
    void set_focal_length(double fx, double fy);
    void set_gaps(double gap_x, double gap_y);
    void set_number_panels(uint_fast64_t num_x, uint_fast64_t num_y);
    void set_mirror_optics(const OpticalProperties& optics);
    void set_back_optics(const OpticalProperties& optics);
    void set_optics(const OpticalProperties& mirror_optics, const OpticalProperties& back_optics);
    // void set_onaxis_canting_distance(double dist);
    // void set_offaxis_canting_sun_position(double azimuth, double zenith);
    void set_canting(CantingType ct, double val1, double val2);
    void set_target_position(const Vector3d &pos);
    void set_tracking_limits(double az_lower, double az_upper,
                             double el_lower, double el_upper);

    Vector3d get_elevation_axis() const
    {
        return this->elevation_axis;
    }

    double get_area() const
    {
        return this->heliostat_area;
    }

    virtual void enforce_user_fields_set() const override;

private:
    bool initialized;

    double aperture_size_x;
    double aperture_size_y;
    double focal_length_x;
    double focal_length_y;
    double gap_x;
    double gap_y;

    uint_fast64_t num_panels_x;
    uint_fast64_t num_panels_y;
    OpticalProperties optics_mirror;
    OpticalProperties optics_back;

    CantingType canting_method;
    double onaxis_canting_distance;
    double offaxis_canting_sun_position_azimuth;
    double offaxis_canting_sun_position_zenith;

    double heliostat_area;
    double tracking_azimuth;
    double tracking_elevation;
    Vector3d elevation_axis;
    // Sun position in global coordinates
    Vector3d sun_position;
    // Target position in global coordinates
    Vector3d target_pos;
    bool target_set;

    std::vector<single_element_ptr> facets;
};

} // namespace SolTrace::Data

/**
 * @}
 */

#endif

/**
 * @file optical_properties.hpp
 * @brief Optical properties definitions for materials
 *
 * Defines optical properties (reflectivity, transmissivity, refractive index)
 * and interaction types for optical surfaces and materials. Includes error
 * distribution parameters for modeling surface imperfections and optical errors.
 */

#ifndef SOLTRACE_OPTICAL_PROPERTIES_H
#define SOLTRACE_OPTICAL_PROPERTIES_H

#include "error_distributions.hpp"

namespace SolTrace::Data {

enum InteractionType
{
    REFLECTION,
    REFRACTION
};

struct OpticalProperties
{
    InteractionType my_type;
    DistributionType error_distribution_type;
    double transmitivity;
    double reflectivity;
    double slope_error;
    double specularity_error;
    double refraction_index_front;
    double refraction_index_back;

    OpticalProperties() : my_type(REFLECTION),
                          error_distribution_type(GAUSSIAN),
                          transmitivity(0.0),
                          reflectivity(0.0),
                          slope_error(0.0),
                          specularity_error(0.0),
                          refraction_index_front(0.0),
                          refraction_index_back(0.0)
    {
    }

    OpticalProperties(InteractionType itype,
                      DistributionType dtype,
                      double trans, double refl,
                      double slope_err, double spec_err,
                      double ri_front, double ri_back)
        : my_type(itype),
          error_distribution_type(dtype),
          transmitivity(trans),
          reflectivity(refl),
          slope_error(slope_err),
          specularity_error(spec_err),
          refraction_index_front(ri_front),
          refraction_index_back(ri_back)
    {
    }

    // TODO: What should the error settings be with the below?

    void set_ideal_absorption()
    {
        this->my_type = REFLECTION;
        this->transmitivity = 0.0;
        this->reflectivity = 0.0;
        return;
    }
    void set_ideal_reflection()
    {
        this->my_type = REFLECTION;
        this->transmitivity = 0.0;
        this->reflectivity = 1.0;
        return;
    }
    void set_ideal_transmission()
    {
        this->my_type = REFRACTION;
        this->transmitivity = 1.0;
        this->reflectivity = 0.0;
        return;
    }
    void set_ideal_transmission(double refraction_index_front,
                                double refraction_index_back)
    {
        this->set_ideal_transmission();
        this->refraction_index_front = refraction_index_front;
        this->refraction_index_back = refraction_index_back;
        return;
    }

    // OpticalProperties &operator=(const OpticalProperties &rhs)
    // {
    //     this->my_type = rhs.my_type;
    //     this->transmitivity = rhs.transmitivity;
    //     this->reflectivity = rhs.reflectivity;
    //     this->slope_error = rhs.slope_error;
    //     this->specularity_error = rhs.specularity_error;
    //     this->refraction_index_front = rhs.refraction_index_front;
    //     this->refraction_index_back = rhs.refraction_index_back;
    //     return *this;
    // }
};

} // namespace SolTrace::Data

#endif

#include "sun.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace SolTrace::Data {

void Sun::set_gaussian_distribution(double _sigma)
{
    if (_sigma <= 0.0)
    {
        throw std::invalid_argument("Gaussian distribution sigma must be positive");
    }
    if (std::isnan(_sigma) || std::isinf(_sigma))
    {
        throw std::invalid_argument("Gaussian distribution sigma must be finite");
    }

    sigma = _sigma;
}

void Sun::set_pillbox_distribution(double _half_width)
{
    if (_half_width <= 0.0)
    {
        throw std::invalid_argument("Pillbox distribution half_width must be positive");
    }
    if (std::isnan(_half_width) || std::isinf(_half_width))
    {
        throw std::invalid_argument("Pillbox distribution half_width must be finite");
    }

    half_width = _half_width;
}

void Sun::set_user_defined_distribution(std::vector<double> _user_angle,
                                        std::vector<double> _user_intensity)
{
    if (_user_angle.empty())
    {
        throw std::invalid_argument("User-defined distribution requires non-empty angle vector");
    }
    if (_user_intensity.empty())
    {
        throw std::invalid_argument("User-defined distribution requires non-empty intensity vector");
    }
    if (_user_angle.size() != _user_intensity.size())
    {
        throw std::invalid_argument("User-defined distribution angle and intensity vectors must have same size");
    }

    // Check for valid angle values
    for (const auto &angle : _user_angle)
    {
        if (std::isnan(angle) || std::isinf(angle))
        {
            throw std::invalid_argument("User-defined distribution angles must be finite");
        }
        if (angle < 0.0)
        {
            throw std::invalid_argument("User-defined distribution angles must be non-negative");
        }
    }

    // Check for valid intensity values
    for (const auto &intensity : _user_intensity)
    {
        if (std::isnan(intensity) || std::isinf(intensity))
        {
            throw std::invalid_argument("User-defined distribution intensities must be finite");
        }
        if (intensity < 0.0)
        {
            throw std::invalid_argument("User-defined distribution intensities must be non-negative");
        }
    }

    // Check that angles are in ascending order
    if (!std::is_sorted(_user_angle.begin(), _user_angle.end()))
    {
        throw std::invalid_argument("User-defined distribution angles must be in ascending order");
    }

    user_angle = std::move(_user_angle);
    user_intensity = std::move(_user_intensity);
}

void Sun::set_shape(SunShape shape,
                    double _sigma,
                    double _half_width,
                    std::vector<double> _user_angle,
                    std::vector<double> _user_intensity)
{
    this->my_shape = shape;

    // Clear arguments
    sigma = std::numeric_limits<double>::quiet_NaN();
    half_width = std::numeric_limits<double>::quiet_NaN();
    user_angle.clear();
    user_intensity.clear();

    switch (shape)
    {
    case (SunShape::GAUSSIAN):
        set_gaussian_distribution(_sigma);
        break;
    case (SunShape::PILLBOX):
        set_pillbox_distribution(_half_width);
        break;
    case (SunShape::USER_DEFINED):
        set_user_defined_distribution(std::move(_user_angle), std::move(_user_intensity));
        break;
    default:
        throw std::invalid_argument("Unknown distribution type");
        break;
    }

    return;
}

} // namespace SolTrace::Data

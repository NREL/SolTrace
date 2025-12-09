/**
 * @file error_distributions.hpp
 * @brief Error distribution types for optical modeling
 *
 * Defines various error distribution types (Gaussian, Pillbox, User-defined)
 * used for modeling optical errors in ray tracing simulations. These distributions
 * are applied to surface slopes, specularity errors, and other optical imperfections.
 */

#ifndef SOLTRACE_ERROR_DISTRIBUTIONS_H
#define SOLTRACE_ERROR_DISTRIBUTIONS_H

#include <map>

namespace SolTrace::Data {

// Add to DistributionTypeMap when adding to DistributionType
enum class DistributionType
{
    GAUSSIAN,
    PILLBOX,
    DIFFUSE,
    USER_DEFINED,        // Not supported in legacy but should be added.
    UNKNOWN
};

inline const std::map<DistributionType, std::string> DistributionTypeMap =
{
    {DistributionType::GAUSSIAN, "GAUSSIAN"},
    {DistributionType::PILLBOX, "PILLBOX"},
    {DistributionType::DIFFUSE, "DIFFUSE"},
    {DistributionType::USER_DEFINED, "USER_DEFINED"},
    {DistributionType::UNKNOWN, "UNKNOWN"}
};

} // namespace SolTrace::Data

#endif

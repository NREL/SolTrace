#ifndef SOLTRACE_TRACING_ERRORS_H
#define SOLTRACE_TRACING_ERRORS_H

#include "optical_properties.hpp"

#include "mtrand.hpp"
#include "native_runner_types.hpp"

namespace SolTrace::NativeRunner {

// Perturbs a surface normal by the slope error.
glm::dvec3 ApplySlopeError(MTRand& myrng,
                           const glm::dvec3& CosIn,
                           const SolTrace::Data::OpticalPropertySet& OptProperties,
                           const bool LastHitBackSide);

// Perturbs an incoming ray by the sun shape.
glm::dvec3 ApplySunShape(MTRand& myrng, const glm::dvec3& CosIn, const TSun& Sun);

// Perturbs an outgoing ray by the specularity error, rejecting directions that
// would pass through an opaque reflecting surface.
glm::dvec3 ApplySpecularityError(MTRand& myrng,
                                 const glm::dvec3& CosIn,
                                 const SolTrace::Data::OpticalPropertySet& OptProperties,
                                 const bool LastHitBackSide,
                                 const glm::dvec3& DFXYZ);


} // namespace SolTrace::NativeRunner

#endif

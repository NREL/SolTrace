#ifndef SOLTRACE_TRACING_ERRORS_H
#define SOLTRACE_TRACING_ERRORS_H

#include "optical_properties.hpp"

#include "mtrand.hpp"
#include "native_runner_types.hpp"

namespace SolTrace::NativeRunner {

void SampleSunShape(MTRand& myrng,
                    const glm::dvec3& CosIn,
                    const TSun* Sun,
                    glm::dvec3& CosOut);

void ApplySurfaceError(MTRand& myrng,
                       const glm::dvec3& CosIn,
                       const SolTrace::Data::OpticalPropertySet* OptProperties,
                       const bool LastHitBackSide,
                       const glm::dvec3& DFXYZ,
                       glm::dvec3& CosOut);

void SurfaceNormalErrors(MTRand& myrng,
                         const glm::dvec3& CosIn,
                         const SolTrace::Data::OpticalPropertySet* OptProperties,
                         const bool LastHitBackSide,
                         glm::dvec3& CosOut);


} // namespace SolTrace::NativeRunner

#endif

#ifndef SOLTRACE_TRACING_ERRORS_H
#define SOLTRACE_TRACING_ERRORS_H

#include "optical_properties.hpp"

#include "mtrand.hpp"
#include "native_runner_types.hpp"

namespace SolTrace::NativeRunner {

void Errors(MTRand& myrng,
            const glm::dvec3& CosIn,
            int Source,
            TSun* Sun,
            const SolTrace::Data::OpticalPropertySet* OptProperties,
            const bool LastHitBackSide,
            glm::dvec3& CosOut,
            const glm::dvec3& DFXYZ);

void SurfaceNormalErrors(MTRand& myrng,
                         const glm::dvec3& CosIn,
                         const SolTrace::Data::OpticalPropertySet* OptProperties,
                         const bool LastHitBackSide,
                         glm::dvec3& CosOut);


} // namespace SolTrace::NativeRunner

#endif

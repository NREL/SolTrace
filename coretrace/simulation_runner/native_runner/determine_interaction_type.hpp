#ifndef SOLTRACE_DETERMINE_INTERACTION_TYPE_H
#define SOLTRACE_DETERMINE_INTERACTION_TYPE_H

#include <optical_properties.hpp>
#include <simulation_result.hpp>

#include "mtrand.hpp"
#include "native_runner_types.hpp"
#include "trace_logger.hpp"

namespace SolTrace::NativeRunner
{

    bool determine_interaction_type(
        trace_logger_ptr logger,
        int_fast64_t stage,
        unsigned thread_id,
        MTRand &myrng,
        const SolTrace::Data::OpticalPropertySet* optics,
        glm::dvec3 const& LastDFXYZ,
        glm::dvec3 const& LastCosRaySurfElement,
        bool LastHitBackSide,
        SolTrace::Result::RayEvent &rev);

} // namespace SolTrace::NativeRunner

#endif

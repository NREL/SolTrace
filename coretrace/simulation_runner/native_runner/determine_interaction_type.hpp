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
        const SolTrace::Data::OpticalProperties *optics,
        const double (&LastDFXYZ)[3],
        const double (&LastCosRaySurfElement)[3],
        // bool LastHitBackSide,
        SolTrace::Result::RayEvent &rev);

} // namespace SolTrace::NativeRunner

#endif

#ifndef SOLTRACE_PROCESS_STAGE_HIT_H
#define SOLTRACE_PROCESS_STAGE_HIT_H

#include "mtrand.hpp"
#include "native_runner_types.hpp"
#include "trace_logger.hpp"

namespace SolTrace::NativeRunner {

enum class StageHitOutcome
{
    CONTINUE,
    ABSORBED,
    ERROR
};

// Classifies and applies the optical interaction for a single element hit
// (errors, angular-table lookup, interaction, transform back to stage/global
// coordinates, and ray-data logging). Shared by the native and Embree
// runners so the interaction ordering can't drift between them again.
StageHitOutcome ProcessStageHit(
    trace_logger_ptr logger,
    TSystem* System,
    MTRand& myrng,
    unsigned thread_id,
    int_fast64_t stage_index,
    const tstage_ptr& Stage,
    bool IncludeSunShape,
    bool IncludeErrors,
    uint_fast64_t MultipleHitCount,
    uint_fast64_t LastElementNumber,
    uint_fast64_t LastRayNumber,
    int LastHitBackSide,
    glm::dvec3& LastDFXYZ,
    glm::dvec3& LastCosRaySurfElement,
    glm::dvec3& LastPosRaySurfElement,
    int& ErrorFlag,
    // Outputs
    glm::dvec3& PosRayStage,
    glm::dvec3& CosRayStage,
    glm::dvec3& PosRayGlob,
    glm::dvec3& CosRayGlob);

} // namespace SolTrace::NativeRunner

#endif

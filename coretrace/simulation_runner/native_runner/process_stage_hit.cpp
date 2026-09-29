#include "process_stage_hit.hpp"

#include "determine_interaction_type.hpp"
#include "process_interaction.hpp"
#include "tracing_errors.hpp"

#include <simulation_data_export.hpp>

namespace SolTrace::NativeRunner
{
    using SolTrace::Result::RayEvent;

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
        glm::dvec3& PosRayStage,
        glm::dvec3& CosRayStage,
        glm::dvec3& PosRayGlob,
        glm::dvec3& CosRayGlob)
    {
        glm::dvec3 PosRayOutElement(0.0);
        glm::dvec3 CosRayOutElement(0.0);

        const SolTrace::Data::OpticalPropertySet* optics_set = nullptr;
        RayEvent rev = RayEvent::VIRTUAL;

        if (Stage->Virtual)
        {
            // If stage is virtual, there is no interaction
            PosRayOutElement = LastPosRaySurfElement;
            CosRayOutElement = LastCosRaySurfElement;
        }
        else
        {
            telement_ptr const& optelm =
                Stage->ElementList[LastElementNumber - 1];
            optics_set = &optelm->Optics;

            // Apply sunshape and slope error before classifying the
            // interaction so the angular-table lookup uses the same
            // perturbed ray/normal as the interaction itself (matches the
            // OptiX runner's ordering).
            if (IncludeSunShape && stage_index == 0 && MultipleHitCount == 1)
            {
                LastCosRaySurfElement = ApplySunShape(
                    myrng, LastCosRaySurfElement, System->Sun);
            }

            if (IncludeErrors)
            {
                LastDFXYZ = ApplySlopeError(
                    myrng, LastDFXYZ, *optics_set, LastHitBackSide);
            }

            bool good = determine_interaction_type(
                logger,
                stage_index,
                thread_id,
                myrng,
                optics_set,
                LastDFXYZ,
                LastCosRaySurfElement,
                LastHitBackSide,
                rev);

            if (!good)
            {
                return StageHitOutcome::ERROR;
            }

            if (rev == RayEvent::ABSORB)
            {
                return StageHitOutcome::ABSORBED;
            }
        }

        // Process Interaction
        int_fast64_t k = LastElementNumber - 1;
        ProcessInteraction(myrng,
                           optics_set,
                           LastHitBackSide,
                           IncludeErrors,
                           stage_index,
                           Stage,
                           MultipleHitCount,
                           LastDFXYZ,
                           LastCosRaySurfElement,
                           ErrorFlag,
                           CosRayOutElement,
                           LastPosRaySurfElement,
                           PosRayOutElement);

        // Transform ray back to stage coordinate system
        TransformToReference(PosRayOutElement,
                             CosRayOutElement,
                             Stage->ElementList[k]->Origin,
                             Stage->ElementList[k]->RLocToRef,
                             PosRayStage,
                             CosRayStage);
        TransformToReference(PosRayStage,
                             CosRayStage,
                             Stage->Origin,
                             Stage->RLocToRef,
                             PosRayGlob,
                             CosRayGlob);

        System->RayData.Append(thread_id,
                               PosRayGlob,
                               CosRayGlob,
                               LastElementNumber,
                               stage_index + 1,
                               LastRayNumber,
                               rev);

        return StageHitOutcome::CONTINUE;
    }

} // namespace SolTrace::NativeRunner

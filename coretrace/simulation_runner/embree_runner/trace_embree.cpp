#include "trace_embree.hpp"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <string>
#include <vector>

// Embree header
#include <embree4/rtcore.h>

// SimulationData header
#include <simulation_data_export.hpp>

// SimulationRunner header
#include <simulation_runner.hpp>

// NativeRunner header(s)
#include <determine_interaction_type.hpp>
#include <generate_ray.hpp>
#include <mtrand.hpp>
#include <native_runner_types.hpp>
#include <process_interaction.hpp>
#include <sun_to_primary_stage.hpp>
#include <thread_manager.hpp>
#include <trace_logger.hpp>

#include "embree_helper.hpp"
#include "find_element_hit_embree.hpp"

namespace SolTrace::EmbreeRunner
{
    using SolTrace::NativeRunner::GlobalRay_refactored;
    using SolTrace::NativeRunner::MTRand;
    using SolTrace::NativeRunner::TElement;
    using SolTrace::NativeRunner::telement_ptr;
    using SolTrace::NativeRunner::thread_manager_ptr;
    using SolTrace::NativeRunner::ThreadManager;
    using SolTrace::NativeRunner::trace_logger_ptr;
    using SolTrace::NativeRunner::tstage_ptr;
    using SolTrace::NativeRunner::TSystem;

    using SolTrace::Runner::RunnerStatus;

    using SolTrace::Result::RayEvent;

    RunnerStatus make_embree_scene(trace_logger_ptr logger,
                                   TSystem *System,
                                   RTCDevice &embree_device,
                                   RTCScene &embree_scene,
                                   unsigned nthreads)
    {
        RunnerStatus sts = RunnerStatus::SUCCESS;
        // // Initialize Embree vars
        // RTCDevice embree_device = nullptr;
        // RTCScene embree_scene = nullptr;
        // bool use_shared_embree = false;

        // Make device
        // std::cout << "Making embree device..." << std::endl;
        std::stringstream ss;
        ss << "threads=" << nthreads;
        embree_device = rtcNewDevice(ss.str().c_str());

        // std::cout << "Setting error function..." << std::endl;
        rtcSetDeviceErrorFunction(embree_device, error_function, NULL);

        // Convert st stages into scene
        // std::cout << "Making scene..." << std::endl;
        embree_scene = make_scene(embree_device, *System);

        // std::cout << "Committing scene..." << std::endl;
        rtcCommitScene(embree_scene);

        // Validate bounds
        RTCError err = rtcGetDeviceError(embree_device);
        if (err != RTC_ERROR_NONE)
        {
            // int asdg = 0;
            // return RunnerStatus::ERROR;
            sts = RunnerStatus::ERROR;
            logger->error_log("Error setting up Embree scene");
        }

        return sts;
    }

    RunnerStatus trace_embree(
        thread_manager_ptr manager,
        trace_logger_ptr logger,
        TSystem *System,
        const std::vector<unsigned> &seeds,
        unsigned nthreads,
        uint_fast64_t NumberOfRays,
        uint_fast64_t MaxNumberOfRays,
        bool IncludeSunShape,
        bool IncludeErrors,
        const RTCScene &embree_scene)
    {

        System->RayData.SetUp(nthreads, NumberOfRays);
        System->SunRayCount = 0;

        // Initialize Sun
        Vector3d PosSunStage;
        bool status = SolTrace::NativeRunner::SunToPrimaryStage(
            logger, System, System->StageList[0].get(),
            &System->Sun, PosSunStage.data);

        if (!status)
            return RunnerStatus::ERROR;

        uint_fast64_t rem = NumberOfRays % nthreads;
        uint_fast64_t nrays_per_thread = NumberOfRays / nthreads;
        uint_fast64_t nrays;

        for (unsigned k = 0; k < nthreads; ++k)
        {
            nrays = k < rem ? nrays_per_thread + 1 : nrays_per_thread;
            ThreadManager::future my_future = std::async(
                std::launch::async,
                trace_embree_single_thread,
                k,
                manager,
                logger,
                System,
                seeds[k],
                nrays,
                MaxNumberOfRays / nthreads + 1,
                IncludeSunShape,
                IncludeErrors,
                PosSunStage,
                embree_scene);
            manager->manage(k, std::move(my_future));
        }

        return manager->monitor_until_completion();
    }

    RunnerStatus trace_embree_single_thread(
        unsigned thread_id,
        thread_manager_ptr manager,
        trace_logger_ptr logger,
        TSystem *System,
        unsigned seed,
        uint_fast64_t NumberOfRays,
        uint_fast64_t MaxNumberOfRays,
        bool IncludeSunShape,
        bool IncludeErrors,
        const SolTrace::Data::Vector3d &PosSunStage,
        const RTCScene &embree_scene)
    {
        // std::cout << "Thread " << thread_id << " with seed " << seed
        //           << std::endl;
        // Initialize Internal State Variables
        MTRand myrng(seed);

        uint_fast64_t update_rate = std::min(
            std::max(static_cast<uint_fast64_t>(1), NumberOfRays / 10),
            static_cast<uint_fast64_t>(1000));
        uint_fast64_t update_count = 0;
        double total_work = System->StageList.size() * NumberOfRays;

        uint_fast64_t RayNumber = 1; // Ray Number of current ray
        bool PreviousStageHasRays = false;
        uint_fast64_t LastRayNumberInPreviousStage = NumberOfRays;

        // Define IncomingRays
        std::vector<GlobalRay_refactored> IncomingRays;
        IncomingRays.resize(NumberOfRays);

        // Initialize stage variables
        uint_fast64_t StageDataArrayIndex = 0;
        uint_fast64_t PreviousStageDataArrayIndex = 0;
        uint_fast64_t n_rays_active = NumberOfRays;
        uint_fast64_t sun_ray_count_local = 0;

        // Loop through stages
        for (uint_fast64_t i = 0; i < System->StageList.size(); i++)
        {
            // Check if previous stage has rays
            bool StageHasRays = true;
            if (i > 0 && PreviousStageHasRays == false)
            {
                StageHasRays = false;
            }

            // Get Current Stage
            tstage_ptr Stage = System->StageList[i];

            // Initialize stage variables
            StageDataArrayIndex = 0;
            PreviousStageDataArrayIndex = 0;

            // Loop through rays
            while (StageHasRays)
            {
                // Initialize Global Coordinates
                double PosRayGlob[3] = {0.0, 0.0, 0.0};
                double CosRayGlob[3] = {0.0, 0.0, 0.0};

                // Initialize Stage Coordinates
                double PosRayStage[3] = {0.0, 0.0, 0.0};
                double CosRayStage[3] = {0.0, 0.0, 0.0};

                // Get Ray
                if (i == 0)
                {
                    // Make ray (if first stage)
                    double PosRaySun[3];
                    SolTrace::NativeRunner::GenerateRay(
                        myrng, PosSunStage.data, Stage->Origin,
                        Stage->RLocToRef, &System->Sun,
                        PosRayGlob, CosRayGlob, PosRaySun);
                    System->SunRayCount++;
                }
                else
                {
                    // Get ray from previous stage
                    RayNumber = IncomingRays[StageDataArrayIndex].Num;
                    CopyVec3(PosRayGlob, IncomingRays[StageDataArrayIndex].Pos);
                    CopyVec3(CosRayGlob, IncomingRays[StageDataArrayIndex].Cos);
                    StageDataArrayIndex++;
                }

                // transform the global incoming ray to local stage coordinates
                TransformToLocal(PosRayGlob, CosRayGlob,
                                 Stage->Origin, Stage->RRefToLoc,
                                 PosRayStage, CosRayStage);

                // Initialize internal variables for ray intersection tracing
                bool RayInStage = true;
                bool in_multi_hit_loop = false;
                double LastPosRaySurfElement[3] = {0.0, 0.0, 0.0};
                double LastCosRaySurfElement[3] = {0.0, 0.0, 0.0};
                double LastPosRaySurfStage[3] = {0.0, 0.0, 0.0};
                double LastCosRaySurfStage[3] = {0.0, 0.0, 0.0};
                double LastDFXYZ[3] = {0.0, 0.0, 0.0};
                uint_fast64_t LastElementNumber = 0;
                uint_fast64_t LastRayNumber = 0;
                int ErrorFlag;
                int LastHitBackSide;
                bool StageHit;
                int MultipleHitCount = 0;
                double PosRayOutElement[3] = {0.0, 0.0, 0.0};
                double CosRayOutElement[3] = {0.0, 0.0, 0.0};

                // Start Loop to trace ray until it leaves stage
                bool RayIsAbsorbed = false;
                while (RayInStage)
                {

                    FindElementHit_embree(embree_scene, i, RayNumber,
                                          PosRayGlob, CosRayGlob,
                                          LastPosRaySurfElement, LastCosRaySurfElement,
                                          LastDFXYZ, LastElementNumber, LastRayNumber,
                                          LastPosRaySurfStage, LastCosRaySurfStage,
                                          ErrorFlag, LastHitBackSide, StageHit);

                    // Breakout if ray left stage
                    if (!StageHit)
                    {
                        RayInStage = false;
                        break;
                    }

                    // Increment MultipleHitCount?
                    MultipleHitCount++;

                    if (i == 0 && MultipleHitCount == 1)
                    {
                        // Add ray to Stage RayData
                        auto r = System->RayData.Append(thread_id,
                                                        PosRayGlob,
                                                        CosRayGlob,
                                                        ELEMENT_NULL,
                                                        i + 1,
                                                        LastRayNumber,
                                                        RayEvent::CREATE);

                        if (r == nullptr)
                        {
                            std::stringstream ss;
                            ss << "Failed to record ray data.\n";
                            logger->error_log(ss.str());
                        }
                    }

                    // Get optics and check for absorption
                    const OpticalProperties *optics = 0;
                    RayEvent rev = RayEvent::VIRTUAL;
                    if (Stage->Virtual)
                    {
                        // If stage is virtual, there is no interaction
                        CopyVec3(PosRayOutElement, LastPosRaySurfElement);
                        CopyVec3(CosRayOutElement, LastCosRaySurfElement);
                    }
                    else
                    {
                        telement_ptr optelm =
                            Stage->ElementList[LastElementNumber - 1];
                        if (LastHitBackSide)
                            optics = &optelm->Optics.Back;
                        else
                            optics = &optelm->Optics.Front;

                        bool good =
                            SolTrace::NativeRunner::determine_interaction_type(
                                logger,
                                i,
                                thread_id,
                                myrng,
                                optics,
                                LastDFXYZ,
                                LastCosRaySurfElement,
                                rev);

                        if (!good)
                        {
                            return RunnerStatus::ERROR;
                        }

                        if (rev == RayEvent::ABSORB)
                        {
                            RayIsAbsorbed = true;
                            break;
                        }
                    }

                    // Process Interaction
                    int k = LastElementNumber - 1;
                    SolTrace::NativeRunner::ProcessInteraction(
                        System,
                        myrng,
                        IncludeSunShape,
                        optics,
                        IncludeErrors,
                        i,
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
                                           i + 1,
                                           LastRayNumber,
                                           rev);

                    // Break out if multiple hits are not allowed
                    if (!Stage->MultiHitsPerRay)
                    {
                        StageHit = false;
                        break;
                    }
                    else
                    {
                        in_multi_hit_loop = true;
                    }
                }

                if (MultipleHitCount > 0)
                    ++update_count;

                if (update_count % update_rate == 0)
                {
                    double progress = update_count / total_work;
                    manager->progress_update(thread_id, progress);
                    if (manager->terminate(thread_id))
                        return RunnerStatus::CANCEL;
                }

                // Handle if Ray was absorbed
                if (RayIsAbsorbed)
                {
                    TransformToReference(LastPosRaySurfStage,
                                         LastCosRaySurfStage,
                                         Stage->Origin,
                                         Stage->RLocToRef,
                                         PosRayGlob,
                                         CosRayGlob);

                    System->RayData.Append(thread_id,
                                           PosRayGlob,
                                           CosRayGlob,
                                           LastElementNumber,
                                           i + 1,
                                           LastRayNumber,
                                           RayEvent::ABSORB);

                    n_rays_active--;

                    // ray was fully absorbed
                    if (RayNumber == LastRayNumberInPreviousStage)
                    {
                        PreviousStageHasRays = false;
                        if (PreviousStageDataArrayIndex > 0)
                        {
                            PreviousStageDataArrayIndex--;
                            PreviousStageHasRays = true;
                        }
                        // goto Label_EndStageLoop;
                        break;
                    }
                    else
                    {
                        if (i == 0)
                        {
                            if (RayNumber == NumberOfRays)
                                // goto Label_EndStageLoop;
                                break;
                            else
                                RayNumber++;
                        }

                        // Next ray in loop
                        continue;
                    }
                }

                // Ray has left the stage
                bool FlagMiss = false;
                if (i == 0)
                {
                    if (MultipleHitCount == 0)
                    {
                        // Ray in first stage missed stage entirely
                        // Generate new ray
                        continue;
                    }
                    else
                    {
                        // Ray hit an element, so save it for next stage
                        CopyVec3(IncomingRays[PreviousStageDataArrayIndex].Pos, PosRayGlob);
                        CopyVec3(IncomingRays[PreviousStageDataArrayIndex].Cos, CosRayGlob);
                        IncomingRays[PreviousStageDataArrayIndex].Num = RayNumber;

                        // Is Ray the last in the stage?
                        if (RayNumber == NumberOfRays)
                        {
                            StageHasRays = false;
                            break;
                        }

                        PreviousStageDataArrayIndex++;
                        PreviousStageHasRays = true;

                        // Move on to next ray
                        RayNumber++;
                        continue;
                    }
                }
                else
                {
                    // After the first stage
                    // Ray hit element OR is traced through stage
                    if (Stage->TraceThrough || MultipleHitCount > 0)
                    {
                        // Ray is saved for the next stage
                        CopyVec3(IncomingRays[PreviousStageDataArrayIndex].Pos, PosRayGlob);
                        CopyVec3(IncomingRays[PreviousStageDataArrayIndex].Cos, CosRayGlob);
                        IncomingRays[PreviousStageDataArrayIndex].Num = RayNumber;

                        // Check if ray is last in stage
                        if (RayNumber == LastRayNumberInPreviousStage)
                        {
                            StageHasRays = false;
                            break;
                        }

                        PreviousStageDataArrayIndex++;
                        PreviousStageHasRays = true;

                        if (MultipleHitCount == 0)
                        {
                            FlagMiss = true;
                        }

                        // Go to next ray
                        continue;
                    }
                    // Ray missed stage entirely and is not traced
                    else
                    {
                        FlagMiss = true;
                    }

                    // Handle FlagMiss condition (
                    if (FlagMiss == true)
                    {
                        LastRayNumber = RayNumber;

                        System->RayData.Append(thread_id,
                                               PosRayGlob,
                                               CosRayGlob,
                                               ELEMENT_NULL,
                                               i + 1,
                                               LastRayNumber,
                                               RayEvent::EXIT);

                        n_rays_active--;

                        if (RayNumber == LastRayNumberInPreviousStage)
                        {
                            if (!Stage->TraceThrough)
                            {
                                PreviousStageHasRays = false;
                                if (PreviousStageDataArrayIndex > 0)
                                {
                                    PreviousStageHasRays = true;
                                    PreviousStageDataArrayIndex--; // last ray was previous one
                                }
                            }

                            // Exit stage
                            StageHasRays = false;
                            break;
                        }
                        else
                        {
                            if (i == 0)
                                RayNumber++; // generate new sun ray

                            // Start new ray
                            continue;
                        }
                    }
                }
            }

            // EndStage section...

            // skipping save_st_data logic

            if (!PreviousStageHasRays)
            {
                LastRayNumberInPreviousStage = 0;
                continue; // No rays to carry forward
            }

            if (PreviousStageDataArrayIndex < IncomingRays.size())
            {
                LastRayNumberInPreviousStage = IncomingRays[PreviousStageDataArrayIndex].Num;
                if (LastRayNumberInPreviousStage == 0)
                {
                    // size_t pp = IncomingRays[PreviousStageDataArrayIndex - 1].Num;
                    // System->errlog("LastRayNumberInPreviousStage=0, stage %d, PrevIdx=%d, CurIdx=%d, pp=%d", i + 1,
                    //                PreviousStageDataArrayIndex, StageDataArrayIndex, pp);
                    return RunnerStatus::ERROR;
                }
            }
            else
            {
                // System->errlog("Invalid PreviousStageDataArrayIndex: %u, @ stage %d",
                //                PreviousStageDataArrayIndex, i + 1);
                return RunnerStatus::ERROR;
            }
        }

        // Close out any remaining rays as misses
        unsigned idx = System->StageList.size() - 1;
        tstage_ptr Stage = System->StageList[idx];
        for (uint_fast64_t k = 0; k < n_rays_active; ++k)
        {
            GlobalRay_refactored ray = IncomingRays[k];
            System->RayData.Append(thread_id,
                                   ray.Pos,
                                   ray.Cos,
                                   ELEMENT_NULL,
                                   idx + 1,
                                   ray.Num,
                                   RayEvent::EXIT);
        }

        // System->SunRayCount is atomic so this is thread safe
        System->SunRayCount += sun_ray_count_local;

        return RunnerStatus::SUCCESS;
    }

} // namespace SolTrace::EmbreeRunner

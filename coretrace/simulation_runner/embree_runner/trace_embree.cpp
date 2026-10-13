#include "trace_embree.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
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
#include <generate_ray.hpp>
#include <mtrand.hpp>
#include <native_runner_types.hpp>
#include <process_stage_hit.hpp>
#include <sun_to_primary_stage.hpp>
#include <thread_manager.hpp>
#include <trace_logger.hpp>

#include "embree_helper.hpp"
#include "find_element_hit_embree.hpp"
#include "ftz_daz.hpp"

namespace SolTrace::EmbreeRunner
{
    using SolTrace::NativeRunner::GlobalRay_refactored;
    using SolTrace::NativeRunner::MTRand;
    using SolTrace::NativeRunner::StageHitOutcome;
    using SolTrace::NativeRunner::TElement;
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

        // // Make device
        // // std::cout << "Making embree device..." << std::endl;
        // std::stringstream ss;
        // ss << "threads=" << nthreads;
        // embree_device = rtcNewDevice(ss.str().c_str());
        
        // TODO: Need to test this on largest scenes we expect to
        // trace. Adding more threads does not significantly help
        // when there are only ~6000 elements.
        embree_device = rtcNewDevice("threads=1");

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

        // using Clock = std::chrono::steady_clock;
        // using Seconds = std::chrono::duration<double>;

        // auto t_start = Clock::now();

        System->RayData.SetUp(nthreads, NumberOfRays);
        System->SunRayCount = 0;

        // auto t_after_setup = Clock::now();

        // Initialize Sun
        glm::dvec3 PosSunStage;
        bool status = SolTrace::NativeRunner::SunToPrimaryStage(logger,
                                                                System,
                                                                System->StageList[0].get(),
                                                                &System->Sun,
                                                                PosSunStage);

        if (!status)
            return RunnerStatus::ERROR;

        // auto t_after_sun_init = Clock::now();

        uint_fast64_t rem = NumberOfRays % nthreads;
        uint_fast64_t nrays_per_thread = NumberOfRays / nthreads;
        uint_fast64_t nrays;

        for (unsigned k = 0; k < nthreads; ++k)
        {
            nrays = k < rem ? nrays_per_thread + 1 : nrays_per_thread;
            const uint_fast64_t ray_index_offset = k * nrays_per_thread + std::min(static_cast<uint_fast64_t>(k), rem);
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
                ray_index_offset,
                IncludeSunShape,
                IncludeErrors,
                PosSunStage,
                embree_scene);
            manager->manage(k, std::move(my_future));
        }

        // auto t_after_launch = Clock::now();

        RunnerStatus result = manager->monitor_until_completion();

        // auto t_done = Clock::now();

        // double s_setup     = Seconds(t_after_setup    - t_start).count();
        // double s_sun_init  = Seconds(t_after_sun_init - t_after_setup).count();
        // double s_launch    = Seconds(t_after_launch   - t_after_sun_init).count();
        // double s_ray_trace = Seconds(t_done           - t_after_launch).count();
        // double s_total     = Seconds(t_done           - t_start).count();

        // std::cout << "[trace_embree] timing (nthreads=" << nthreads
        //           << ", rays=" << NumberOfRays << ")\n"
        //           << "  ray_data_setup : " << s_setup     << " s\n"
        //           << "  sun_init       : " << s_sun_init  << " s\n"
        //           << "  thread_launch  : " << s_launch    << " s\n"
        //           << "  ray_trace      : " << s_ray_trace << " s\n"
        //           << "  total          : " << s_total     << " s\n";

        return result;
    }

    RunnerStatus trace_embree_single_thread(
        unsigned thread_id,
        thread_manager_ptr manager,
        trace_logger_ptr logger,
        TSystem *System,
        unsigned seed,
        uint_fast64_t NumberOfRays,
        uint_fast64_t MaxNumberOfRays,
        uint_fast64_t ray_index_offset,
        bool IncludeSunShape,
        bool IncludeErrors,
        const glm::dvec3 &PosSunStage,
        const RTCScene &embree_scene)
    {
        // std::cout << "Thread " << thread_id << " with seed " << seed
        //           << std::endl;
        // Set flush-to-zero and denormals-are-zero for this thread to avoid
        // slow denormal handling in the FPU (recommended by Embree docs).
        SOLTRACE_SET_FTZ_DAZ();

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
            tstage_ptr const& Stage = System->StageList[i];

            // Initialize stage variables
            StageDataArrayIndex = 0;
            PreviousStageDataArrayIndex = 0;
            int ErrorFlag = 0;

            // Loop through rays
            while (StageHasRays)
            {
                // Initialize Global Coordinates
                glm::dvec3 PosRayGlob(0.0);
                glm::dvec3 CosRayGlob(0.0);

                // Initialize Stage Coordinates
                glm::dvec3 PosRayStage(0.0);
                glm::dvec3 CosRayStage(0.0);


                // Get Ray
                if (i == 0)
                {
                    const uint_fast64_t sample_index = ray_index_offset + sun_ray_count_local + 1;

                    // Make ray (if first stage)
                    glm::dvec3 PosRaySun;
                    SolTrace::NativeRunner::GenerateRay(
                        myrng, PosSunStage, Stage->Origin,
                        Stage->RLocToRef, &System->Sun,
                        sample_index,
                        PosRayGlob, CosRayGlob, PosRaySun,
                        ErrorFlag);

                    if (ErrorFlag != 0)
                    {
                        return RunnerStatus::ERROR;
                    }

                    sun_ray_count_local++;
                    
                }
                else
                {
                    // Get ray from previous stage
                    RayNumber = IncomingRays[StageDataArrayIndex].Num;
                    PosRayGlob = IncomingRays[StageDataArrayIndex].Pos;
                    CosRayGlob = IncomingRays[StageDataArrayIndex].Cos;
                    StageDataArrayIndex++;
                }

                // transform the global incoming ray to local stage coordinates
                {
                    // auto _t0 = Clock::now();
                    TransformToLocal(PosRayGlob, CosRayGlob,
                                     Stage->Origin, Stage->RRefToLoc,
                                     PosRayStage, CosRayStage);
                    // t_transform_to_local += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    //     Clock::now() - _t0).count();
                }

                // Initialize internal variables for ray intersection tracing
                bool RayInStage = true;
                bool in_multi_hit_loop = false;

                glm::dvec3 LastPosRaySurfElement(0.0);
                glm::dvec3 LastCosRaySurfElement(0.0);
                glm::dvec3 LastPosRaySurfStage(0.0);
                glm::dvec3 LastCosRaySurfStage(0.0);
                glm::dvec3 LastDFXYZ(0.0);

                uint_fast64_t LastElementNumber = 0;
                uint_fast64_t LastRayNumber = 0;
                int LastHitBackSide;
                bool StageHit;
                int MultipleHitCount = 0;

                // Start Loop to trace ray until it leaves stage
                bool RayIsAbsorbed = false;
                while (RayInStage)
                {

                    {
                        // auto _t0 = Clock::now();
                        FindElementHit_embree(embree_scene, i, RayNumber,
                                              PosRayGlob, CosRayGlob,
                                              LastPosRaySurfElement, LastCosRaySurfElement,
                                              LastDFXYZ, LastElementNumber, LastRayNumber,
                                              LastPosRaySurfStage, LastCosRaySurfStage,
                                              ErrorFlag, LastHitBackSide, StageHit);
                        // t_find_element_hit += std::chrono::duration_cast<std::chrono::nanoseconds>(
                        //     Clock::now() - _t0).count();
                        // ++n_find_element_hit;
                    }

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
                        // auto _t0_append = Clock::now();
                        auto r = System->RayData.Append(thread_id,
                                                        PosRayGlob,
                                                        CosRayGlob,
                                                        ELEMENT_NULL,
                                                        i + 1,
                                                        LastRayNumber,
                                                        RayEvent::CREATE);
                        // t_ray_data_append += std::chrono::duration_cast<std::chrono::nanoseconds>(
                        //     Clock::now() - _t0_append).count();
                        // ++n_ray_data_append;

                        if (r == nullptr)
                        {
                            std::stringstream ss;
                            ss << "Failed to record ray data.\n";
                            logger->error_log(ss.str());
                        }
                    }

                    // Get optics, apply errors, classify, and interact -
                    // shared with the native runner to keep ordering in sync.
                    StageHitOutcome outcome = SolTrace::NativeRunner::ProcessStageHit(
                        logger,
                        System,
                        myrng,
                        thread_id,
                        i,
                        Stage,
                        IncludeSunShape,
                        IncludeErrors,
                        MultipleHitCount,
                        LastElementNumber,
                        LastRayNumber,
                        LastHitBackSide,
                        LastDFXYZ,
                        LastCosRaySurfElement,
                        LastPosRaySurfElement,
                        ErrorFlag,
                        PosRayStage,
                        CosRayStage,
                        PosRayGlob,
                        CosRayGlob);

                    if (outcome == StageHitOutcome::ERROR)
                    {
                        return RunnerStatus::ERROR;
                    }

                    if (outcome == StageHitOutcome::ABSORBED)
                    {
                        RayIsAbsorbed = true;
                        break;
                    }

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
                    // auto _t0 = Clock::now();
                    double progress = update_count / total_work;
                    manager->progress_update(thread_id, progress);
                    bool should_cancel = manager->terminate(thread_id);
                    // t_progress_update += std::chrono::duration_cast<std::chrono::nanoseconds>(
                    //     Clock::now() - _t0).count();
                    if (should_cancel)
                    {
                        // write_timing();
                        return RunnerStatus::CANCEL;
                    }
                }

                // Handle if Ray was absorbed
                if (RayIsAbsorbed)
                {
                    {
                        // auto _t0 = Clock::now();
                        TransformToReference(LastPosRaySurfStage,
                                             LastCosRaySurfStage,
                                             Stage->Origin,
                                             Stage->RLocToRef,
                                             PosRayGlob,
                                             CosRayGlob);
                        // t_transform_to_reference += std::chrono::duration_cast<std::chrono::nanoseconds>(
                        //     Clock::now() - _t0).count();
                    }

                    {
                        // auto _t0 = Clock::now();
                        System->RayData.Append(thread_id,
                                               PosRayGlob,
                                               CosRayGlob,
                                               LastElementNumber,
                                               i + 1,
                                               LastRayNumber,
                                               RayEvent::ABSORB);
                        // t_ray_data_append += std::chrono::duration_cast<std::chrono::nanoseconds>(
                        //     Clock::now() - _t0).count();
                        // ++n_ray_data_append;
                    }

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
                        IncomingRays[PreviousStageDataArrayIndex].Pos = PosRayGlob;
                        IncomingRays[PreviousStageDataArrayIndex].Cos = CosRayGlob;
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
                        IncomingRays[PreviousStageDataArrayIndex].Pos = PosRayGlob;
                        IncomingRays[PreviousStageDataArrayIndex].Cos = CosRayGlob;
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

                        {
                            // auto _t0 = Clock::now();
                            System->RayData.Append(thread_id,
                                                   PosRayGlob,
                                                   CosRayGlob,
                                                   ELEMENT_NULL,
                                                   i + 1,
                                                   LastRayNumber,
                                                   RayEvent::EXIT);
                            // t_ray_data_append += std::chrono::duration_cast<std::chrono::nanoseconds>(
                            //     Clock::now() - _t0).count();
                            // ++n_ray_data_append;
                        }

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
                    // write_timing();
                    return RunnerStatus::ERROR;
                }
            }
            else
            {
                // System->errlog("Invalid PreviousStageDataArrayIndex: %u, @ stage %d",
                //                PreviousStageDataArrayIndex, i + 1);
                // write_timing();
                return RunnerStatus::ERROR;
            }
        }

        // Close out any remaining rays as misses
        unsigned idx = System->StageList.size() - 1;
        tstage_ptr const& Stage = System->StageList[idx];
        for (uint_fast64_t k = 0; k < n_rays_active; ++k)
        {
            GlobalRay_refactored ray = IncomingRays[k];
            {
                // auto _t0 = Clock::now();
                System->RayData.Append(thread_id,
                                       ray.Pos,
                                       ray.Cos,
                                       ELEMENT_NULL,
                                       idx + 1,
                                       ray.Num,
                                       RayEvent::EXIT);
                // t_ray_data_append += std::chrono::duration_cast<std::chrono::nanoseconds>(
                //     Clock::now() - _t0).count();
                // ++n_ray_data_append;
            }
        }

        // System->SunRayCount is atomic so this is thread safe
        System->SunRayCount += sun_ray_count_local;

        if (manager->terminate(thread_id))
            return RunnerStatus::CANCEL;

        return RunnerStatus::SUCCESS;
    }

} // namespace SolTrace::EmbreeRunner

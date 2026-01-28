#ifndef SOLTRACE_TRACE_H
#define SOLTRACE_TRACE_H

#include <cstdint>
#include <vector>

#include <simulation_data.hpp>
#include <simulation_runner.hpp>

#include "mtrand.hpp"
#include "native_runner_types.hpp"
#include "thread_manager.hpp"
#include "trace_logger.hpp"
#include "treemesh.hpp"

namespace SolTrace::NativeRunner
{
    SolTrace::Runner::RunnerStatus trace_native(
        thread_manager_ptr manager,
        trace_logger_ptr logger,
        TSystem *System,
        const std::vector<unsigned int> &seeds,
        unsigned nthreads,
        uint_fast64_t NumberOfRays,
        uint_fast64_t MaxNumberOfRays,
        bool IncludeSunShape,
        bool IncludeErrors,
        bool AsPowerTower);

    SolTrace::Runner::RunnerStatus trace_single_thread(
        unsigned thread_id,
        thread_manager_ptr manager,
        trace_logger_ptr logger,
        TSystem *System,
        unsigned seed,
        uint_fast64_t NumberOfRays,
        uint_fast64_t MaxNumberOfRays,
        bool IncludeSunShape,
        bool IncludeErrors,
        bool AsPowerTower,
        const SolTrace::Data::Vector3d &PosSunStage,
        st_hash_tree *sun_hash,
        st_hash_tree *rec_hash,
        const SolTrace::Data::Vector3d &reccm_helio);

    struct ThreadInfo
    {
        thread_manager_ptr manager;
        trace_logger_ptr logger;
        TSystem *System;
        // unsigned int seed;
        uint_fast64_t NumberOfRays;
        uint_fast64_t MaxNumberOfRays;
        bool IncludeSunShape;
        bool IncludeErrors;
        bool AsPowerTower;
        SolTrace::Data::Vector3d PosSunStage;
        st_hash_tree *sun_hash;
        st_hash_tree *rec_hash;
        SolTrace::Data::Vector3d reccm_helio;
    };

    // Hack to get around stupid compiler issue
    inline SolTrace::Runner::RunnerStatus trace_single_compact(
        unsigned thread_id,
        unsigned seed,
        ThreadInfo info)
    {
        return trace_single_thread(
            thread_id,
            info.manager,
            info.logger,
            info.System,
            seed,
            info.NumberOfRays,
            info.MaxNumberOfRays,
            info.IncludeSunShape,
            info.IncludeErrors,
            info.AsPowerTower,
            info.PosSunStage,
            info.sun_hash,
            info.rec_hash,
            info.reccm_helio);
    }

} // namespace SolTrace::NativeRunner

#endif

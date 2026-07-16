#pragma once

#include <cstddef>
#include <cstdio>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "core/soltrace_state.h" // SoltraceState
#include "core/vec3d.h"          // Vec3d
#include "core/timer.h"
#include "core/CspElement.h"  // CspElement
#include "core/Surface.h"     // Surface and derived classes
#include "shaders/Soltrace.h" // HitRecord, HitType
#include "ray_utils.h"        // CompactionScratch

#include "../../../../../simulation_data/simulation_data_export.hpp"

namespace OptixCSP
{

    class GeometryManager;
    class pipelineManager;
    class dataManager;
    class CspElement;
    class Vec3d;
    class Surface;

    static constexpr SolTrace::Data::SunShape kSupportedSunshapes[] = {
        SolTrace::Data::SunShape::GAUSSIAN,
        SolTrace::Data::SunShape::PILLBOX,
        SolTrace::Data::SunShape::BUIE_CSR,
        SolTrace::Data::SunShape::LIMBDARKENED,
        SolTrace::Data::SunShape::USER_DEFINED};

    class SolTraceSystem
    {

    public:
        SolTraceSystem();
        ~SolTraceSystem() noexcept;

        /// Call to this function mark the completion of the simulation setup
        void initialize();

        /// Execute the ray tracing simulation
        void run();

        /// Update launch params
        void update();

        // Get all hit points
        void get_hp_output(std::vector<float4> &hp_vec, std::vector<uint_fast64_t> &raynumber_vec, std::vector<int32_t> &element_id_vec,
                           std::vector<uint8_t> &hit_type_vec);

        /// Explicit cleanup
        void clean_up();

        // Reset sys
        void reset();

        void set_verbose(bool verbose); // Set verbosity for debugging
        bool is_verbose() const { return m_verbose; }
        /// <summary>
        /// set the number of rays launched
        /// </summary>
        /// <param name="numSunPoints"></param>
        void set_number_of_rays(uint_fast64_t nrays, uint_fast64_t maxrays)
        {
            m_number_of_rays = nrays;
            m_max_number_of_rays = maxrays;
        }

        /// Set the number of rays launched per iteration.
        /// Use 0 (default) to let determine_batch_size() automatically compute a
        /// batch size that fits the ray-data buffers in available GPU memory.
        /// Throws std::out_of_range if batch_size exceeds the maximum int value.
        void set_batch_size(uint_fast64_t batch_size)
        {
            if (batch_size > static_cast<uint_fast64_t>(std::numeric_limits<int>::max()))
                throw std::out_of_range("batch_size exceeds std::numeric_limits<int>::max()");
            m_batch_size = batch_size;
        }
        uint_fast64_t get_batch_size() const { return m_batch_size; }

        /// Set the maximum ray interaction depth. Must be called before initialize().
        /// Values are clamped to [2, 255]. Defaults to DEFAULT_MAX_TRACE_DEPTH = 5.
        void set_max_ray_depth(uint64_t depth);
        uint8_t get_max_ray_depth() const { return m_max_ray_depth; }

        void set_sun(SolTrace::Data::Sun *sun) { m_sun = sun; }

        void set_seed(uint64_t seed) { m_seed = seed; } // Set sun seed

        void set_optical_errors(bool include_optical_errors)
        {
            m_optical_errors = include_optical_errors;
        }

        /// <summary>
        /// add element
        /// /// </summary>
        void add_element(std::shared_ptr<CspElement> element);

        double get_time_trace();
        double get_time_setup();

        /// Print a formatted summary of all timing information collected during
        /// the last initialize() and run() calls.
        void print_timing() const;

        void print_launch_params();

        /// <summary>
        /// Return sun plane area (parallelogram spanned by sun_v0..sun_v3) in world units.
        /// Computed as |(v0 - v1) x (v1 - v2)|.
        /// </summary>
        double get_sun_plane_area() const;

        uint_fast64_t get_N_sun_rays() const { return m_n_sun_rays; }

        /// Returns the number of run() iterations executed during the last run() call.
        uint64_t get_N_run_iterations() const { return m_n_run_iterations; }

        /// Returns the compacted hit records (CREATE + hits, misses excluded).
        const std::vector<HitRecord> &get_hit_records() const { return m_hit_records; }

        /// Returns the number of rays that hit at least one element.
        uint_fast64_t get_N_hit_rays() const { return m_n_hit_rays; }

        /// Returns the number of rays terminated by max depth (excludes absorption at max depth).
        uint_fast64_t get_N_depth_exceeded_rays() const { return m_n_depth_exceeded_rays; }
        void set_sun_shape_errors(bool flag) { this->m_include_sun_shape_errors = flag; }

        /// Enable or disable trimming excess rays at the end of run() so that
        /// exactly m_number_of_rays hit rays are returned.  Enabled by default.
        void set_trim_excess_rays(bool trim) { m_trim_excess_rays = trim; }
        bool get_trim_excess_rays() const { return m_trim_excess_rays; }

    private:
        // m_verbose and m_state must be declared before the shared_ptr managers so
        // that they are initialized first (C++ initializes members in declaration order).
        // GeometryManager and pipelineManager store references/copies of these at
        // construction time, so they must be valid when the shared_ptrs are built.
        bool m_verbose = false;
        OptixCSP::SoltraceState m_state;

        std::shared_ptr<GeometryManager> geometry_manager;
        std::shared_ptr<pipelineManager> pipeline_manager;
        std::shared_ptr<dataManager> data_manager;

        uint_fast64_t m_number_of_rays;
        uint_fast64_t m_max_number_of_rays;
        uint_fast64_t m_batch_size = 0; // 0 means auto-size: determine_batch_size() calls automatic_batch_size()
        uint8_t m_max_ray_depth = DEFAULT_MAX_TRACE_DEPTH;
        uint_fast64_t m_n_depth_exceeded_rays = 0; // rays stopped by max depth, not absorption

        // Sun
        // OptixCSP::Vec3d m_sun_vector;
        // double m_sun_angle;

        SolTrace::Data::Sun *m_sun;
        bool m_include_sun_shape_errors = false;
        bool m_trim_excess_rays = true;

        uint64_t m_seed = 123456ULL;
        bool m_optical_errors = false;

        // Results

        // Compacted hit records: one contiguous array of HitRecord.
        // Each ray group starts with a HIT_CREATE record followed by its hits.
        // Rays that produced no hits (CREATE-only) are excluded.
        std::vector<HitRecord> m_hit_records;

        // Global ray index (ray_offset + local_index) for each logical hit ray in m_hit_records.
        // Parallel to the logical rays (not records): m_hit_ray_ids.size() == m_n_hit_rays.
        std::vector<uint64_t> m_hit_ray_ids;

        // Count of rays that produced at least one non-CREATE hit.
        uint_fast64_t m_n_hit_rays = 0;

        // Total rays generated (launched from the sun plane) across all run() iterations.
        uint_fast64_t m_n_sun_rays = 0;

        // Current allocated device launch buffer sizes.
        size_t m_hit_buffer_size_allocated = 0;
        size_t m_sun_dir_buffer_size_allocated = 0;

        // Pre-allocated device scratch buffers for GPU stream compaction.
        CompactionScratch m_compaction_scratch;
        CompactionTimings m_compaction_timings;

        std::vector<std::shared_ptr<CspElement>> m_element_list;
        void create_shader_binding_table();
        void allocate_device_buffers();
        void setup_device_buffer();
        // GPU-side compaction: count hits, compact buffer on device, copy result to m_hit_records.
        // Increments m_n_hit_rays by the number of newly collected hit rays.
        void get_buffer_results();
        /// Computes the maximum rays-per-batch that fit in 80 % of current free
        /// GPU memory, accounting for all per-ray device buffers and compaction
        /// scratch. Returns 0 if memory cannot be queried.
        uint_fast64_t automatic_batch_size() const;
        /// Returns the effective batch size for a run() call.
        /// If m_batch_size > 0 the user-supplied value is used as-is.
        /// Otherwise automatic_batch_size() is called and the result is capped
        /// at m_number_of_rays.
        uint_fast64_t determine_batch_size() const;

        Timer m_timer_setup;
        Timer m_timer_trace;

        // initialize() sub-timers
        Timer m_timer_aabb;
        Timer m_timer_geometry;
        Timer m_timer_pipeline;
        Timer m_timer_sbt;

        // run() sub-timers
        Timer m_timer_setup_buffer;
        Timer m_timer_optix_launch;
        Timer m_timer_collect_results;
        uint64_t m_n_run_iterations;

        // memory usage
        size_t m_mem_free_before;     ///< Free GPU memory at the start of initialize(), before any setup allocations.
        size_t m_mem_free_post_setup; ///< Free GPU memory at the end of initialize(), after all setup allocations (BVH,
                                      ///  pipeline, SBT, geometry/material arrays). Used as the baseline in
                                      ///  automatic_batch_size() so batch sizing is stable across run() calls.
        size_t m_mem_free_after;      ///< Free GPU memory sampled during run() for per-launch memory reporting.
    };
}

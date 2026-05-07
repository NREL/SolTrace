#pragma once

#include <string>
#include <vector>
#include <memory>                 
#include <cstddef>                
#include <cstdio>                 



#include "core/soltrace_state.h" // SoltraceState
#include "core/vec3d.h"      // Vec3d
#include "core/timer.h"
#include "core/CspElement.h" // CspElement
#include "core/Surface.h"    // Surface and derived classes

#include "../../../../../simulation_data/simulation_data_export.hpp"

namespace OptixCSP {

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
            SolTrace::Data::SunShape::USER_DEFINED
    };

    class SolTraceSystem {

    public:
        
        SolTraceSystem();
        ~SolTraceSystem();

        /// Call to this function mark the completion of the simulation setup
        void initialize();

        /// Execute the ray tracing simulation
        void run();

        /// Update launch params
        void update();

        // Get all hit points
        void get_hp_output(std::vector<float4>& hp_vec, std::vector<uint_fast64_t>& raynumber_vec, std::vector<int32_t>& element_id_vec,
            std::vector<uint8_t>& hit_type_vec);

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

        void set_sun(SolTrace::Data::Sun* sun) { m_sun = sun; }

        void set_seed(uint64_t seed) { m_seed = seed; }  // Set sun seed

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

        void print_launch_params();

        /// <summary>
        /// Return sun plane area (parallelogram spanned by sun_v0..sun_v3) in world units.
        /// Computed as |(v0 - v1) x (v1 - v2)|.
        /// </summary>
        double get_sun_plane_area() const;

        uint_fast64_t get_N_sun_rays() 
        { 
            if (m_sunraynumber_vec.empty())
                return 0;
            return m_sunraynumber_vec.back(); 
        }

        std::vector<uint_fast64_t> get_sunraynumber_vec() const { return m_sunraynumber_vec; }
        void set_sun_shape_errors(bool flag) { this->m_include_sun_shape_errors = flag; }

        

    private:

        std::shared_ptr<GeometryManager> geometry_manager;
        std::shared_ptr<pipelineManager> pipeline_manager;
        std::shared_ptr<dataManager>     data_manager;

        uint_fast64_t m_number_of_rays;
        uint_fast64_t m_max_number_of_rays;

        bool m_verbose;

        // Sun
        //OptixCSP::Vec3d m_sun_vector;
        //double m_sun_angle;
        
        SolTrace::Data::Sun* m_sun;
        bool m_include_sun_shape_errors = false;


        uint64_t m_seed = 123456ULL;
        bool m_optical_errors;
        OptixCSP::SoltraceState m_state;

        // Results

        // Contains information on rays that hit objects
        std::vector<float4> m_hp_vec;
        std::vector<uint_fast64_t> m_raynumber_vec;
        std::vector<int32_t> m_element_id_vec;
        std::vector<uint8_t> m_hit_type_vec;
        std::vector<uint_fast64_t> m_sunraynumber_vec;    // This is ID of hit rays out of all generated rays

        // Reused host-side scratch buffers for copying launch results back from device.
        std::vector<float4> m_hp_output_buffer_host;
        std::vector<int32_t> m_element_id_buffer_host;
        std::vector<uint8_t> m_hit_type_buffer_host;

        // Current allocated device launch buffer sizes.
        size_t m_hit_point_buffer_size_allocated = 0;
        size_t m_element_id_buffer_size_allocated = 0;
        size_t m_hit_type_buffer_size_allocated = 0;
        size_t m_sun_dir_buffer_size_allocated = 0;

        std::vector<std::shared_ptr<CspElement>> m_element_list;
        void create_shader_binding_table();
        void setup_device_buffer();
        void get_buffer_results(std::vector<float4>& hp_vec, std::vector<uint_fast64_t>& raynumber_vec, 
            std::vector<int32_t>& element_id_vec, std::vector<uint8_t>& hit_type_vec, 
            std::vector<uint_fast64_t>& sunraynumber_vec);

        Timer m_timer_setup;
        Timer m_timer_trace;

        // memory usage
        size_t m_mem_free_before;
        size_t m_mem_free_after;


    };
}

#pragma once

#include "GeometryDataST.h"
#include "MaterialDataST.h"
#include "soltrace_constants.h"

#include <cstdint>
#include <vector_types.h>
#include <optix.h>
#include <curand_kernel.h>

namespace OptixCSP{

    const unsigned int NUM_ATTRIBUTE_VALUES = 4u;
    const unsigned int NUM_PAYLOAD_VALUES   = 2u;
    // NOTE: Maximum number of ray interactions in tracing with the geometry is
    // DEFAULT_MAX_TRACE_DEPTH - 1 (so currently 4). See the end of the function
    // __closesthit__element in materials.cu. Note the type. Limited to 255.
    const uint8_t DEFAULT_MAX_TRACE_DEPTH = 5u;

    struct HitGroupData
    {
        MaterialData material_data;
    };

    struct HitRecord {
        float4 hit_point;
        int32_t element_id;
        uint8_t hit_type;
    };

    enum RayType
    {
        RAY_TYPE_RADIANCE = 0,
        RAY_TYPE_COUNT = 1         // not using occlusion/shadow rays atm
    };

    enum OpticalEntityType : unsigned int {
        RECTANGLE_FLAT          = 0,
        RECTANGLE_PARABOLIC     = 1,
        CYLINDRICAL             = 2,
        TRIANGLE_FLAT           = 3,
        QUADRILATERAL_FLAT      = 4,
        CIRCLE_FLAT             = 5,
        HEXAGON_FLAT            = 6,
        ANNULUS_FLAT            = 7,
        CIRCLE_PARABOLIC        = 8,
        HEXAGON_PARABOLIC       = 9,
        TRIANGLE_PARABOLIC      = 10,
        ANNULUS_PARABOLIC       = 11,
        QUADRILATERAL_PARABOLIC = 12,
        RECTANGLE_SPHERICAL     = 13,
        CIRCLE_SPHERICAL        = 14,
        HEXAGON_SPHERICAL       = 15,
        ANNULUS_SPHERICAL       = 16,
        TRIANGLE_SPHERICAL      = 17,
        QUADRILATERAL_SPHERICAL = 18,
	    NUM_OPTICAL_ENTITY_TYPES
    };

    struct LaunchParams
    {
        bool                        optical_errors;

        unsigned int                width;   // essentially number of rays launched and sun points 
        unsigned int                height;
        unsigned int                max_depth;
        unsigned long long          ray_offset; // Global offset for current branch

        // float4*                     hit_point_buffer;
        HitRecord*                  hit_buffer;
        float3*                     sun_dir_buffer;
        curandState*                rng_states;
        OptixTraversableHandle      handle;
        // int32_t*                    element_id_buffer;
        // uint8_t*                    hit_type_buffer;
        uint64_t*                   d_depth_exceeded_count; // Atomic counter: rays stopped by max depth, not absorption


        float3                      sun_vector;
        bool                        include_sun_shape_errors;
        SunShape                    sun_shape;      // OptixCSP::SunShape (mirrors SolTrace::Data::SunShape)
        float                       sigma;          // [mrad] for GAUSSIAN
        float                       half_width;     // [mrad] For PILLBOX
        float                       buie_kappa;     // Used by buie csr
        float                       buie_gamma;     // Used by buie csr
        float                       sun_max_angle;  // Calculated based on sunshape within SimulationData
        float                       sun_max_intensity;  // ^
		unsigned long long          sun_dir_seed;   // seed for the sun direction randomization
        GenType                     sun_gen_type;

        float*                      sun_user_angle; // User defined sun angle
        float*                      sun_user_intensity; // User defined sun intensity
        int                         sun_user_capacity; // Number of user defined values

        float3                      sun_v0;
        float3                      sun_v1;
        float3                      sun_v2;
        float3                      sun_v3;

	    GeometryDataST*             geometry_data_array;
		MaterialData*               material_data_array_front;
        MaterialData*               material_data_array_back;
    };

    struct PerRayData
    {
        unsigned int ray_path_index;  // Index of the ray in the ray path buffer
        unsigned int depth;           // Trace depth
    };

} // end namespace OptixCSP

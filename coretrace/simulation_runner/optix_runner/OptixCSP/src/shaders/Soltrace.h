#pragma once

#include "GeometryDataST.h"
#include "MaterialDataST.h"
#include "soltrace_constants.h"

#include <vector_types.h>
#include <optix.h>
#include <curand_kernel.h>

namespace OptixCSP{

    const unsigned int NUM_ATTRIBUTE_VALUES = 4u;
    const unsigned int NUM_PAYLOAD_VALUES   = 2u;
    const unsigned int MAX_TRACE_DEPTH      = 5u;

    struct HitGroupData
    {
        MaterialData material_data;
    };

    enum RayType
    {
        RAY_TYPE_RADIANCE = 0,
        RAY_TYPE_COUNT = 1         // not using occlusion/shadow rays atm
    };

    enum OpticalEntityType : unsigned int {
        RECTANGLE_FLAT       = 0,
        RECTANGLE_PARABOLIC  = 1,
        CYLINDRICAL          = 2,
        TRIANGLE_FLAT        = 3,
        QUADRILATERAL_FLAT   = 4,
	    NUM_OPTICAL_ENTITY_TYPES
    };

    struct LaunchParams
    {
        bool                        optical_errors;

        unsigned int                width;   // essentially number of rays launched and sun points 
        unsigned int                height;
        int                         max_depth;
        unsigned int                ray_offset; // Global offset for current branch

        float4*                     hit_point_buffer;
        float3*                     sun_dir_buffer;
        curandState*                rng_states;
        OptixTraversableHandle      handle;
        int32_t*                    element_id_buffer;
        uint8_t*                    hit_type_buffer;

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

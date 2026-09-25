#pragma once

#include <cstdint>

#include <cuda_runtime.h>

#include "soltrace_constants.h"

namespace OptixCSP
{
    struct MaterialData
    {
        float reflectivity;
        float transmissivity;
        float mu; // relative refractive index: incident / transmitted
        float slope_error;
        float specularity_error;
        OpticalDistribution optical_dist; // uint8_t
		bool  use_refraction;  // todo: for now, the ray goes through the object if true, otherwise it reflects

        // Angle-dependent reflectivity or transmissivity (whichever use_refraction
        // selects), looked up by cosine of the incidence angle. offset/count index
        // into the shared, de-duplicated pool held in
        // LaunchParams::angular_table_cos/value.
        bool     use_angular_table;
        uint32_t angular_table_offset;
        uint32_t angular_table_count;
    };
}

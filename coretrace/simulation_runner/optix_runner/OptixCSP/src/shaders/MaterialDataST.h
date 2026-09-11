#pragma once

#include <cuda_runtime.h>

#include "soltrace_constants.h"

namespace OptixCSP
{
    struct MaterialData
    {
        float reflectivity;
        float transmissivity;
        float refractive_index_incident;
        float refractive_index_transmitted;
        float slope_error;
        float specularity_error;
        OpticalDistribution optical_dist; // uint8_t
		bool  use_refraction;  // todo: for now, the ray goes through the object if true, otherwise it reflects
    };
}

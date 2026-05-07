#include <optix_device.h>
#include <curand_kernel.h>
#include <vector_types.h>

// todo: move curand initializatin to global function

//#include <cuda/helpers.h>
//#include <cuda/random.h>
#include "Soltrace.h"
#include "soltrace_constants.h"

#include <cstdio>

// Launch parameters for soltrace
extern "C" {
    __constant__ OptixCSP::LaunchParams params;
}

namespace OptixCSP {

    // Halton sequence generator, used for quasi-random sampling
    // Generates a Halton sequence value for a given index and base
    __device__ float halton(int index, int base) {
        float f = 1.0f, result = 0.0f;
        while (index > 0) {
            f = f / base;
            result = result + f * (index % base);
            index = index / base;
        }
        return result;
    }

    // Generate a sample point within a parallelogram defined by the AABB (Axis-Aligned Bounding Box)
    // Uses the Halton sequence for sampling
    __device__ float3 haltonSampleInParallelogram(unsigned int sample_index) {
        // Generate Halton sequence values
        float u = halton(sample_index, 2); // Base 2 for x
        float v = halton(sample_index, 3); // Base 3 for y

        // Compute the two edge vectors of the parallelogram
        float3 edge1 = params.sun_v1 - params.sun_v0; // First edge vector
        float3 edge2 = params.sun_v3 - params.sun_v0; // Second edge vector

        return params.sun_v0 + u * edge1 + v * edge2;
    }

    __device__ float3 randomSampleInParallelogram(unsigned int ray_number)
    {
        curandState rng_state = params.rng_states[ray_number];

        const float u = curand_uniform(&rng_state);
        const float v = curand_uniform(&rng_state);

        params.rng_states[ray_number] = rng_state;

        const float3 edge1 = params.sun_v1 - params.sun_v0;
        const float3 edge2 = params.sun_v3 - params.sun_v0;

        return params.sun_v0 + u * edge1 + v * edge2;
    }

    // Sample a random ray direction within a cone defined by a maximum angle
    __device__ float3 sampleRayDirectionInCone_Pillbox(float3 dir, float half_angle, unsigned int ray_number) {
        curandState rng_state = params.rng_states[ray_number];

        const float half_angle_mrad = half_angle;

        // Build an orthonormal basis
        float3 w = normalize(dir);
        float3 u = normalize(cross(fabs(w.x) > 0.99f ? make_float3(0, 1, 0) : make_float3(1, 0, 0), w));
        float3 v = cross(w, u);

        float thetax = 0.0f;
        float thetay = 0.0f;
        float theta2 = 0.0f;
        do
        {
            thetax = 2.0f * half_angle_mrad * curand_uniform(&rng_state) - half_angle_mrad;
            thetay = 2.0f * half_angle_mrad * curand_uniform(&rng_state) - half_angle_mrad;
            theta2 = thetax * thetax + thetay * thetay;
        } while (theta2 > (half_angle_mrad * half_angle_mrad));

        const float theta_rad = 0.001f * sqrtf(theta2); // mrad -> rad
        const float phi = 2.0f * M_PIf * curand_uniform(&rng_state);
        const float sin_t = sinf(theta_rad);
        const float cos_t = cosf(theta_rad);

        params.rng_states[ray_number] = rng_state;

        // Transform to world space
        return normalize(sin_t * (cosf(phi) * u + sinf(phi) * v) + cos_t * w);
    }

    __device__ float3 sampleRayDirectionInCone_Gaussian(float3 dir, float sigma, unsigned int ray_number) {
        curandState rng = params.rng_states[ray_number];

        const float sigma_rad = sigma * 0.001f;   // Convert to rad

        // Build an orthonormal basis
        float3 w = normalize(dir);
        float3 u = normalize(cross(fabs(w.x) > 0.99f ? make_float3(0, 1, 0) : make_float3(1, 0, 0), w));
        float3 v = cross(w, u);

        float gx = curand_normal(&rng);
        float gy = curand_normal(&rng);

        float thetax = sigma_rad * gx;
        float thetay = sigma_rad * gy;
        float theta2 = thetax * thetax + thetay * thetay;
        float z = sqrtf(1.0f - theta2);

        params.rng_states[ray_number] = rng;

        // Transform to world space
        return normalize(thetax * u + thetay * v + z * w);
    }

    __device__ float3 sampleRayDirectionInCone_BuieCSR(float3 dir, float buie_kappa, float buie_gamma, unsigned int ray_number)
    {
        curandState rng = params.rng_states[ray_number];

        const float max_angle_mrad = params.sun_max_angle;      // [mrad]

        // Orthonormal basis about dir
        float3 w = normalize(dir);
        float3 u = normalize(cross(fabsf(w.x) > 0.99f ? make_float3(0, 1, 0) : make_float3(1, 0, 0), w));
        float3 v = cross(w, u);

        float theta = 0.0f;
        float theta2 = 0.0f;
        float stest = 0.0f;
        float max_int = params.sun_max_intensity;

        // Rejection sampling in theta-space, matching CPU logic
        do
        {
            // Uniform sample in square [-max_angle, max_angle]^2 in mrad
            float thetax = (2.0f * curand_uniform(&rng) - 1.0f) * max_angle_mrad;
            float thetay = (2.0f * curand_uniform(&rng) - 1.0f) * max_angle_mrad;

            theta2 = thetax * thetax + thetay * thetay;
            theta = sqrtf(theta2);

            if (theta <= 4.65f)  // within solar disc (mrad, as in CPU code)
            {
                // stest = cos(0.326 * theta) / cos(0.308 * theta);
                float t = theta;
                stest = cosf(0.326f * t) / cosf(0.308f * t);
            }
            else // within circumsolar region
            {
                // stest = exp(kappa) * pow(|theta|, gamma);
                stest = expf(buie_kappa) * powf(fabsf(theta), buie_gamma);
            }

            // Reject if outside max_angle or fails intensity test
        } while ((curand_uniform(&rng) > (stest / max_int)) || (theta2 > (max_angle_mrad * max_angle_mrad)));

        // Convert theta (mrad) to radians
        float theta_rad = theta * 0.001f;

        // Random azimuth
        float phi = 2.0f * M_PIf * curand_uniform(&rng);

        // Local direction in cone coordinates
        float sin_t = sinf(theta_rad);
        float cos_t = cosf(theta_rad);
        float3 local_dir = make_float3(
            sin_t * cosf(phi),
            sin_t * sinf(phi),
            cos_t
        );

        params.rng_states[ray_number] = rng;

        // Transform to world space
        float3 world_dir = normalize(local_dir.x * u + local_dir.y * v + local_dir.z * w);
        return world_dir;
    }

    __device__ float3 sampleRayDirectionInCone_LimbDarkened(float3 dir, unsigned int ray_number)
    {
        curandState rng = params.rng_states[ray_number];

        const float max_angle_mrad = params.sun_max_angle; // [mrad]
        const float max_int = params.sun_max_intensity;

        // Orthonormal basis about dir
        float3 w = normalize(dir);
        float3 u = normalize(cross(fabsf(w.x) > 0.99f ? make_float3(0, 1, 0) : make_float3(1, 0, 0), w));
        float3 v = cross(w, u);

        float theta = 0.f;
        float theta2 = 0.f;
        float stest = 0.f;

        do
        {
            float thetax = (2.f * curand_uniform(&rng) - 1.f) * max_angle_mrad;
            float thetay = (2.f * curand_uniform(&rng) - 1.f) * max_angle_mrad;
            theta2 = thetax * thetax + thetay * thetay;
            theta = sqrtf(theta2);

            stest = 1.f - 0.5138f * powf(theta / max_angle_mrad, 4.f);
        } while ((curand_uniform(&rng) > (stest / max_int)) || (theta2 > (max_angle_mrad * max_angle_mrad)));

        // Convert theta (mrad) to radians
        float theta_rad = theta * 0.001f;

        // Random azimuth
        float phi = 2.f * M_PIf * curand_uniform(&rng);

        // Local direction in cone coordinates
        float sin_t = sinf(theta_rad);
        float cos_t = cosf(theta_rad);
        float3 local_dir = make_float3(
            sin_t * cosf(phi),
            sin_t * sinf(phi),
            cos_t
        );

        params.rng_states[ray_number] = rng;

        // Transform to world space
        float3 world_dir = normalize(local_dir.x * u + local_dir.y * v + local_dir.z * w);
        return world_dir;
    }

    __device__ float3 sampleRayDirectionInCone_UserDefined(float3 dir, int user_capacity, float* user_angle, 
        float* user_intensity, unsigned int ray_number)
    {
        curandState rng = params.rng_states[ray_number];

        if (user_capacity <= 0 || user_angle == nullptr || user_intensity == nullptr)
        {
            return normalize(dir);
        }

        const float max_angle_mrad = params.sun_max_angle;  // [mrad]
        const float max_int = params.sun_max_intensity;

        // Orthonormal basis about dir
        float3 w = normalize(dir);
        float3 u = normalize(cross(fabsf(w.x) > 0.99f ? make_float3(0, 1, 0) : make_float3(1, 0, 0), w));
        float3 v = cross(w, u);

        float theta = 0.f;
        float theta2 = 0.f;
        float stest = 0.f;

        do
        {
            float thetax = (2.f * curand_uniform(&rng) - 1.f) * max_angle_mrad;
            float thetay = (2.f * curand_uniform(&rng) - 1.f) * max_angle_mrad;
            theta2 = thetax * thetax + thetay * thetay;
            theta = sqrtf(theta2);

            int i = 0;
            while (i < user_capacity - 1 && user_angle[i] < theta)
                i++;

            if (i == 0)
                stest = user_intensity[0];
            else
            {
                const float denom = user_angle[i] - user_angle[i - 1];
                if (fabsf(denom) <= 1.0e-7f)
                    stest = user_intensity[i];
                else
                {
                    stest = user_intensity[i - 1] + (user_intensity[i] - user_intensity[i - 1]) * (theta - user_angle[i - 1])
                        / denom;
                }
            }

        } while ((curand_uniform(&rng) > (stest / max_int)) || (theta2 > (max_angle_mrad * max_angle_mrad)));

        // Convert theta (mrad) to radians
        float theta_rad = theta * 0.001f;

        // Random azimuth
        float phi = 2.f * M_PIf * curand_uniform(&rng);

        // Local direction in cone coordinates
        float sin_t = sinf(theta_rad);
        float cos_t = cosf(theta_rad);
        float3 local_dir = make_float3(
            sin_t * cosf(phi),
            sin_t * sinf(phi),
            cos_t
        );

        params.rng_states[ray_number] = rng;

        // Transform to world space
        float3 world_dir = normalize(local_dir.x * u + local_dir.y * v + local_dir.z * w);
        return world_dir;
    }

}
// == Ray Generation Program - Sun Source (Parallelogram Sampling)
extern "C" __global__ void __raygen__sun_source()
{
    // Lookup location in launch grid
    const uint3 launch_idx = optixGetLaunchIndex();         // Index of the current launch thread
    const uint3 launch_dims = optixGetLaunchDimensions();   // Dimensions of the launch grid
    const unsigned int ray_number = launch_idx.y * launch_dims.x + launch_idx.x;  // Unique ray ID
    const unsigned int ray_number_global = ray_number + params.ray_offset;  // Global unique ray ID

    float3 sun_sample_pos;
    switch (params.sun_gen_type)
    {
        case(OptixCSP::GenType::RANDOM):
            sun_sample_pos = OptixCSP::randomSampleInParallelogram(ray_number);
            break;
        case(OptixCSP::GenType::HALTON):
            sun_sample_pos = OptixCSP::haltonSampleInParallelogram(ray_number_global);
            break;
        default:          
            return;
    }

    // Sample emission angle here - capturing sun distribution
    const float3 ray_gen_pos = sun_sample_pos;

    float3 init_ray_dir = -normalize(params.sun_vector);

    // Apply sun shape errors
    float3 ray_dir;
    if (params.include_sun_shape_errors)
    {
        switch (params.sun_shape)
        {
            case(OptixCSP::SunShape::PILLBOX):
                ray_dir = OptixCSP::sampleRayDirectionInCone_Pillbox(init_ray_dir, params.half_width, ray_number);
                break;
            case(OptixCSP::SunShape::GAUSSIAN):
                ray_dir = OptixCSP::sampleRayDirectionInCone_Gaussian(init_ray_dir, params.sigma, ray_number);
                break;
            case(OptixCSP::SunShape::BUIE_CSR):
                ray_dir = OptixCSP::sampleRayDirectionInCone_BuieCSR(init_ray_dir, params.buie_kappa, params.buie_gamma, ray_number);
                break;
            case(OptixCSP::SunShape::LIMBDARKENED):
                ray_dir = OptixCSP::sampleRayDirectionInCone_LimbDarkened(init_ray_dir, ray_number);
                break;
            case(OptixCSP::SunShape::USER_DEFINED):
                ray_dir = OptixCSP::sampleRayDirectionInCone_UserDefined(init_ray_dir, params.sun_user_capacity, 
                    params.sun_user_angle, params.sun_user_intensity, ray_number);
                break;
            default:
                assert(false);
                // Just return since the sun shape is not supported
                return;
        }
    }
    else
    {
        ray_dir = init_ray_dir;
    }
    
    
    //float3 ray_dir = OptixCSP::sampleRayDirectionInCone_Gaussian(init_ray_dir, params.max_sun_angle, ray_number);

    // Create the PerRayData structure to track ray state (e.g., path index and recursion depth)
    OptixCSP::PerRayData prd;
    prd.ray_path_index = ray_number;
    prd.depth = 0;

    // TODO make this a launch parameter
    params.hit_point_buffer[params.max_depth * prd.ray_path_index] = make_float4(0.0f, ray_gen_pos);
    params.element_id_buffer[params.max_depth * prd.ray_path_index] = OptixCSP::kElementIdRayGen;
    params.hit_type_buffer[params.max_depth * prd.ray_path_index] = OptixCSP::HitType::HIT_CREATE;
    params.sun_dir_buffer[prd.ray_path_index] = ray_dir;
    

    // Cast and trace the ray through the scene
    optixTrace(
        params.handle,               // Acceleration structure handle
        ray_gen_pos,                 // Ray origin
        ray_dir,                     // Ray direction
        0.001f,                      // Minimum ray distance (near hit distance)
        1e16f,                       // Maximum ray distance (far hit distance)
        0.0f,                        // Time parameter (static for now)
        OptixVisibilityMask(1),      // Visibility mask (e.g., to restrict ray interactions)
        OPTIX_RAY_FLAG_NONE,         // Ray flags (no special flags)
        OptixCSP::RAY_TYPE_RADIANCE, // Ray type (radiance for sunlight)
        OptixCSP::RAY_TYPE_COUNT,    // Number of ray types
        OptixCSP::RAY_TYPE_RADIANCE, // SBT offset (ray type to launch)
        reinterpret_cast<unsigned int&>(prd.ray_path_index),
        reinterpret_cast<unsigned int&>(prd.depth)  
    );
}
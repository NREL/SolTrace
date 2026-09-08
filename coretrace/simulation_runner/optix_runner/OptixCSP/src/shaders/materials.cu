#include "MaterialDataST.h"
#include "Soltrace.h"
#include "curand_kernel.h"
#include "soltrace_constants.h"
#include <optix_device.h>
#include <stdio.h>
#include <vector_types.h>

// Launch parameters for soltrace
extern "C"
{
    __constant__ OptixCSP::LaunchParams params;
}

namespace OptixCSP
{
static __device__ __inline__ OptixCSP::PerRayData getPayload()
{
    OptixCSP::PerRayData prd;
    prd.ray_path_index = optixGetPayload_0();
    prd.depth          = optixGetPayload_1();
    return prd;
}

static __device__ __inline__ void setPayload(const OptixCSP::PerRayData& prd)
{
    optixSetPayload_0(prd.ray_path_index);
    optixSetPayload_1(prd.depth);
}

// // 32-bit avalanche mix (fast, good diffusion)
// static __device__ __inline__ float rng_uniform(uint32_t x)
// {
//     x ^= x >> 16;
//     x *= 0x85EBCA6Bu;
//     x ^= x >> 13;
//     x *= 0xC2B2AE35u;
//     x ^= x >> 16;
//     return float(x >> 8) * (1.0f / 16777216.0f); // Scale to [0, 1)
// }
} // namespace OptixCSP

// Assumes that v is a unit vector
extern "C" __device__ __inline__ float3 orthonormal_vector(float3 v)
{
    // TODO: Need to handle w = c * v case...
    float3 u;
    if (fabs(v.x) < 0.9f)
    {
        // Code does the following:
        // w = make_float3(1.0f, 0.0f, 0.0f);
        // u = cross(v, w);
        u = make_float3(0.0f, v.z, -v.y);
    }
    else
    {
        // Code does the following:
        // w = make_float3(0.0f, 1.0f, 0.0f);
        // u = cross(v, w);
        u = make_float3(-v.z, 0.0f, v.x);
    }
    return normalize(u);
}

// Add perturbation ortogonal to given vector. Perturbation is uniform over
// a disk of radius a centered at the vector n. Returned vector is a
// unit vector.
extern "C" __device__ float3 apply_uniform_errors(float                 a,
                                                  float3                v,
                                                  OptixCSP::PerRayData& prd)
{
    curandState local_rng = params.rng_states[prd.ray_path_index];
    float3      eta       = orthonormal_vector(v);
    float3      xi        = cross(v, eta);
    float       phi       = 2.0f * M_PIf * curand_uniform(&local_rng);
    float       r         = a * sqrtf(curand_uniform(&local_rng));
    eta                   = r * cosf(phi) * eta + r * sinf(phi) * xi;
    params.rng_states[prd.ray_path_index] = local_rng;
    return normalize(v + eta);
}

// Add perturbation orthogonal to given vector. Magnitude is Gaussian with
// standard deviation sigma and the angle is uniform over [-pi, pi).
// Returned vector is a unit vector.
extern "C" __device__ float3 apply_gaussian_errors(float                 sigma,
                                                   float3                n,
                                                   OptixCSP::PerRayData& prd)
{
    curandState  local_rng = params.rng_states[prd.ray_path_index];
    float3       eta       = orthonormal_vector(n);
    const float3 xi        = cross(n, eta);
    float        phi       = curand_normal(&local_rng);
    float        r         = curand_normal(&local_rng);
    r                      = sigma * sqrtf(r * r + phi * phi);
    phi                    = 2.0f * M_PIf * curand_uniform(&local_rng);
    eta                    = r * cosf(phi) * eta + r * sinf(phi) * xi;
    params.rng_states[prd.ray_path_index] = local_rng;
    return normalize(n + eta);
}

extern "C" __device__ float3 apply_lambertian_errors(float,
                                                     float3                n,
                                                     OptixCSP::PerRayData& prd)
{
    curandState  local_rng = params.rng_states[prd.ray_path_index];
    const float3 eta       = orthonormal_vector(n);
    const float3 xi        = cross(n, eta);
    float        theta     = asinf(sqrtf(curand_uniform(&local_rng)));
    float        phi       = 2.0 * M_PIf * curand_uniform(&local_rng);
    float3 d = sinf(theta) * cosf(phi) * eta + sinf(theta) * sinf(phi) * xi +
               cosf(theta) * n;
    params.rng_states[prd.ray_path_index] = local_rng;
    return normalize(d);
}

extern "C" __global__ void __closesthit__element()
{
    // Determine if we are using optical errors
    const bool optical_errors = params.optical_errors;

    // Fetch the normal vector from the hit attributes passed by OptiX
    float3 object_normal = make_float3(__uint_as_float(optixGetAttribute_0()),
                                       __uint_as_float(optixGetAttribute_1()),
                                       __uint_as_float(optixGetAttribute_2()));
    // // Transform the object-space normal to world space using OptiX built-in
    // function float3 world_normal =
    // normalize(optixTransformNormalFromObjectToWorldSpace(object_normal)); All
    // reports from the intersection functions are in world coordinates
    float3 world_normal = normalize(object_normal);

    // Compute the facing normal, which handles the direction of the normal
    // based on the incoming ray direction
    const float3 ray_dir  = optixGetWorldRayDirection();
    float3       ffnormal = faceforward(world_normal, -ray_dir, world_normal);

    // Get the incoming ray's origin, direction, and max t (intersection
    // distance)
    const float3 ray_orig = optixGetWorldRayOrigin();
    const float  ray_t    = optixGetRayTmax();

    // Determine hit direction (front or back)
    const float dot_nd         = dot(ray_dir, world_normal);
    const bool  hit_front_face = (dot_nd < 0.0f);

    // Compute the hit point of the ray using its origin and direction, scaled
    // by the intersection distance (ray_t)
    const float3 hit_point = ray_orig + ray_t * ray_dir;

    OptixCSP::PerRayData prd = OptixCSP::getPayload();
    const unsigned int   new_depth =
        prd.depth + 1; // Increment the ray depth for recursive tracing

    // we have two scenarios here
    // if we use refraction, then we look at transmissivity to determine if the
    // ray will refract or get obsorbed. otherwise, it will get reflected.
    float3 new_dir;
    bool absorbed = false; // determine whether the ray is absorbed or not, this
                           // is montecarlo based, should be applied to
                           // reflection and refraction
    uint8_t hit_type = OptixCSP::HitType::HIT_UNASSIGNED;

    // Get optical properties
    OptixCSP::MaterialData material =
        hit_front_face
            ? params.material_data_array_front[optixGetPrimitiveIndex()]
            : params.material_data_array_back[optixGetPrimitiveIndex()];
    const float transmissivity     = material.transmissivity;
    const bool  use_transmissivity = material.use_refraction;
    const float reflectivity       = material.reflectivity;
    const float normal_sigma =
        1e-3f * material.slope_error; // Convert mrad to rad
    const float spec_sigma =
        1e-3f * material.specularity_error; // Convert mrad to rad
    const uint8_t error_type = material.optical_dist;

    // Surface normal (macro-surface) errors
    if (optical_errors)
    {
        if (error_type == OptixCSP::OpticalDistribution::OPT_GAUSSIAN)
        {
            ffnormal = apply_gaussian_errors(normal_sigma, ffnormal, prd);
        }
        else if (error_type == OptixCSP::OpticalDistribution::OPT_PILLBOX)
        {
            ffnormal = apply_uniform_errors(normal_sigma, ffnormal, prd);
        }
        else
        {
            ; // Intentional no-op
        }
    }

    // now we figure out the random number to determine if the ray is absorbed
    // or refracted float xi = OptixCSP::rng_uniform(prd); // random number in
    // [0,1)
    const float xi = curand_uniform(&params.rng_states[prd.ray_path_index]);

    if (use_transmissivity)
    {
        if (xi > transmissivity)
        {
            absorbed = true;
            hit_type = OptixCSP::HitType::HIT_ABSORB;
            // printf("ray is absorbed! ray index is %d, depth %d\n",
            // prd.ray_path_index, prd.depth);
        } // ray is absorbed
        else
        {
            new_dir  = refract(ray_dir, ffnormal);
            hit_type = OptixCSP::HitType::HIT_TRANSMIT;
        }
    }
    else
    {
        if (xi > reflectivity)
        {
            absorbed = true;
            hit_type = OptixCSP::HitType::HIT_ABSORB;
            // printf("ray is absorbed! ray index is %d, depth %d\n",
            // prd.ray_path_index, prd.depth);
        } // ray is absorbed
        else
        {
            new_dir  = reflect(ray_dir, ffnormal);
            hit_type = OptixCSP::HitType::HIT_REFLECT;
        }
    }

    if (optical_errors && !absorbed)
    {
        // Optical (micro-surface) errors
        if (error_type == OptixCSP::OpticalDistribution::OPT_GAUSSIAN)
        {
            new_dir = apply_gaussian_errors(spec_sigma, new_dir, prd);
        }
        else if (error_type == OptixCSP::OpticalDistribution::OPT_PILLBOX)
        {
            new_dir = apply_uniform_errors(spec_sigma, new_dir, prd);
        }
        else if (error_type == OptixCSP::OpticalDistribution::OPT_DIFFUSE)
            new_dir = apply_lambertian_errors(spec_sigma, ffnormal, prd);
        else
        {
            ; // Intentional no-op
        }
    }

    ////////////////////////////////////////////////////////////////
    // // Alternative to the above doing this twice is to combine it
    // // into a single spot...may not work as intended for pillbox?
    // // Would be more efficient since we could combine the two
    // // error values during setup and then we would only need one
    // // field and we would avoid a whole if clause.
    // const float sigma = sqrt(4.0f * normal_sigma + spec_sigma);
    // if (optical_errors)
    // {
    //     // Optical (micro-surface) errors
    //     if (error_type == OptixCSP::OpticalDistribution::OPT_GAUSSIAN)
    //     {
    //         new_dir = apply_gaussian_errors(sigma, new_dir, prd);
    //     }
    //     else if (error_type == OptixCSP::OpticalDistribution::OPT_PILLBOX)
    //     {
    //         new_dir = apply_uniform_errors(sigma, new_dir, prd);
    //     }
    //     else
    //     {
    //         ; // Intentional no-op
    //     }
    // }
    /////////////////////////////////////////////////////////////////

    // Check if the maximum recursion depth has not been reached
    if (new_depth < params.max_depth)
    {

        // Get buffer slot
        const unsigned int slot =
            params.max_depth * prd.ray_path_index + new_depth;

        // Store the hit point in the hit point buffer (used for visualization
        // or further calculations)
        params.hit_buffer[slot].hit_point = make_float4(new_depth, hit_point);

        // Store element id
        const int32_t elementId =
            params.geometry_data_array[optixGetPrimitiveIndex()].id;
        params.hit_buffer[slot].element_id = elementId;

        // Store hit type
        params.hit_buffer[slot].hit_type = hit_type;

        // Store the reflected direction in its buffer (used for visualization
        // or further calculations)
        /*
        params.reflected_dir_buffer[params.max_depth * prd.ray_path_index +
        new_depth] = make_float4(1.0f, reflected_dir);
        */

        // Trace the reflected ray
        prd.depth = new_depth;
        if (!absorbed)
        {
            // printf("launching rays that are not absorbed: %d, depth, %d\n",
            // prd.ray_path_index, prd.depth - 1);
            optixTrace(
                params.handle, // The handle to the acceleration structure
                hit_point,     // The starting point of the reflected ray
                new_dir,       // The direction of the reflected ray
                0.01f, // A small offset to avoid self-intersection (shadow
                       // acne)
                1e16f, // Maximum distance the ray can travel
                0.0f,  // Ray time (used for time-dependent effects)
                OptixVisibilityMask(1), // Visibility mask (defines what the ray
                                        // can interact with)
                OPTIX_RAY_FLAG_NONE,    // Ray flags (no special flags for now)
                OptixCSP::RAY_TYPE_RADIANCE, // Use the radiance ray type
                OptixCSP::RAY_TYPE_COUNT,    // Total number of ray types
                OptixCSP::RAY_TYPE_RADIANCE, // The ray type's offset into the
                                             // SBT
                reinterpret_cast<unsigned int&>(
                    prd.ray_path_index), // Pass the ray path index
                reinterpret_cast<unsigned int&>(
                    prd.depth) // Pass the updated depth
            );
        }
        else
        {
            // printf("ray %d is absorbed, terminate! depth = %d\n",
            // prd.ray_path_index, prd.depth - 1);
            prd.depth = params.max_depth; // terminate the ray by setting depth
                                          // to max depth
        }
    }
    else if (!absorbed)
    {
        // Ray hit an element but max depth was reached; count it (absorption at
        // this depth does not count).
        atomicAdd(reinterpret_cast<unsigned long long*>(
                      params.d_depth_exceeded_count),
                  1ULL);
    }

    setPayload(prd);
}

extern "C" __global__ void __miss__ms()
{
    OptixCSP::PerRayData prd = OptixCSP::getPayload();

    // Do not record direct solar misses; no intersection exists to exit from.
    const unsigned int exit_depth = prd.depth + 1;
    if (prd.depth > 0 && exit_depth < params.max_depth)
    {
        const float3 exit_direction = normalize(optixGetWorldRayDirection());
        const float3 exit_point = optixGetWorldRayOrigin() + exit_direction;
        const unsigned int slot =
            params.max_depth * prd.ray_path_index + exit_depth;

        params.hit_buffer[slot].hit_point = make_float4(exit_depth, exit_point);
        params.hit_buffer[slot].element_id = OptixCSP::kElementIdUnassigned;
        params.hit_buffer[slot].hit_type = OptixCSP::HitType::HIT_EXIT;
    }

    // Set the payload values to 0, indicating that the ray missed all geometry.
    optixSetPayload_0(0); // Default value
    optixSetPayload_1(0);
}

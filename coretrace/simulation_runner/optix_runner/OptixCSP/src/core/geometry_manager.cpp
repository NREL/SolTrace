#include "geometry_manager.h"
#include "shaders/GeometryDataST.h"
#include "sun_utils.h"
#include "soltrace_state.h"
#include "utils/util_check.hpp"

#include <optix_stubs.h>

#include <iostream>
#include <stdexcept>
#include <sstream>
#include <vector>

using namespace OptixCSP;

// Collect geometry and material information from the scene elements
void GeometryManager::collect_geometry_info(const std::vector<std::shared_ptr<CspElement>> &element_list,
                                            LaunchParams &params)
{
    m_aabb_list_H.clear();           // Clear the existing AABB list
    m_sbt_index_H.clear();           // Clear the existing SBT index list
    m_geometry_data_array_H.clear(); // Clear the existing geometry data array
    m_material_data_array_front_H.clear();
    m_material_data_array_back_H.clear();

    m_obj_counts = static_cast<uint32_t>(element_list.size()); // Number of objects in the scene

    // Resize
    m_aabb_list_H.resize(m_obj_counts);
    m_geometry_data_array_H.resize(m_obj_counts);
    m_sbt_index_H.resize(m_obj_counts);
    m_material_data_array_front_H.resize(m_obj_counts);
    m_material_data_array_back_H.resize(m_obj_counts);

    for (uint32_t i = 0; i < m_obj_counts; i++)
    {

        std::shared_ptr<CspElement> element = element_list[i];

        // Create an OptixAabb from the geometry data
        OptixAabb aabb;
        uint32_t sbt_offset = 0;

        if (element->get_aperture_type() == ApertureType::ANNULUS)
        {
            if (element->get_surface_type() == SurfaceType::FLAT)
            {
                sbt_offset = static_cast<uint32_t>(OpticalEntityType::ANNULUS_FLAT);
            }
            else if (element->get_surface_type() == SurfaceType::PARABOLIC)
            {
                sbt_offset = static_cast<uint32_t>(OpticalEntityType::ANNULUS_PARABOLIC);
            }
            else
            {
                std::stringstream ss;
                ss << "Unimplemented surface type ("
                   << static_cast<int>(element->get_surface_type())
                   << ") for annular aperture ("
                   << static_cast<int>(element->get_aperture_type());
                throw std::runtime_error(ss.str());
            }
        }

        if (element->get_aperture_type() == ApertureType::CIRCLE)
        {
            if (element->get_surface_type() == SurfaceType::FLAT)
            {
                sbt_offset = static_cast<uint32_t>(OpticalEntityType::CIRCLE_FLAT);
            }
            else if (element->get_surface_type() == SurfaceType::PARABOLIC)
            {
                sbt_offset = static_cast<uint32_t>(OpticalEntityType::CIRCLE_PARABOLIC);
            }
            else
            {
                std::stringstream ss;
                ss << "Unimplemented surface type ("
                   << static_cast<int>(element->get_surface_type())
                   << ") for circular aperture ("
                   << static_cast<int>(element->get_aperture_type());
                throw std::runtime_error(ss.str());
            }
        }

        if (element->get_aperture_type() == ApertureType::HEXAGON)
        {
            if (element->get_surface_type() == SurfaceType::FLAT)
            {
                sbt_offset = static_cast<uint32_t>(OpticalEntityType::HEXAGON_FLAT);
            }
            else if (element->get_surface_type() == SurfaceType::PARABOLIC)
            {
                sbt_offset = static_cast<uint32_t>(OpticalEntityType::HEXAGON_PARABOLIC);
            }
            else
            {
                std::stringstream ss;
                ss << "Unimplemented surface type ("
                   << static_cast<int>(element->get_surface_type())
                   << ") for hexagon aperture ("
                   << static_cast<int>(element->get_aperture_type());
                throw std::runtime_error(ss.str());
            }
        }

        if (element->get_aperture_type() == ApertureType::RECTANGLE)
        {
            if (element->get_surface_type() == SurfaceType::PARABOLIC)
            {
                sbt_offset = static_cast<uint32_t>(OpticalEntityType::RECTANGLE_PARABOLIC);
            }
            else if (element->get_surface_type() == SurfaceType::FLAT)
            {
                sbt_offset = static_cast<uint32_t>(OpticalEntityType::RECTANGLE_FLAT);
            }
            else if (element->get_surface_type() == SurfaceType::CYLINDER)
            {
                sbt_offset = static_cast<uint32_t>(OpticalEntityType::CYLINDRICAL);
            }
            else
            {
                std::stringstream ss;
                ss << "Unimplemented surface type ("
                   << static_cast<int>(element->get_surface_type())
                   << ") for rectangle aperture ("
                   << static_cast<int>(element->get_aperture_type());
                throw std::runtime_error(ss.str());
            }
        }

        if (element->get_aperture_type() == ApertureType::TRIANGLE)
        {
            if (element->get_surface_type() == SurfaceType::FLAT)
            {
                sbt_offset = static_cast<uint32_t>(OpticalEntityType::TRIANGLE_FLAT);
            }
            else if (element->get_surface_type() == SurfaceType::PARABOLIC)
            {
                sbt_offset = static_cast<uint32_t>(OpticalEntityType::TRIANGLE_PARABOLIC);
            }
            else
            {
                std::stringstream ss;
                ss << "Unimplemented surface type ("
                   << static_cast<int>(element->get_surface_type())
                   << ") for triangle aperture ("
                   << static_cast<int>(element->get_aperture_type());
                throw std::runtime_error(ss.str());
            }
        }

        if (element->get_aperture_type() == ApertureType::QUADRILATERAL)
        {
            if (element->get_surface_type() == SurfaceType::FLAT)
            {
                sbt_offset = static_cast<uint32_t>(OpticalEntityType::QUADRILATERAL_FLAT);
            }
            else if (element->get_surface_type() == SurfaceType::PARABOLIC)
            {
                sbt_offset = static_cast<uint32_t>(OpticalEntityType::QUADRILATERAL_PARABOLIC);
            }
            else
            {
                std::stringstream ss;
                ss << "Unimplemented surface type ("
                   << static_cast<int>(element->get_surface_type())
                   << ") for quadrilateral aperture ("
                   << static_cast<int>(element->get_aperture_type());
                throw std::runtime_error(ss.str());
            }
        }

        float buffer = 0;

        aabb.minX = (float)(element->get_lower_bounding_box()[0]) - buffer;
        aabb.minY = (float)(element->get_lower_bounding_box()[1]) - buffer;
        aabb.minZ = (float)(element->get_lower_bounding_box()[2]) - buffer;

        aabb.maxX = (float)(element->get_upper_bounding_box()[0]) + buffer;
        aabb.maxY = (float)(element->get_upper_bounding_box()[1]) + buffer;
        aabb.maxZ = (float)(element->get_upper_bounding_box()[2]) + buffer;

        m_aabb_list_H[i] = aabb;       // Store the AABB in the list
        m_sbt_index_H[i] = sbt_offset; // Store the SBT index
        m_geometry_data_array_H[i] = element->toDeviceGeometryData();

        // now we set the material data for each element, use placeholder values for now
        m_material_data_array_front_H[i] = element->toDeviceMaterialDataFront();
        m_material_data_array_back_H[i] = element->toDeviceMaterialDataBack();
    }

    // print out computed minimum distance
    if(m_verbose)
        std::cout << "Minimum distance to sun plane: " << m_sun_plane_distance << std::endl;
}

void GeometryManager::compute_sun_plane_H(LaunchParams &params)
{

    m_sun_plane_distance = -1;
    float3 sun_vector = params.sun_vector;

    float *sun_plane_dist_D;
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&sun_plane_dist_D), sizeof(float)));
    CUDA_CHECK(cudaMemset(sun_plane_dist_D, 0, sizeof(float)));

    compute_d_on_gpu(reinterpret_cast<const OptixAabb *>(m_aabb_list_D), m_obj_counts, sun_vector, sun_plane_dist_D);

    CUDA_CHECK(cudaMemcpy(reinterpret_cast<void *>(&m_sun_plane_distance),
                          reinterpret_cast<void *>(sun_plane_dist_D),
                          sizeof(float),
                          cudaMemcpyDeviceToHost));

    // Place the sun plane just above the scene (1 world-unit margin) so that
    // rays always travel at least tmin before reaching any element, without
    // amplifying angular errors from sun shape perturbations.
    m_sun_plane_distance = m_sun_plane_distance + 1.0f;

    float3 sun_u, sun_v;
    float3 axis = (abs(sun_vector.x) < 0.9f) ? make_float3(1.0f, 0.0f, 0.0f) : make_float3(0.0f, 1.0f, 0.0f);
    sun_u = normalize(cross(axis, sun_vector));
    sun_v = normalize(cross(sun_vector, sun_u));
    // allocate uv bounds on device, float array of size four
    float sun_uv_bounds[4] = {0.0f, 0.0f, 0.0f, 0.0f};
    // allocate memory for the sun_u_bounds on device
    float *sun_uv_bounds_D;
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&sun_uv_bounds_D), 4 * sizeof(float)));

    compute_uv_bounds_on_gpu(reinterpret_cast<const OptixAabb *>(m_aabb_list_D),
                             m_obj_counts,
                             m_sun_plane_distance,
                             sun_vector,
                             sun_u,
                             sun_v,
                             tan(params.sun_max_angle * 1e-3), // Convert from mrad to rad
                             sun_uv_bounds_D);

    // Copy the computed bounds back to the host
    cudaMemcpy(sun_uv_bounds, sun_uv_bounds_D, 4 * sizeof(float), cudaMemcpyDeviceToHost);

    float u_min = sun_uv_bounds[0];
    float u_max = sun_uv_bounds[1];
    float v_min = sun_uv_bounds[2];
    float v_max = sun_uv_bounds[3];

    params.sun_v0 = u_min * sun_u + v_min * sun_v + m_sun_plane_distance * sun_vector;
    params.sun_v1 = u_max * sun_u + v_min * sun_v + m_sun_plane_distance * sun_vector;
    params.sun_v2 = u_max * sun_u + v_max * sun_v + m_sun_plane_distance * sun_vector;
    params.sun_v3 = u_min * sun_u + v_max * sun_v + m_sun_plane_distance * sun_vector;

    CUDA_CHECK(cudaFree(reinterpret_cast<void *>(sun_plane_dist_D)));
    CUDA_CHECK(cudaFree(reinterpret_cast<void *>(sun_uv_bounds_D)));
}

void GeometryManager::clean_up()
{
    if (m_aabb_list_D)
    {
        CUDA_CHECK(cudaFree(reinterpret_cast<void *>(m_aabb_list_D)));
        m_aabb_list_D = 0;
    }

    if (m_sbt_index_D)
    {
        CUDA_CHECK(cudaFree(reinterpret_cast<void *>(m_sbt_index_D)));
        m_sbt_index_D = 0;
    }

    if (m_temp_buffer)
    {
        CUDA_CHECK(cudaFree(reinterpret_cast<void *>(m_temp_buffer)));
        m_temp_buffer = 0;
    }

    if (m_output_buffer)
    {
        CUDA_CHECK(cudaFree(reinterpret_cast<void *>(m_output_buffer)));
        m_output_buffer = 0;
    }

    m_temp_buffer_size = 0;
    m_output_buffer_size = 0;
    m_obj_counts = 0;
    m_sun_plane_distance = -1.0f;
    m_aabb_input = {};
    m_accel_build_options = {};
}

GeometryManager::~GeometryManager() noexcept
{
    try
    {
        clean_up();
    }
    catch (const std::exception &error)
    {
        std::cerr << "[OptixRunner] Geometry cleanup failed during destruction: "
                  << error.what() << '\n';
    }
    catch (...)
    {
        std::cerr << "[OptixRunner] Geometry cleanup failed during destruction with an unknown error.\n";
    }
}

void GeometryManager::create_geometries(LaunchParams &params)
{

    // Allocate memory on the device for the AABB array.
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&m_aabb_list_D), m_obj_counts * sizeof(OptixAabb)));
    CUDA_CHECK(cudaMemcpy(reinterpret_cast<void *>(m_aabb_list_D),
                          m_aabb_list_H.data(),
                          m_obj_counts * sizeof(OptixAabb),
                          cudaMemcpyHostToDevice));

    compute_sun_plane_H(params);

    // populate aabb_input_flags vector, size of types, no rebuild
    std::vector<uint32_t> aabb_input_flags(NUM_OPTICAL_ENTITY_TYPES);
    for (int i = 0; i < NUM_OPTICAL_ENTITY_TYPES; i++)
    {
        aabb_input_flags[i] = OPTIX_GEOMETRY_FLAG_DISABLE_ANYHIT;
    }

    // device vector for SBT index, no need to rebuild
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&m_sbt_index_D), m_obj_counts * sizeof(uint32_t)));
    CUDA_CHECK(cudaMemcpy(reinterpret_cast<void *>(m_sbt_index_D),
                          m_sbt_index_H.data(),
                          m_obj_counts * sizeof(uint32_t),
                          cudaMemcpyHostToDevice));

    // Configure the input for the GAS build process.
    m_aabb_input.type = OPTIX_BUILD_INPUT_TYPE_CUSTOM_PRIMITIVES;
    m_aabb_input.customPrimitiveArray.aabbBuffers = &m_aabb_list_D;
    m_aabb_input.customPrimitiveArray.flags = aabb_input_flags.data();
    m_aabb_input.customPrimitiveArray.numSbtRecords = NUM_OPTICAL_ENTITY_TYPES;
    m_aabb_input.customPrimitiveArray.numPrimitives = m_obj_counts;
    m_aabb_input.customPrimitiveArray.sbtIndexOffsetBuffer = m_sbt_index_D;
    m_aabb_input.customPrimitiveArray.sbtIndexOffsetSizeInBytes = sizeof(uint32_t);
    m_aabb_input.customPrimitiveArray.primitiveIndexOffset = 0;

    // Set up acceleration structure (AS) build options.
    m_accel_build_options = {
        OPTIX_BUILD_FLAG_ALLOW_UPDATE, // allow update
        OPTIX_BUILD_OPERATION_BUILD    // operation type, build a new aceleration structure
    };

    OptixAccelBufferSizes gas_buffer_sizes; // sizes for temp and output buffers.

    // Query the memory usage required for building the GAS.
    OPTIX_CHECK(optixAccelComputeMemoryUsage(m_state.context,
                                             &m_accel_build_options,
                                             &m_aabb_input,
                                             1,
                                             &gas_buffer_sizes));

    m_temp_buffer_size = gas_buffer_sizes.tempSizeInBytes;
    m_output_buffer_size = gas_buffer_sizes.outputSizeInBytes;

    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&m_temp_buffer), m_temp_buffer_size));
    CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&m_output_buffer), m_output_buffer_size));

    // Build the GAS.
    OPTIX_CHECK(optixAccelBuild(m_state.context, // OptiX context
                                m_state.stream,  // CUDA stream (default is 0)
                                &m_accel_build_options,
                                &m_aabb_input,
                                1,
                                m_temp_buffer,
                                m_temp_buffer_size,
                                m_output_buffer,
                                m_output_buffer_size,
                                &m_state.gas_handle, // Output handle for the GAS
                                nullptr,             // Emitted properties (not used here)
                                0));                 // Number of emitted properties
}

void GeometryManager::update_geometry_info(const std::vector<std::shared_ptr<CspElement>> &element_list,
                                           LaunchParams &params)
{
    // Recollect geometry info
    collect_geometry_info(element_list, params);

    // update device aabb list
    CUDA_CHECK(cudaMemcpyAsync(
        reinterpret_cast<void *>(m_aabb_list_D),
        m_aabb_list_H.data(),
        m_aabb_list_H.size() * sizeof(OptixAabb),
        cudaMemcpyHostToDevice, m_state.stream));

    std::vector<uint32_t> aabb_input_flags(NUM_OPTICAL_ENTITY_TYPES);
    for (int i = 0; i < NUM_OPTICAL_ENTITY_TYPES; i++)
    {
        aabb_input_flags[i] = OPTIX_GEOMETRY_FLAG_DISABLE_ANYHIT;
    }

    m_aabb_input.customPrimitiveArray.flags = aabb_input_flags.data();

    CUDA_CHECK(cudaMemcpyAsync(
        reinterpret_cast<void *>(m_aabb_input.customPrimitiveArray.aabbBuffers[0]),
        m_aabb_list_H.data(),
        m_aabb_list_H.size() * sizeof(OptixAabb),
        cudaMemcpyHostToDevice, m_state.stream));

    m_accel_build_options.operation = OPTIX_BUILD_OPERATION_UPDATE; // set to update

    OPTIX_CHECK(optixAccelBuild(m_state.context, // OptiX context
                                m_state.stream,  // CUDA stream (default is 0)
                                &m_accel_build_options,
                                &m_aabb_input,
                                1,
                                m_temp_buffer,
                                m_temp_buffer_size,
                                m_output_buffer,
                                m_output_buffer_size,
                                &m_state.gas_handle, // Output handle for the GAS
                                nullptr,             // Emitted properties (not used here)
                                0));                 // Number of emitted properties

    compute_sun_plane_H(params);
}

#ifndef NOMINMAX
#define NOMINMAX
#endif

#include "soltrace_system.h"
#include "ray_utils.h"

#include "CspElement.h"
#include "data_manager.h"
#include "geometry_manager.h"
#include "pipeline_manager.h"
#include "soltrace_constants.h"
#include "soltrace_type.h"
#include "timer.h"

#include "shaders/Soltrace.h"

#include "utils/util_record.hpp"
#include "utils/util_check.hpp"
#include "utils/math_util.h"

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <limits>

#include <optix_function_table_definition.h>
#include <optix_stubs.h>

using namespace OptixCSP;

// TODO: optix related type should go into one header file
// i don't know what we should put here in hit group record, material is handled through a global array
// we can leave this empty for now ...
// note that this is has to be per optical entity type.
typedef Record<OptixCSP::HitGroupData> HitGroupRecord;

SolTraceSystem::SolTraceSystem()
    : m_number_of_rays(0),
      m_max_number_of_rays(0),
      m_batch_size(0),
      m_max_ray_depth(DEFAULT_MAX_TRACE_DEPTH),
      m_verbose(false),
      m_mem_free_before(0),
      m_mem_free_after(0),
      m_optical_errors(false),
      m_n_hit_rays(0),
      m_n_sun_rays(0),
      m_n_depth_exceeded_rays(0),
      m_include_sun_shape_errors(false),
      m_timer_setup(),
      m_timer_trace(),
      m_timer_aabb(),
      m_timer_geometry(),
      m_timer_pipeline(),
      m_timer_sbt(),
      m_timer_setup_buffer(),
      m_timer_optix_launch(),
      m_timer_collect_results(),
      m_n_run_iterations(0),
      m_mem_free_post_setup(0),
      geometry_manager(std::make_shared<GeometryManager>(m_state, m_verbose)),
      data_manager(std::make_shared<dataManager>()),
      pipeline_manager(std::make_shared<pipelineManager>(m_state)),
      m_sun(nullptr)
{
    unsigned int major = OPTIX_VERSION / 10000;
    unsigned int minor = (OPTIX_VERSION % 10000) / 100;
    unsigned int micro = OPTIX_VERSION % 100;
    if (m_verbose)
    {
        std::cout << "Using OPTIX Version: " << major
                  << "." << minor
                  << "." << micro
                  << std::endl;
    }

    CUDA_CHECK(cudaFree(0));
    CUcontext cuCtx = 0;
    OPTIX_CHECK(optixInit());
    OptixDeviceContextOptions options = {};
    options.logCallbackFunction = [](unsigned int level, const char *tag, const char *message, void *)
    {
        std::cerr << "[" << std::setw(2) << level << "][" << std::setw(12) << tag << "]: " << message << "\n";
    };
    options.logCallbackLevel = m_verbose ? 4 : 0;
    m_state.context = nullptr;
    m_state.stream = nullptr;
    m_state.sbt = {};
    m_state.d_gas_output_buffer = 0;
}

SolTraceSystem::~SolTraceSystem() noexcept
{
    try
    {
        clean_up();
    }
    catch (const std::exception &error)
    {
        std::cerr << "[OptixRunner] Cleanup failed during destruction: "
                  << error.what() << '\n';
    }
    catch (...)
    {
        std::cerr << "[OptixRunner] Cleanup failed during destruction with an unknown error.\n";
    }
}

void SolTraceSystem::set_max_ray_depth(uint64_t depth)
{
    if (depth < 2)
    {
        std::cerr << "[OptixRunner] WARNING: max_ray_depth (" << depth
                  << ") is below the minimum of 2. Clamping to 2.\n";
        depth = 2;
    }
    else if (depth > 255)
    {
        std::cerr << "[OptixRunner] WARNING: max_ray_depth (" << depth
                  << ") exceeds the maximum of 255. Clamping to 255.\n";
        depth = 255;
    }
    m_max_ray_depth = depth;
}

void SolTraceSystem::set_verbose(bool verbose)
{
    m_verbose = verbose;
    geometry_manager->set_verbose(verbose);
    pipeline_manager->set_verbose(verbose);
}

void SolTraceSystem::print_launch_params()
{
    if (!m_verbose)
    {
        return;
    }

    LaunchParams params = data_manager->launch_params_H;

    float3 sun_box_a = params.sun_v0 - params.sun_v1;
    float3 sun_box_b = params.sun_v1 - params.sun_v2;

    float sun_box_edge_a = sqrtf(sun_box_a.x * sun_box_a.x + sun_box_a.y * sun_box_a.y + sun_box_a.z * sun_box_a.z);
    float sun_box_edge_b = sqrtf(sun_box_b.x * sun_box_b.x + sun_box_b.y * sun_box_b.y + sun_box_b.z * sun_box_b.z);

    std::cout << "print launch params: " << std::endl;
    std::cout << "width              : " << params.width << std::endl;
    std::cout << "height             : " << params.height << std::endl;
    std::cout << "max_depth          : " << params.max_depth << std::endl;
    // std::cout << "hit_point_buffer   : " << params.hit_point_buffer << std::endl;
    std::cout << "hit_buffer         : " << params.hit_buffer << std::endl;
    std::cout << "sun_dir_buffer     : " << params.sun_dir_buffer << std::endl;
    std::cout << "sun_vector         : " << params.sun_vector.x << " " << params.sun_vector.y << " " << params.sun_vector.z << std::endl;
    // std::cout << "max_sun_angle      : " << params.max_sun_angle << std::endl;
    std::cout << "sun_v0             : " << params.sun_v0.x << " " << params.sun_v0.y << " " << params.sun_v0.z << std::endl;
    std::cout << "sun_v1             : " << params.sun_v1.x << " " << params.sun_v1.y << " " << params.sun_v1.z << std::endl;
    std::cout << "sun_v2             : " << params.sun_v2.x << " " << params.sun_v2.y << " " << params.sun_v2.z << std::endl;
    std::cout << "sun_v3             : " << params.sun_v3.x << " " << params.sun_v3.y << " " << params.sun_v3.z << std::endl;
    std::cout << "sun_box_edge_a     : " << sun_box_edge_a << std::endl;
    std::cout << "sun_box_edge_b     : " << sun_box_edge_b << std::endl;
}

void SolTraceSystem::initialize()
{

    // Create OptiX context on first initialize
    if (!m_state.context)
    {
        CUDA_CHECK(cudaFree(0));
        CUcontext cuCtx = 0;
        OPTIX_CHECK(optixInit());
        OptixDeviceContextOptions options = {};
        options.logCallbackFunction = [](unsigned int level, const char *tag, const char *message, void *)
        {
            std::cerr << "[" << std::setw(2) << level << "][" << std::setw(12) << tag << "]: " << message << "\n";
        };
        options.logCallbackLevel = m_verbose ? 4 : 0;
        OPTIX_CHECK(optixDeviceContextCreate(cuCtx, &options, &m_state.context));
    }

    // Create the CUDA stream immediately after the context so that all subsequent
    // GPU work (GAS build, kernel launches, optixLaunch) uses the same named stream.
    // Doing this before create_geometries() ensures optixAccelBuild and optixLaunch
    // share a single stream, giving explicit serial ordering without relying solely
    // on legacy null-stream synchronization semantics.
    if (!m_state.stream)
        CUDA_CHECK(cudaStreamCreate(&m_state.stream));

    {
        size_t mem_total;
        CUDA_CHECK(cudaMemGetInfo(&m_mem_free_before, &mem_total));
    }
    m_timer_setup.start();

    // set up input related to sun
    glm::dvec3 sun_vec = m_sun->get_position();
    glm::dvec3 sun_vec_norm = sun_vec;
    sun_vec_norm = glm::normalize(sun_vec_norm);

    data_manager->launch_params_H.sun_vector = make_float3(static_cast<float>(sun_vec_norm[0]),
                                                           static_cast<float>(sun_vec_norm[1]),
                                                           static_cast<float>(sun_vec_norm[2]));

    // Set generation type
    switch (m_sun->get_gen_type())
    {
    case (SolTrace::Data::GenType::RANDOM):
        data_manager->launch_params_H.sun_gen_type = OptixCSP::GenType::RANDOM;
        break;
    case (SolTrace::Data::GenType::HALTON):
        data_manager->launch_params_H.sun_gen_type = OptixCSP::GenType::HALTON;
        break;
    default:
        data_manager->launch_params_H.sun_gen_type = OptixCSP::GenType::UNKNOWN;
    }

    // Assign sun shape parameters (if necessary)
    data_manager->launch_params_H.include_sun_shape_errors = this->m_include_sun_shape_errors;
    data_manager->allocateSunUserData({}, {}); // Clear sun user data
    if (this->m_include_sun_shape_errors)
    {
        // Map SolTrace::Data::SunShape to OptixCSP::SunShape for device code
        switch (m_sun->get_shape())
        {
        case SolTrace::Data::SunShape::GAUSSIAN:
            data_manager->launch_params_H.sun_shape = SunShape::GAUSSIAN;
            data_manager->launch_params_H.sigma = static_cast<float>(m_sun->get_sigma());
            break;
        case SolTrace::Data::SunShape::PILLBOX:
            data_manager->launch_params_H.sun_shape = SunShape::PILLBOX;
            data_manager->launch_params_H.half_width = static_cast<float>(m_sun->get_half_width());
            break;
        case SolTrace::Data::SunShape::LIMBDARKENED:
            data_manager->launch_params_H.sun_shape = SunShape::LIMBDARKENED;
            break;
        case SolTrace::Data::SunShape::BUIE_CSR:
        {
            data_manager->launch_params_H.sun_shape = SunShape::BUIE_CSR;
            double kappa, gamma;
            m_sun->calculate_buie_parameters(kappa, gamma);
            data_manager->launch_params_H.buie_kappa = kappa;
            data_manager->launch_params_H.buie_gamma = gamma;
            break;
        }
        case SolTrace::Data::SunShape::USER_DEFINED:
        {
            data_manager->launch_params_H.sun_shape = SunShape::USER_DEFINED;
            std::vector<double> user_angle, user_intensity;
            m_sun->get_user_data(user_angle, user_intensity);

            if (user_angle.empty() || user_intensity.empty())
            {
                throw std::runtime_error("User-defined sun shape requires non-empty angle and intensity arrays.");
            }

            std::vector<float> user_angle_float(user_angle.begin(), user_angle.end());
            std::vector<float> user_intensity_float(user_intensity.begin(), user_intensity.end());
            data_manager->allocateSunUserData(user_angle_float, user_intensity_float);
            break;
        }
        case SolTrace::Data::SunShape::UNKNOWN:
        default:
            data_manager->launch_params_H.sun_shape = SunShape::UNKNOWN;
            break;
        }

        data_manager->launch_params_H.sun_max_angle = static_cast<float>(m_sun->get_max_sun_angle());
        data_manager->launch_params_H.sun_max_intensity = static_cast<float>(m_sun->get_max_intensity());
    }

    m_timer_aabb.reset();
    m_timer_aabb.start();
    geometry_manager->collect_geometry_info(m_element_list, data_manager->launch_params_H);
    m_timer_aabb.stop();

    m_timer_geometry.reset();
    m_timer_geometry.start();
    geometry_manager->create_geometries(data_manager->launch_params_H);
    m_timer_geometry.stop();

    // Pipeline setup.
    m_timer_pipeline.reset();
    m_timer_pipeline.start();
    pipeline_manager->set_max_trace_depth(m_max_ray_depth);
    pipeline_manager->createPipeline();
    m_timer_pipeline.stop();

    m_timer_sbt.reset();
    m_timer_sbt.start();
    create_shader_binding_table();
    m_timer_sbt.stop();

    // seed for randomization
    data_manager->launch_params_H.sun_dir_seed = m_seed;
    data_manager->launch_params_H.optical_errors = m_optical_errors;

    // Link the GAS handle.
    data_manager->launch_params_H.handle = m_state.gas_handle;
    data_manager->allocateGeometryDataArray(geometry_manager->get_geometry_data_array());
    data_manager->allocateMaterialDataArray(geometry_manager->get_material_data_array_front(),
                                            geometry_manager->get_material_data_array_back());

    if (m_verbose)
    {
        std::cout << "Time to compute AABB: " << m_timer_aabb.get_time_sec() << " seconds" << std::endl;
        std::cout << "Time to create geometries: " << m_timer_geometry.get_time_sec() << " seconds" << std::endl;
        std::cout << "Time to create pipeline: " << m_timer_pipeline.get_time_sec() << " seconds" << std::endl;
        std::cout << "Time to create SBT: " << m_timer_sbt.get_time_sec() << " seconds" << std::endl;

        print_launch_params();
    }

    data_manager->allocateLaunchParams();

    // Snapshot free GPU memory now that all setup allocations (BVH, pipeline,
    // SBT, geometry/material arrays, launch params) are complete but before any
    // ray buffers exist.  automatic_batch_size() uses this as a stable baseline
    // so that batch sizing is consistent across every run() call.
    // Memory used by setup = m_mem_free_before - m_mem_free_post_setup.
    {
        size_t mem_total;
        CUDA_CHECK(cudaMemGetInfo(&m_mem_free_post_setup, &mem_total));
    }

    m_timer_setup.stop();
}

void SolTraceSystem::run()
{
    // Initialize results
    m_hit_records.clear();
    m_hit_ray_ids.clear();
    m_n_hit_rays = 0;
    m_n_sun_rays = 0;
    m_n_depth_exceeded_rays = 0;
    uint_fast64_t N_ray_hit = 0;
    uint_fast64_t N_ray_gen = 0;

    m_timer_trace.reset();
    m_timer_trace.start();
    m_timer_setup_buffer.reset();
    m_timer_optix_launch.reset();
    m_timer_collect_results.reset();
    m_n_run_iterations = 0;
    m_compaction_timings = CompactionTimings{};

    // Allocate device buffers and initialize RNG states once (sizes are constant across the while loop).
    allocate_device_buffers();

    while (N_ray_hit < m_number_of_rays && N_ray_gen < m_max_number_of_rays)
    {
        ++m_n_run_iterations;

        // Update ray offset (pushed to device in setup_device_buffer)
        data_manager->launch_params_H.ray_offset = N_ray_gen;

        // Allocate buffer (sets data_manager->launch_params_H buffer)
        m_timer_setup_buffer.start();
        {
            setup_device_buffer();
        }
        m_timer_setup_buffer.stop();

        int width = data_manager->launch_params_H.width;
        int height = data_manager->launch_params_H.height;

        size_t mem_total;
        cudaMemGetInfo(&m_mem_free_after, &mem_total);

        if (m_verbose)
            std::cout << "Memory used by launch: " << (m_mem_free_before - m_mem_free_after) / (1024.0 * 1024.0) << " MB\n";

        // Launch the simulation.
        m_timer_optix_launch.start();
        {
            OPTIX_CHECK(optixLaunch(
                m_state.pipeline,
                m_state.stream, // Assume this stream is properly created.
                reinterpret_cast<CUdeviceptr>(data_manager->getDeviceLaunchParams()),
                sizeof(OptixCSP::LaunchParams),
                &m_state.sbt, // Shader Binding Table.
                width,        // Launch dimensions
                height,
                1));
            CUDA_SYNC_CHECK();
        }
        m_timer_optix_launch.stop();

        // Collect results
        m_timer_collect_results.start();
        {
            get_buffer_results();
        }
        m_timer_collect_results.stop();

        // Read back depth-exceeded count for this batch
        uint64_t iter_depth_exceeded = 0;
        CUDA_CHECK(cudaMemcpy(&iter_depth_exceeded,
                              data_manager->launch_params_H.d_depth_exceeded_count,
                              sizeof(uint64_t), cudaMemcpyDeviceToHost));
        m_n_depth_exceeded_rays += iter_depth_exceeded;

        N_ray_hit = m_n_hit_rays;
        N_ray_gen += width;
        m_n_sun_rays = N_ray_gen;
    }

    // Trim excess rays: remove ray groups from the tail until m_n_hit_rays == m_number_of_rays.
    // Each group starts at the last HIT_CREATE record in m_hit_records.
    while (m_trim_excess_rays && m_n_hit_rays > m_number_of_rays && !m_hit_records.empty())
    {
        // Walk backwards to find the last CREATE record
        auto rit = std::find_if(m_hit_records.rbegin(), m_hit_records.rend(),
                                [](const HitRecord &r)
                                { return r.hit_type == HitType::HIT_CREATE; });
        if (rit == m_hit_records.rend())
            break;
        m_hit_records.erase(std::prev(rit.base()), m_hit_records.end());
        m_hit_ray_ids.pop_back();
        --m_n_hit_rays;
    }
    // m_n_sun_rays = rays generated up to and including the last retained hit ray.
    if (!m_hit_ray_ids.empty())
        m_n_sun_rays = m_hit_ray_ids.back() + 1;

    m_timer_trace.stop();

    if (m_n_depth_exceeded_rays > 0)
        std::cout << "[SolTraceSystem] " << m_n_depth_exceeded_rays
                  << " ray(s) were terminated due to reaching max_depth ("
                  << static_cast<unsigned>(data_manager->launch_params_H.max_depth) << ").\n";

    if (m_verbose)
    {
        const double t_setup = m_timer_setup_buffer.get_time_sec();
        const double t_launch = m_timer_optix_launch.get_time_sec();
        const double t_collect = m_timer_collect_results.get_time_sec();
        const double t_total = t_setup + t_launch + t_collect;
        const double inv_n = m_n_run_iterations > 0 ? 1.0 / static_cast<double>(m_n_run_iterations) : 0.0;

        std::cout << "\n--- SolTraceSystem::run() timing (" << m_n_run_iterations << " iteration"
                  << (m_n_run_iterations == 1 ? "" : "s") << ") ---\n";
        std::cout << std::fixed << std::setprecision(6);
        std::cout << "  setup_device_buffer : total = " << t_setup << " s"
                  << "  avg = " << t_setup * inv_n << " s"
                  << "  fraction = " << (t_total > 0.0 ? 100.0 * t_setup / t_total : 0.0) << " %\n";
        std::cout << "  optixLaunch         : total = " << t_launch << " s"
                  << "  avg = " << t_launch * inv_n << " s"
                  << "  fraction = " << (t_total > 0.0 ? 100.0 * t_launch / t_total : 0.0) << " %\n";
        std::cout << "  get_buffer_results  : total = " << t_collect << " s"
                  << "  avg = " << t_collect * inv_n << " s"
                  << "  fraction = " << (t_total > 0.0 ? 100.0 * t_collect / t_total : 0.0) << " %\n";
        std::cout << "  total (3 sections)  : " << t_total << " s\n";
        std::cout << "----------------------------------------------\n";
    }
}

void SolTraceSystem::update()
{

    const size_t N_slots = static_cast<size_t>(data_manager->launch_params_H.width) * static_cast<size_t>(data_manager->launch_params_H.height) * static_cast<size_t>(data_manager->launch_params_H.max_depth);
    const size_t hit_buffer_size = N_slots * sizeof(HitRecord);

    // update aabb and sun plane accordingly
    geometry_manager->update_geometry_info(m_element_list, data_manager->launch_params_H);

    // update data on the device
    data_manager->updateGeometryDataArray(geometry_manager->get_geometry_data_array());
    CUDA_CHECK(cudaMemset(data_manager->launch_params_H.hit_buffer, 0, hit_buffer_size));

    data_manager->updateLaunchParams();
}

void SolTraceSystem::get_hp_output(std::vector<float4> &hp_vec,
                                   std::vector<uint_fast64_t> &raynumber_vec,
                                   std::vector<int32_t> &element_id_vec,
                                   std::vector<uint8_t> &hit_type_vec)
{
    hp_vec.clear();
    raynumber_vec.clear();
    element_id_vec.clear();
    hit_type_vec.clear();
    hp_vec.reserve(m_hit_records.size());
    raynumber_vec.reserve(m_hit_records.size());
    element_id_vec.reserve(m_hit_records.size());
    hit_type_vec.reserve(m_hit_records.size());

    uint_fast64_t ray_number = 0;
    for (const HitRecord &r : m_hit_records)
    {
        if (r.hit_type == HitType::HIT_CREATE)
            ++ray_number;
        hp_vec.push_back(r.hit_point);
        raynumber_vec.push_back(ray_number);
        element_id_vec.push_back(r.element_id);
        hit_type_vec.push_back(r.hit_type);
    }
}

void SolTraceSystem::clean_up()
{

    // Nothing to clean if device not initialized
    if (!m_state.context)
    {
        return;
    }

    CUDA_CHECK(cudaDeviceSynchronize());

    geometry_manager->clean_up();

    // destroy pipeline related resources
    pipeline_manager->cleanup();

    // destory CUDA stream
    if (m_state.stream)
    {
        CUDA_CHECK(cudaStreamDestroy(m_state.stream));
    }

    OPTIX_CHECK(optixDeviceContextDestroy(m_state.context));

    // Free OptiX shader binding table (SBT) memory
    CUDA_CHECK(cudaFree(reinterpret_cast<void *>(m_state.sbt.raygenRecord)));
    CUDA_CHECK(cudaFree(reinterpret_cast<void *>(m_state.sbt.missRecordBase)));
    CUDA_CHECK(cudaFree(reinterpret_cast<void *>(m_state.sbt.hitgroupRecordBase)));

    // Free OptiX GAS output buffer
    CUDA_CHECK(cudaFree(reinterpret_cast<void *>(m_state.d_gas_output_buffer)));

    // Free device-side launch parameter memory
    CUDA_CHECK(cudaFree(reinterpret_cast<void *>(data_manager->launch_params_H.hit_buffer)));
    CUDA_CHECK(cudaFree(reinterpret_cast<void *>(data_manager->launch_params_H.sun_dir_buffer)));

    data_manager->launch_params_H.hit_buffer = nullptr;
    data_manager->launch_params_H.sun_dir_buffer = nullptr;
    m_hit_buffer_size_allocated = 0;
    m_sun_dir_buffer_size_allocated = 0;

    free_compaction_scratch(m_compaction_scratch);

    data_manager->cleanup();

    m_state.context = nullptr;
    m_state.stream = nullptr;
    m_state.pipeline = nullptr;
    m_state.raygen_prog_group = nullptr;
    m_state.radiance_miss_prog_group = nullptr;
    m_state.geometry_module = nullptr;
    m_state.shading_module = nullptr;
    m_state.sun_module = nullptr;
    m_state.gas_handle = 0;
    m_state.sbt = {};
    m_state.d_gas_output_buffer = 0;

    m_mem_free_before = 0;
    m_mem_free_post_setup = 0;
    m_mem_free_after = 0;
}

void SolTraceSystem::reset()
{
    clean_up();

    m_element_list.clear();
    m_hit_records.clear();
    m_hit_ray_ids.clear();
    m_n_hit_rays = 0;
    m_n_sun_rays = 0;
    m_n_depth_exceeded_rays = 0;

    m_sun = nullptr;
    m_number_of_rays = 0;
    m_max_number_of_rays = 0;
}

// Create and configure the Shader Binding Table (SBT).
// The SBT is a crucial data structure in OptiX that links geometry and ray types
// with their corresponding programs (ray generation, miss, and hit group).
void SolTraceSystem::create_shader_binding_table()
{
    // Free any previously allocated SBT records to avoid leaks on re-initialization.
    CUDA_CHECK(cudaFree(reinterpret_cast<void *>(m_state.sbt.raygenRecord)));
    CUDA_CHECK(cudaFree(reinterpret_cast<void *>(m_state.sbt.missRecordBase)));
    CUDA_CHECK(cudaFree(reinterpret_cast<void *>(m_state.sbt.hitgroupRecordBase)));
    m_state.sbt = {};

    // Ray generation program record
    {
        CUdeviceptr d_raygen_record; // Device pointer to hold the raygen SBT record.
        size_t sizeof_raygen_record = sizeof(EmptyRecord);

        CUDA_CHECK(cudaMalloc(
            reinterpret_cast<void **>(&d_raygen_record),
            sizeof_raygen_record));

        EmptyRecord rg_sbt; // host

        optixSbtRecordPackHeader(m_state.raygen_prog_group, &rg_sbt);

        // Copy the raygen SBT record from host to device.
        CUDA_CHECK(cudaMemcpy(
            reinterpret_cast<void *>(d_raygen_record),
            &rg_sbt,
            sizeof_raygen_record,
            cudaMemcpyHostToDevice));

        // Assign the device pointer to the raygenRecord field in the SBT.
        m_state.sbt.raygenRecord = d_raygen_record;
    }

    // Miss program record
    {
        CUdeviceptr d_miss_record;
        size_t sizeof_miss_record = sizeof(EmptyRecord);

        CUDA_CHECK(cudaMalloc(
            reinterpret_cast<void **>(&d_miss_record),
            sizeof_miss_record * OptixCSP::RAY_TYPE_COUNT));

        EmptyRecord ms_sbt[OptixCSP::RAY_TYPE_COUNT];
        // Pack the program header into the first miss SBT record.
        optixSbtRecordPackHeader(m_state.radiance_miss_prog_group, &ms_sbt[0]);

        CUDA_CHECK(cudaMemcpy(
            reinterpret_cast<void *>(d_miss_record),
            ms_sbt,
            sizeof_miss_record * OptixCSP::RAY_TYPE_COUNT,
            cudaMemcpyHostToDevice));

        // Configure the SBT miss program fields.
        m_state.sbt.missRecordBase = d_miss_record;                                      // Base address of the miss records.
        m_state.sbt.missRecordCount = OptixCSP::RAY_TYPE_COUNT;                          // Number of miss records.
        m_state.sbt.missRecordStrideInBytes = static_cast<uint32_t>(sizeof_miss_record); // Stride between miss records.
    }

    // Hitgroup program record
    {
        // Total number of hitgroup records is the number of optical entity types
        const unsigned int count_records = OptixCSP::NUM_OPTICAL_ENTITY_TYPES;
        std::vector<HitGroupRecord> hitgroup_records_list(count_records);

        // now we need to populate hitgroup_records_list, basically match the
        // OpticalEntityType with the corresponding m_program_group
        for (unsigned int i = 0; i < count_records; i++)
        {

            OptixCSP::OpticalEntityType my_type = static_cast<OptixCSP::OpticalEntityType>(i);
            // initialize program handle and data
            OptixProgramGroup program_group_handle = pipeline_manager->getElementProgram(my_type);
            hitgroup_records_list[i].data.material_data = {0.875425, 0, 0, 0};
            OPTIX_CHECK(optixSbtRecordPackHeader(program_group_handle, &hitgroup_records_list[i].header));
        }

        // Allocate memory for hitgroup records on the device.
        CUdeviceptr hitgroup_records_D;
        size_t sizeof_hitgroup_record = sizeof(HitGroupRecord);
        CUDA_CHECK(cudaMalloc(
            reinterpret_cast<void **>(&hitgroup_records_D),
            sizeof_hitgroup_record * count_records));

        // Copy hitgroup records from host to device.
        CUDA_CHECK(cudaMemcpy(
            reinterpret_cast<void *>(hitgroup_records_D),
            hitgroup_records_list.data(),
            sizeof_hitgroup_record * count_records,
            cudaMemcpyHostToDevice));

        // Configure the SBT hitgroup fields.
        m_state.sbt.hitgroupRecordBase = hitgroup_records_D;                                     // Base address of hitgroup records.
        m_state.sbt.hitgroupRecordCount = count_records;                                         // Total number of hitgroup records.
        m_state.sbt.hitgroupRecordStrideInBytes = static_cast<uint32_t>(sizeof_hitgroup_record); // Stride size.
    }
}

void SolTraceSystem::allocate_device_buffers()
{
    // Set constant launch params (unchanged across the while loop).
    const uint_fast64_t effective_batch = determine_batch_size();
    data_manager->launch_params_H.width = static_cast<int>(effective_batch);
    data_manager->launch_params_H.height = 1;
    data_manager->launch_params_H.max_depth = m_max_ray_depth;

    const size_t hit_buffer_size = static_cast<size_t>(data_manager->launch_params_H.width) * static_cast<size_t>(data_manager->launch_params_H.height) * static_cast<size_t>(data_manager->launch_params_H.max_depth) * sizeof(HitRecord);
    const size_t sun_dir_size = static_cast<size_t>(data_manager->launch_params_H.width) * static_cast<size_t>(data_manager->launch_params_H.height) * sizeof(float3);

    // NOTE: cudaFree is nullptr safe

    if (data_manager->launch_params_H.hit_buffer == nullptr || m_hit_buffer_size_allocated != hit_buffer_size)
    {
        CUDA_CHECK(cudaFree(reinterpret_cast<void *>(data_manager->launch_params_H.hit_buffer)));
        CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&data_manager->launch_params_H.hit_buffer), hit_buffer_size));
        m_hit_buffer_size_allocated = hit_buffer_size;

        // Reallocate compaction scratch whenever ray-buffer dimensions change
        const uint64_t num_rays = data_manager->launch_params_H.width * data_manager->launch_params_H.height;
        const uint64_t max_depth = static_cast<uint64_t>(data_manager->launch_params_H.max_depth);
        allocate_compaction_scratch(m_compaction_scratch, num_rays, max_depth);
    }

    if (data_manager->launch_params_H.sun_dir_buffer == nullptr || m_sun_dir_buffer_size_allocated != sun_dir_size)
    {
        CUDA_CHECK(cudaFree(reinterpret_cast<void *>(data_manager->launch_params_H.sun_dir_buffer)));
        CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&data_manager->launch_params_H.sun_dir_buffer), sun_dir_size));
        m_sun_dir_buffer_size_allocated = sun_dir_size;
    }

    // Initialize RNG states once (sizes are constant across the while loop).
    // curand states are persistent on the device and advance naturally across kernel launches.
    const unsigned int num_rng_states = static_cast<unsigned int>(
        data_manager->launch_params_H.width * data_manager->launch_params_H.height);
    data_manager->ensureCurandStates(
        num_rng_states,
        data_manager->launch_params_H.sun_dir_seed,
        0,
        m_state.stream);

    data_manager->ensureDepthExceededCounter();
}

void SolTraceSystem::setup_device_buffer()
{
    const size_t hit_buffer_size = static_cast<size_t>(data_manager->launch_params_H.width) * static_cast<size_t>(data_manager->launch_params_H.height) * static_cast<size_t>(data_manager->launch_params_H.max_depth) * sizeof(HitRecord);
    const size_t sun_dir_size = static_cast<size_t>(data_manager->launch_params_H.width) * static_cast<size_t>(data_manager->launch_params_H.height) * sizeof(float3);

    CUDA_CHECK(cudaMemset(data_manager->launch_params_H.hit_buffer, 0, hit_buffer_size));
    CUDA_CHECK(cudaMemset(data_manager->launch_params_H.sun_dir_buffer, 0, sun_dir_size));
    CUDA_CHECK(cudaMemset(data_manager->launch_params_H.d_depth_exceeded_count, 0, sizeof(uint64_t)));

    data_manager->updateLaunchParams();
}

// Compacts the device hit buffer on the GPU, then copies only qualifying records to host.
// Rays that produced only a HIT_CREATE event (missed all elements) are discarded.
// Empty depth slots are discarded. The compacted HitRecord array is appended to
// m_hit_records and m_n_hit_rays is incremented by the number of newly collected hit rays.
void SolTraceSystem::get_buffer_results()
{
    const uint32_t num_rays = static_cast<uint32_t>(data_manager->launch_params_H.width *
                                                    data_manager->launch_params_H.height);
    const uint32_t max_depth = static_cast<uint32_t>(data_manager->launch_params_H.max_depth);

    const uint32_t n_new_hits = gpu_compact_hit_buffer(
        data_manager->launch_params_H.hit_buffer,
        num_rays,
        max_depth,
        data_manager->launch_params_H.ray_offset,
        m_hit_records,
        m_hit_ray_ids,
        m_state.stream,
        m_compaction_scratch,
        &m_compaction_timings);
    m_n_hit_rays += n_new_hits;
}

void SolTraceSystem::add_element(std::shared_ptr<CspElement> e)
{
    m_element_list.push_back(e);
}

double SolTraceSystem::get_time_trace()
{
    return m_timer_trace.get_time_sec();
}

double SolTraceSystem::get_time_setup()
{
    return m_timer_setup.get_time_sec();
}

void SolTraceSystem::print_timing() const
{
    const double t_setup = m_timer_setup.get_time_sec();
    const double t_aabb = m_timer_aabb.get_time_sec();
    const double t_geometry = m_timer_geometry.get_time_sec();
    const double t_pipeline = m_timer_pipeline.get_time_sec();
    const double t_sbt = m_timer_sbt.get_time_sec();

    const double t_trace = m_timer_trace.get_time_sec();
    const double t_buf_setup = m_timer_setup_buffer.get_time_sec();
    const double t_launch = m_timer_optix_launch.get_time_sec();
    const double t_collect = m_timer_collect_results.get_time_sec();

    const double inv_n = m_n_run_iterations > 0
                             ? 1.0 / static_cast<double>(m_n_run_iterations)
                             : 0.0;

    const auto pct = [](double num, double denom) -> double
    {
        return denom > 0.0 ? 100.0 * num / denom : 0.0;
    };

    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\n=== SolTraceSystem Timing Summary ===\n";

    std::cout << "\n--- initialize() ---\n";
    std::cout << "  AABB computation    : " << t_aabb << " s  (" << pct(t_aabb, t_setup) << " %)\n";
    std::cout << "  Geometry creation   : " << t_geometry << " s  (" << pct(t_geometry, t_setup) << " %)\n";
    std::cout << "  Pipeline creation   : " << t_pipeline << " s  (" << pct(t_pipeline, t_setup) << " %)\n";
    std::cout << "  SBT creation        : " << t_sbt << " s  (" << pct(t_sbt, t_setup) << " %)\n";
    std::cout << "  Total setup         : " << t_setup << " s\n";

    std::cout << "\n--- run() [" << m_n_run_iterations
              << " iteration" << (m_n_run_iterations == 1 ? "" : "s") << "] ---\n";
    std::cout << "  Setup device buffer : total = " << t_buf_setup << " s"
              << "  avg/iter = " << t_buf_setup * inv_n << " s"
              << "  (" << pct(t_buf_setup, t_trace) << " %)\n";
    std::cout << "  OptiX launch        : total = " << t_launch << " s"
              << "  avg/iter = " << t_launch * inv_n << " s"
              << "  (" << pct(t_launch, t_trace) << " %)\n";
    std::cout << "  Collect results     : total = " << t_collect << " s"
              << "  avg/iter = " << t_collect * inv_n << " s"
              << "  (" << pct(t_collect, t_trace) << " %)\n";
    if (m_compaction_timings.n_calls > 0)
    {
        const float inv_c = 1.0f / static_cast<float>(m_compaction_timings.n_calls);
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "    GPU pass 1 (count/scan/reduce) : total = " << m_compaction_timings.gpu_phase1_ms << " ms"
                  << "  avg/call = " << m_compaction_timings.gpu_phase1_ms * inv_c << " ms\n";
        std::cout << "    D->H scalars (3x memcpy)       : total = " << m_compaction_timings.scalar_dth_ms << " ms"
                  << "  avg/call = " << m_compaction_timings.scalar_dth_ms * inv_c << " ms\n";
        std::cout << "    GPU pass 2 (compact/select)    : total = " << m_compaction_timings.gpu_phase2_ms << " ms"
                  << "  avg/call = " << m_compaction_timings.gpu_phase2_ms * inv_c << " ms\n";
        std::cout << "    D->H bulk (records+ids)        : total = " << m_compaction_timings.bulk_dth_ms << " ms"
                  << "  avg/call = " << m_compaction_timings.bulk_dth_ms * inv_c << " ms\n";
        std::cout << std::fixed << std::setprecision(6);
    }
    std::cout << "  Total trace         : " << t_trace << " s\n";

    std::cout << "\n--- Grand Total ---\n";
    std::cout << "  Setup + Trace       : " << (t_setup + t_trace) << " s\n";

    std::cout << "\n--- GPU Memory Usage ---\n";
    constexpr double kMB = 1.0 / (1024.0 * 1024.0);
    if (m_mem_free_before > 0)
    {
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "  Free before setup   : " << m_mem_free_before * kMB << " MB\n";
        if (m_mem_free_post_setup > 0)
        {
            std::cout << "  Free after setup    : " << m_mem_free_post_setup * kMB << " MB\n";
            std::cout << "  Setup structures    : " << (m_mem_free_before - m_mem_free_post_setup) * kMB << " MB\n";
            if (m_mem_free_after > 0)
            {
                std::cout << "  Ray buffers         : " << (m_mem_free_post_setup - m_mem_free_after) * kMB << " MB\n";
                std::cout << "  Total used          : " << (m_mem_free_before - m_mem_free_after) * kMB << " MB\n";
            }
        }
        std::cout << std::fixed << std::setprecision(6);
    }
    std::cout << "=====================================\n";
}

double SolTraceSystem::get_sun_plane_area() const
{
    const LaunchParams &lp = data_manager->launch_params_H;
    // Parallelogram area spanned by two adjacent edges of the sun box
    const float3 a = make_float3(lp.sun_v0.x - lp.sun_v1.x, lp.sun_v0.y - lp.sun_v1.y, lp.sun_v0.z - lp.sun_v1.z);
    const float3 b = make_float3(lp.sun_v1.x - lp.sun_v2.x, lp.sun_v1.y - lp.sun_v2.y, lp.sun_v1.z - lp.sun_v2.z);
    const float3 cross = make_float3(
        a.y * b.z - a.z * b.y,
        a.z * b.x - a.x * b.z,
        a.x * b.y - a.y * b.x);
    return static_cast<double>(sqrtf(cross.x * cross.x + cross.y * cross.y + cross.z * cross.z));
}

uint_fast64_t SolTraceSystem::automatic_batch_size() const
{
    // Use the free-memory snapshot taken at the end of initialize(), after all
    // setup allocations (BVH, pipeline, SBT, etc.) but before any ray buffers.
    // This gives a stable baseline that does not shrink on subsequent run() calls
    // due to the already-allocated (and reused) ray buffers being counted as used.
    const size_t mem_free = m_mem_free_post_setup;

    // Reserve 20 % headroom for OptiX internal allocations, memory
    // fragmentation, and any other transient allocations during launch.
    constexpr double kUsableFraction = 0.80;
    const size_t usable_bytes = static_cast<size_t>(
        static_cast<double>(mem_free) * kUsableFraction);

    // Per-ray device memory charged by allocate_device_buffers() and
    // allocate_compaction_scratch():
    //   hit_buffer      DEFAULT_MAX_TRACE_DEPTH * sizeof(HitRecord)  -- trace output
    //   d_compacted     DEFAULT_MAX_TRACE_DEPTH * sizeof(HitRecord)  -- worst-case compacted copy
    //   sun_dir_buffer  sizeof(float3)                       -- sun ray direction
    //   curand states   sizeof(curandState)                  -- RNG state
    //   d_offsets       sizeof(uint64_t)                     -- compaction prefix sum / global ray IDs
    //   d_count         sizeof(uint8_t)                      -- compaction hit count (bounded by DEFAULT_MAX_TRACE_DEPTH <= 255)
    //   d_has_hit       sizeof(uint8_t)                      -- per-ray hit flag
    const size_t bytes_per_ray =
        2u * m_max_ray_depth * sizeof(HitRecord) + sizeof(float3) + sizeof(curandState) + sizeof(uint64_t) + 2u * sizeof(uint8_t);

    const uint_fast64_t computed =
        (bytes_per_ray > 0) ? static_cast<uint_fast64_t>(usable_bytes / bytes_per_ray) : 0u;

    // Cap at int max / m_max_ray_depth (OptiX launch width is signed int).
    uint_fast64_t batch_size = std::min(
        computed,
        static_cast<uint_fast64_t>(std::numeric_limits<int>::max() / m_max_ray_depth));

    if (m_verbose)
    {
        std::cout << "automatic_batch_size:"
                  << " free=" << mem_free / (1024.0 * 1024.0) << " MB"
                  << ", usable=" << usable_bytes / (1024.0 * 1024.0) << " MB"
                  << ", bytes_per_ray=" << bytes_per_ray
                  << ", batch_size=" << batch_size << "\n";
    }

    return batch_size;
}

uint_fast64_t SolTraceSystem::determine_batch_size() const
{
    // Estimates number of rays that can be traced in a single batch based on
    // available GPU memory.
    uint_fast64_t batch_size = automatic_batch_size();

    if (m_batch_size > 0)
    {
        if (m_batch_size > batch_size && batch_size > 0)
        {
            std::cerr << "[SolTraceSystem] WARNING: user-supplied batch_size ("
                      << m_batch_size
                      << ") exceeds the GPU-memory-safe automatic batch size ("
                      << batch_size
                      << "). This may cause device out-of-memory errors or "
                         "degraded GPU performance.\n";
        }
        batch_size = m_batch_size;
    }
    else
    {
        // Take the smaller of the automatic batch_size and number of rays?
        batch_size = batch_size > 0 ? std::min(batch_size, m_number_of_rays) : m_number_of_rays;
    }

    return batch_size;
}

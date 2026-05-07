#include "soltrace_system.h"
#include "geometry_manager.h"
#include "data_manager.h"
#include "pipeline_manager.h"
#include "soltrace_type.h"
#include "CspElement.h"
#include "timer.h"
#include "soltrace_constants.h"
#include "../../../../../simulation_data/simdata_io.hpp"
#include "../../../../../simulation_data/solar_position_calculators/basic_sun_position.hpp"

#include "utils/util_record.hpp"
#include "utils/util_check.hpp"
#include "utils/math_util.h"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cstring>

#include <optix_function_table_definition.h>
#include <optix_stubs.h>

using namespace OptixCSP;

// TODO: optix related type should go into one header file
// i don't know what we should put here in hit group record, material is handled through a global array
// we can leave this empty for now ...
// note that this is has to be per optical entity type.
typedef Record<OptixCSP::HitGroupData> HitGroupRecord;

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
    std::cout << "hit_point_buffer   : " << params.hit_point_buffer << std::endl;
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

SolTraceSystem::SolTraceSystem()
    : m_number_of_rays(0),
      m_max_number_of_rays(0),
      m_verbose(false),
      m_mem_free_before(0),
      m_mem_free_after(0),
      m_optical_errors(false),
      m_include_sun_shape_errors(false),
      m_timer_setup(),
      m_timer_trace(),
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

SolTraceSystem::~SolTraceSystem()
{
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

    size_t mem_total;
    cudaMemGetInfo(&m_mem_free_before, &mem_total);
    m_timer_setup.start();

    // set up input related to sun
    glm::dvec3 sun_vec = m_sun->get_position();
    glm::dvec3 sun_vec_norm = sun_vec;
    sun_vec_norm = glm::normalize(sun_vec_norm);

    data_manager->launch_params_H.sun_vector = make_float3(static_cast<float>(sun_vec_norm[0]),
                                                           static_cast<float>(sun_vec_norm[1]), static_cast<float>(sun_vec_norm[2]));

    // Set generation type
    switch (m_sun->get_gen_type())
    {
        case(SolTrace::Data::GenType::RANDOM):
            data_manager->launch_params_H.sun_gen_type = OptixCSP::GenType::RANDOM;
            break;
        case(SolTrace::Data::GenType::HALTON):
            data_manager->launch_params_H.sun_gen_type = OptixCSP::GenType::HALTON;
            break;
        default:
            data_manager->launch_params_H.sun_gen_type = OptixCSP::GenType::UNKNOWN;
    }

    // Assign sun shape parameters (if necessary)
    data_manager->launch_params_H.include_sun_shape_errors = this->m_include_sun_shape_errors;
    data_manager->allocateSunUserData({}, {});  // Clear sun user data
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

    Timer AABB_timer;
    AABB_timer.start();
    geometry_manager->collect_geometry_info(m_element_list, data_manager->launch_params_H);
    AABB_timer.stop();

    Timer geometry_timer;
    geometry_timer.start();
    geometry_manager->create_geometries(data_manager->launch_params_H);
    geometry_timer.stop();

    // Pipeline setup.
    Timer pipeline_timer;
    pipeline_timer.start();
    pipeline_manager->createPipeline();
    pipeline_timer.stop();

    Timer sbt_timer;
    sbt_timer.start();
    create_shader_binding_table();
    sbt_timer.stop();

    // seed for randomization
    data_manager->launch_params_H.sun_dir_seed = m_seed;
    data_manager->launch_params_H.optical_errors = m_optical_errors;

    // Create a CUDA stream for asynchronous operations.
    CUDA_CHECK(cudaStreamCreate(&m_state.stream));

    // Link the GAS handle.
    data_manager->launch_params_H.handle = m_state.gas_handle;
    data_manager->allocateGeometryDataArray(geometry_manager->get_geometry_data_array());
    data_manager->allocateMaterialDataArray(geometry_manager->get_material_data_array_front(),
                                            geometry_manager->get_material_data_array_back());
    
    if (m_verbose)
    {
        std::cout << "Time to compute AABB: " << AABB_timer.get_time_sec() << " seconds" << std::endl;
        std::cout << "Time to create geometries: " << geometry_timer.get_time_sec() << " seconds" << std::endl;
        std::cout << "Time to create pipeline: " << pipeline_timer.get_time_sec() << " seconds" << std::endl;
        std::cout << "Time to create SBT: " << sbt_timer.get_time_sec() << " seconds" << std::endl;

        print_launch_params();
    }
        

    data_manager->allocateLaunchParams();
    m_timer_setup.stop();
}

void SolTraceSystem::run()
{

    // Initialize results vectors
    m_hp_vec.clear();
    m_raynumber_vec.clear();
    m_element_id_vec.clear();
    m_hit_type_vec.clear();
    uint_fast64_t N_ray_hit = 0;
    uint_fast64_t N_ray_gen = 0;

    while (N_ray_hit < m_number_of_rays && N_ray_gen < m_max_number_of_rays)
    {
        // Update ray offset (pushed to device in setup_device_buffer)
        data_manager->launch_params_H.ray_offset = N_ray_gen;

        // Allocate buffer (sets data_manager->launch_params_H buffer)
        setup_device_buffer();

        int width = data_manager->launch_params_H.width;
        int height = data_manager->launch_params_H.height;

        size_t m_mem_free_after;
	    size_t mem_total;
        cudaMemGetInfo(&m_mem_free_after, &mem_total);

        if(m_verbose)
            std::cout << "Memory used by launch: " << (m_mem_free_before - m_mem_free_after) / (1024.0 * 1024.0) << " MB\n";

        m_timer_trace.start();
        // Launch the simulation.
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

        // Collect results
        get_buffer_results(m_hp_vec, m_raynumber_vec, m_element_id_vec, m_hit_type_vec,
                           m_sunraynumber_vec);
        N_ray_hit = m_raynumber_vec.empty() ? 0 : m_raynumber_vec.back();
        N_ray_gen += width;
    }

    // Trim excess rays
    if (N_ray_hit > m_number_of_rays)
    {
        while (m_raynumber_vec.back() > m_number_of_rays)
        {
            m_hp_vec.pop_back();
            m_raynumber_vec.pop_back();
            m_element_id_vec.pop_back();
            m_hit_type_vec.pop_back();
            m_sunraynumber_vec.pop_back();
        }
    }

    m_timer_trace.stop();
}

void SolTraceSystem::update()
{

    const int N_slots = data_manager->launch_params_H.width * data_manager->launch_params_H.height * data_manager->launch_params_H.max_depth;
    const size_t hit_point_buffer_size = N_slots * sizeof(float4);
    const size_t element_id_size = N_slots * sizeof(int32_t);
    const size_t hit_type_buffer_size = N_slots * sizeof(uint8_t);

    // update aabb and sun plane accordingly
    geometry_manager->update_geometry_info(m_element_list, data_manager->launch_params_H);

    // update data on the device
    data_manager->updateGeometryDataArray(geometry_manager->get_geometry_data_array());
    CUDA_CHECK(cudaMemset(data_manager->launch_params_H.hit_point_buffer, 0, hit_point_buffer_size));
    CUDA_CHECK(cudaMemset(data_manager->launch_params_H.element_id_buffer, kElementIdBuffer, element_id_size));
    CUDA_CHECK(cudaMemset(data_manager->launch_params_H.hit_type_buffer, HitType::HIT_UNASSIGNED, hit_type_buffer_size));

    data_manager->updateLaunchParams();
}

void SolTraceSystem::get_hp_output(std::vector<float4> &hp_vec,
                                   std::vector<uint_fast64_t> &raynumber_vec,
                                   std::vector<int32_t> &element_id_vec,
                                   std::vector<uint8_t> &hit_type_vec)
{
    hp_vec = m_hp_vec;
    raynumber_vec = m_raynumber_vec;
    element_id_vec = m_element_id_vec;
    hit_type_vec = m_hit_type_vec;
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
    CUDA_CHECK(cudaFree(reinterpret_cast<void *>(data_manager->launch_params_H.hit_point_buffer)));
    CUDA_CHECK(cudaFree(reinterpret_cast<void *>(data_manager->launch_params_H.element_id_buffer)));
    CUDA_CHECK(cudaFree(reinterpret_cast<void *>(data_manager->launch_params_H.hit_type_buffer)));
    CUDA_CHECK(cudaFree(reinterpret_cast<void *>(data_manager->launch_params_H.sun_dir_buffer)));

    data_manager->launch_params_H.hit_point_buffer = nullptr;
    data_manager->launch_params_H.element_id_buffer = nullptr;
    data_manager->launch_params_H.hit_type_buffer = nullptr;
    data_manager->launch_params_H.sun_dir_buffer = nullptr;
    m_hit_point_buffer_size_allocated = 0;
    m_element_id_buffer_size_allocated = 0;
    m_hit_type_buffer_size_allocated = 0;
    m_sun_dir_buffer_size_allocated = 0;

    data_manager->cleanup();

    m_hp_output_buffer_host.clear();
    m_hp_output_buffer_host.shrink_to_fit();
    m_element_id_buffer_host.clear();
    m_element_id_buffer_host.shrink_to_fit();
    m_hit_type_buffer_host.clear();
    m_hit_type_buffer_host.shrink_to_fit();

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
}

void SolTraceSystem::reset()
{
    clean_up();

    m_element_list.clear();
    m_hp_vec.clear();
    m_raynumber_vec.clear();
    m_element_id_vec.clear();
    m_hit_type_vec.clear();
    m_sunraynumber_vec.clear();

    m_hp_output_buffer_host.clear();
    m_element_id_buffer_host.clear();
    m_hit_type_buffer_host.clear();

    m_sun = nullptr;
    m_number_of_rays = 0;
    m_max_number_of_rays = 0;
}

// Create and configure the Shader Binding Table (SBT).
// The SBT is a crucial data structure in OptiX that links geometry and ray types
// with their corresponding programs (ray generation, miss, and hit group).
void SolTraceSystem::create_shader_binding_table()
{

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
            // OptixProgramGroup program_group_handle = nullptr;
            // SurfaceApertureMap map = {};

            // switch (my_type)
            // {
            // case OptixCSP::OpticalEntityType::RECTANGLE_FLAT:
            //     map = {SurfaceType::FLAT, ApertureType::RECTANGLE};
            //     program_group_handle = pipeline_manager->getElementProgram(map);
            //     hitgroup_records_list[i].data.material_data = {0.875425, 0, 0, 0};
            //     printf("RECTANGLE_FLAT, program group address: %p \n", program_group_handle);

            //     break;

            // case OptixCSP::OpticalEntityType::RECTANGLE_PARABOLIC:
            //     map = {SurfaceType::PARABOLIC, ApertureType::RECTANGLE};
            //     program_group_handle = pipeline_manager->getElementProgram(map);
            //     hitgroup_records_list[i].data.material_data = {0.875425, 0, 0, 0};
            //     printf("RECTANGLE_PARABOLIC, program group address: %p \n", program_group_handle);

            //     break;

            // case OptixCSP::OpticalEntityType::CYLINDRICAL:
            //     map = {SurfaceType::CYLINDER, ApertureType::RECTANGLE};
            //     program_group_handle = pipeline_manager->getElementProgram(map);
            //     hitgroup_records_list[i].data.material_data = {0.95, 0, 0, 0};
            //     printf("CYLINDRICAL, program group address: %p \n", program_group_handle);

            //     break;

            // case OptixCSP::OpticalEntityType::TRIANGLE_FLAT:
            //     map = {SurfaceType::FLAT, ApertureType::TRIANGLE};
            //     program_group_handle = pipeline_manager->getElementProgram(map);
            //     hitgroup_records_list[i].data.material_data = {0.95, 0, 0, 0};
            //     printf("FLAT_TRIANGLE, program group address: %p \n", program_group_handle);

            //     break;

            // case OptixCSP::OpticalEntityType::QUADRILATERAL_FLAT:
            //     ma = {SurfaceType::FLAT, ApertureType::QUADRILATERAL};
            //     program_group_handle = pipeline_manager->getElementProgram(map);
            //     hitgroup_records_list[i].data.material_data = {0.875425, 0, 0, 0};
            //     printf("FLAT_QUADRILATERAL, program group address: %p \n", program_group_handle);

            //     break;

            // default:
            //     std::cerr << "Unknown OpticalEntityType: " << my_type << std::endl;
            // }

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

void SolTraceSystem::setup_device_buffer()
{
    // Initialize launch params
    data_manager->launch_params_H.width = m_number_of_rays;
    data_manager->launch_params_H.height = 1;
    data_manager->launch_params_H.max_depth = MAX_TRACE_DEPTH;

    const size_t hit_point_buffer_size = data_manager->launch_params_H.width * data_manager->launch_params_H.height * sizeof(float4) * data_manager->launch_params_H.max_depth;
    const size_t element_id_size = data_manager->launch_params_H.width * data_manager->launch_params_H.height * sizeof(int32_t) * data_manager->launch_params_H.max_depth;
    const size_t hit_type_size = data_manager->launch_params_H.width * data_manager->launch_params_H.height * sizeof(uint8_t) * data_manager->launch_params_H.max_depth;
    const size_t sun_dir_size = data_manager->launch_params_H.width * data_manager->launch_params_H.height * sizeof(float3);

    if (data_manager->launch_params_H.hit_point_buffer == nullptr || m_hit_point_buffer_size_allocated != hit_point_buffer_size)
    {
        if (data_manager->launch_params_H.hit_point_buffer != nullptr)
            CUDA_CHECK(cudaFree(reinterpret_cast<void *>(data_manager->launch_params_H.hit_point_buffer)));
        CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&data_manager->launch_params_H.hit_point_buffer), hit_point_buffer_size));
        m_hit_point_buffer_size_allocated = hit_point_buffer_size;
    }
    CUDA_CHECK(cudaMemset(data_manager->launch_params_H.hit_point_buffer, 0, hit_point_buffer_size));

    if (data_manager->launch_params_H.element_id_buffer == nullptr || m_element_id_buffer_size_allocated != element_id_size)
    {
        if (data_manager->launch_params_H.element_id_buffer != nullptr)
            CUDA_CHECK(cudaFree(reinterpret_cast<void *>(data_manager->launch_params_H.element_id_buffer)));
        CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&data_manager->launch_params_H.element_id_buffer), element_id_size));
        m_element_id_buffer_size_allocated = element_id_size;
    }
    CUDA_CHECK(cudaMemset(data_manager->launch_params_H.element_id_buffer, kElementIdBuffer, element_id_size));

    if (data_manager->launch_params_H.hit_type_buffer == nullptr || m_hit_type_buffer_size_allocated != hit_type_size)
    {
        if (data_manager->launch_params_H.hit_type_buffer != nullptr)
            CUDA_CHECK(cudaFree(reinterpret_cast<void *>(data_manager->launch_params_H.hit_type_buffer)));
        CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&data_manager->launch_params_H.hit_type_buffer), hit_type_size));
        m_hit_type_buffer_size_allocated = hit_type_size;
    }
    CUDA_CHECK(cudaMemset(data_manager->launch_params_H.hit_type_buffer, HitType::HIT_UNASSIGNED, hit_type_size));

    if (data_manager->launch_params_H.sun_dir_buffer == nullptr || m_sun_dir_buffer_size_allocated != sun_dir_size)
    {
        if (data_manager->launch_params_H.sun_dir_buffer != nullptr)
            CUDA_CHECK(cudaFree(reinterpret_cast<void *>(data_manager->launch_params_H.sun_dir_buffer)));
        CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&data_manager->launch_params_H.sun_dir_buffer), sun_dir_size));
        m_sun_dir_buffer_size_allocated = sun_dir_size;
    }
    CUDA_CHECK(cudaMemset(data_manager->launch_params_H.sun_dir_buffer, 0, sun_dir_size));

    const unsigned int num_rng_states = static_cast<unsigned int>(data_manager->launch_params_H.width * data_manager->launch_params_H.height);
    data_manager->ensureCurandStates(
        num_rng_states,
        data_manager->launch_params_H.sun_dir_seed,
        data_manager->launch_params_H.ray_offset,
        m_state.stream);

    data_manager->updateLaunchParams();
}

// Collects results from device buffer
// only keeps rays that hit elements
void SolTraceSystem::get_buffer_results(std::vector<float4> &hp_vec, std::vector<uint_fast64_t> &raynumber_vec,
                                        std::vector<int32_t> &element_id_vec, std::vector<uint8_t> &hit_type_vec,
                                        std::vector<uint_fast64_t> &sunraynumber_vec)
{
    const int max_depth = data_manager->launch_params_H.max_depth;
    const int num_rays = data_manager->launch_params_H.width * data_manager->launch_params_H.height;
    const int output_size = data_manager->launch_params_H.width * data_manager->launch_params_H.height * data_manager->launch_params_H.max_depth;

    if (static_cast<int>(m_hp_output_buffer_host.size()) != output_size)
        m_hp_output_buffer_host.resize(output_size);
    if (static_cast<int>(m_element_id_buffer_host.size()) != output_size)
        m_element_id_buffer_host.resize(output_size);
    if (static_cast<int>(m_hit_type_buffer_host.size()) != output_size)
        m_hit_type_buffer_host.resize(output_size);

    CUDA_CHECK(cudaMemcpy(m_hp_output_buffer_host.data(), data_manager->launch_params_H.hit_point_buffer, output_size * sizeof(float4), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(m_element_id_buffer_host.data(), data_manager->launch_params_H.element_id_buffer, output_size * sizeof(int32_t), cudaMemcpyDeviceToHost));
    CUDA_CHECK(cudaMemcpy(m_hit_type_buffer_host.data(), data_manager->launch_params_H.hit_type_buffer, output_size * sizeof(uint8_t), cudaMemcpyDeviceToHost));

    // Loop through each buffer slot
    uint_fast64_t ray_number = raynumber_vec.empty() ? 0 : raynumber_vec.back();
    uint_fast64_t sunray_number = sunraynumber_vec.empty() ? 0 : sunraynumber_vec.back();
    for (int i = 0; i < output_size; ++i)
    {

        // Get hit type
        const uint8_t &hit_type = m_hit_type_buffer_host[i];

        // Skip if empty
        if (hit_type < HitType::HIT_CREATE || hit_type > HitType::HIT_EXIT)
        {
            continue;
        }

        // If new ray, check if previous ray hit anything
        if (hit_type == HitType::HIT_CREATE)
        {
            // Remove last ray if it has no hits
            if (!hit_type_vec.empty() && hit_type_vec.back() == HitType::HIT_CREATE)
            {
                hp_vec.pop_back();
                raynumber_vec.pop_back();
                hit_type_vec.pop_back();
                element_id_vec.pop_back();
                sunraynumber_vec.pop_back();
                ray_number--;
            }

            // New ray
            ray_number++;

            // Sun ray number always increments, even if no hit
            sunray_number++;
        }

        // Get hit record, element_id
        const float4 &hit_record = m_hp_output_buffer_host[i]; // [depth, pos x, pos y, pos z]
        const int32_t &element_id = m_element_id_buffer_host[i];

        // Collect results
        hp_vec.push_back(hit_record);
        raynumber_vec.push_back(ray_number);
        hit_type_vec.push_back(hit_type);
        element_id_vec.push_back(element_id);
        sunraynumber_vec.push_back(sunray_number);
    }

    // Remove last ray if it is only CREATE
    if (!hit_type_vec.empty() && hit_type_vec.back() == HitType::HIT_CREATE)
    {
        hp_vec.pop_back();
        raynumber_vec.pop_back();
        element_id_vec.pop_back();
        hit_type_vec.pop_back();
        sunraynumber_vec.pop_back();
    }

    return;
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

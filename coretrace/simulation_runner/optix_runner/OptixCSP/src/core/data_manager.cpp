#include "data_manager.h"
#include "soltrace_system.h"
#include "utils/util_check.hpp"
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include "sun_utils.h"

using namespace OptixCSP;

dataManager::dataManager()
	: launch_params_D(nullptr),
	  geometry_data_array_D(nullptr),
	  material_data_array_front_D(nullptr),
	  material_data_array_back_D(nullptr),
	  sun_user_angle_D(nullptr),
	  sun_user_intensity_D(nullptr),
	  sun_user_capacity(0),
	  rng_states_D(nullptr),
	  rng_states_capacity(0),
	  depth_exceeded_count_D(nullptr)
{

	// Initialize launch parameters with default values
	launch_params_H.width = 10;
	launch_params_H.height = 1;
	launch_params_H.max_depth = 5;
	launch_params_H.ray_offset = 0;

	// launch_params_H.hit_point_buffer = nullptr;
	launch_params_H.hit_buffer = nullptr;
	launch_params_H.sun_dir_buffer = nullptr;
	launch_params_H.rng_states = nullptr;
	// launch_params_H.element_id_buffer = nullptr;
	// launch_params_H.hit_type_buffer = nullptr;
	launch_params_H.sun_vector = make_float3(0.0f, 0.0f, 10.0f);
	launch_params_H.sun_shape = OptixCSP::SunShape::UNKNOWN;
	launch_params_H.include_sun_shape_errors = false;
	launch_params_H.sigma = 0.0f;
	launch_params_H.half_width = 0.0f;
	launch_params_H.buie_kappa = 0.0f;
	launch_params_H.buie_gamma = 0.0f;
	launch_params_H.sun_max_angle = 0.0f;
	launch_params_H.sun_max_intensity = 0.0f;
	launch_params_H.sun_gen_type = OptixCSP::GenType::UNKNOWN;
	launch_params_H.sun_user_angle = nullptr;
	launch_params_H.sun_user_intensity = nullptr;
	launch_params_H.sun_user_capacity = 0;

	launch_params_H.sun_v0 = make_float3(0.0f, 0.0f, 0.0f);
	launch_params_H.sun_v1 = make_float3(0.0f, 0.0f, 0.0f);
	launch_params_H.sun_v2 = make_float3(0.0f, 0.0f, 0.0f);
	launch_params_H.sun_v3 = make_float3(0.0f, 0.0f, 0.0f);

	launch_params_H.optical_errors = false;
	launch_params_H.material_data_array_front = nullptr;
	launch_params_H.material_data_array_back = nullptr;
	launch_params_H.sun_dir_seed = 0ULL;

	launch_params_H.d_depth_exceeded_count = nullptr;

	launch_params_H.geometry_data_array = nullptr;
	launch_params_H.handle = OptixTraversableHandle{};

	
}

dataManager::~dataManager() noexcept
{
	try
	{
		cleanup();
	}
	catch (const std::exception &error)
	{
		std::cerr << "[OptixRunner] Data cleanup failed during destruction: "
		          << error.what() << '\n';
	}
	catch (...)
	{
		std::cerr << "[OptixRunner] Data cleanup failed during destruction with an unknown error.\n";
	}
}

OptixCSP::LaunchParams *dataManager::getDeviceLaunchParams() const { return launch_params_D; }

void dataManager::allocateLaunchParams()
{
	if (launch_params_D)
	{
		CUDA_CHECK(cudaFree(launch_params_D));
		launch_params_D = nullptr;
	}
	CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&launch_params_D), sizeof(LaunchParams)));
}

void dataManager::updateLaunchParams()
{
	CUDA_CHECK(cudaMemcpy(launch_params_D, &launch_params_H, sizeof(LaunchParams), cudaMemcpyHostToDevice));
}

void dataManager::ensureCurandStates(
	unsigned int num_states,
	unsigned long long seed,
	unsigned int sequence_offset,
	cudaStream_t stream)
{

	if (num_states == 0)
	{
		launch_params_H.rng_states = nullptr;
		return;
	}

	if (rng_states_capacity < num_states)
	{
		if (rng_states_D != nullptr)
		{
			CUDA_CHECK(cudaFree(rng_states_D));
		}

		CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&rng_states_D), num_states * sizeof(curandState)));
		rng_states_capacity = num_states;
	}

	initialize_curand_states_on_gpu(rng_states_D, num_states, seed, sequence_offset, stream);
	launch_params_H.rng_states = rng_states_D;
}

void dataManager::ensureDepthExceededCounter()
{
	if (depth_exceeded_count_D == nullptr)
	{
		CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&depth_exceeded_count_D), sizeof(uint64_t)));
		CUDA_CHECK(cudaMemset(depth_exceeded_count_D, 0, sizeof(uint64_t)));
		launch_params_H.d_depth_exceeded_count = depth_exceeded_count_D;
	}
}

void dataManager::allocateGeometryDataArray(std::vector<GeometryDataST> geometry_data_array_H)
{
	if (geometry_data_array_D)
	{
		CUDA_CHECK(cudaFree(geometry_data_array_D));
		geometry_data_array_D = nullptr;
	}

	CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&geometry_data_array_D),
						  geometry_data_array_H.size() * sizeof(GeometryDataST)));

	CUDA_CHECK(cudaMemcpy(geometry_data_array_D, geometry_data_array_H.data(),
						  geometry_data_array_H.size() * sizeof(GeometryDataST), cudaMemcpyHostToDevice));
	// make sure launch_params_H is updated with the new geometry data array
	launch_params_H.geometry_data_array = geometry_data_array_D;
}

void dataManager::updateGeometryDataArray(std::vector<GeometryDataST> geometry_data_array_H)
{

	if (geometry_data_array_D == nullptr)
	{
		throw std::runtime_error("Geometry data array is not allocated.");
	}

	CUDA_CHECK(cudaMemcpy(geometry_data_array_D, geometry_data_array_H.data(),
						  geometry_data_array_H.size() * sizeof(GeometryDataST), cudaMemcpyHostToDevice));

	// launch_params_H.geometry_data_array = geometry_data_array_D;
}

void dataManager::allocateMaterialDataArray(std::vector<MaterialData> material_data_array_front_H,
											std::vector<MaterialData> material_data_array_back_H)
{
	if (material_data_array_front_D)
	{
		CUDA_CHECK(cudaFree(material_data_array_front_D));
		material_data_array_front_D = nullptr;
	}
	if (material_data_array_back_D)
	{
		CUDA_CHECK(cudaFree(material_data_array_back_D));
		material_data_array_back_D = nullptr;
	}

	CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&material_data_array_front_D),
						  material_data_array_front_H.size() * sizeof(MaterialData)));

	CUDA_CHECK(cudaMemcpy(material_data_array_front_D, material_data_array_front_H.data(),
						  material_data_array_front_H.size() * sizeof(MaterialData), cudaMemcpyHostToDevice));
	// make sure launch_params_H is updated with the new geometry data array
	launch_params_H.material_data_array_front = material_data_array_front_D;

	CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&material_data_array_back_D),
						  material_data_array_back_H.size() * sizeof(MaterialData)));

	CUDA_CHECK(cudaMemcpy(material_data_array_back_D, material_data_array_back_H.data(),
						  material_data_array_back_H.size() * sizeof(MaterialData), cudaMemcpyHostToDevice));
	// make sure launch_params_H is updated with the new geometry data array
	launch_params_H.material_data_array_back = material_data_array_back_D;
}

void dataManager::updateMaterialDataArray(std::vector<MaterialData> material_data_array_H)
{
	if (material_data_array_front_D == nullptr)
	{
		throw std::runtime_error("Not implemented yet ... does material data change??");
	}

	if (material_data_array_back_D == nullptr)
	{
		throw std::runtime_error("Not implemented yet ... does material data change??");
	}
}

void dataManager::allocateSunUserData(const std::vector<float>& user_angle,
	                                  const std::vector<float>& user_intensity)
{
	if (user_angle.size() != user_intensity.size())
	{
		throw std::runtime_error("User-defined sun angle/intensity vectors must have the same size.");
	}

	if (sun_user_angle_D != nullptr)
	{
		CUDA_CHECK(cudaFree(sun_user_angle_D));
		sun_user_angle_D = nullptr;
	}

	if (sun_user_intensity_D != nullptr)
	{
		CUDA_CHECK(cudaFree(sun_user_intensity_D));
		sun_user_intensity_D = nullptr;
	}

	launch_params_H.sun_user_angle = nullptr;
	launch_params_H.sun_user_intensity = nullptr;
	launch_params_H.sun_user_capacity = 0;
	sun_user_capacity = 0;

	if (user_angle.empty())
	{
		return;
	}

	CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&sun_user_angle_D),
		                  user_angle.size() * sizeof(float)));
	CUDA_CHECK(cudaMalloc(reinterpret_cast<void **>(&sun_user_intensity_D),
		                  user_intensity.size() * sizeof(float)));

	CUDA_CHECK(cudaMemcpy(sun_user_angle_D, user_angle.data(),
		                  user_angle.size() * sizeof(float), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaMemcpy(sun_user_intensity_D, user_intensity.data(),
		                  user_intensity.size() * sizeof(float), cudaMemcpyHostToDevice));

	launch_params_H.sun_user_angle = sun_user_angle_D;
	launch_params_H.sun_user_intensity = sun_user_intensity_D;
	launch_params_H.sun_user_capacity = static_cast<int>(user_angle.size());
	sun_user_capacity = user_angle.size();
}

void dataManager::cleanup() {
	if (launch_params_D) {
		CUDA_CHECK(cudaFree(launch_params_D));
		launch_params_D = nullptr;
	}

	if (geometry_data_array_D) {
		CUDA_CHECK(cudaFree(geometry_data_array_D));
		geometry_data_array_D = nullptr;
	}
	launch_params_H.geometry_data_array = nullptr;

	if (material_data_array_front_D) {
		CUDA_CHECK(cudaFree(material_data_array_front_D));
		material_data_array_front_D = nullptr;
	}
	launch_params_H.material_data_array_front = nullptr;

	if (material_data_array_back_D) {
		CUDA_CHECK(cudaFree(material_data_array_back_D));
		material_data_array_back_D = nullptr;
	}
	launch_params_H.material_data_array_back = nullptr;

	if (sun_user_angle_D != nullptr) {
		CUDA_CHECK(cudaFree(sun_user_angle_D));
		sun_user_angle_D = nullptr;
	}
	launch_params_H.sun_user_angle = nullptr;

	if (sun_user_intensity_D != nullptr) {
		CUDA_CHECK(cudaFree(sun_user_intensity_D));
		sun_user_intensity_D = nullptr;
	}
	launch_params_H.sun_user_intensity = nullptr;
	launch_params_H.sun_user_capacity = 0;
	sun_user_capacity = 0;

	if (rng_states_D != nullptr) {
		CUDA_CHECK(cudaFree(rng_states_D));
		rng_states_D = nullptr;
		rng_states_capacity = 0;
	}
	launch_params_H.rng_states = nullptr;

	if (depth_exceeded_count_D != nullptr) {
		CUDA_CHECK(cudaFree(depth_exceeded_count_D));
		depth_exceeded_count_D = nullptr;
	}
	launch_params_H.d_depth_exceeded_count = nullptr;

	launch_params_H.handle = OptixTraversableHandle{};
}

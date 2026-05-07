#pragma once

#include "shaders/Soltrace.h"
#include "CspElement.h"
#include <vector>
#include <cuda_runtime_api.h>

namespace OptixCSP
{

    // Class to manage data on the host and device.
    class dataManager
    {
    public:
        // TODO: move this to private member for best practice
        // Host copy of launch parameters.
        OptixCSP::LaunchParams launch_params_H;
        // Device pointer to launch parameters.
        OptixCSP::LaunchParams *launch_params_D;

        // device pointer to geometry data
        GeometryDataST *geometry_data_array_D;

        // device pointer to material data
        MaterialData *material_data_array_front_D;
        MaterialData *material_data_array_back_D;

        // CURAND state memory on device
        curandState *rng_states_D;
        size_t rng_states_capacity;

        float *sun_user_angle_D;
        float *sun_user_intensity_D;
        size_t sun_user_capacity;

        dataManager();
        ~dataManager();

        void cleanup();

        OptixCSP::LaunchParams *getDeviceLaunchParams() const;

        void allocateLaunchParams();

        void updateLaunchParams();

        // create geometry_data_array_D on the device
        // then launch_params_D.geometry_data_array = geometry_data_array_D gets a copy.
        void allocateGeometryDataArray(std::vector<GeometryDataST> geometry_data_array);

        // update geometry_data_array_D on the device
        // then launch_params_D.geometry_data_array = geometry_data_array_D gets a copy.
        void updateGeometryDataArray(std::vector<GeometryDataST> geometry_data_array_H);

        // create material_data_array_D on the device
        // then launch_params_D.material_data_array = material_data_array_D gets a copy.
        void allocateMaterialDataArray(std::vector<MaterialData> material_data_array_front,
                                       std::vector<MaterialData> material_data_array_back);

        // update material_data_array_D on the device
        // then launch_params_D.material_data_array = material_data_array_D gets a copy.
        void updateMaterialDataArray(std::vector<MaterialData> material_data_array_H);

        void allocateSunUserData(const std::vector<float>& user_angle,
                                 const std::vector<float>& user_intensity);

        void ensureCurandStates(unsigned int num_states,
                                unsigned long long seed,
                                unsigned int sequence_offset,
                                cudaStream_t stream);
    };
}

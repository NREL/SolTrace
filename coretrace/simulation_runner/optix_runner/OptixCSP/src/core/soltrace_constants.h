#pragma once
                
#include <cstdint>

namespace OptixCSP {

    // Mirror of SolTrace::Data::SunShape for OptiX/device code.
    // Keep values in sync with SolTrace::Data::SunShape in ray_source.hpp.
    enum class SunShape : uint8_t {
        GAUSSIAN = 0,
        PILLBOX = 1,
        LIMBDARKENED = 2,
        BUIE_CSR = 3,
        USER_DEFINED = 4,
        UNKNOWN = 5
    };

    // Hit type
    // TODO: Replace this with RayEvent from simulation_result.hpp?
    enum HitType : uint8_t {
        HIT_UNASSIGNED = 0,
        HIT_CREATE = 1,
        HIT_ABSORB = 2,
        HIT_REFLECT = 3,
        HIT_TRANSMIT = 4,
        HIT_VIRTUAL = 5,
        HIT_EXIT = 6,
        HIT_UNKNOWN = 255
    };

    enum OpticalDistribution : uint8_t {
        OPT_NONE = 0,
        OPT_GAUSSIAN = 1,
        OPT_PILLBOX = 2,
        OPT_DIFFUSE = 3,
        OPT_UNKNOWN = 255
    };

    enum class GenType : uint8_t {
        RANDOM,
        HALTON,
        UNKNOWN
    };

    // Element ID consts
    constexpr int32_t kElementIdBuffer = 0;
    constexpr int32_t kElementIdRayGen = -1;
    constexpr int32_t kElementIdUnassigned = -2;

}
# ──────────────────────────────────────────────────────────────
#  src/CMakeLists.txt   – builds the OptixCSP core library
# ──────────────────────────────────────────────────────────────

include(GNUInstallDirs)

# ---------------------------------------------------------------------------
# core target
# ---------------------------------------------------------------------------
add_library(OptixCSP_core STATIC)

# Link against CUDA and simulation_data (simdata)
target_link_libraries(OptixCSP_core PRIVATE
    CUDA::cuda_driver
    CUDA::cudart
    simdata
)

find_package(CUDAToolkit REQUIRED)
find_package(OptiX REQUIRED)

# OptixCSP's FindOptiX sets OptiX_INCLUDE, normalize to OptiX_INCLUDE for consistency
if(NOT OptiX_INCLUDE AND OptiX_INCLUDE_DIR)
    set(OptiX_INCLUDE ${OptiX_INCLUDE_DIR})
endif()
message(STATUS "OptiX include directory: ${OptiX_INCLUDE}")

# Build for a declared minimum GPU rather than the build host. Release builders
# are often headless, and host detection would make artifacts depend on whatever
# GPU happened to be installed on the build machine.
set(
    SOLTRACE_OPTIX_COMPUTE_CAPABILITY "75"
    CACHE STRING
    "Minimum CUDA compute capability for the OptiX runner, without a decimal point"
)
if(NOT SOLTRACE_OPTIX_COMPUTE_CAPABILITY MATCHES "^[0-9]+$")
    message(FATAL_ERROR
        "SOLTRACE_OPTIX_COMPUTE_CAPABILITY must be numeric, for example 75"
    )
endif()
message(STATUS
    "OptiX minimum CUDA compute capability: ${SOLTRACE_OPTIX_COMPUTE_CAPABILITY}"
)


file(GLOB_RECURSE CORE_SRC  CONFIGURE_DEPENDS  core/*.cpp core/*.cu)
file(GLOB_RECURSE CORE_HDR  CONFIGURE_DEPENDS  core/*.hpp core/*.h utils/*.hpp utils/*.h shaders/*.h)

target_sources(OptixCSP_core PRIVATE ${CORE_SRC} ${CORE_HDR})

# A numeric CMake CUDA architecture emits both native code for the minimum GPU
# and virtual PTX that NVIDIA drivers can JIT for newer GPU architectures.
set_property(
    TARGET OptixCSP_core
    PROPERTY CUDA_ARCHITECTURES "${SOLTRACE_OPTIX_COMPUTE_CAPABILITY}"
)

target_include_directories(OptixCSP_core
   PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}
           ${CMAKE_CURRENT_SOURCE_DIR}/core
           ${OptiX_INCLUDE}
           ${CUDAToolkit_INCLUDE_DIRS}
   PRIVATE ${CMAKE_CURRENT_BINARY_DIR})

target_link_libraries(OptixCSP_core
    PRIVATE CUDA::cuda_driver
            CUDA::cudart)

if (NOT APPLE AND ENABLE_COVERAGE)
    SET(CMAKE_CXX_FLAGS "-O0 -coverage")
endif()

# have to specify the Optix_INCLUDE directory for the CUDA compiler 
target_compile_options(OptixCSP_core PRIVATE
    $<$<COMPILE_LANGUAGE:CUDA>:--use_fast_math -lineinfo -I"${OptiX_INCLUDE}">
)
# ---------------------------------------------------------------------------
# shaders target
# ---------------------------------------------------------------------------

set(SHADERS
    intersection
    materials
    sun
)

# Use a reliable PTX output directory
set(PTX_OUTPUT_DIR "${CMAKE_BINARY_DIR}/lib/ptx")
file(MAKE_DIRECTORY "${PTX_OUTPUT_DIR}")

message(STATUS "PTX output directory: ${PTX_OUTPUT_DIR}")

foreach(shader ${SHADERS})
    set(SHADER_SRC ${CMAKE_CURRENT_SOURCE_DIR}/shaders/${shader}.cu)
    set(SHADER_PTX ${PTX_OUTPUT_DIR}/${shader}.ptx)

    add_custom_command(
        OUTPUT ${SHADER_PTX}
        COMMAND ${CMAKE_CUDA_COMPILER}
            --use_fast_math
            -lineinfo
            -std=c++17
            -arch=compute_${SOLTRACE_OPTIX_COMPUTE_CAPABILITY}
            -I"${OptiX_INCLUDE}"
            -I${CMAKE_CURRENT_SOURCE_DIR}/core
            -DOPTIX_INCLUDE_INTERNAL_FUNCTIONS
            -ptx
            -o ${SHADER_PTX}
            ${SHADER_SRC}
        DEPENDS ${SHADER_SRC}
        COMMENT
            "Generating PTX for compute_${SOLTRACE_OPTIX_COMPUTE_CAPABILITY}: ${shader}.ptx"
        VERBATIM
    )

    list(APPEND PTX_OUTPUTS ${SHADER_PTX})
    list(APPEND SHADER_FILES ${SHADER_SRC})

endforeach()

add_custom_target(OptixCSP_shaders DEPENDS ${PTX_OUTPUTS} ${SHADER_FILES}
                                   SOURCES ${SHADER_FILES}
)

install(FILES ${PTX_OUTPUTS}
    DESTINATION ${CMAKE_INSTALL_BINDIR}/ptx
)

# ---------------------------------------------------------------------------
# core depends on shaders target
# ---------------------------------------------------------------------------
add_dependencies(OptixCSP_core OptixCSP_shaders)


# ---------------------------------------------------------------------------
# generate sampleConfig.h
# ---------------------------------------------------------------------------
configure_file(
    ../sampleConfig.h.in               # contains:  #define SAMPLES_PTX_DIR "@PTX_OUTPUT_DIR@/"
    ${CMAKE_CURRENT_BINARY_DIR}/sampleConfig.h
    @ONLY)


# ---------------------------------------------------------------------------
# export targets for linking to third-party projects
# ---------------------------------------------------------------------------
#set(OptixCSP_EXPORT_CONFIG_DIR "${CMAKE_BINARY_DIR}/cmake/")

# simdata is a PRIVATE dependency of OptixCSP_core; it must be included in the
# build-tree export so consumers can resolve it. Export it without a namespace
# prefix first, then append OptixCSP_core with its own namespace.
#export(TARGETS simdata
#       FILE ${OptixCSP_EXPORT_CONFIG_DIR}/OptixCSPTargets.cmake
#)
#export(TARGETS OptixCSP_core
#       NAMESPACE OptixCSP::
#       APPEND FILE ${OptixCSP_EXPORT_CONFIG_DIR}/OptixCSPTargets.cmake
#)

# Generate Config file
#include(CMakePackageConfigHelpers)

#configure_package_config_file(
#    ${CMAKE_CURRENT_SOURCE_DIR}/../cmake/OptixCSPConfig.cmake.in
#    ${OptixCSP_EXPORT_CONFIG_DIR}/OptixCSPConfig.cmake
#    INSTALL_DESTINATION lib/cmake/OptixCSP
#)

# ---------------------------------------------------------------------------
# standalone package export
# ---------------------------------------------------------------------------
#
# OptixCSP is built as part of the parent SolTrace project. Do not export
# `simdata` from here:
#
# - `simdata` belongs to the parent SolTrace package
# - `simdata` has PUBLIC dependencies on `nlohmann_json::nlohmann_json`
#   and `glm::glm`
# - exporting it here causes CMake generate-time failures because those
#   dependencies are not part of the OptixCSP export set
#
# The parent project already owns installation/export of SolTrace targets,
# so skip nested OptixCSP package export in this build.

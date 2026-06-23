Project(optix_runner)

add_subdirectory(OptixCSP)

# OptixCSP_core's INTERFACE_INCLUDE_DIRECTORIES were set with plain absolute paths
# inside the submodule. Wrap them all in BUILD_INTERFACE so they don't leak into
# the installed targets file, and add an INSTALL_INTERFACE path for consumers.
get_target_property(_optixcsp_raw_includes OptixCSP_core INTERFACE_INCLUDE_DIRECTORIES)
set(_optixcsp_wrapped_includes)
foreach(_inc IN LISTS _optixcsp_raw_includes)
    list(APPEND _optixcsp_wrapped_includes "$<BUILD_INTERFACE:${_inc}>")
endforeach()
list(APPEND _optixcsp_wrapped_includes
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SolTrace/optix>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SolTrace/optix/core>"
)
set_property(TARGET OptixCSP_core PROPERTY INTERFACE_INCLUDE_DIRECTORIES
    ${_optixcsp_wrapped_includes}
)

set(OPTIX_RUNNER_SRC
    optix_runner.cpp
)

set(POSITION_INDEPENDENT_CODE ON)

add_library(optix_runner STATIC ${OPTIX_RUNNER_SRC})

set_target_properties(
    optix_runner
    PROPERTIES
    DEBUG_POSTFIX "d"
	PREFIX ""
)

# Publish include directories so dependents inherit them.
# optix_runner.hpp includes core/soltrace_system.h -> <optix.h>, so OptiX and
# CUDA include paths must be PUBLIC so any consumer that includes optix_runner.hpp
# can find those system headers.
target_include_directories(optix_runner
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../..>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../simulation_data>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../simulation_runner>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../simulation_results>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../../../app/deploy/samples>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/OptixCSP/src>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/OptixCSP/src/core>
        $<BUILD_INTERFACE:${OptiX_ROOT_DIR}/include>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SolTrace/simulation_runner/optix_runner>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SolTrace/simulation_runner>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SolTrace/optix>
)

target_link_libraries(optix_runner
    PRIVATE OptixCSP::OptixCSP_core
    PUBLIC  CUDA::cuda_driver
            simdata
)

install(TARGETS optix_runner OptixCSP_core
    EXPORT SolTraceTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/optix_runner.hpp"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SolTrace/simulation_runner/optix_runner
)

install(DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/OptixCSP/src/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SolTrace/optix
    FILES_MATCHING
        PATTERN "*.hpp"
        PATTERN "*.h"
        PATTERN "CMakeLists.txt" EXCLUDE
)

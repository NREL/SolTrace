Project(simdata)

include_directories(.)
include_directories(cst_templates)
include_directories(solar_position_calculators)

CPMAddPackage(
  GITHUB_REPOSITORY "nlohmann/json"
  GIT_TAG "v3.11.3"
  GIT_SHALLOW TRUE
  OPTIONS "JSON_BuildTests OFF"
)

CPMAddPackage(
  GITHUB_REPOSITORY "g-truc/glm"
  GIT_TAG "1.0.3"
  GIT_SHALLOW TRUE
  OPTIONS
    "GLM_ENABLE_CXX_17 ON"
)

set(SIMDATA_SRC
    aperture.cpp
    composite_element.cpp
    element.cpp
    json_schema.cpp
    matvec.cpp
    optical_properties.cpp
    simdata_io.cpp
    simulation_data.cpp
    single_element.cpp
    stage_element.cpp
    sun.cpp
    surface.cpp
    utilities.cpp
    virtual_element.cpp
    cst_templates/arclength.cpp
    cst_templates/heliostat.cpp
    cst_templates/linear_fresnel.cpp
    cst_templates/parabolic_dish.cpp
    cst_templates/parabolic_trough.cpp
    solar_position_calculators/solar_position_calculator.cpp
    solar_position_calculators/lib_irradproc.cpp
    solar_position_calculators/solpos.cpp
    solar_position_calculators/basic_sun_position.cpp
)

# Only headers in simulation_data (including cst_templates subfolder)
set(SIMDATA_HDRS
    aperture.hpp
    composite_element.hpp
    constants.hpp
    container.hpp
    datetime.hpp
    element.hpp
    error_distributions.hpp
    json_helpers.hpp
    json_schema.hpp
    simulation_data_export.hpp
    matvec.hpp
    optical_properties.hpp
    ray_source.hpp
    simdata_io.hpp
    simulation_data.hpp
    simulation_parameters.hpp
    single_element.hpp
    stage_element.hpp
    sun.hpp
    surface.hpp
    vector_utility.hpp
    utilities.hpp
    virtual_element.hpp
    cst_templates/arclength.hpp
    cst_templates/heliostat.hpp
    cst_templates/linear_fresnel.hpp
    cst_templates/parabolic_dish.hpp
    cst_templates/parabolic_trough.hpp
    solar_position_calculators/solar_position_calculator.hpp
    solar_position_calculators/lib_irradproc.h
    solar_position_calculators/solpos00.h
    solar_position_calculators/basic_sun_position.hpp
)

source_group("Header Files" FILES ${SIMDATA_HDRS})

set(POSITION_INDEPENDENT_CODE ON)

add_library(simdata STATIC ${SIMDATA_SRC} ${SIMDATA_HDRS})
set_target_properties(simdata
    PROPERTIES
    DEBUG_POSTFIX "d"
    PREFIX ""
)

target_include_directories(simdata
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cst_templates>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/solar_position_calculators>
        $<BUILD_INTERFACE:${nlohmann_json_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${glm_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SolTrace/simulation_data>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SolTrace/simulation_data/cst_templates>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SolTrace/simulation_data/solar_position_calculators>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# Link dependencies via BUILD_INTERFACE: nlohmann_json and glm are header-only;
# their include directories are propagated above via target_include_directories,
# and consumers find them through find_dependency() in SolTraceConfig.cmake.
# Wrapping in $<BUILD_INTERFACE:...> excludes these CPM-fetched (non-IMPORTED)
# targets from the install interface, avoiding the CMake export error:
#   install(EXPORT "SolTraceTargets") includes target "simdata" which requires
#   target "nlohmann_json"/"glm" that is not in any export set.
target_link_libraries(simdata
    PRIVATE
    $<BUILD_INTERFACE:nlohmann_json::nlohmann_json>
    $<BUILD_INTERFACE:glm::glm>
)

if (NOT APPLE AND ENABLE_COVERAGE)
    # SET(CMAKE_CXX_FLAGS "-O0 -coverage -fkeep-inline-functions")
    SET(CMAKE_CXX_FLAGS "-O0 -coverage")
endif()

install(TARGETS simdata
    EXPORT SolTraceTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SolTrace/simulation_data
    FILES_MATCHING
        PATTERN "*.hpp"
        PATTERN "*.h"
        PATTERN "CMakeLists.txt" EXCLUDE
)

# Install nlohmann_json headers for consumers of simdata
install(DIRECTORY "${nlohmann_json_SOURCE_DIR}/include/"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install glm headers for consumers of simdata
install(DIRECTORY "${glm_SOURCE_DIR}/glm"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

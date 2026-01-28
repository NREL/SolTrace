Project(simdata)

include_directories(.)
include_directories(cst_templates)
include_directories(solar_position_calculators)

include(FetchContent)
    FetchContent_Declare(
        nlohmann_json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG v3.11.3
        GIT_SHALLOW TRUE
    )
    FetchContent_MakeAvailable(nlohmann_json)

set(SIMDATA_SRC
    aperture.cpp
    composite_element.cpp
    element.cpp
    matvec.cpp
    optical_properties.cpp
    simdata_io.cpp
    simulation_data.cpp
    single_element.cpp
    stage_element.cpp
    sun.cpp
    surface.cpp
    utilities.cpp
    vector3d.cpp
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
    utilities.hpp
    vector3d.hpp
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

if(WIN32)
    # Hack to build on Windows
    add_library(simdata STATIC ${SIMDATA_SRC} ${SIMDATA_HDRS})
else(WIN32)
    add_library(simdata SHARED ${SIMDATA_SRC} ${SIMDATA_HDRS})
endif()
set_target_properties(simdata
    PROPERTIES
    DEBUG_POSTFIX "d"
    PREFIX ""
)

# Link JSON (propagates its include directory)
target_link_libraries(simdata PUBLIC nlohmann_json::nlohmann_json)

if (NOT APPLE AND ENABLE_COVERAGE)
    # SET(CMAKE_CXX_FLAGS "-O0 -coverage -fkeep-inline-functions")
    SET(CMAKE_CXX_FLAGS "-O0 -coverage")
endif()

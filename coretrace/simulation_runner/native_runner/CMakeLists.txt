Project(native_runner)

include_directories(
    .
    ../../simulation_data
    ../../simulation_results
    ../../simulation_runner
)

set(NATIVE_RUNNER_SRC
    calculator_factory.cpp
    cylinder_calculator.cpp
    determine_element_intersection_new.cpp
    determine_interaction_type.cpp
    find_element_hit.cpp
    flat_calculator.cpp
    generate_ray.cpp
    native_runner.cpp
    native_runner_types.cpp
    newton_calculator.cpp
    parabola_calculator.cpp
    process_interaction.cpp
    pt_optimizations.cpp
    sphere_calculator.cpp
    sun_to_primary_stage.cpp
    thread_manager.cpp
    trace.cpp
    tracing_errors.cpp
    trace_logger.cpp
    treemesh.cpp
)

set(POSITION_INDEPENDENT_CODE ON)

if(WIN32)
    # Hack to build on Windows
    add_library(native_runner STATIC ${NATIVE_RUNNER_SRC})
else(WIN32)
    add_library(native_runner SHARED ${NATIVE_RUNNER_SRC})
endif()

set_target_properties(
    native_runner
    PROPERTIES
    DEBUG_POSTFIX "d"
	PREFIX ""
)

target_link_libraries(
    native_runner
    simdata
    simresult
)

if (NOT APPLE AND ENABLE_COVERAGE)
    SET(CMAKE_CXX_FLAGS "-O0 -coverage")
endif()

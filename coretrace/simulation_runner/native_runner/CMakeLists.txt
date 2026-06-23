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

add_library(native_runner STATIC ${NATIVE_RUNNER_SRC})

set_target_properties(
    native_runner
    PROPERTIES
    DEBUG_POSTFIX "d"
	PREFIX ""
)

target_include_directories(native_runner
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/..>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SolTrace/simulation_runner/native_runner>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SolTrace/simulation_runner>
)

target_link_libraries(
    native_runner
    simdata
    simresult
)

if (NOT APPLE AND ENABLE_COVERAGE)
    SET(CMAKE_CXX_FLAGS "-O0 -coverage")
endif()

install(TARGETS native_runner
    EXPORT SolTraceTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/../simulation_runner.hpp"
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SolTrace/simulation_runner
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SolTrace/simulation_runner/native_runner
    FILES_MATCHING
        PATTERN "*.hpp"
        PATTERN "*.h"
        PATTERN "CMakeLists.txt" EXCLUDE
)

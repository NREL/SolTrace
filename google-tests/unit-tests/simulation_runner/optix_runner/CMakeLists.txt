Project(OptixRunnerUnitTests)

include_directories(
    .
    ../../common
    ${SAMPLE_DIR}
    ${TEST_INPUT_DIR}
    ../../../../coretrace/simulation_data
    ../../../../coretrace/simulation_results
    ../../../../coretrace/simulation_runner
    ../../../../coretrace/simulation_runner/native_runner
    ../../../../coretrace/simulation_runner/optix_runner
)

# Add test files for optix_runner
set(OPTIX_RUNNER_TEST_SRC
    ../../common/common.cpp
    vector_error_test.cpp
    geometry_intersection_test.cpp
    gpu_tower_demo.cpp
    flat_optical_test.cpp
    two_plate_test.cpp
    sun_test.cpp
    ray_position_sampling_test.cpp
    batch_size_test.cpp
    max_ray_depth_test.cpp
    high_flux_solar_furnace_test.cpp
)

add_executable(OptixRunnerUnitTests
    ${OPTIX_RUNNER_TEST_SRC}
)

target_link_libraries(
    OptixRunnerUnitTests
    simdata
    simresult
    native_runner
    optix_runner
    OptixCSP_core
    GTest::gtest_main
)

if (ENABLE_COVERAGE AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_link_options(OptixRunnerUnitTests PRIVATE --coverage)
    target_link_libraries(OptixRunnerUnitTests gcov)
endif()

if(UNIX)
    target_link_libraries(OptixRunnerUnitTests pthread)
endif()

add_test(NAME OptixRunnerUnitTests COMMAND OptixRunnerUnitTests)

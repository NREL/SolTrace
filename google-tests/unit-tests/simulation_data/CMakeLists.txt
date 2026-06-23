Project(SimulationDataUnitTests)

include_directories(
    .
    ../common
    ${SAMPLE_DIR}
    ${TEST_INPUT_DIR}
    ../../../coretrace/simulation_data
    ../../../coretrace/simulation_results
    ../../../coretrace/simulation_runner
    ../../../coretrace/simulation_runner/native_runner
    ../../../coretrace/simulation_data/solar_position_calculators
)

# Add CST template tests subdirectory
add_subdirectory(cst-templates)

# Add test files for simulation_data
set(SIMULATION_DATA_TEST_SRC
    ../common/common.cpp
    aperture_test.cpp
    surface_test.cpp
    container_test.cpp
    element_test.cpp
    file_io_test.cpp
    linear_algebra_test.cpp
    optical_properties_test.cpp
    simulation_data_test.cpp
    sun_position_test.cpp
    sun_test.cpp
    utilities_test.cpp
)

add_executable(SimulationDataUnitTests
    ${SIMULATION_DATA_TEST_SRC}
)

target_link_libraries(
    SimulationDataUnitTests
    simdata
    simresult
    native_runner
    testtools
    GTest::gtest_main
)

if (ENABLE_COVERAGE AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_link_options(SimulationDataUnitTests PRIVATE --coverage)
    target_link_libraries(SimulationDataUnitTests gcov)
endif()

if(UNIX)
    target_link_libraries(SimulationDataUnitTests pthread)
endif()

add_test(NAME SimulationDataUnitTests COMMAND SimulationDataUnitTests)

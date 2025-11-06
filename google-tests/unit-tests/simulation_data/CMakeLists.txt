Project(SimulationDataUnitTests)

include_directories(
    .
    ../common
    ../../app/deploy/samples
    ../../../coretrace/simulation_data
    ../../../coretrace/simulation_data/solar_position_calculators
)

# Add CST template tests subdirectory
add_subdirectory(cst-templates)

# Add test files for simulation_data
set(SIMULATION_DATA_TEST_SRC
    ../common/common.cpp
    linear_algebra_test.cpp
    aperture_test.cpp
    surface_test.cpp
    container_test.cpp
    element_test.cpp
    simulation_data_test.cpp
    sun_position_test.cpp
    sun_test.cpp
)

add_executable(SimulationDataUnitTests
    ${SIMULATION_DATA_TEST_SRC}
)

target_link_libraries(
    SimulationDataUnitTests
    simdata
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

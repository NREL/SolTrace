Project(SimulationResultsUnitTests)

include_directories(
    .
    ../common
    ${SAMPLE_DIR}
    ${TEST_INPUT_DIR}
    ../../../coretrace/simulation_results
    ../../../coretrace/simulation_data
)

# Add test files for simulation_results
set(SIMULATION_RESULTS_TEST_SRC
    ../common/common.cpp
    simulation_result_test.cpp
)

add_executable(SimulationResultsUnitTests
    ${SIMULATION_RESULTS_TEST_SRC}
)

target_link_libraries(
    SimulationResultsUnitTests
    simdata
    simresult
    GTest::gtest_main
)

if (ENABLE_COVERAGE AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_link_options(SimulationResultsUnitTests PRIVATE --coverage)
    target_link_libraries(SimulationResultsUnitTests gcov)
endif()

if(UNIX)
    target_link_libraries(SimulationResultsUnitTests pthread)
endif()

add_test(NAME SimulationResultsUnitTests COMMAND SimulationResultsUnitTests)

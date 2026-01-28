Project(NativeRunnerUnitTests)

include_directories(
    .
    ../../common
    ../../../test_tools
    ../../../app/deploy/samples
    ../../../../coretrace/simulation_data
    ../../../../coretrace/simulation_results
    ../../../../coretrace/simulation_runner
    ../../../../coretrace/simulation_runner/native_runner
)

# Add test files for native_runner
set(NATIVE_RUNNER_TEST_SRC
    ../../common/common.cpp
    flat_calculator_test.cpp
    newton_calculator_test.cpp
    parabola_calculator_test.cpp
    sphere_calculator_test.cpp
    cylinder_calculator_test.cpp
    calculator_factory_test.cpp
    thread_manager_test.cpp
    native_runner_test.cpp
    tower_demo.cpp
)

if (NOT ENABLE_COVERAGE)
    list(APPEND
    NATIVE_RUNNER_TEST_SRC
    native_runner_multithreading_test.cpp
    )
endif()

add_executable(NativeRunnerUnitTests
    ${NATIVE_RUNNER_TEST_SRC}
)

target_link_libraries(
    NativeRunnerUnitTests
    simdata
    simresult
    native_runner
    testtools
    GTest::gtest_main
)

if (ENABLE_COVERAGE AND CMAKE_BUILD_TYPE STREQUAL "Debug")
    target_link_options(NativeRunnerUnitTests PRIVATE --coverage)
    target_link_libraries(NativeRunnerUnitTests gcov)
endif()

if(UNIX)
    target_link_libraries(NativeRunnerUnitTests pthread)
endif()

add_test(NAME NativeRunnerUnitTests COMMAND NativeRunnerUnitTests)

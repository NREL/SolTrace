
Project(RegressionTests)

include_directories(
	.
	# ../app/deploy/samples
	../test_tools
	../../coretrace/simulation_data
	../../coretrace/simulation_data/cst_templates
	../../coretrace/simulation_results
	../../coretrace/simulation_runner
	../../coretrace/simulation_runner/native_runner
)

# file(GLOB UNIT_TEST_SRC *.cpp)

# Add test files below as needed. Add more basic ones earlier in the list
# so those tests run prior to the more complicated ones.
# common.hpp and common.cpp contain helper methods used in multiple files.
set(REGRESSION_TEST_SRC
	native_runner_validation_test.cpp
)

if(CMAKE_BUILD_TYPE STREQUAL "Release" OR SOLTRACE_BUILD_PERF_TEST)
	list(APPEND 
	REGRESSION_TEST_SRC 
	native_runner_multithreading_performance_test.cpp
	native_runner_performance_test.cpp
	)
endif()

add_executable(RegressionTests
	${CSV_TOOL}
    ${REGRESSION_TEST_SRC}
)

target_link_libraries(
    RegressionTests
    simdata
	simresult
	native_runner
	testtools
	GTest::gtest_main
)

if (ENABLE_COVERAGE AND CMAKE_BUILD_TYPE STREQUAL "Debug")
	target_link_options(RegressionTests PRIVATE --coverage)
	target_link_libraries(RegressionTests gcov)
# else()
# 	target_link_libraries(RegressionTests coretrace GTest::gtest_main)
endif()

if(UNIX)
	target_link_libraries(RegressionTests pthread)
endif()

add_test(NAME RegressionTests COMMAND RegressionTests)

Project(CSTUnitTests)

include_directories(
	.
	../
	../../common
	../../../test_tools
	${SAMPLE_DIR}
	${TEST_INPUT_DIR}
	../../../../coretrace/simulation_data
	../../../../coretrace/simulation_results
	../../../../coretrace/simulation_runner
	../../../../coretrace/simulation_runner/native_runner
	../../../../coretrace/simulation_runner/optix_runner
	../../../../coretrace/simulation_runner/embree_runner
)

# CST Template test files
set(CST_TEST_SRC
	# utilities_test.cpp
	linear_fresnel_test.cpp
	parabolic_trough_test.cpp
	parabolic_dish_test.cpp
	heliostat_test.cpp
)

add_executable(CSTUnitTests
	${CST_TEST_SRC}
	${CST_TEST_HEADERS}
)

# Always link core simdata and native runner
set(CST_LIBS
	simdata
	native_runner
	testtools
	GTest::gtest_main
)

# Conditionally link OptiX and Embree runners
if(SOLTRACE_BUILD_OPTIX_SUPPORT)
	list(APPEND CST_LIBS optix_runner)
endif()

if(SOLTRACE_BUILD_EMBREE_SUPPORT)
	list(APPEND CST_LIBS embree_runner)
endif()

target_link_libraries(CSTUnitTests ${CST_LIBS})

if (ENABLE_COVERAGE AND CMAKE_BUILD_TYPE STREQUAL "Debug")
	target_link_options(CSTUnitTests PRIVATE --coverage)
	target_link_libraries(CSTUnitTests gcov)
endif()

if(UNIX)
	target_link_libraries(CSTUnitTests pthread)
endif()

add_test(NAME CSTUnitTests COMMAND CSTUnitTests)

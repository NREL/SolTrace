Project(CSTUnitTests)

include_directories(
	.
	../
	../../common
	../../../test_tools
	../../../app/deploy/samples
	../../../../coretrace/simulation_data
	../../../../coretrace/simulation_results
	../../../../coretrace/simulation_runner
	../../../../coretrace/simulation_runner/native_runner
)

# CST Template test files
set(CST_TEST_SRC
	utilities_test.cpp
	linear_fresnel_test.cpp
	parabolic_trough_test.cpp
	parabolic_dish_test.cpp
	heliostat_test.cpp
)

if(CMAKE_BUILD_TYPE STREQUAL "Release" OR FULL_FIELD_VALIDATION_TEST)
	list(APPEND 
	CST_TEST_SRC
	single_heliostat_test.cpp
	heliostat_field_test.cpp
	)
endif()

add_executable(CSTUnitTests
	${CST_TEST_SRC}
)

target_link_libraries(
	CSTUnitTests
	simdata
	native_runner
	testtools
	GTest::gtest_main
)

if (ENABLE_COVERAGE AND CMAKE_BUILD_TYPE STREQUAL "Debug")
	target_link_options(CSTUnitTests PRIVATE --coverage)
	target_link_libraries(CSTUnitTests gcov)
endif()

if(UNIX)
	target_link_libraries(CSTUnitTests pthread)
endif()

add_test(NAME CSTUnitTests COMMAND CSTUnitTests)


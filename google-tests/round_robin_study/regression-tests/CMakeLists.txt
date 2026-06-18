Project(CSTRegressionTests)

include_directories(
	.
	../../unit-tests/common
	../../test_tools
	../../app/deploy/samples
	../../../coretrace/simulation_data
	../../../coretrace/simulation_results
	../../../coretrace/simulation_runner
	../../../coretrace/simulation_runner/native_runner
	../../../coretrace/simulation_runner/optix_runner
	../../../coretrace/simulation_runner/embree_runner
)

# Add headers for Visual Studio solution
set(CST_TEST_HEADERS
	heliostat_field_test_template.hpp
	single_heliostat_test_template.hpp
)

if(CMAKE_BUILD_TYPE STREQUAL "Release" OR FULL_FIELD_VALIDATION_TEST)
	# In release or full-field validation builds, also add backend-specific tests
	if(SOLTRACE_BUILD_OPTIX_SUPPORT)
		list(APPEND CST_TEST_SRC
			single_heliostat_test_optix.cpp
			heliostat_field_test_optix.cpp
		)
	endif()
	if(SOLTRACE_BUILD_EMBREE_SUPPORT)
		list(APPEND CST_TEST_SRC
			single_heliostat_test_embree.cpp
			heliostat_field_test_embree.cpp
		)
	endif()
	if(FULL_FIELD_VALIDATION_TEST)
		list(APPEND CST_TEST_SRC
			single_heliostat_test_native.cpp
			heliostat_field_test_native.cpp
		)
	endif()
endif()

if(CST_TEST_SRC)
	add_executable(CSTRegressionTests
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

	target_link_libraries(CSTRegressionTests PRIVATE ${CST_LIBS})

	if(ENABLE_COVERAGE AND CMAKE_BUILD_TYPE STREQUAL "Debug")
		target_link_options(CSTRegressionTests PRIVATE --coverage)
		target_link_libraries(CSTRegressionTests PRIVATE gcov)
	endif()

	if(UNIX)
		target_link_libraries(CSTRegressionTests PRIVATE pthread)
	endif()

	add_test(NAME CSTRegressionTests COMMAND CSTRegressionTests)
else()
	message(STATUS "Skipping CSTRegressionTests because no test source files were selected.")
endif()


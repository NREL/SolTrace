
include(FetchContent)
FetchContent_Declare(
	googletest
	GIT_REPOSITORY https://github.com/google/googletest.git
	GIT_TAG main
)

# Common for all tests
add_definitions(-DPROJECT_DIR="${CMAKE_CURRENT_BINARY_DIR}")
# file(GLOB CSV_TOOL test_tools/*.cpp test_tools/*.h)

add_subdirectory(test_tools)
add_subdirectory(unit-tests)
add_subdirectory(regression-tests)

set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

file(GLOB SAMPLE_FILES ../app/deploy/samples/*.stinput ../app/deploy/samples/*.mon)
file(GLOB GROUND_FILES ./ground_results/*.csv)

file(COPY ${SAMPLE_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${GROUND_FILES} DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/round_robin_study" DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if(SOLTRACE_BUILD_CORETRACE)

	include_directories(
		soltrace_example_tests
		test_tools
		../coretrace
		../coretrace/strace
		../app/deploy/samples
		../coretrace/simulation_data
	)

	file(GLOB SOLTRACE_EXAMPLE_TESTS soltrace_example_tests/*.cpp)

	set(STRACE_SRC ../coretrace/strace)
	set(CORETRACE ../coretrace/stapi.cpp)
	set(STRACE ../coretrace/strace/strace.cpp)
	set(ENV{SOLTRACE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}")

	message("Binary dir: ${CMAKE_CURRENT_BINARY_DIR}")
	message("Source dir: ${CMAKE_CURRENT_SOURCE_DIR}")
	message("Build type: ${CMAKE_BUILD_TYPE}")

	add_executable( Test
		# ${CSV_TOOL}
		${SOLTRACE_EXAMPLE_TESTS}
		${CORETRACE}
		${STRACE}
	)

	target_link_libraries(
		Test
		coretrace
		testtools
		GTest::gtest_main
	)

	if (ENABLE_COVERAGE AND CMAKE_BUILD_TYPE STREQUAL "Debug")
		target_link_options(Test PRIVATE --coverage)
		target_link_libraries(Test gcov)
	endif()

	if(UNIX)
		target_link_libraries(Test pthread)
	endif()

endif()

enable_testing()

if(SOLTRACE_BUILD_CORETRACE)
	add_test(NAME Test COMMAND Test)
endif()

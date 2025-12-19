cmake_minimum_required(VERSION 3.19)

Project(soltrace_ui VERSION 1.0.0)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(SOLTRACE_BUILD_GUI "Build the GUI parts of Soltrace" OFF)
option(SOLTRACE_BUILD_CORETRACE "Build the coretrace library (ray tracing core)" ON)
option(SOLTRACE_BUILD_OPTIX_SUPPORT "Build the OptiX support for ray tracing" OFF)
option(SOLTRACE_BUILD_PERF_TEST "Build performance tests (not built for debug)" OFF)

option(FULL_FIELD_VALIDATION_TEST "Include full field validation tests" OFF)


enable_testing()

add_subdirectory(coretrace)
add_subdirectory(google-tests)

if(SOLTRACE_BUILD_GUI)
	add_subdirectory(app)
endif()


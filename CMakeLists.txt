cmake_minimum_required(VERSION 3.19)

Project(soltrace_ui VERSION 4.0.0)

include(cmake/CPM.cmake)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(MSVC)
	add_compile_options("$<$<AND:$<CONFIG:Debug>,$<COMPILE_LANG_AND_ID:CXX,MSVC>>:/bigobj>")
endif()

option(SOLTRACE_BUILD_GUI "Build the GUI parts of Soltrace" ON)
option(
	SOLTRACE_BUILD_CORETRACE 
	"Legacy ray tracing backend including the strace executable"
	OFF)

option(SOLTRACE_BUILD_EMBREE_SUPPORT "Build Embree support for ray tracing" OFF)
option(SOLTRACE_BUILD_OPTIX_SUPPORT "Build the OptiX support for ray tracing" OFF)

option(SOLTRACE_BUILD_PERF_TEST "Build performance tests (not built for debug)" OFF)
option(SOLTRACE_BUILD_TESTS "Build SolTrace test targets" ON)
option(SOLTRACE_ENABLE_ASAN "Enable AddressSanitizer for supported builds" OFF)

option(FULL_FIELD_VALIDATION_TEST "Include full field validation tests" OFF)

if(SOLTRACE_ENABLE_ASAN)
	if(MSVC)
		message(WARNING "SOLTRACE_ENABLE_ASAN is not yet supported with MSVC or GCC")
	elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang|AppleClang")
		message(STATUS "Enabling address sanitizer")
		add_compile_options(
			$<$<COMPILE_LANGUAGE:C,CXX>:-fsanitize=address>
			$<$<COMPILE_LANGUAGE:C,CXX>:-fno-omit-frame-pointer>
		)
		add_link_options(-fsanitize=address)
	else()
		message(WARNING "SOLTRACE_ENABLE_ASAN is only supported with Clang")
	endif()
endif()


if(SOLTRACE_BUILD_TESTS)
	enable_testing()
endif()

include(GNUInstallDirs)
include(CMakePackageConfigHelpers)

configure_package_config_file(
	"${CMAKE_CURRENT_SOURCE_DIR}/cmake/SolTraceConfig.cmake.in"
	"${CMAKE_CURRENT_BINARY_DIR}/SolTraceConfig.cmake"
	INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/SolTrace"
)

write_basic_package_version_file(
	"${CMAKE_CURRENT_BINARY_DIR}/SolTraceConfigVersion.cmake"
	VERSION ${PROJECT_VERSION}
	COMPATIBILITY SameMajorVersion
)

install(FILES
	"${CMAKE_CURRENT_BINARY_DIR}/SolTraceConfig.cmake"
	"${CMAKE_CURRENT_BINARY_DIR}/SolTraceConfigVersion.cmake"
	"${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindOptiX.cmake"
	DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/SolTrace"
)

install(EXPORT SolTraceTargets
	FILE SolTraceTargets.cmake
	NAMESPACE SolTrace::
	DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/SolTrace"
)

add_subdirectory(coretrace)

if(SOLTRACE_BUILD_TESTS)
	add_subdirectory(google-tests)
endif()

if(SOLTRACE_BUILD_GUI)
	add_subdirectory(gui)
endif()

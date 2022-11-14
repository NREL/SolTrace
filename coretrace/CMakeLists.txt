#####################################################################################################################
#
# CMake and System Settings
#
#####################################################################################################################

set(CMAKE_VERBOSE_MAKEFILE ON)

if(APPLE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9" CACHE STRING "Minimum OS X deployment version")
endif()


#####################################################################################################################
#
# Project Settings 
#
#####################################################################################################################

Project(coretrace)

include_directories(.)

set(CORETRACE_SRC
	interpolate.cpp
	apertureplane.cpp
	determineelementintersectionnew.cpp
	dumpsys.cpp
	errors.cpp
	generateray.cpp
	input.cpp
	interaction.cpp
	intersect.cpp
	mathproc.cpp
	newzstartforcubicsplinesurf.cpp
	quadricsurfaceclosedform.cpp
	raytrace.cpp
	root432.cpp
	spencerandmurtysurfaceclosedform.cpp
	stapi.cpp
	suntoprimarystage.cpp
	surface.cpp
	surfacenormalerrors.cpp
	surfacezatxypair.cpp
	torusclosedform.cpp
	treemesh.cpp
	types.cpp
	vshot.cpp)


#####################################################################################################################
#
# Compile Options per Platform
#
#####################################################################################################################

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if(MSVC)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
	add_compile_options(/MP /W0)
    add_compile_definitions( __64BIT__ __WINDOWS__ _WINDOWS _CRT_SECURE_NO_WARNINGS WIN32 _MBCS )
    foreach (flag_var CMAKE_C_FLAGS_DEBUG CMAKE_CXX_FLAGS_DEBUG)
		set(${flag_var} "${${flag_var}} /D_DEBUG" CACHE STRING "compile flags" FORCE)
	endforeach()
else(MSVC)
	if (APPLE)
		add_compile_options(-arch x86_64 -fno-common)
  		add_definitions(-DWX_PRECOMP)
	else()
		add_definitions(-D__64BIT__)
	endif()
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unknown-pragmas -Wno-sign-compare -Wno-unused-variable -Wno-uninitialized -Wno-int-to-pointer-cast ")
	if(CMAKE_BUILD_TYPE STREQUAL "Debug")
		add_compile_definitions(_DEBUG)
	else()
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -DNDEBUG" )
	endif()
endif(MSVC)


#####################################################################################################################
# 
# CMake Targets
#
#####################################################################################################################

add_library(coretrace_api SHARED ${CORETRACE_SRC})
set_target_properties(coretrace_api
    PROPERTIES
    DEBUG_POSTFIX "d"
	PREFIX ""
	)
if(MSVC)
	target_compile_definitions(coretrace_api PRIVATE STCORE_API_EXPORTS _STCOREDLL_)
endif()
add_custom_command(
	TARGET coretrace_api
	POST_BUILD
	COMMAND ${CMAKE_COMMAND}
	ARGS -E copy $<TARGET_FILE:coretrace_api> ${CMAKE_CURRENT_SOURCE_DIR}/../app/deploy/api
	)
	
add_library(coretrace STATIC ${CORETRACE_SRC})
set_target_properties(coretrace
    PROPERTIES
    DEBUG_POSTFIX "d"
	PREFIX ""
	)


#####################################################################################################################
#
# Link Libraries and Options
#
#####################################################################################################################

if (MSVC)
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /SUBSYSTEM:WINDOWS")
endif()

add_subdirectory(strace)
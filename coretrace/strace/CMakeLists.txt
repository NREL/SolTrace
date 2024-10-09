set(CORETRACEDIR ${CMAKE_CURRENT_SOURCE_DIR}/../coretrace)


#####################################################################################################################
#
# CMake and System Settings
#
#####################################################################################################################
cmake_minimum_required(VERSION 3.19)

set(CMAKE_VERBOSE_MAKEFILE ON)

if(APPLE)
    set(CMAKE_OSX_DEPLOYMENT_TARGET "10.9" CACHE STRING "Minimum OS X deployment version")
endif()

if (UNIX AND NOT CMAKE_C_COMPILER)
	set(CMAKE_C_COMPILER gcc)
	set(CMAKE_CXX_COMPILER g++)
endif()

#####################################################################################################################
#
# Project Settings 
#
#####################################################################################################################
project(strace)

include_directories(. ..)

set(STRACE_SRC
	main.cpp
)


#####################################################################################################################
#
# Compile Options per Platform
#
#####################################################################################################################

cmake_policy(SET CMP0043 NEW)

if(MSVC)
	add_compile_options(/W4 /MP /GS)
	
	#add_compile_definitions(WIN32)
	add_compile_definitions(_CRT_SECURE_NO_DEPRECATE=1)
	add_compile_definitions(_CRT_NON_CONFORMING_SWPRINTFS=1)
	add_compile_definitions(_SCL_SECURE_NO_WARNINGS=1)
	add_compile_definitions(NDEBUG)
	add_compile_definitions(_UNICODE)
	add_compile_definitions(_WINDOWS)
	add_compile_definitions(NOPCH)
	
	foreach (flag_var CMAKE_C_FLAGS_DEBUG CMAKE_CXX_FLAGS_DEBUG)
		set(${flag_var} "${${flag_var}} /D_DEBUG" CACHE STRING "compile flags" FORCE)
	endforeach()
else(MSVC)
	if (APPLE)
		set(CMAKE_POSITION_INDEPENDENT_CODE ON)
		add_compile_options(-arch x86_64)
	else()
		add_definitions(-D__64BIT__)
	endif()
	add_compile_options(-Wall -Wno-deprecated -Wno-unknown-pragmas -Wno-overloaded-virtual -fno-common -std=c++11)
	add_definitions( -DWX_PRECOMP )
	if(CMAKE_BUILD_TYPE STREQUAL "Debug")
		add_compile_definitions(_DEBUG)
	else()
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3" )
	endif()
endif(MSVC)


#####################################################################################################################
#
# Packages
#
#####################################################################################################################

#####################################################################################################################
# 
# CMake Targets
#
#####################################################################################################################

# Setup Target
add_executable(strace ${STRACE_SRC} )

set_target_properties(strace
    PROPERTIES
    DEBUG_POSTFIX "d"
	PREFIX ""
	RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_SOURCE_DIR}
	RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_SOURCE_DIR}
	)

if (APPLE)
	set_target_properties(strace PROPERTIES
		MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/strace.plist
	)
endif()

if (MSVC)
	set_target_properties(strace PROPERTIES
        LINK_FLAGS /SUBSYSTEM:CONSOLE)
endif()

#####################################################################################################################
#
# Link Libraries and Options
#
#####################################################################################################################

target_link_libraries(strace coretrace)

if (UNIX AND NOT APPLE)
	target_link_libraries(strace -ldl -lfontconfig)
endif()

if (UNIX)
	target_link_libraries(strace -lm )
elseif (MSVC)
	target_link_libraries(strace Winhttp)
endif()

#####################################################################################################################
#
# Target Installation
#
#####################################################################################################################





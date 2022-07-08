
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
set(CORETRACEDIR ${CMAKE_CURRENT_SOURCE_DIR}/../coretrace)

include_directories(. ${CORETRACEDIR} $ENV{WEXDIR}/include $ENV{LKDIR}/include )

set(SOLTRACE_SRC
	src/console.cpp
	src/elementlist.cpp
	src/fluxmap.cpp
	src/geometry.cpp
	src/intersections.cpp
	src/optics.cpp
	src/project.cpp
	src/raydata.cpp
	src/script.cpp
	src/soltrace.cpp
	src/sunshape.cpp
	src/trace.cpp
)


#####################################################################################################################
#
# Compile Options per Platform
#
#####################################################################################################################

cmake_policy(SET CMP0043 NEW)

if(MSVC)
	add_compile_options(/W4 /MP /GS)
	
	add_compile_definitions(LK_USE_WXWIDGETS)
	add_compile_definitions(WIN32)
	add_compile_definitions(_CRT_SECURE_NO_DEPRECATE=1)
	add_compile_definitions(_CRT_NON_CONFORMING_SWPRINTFS=1)
	add_compile_definitions(_SCL_SECURE_NO_WARNINGS=1)
	add_compile_definitions(__WXMSW__)
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
		add_definitions(-D__WXOSX__)
	else()
		add_definitions(-D__64BIT__)
	endif()
	add_compile_options(-Wall -Wno-deprecated -Wno-unknown-pragmas -Wno-overloaded-virtual -fno-common -std=c++11)
	add_definitions(-DLK_USE_WXWIDGETS -DWX_PRECOMP )
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

if (UNIX)
	set(wxWidgets_CONFIG_EXECUTABLE /usr/local/bin/wx-config-3)
	find_package(wxWidgets COMPONENTS qa aui richtext html propgrid adv net stc core xml base REQUIRED)
else()
	set(wxWidgets_ROOT_DIR $ENV{WXMSW3})
	find_package(wxWidgets COMPONENTS qa webview aui richtext html propgrid adv net stc core xml base scintilla REQUIRED)
endif()

include(${wxWidgets_USE_FILE})


#####################################################################################################################
# 
# CMake Targets
#
#####################################################################################################################

# Setup Application
if (APPLE)
	set(MACOSX_BUNDLE_ICON_FILE ${CMAKE_CURRENT_SOURCE_DIR}/resource/soltrace.ico)
	set(SOLTRACE_ICON ${CMAKE_CURRENT_SOURCE_DIR}/resource/soltrace.ico)
	set_source_files_properties(${SOLTRACE_ICON} PROPERTIES
		MACOSX_PACKAGE_LOCATION "Resources")
elseif(WIN32)
	set(SOLTRACE_ICON ${CMAKE_CURRENT_SOURCE_DIR}/resource/soltrace.rc)
endif()

# Setup Target
add_executable(SolTrace
		WIN32
		MACOSX_BUNDLE 
		${SOLTRACE_SRC}
		${SOLTRACE_ICON})
set_target_properties(SolTrace
    PROPERTIES
    DEBUG_POSTFIX "d"
	PREFIX ""
	RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_SOURCE_DIR}/deploy/x64
	RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_SOURCE_DIR}/deploy/x64
	)

if (APPLE)
	set_target_properties(SolTrace PROPERTIES
		MACOSX_BUNDLE_INFO_PLIST ${CMAKE_CURRENT_SOURCE_DIR}/../resource/Info-SolTrace.plist
	)
endif()

if (MSVC)
	set_target_properties(SolTrace PROPERTIES
        LINK_FLAGS /SUBSYSTEM:WINDOWS)
endif()

#####################################################################################################################
#
# Link Libraries and Options
#
#####################################################################################################################

target_link_libraries(SolTrace coretrace wex lk)
target_link_libraries(SolTrace ${wxWidgets_LIBRARIES})

if (UNIX AND NOT APPLE)
	target_link_libraries(SolTrace -ldl -lfontconfig)
endif()

if (UNIX)
	target_link_libraries(SolTrace -lm -lcurl)
elseif (MSVC)
	find_library( CURL_LIB
			NAMES libcurl.lib
			PATHS $ENV{WEXDIR}/build_resources/libcurl_ssl_x64/lib)
	target_link_libraries(SolTrace ${CURL_LIB})
	target_link_libraries(SolTrace Winhttp)
endif()

#####################################################################################################################
#
# Target Installation
#
#####################################################################################################################

if (MSVC)
	set( LIBCURL_FILES
	   dbghelp.dll
	   libcurl.dll
	   libeay32.dll
	   libssh2.dll
	   msvcp140.dll
	   openssl.exe
	   ssleay32.dll
	   vcruntime140.dll
	)

	set(SOLTRACE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/deploy/x64)

	foreach( file_i ${LIBCURL_FILES})
		add_custom_command(
		TARGET SolTrace
		POST_BUILD
		COMMAND ${CMAKE_COMMAND}
		ARGS -E copy $ENV{WEXDIR}/build_resources/libcurl_ssl_x64/bin/${file_i} ${SOLTRACE_DIR}
	)
	endforeach( file_i )
endif()



cmake_minimum_required(VERSION 3.8)

project(OptixCSP)

# When built as subdirectory, OptiX is already found by parent
if(NOT OptiX_INCLUDE)
    if(OptiX_INCLUDE_DIR)
        set(OptiX_INCLUDE ${OptiX_INCLUDE_DIR})
    else()
        find_package(OptiX REQUIRED)
    endif()
endif()

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

add_library(OptixCSP_interface_soltrace STATIC
    soltrace_optix_interface.cpp
)

# add interface folder to the include path
target_include_directories(OptixCSP_interface_soltrace
    PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}       # expose interface folder to public
    PRIVATE ${PROJECT_SOURCE_DIR}/src/shaders
            ${PROJECT_SOURCE_DIR}/src/core    
)

target_link_libraries(OptixCSP_interface_soltrace
    PRIVATE OptixCSP_core)
cmake_minimum_required(VERSION 3.20)

project(OptixCSP LANGUAGES CXX CUDA)

if(MSVC)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -Xcompiler=/Zc:preprocessor")
endif()

# C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# CUDA architectures 
set(CMAKE_CUDA_ARCHITECTURES 75;86;89;90)

# Set CMAKE_MODULE_PATH BEFORE add_subdirectory so child can find FindOptiX.cmake
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/")

# Set output directories - needed for PTX output
# Use project-specific output dir to avoid conflicts when built as subdirectory
if(NOT CMAKE_LIBRARY_OUTPUT_DIRECTORY)
    set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
endif()
if(NOT CMAKE_RUNTIME_OUTPUT_DIRECTORY)
    set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")
endif()
if(NOT CMAKE_ARCHIVE_OUTPUT_DIRECTORY)
    set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
endif()

# Only copy data when building standalone
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
    file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data/ DESTINATION ${CMAKE_BINARY_DIR}/bin/data/)
    message(STATUS "Data directory copied to: ${CMAKE_BINARY_DIR}/bin/data/")
endif()

# Add subdirectories
add_subdirectory(src)
#add_subdirectory(demos)

# Create an alias target for consistent naming (works for both in-tree and exported)
if(NOT TARGET OptixCSP::OptixCSP_core)
    add_library(OptixCSP::OptixCSP_core ALIAS OptixCSP_core)
endif()

message(STATUS "Configured OptixCSP project.")

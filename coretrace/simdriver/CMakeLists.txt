project(simdriver)

cmake_minimum_required(VERSION 3.19)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

include_directories(
    .
    ../simulation_data
    ../simulation_results
    ../simulation_runner/native_runner
)

add_executable(simdriver main.cpp)

if(SOLTRACE_BUILD_EMBREE_SUPPORT)
    target_include_directories(simdriver PRIVATE
        ../simulation_runner/embree_runner
    )
    target_compile_definitions(simdriver PRIVATE SOLTRACE_EMBREE_SUPPORT)
    target_link_libraries(simdriver PRIVATE embree_runner)
endif()

if(SOLTRACE_BUILD_OPTIX_SUPPORT)
    target_include_directories(simdriver PRIVATE
        ../simulation_runner/optix_runner
    )
    target_compile_definitions(simdriver PRIVATE SOLTRACE_OPTIX_SUPPORT)
    target_link_libraries(simdriver PRIVATE optix_runner)
endif()

target_link_libraries(simdriver PRIVATE
    simdata
    simresult
    native_runner
)

if(UNIX)
    target_link_libraries(simdriver PRIVATE pthread)
endif()

if (NOT APPLE AND ENABLE_COVERAGE)
    target_link_options(simdriver PRIVATE --coverage)
    target_link_libraries(simdriver PRIVATE gcov)
endif()

install(TARGETS simdriver
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

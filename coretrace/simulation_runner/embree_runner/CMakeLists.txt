Project(embree_runner)

# include(FetchContent)
# FetchContent_Declare(
# 	embree
# 	GIT_REPOSITORY https://github.com/RenderKit/embree.git
# 	GIT_TAG v4.4.0
# )

include_directories(
    .
    ../../simulation_data
    ../../simulation_runner
    ../../simulation_runner/native_runner
    ../../simulation_results
)

set(EMBREE_RUNNER_SRC
    bbox_calculator.cpp
    find_element_hit_embree.cpp
    embree_helper.cpp
    embree_runner.cpp
    trace_embree.cpp
)

set(POSITION_INDEPENDENT_CODE ON)

FIND_PACKAGE(embree 4 REQUIRED)

if(WIN32)
    # Hack to build on Windows
    add_library(embree_runner STATIC ${EMBREE_RUNNER_SRC})
else(WIN32)
    add_library(embree_runner SHARED ${EMBREE_RUNNER_SRC})
endif()

set_target_properties(
    embree_runner
    PROPERTIES
    DEBUG_POSTFIX "d"
	PREFIX ""
)

target_link_libraries(
    embree_runner
    PUBLIC
    simdata
    simresult
    native_runner
    embree
)

if (NOT APPLE AND ENABLE_COVERAGE)
    SET(CMAKE_CXX_FLAGS "-O0 -coverage")
endif()

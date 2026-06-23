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

add_library(embree_runner STATIC ${EMBREE_RUNNER_SRC})

set_target_properties(
    embree_runner
    PROPERTIES
    DEBUG_POSTFIX "d"
	PREFIX ""
)

target_include_directories(embree_runner
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/../../simulation_runner>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SolTrace/simulation_runner/embree_runner>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SolTrace/simulation_runner>
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

install(TARGETS embree_runner
    EXPORT SolTraceTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SolTrace/simulation_runner/embree_runner
    FILES_MATCHING
        PATTERN "*.hpp"
        PATTERN "*.h"
        PATTERN "CMakeLists.txt" EXCLUDE
)


Project(simresult)

include_directories(
    . 
    ../simulation_data
)
# include_directories(cst_templates)

set(SIMRES_SRC
    records.cpp
    simulation_result.cpp
)

set(SIMRES_HDRS
    records.hpp
    simulation_result.hpp
    simulation_result_api.hpp
    simulation_result_export.hpp
)

set(POSITION_INDEPENDENT_CODE ON)

add_library(simresult STATIC ${SIMRES_SRC})

set_target_properties(simresult
    PROPERTIES
    DEBUG_POSTFIX "d"
	PREFIX ""
)

target_include_directories(simresult
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/SolTrace/simulation_results>
)

target_link_libraries(
    simresult
    simdata
)

if (NOT APPLE AND ENABLE_COVERAGE)
	# SET(CMAKE_CXX_FLAGS "-O0 -coverage -fkeep-inline-functions")
    SET(CMAKE_CXX_FLAGS "-O0 -coverage")
endif()

install(TARGETS simresult
    EXPORT SolTraceTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/SolTrace/simulation_results
    FILES_MATCHING
        PATTERN "*.hpp"
        PATTERN "*.h"
        PATTERN "CMakeLists.txt" EXCLUDE
)

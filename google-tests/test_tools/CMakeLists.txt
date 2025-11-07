Project(TestTools)

include_directories(
	.
	../../coretrace
    ../../coretrace/simulation_data
    ../../coretrace/simulation_results
    ../../coretrace/simulation_runner
    ../../coretrace/simulation_runner/native_runner
)

# Test tools files
set(TEST_TOOLS_SRC
	count_absorbed_native.cpp
    split_csv.cpp
)

set(POSITION_INDEPENDENT_CODE ON)

if(WIN32)
    # Hack to build on Windows
    add_library(testtools STATIC ${TEST_TOOLS_SRC})
else(WIN32)
    add_library(testtools SHARED ${TEST_TOOLS_SRC})
endif()

set_target_properties(testtools
    PROPERTIES
    DEBUG_POSTFIX "d"
    PREFIX ""
)

target_link_libraries(
    testtools
    simdata
    simresult
    native_runner
)

if(UNIX)
	target_link_libraries(testtools pthread)
endif()

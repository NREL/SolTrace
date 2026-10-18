
Project(UnitTests)

# Add the separate test directories
add_subdirectory(simulation_data)
add_subdirectory(simulation_results)  
add_subdirectory(simulation_runner)

if(SOLTRACE_BUILD_GUI)
    add_subdirectory(gui)
endif()

if(SOLTRACE_BUILD_API)
    add_subdirectory(stapi_v2)
endif()

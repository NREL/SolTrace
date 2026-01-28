Project(SimulationRunnerUnitTests)

# Add native runner tests (always built)
add_subdirectory(native_runner)

# Add embree runner tests only if Embree is enabled
if(SOLTRACE_BUILD_EMBREE_SUPPORT)
    add_subdirectory(embree_runner)
endif()

# Add optix runner tests only if OptiX is enabled
if(SOLTRACE_BUILD_OPTIX_SUPPORT)
    add_subdirectory(optix_runner)
endif()

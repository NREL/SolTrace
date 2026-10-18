Project(STAPIUnitTests)

# Add native runner tests (always built)
add_subdirectory(native)
message(STATUS "[stapi_v2_tests] - Added tests for stapi_v2 with only native runner")

# Add embree runner tests only if Embree is enabled
if(SOLTRACE_BUILD_EMBREE_SUPPORT)
    add_subdirectory(embree)
    message(STATUS "[stapi_v2_tests] - Added tests for stapi_v2 with embree runner")
endif()

# Add optix runner tests only if OptiX is enabled
if(SOLTRACE_BUILD_OPTIX_SUPPORT)
    add_subdirectory(optix)
    message(STATUS "[stapi_v2_tests] - Added tests for stapi_v2 with optix runner")
endif()

# Add all runner tests if both are enabled
if(SOLTRACE_BUILD_EMBREE_SUPPORT AND SOLTRACE_BUILD_OPTIX_SUPPORT)
    add_subdirectory(all_runners)
    message(STATUS "[stapi_v2_tests] - Added tests for stapi_v2 with all runners")
endif()

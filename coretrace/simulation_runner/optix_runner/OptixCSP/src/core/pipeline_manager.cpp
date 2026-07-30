#include "pipeline_manager.h"

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <cuda_runtime.h>
#include <optix.h>
#include <sampleConfig.h>
#include <optix_stack_size.h>
#include <optix_stubs.h>

#include "utils/util_check.hpp"
#include "shaders/Soltrace.h"

#include "data_manager.h"
#include "soltrace_state.h"

#if defined(_WIN32)
#include <windows.h>
#elif defined(__APPLE__)
#include <mach-o/dyld.h>
#else
#include <limits.h>
#include <unistd.h>
#endif

using namespace OptixCSP;

char LOG[2048] = {}; // A mutable log buffer.
size_t LOG_SIZE = sizeof(LOG);

// const char* intersectionFuncs[] = {
//     "__intersection__rectangle_parabolic",
//     "__intersection__rectangle_flat",
//     "__intersection__triangle_flat",
//     // "__intersection__cylinder_y_capped"
//     "__intersection_cylinder_y",
// };

// const std::map<SurfaceApertureMap, std::string> IntersectionKernelMap = {
//     {SurfaceApertureMap(SurfaceType::PARABOLIC, ApertureType::RECTANGLE), "__intersection__rectangle_parabolic"},
//     {SurfaceApertureMap(SurfaceType::FLAT, ApertureType::RECTANGLE), "__intersection__rectangle_flat"},
//     {SurfaceApertureMap(SurfaceType::FLAT, ApertureType::TRIANGLE), "__intersection__triangle_flat"},
//     {SurfaceApertureMap(SurfaceType::CYLINDER, ApertureType::RECTANGLE), "__intersection__cylinder_y"},
//     {SurfaceApertureMap(SurfaceType::FLAT, ApertureType::QUADRILATERAL), "__intersection__flat_quadrilateral"}};

const std::map<OpticalEntityType, std::string> IntersectionKernelMap = {
    {OpticalEntityType::RECTANGLE_PARABOLIC, "__intersection__rectangle_parabolic"},
    {OpticalEntityType::RECTANGLE_FLAT, "__intersection__rectangle_flat"},
    {OpticalEntityType::TRIANGLE_FLAT, "__intersection__triangle_flat"},
    {OpticalEntityType::CYLINDRICAL, "__intersection__cylinder_y"},
    {OpticalEntityType::QUADRILATERAL_FLAT, "__intersection__quadrilateral_flat"},
    {OpticalEntityType::CIRCLE_FLAT, "__intersection__circle_flat"},
    {OpticalEntityType::HEXAGON_FLAT, "__intersection__hexagon_flat"},
    {OpticalEntityType::ANNULUS_FLAT, "__intersection__annulus_flat"},
    {OpticalEntityType::CIRCLE_PARABOLIC, "__intersection__circle_parabolic"},
    {OpticalEntityType::HEXAGON_PARABOLIC, "__intersection__hexagon_parabolic"},
    {OpticalEntityType::TRIANGLE_PARABOLIC, "__intersection__triangle_parabolic"},
    {OpticalEntityType::ANNULUS_PARABOLIC, "__intersection__annulus_parabolic"},
    {OpticalEntityType::QUADRILATERAL_PARABOLIC, "__intersection__quadrilateral_parabolic"},
    {OpticalEntityType::RECTANGLE_SPHERICAL, "__intersection__rectangle_spherical"},
    {OpticalEntityType::CIRCLE_SPHERICAL, "__intersection__circle_spherical"},
    {OpticalEntityType::HEXAGON_SPHERICAL, "__intersection__hexagon_spherical"},
    {OpticalEntityType::ANNULUS_SPHERICAL, "__intersection__annulus_spherical"},
    {OpticalEntityType::TRIANGLE_SPHERICAL, "__intersection__triangle_spherical"},
    {OpticalEntityType::QUADRILATERAL_SPHERICAL, "__intersection__quadrilateral_spherical"}};

pipelineManager::pipelineManager(SoltraceState &state) : m_state(state) {}

pipelineManager::~pipelineManager() noexcept
{
    try
    {
        cleanup();
    }
    catch (const std::exception &error)
    {
        std::cerr << "[OptixRunner] Pipeline cleanup failed during destruction: "
                  << error.what() << '\n';
    }
    catch (...)
    {
        std::cerr << "[OptixRunner] Pipeline cleanup failed during destruction with an unknown error.\n";
    }
}

namespace
{

std::filesystem::path executableDirectory()
{
#if defined(_WIN32)
    std::vector<char> buffer(MAX_PATH);
    DWORD length = 0;
    while (true)
    {
        length = GetModuleFileNameA(nullptr, buffer.data(), static_cast<DWORD>(buffer.size()));
        if (length == 0)
            return {};
        if (length < buffer.size())
            break;
        buffer.resize(buffer.size() * 2);
    }
    return std::filesystem::path(buffer.data()).parent_path();
#elif defined(__APPLE__)
    uint32_t size = 0;
    _NSGetExecutablePath(nullptr, &size);
    std::vector<char> buffer(size);
    if (_NSGetExecutablePath(buffer.data(), &size) != 0)
        return {};
    return std::filesystem::weakly_canonical(buffer.data()).parent_path();
#else
    std::vector<char> buffer(PATH_MAX);
    ssize_t length = readlink("/proc/self/exe", buffer.data(), buffer.size() - 1);
    if (length <= 0)
        return {};
    buffer[length] = '\0';
    return std::filesystem::path(buffer.data()).parent_path();
#endif
}

std::vector<std::filesystem::path> ptxSearchDirectories()
{
    std::vector<std::filesystem::path> directories;

    if (const char *env_path = std::getenv("SOLTRACE_PTX_DIR"))
    {
        if (*env_path != '\0')
            directories.emplace_back(env_path);
    }

    directories.emplace_back(SAMPLES_PTX_DIR);

    const auto cwd = std::filesystem::current_path();
    directories.emplace_back(cwd / "ptx");
    directories.emplace_back(cwd / "lib" / "ptx");

    const auto exe_dir = executableDirectory();
    if (!exe_dir.empty())
    {
        directories.emplace_back(exe_dir / "ptx");
        directories.emplace_back(exe_dir / "lib" / "ptx");
        directories.emplace_back(exe_dir / ".." / "lib" / "ptx");
        directories.emplace_back(exe_dir / ".." / "share" / "SolTrace" / "ptx");
    }

    return directories;
}

} // namespace

void pipelineManager::cleanup()
{
    if (m_state.pipeline)
    {
        OPTIX_CHECK(optixPipelineDestroy(m_state.pipeline));
        m_state.pipeline = nullptr;
    }

    // destroy all program groups
    for (auto &prog_group : m_program_groups)
    {
        if (prog_group)
        {
            OPTIX_CHECK(optixProgramGroupDestroy(prog_group));
        }
    }
    m_program_groups.clear();
    m_intersection_program_group_map.clear();

    m_state.raygen_prog_group = nullptr;
    m_state.radiance_miss_prog_group = nullptr;
    m_state.radiance_receiver_prog_group = nullptr;

    if (m_state.geometry_module)
    {
        OPTIX_CHECK(optixModuleDestroy(m_state.geometry_module));
        m_state.geometry_module = nullptr;
    }

    if (m_state.shading_module)
    {
        OPTIX_CHECK(optixModuleDestroy(m_state.shading_module));
        m_state.shading_module = nullptr;
    }

    if (m_state.sun_module)
    {
        OPTIX_CHECK(optixModuleDestroy(m_state.sun_module));
        m_state.sun_module = nullptr;
    }
}

std::string pipelineManager::loadPtxFromFile(const std::string &kernelName)
{
    const auto filename = kernelName + ".ptx";
    std::stringstream attempted_paths;

    for (const auto &directory : ptxSearchDirectories())
    {
        const auto ptx_file = directory / filename;
        std::ifstream file(ptx_file);
        if (!file.good())
        {
            attempted_paths << "\n  " << ptx_file.string();
            continue;
        }

        std::stringstream buffer;
        buffer << file.rdbuf();
        return buffer.str();
    }

    throw std::runtime_error("PTX file not found. Tried:" + attempted_paths.str());
}

void pipelineManager::loadModules()
{
    OptixModuleCompileOptions moduleCompileOptions = {};
#if !defined(NDEBUG)
    moduleCompileOptions.optLevel = OPTIX_COMPILE_OPTIMIZATION_LEVEL_0;
    moduleCompileOptions.debugLevel = OPTIX_COMPILE_DEBUG_LEVEL_FULL;
#endif

    {
        std::string ptx = loadPtxFromFile("intersection");
        LOG_SIZE = sizeof(LOG);
        // We are temporarily replacing OPTIX_CHECK with a manual check
        OptixResult result = optixModuleCreate(
            m_state.context,
            &moduleCompileOptions,
            &m_state.pipeline_compile_options,
            ptx.c_str(),
            ptx.size(),
            LOG, &LOG_SIZE,
            &m_state.geometry_module);

        // If it fails, print the REAL error message from the LOG buffer
        if (result != OPTIX_SUCCESS)
        {
            std::cerr << "--- OPTIX COMPILATION LOG ---\n"
                      << std::string(LOG, LOG_SIZE)
                      << "\n--- END LOG ---\n";
            // Now, re-throw the error so the test still fails
            throw std::runtime_error("optixModuleCreate failed for intersection.ptx");
        }
    }
    // Shading/materials module.
    {
        std::string ptx = loadPtxFromFile("materials");
        LOG_SIZE = sizeof(LOG);
        OPTIX_CHECK(optixModuleCreate(
            m_state.context,
            &moduleCompileOptions,
            &m_state.pipeline_compile_options,
            ptx.c_str(),
            ptx.size(),
            LOG, &LOG_SIZE,
            &m_state.shading_module));
    }
    // Sun module.
    {
        std::string ptx = loadPtxFromFile("sun");
        LOG_SIZE = sizeof(LOG);
        OPTIX_CHECK(optixModuleCreate(
            m_state.context,
            &moduleCompileOptions,
            &m_state.pipeline_compile_options,
            ptx.c_str(),
            ptx.size(),
            LOG, &LOG_SIZE,
            &m_state.sun_module));
    }
}

void pipelineManager::createPipeline()
{
    m_state.pipeline_compile_options = {
        false,                                           // usesMotionBlur: Disable motion blur.
        OPTIX_TRAVERSABLE_GRAPH_FLAG_ALLOW_SINGLE_GAS,   // traversableGraphFlags: Allow only a single GAS.
        2, /* RadiancePRD uses 5 payloads */             // numPayloadValues
        5, /* Parallelogram intersection uses 5 attrs */ // numAttributeValues
        OPTIX_EXCEPTION_FLAG_NONE,                       // exceptionFlags
        "params"                                         // pipelineLaunchParamsVariableName
    };

    // Prepare modules and program groups
    loadModules();
    createSunProgram();
    createElementPrograms();
    createMissProgram();

    // Link program groups to pipeline
    OptixPipelineLinkOptions pipeline_link_options = {};
    pipeline_link_options.maxTraceDepth = m_max_trace_depth; // Maximum recursion depth for ray tracing.

    // Create the OptiX pipeline by linking the program groups.
    LOG_SIZE = sizeof(LOG);
    OPTIX_CHECK(optixPipelineCreate(
        m_state.context,                                    // OptiX context.
        &m_state.pipeline_compile_options,                  // Compile options for the pipeline.
        &pipeline_link_options,                             // Link options for the pipeline.
        m_program_groups.data(),                            // Array of program groups.
        static_cast<unsigned int>(m_program_groups.size()), // Number of program groups.
        LOG, &LOG_SIZE,                                     // Logs for diagnostics.
        &m_state.pipeline                                   // Output: Handle for the created pipeline.
        ));

    // Compute and configure the stack sizes for the pipeline.
    OptixStackSizes stack_sizes = {};
    for (auto &prog_group : m_program_groups)
    {
        OPTIX_CHECK(optixUtilAccumulateStackSizes(prog_group, &stack_sizes, m_state.pipeline));
    }

    uint32_t direct_callable_stack_size_from_traversal;
    uint32_t direct_callable_stack_size_from_state;
    uint32_t continuation_stack_size;

    // Compute stack sizes based on the maximum trace depth and other settings.
    OPTIX_CHECK(optixUtilComputeStackSizes(
        &stack_sizes,                               // Input stack sizes.
        m_max_trace_depth,                          // Maximum trace depth.
        0,                                          // maxCCDepth: Maximum depth of continuation callables (none in this case).
        0,                                          // maxDCDepth: Maximum depth of direct callables (none in this case).
        &direct_callable_stack_size_from_traversal, // Output: Stack size for callable traversal.
        &direct_callable_stack_size_from_state,     // Output: Stack size for callable state.
        &continuation_stack_size                    // Output: Stack size for continuation stack.
        ));
    // Set the computed stack sizes for the pipeline.
    OPTIX_CHECK(optixPipelineSetStackSize(
        m_state.pipeline,                          // Pipeline to configure.
        direct_callable_stack_size_from_traversal, // Stack size for direct callable traversal.
        direct_callable_stack_size_from_state,     // Stack size for direct callable state.
        continuation_stack_size,                   // Stack size for continuation stack.
        1                                          // maxTraversableDepth: Maximum depth of traversable hierarchy.
        ));
}

OptixPipeline pipelineManager::getPipeline() const
{
    return m_state.pipeline;
}

void pipelineManager::createHitGroupProgram(OptixProgramGroup &group,
                                            OptixModule intersectionModule, const char *intersectionFunc,
                                            OptixModule closestHitModule, const char *closestHitFunc)
{
    OptixProgramGroupOptions options = {};
    OptixProgramGroupDesc desc = {};
    desc.kind = OPTIX_PROGRAM_GROUP_KIND_HITGROUP;
    desc.hitgroup.moduleIS = intersectionModule;
    desc.hitgroup.entryFunctionNameIS = intersectionFunc;
    desc.hitgroup.moduleCH = closestHitModule;
    desc.hitgroup.entryFunctionNameCH = closestHitFunc;
    desc.hitgroup.moduleAH = nullptr;
    desc.hitgroup.entryFunctionNameAH = nullptr;

    LOG_SIZE = sizeof(LOG);
    OPTIX_CHECK(optixProgramGroupCreate(
        m_state.context,
        &desc,
        1,
        &options,
        LOG, &LOG_SIZE,
        &group));
}

// TODO: simplify
void pipelineManager::createSunProgram()
{
    OptixProgramGroup group;               // Handle for the sun program group.
    OptixProgramGroupOptions options = {}; // Options for the program group
    OptixProgramGroupDesc desc = {};       // Descriptor to define the program group.

    // Specify the kind of program group (Ray Generation).
    desc.kind = OPTIX_PROGRAM_GROUP_KIND_RAYGEN;
    // Link the ray generation program to the sun module and specify the function name.
    desc.raygen.module = m_state.sun_module;
    desc.raygen.entryFunctionName = "__raygen__sun_source";

    // Create the program group
    // Note: OPTIX_CHECK_LOG is not used here because the macro creates its own
    // local LOG_ buffer (not the global LOG), causing it to always print 2048
    // null bytes to stderr. Instead we use OPTIX_CHECK and manually print any
    // non-empty log content when verbose mode is enabled.
    LOG_SIZE = sizeof(LOG);
    OPTIX_CHECK(optixProgramGroupCreate(
        m_state.context, // OptiX context.
        &desc,           // Descriptor defining the program group.
        1,               // Number of program groups to create (1 in this case).
        &options,        // Options for the program group.
        LOG, &LOG_SIZE,  // Logs to capture diagnostic information.
        &group           // Output: Handle for the created program group.
        ));
    if (LOG_SIZE > 1 && LOG[0] != '\0')
    {
        std::cerr << "OptiX log for optixProgramGroupCreate (sun):\n"
                  << std::string(LOG, LOG + LOG_SIZE) << std::endl;
    }

    m_program_groups.push_back(group);
    m_state.raygen_prog_group = group;
}

// Create program group for handling rays interacting with elements.
void pipelineManager::createElementPrograms()
{

    // number of element programs
    // size_t numElementPrograms = sizeof(intersectionFuncs) / sizeof(intersectionFuncs[0]);

    // for (size_t i = 0; i < numElementPrograms; i++) {
    // 	OptixProgramGroup group;

    // 	createHitGroupProgram(group,
    //         			      m_state.geometry_module,
    // 			              intersectionFuncs[i],
    //         			      m_state.shading_module,
    // 			              "__closesthit__element");

    // 	m_program_groups.push_back(group);
    // }

    m_intersection_program_group_map.clear();
    size_t idx;
    for (const auto &[optype, kernel_name] : IntersectionKernelMap)
    {
        OptixProgramGroup group;

        createHitGroupProgram(group,
                              m_state.geometry_module,
                              kernel_name.c_str(),
                              m_state.shading_module,
                              "__closesthit__element");

        idx = m_program_groups.size();

        // std::cout << "Key: " << optype
        // 	  << " Kernel: " << kernel_name
        // 	  << " Program Group Index: " << idx
        // 	  << std::endl;

        auto sts = m_intersection_program_group_map.insert_or_assign(optype, idx);
        if (!sts.second)
        {
            // This should be impossible since we are iterating over a map with the same
            // key as we are inserting here.
            throw std::runtime_error("Duplicate optical entity!");
        }

        m_program_groups.push_back(group);
    }
}

// Create program group for handling rays that miss all geometry.
void pipelineManager::createMissProgram()
{
    OptixProgramGroup group;               // Handle for the miss program group.
    OptixProgramGroupOptions options = {}; // Options for the program group (none needed).
    OptixProgramGroupDesc desc = {};       // Descriptor for the program group.

    // Specify the kind of program group (Miss Program for handling missed rays).
    desc.kind = OPTIX_PROGRAM_GROUP_KIND_MISS;
    // Link the miss shader (background or environment shading) to the shading module.
    desc.miss.module = m_state.shading_module;
    desc.miss.entryFunctionName = "__miss__ms";

    // Create the program grou
    LOG_SIZE = sizeof(LOG);
    OPTIX_CHECK(optixProgramGroupCreate(
        m_state.context,
        &desc,
        1,
        &options,
        LOG, &LOG_SIZE,
        &group));

    m_program_groups.push_back(group);
    m_state.radiance_miss_prog_group = group;

    // now let's print out program groups address
    for (int i = 0; i < m_program_groups.size(); i++)
    {
        OptixProgramGroup group = m_program_groups[i];
        if (m_verbose)
        {
            std::cout << "Program group " << i << " address: " << group << ", kind: " << desc.kind << std::endl;
        }
    }
}

OptixProgramGroup pipelineManager::getElementProgram(OpticalEntityType map) const
{
    auto iter = m_intersection_program_group_map.find(map);
    if (iter != m_intersection_program_group_map.cend())
    {
        return m_program_groups[iter->second];
    }
    throw std::runtime_error("Unsupported surface or aperture type in getElementProgram");
}

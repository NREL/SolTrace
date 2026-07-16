#pragma once
#include <map>
#include <string>
#include <vector>
#include "optix.h"
#include "shaders/Soltrace.h"
#include "soltrace_type.h"
#include "soltrace_state.h"

namespace OptixCSP
{
    /**
     * @class pipelineManager
     * @brief Manages the OptiX pipeline, including loading modules, creating programs, and handling cleanup.
     */
    class pipelineManager
    {
    public:
        /**
         * @brief Constructs a pipelineManager instance.
         * @param state Reference to SoltraceState for managing OptiX state.
         */
        pipelineManager(SoltraceState &state);

        /**
         * @brief Destroys the pipelineManager and cleans up resources.
         */
        ~pipelineManager() noexcept;

        /**
         * @brief Cleans up all allocated resources including OptiX programs and modules.
         */
        void cleanup();

        void set_verbose(bool verbose) { m_verbose = verbose; }

        void set_max_trace_depth(uint8_t depth) { m_max_trace_depth = depth; }

        /**
         * @brief Loads PTX code from a file.
         * @param kernelName The name of the PTX kernel file to load.
         * @return The loaded PTX code as a string.
         */
        std::string loadPtxFromFile(const std::string &kernelName);

        /**
         * @brief Loads all required OptiX modules.
         */
        void loadModules();

        /**
         * @brief Creates the OptiX pipeline with associated program groups.
         */
        void createPipeline();

        /**
         * @brief Retrieves the OptiX pipeline object.
         * @return The OptiX pipeline.
         */
        OptixPipeline getPipeline() const;

        /**
         * @brief Creates the ray generation program for the sun.
         *
         * Currently, there is only one ray generation source; this function can be extended for multiple sources.
         */
        void createSunProgram();

        /**
         * @brief Creates programs for element interactions.
         */
        void createElementPrograms();

        /**
         * @brief Creates the miss program, handling rays that do not hit any geometry.
         */
        void createMissProgram();

        /**
         * @brief Helper function to create a hit group program given intersection and closest hit functions.
         * @param group Reference to an OptixProgramGroup to be created.
         * @param intersectionModule OptiX module containing the intersection function.
         * @param intersectionFunc Name of the intersection function.
         * @param closestHitModule OptiX module containing the closest hit function.
         * @param closestHitFunc Name of the closest hit function.
         */
        void createHitGroupProgram(OptixProgramGroup &group,
                                   OptixModule intersectionModule,
                                   const char *intersectionFunc,
                                   OptixModule closestHitModule,
                                   const char *closestHitFunc);

        // /**
        //  * @brief Retrieves the ray generation program group.
        //  * @return The OptixProgramGroup for ray generation.
        //  */
        // OptixProgramGroup getRaygenProgram() const;

        // /**
        //  * @brief Retrieves the miss program group.
        //  * @return The OptixProgramGroup for the miss shader.
        //  */
        // OptixProgramGroup getMissProgram() const;

        /**
         * @brief Retrieves the appropriate element program group based on surface and aperture type.
         * @param map OpticalEntityType specifying the surface and aperture combination.
         * @return The corresponding OptixProgramGroup.
         */
        OptixProgramGroup getElementProgram(OpticalEntityType map) const;

    private:
        SoltraceState &m_state;                                               ///< Reference to the simulation's OptiX state.
        std::vector<OptixProgramGroup> m_program_groups;                      ///< Stores all created OptiX program groups.
        std::map<OpticalEntityType, size_t> m_intersection_program_group_map; ///< Map surface-aperture combinations to index in m_program_groups
        bool m_verbose = false;
        uint8_t m_max_trace_depth = DEFAULT_MAX_TRACE_DEPTH;

        // // Number of program groups categorized by type.
        // int num_raygen_programs = 1; ///< Number of ray generation programs.
        // int num_heliostat_programs = 4; ///< Number of heliostat-related programs.
        // int num_miss_programs = 1; ///< Number of miss programs.
    };
}

#include <gtest/gtest.h>

#include <fstream>
#include <string>

#include <native_runner.hpp>
#include <optix_runner.hpp>
#include <simulation_data.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result_export.hpp>

using SolTrace::Result::RayEvent;
using SolTrace::Runner::RunnerStatus;

// -----------------------------------------------------------------------
// Tests that the GPU (OptiX) runner can load and execute the High Flux
// Solar Furnace example file without errors.
//
// The HFSF file uses spherical surfaces with hexagonal apertures, which
// requires the SPHERICAL surface type to be supported by the GPU runner.
// -----------------------------------------------------------------------

TEST(HighFluxSolarFurnace, OptixRunnerCanRun)
{
    // Load input file
    std::string sample_path =
        std::string(PROJECT_DIR) +
        std::string("/high_flux_solar_furnace_test.stinput");

    SimulationData sd;
    bool           ok = sd.import_from_file(sample_path);
    ASSERT_TRUE(ok) << "Failed to load HFSF stinput file: " << sample_path;
    EXPECT_GT(sd.get_number_of_elements(), 0u);
    EXPECT_GT(sd.get_number_of_ray_sources(), 0u);

    // Configure simulation
    const uint_fast64_t   NRAYS     = 10000;
    SimulationParameters& params    = sd.get_simulation_parameters();
    params.include_optical_errors   = false;
    params.include_sun_shape_errors = false;
    params.max_number_of_rays       = NRAYS * 100;
    params.number_of_rays           = NRAYS;
    params.seed                     = 1;

    // Run GPU simulation
    OptixRunner  runner;
    RunnerStatus sts = runner.initialize();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    sts = runner.setup_simulation(&sd);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    sts = runner.run_simulation();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    // Collect results
    SimulationResult result;
    sts = runner.report_simulation(&result, 0);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    // Verify that rays were traced and results are available
    EXPECT_EQ(result.get_number_of_records(), static_cast<int>(NRAYS));

    // Verify some rays hit something (absorbed or reflected)
    int total_hits = 0;
    int n_records  = result.get_number_of_records();
    for (int i = 0; i < n_records; ++i)
    {
        ray_record_ptr rec            = result[i];
        int            n_interactions = rec->get_number_of_interactions();
        for (int j = 0; j < n_interactions; ++j)
        {
            RayEvent rev = rec->get_event(j);
            if (rev == RayEvent::ABSORB || rev == RayEvent::REFLECT)
                ++total_hits;
        }
    }
    EXPECT_GT(total_hits, 0)
        << "No ray interactions recorded – GPU trace likely failed silently";
}

TEST(HighFluxSolarFurnace, OptixRunnerResultsConsistentWithNativeRunner)
{
    // Load input file
    std::string sample_path =
        std::string(PROJECT_DIR) +
        std::string("/high_flux_solar_furnace_test.stinput");

    const uint_fast64_t NRAYS = 10000;

    // Run native runner first to get reference hit counts
    int native_total_hits = 0;
    {
        SimulationData sd;
        ASSERT_TRUE(sd.import_from_file(sample_path));

        SimulationParameters& params    = sd.get_simulation_parameters();
        params.include_optical_errors   = false;
        params.include_sun_shape_errors = false;
        params.max_number_of_rays       = NRAYS * 100;
        params.number_of_rays           = NRAYS;
        params.seed                     = 1;

        SolTrace::NativeRunner::NativeRunner native_runner;
        native_runner.set_number_of_threads(1);
        native_runner.disable_stages();
        RunnerStatus sts = native_runner.initialize();
        ASSERT_EQ(sts, RunnerStatus::SUCCESS);
        sts = native_runner.setup_simulation(&sd);
        ASSERT_EQ(sts, RunnerStatus::SUCCESS);
        sts = native_runner.run_simulation();
        ASSERT_EQ(sts, RunnerStatus::SUCCESS);

        SimulationResult result;
        sts = native_runner.report_simulation(&result, 0);
        ASSERT_EQ(sts, RunnerStatus::SUCCESS);
        ASSERT_EQ(result.get_number_of_records(), static_cast<int>(NRAYS));

        for (int i = 0; i < result.get_number_of_records(); ++i)
        {
            ray_record_ptr rec = result[i];
            for (int j = 0; j < rec->get_number_of_interactions(); ++j)
            {
                RayEvent rev = rec->get_event(j);
                if (rev == RayEvent::ABSORB || rev == RayEvent::REFLECT)
                    ++native_total_hits;
            }
        }
    }

    // Run GPU runner
    int optix_total_hits = 0;
    {
        SimulationData sd;
        ASSERT_TRUE(sd.import_from_file(sample_path));

        SimulationParameters& params    = sd.get_simulation_parameters();
        params.include_optical_errors   = false;
        params.include_sun_shape_errors = false;
        params.max_number_of_rays       = NRAYS * 100;
        params.number_of_rays           = NRAYS;
        params.seed                     = 1;

        OptixRunner  runner;
        RunnerStatus sts = runner.initialize();
        ASSERT_EQ(sts, RunnerStatus::SUCCESS);
        sts = runner.setup_simulation(&sd);
        ASSERT_EQ(sts, RunnerStatus::SUCCESS);
        sts = runner.run_simulation();
        ASSERT_EQ(sts, RunnerStatus::SUCCESS);

        SimulationResult result;
        sts = runner.report_simulation(&result, 0);
        ASSERT_EQ(sts, RunnerStatus::SUCCESS);
        ASSERT_EQ(result.get_number_of_records(), static_cast<int>(NRAYS));

        for (int i = 0; i < result.get_number_of_records(); ++i)
        {
            ray_record_ptr rec = result[i];
            for (int j = 0; j < rec->get_number_of_interactions(); ++j)
            {
                RayEvent rev = rec->get_event(j);
                if (rev == RayEvent::ABSORB || rev == RayEvent::REFLECT)
                    ++optix_total_hits;
            }
        }
    }

    // Both runners should see a similar number of total interactions.
    // Allow a 5% relative tolerance due to floating-point and stochastic
    // differences.
    EXPECT_GT(native_total_hits, 0);
    EXPECT_GT(optix_total_hits, 0);

    const double tolerance = 0.1;
    const double relative_diff =
        fabs(static_cast<double>(optix_total_hits - native_total_hits) /
             static_cast<double>(native_total_hits));
    EXPECT_NEAR(relative_diff, 0.0, tolerance)
        << "GPU hit count (" << optix_total_hits
        << ") diverges from native hit count (" << native_total_hits
        << ") by more than " << (tolerance * 100) << "%";
}

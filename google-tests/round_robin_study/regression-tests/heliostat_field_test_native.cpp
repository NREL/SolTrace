#include "heliostat_field_test_template.hpp"

#include <native_runner.hpp>

using NativeRunnerType = SolTrace::NativeRunner::NativeRunner;

using HeliostatFieldSimulationNative = HeliostatFieldSimulation<NativeRunnerType>;

static const int N_threads = static_cast<int>(std::max(1u, std::min(std::thread::hardware_concurrency(), 10u)));

TEST_F(HeliostatFieldSimulationNative, singleFacet_SlantFocused)
{
    this->runner.set_number_of_threads(N_threads);

    // Centerline aimpoints
    create_heliostat_field();
    setup_simData();
    simulate_check_outputs("1a", "1", "8");
    simulate_check_outputs("1a", "1", "12");

    // Scatter aimpoints
    set_scatter_aimpoints();
    simulate_check_outputs("1a", "2", "8");
    simulate_check_outputs("1a", "2", "12");
}

TEST_F(HeliostatFieldSimulationNative, singleFacet_BandFocused)
{
    this->runner.set_number_of_threads(N_threads);

    // Centerline aimpoints
    create_heliostat_field();
    assign_focal_lengths_canting_banded();
    setup_simData();
    simulate_check_outputs("1b", "1", "8");
    simulate_check_outputs("1b", "1", "12");

    // Scatter aimpoints
    set_scatter_aimpoints();
    simulate_check_outputs("1b", "2", "8");
    simulate_check_outputs("1b", "2", "12");
}

TEST_F(HeliostatFieldSimulationNative, multiFacet_SlantCanted)
{
    this->runner.set_number_of_threads(N_threads);

    // Centerline aimpoints
    create_heliostat_field();
    assign_canted_slant(true);      // Flat facets

    setup_simData();
    simulate_check_outputs("2a", "1", "8");
    simulate_check_outputs("2a", "1", "12");

    // Scatter aimpoints
    set_scatter_aimpoints();
    simulate_check_outputs("2a", "2", "8");
    simulate_check_outputs("2a", "2", "12");
}

TEST_F(HeliostatFieldSimulationNative, multiFacet_BandCanted)
{
    this->runner.set_number_of_threads(N_threads);

    // Centerline aimpoints
    create_heliostat_field();
    assign_canted_banded(true);     // Flat facets

    setup_simData();
    simulate_check_outputs("2b", "1", "8");
    simulate_check_outputs("2b", "1", "12");

    // Scatter aimpoints
    set_scatter_aimpoints();
    simulate_check_outputs("2b", "2", "8");
    simulate_check_outputs("2b", "2", "12");
}

TEST_F(HeliostatFieldSimulationNative, multiFacet_SlantFocused_SlantCanted)
{
    this->runner.set_number_of_threads(N_threads);

    // Center aimpoints;
    create_heliostat_field();
    assign_canted_slant(false);     // Slant focused (default)

    setup_simData();
    simulate_check_outputs("3a", "1", "8");
    simulate_check_outputs("3a", "1", "12");

    // Scatter aimpoints
    set_scatter_aimpoints();
    simulate_check_outputs("3a", "2", "8");
    simulate_check_outputs("3a", "2", "12");
}

TEST_F(HeliostatFieldSimulationNative, multiFacet_BandFocused_BandCanted)
{
    this->runner.set_number_of_threads(N_threads);

    // Center aimpoints
    create_heliostat_field();
    assign_canted_banded(false);    // Canted by band
    assign_focal_lengths_banded();  // Focused by band

    setup_simData();
    simulate_check_outputs("3b", "1", "8");
    simulate_check_outputs("3b", "1", "12");

    // Scatter aimpoints
    set_scatter_aimpoints();
    simulate_check_outputs("3b", "2", "8");
    simulate_check_outputs("3b", "2", "12");
}

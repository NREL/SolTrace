#include <gtest/gtest.h>

#include <iostream>
#include <string>

#include <native_runner.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result_export.hpp>
#include <simulation_runner.hpp>

TEST(NativeRunner, Multithreading)
{
    const uint_fast64_t NRAYS = 50000;
    const uint_fast64_t NTHREADS = 8;

    // TODO: Make path work with cmake system
    // std::string input_file = "../../input-files/ANU_NREL_roundRobin_fullField_Task_2a_12_longer_aimpoint.stinput";
    std::string input_file = std::string(PROJECT_DIR) + std::string("/Power-tower-surround_singlefacet.stinput");

    std::cout << "Loading file: " << input_file << std::endl;
    SimulationData sd;
    ASSERT_TRUE(sd.import_from_file(input_file));

    std::cout << "Number of elements: " << sd.get_number_of_elements()
              << "\nNumber of rays: " << NRAYS
              << std::endl;

    SimulationParameters &params = sd.get_simulation_parameters();
    params.include_optical_errors = false;
    params.include_sun_shape_errors = false;
    params.max_number_of_rays = NRAYS * 100;
    params.number_of_rays = NRAYS;
    params.seed = 1;

    SolTrace::Runner::RunnerStatus sts;
    SolTrace::NativeRunner::NativeRunner nr;
    nr.disable_point_focus();
    nr.disable_power_tower();
    nr.set_number_of_threads(NTHREADS);

    std::cout << "Initializing NativeRunner..." << std::endl;
    sts = nr.initialize();
    ASSERT_EQ(sts, SolTrace::Runner::RunnerStatus::SUCCESS);

    std::cout << "Setting up simulation..." << std::endl;
    sts = nr.setup_simulation(&sd);
    ASSERT_EQ(sts, SolTrace::Runner::RunnerStatus::SUCCESS);

    std::cout << "Running simulation..." << std::endl;
    auto t0 = std::chrono::high_resolution_clock::now();
    sts = nr.run_simulation();
    auto t1 = std::chrono::high_resolution_clock::now();
    ASSERT_EQ(sts, SolTrace::Runner::RunnerStatus::SUCCESS);

    std::cout << "Done." << std::endl;

    std::chrono::duration<double, std::milli> dur = t1 - t0;
    std::cout << "Time: " << dur.count() << " ms" << std::endl;

    EXPECT_TRUE(dur.count() < 20000.0);

    // element_id absorber_id = 6285;
    // int_fast64_t nabsorbed = count_element_event(result, absorber_id, RayEvent::ABSORB);
    // int_fast64_t nreflect = count_element_event(result, absorber_id, RayEvent::REFLECT);
    // int_fast64_t nevents = nabsorbed + nreflect;

    // std::cout << "Total: " << nevents
    //           << "\nAbsorb: " << nabsorbed << " ("
    //           << static_cast<double>(nabsorbed) / nevents << ")"
    //           << "\nReflect: " << nreflect << " ("
    //           << static_cast<double>(nreflect) / nevents << ")"
    //           << std::endl;
}

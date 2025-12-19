#include <gtest/gtest.h>

#include <chrono>
#include <future>
#include <thread>

#include <native_runner.hpp>
#include <simulation_data_export.hpp>
#include <simulation_runner.hpp>
#include <simulation_result_export.hpp>

using SolTrace::Runner::RunnerStatus;
using SolTrace::NativeRunner::NativeRunner;

TEST(NativeRunner, StatusAndCancelMultiThread)
{
    std::string sample_path = std::string(PROJECT_DIR) +
                              std::string("/Power-tower-surround_singlefacet.stinput");

    SimulationData sd;
    EXPECT_TRUE(sd.import_from_file(sample_path));
    sd.set_number_of_rays(100000);

    NativeRunner runner;
    runner.disable_point_focus();
    runner.disable_power_tower();
    runner.set_number_of_threads(2);
    RunnerStatus sts;
    sts = runner.setup_simulation(&sd);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);

    auto t0 = std::chrono::high_resolution_clock::now();

    auto fsts = std::async(&NativeRunner::run_simulation, &runner);

    std::this_thread::sleep_for(std::chrono::milliseconds(200));
    sts = runner.status_simulation();
    EXPECT_EQ(sts, RunnerStatus::RUNNING);

    double prog;
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    sts = runner.status_simulation(&prog);
    EXPECT_EQ(sts, RunnerStatus::RUNNING);
    EXPECT_LE(prog, 1.0);
    EXPECT_GE(prog, 0.0);

    runner.cancel_simulation();
    fsts.wait();

    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> dur = t1 - t0;

    EXPECT_EQ(fsts.get(), RunnerStatus::CANCEL);
    EXPECT_LT(dur.count(), 10000.0);

    std::cout << "Time for run: " << dur.count() << std::endl;
    std::cout << "Progress before cancel: " << prog << std::endl;
}

TEST(NativeRunner, CancelMultithread)
{
    std::string sample_path = std::string(PROJECT_DIR) +
                              std::string("/Power-tower-surround_singlefacet.stinput");

    SimulationData sd;
    EXPECT_TRUE(sd.import_from_file(sample_path));
    sd.set_number_of_rays(1000000);

    NativeRunner runner;
    runner.disable_point_focus();
    runner.disable_power_tower();
    runner.set_number_of_threads(4);
    RunnerStatus sts;
    sts = runner.setup_simulation(&sd);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    auto fsts = std::async(&NativeRunner::run_simulation, &runner);

    // Give time to start processing
    std::this_thread::sleep_for(std::chrono::milliseconds(500));
    sts = runner.status_simulation();
    EXPECT_EQ(sts, RunnerStatus::RUNNING);

    // Shut everything down to make sure it doesn't hang
    auto t0 = std::chrono::high_resolution_clock::now();
    runner.cancel_simulation();
    ASSERT_EQ(fsts.wait_for(std::chrono::seconds(10)), std::future_status::ready);
    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> dur = t1 - t0;

    EXPECT_EQ(fsts.get(), RunnerStatus::CANCEL);

    std::cout << "Time to cancel: " << dur.count() << std::endl;
}

TEST(NativeRunner, RayIdAssignmentMultiThread)
{
    const uint_fast64_t NRAYS = 50000;
    const unsigned NTHREADS = 12;

    std::string project_path = std::string(PROJECT_DIR);
    std::string sample_path = project_path +
                              std::string("/simple_test_case.stinput");

    // Load simulation data from file
    SimulationData sd;
    ASSERT_TRUE(sd.import_from_file(sample_path));

    sd.set_number_of_rays(NRAYS);

    // Create and run the native runner
    NativeRunner runner;
    runner.set_number_of_threads(NTHREADS);
    RunnerStatus sts = runner.initialize();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.setup_simulation(&sd);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.run_simulation();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    SimulationResult result;
    sts = runner.report_simulation(&result, 0);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    // std::cout << result << std::endl;

    uint_fast64_t nrec = result.get_number_of_records();

    std::cout << "Number of rays: " << NRAYS
              << "\nNumber of records: " << nrec
              << std::endl;

    if (nrec != NRAYS)
    {
        runner.print_log(std::cout);
    }

    ASSERT_EQ(nrec, NRAYS);

    for (uint_fast64_t k = 0; k < nrec; ++k)
    {
        auto ray_rec = result[k];
        ASSERT_EQ(ray_rec->id, k + 1);
    }
}

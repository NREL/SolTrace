#include <gtest/gtest.h>

#include <limits>
#include <stdexcept>

#include <optix_runner.hpp>
#include <simulation_data.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result_export.hpp>

using SolTrace::Runner::RunnerStatus;

// Reuse the two-plate scene defined in two_plate_test.cpp
void make_two_plate_sd(SimulationData &sd, element_ptr &plate1, element_ptr &plate2);

// ---------------------------------------------------------------------------
// set_batch_size / get_batch_size accessor tests (no GPU required)
// ---------------------------------------------------------------------------

// Default value of 0 means automatic batch sizing: determine_batch_size() will
// call automatic_batch_size() to compute a GPU-memory-safe batch size at run
// time. It does NOT mean "launch all rays in a single batch".
TEST(OptixRunnerBatchSize, DefaultIsZeroMeansAutoSize)
{
    OptixRunner runner;
    EXPECT_EQ(runner.get_batch_size(), 0u);
}

TEST(OptixRunnerBatchSize, SetAndGet)
{
    OptixRunner runner;
    runner.set_batch_size(500);
    EXPECT_EQ(runner.get_batch_size(), 500u);
}

// Setting back to 0 restores automatic GPU-memory-based sizing.
TEST(OptixRunnerBatchSize, SetZeroRestoresAutoSize)
{
    OptixRunner runner;
    runner.set_batch_size(1000);
    runner.set_batch_size(0);
    EXPECT_EQ(runner.get_batch_size(), 0u);
}

TEST(OptixRunnerBatchSize, ThrowsOnOverflow)
{
    OptixRunner runner;
    const uint_fast64_t too_large =
        static_cast<uint_fast64_t>(std::numeric_limits<int>::max()) + 1ULL;
    EXPECT_THROW(runner.set_batch_size(too_large), std::out_of_range);
    // Value should be unchanged after the throw
    EXPECT_EQ(runner.get_batch_size(), 0u);
}

// INT_MAX itself is the largest valid batch size; it must not throw.
TEST(OptixRunnerBatchSize, MaxIntBoundaryDoesNotThrow)
{
    OptixRunner runner;
    const uint_fast64_t max_valid =
        static_cast<uint_fast64_t>(std::numeric_limits<int>::max());
    EXPECT_NO_THROW(runner.set_batch_size(max_valid));
    EXPECT_EQ(runner.get_batch_size(), max_valid);
}

// A failed set_batch_size must not corrupt a previously stored non-zero value.
TEST(OptixRunnerBatchSize, ThrowPreservesExistingValue)
{
    OptixRunner runner;
    runner.set_batch_size(999);
    const uint_fast64_t too_large =
        static_cast<uint_fast64_t>(std::numeric_limits<int>::max()) + 1ULL;
    EXPECT_THROW(runner.set_batch_size(too_large), std::out_of_range);
    EXPECT_EQ(runner.get_batch_size(), 999u);
}

// ---------------------------------------------------------------------------
// Simulation correctness: batched run should yield the same hit count as the
// default single-batch run.
// ---------------------------------------------------------------------------

TEST(OptixRunnerBatchSize, BatchedRunMatchesHitCount)
{
    const int N_rays = 10000;

    // --- reference run (default auto-sized batch: m_batch_size == 0 defers to
    //     determine_batch_size() / automatic_batch_size()) ---
    SimulationData sd_ref;
    element_ptr p1_ref, p2_ref;
    make_two_plate_sd(sd_ref, p1_ref, p2_ref);
    sd_ref.get_simulation_parameters().number_of_rays = N_rays;
    sd_ref.get_simulation_parameters().max_number_of_rays = N_rays * 100;

    OptixRunner ref_runner;
    ASSERT_EQ(ref_runner.initialize(), RunnerStatus::SUCCESS);
    ASSERT_EQ(ref_runner.setup_simulation(&sd_ref), RunnerStatus::SUCCESS);
    ASSERT_EQ(ref_runner.run_simulation(), RunnerStatus::SUCCESS);

    RayHistoryResult ref_result;
    ASSERT_EQ(ref_runner.report_simulation(&ref_result, 0), RunnerStatus::SUCCESS);
    const int ref_hits = ref_result.get_number_of_records();

    // --- explicit batched run (user-supplied batch_size = 1000, ~10 iterations) ---
    SimulationData sd_batch;
    element_ptr p1_batch, p2_batch;
    make_two_plate_sd(sd_batch, p1_batch, p2_batch);
    sd_batch.get_simulation_parameters().number_of_rays = N_rays;
    sd_batch.get_simulation_parameters().max_number_of_rays = N_rays * 100;

    OptixRunner batch_runner;
    batch_runner.set_batch_size(1000);
    ASSERT_EQ(batch_runner.initialize(), RunnerStatus::SUCCESS);
    ASSERT_EQ(batch_runner.setup_simulation(&sd_batch), RunnerStatus::SUCCESS);
    ASSERT_EQ(batch_runner.run_simulation(), RunnerStatus::SUCCESS);

    RayHistoryResult batch_result;
    ASSERT_EQ(batch_runner.report_simulation(&batch_result, 0), RunnerStatus::SUCCESS);
    const int batch_hits = batch_result.get_number_of_records();

    EXPECT_EQ(ref_hits, N_rays);
    EXPECT_EQ(batch_hits, N_rays);

    // With a user-supplied batch of 1000 and 10000 rays, at least 10 iterations
    // are required regardless of available GPU memory.
    EXPECT_GE(batch_runner.get_N_run_iterations(), 10u);
}

// ---------------------------------------------------------------------------
// Batch size smaller than total rays forces multiple iterations.
// ---------------------------------------------------------------------------

TEST(OptixRunnerBatchSize, SmallBatchMultipleIterations)
{
    const int N_rays = 5000;
    const int batch = 500; // 10+ iterations needed

    SimulationData sd;
    element_ptr plate1, plate2;
    make_two_plate_sd(sd, plate1, plate2);
    sd.get_simulation_parameters().number_of_rays = N_rays;
    sd.get_simulation_parameters().max_number_of_rays = N_rays * 100;

    OptixRunner runner;
    runner.set_batch_size(batch);
    ASSERT_EQ(runner.initialize(), RunnerStatus::SUCCESS);
    ASSERT_EQ(runner.setup_simulation(&sd), RunnerStatus::SUCCESS);
    ASSERT_EQ(runner.run_simulation(), RunnerStatus::SUCCESS);

    RayHistoryResult result;
    ASSERT_EQ(runner.report_simulation(&result, 0), RunnerStatus::SUCCESS);

    EXPECT_EQ(result.get_number_of_records(), N_rays);

    // With a small batch the runner must have generated at least N_rays sun rays
    EXPECT_GE(runner.get_N_sun_rays(), static_cast<uint_fast64_t>(N_rays));

    // Multiple iterations must have been needed
    EXPECT_GE(runner.get_N_run_iterations(), 10u);
}

// ---------------------------------------------------------------------------
// Batch size >= N_rays: should complete in exactly one iteration.
// ---------------------------------------------------------------------------

TEST(OptixRunnerBatchSize, BatchSizeExceedingRaysCompletesInOneIteration)
{
    const int N_rays = 1000;

    SimulationData sd;
    element_ptr plate1, plate2;
    make_two_plate_sd(sd, plate1, plate2);
    sd.get_simulation_parameters().number_of_rays = N_rays;
    sd.get_simulation_parameters().max_number_of_rays = N_rays * 100;

    OptixRunner runner;
    runner.set_batch_size(N_rays * 8); // larger than N_rays
    ASSERT_EQ(runner.initialize(), RunnerStatus::SUCCESS);
    ASSERT_EQ(runner.setup_simulation(&sd), RunnerStatus::SUCCESS);
    ASSERT_EQ(runner.run_simulation(), RunnerStatus::SUCCESS);

    RayHistoryResult result;
    ASSERT_EQ(runner.report_simulation(&result, 0), RunnerStatus::SUCCESS);

    EXPECT_EQ(result.get_number_of_records(), N_rays);
    EXPECT_EQ(runner.get_N_run_iterations(), 1u);
}

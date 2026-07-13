#include <gtest/gtest.h>

#include <chrono>
#include <future>
#include <mutex>
#include <set>
#include <sstream>
#include <thread>
#include <tuple>
#include <vector>

#include <generate_ray.hpp>
#include <native_runner.hpp>
#include <simulation_data_export.hpp>
#include <simulation_runner.hpp>
#include <simulation_result_export.hpp>

using SolTrace::Runner::RunnerStatus;
using SolTrace::NativeRunner::GenerateRay;
using SolTrace::NativeRunner::MTRand;
using SolTrace::NativeRunner::NativeRunner;
using SolTrace::NativeRunner::TSun;
using SolTrace::Result::RayEvent;

namespace
{
    SimulationData make_halton_test_simulation(uint_fast64_t nrays)
    {
        SimulationData sd;

        auto sun = SolTrace::Data::make_ray_source<Sun>();
        sun->set_position(0.0, 0.0, 100.0);
        sun->set_shape(SolTrace::Data::SunShape::GAUSSIAN, 1.0, -5.0, 0.0);
        sun->set_gen_type(SolTrace::Data::GenType::HALTON);
        sd.add_ray_source(sun);

        auto stage = SolTrace::Data::make_stage(0);

        auto mirror = SolTrace::Data::make_element<SingleElement>();
        mirror->set_aperture(SolTrace::Data::make_aperture<Rectangle>(20.0, 20.0));
        mirror->set_surface(SolTrace::Data::make_surface<Flat>());
        mirror->set_reference_frame_geometry(glm::dvec3(0.0, 0.0, 0.0),
            glm::dvec3(0.0, 0.0, 1.0),
            0.0);

        SolTrace::Data::OpticalPropertySet mirror_optics_set(SolTrace::Data::InteractionType::REFLECTION, "Mirror");
        mirror_optics_set.set_ideal_reflection(SolTrace::Data::OpticalSide::Both);
        mirror_optics_set.set_errors(SolTrace::Data::OpticalSide::Both, SolTrace::Data::DistributionType::NONE, 0.0, 0.0);
        auto ref = sd.add_optical_property_set(mirror_optics_set);
        mirror->set_optical_property_set(ref);

        stage->add_element(mirror);

        sd.add_stage(stage);

        auto &params = sd.get_simulation_parameters();
        params.number_of_rays = nrays;
        params.max_number_of_rays = nrays * 10;
        params.include_optical_errors = false;
        params.include_sun_shape_errors = false;
        params.seed = 12345;

        return sd;
    }

    void make_identity(double R[3][3])
    {
        for (int i = 0; i < 3; ++i)
        {
            for (int j = 0; j < 3; ++j)
            {
                R[i][j] = (i == j ? 1.0 : 0.0);
            }
        }
    }
}

TEST(GenerateRay, HaltonCommonSenseSingleThread)
{
    const uint_fast64_t NRAYS = 8;

    SimulationData sd = make_halton_test_simulation(NRAYS);

    NativeRunner runner;
    runner.set_number_of_threads(1);

    RunnerStatus sts = runner.initialize();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.setup_simulation(&sd);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.run_simulation();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    RayHistoryResult result;
    sts = runner.report_simulation(&result, 0);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    ASSERT_EQ(result.get_number_of_records(), NRAYS);

    std::set<std::tuple<double, double, double>> create_points;
    for (int_fast64_t k = 0; k < result.get_number_of_records(); ++k)
    {
        auto rr = result[k];
        ASSERT_GT(rr->get_number_of_interactions(), 0);
        EXPECT_EQ(rr->get_event(0), RayEvent::CREATE);

        glm::dvec3 pos;
        rr->get_position(0, pos);
        EXPECT_DOUBLE_EQ(pos[2], 10000.0);
        create_points.emplace(pos[0], pos[1], pos[2]);
    }

    EXPECT_EQ(create_points.size(), NRAYS);
}

TEST(GenerateRay, HaltonUniqueMultiThread)
{
    const uint_fast64_t NRAYS = 64;
    const unsigned NTHREADS = 8;
    const uint_fast64_t rem = NRAYS % NTHREADS;
    const uint_fast64_t nrays_per_thread = NRAYS / NTHREADS;

    TSun sun;
    sun.PointSource = false;
    sun.GenTypeIndex = SolTrace::Data::GenType::HALTON;
    sun.MinXSun = -10.0;
    sun.MaxXSun = 10.0;
    sun.MinYSun = -10.0;
    sun.MaxYSun = 10.0;
    sun.RLocToRef = glm::dmat3(1.0);    // Set identity

    const glm::dvec3 PosSunStage = {0.0, 0.0, 0.0};
    glm::dvec3 Origin = {0.0, 0.0, 0.0};
    glm::dmat3 RLocToRef(1.0);

    struct ThreadResult
    {
        uint_fast64_t expected_count;
        std::set<std::tuple<double, double, double>> points;
    };

    std::set<std::tuple<double, double, double>> all_points;
    std::vector<std::future<ThreadResult>> futures;

    for (unsigned thread_id = 0; thread_id < NTHREADS; ++thread_id)
    {
        futures.push_back(std::async(std::launch::async, [&, thread_id]() {
            MTRand myrng(12345 + thread_id);
            const uint_fast64_t local_count = thread_id < rem ? nrays_per_thread + 1 : nrays_per_thread;
            const uint_fast64_t ray_index_offset = thread_id * nrays_per_thread + std::min(static_cast<uint_fast64_t>(thread_id), rem);

            ThreadResult result;
            result.expected_count = local_count;
            for (uint_fast64_t k = 0; k < local_count; ++k)
            {
                const uint_fast64_t sample_index = ray_index_offset + k + 1;
                glm::dvec3 PosRayGlobal = {0.0, 0.0, 0.0};
                glm::dvec3 CosRayGlobal = {0.0, 0.0, 0.0};
                glm::dvec3 PosRaySun = {0.0, 0.0, 0.0};
                int ErrorFlag = 0;

                GenerateRay(myrng, PosSunStage, Origin, RLocToRef, &sun,
                            sample_index, PosRayGlobal, CosRayGlobal, PosRaySun,
                            ErrorFlag);

                result.points.emplace(PosRaySun[0], PosRaySun[1], PosRaySun[2]);
            }

            return result;
        }));
    }

    for (auto &f : futures)
    {
        const ThreadResult result = f.get();
        EXPECT_EQ(result.points.size(), result.expected_count);
        all_points.insert(result.points.begin(), result.points.end());
    }

    EXPECT_EQ(all_points.size(), NRAYS);
}

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

    double prog = -1.0;
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

    RayHistoryResult result;
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

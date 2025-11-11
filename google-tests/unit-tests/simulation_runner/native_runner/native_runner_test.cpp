#include <gtest/gtest.h>

#include <chrono>
#include <future>
#include <memory>
#include <string>
#include <thread>

#include <constants.hpp>
#include <error_distributions.hpp>
#include <mtrand.hpp>
#include <native_runner.hpp>
#include <native_runner_types.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result_export.hpp>
#include <simulation_runner.hpp>

#include "common.hpp"
#include "count_absorbed_native.h"

using SolTrace::Runner::RunnerStatus;

using SolTrace::NativeRunner::MTRand;
using SolTrace::NativeRunner::NativeRunner;
using SolTrace::NativeRunner::TRayData;
using SolTrace::NativeRunner::TSystem;

TEST(RandomNumberGenerator, SingleNumberMersenneTwister)
{
    MTRand myrng(1);
    double random_number = myrng.rand();

    EXPECT_NEAR(random_number, 0.13387664401253274, 1e-7);
}

TEST(NativeRunnerTypes, TSun)
{

    SimulationData my_sim;
    auto sun = SolTrace::Data::make_ray_source<Sun>();
    Vector3d spos(1.0, 2.0, 3.0);
    sun->set_position(spos);
    sun->set_shape(SunShape::PILLBOX, -1.0, 1.0, 0.0);
    my_sim.add_ray_source(sun);

    NativeRunner runner;
    runner.setup_sun(&my_sim);
    auto sys = runner.get_system();
    EXPECT_TRUE(is_identical(sys->Sun.Origin, sun->get_position()));
    EXPECT_EQ(sys->Sun.ShapeIndex, SunShape::PILLBOX);
}

TEST(NativeRunnerTypes, MakeElement)
{
}

TEST(NativeRunnerTypes, MakeStage)
{
}

// TEST(NativeRunnerTypes, TElement)
// {
//     // TODO: Implement a test here...

//     // SimulationData my_sim;
//     // // **** Setup Answers **** //
//     // // Origin
//     // Vector3d Origin1(1.0, 2.0, 3.0);
//     // // Corresponding Euler angles in radians
//     // const double a1 = 0.0;
//     // const double b1 = asin(-1.0 / sqrt(3.0));
//     // const double g1 = acos(1.0 / cos(b1) * 1.0 / sqrt(6.0)); // approximately 0.615
//     // // Corresponding aim vector (local z-axis in reference coordinates)
//     // Vector3d aim1(0.0, -1.0 / sqrt(3.0), sqrt(2.0 / 3.0));
//     // vector_add(1.0, Origin1, 1.0, aim1);

//     // // Z-Rotation is the last of the Euler angles but in degrees
//     // const double zrot1 = g1 * 180.0 / PI;

//     // // Origin
//     // Vector3d Origin2(-3.0, 1.0, -5.0);
//     // const double a2 = PI / 4.0;
//     // const double b2 = PI / 6.0;
//     // const double g2 = PI / 3.0;
//     // // Corresponding aim vector (local z-axis in reference coordinates)
//     // Vector3d aim2(sqrt(3.0 / 8.0), 0.5, sqrt(3.0 / 8.0));
//     // vector_add(1.0, Origin2, 1.0, aim2);

//     // // Z-Rotation is the last of the Euler angles but in degrees
//     // const double zrot2 = 60.0;

//     // // **** Setup Elements **** //
//     // auto el = make_element<SingleElement>();
//     // el->set_aperture(make_aperture<Circle>(2.0));
//     // el->set_surface(make_surface<Flat>());
//     // el->set_reference_frame_geometry(Origin1, aim1, zrot1);

//     // auto st = make_stage(0);
//     // st->set_reference_frame_geometry(Origin2, aim2, zrot2);
//     // st->add_element(el);
// }

// TEST(NativeRunnerTypes, TStage)
// {
//     // TODO: Implement test
// }

TEST(NativeRunner, SmokeTest)
{
    const unsigned NRAYS = 10;
    NativeRunner runner;
    SimulationData my_sim;

    SimulationParameters &params = my_sim.get_simulation_parameters();
    params.include_optical_errors = false;
    params.include_sun_shape_errors = false;
    params.number_of_rays = NRAYS;
    params.max_number_of_rays = 10 * NRAYS;
    // my_sim.set_number_of_rays(10);
    // my_sim.set_max_rays_traced(100);

    auto sun = SolTrace::Data::make_ray_source<Sun>();
    sun->set_position(0.0, 0.0, 100.0);
    sun->set_shape(SolTrace::Data::SunShape::GAUSSIAN, 1.0, -5.0, 0.0);
    my_sim.add_ray_source(sun);

    auto my_st = SolTrace::Data::make_stage(0);
    const int NUM_ELEMENTS = 4;
    double x[NUM_ELEMENTS] = {1.0, 0.0, -1.0, 0.0};
    double y[NUM_ELEMENTS] = {0.0, 1.0, 0.0, -1.0};
    OpticalProperties optics(SolTrace::Data::InteractionType::REFLECTION,
                             SolTrace::Data::DistributionType::GAUSSIAN,
                             0.0, 1.0, 0.0, 0.0, 1.0, 1.0);

    for (int k = 0; k < NUM_ELEMENTS; ++k)
    {
        element_ptr el = SolTrace::Data::make_element<SingleElement>();
        el->set_aperture(SolTrace::Data::make_aperture<Circle>(2.0));
        el->set_surface(SolTrace::Data::make_surface<Flat>());
        el->set_reference_frame_geometry(Vector3d(x[k], y[k], 0.0),
                                         Vector3d(-x[k], -y[k], 1.0),
                                         0.0);
        el->set_front_optical_properties(optics);
        el->set_back_optical_properties(optics);
        my_st->add_element(el);
    }

    EXPECT_EQ(my_st->get_number_of_elements(), NUM_ELEMENTS);
    my_sim.add_stage(my_st);
    EXPECT_EQ(my_sim.get_number_of_elements(), NUM_ELEMENTS);

    RunnerStatus sts;
    sts = runner.initialize();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.setup_simulation(&my_sim);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.run_simulation();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);

    const TSystem *sys = runner.get_system();
    // sys->AllRayData.Print();
    // const TRayData *ray_data = &(sys->AllRayData);
    sys->RayData.Print();
    const TRayData *ray_data = &(sys->RayData);
    uint_fast64_t num_absorbed = count_absorbed_native(ray_data);
    uint_fast64_t ncreate = count_event_native(ray_data, RayEvent::CREATE);
    uint_fast64_t nexit = count_event_native(ray_data, RayEvent::EXIT);
    size_t n = ray_data->Count();

    std::cout << "Number Created: " << ncreate << std::endl;
    std::cout << "Number Absorbed: " << num_absorbed << std::endl;
    std::cout << "Number Interactions: " << n << std::endl;
    std::cout << "Number Exit: " << nexit << std::endl;

    EXPECT_EQ(ncreate, NRAYS);
    EXPECT_EQ(num_absorbed + nexit, ncreate);
    EXPECT_EQ(num_absorbed + nexit, NRAYS);

    ray_data->Print();

    EXPECT_TRUE(n >= NRAYS);

    SimulationResult result;
    sts = runner.report_simulation(&result, 0);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);

    EXPECT_EQ(result.get_number_of_records(), NRAYS);
    for (int_fast64_t idx = 0; idx < result.get_number_of_records(); ++idx)
    {
        const ray_record_ptr rr = result[idx];
        EXPECT_TRUE(rr->get_number_of_interactions() > 0);
    }

    std::cout << "Number of ray records: "
              << result.get_number_of_records()
              << std::endl;
}

TEST(NativeRunner, PowerTowerSmokeTest)
{
    SimulationData sd;

    // Sun
    auto sun = SolTrace::Data::make_ray_source<Sun>();
    sun->set_position(0.0, 0.0, 100.0);
    sd.add_ray_source(sun);

    // Absorber -- Flat
    auto absorber = SolTrace::Data::make_element<SingleElement>();
    absorber->set_origin(0.0, 0.0, 10.0);
    absorber->set_aim_vector(0.0, 0.0, 5.0);
    absorber->set_zrot(0.0);
    absorber->compute_coordinate_rotations();
    absorber->set_surface(SolTrace::Data::make_surface<Flat>()); // surface(nullptr)
    absorber->set_aperture(SolTrace::Data::make_aperture<Rectangle>(2.0, 2.0));
    OpticalProperties *foptics = absorber->get_front_optical_properties();
    foptics->my_type = InteractionType::REFLECTION;
    foptics->reflectivity = 0.0;

    // Make stage 1 -- second stage -- these can be added to SimulationData
    // in any order but should be numbered in the desired order
    auto st1 = SolTrace::Data::make_stage(1);
    // Origin is initialized to zero but set it explicitly
    st1->set_origin(0.0, 0.0, 0.0);
    // Set aim vector so stage and global coordinates are identical
    st1->set_aim_vector(0.0, 0.0, 1.0);
    // Set z rotation so stage and global coordinates are identical
    st1->set_zrot(0.0);
    // Compute coordinate rotations
    st1->compute_coordinate_rotations();
    st1->add_element(absorber);
    // Optional -- to help the user identify things

    // Make stage 0 -- this will be the first stage if the runner uses stages
    auto st0 = SolTrace::Data::make_stage(0);
    st0->set_origin(0.0, 0.0, 0.0);
    st0->set_aim_vector(0.0, 0.0, 1.0);

    Vector3d rvec, svec, avec;
    Vector3d aim, pos;

    const int NUM_ELEMENTS = 10;
    for (int k = 0; k < NUM_ELEMENTS; ++k)
    {
        auto el = SolTrace::Data::make_element<SingleElement>();
        foptics = el->get_front_optical_properties();
        foptics->reflectivity = 1.0;

        pos.set_values(5 * sin(k * PI * 2.0 / NUM_ELEMENTS),
                       5 * cos(k * PI * 2.0 / NUM_ELEMENTS),
                       0.0);
        vector_add(1.0, absorber->get_origin_global(),
                   -1.0, pos,
                   rvec);
        make_unit_vector(rvec);
        svec = sun->get_position();
        make_unit_vector(svec);
        vector_add(0.5, rvec, 0.5, svec, avec);
        vector_add(1.0, pos, 100.0, avec, aim);

        el->set_reference_frame_geometry(pos, aim, 0.0);

        el->set_surface(SolTrace::Data::make_surface<Flat>());
        el->set_aperture(SolTrace::Data::make_aperture<Circle>(2.0));

        st0->add_element(el);
    }
    EXPECT_EQ(st0->get_number_of_elements(), NUM_ELEMENTS);
    sd.add_stage(st0);
    EXPECT_EQ(sd.get_number_of_elements(), NUM_ELEMENTS);
    sd.add_stage(st1);
    EXPECT_EQ(sd.get_number_of_elements(), NUM_ELEMENTS + 1);

    // Set parameters
    const uint_fast64_t NRAYS = 10000;
    SimulationParameters &params = sd.get_simulation_parameters();
    params.number_of_rays = NRAYS;
    params.max_number_of_rays = params.number_of_rays * 100;
    params.include_optical_errors = false;
    params.include_sun_shape_errors = false;
    params.seed = 12345;

    NativeRunner runner;
    RunnerStatus sts = runner.initialize();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.setup_simulation(&sd);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.run_simulation();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);

    const TSystem *sys = runner.get_system();
    // sys->RayData.Print();
    const TRayData *ray_data = &(sys->RayData);
    uint_fast64_t num_absorbed = count_absorbed_native(ray_data);
    uint_fast64_t ncreate = count_event_native(ray_data, RayEvent::CREATE);
    uint_fast64_t nexit = count_event_native(ray_data, RayEvent::EXIT);
    size_t n = ray_data->Count();

    std::cout << "Number Created: " << ncreate << std::endl;
    std::cout << "Number Absorbed: " << num_absorbed << std::endl;
    std::cout << "Number Interactions: " << n << std::endl;
    std::cout << "Number Exit: " << nexit << std::endl;

    scan_events_native(ray_data);

    EXPECT_EQ(ncreate, NRAYS);
    EXPECT_EQ(num_absorbed + nexit, ncreate);
    EXPECT_EQ(num_absorbed + nexit, NRAYS);

    // ray_data->Print();

    EXPECT_TRUE(n >= NRAYS);
    EXPECT_TRUE(num_absorbed > 0);

    SimulationResult result;
    sts = runner.report_simulation(&result, 0);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);

    EXPECT_EQ(result.get_number_of_records(), NRAYS);
    for (int_fast64_t idx = 0; idx < result.get_number_of_records(); ++idx)
    {
        const ray_record_ptr rr = result[idx];
        EXPECT_TRUE(rr->get_number_of_interactions() > 0);
    }

    std::cout << "Number of ray records: "
              << result.get_number_of_records()
              << std::endl;
}

TEST(NativeRunner, SingleRayValidationTest)
{
    const double TOL = 1e-8;

    constexpr double c = 0.09;
    constexpr double R = 1.0 / c;
    constexpr double x = -3.0621423346154577;
    constexpr double y = 5.9286205128611948;
    constexpr double z0 = 15.0;

    double zref = R - sqrt(R * R - (x * x + y * y));
    double z = z0 - zref;

    Vector3d u(0.0, 0.0, -1.0);
    Vector3d v(2.0 * x, 2.0 * y, -2.0 * (z0 - z - R));
    v.make_unit();
    Vector3d w;
    double alpha = dot_product(u, v);
    vector_add(1.0, u, -2.0 * alpha, v, w);

    SimulationData sd;

    // Set parameters
    const uint_fast64_t NRAYS = 1;
    SimulationParameters &params = sd.get_simulation_parameters();
    params.number_of_rays = NRAYS;
    params.max_number_of_rays = params.number_of_rays * 100;
    params.include_optical_errors = false;
    params.include_sun_shape_errors = false;
    params.seed = 1;

    // Sun
    auto sun = SolTrace::Data::make_ray_source<Sun>();
    sun->set_position(0.0, 0.0, 100.0);
    sun->set_shape(SolTrace::Data::SunShape::PILLBOX, -1.0, 1.0, 0.0);
    sd.add_ray_source(sun);

    auto sph = SolTrace::Data::make_element<SingleElement>();
    Vector3d origin(0.0, 0.0, z0);
    Vector3d aim(0.0, 0.0, -1.0);
    double zrot = 0.0;
    sph->set_reference_frame_geometry(origin, aim, zrot);
    sph->set_aperture(SolTrace::Data::make_aperture<Hexagon>(20.0));
    sph->set_surface(SolTrace::Data::make_surface<Sphere>(c));
    sph->get_front_optical_properties()->set_ideal_reflection();
    sph->get_back_optical_properties()->set_ideal_reflection();
    sph->set_name("Sphere");
    sd.add_element(sph);

    auto para = SolTrace::Data::make_element<SingleElement>();
    origin.set_values(0.0, 0.0, -1.0);
    aim.set_values(0.0, 0.0, 0.0);
    zrot = 0.0;
    para->set_reference_frame_geometry(origin, aim, zrot);
    para->set_aperture(SolTrace::Data::make_aperture<Rectangle>(31.0, 31.0));
    para->set_surface(SolTrace::Data::make_surface<Parabola>(0.5 / 0.03, 0.5 / 0.03));
    para->set_name("Parabola");
    sd.add_element(para);

    NativeRunner runner;
    RunnerStatus sts = runner.initialize();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.setup_simulation(&sd);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.run_simulation();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);

    const TSystem *sys = runner.get_system();
    sys->RayData.Print();
    const TRayData *ray_data = &(sys->RayData);
    size_t n = ray_data->Count();

    // sys->RayData.Print();

    EXPECT_EQ(n, 3);

    Vector3d ipoint, idir;
    int element, stage;
    unsigned int raynum;
    SolTrace::Result::RayEvent rev;
    sys->RayData.Query(0, ipoint.data, idir.data,
                       &element, &stage, &raynum, &rev);

    EXPECT_EQ(raynum, 1);
    EXPECT_EQ(rev, SolTrace::Result::RayEvent::CREATE);

    EXPECT_NEAR(ipoint[0], x, TOL);
    EXPECT_NEAR(ipoint[1], y, TOL);
    EXPECT_NEAR(ipoint[2], 10000.0, TOL);

    EXPECT_NEAR(idir[0], u[0], TOL);
    EXPECT_NEAR(idir[1], u[1], TOL);
    EXPECT_NEAR(idir[2], u[2], TOL);

    sys->RayData.Query(1, ipoint.data, idir.data,
                       &element, &stage, &raynum, &rev);

    EXPECT_EQ(raynum, 1);
    EXPECT_EQ(rev, SolTrace::Result::RayEvent::REFLECT);

    EXPECT_NEAR(ipoint[0], x, TOL);
    EXPECT_NEAR(ipoint[1], y, TOL);
    EXPECT_NEAR(ipoint[2], z, TOL);

    EXPECT_NEAR(idir[0], w[0], TOL);
    EXPECT_NEAR(idir[1], w[1], TOL);
    EXPECT_NEAR(idir[2], w[2], TOL);
}

TEST(NativeRunner, LegacyFileLoadTest)
{
    std::string project_path = std::string(PROJECT_DIR);
    std::string sample_path = project_path +
                              std::string("/simple_test_case.stinput");

    // Load simulation data from file
    SimulationData sd;
    bool success = sd.import_from_file(sample_path);
    EXPECT_TRUE(success);

    EXPECT_EQ(sd.get_number_of_ray_sources(), 1);
    EXPECT_EQ(sd.get_number_of_elements(), 2);

    // Create and run the native runner
    NativeRunner runner;
    RunnerStatus sts = runner.initialize();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.setup_simulation(&sd);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.run_simulation();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
}

TEST(NativeRunner, StatusAndCancel)
{
    std::string sample_path = std::string(PROJECT_DIR) +
                              std::string("/Power-tower-surround_singlefacet.stinput");

    SimulationData sd;
    EXPECT_TRUE(sd.import_from_file(sample_path));
    sd.set_number_of_rays(50000);

    NativeRunner runner;
    runner.disable_point_focus();
    runner.disable_power_tower();
    RunnerStatus sts;
    sts = runner.setup_simulation(&sd);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);

    auto t0 = std::chrono::high_resolution_clock::now();

    auto fsts = std::async(&NativeRunner::run_simulation, &runner);

    std::this_thread::sleep_for(std::chrono::milliseconds(200));
    sts = runner.status_simulation();
    EXPECT_EQ(sts, RunnerStatus::RUNNING);

    double prog;
    std::this_thread::sleep_for(std::chrono::milliseconds(200));
    sts = runner.status_simulation(&prog);
    EXPECT_EQ(sts, RunnerStatus::RUNNING);
    EXPECT_LE(prog, 1.0);
    EXPECT_GE(prog, 0.0);

    runner.cancel_simulation();
    fsts.wait();

    auto t1 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double, std::milli> dur = t1 - t0;

    EXPECT_EQ(fsts.get(), RunnerStatus::CANCEL);
    EXPECT_LT(dur.count(), 2000.0);

    std::cout << "Time for run: " << dur.count() << std::endl;
    std::cout << "Progress before cancel: " << prog << std::endl;
}

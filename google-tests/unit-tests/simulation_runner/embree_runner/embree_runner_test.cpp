#include <gtest/gtest.h>

#include <embree_runner.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result.hpp>

#include "common.hpp"
#include "count_absorbed_native.h"

using SolTrace::EmbreeRunner::EmbreeRunner;
using SolTrace::EmbreeRunner::RunnerStatus;
using SolTrace::EmbreeRunner::TSystem;
using SolTrace::EmbreeRunner::TRayData;

using SolTrace::Result::ray_record_ptr;
using SolTrace::Result::RayEvent;
using SolTrace::Result::SimulationResult;

TEST(EmbreeRunner, SingleRayValidationTest)
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

    // std::cout << "Constructing..." << std::endl;
    EmbreeRunner runner;
    runner.set_number_of_threads(1);
    // std::cout << "Initializing..." << std::endl;
    RunnerStatus sts = runner.initialize();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    // std::cout << "Setting up simulation..." << std::endl;
    sts = runner.setup_simulation(&sd);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    // std::cout << "Running simulation..." << std::endl;
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
    uint_fast64_t raynum;
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

TEST(EmbreeRunner, PowerTowerSmokeTest)
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

    EmbreeRunner runner;
    runner.set_number_of_threads(1);
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

TEST(EmbreeRunner, PowerTowerTest)
{
    // Pulling in path variable from CMake and creating path to .stinput sample file
    std::string sample_path = std::string(PROJECT_DIR) + std::string("/Power-tower-surround_singlefacet.stinput");

    // Create Simuluation Data
    SimulationData sd;

    // Constants
    const uint_fast64_t NRAYS = 10000;
    const double TOL = 1e-4;

    // Read Input File
    bool success = sd.import_from_file(sample_path);
    EXPECT_TRUE(success);
    EXPECT_TRUE(sd.get_number_of_elements() > 0);
    EXPECT_TRUE(sd.get_number_of_ray_sources() > 0);

    std::cout << "Num Elements: " << sd.get_number_of_elements() << std::endl;

    // Parameters
    SimulationParameters &params = sd.get_simulation_parameters();
    params.include_optical_errors = true;
    params.include_sun_shape_errors = true;
    params.max_number_of_rays = NRAYS * 100;
    params.number_of_rays = NRAYS;
    params.seed = 1;

    // Run Ray Trace
    EmbreeRunner runner;
    runner.set_number_of_threads(1);
    RunnerStatus sts;
    sts = runner.initialize();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.setup_simulation(&sd);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);

    auto t0 = std::chrono::high_resolution_clock::now();
    sts = runner.run_simulation();
    auto t1 = std::chrono::high_resolution_clock::now();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);

    std::chrono::duration<double, std::milli> dur = t1 - t0;

    std::cout << "Time: " << dur.count() << " ms" << std::endl;

    SimulationResult result;
    sts = runner.report_simulation(&result, 0);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    EXPECT_EQ(result.get_number_of_records(), NRAYS);

    element_id absorber_id = 6285;
    int_fast64_t nabsorbed = count_element_event(result, absorber_id, RayEvent::ABSORB);
    int_fast64_t nreflect = count_element_event(result, absorber_id, RayEvent::REFLECT);
    int_fast64_t nevents = nabsorbed + nreflect;

    std::cout << "Total: " << nevents
              << "\nAbsorb: " << nabsorbed << " ("
              << static_cast<double>(nabsorbed) / nevents << ")"
              << "\nReflect: " << nreflect << " ("
              << static_cast<double>(nreflect) / nevents << ")"
              << std::endl;
}

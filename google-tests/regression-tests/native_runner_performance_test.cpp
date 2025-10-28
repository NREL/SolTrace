#include <gtest/gtest.h>

#include <chrono>

#include <aperture.hpp>
#include <element.hpp>
#include <heliostat.hpp>
#include <native_runner.hpp>
#include <simulation_data.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result.hpp>
#include <stage_element.hpp>
#include <sun.hpp>
#include <surface.hpp>

// #include "split_csv.h"
#include "count_absorbed_native.h"

using SolTrace::Runner::RunnerStatus;

using SolTrace::NativeRunner::NativeRunner;
using SolTrace::NativeRunner::TRayData;
using SolTrace::NativeRunner::TSystem;

TEST(NativeRunner, PerformanceTest)
{
    const uint_fast64_t NRAYS = 100000;
    const Vector3d zero(0.0, 0.0, 0.0);
    const Vector3d khat(0.0, 0.0, 1.0);

    const uint_fast64_t NX = 4;
    const uint_fast64_t NY = 2;
    const uint_fast64_t NHX = 20;
    const uint_fast64_t NHY = 12;

    const double LX = 10.0;
    const double LY = 5.0;
    const double dx = 2.0 * LX / NHX;
    const double dy = 2.0 * LY / NHY;
    const double ABS_RADIUS = 10.0;

    std::cout << "Num Elements: " << NX * NY * NHX * NHY << std::endl;
    std::cout << "dx = " << dx << "\ndy = " << dy << std::endl;

    SimulationData sdata;
    NativeRunner my_runner;

    SimulationParameters &params = sdata.get_simulation_parameters();
    params.include_optical_errors = false;
    params.include_sun_shape_errors = false;
    params.max_number_of_rays = NRAYS * 100;
    params.number_of_rays = NRAYS;
    params.seed = 12345;

    my_runner.enable_power_tower();
    my_runner.enable_point_focus();

    OpticalProperties mirror;
    mirror.set_ideal_reflection();

    stage_ptr st1 = make_stage(1);
    st1->set_reference_frame_geometry(zero, khat, 0.0);
    stage_ptr st2 = make_stage(2);
    st2->set_reference_frame_geometry(zero, khat, 0.0);

    Vector3d sun_pos(0.0, 0.0, 10000.0);
    Vector3d abs_origin(0.0, 0.0, 10.0);
    Vector3d hs_origin;
    Vector3d v1;
    Vector3d v2;
    Vector3d aim;
    Vector3d aim_point;

    // double xpos = -1.0 * LX;
    // double ypos = -1.0 * LY;
    double xpos, ypos;

    for (auto ix = 0; ix < NHX; ++ix)
    {
        xpos = dx * ix - LX;
        for (auto jy = 0; jy < NHY; ++jy)
        {
            ypos = dy * jy - LY;
            hs_origin.set_values(xpos, ypos, 0.0);

            // vector_add(1.0, sun_pos, -1.0, hs_origin, v1);
            // vector_add(1.0, abs_origin, -1.0, hs_origin, v2);
            // vector_add(0.5, v1, 0.5, v2, aim);
            // std::cout << aim << std::endl;
            // vector_add(1.0, hs_origin, 1.0, aim, aim_point);

            auto hs = make_element<Heliostat>();
            hs->set_optics(mirror);
            // hs->set_reference_frame_geometry(hs_origin, aim, 0.0);
            hs->set_origin(hs_origin);
            hs->set_aperture_size(2.0 * dx, 2.0 * dy);
            hs->set_number_panels(NX, NY);
            hs->set_gaps(0.0, 0.0);
            hs->set_focal_length(0.0);
            hs->set_canting(Heliostat::NONE, 0.0, 0.0);
            hs->set_target_position(abs_origin);
            hs->set_name("Heliostat");
            hs->enable();

            hs->create_geometry();
            auto ret = st1->add_element(hs);
            EXPECT_TRUE(SolTrace::Data::Element::is_success(ret));
            hs->update_geometry(0.0, 90.0);

            // std::cout << "****************"
            //           << "\nix = " << ix << "  jy = " << jy
            //           << "\norigin: " << hs_origin
            //           << "\naim: " << aim
            //           << "\naim_point: " << aim_point
            //           << "\nv1: " << v1
            //           << "\nv2: " << v2
            //           << std::endl;
        }
    }

    auto absorb = make_element<SingleElement>();
    absorb->get_front_optical_properties()->set_ideal_absorption();
    absorb->get_back_optical_properties()->set_ideal_absorption();
    absorb->set_aperture(make_aperture<Circle>(2.0 * ABS_RADIUS));
    absorb->set_surface(make_surface<Sphere>(1.0 / ABS_RADIUS));
    // absorb->set_surface(make_surface<Flat>());

    aim_point = abs_origin;
    aim_point[2] += vector_norm(aim_point);
    absorb->set_reference_frame_geometry(abs_origin, aim_point, 0.0);
    absorb->set_name("Absorber");
    absorb->enable();
    auto ret = st2->add_element(absorb);
    EXPECT_TRUE(SolTrace::Data::Element::is_success(ret));

    sdata.add_stage(st1);
    sdata.add_stage(st2);

    auto sun = make_ray_source<Sun>();
    sun->set_position(sun_pos);
    sun->set_shape(DistributionType::PILLBOX, 0.0, 0.5);
    sdata.add_ray_source(sun);

    RunnerStatus sts = my_runner.initialize();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    // Setup runs but is not complete
    sts = my_runner.setup_simulation(&sdata);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    // Run simulation runs but returns RunnerStatus::ERROR
    auto t0 = std::chrono::high_resolution_clock::now();
    sts = my_runner.run_simulation();
    auto t1 = std::chrono::high_resolution_clock::now();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);

    std::chrono::duration<double, std::milli> dur = t1 - t0;
    EXPECT_TRUE(dur.count() < 20000.0);

    const TSystem *sys = my_runner.get_system();
    // sys->RayData.Print();
    const TRayData *ray_data = &(sys->RayData);
    size_t n = ray_data->Count();
    uint_fast64_t num_absorbed = count_absorbed_native(ray_data);

    std::cout << "Time: " << dur.count() << " ms" << std::endl;
    std::cout << "Number Absorbed: " << num_absorbed << std::endl;
    std::cout << "Number Interactions: " << n << std::endl;

    EXPECT_TRUE(n >= NRAYS);
    EXPECT_TRUE(num_absorbed > 0);

    // sys->AllRayData.Print();
}

// TEST(NativeRunner, LargePerformanceTest)
// {
//     // Pulling in path variable from CMake and creating path to .stinput sample file
//     std::string sample_path = std::string(PROJECT_DIR) +
//                               std::string("/Power-tower-surround_singlefacet.stinput");

//     // Path to .csv exported from Soltrace as ground truth
//     std::string ground_csv_path = PROJECT_DIR +
//                                   std::string("/powertower_example_raydata.csv");

//     // std::ifstream csv_file(ground_csv_path);
//     // std::vector<std::vector<std::string>> ground_raydata = split_csv(
//     //     ground_csv_path);

//     // const char *file = sample_path.data();

//     // // Soltrace case parameters
//     // int nrays = 50000;
//     // int maxrays = 5000000;
//     // int seed = 1; // Any positive integer will produce the same results each time, -1 will be a random seed
//     // int sunshape = 0;
//     // int errors = 0;
//     // int powertower = 0; // Toggles optimizations for power tower cases

//     // Create Simuluation Data
//     SimulationData sd;

//     // Constants
//     const uint_fast64_t NRAYS = 50000;
//     const double TOL = 1e-4;

//     // Read Input File
//     bool success = sd.import_from_file(sample_path);
//     EXPECT_TRUE(success);
//     EXPECT_TRUE(sd.get_number_of_elements() > 0);
//     EXPECT_TRUE(sd.get_number_of_ray_sources() > 0);

//     std::cout << "Num Elements: " << sd.get_number_of_elements() << std::endl;

//     // Parameters
//     SimulationParameters &params = sd.get_simulation_parameters();
//     params.include_optical_errors = false;
//     params.include_sun_shape_errors = false;
//     params.max_number_of_rays = NRAYS * 100;
//     params.number_of_rays = NRAYS;
//     params.seed = 1;

//     // Run Ray Trace
//     NativeRunner runner;
//     runner.disable_point_focus();
//     runner.disable_power_tower();
//     RunnerStatus sts;
//     sts = runner.initialize();
//     EXPECT_EQ(sts, RunnerStatus::SUCCESS);
//     sts = runner.setup_simulation(&sd);
//     EXPECT_EQ(sts, RunnerStatus::SUCCESS);

//     auto t0 = std::chrono::high_resolution_clock::now();
//     sts = runner.run_simulation();
//     auto t1 = std::chrono::high_resolution_clock::now();
//     EXPECT_EQ(sts, RunnerStatus::SUCCESS);

//     std::chrono::duration<double, std::milli> dur = t1 - t0;
//     EXPECT_TRUE(dur.count() < 60000.0);

//     const TSystem *sys = runner.get_system();
//     // sys->RayData.Print();
//     const TRayData *ray_data = &(sys->RayData);
//     size_t n = ray_data->Count();
//     uint_fast64_t num_absorbed = count_absorbed_native(ray_data);

//     std::cout << "Time: " << dur.count() << " ms" << std::endl;
//     std::cout << "Number Absorbed: " << num_absorbed << std::endl;
//     std::cout << "Number Interactions: " << n << std::endl;

//     EXPECT_TRUE(n >= NRAYS);
//     EXPECT_TRUE(num_absorbed > 0);

//     // const TSystem *sys = runner.get_system();
//     // const TRayData *ray_data = &(sys->AllRayData);
//     // size_t nrdata = ray_data->Count();

//     // Vector3d point, cosines;
//     // int element;
//     // int stage;
//     // unsigned int raynum;

//     // // Retrieving data from Runner
//     // for (size_t i = 0; i < nrdata; i++)
//     // {
//     //     EXPECT_TRUE(ray_data->Query(i, point.data, cosines.data,
//     //                                 &element, &stage, &raynum));

//     //     EXPECT_EQ(element, stoi(ground_raydata[6][i + 1]));
//     //     EXPECT_EQ(stage, stoi(ground_raydata[7][i + 1]));
//     //     EXPECT_EQ(raynum, stoul(ground_raydata[8][i + 1]));

//     //     if (element == stoi(ground_raydata[6][i + 1]) &&
//     //         stage == stoi(ground_raydata[7][i + 1]) &&
//     //         raynum == stoul(ground_raydata[8][i + 1]))
//     //     {
//     //         EXPECT_NEAR(point[0], stod(ground_raydata[0][i + 1]), TOL);
//     //         EXPECT_NEAR(point[1], stod(ground_raydata[1][i + 1]), TOL);
//     //         EXPECT_NEAR(point[2], stod(ground_raydata[2][i + 1]), TOL);

//     //         EXPECT_NEAR(cosines[0], stod(ground_raydata[3][i + 1]), TOL);
//     //         EXPECT_NEAR(cosines[1], stod(ground_raydata[4][i + 1]), TOL);
//     //         EXPECT_NEAR(cosines[2], stod(ground_raydata[5][i + 1]), TOL);
//     //     }
//     //     else
//     //     {
//     //         std::cout << "Line: " << i
//     //                   << "\nElement: " << element
//     //                   << "\nStage: " << stage
//     //                   << "\nRay Number: " << raynum
//     //                   << std::endl;
//     //         break;
//     //     }
//     // }
// }

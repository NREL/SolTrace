#include <gtest/gtest.h>

#include <native_runner.hpp>
#include <native_runner_types.hpp>
#include <simulation_data.hpp>
#include <simulation_result_export.hpp>
#include <stage_element.hpp>
#include <sun.hpp>
#include <utilities.hpp>

#include <cst_templates/heliostat.hpp>

#include "common.hpp"
#include "count_absorbed_native.h"

using Heliostat = SolTrace::Data::Heliostat;

using SolTrace::Runner::RunnerStatus;
using SolTrace::NativeRunner::NativeRunner;
using SolTrace::NativeRunner::TRayData;
using SolTrace::NativeRunner::TSystem;
using SolTrace::NativeRunner::TSun;

// Error Checking Tests for Heliostat
TEST(Heliostat, ErrorChecking_SetApertureSize)
{
    auto hs = SolTrace::Data::make_element<Heliostat>();

    // Test negative aperture size
    EXPECT_THROW(hs->set_aperture_size(-12.0, 12.0), std::invalid_argument);
    EXPECT_THROW(hs->set_aperture_size(12.0, -12.0), std::invalid_argument);
    EXPECT_THROW(hs->set_aperture_size(-12.0, -12.0), std::invalid_argument);

    // Test zero aperture size
    EXPECT_THROW(hs->set_aperture_size(0.0, 12.0), std::invalid_argument);
    EXPECT_THROW(hs->set_aperture_size(12.0, 0.0), std::invalid_argument);

    // Test valid aperture sizes
    EXPECT_NO_THROW(hs->set_aperture_size(12.0, 12.0));
    EXPECT_NO_THROW(hs->set_aperture_size(15.5, 10.3));
}

TEST(Heliostat, ErrorChecking_SetFocalLength)
{
    auto hs = SolTrace::Data::make_element<Heliostat>();

    // Test negative focal length
    EXPECT_THROW(hs->set_focal_length(-156.06), std::invalid_argument);
    EXPECT_THROW(hs->set_focal_length(-0.1), std::invalid_argument);

    // Test valid focal lengths (including zero for flat heliostat)
    EXPECT_NO_THROW(hs->set_focal_length(0.0));
    EXPECT_NO_THROW(hs->set_focal_length(156.06));
    EXPECT_NO_THROW(hs->set_focal_length(250.0));
}

TEST(Heliostat, ErrorChecking_SetFocalLengthXY)
{
    auto hs = SolTrace::Data::make_element<Heliostat>();

    // Test negative focal lengths
    EXPECT_THROW(hs->set_focal_length(-156.06, 156.06), std::invalid_argument);
    EXPECT_THROW(hs->set_focal_length(156.06, -156.06), std::invalid_argument);
    EXPECT_THROW(hs->set_focal_length(-156.06, -156.06), std::invalid_argument);

    // Test valid focal lengths (including zero for flat heliostat)
    EXPECT_NO_THROW(hs->set_focal_length(0.0, 0.0));
    EXPECT_NO_THROW(hs->set_focal_length(156.06, 156.06));
    EXPECT_NO_THROW(hs->set_focal_length(200.0, 250.0));
}

TEST(Heliostat, ErrorChecking_SetNumberPanels)
{
    auto hs = SolTrace::Data::make_element<Heliostat>();

    // Test invalid panel counts
    EXPECT_THROW(hs->set_number_panels(0, 4), std::invalid_argument);
    EXPECT_THROW(hs->set_number_panels(3, 0), std::invalid_argument);
    // Note: negative values would be caught by uint_fast64_t type

    // Test valid panel counts
    EXPECT_NO_THROW(hs->set_number_panels(1, 1));
    EXPECT_NO_THROW(hs->set_number_panels(3, 4));
    EXPECT_NO_THROW(hs->set_number_panels(10, 20));
}

TEST(Heliostat, ErrorChecking_SetGaps)
{
    auto hs = SolTrace::Data::make_element<Heliostat>();

    // Test negative gap values
    EXPECT_THROW(hs->set_gaps(-0.1, 0.1), std::invalid_argument);
    EXPECT_THROW(hs->set_gaps(0.1, -0.1), std::invalid_argument);

    // Test valid gap values (including zero)
    EXPECT_NO_THROW(hs->set_gaps(0.0, 0.0));
    EXPECT_NO_THROW(hs->set_gaps(0.1, 0.1));
}

TEST(Heliostat, ErrorChecking_SetCanting)
{
    auto hs = SolTrace::Data::make_element<Heliostat>();

    // Test valid canting types
    EXPECT_NO_THROW(hs->set_canting(Heliostat::NONE, 0.0, 0.0));
    EXPECT_NO_THROW(hs->set_canting(Heliostat::ON_AXIS, 100.0, 0.0));
    EXPECT_NO_THROW(hs->set_canting(Heliostat::OFF_AXIS, 45.0, 30.0));

    // Test invalid canting parameters
    EXPECT_THROW(hs->set_canting(Heliostat::ON_AXIS, -100.0, 0.0),
                 std::invalid_argument);
    EXPECT_THROW(hs->set_canting(Heliostat::OFF_AXIS, -1.0, 30.0),
                 std::invalid_argument);
    EXPECT_THROW(hs->set_canting(Heliostat::OFF_AXIS, 45.0, -1.0),
                 std::invalid_argument);
}

TEST(Heliostat, ErrorChecking_CreateGeometryWithoutParameters)
{
    auto hs = SolTrace::Data::make_element<Heliostat>();

    // Test create_geometry without setting required parameters
    EXPECT_THROW(hs->create_geometry(), std::invalid_argument);

    // Set aperture size and test again
    hs->set_aperture_size(12.0, 12.0);
    EXPECT_THROW(hs->create_geometry(), std::invalid_argument);

    // Set focal length and test again
    hs->set_focal_length(156.06);
    EXPECT_THROW(hs->create_geometry(), std::invalid_argument);

    // Set number of panels and test again
    hs->set_number_panels(3, 4);
    EXPECT_THROW(hs->create_geometry(), std::invalid_argument);

    hs->set_gaps(0.0, 0.0);
    EXPECT_THROW(hs->create_geometry(), std::invalid_argument);

    // Set canting and test again
    hs->set_canting(SolTrace::Data::Heliostat::NONE, 0.0, 0.0);
    EXPECT_THROW(hs->create_geometry(), std::invalid_argument);

    hs->set_target_position(Vector3d(0.0, 0.0, 10.0));
    EXPECT_NO_THROW(hs->create_geometry());
}

TEST(Heliostat, BuildParabolaNone)
{
    OpticalProperties mirror;
    mirror.set_ideal_reflection();

    auto hs = SolTrace::Data::make_element<Heliostat>();
    hs->set_mirror_optics(mirror);
    hs->set_origin(1.0, 1.0, 0.0);
    hs->set_aim_vector(0.0, 0.0, 100.0);
    hs->set_aperture_size(12.0, 12.0);
    hs->set_number_panels(3, 4);
    hs->set_gaps(0.1, 0.1);
    hs->set_focal_length(156.06);
    // hs->set_focal_point(0.0, 0.0, 156.06);
    hs->set_canting(Heliostat::NONE, 0.0, 0.0);
    hs->set_target_position(Vector3d(0.0, 0.0, 1.0));
    hs->create_geometry();

    // TODO: Check that everything ends up in the proper position
}

TEST(Heliostat, BuildFlatOnAxis)
{
    OpticalProperties mirror;
    mirror.set_ideal_reflection();

    auto hs = SolTrace::Data::make_element<Heliostat>();
    hs->set_mirror_optics(mirror);
    hs->set_origin(1.0, 1.0, 0.0);
    hs->set_aim_vector(0.0, 0.0, 100.0);
    hs->set_aperture_size(12.0, 12.0);
    hs->set_number_panels(3, 4);
    hs->set_gaps(0.1, 0.1);
    hs->set_focal_length(0.0);
    // hs->set_focal_point(Vector3d(0.0, 0.0, 10.0));
    hs->set_canting(Heliostat::NONE, 0.0, 0.0);
    hs->set_target_position(Vector3d(0.0, 0.0, 1.0));
    hs->create_geometry();

    // TODO: Check that everything ends up in the proper position
}

TEST(Heliostat, Trace)
{
    constexpr uint_fast64_t NRAYS = 10000;
    constexpr uint_fast64_t N_ABSORBED_THRESH = NRAYS / 10;
    const Vector3d zero(0.0, 0.0, 0.0);
    const Vector3d khat(0.0, 0.0, 1.0);

    SimulationData my_sim;
    // Set parameters
    SimulationParameters &params = my_sim.get_simulation_parameters();
    params.number_of_rays = NRAYS;
    params.max_number_of_rays = params.number_of_rays * 100;
    params.include_optical_errors = false;
    params.include_sun_shape_errors = false;
    params.seed = 12345;

    NativeRunner my_runner;
    my_runner.disable_power_tower();
    my_runner.enable_point_focus();

    OpticalProperties mirror;
    mirror.set_ideal_reflection();

    stage_ptr st1 = SolTrace::Data::make_stage(1);
    st1->set_reference_frame_geometry(zero, khat, 0.0);
    stage_ptr st2 = SolTrace::Data::make_stage(2);
    st2->set_reference_frame_geometry(zero, khat, 0.0);

    Vector3d sun_pos(0.0, 0.0, 1000.0);
    Vector3d hs_origin(1.0, 1.0, 0.0);
    Vector3d abs_origin(0.0, 0.0, 10.0);
    Vector3d v1;
    Vector3d v2;
    Vector3d aim;
    Vector3d aim_point;
    vector_add(1.0, sun_pos, -1.0, hs_origin, v1);
    vector_add(1.0, abs_origin, -1.0, hs_origin, v2);
    vector_add(0.5, v1, 0.5, v2, aim);
    vector_add(1.0, hs_origin, 1.0, aim, aim_point);

    auto hs = SolTrace::Data::make_element<Heliostat>();
    hs->set_mirror_optics(mirror);
    // hs->set_origin(hs_origin);
    // hs->set_aim_vector(0.0, 0.0, 2.0);
    // hs->set_zrot(0.0);
    // hs->compute_coordinate_rotations();
    // hs->convert_global_to_local(aim, abs_origin);
    hs->set_reference_frame_geometry(hs_origin, aim_point, 0.0);
    hs->set_aperture_size(12.0, 12.0);
    hs->set_number_panels(3, 4);
    hs->set_gaps(0.1, 0.1);
    hs->set_focal_length(0.0);
    hs->set_canting(Heliostat::NONE, 0.0, 0.0);
    hs->set_target_position(abs_origin);
    hs->create_geometry();
    hs->set_name("Heliostat");
    hs->enable();

    auto ret = st1->add_element(hs);
    EXPECT_TRUE(SolTrace::Data::Element::is_success(ret));

    auto absorb = SolTrace::Data::make_element<SingleElement>();
    absorb->get_front_optical_properties()->set_ideal_absorption();
    absorb->get_back_optical_properties()->set_ideal_absorption();
    absorb->set_aperture(SolTrace::Data::make_aperture<SolTrace::Data::Rectangle>(5.0, 5.0));
    absorb->set_surface(SolTrace::Data::make_surface<SolTrace::Data::Flat>());
    // absorb->set_origin(abs_origin);
    // absorb->set_aim_vector(0.0, 0.0, -1.0);
    // absorb->set_zrot(0.0);
    // absorb->compute_coordinate_rotations();
    // aim.scalar_mult(-1.0);
    vector_add(1.0, hs_origin, -1.0, abs_origin, aim);
    vector_add(1.0, abs_origin, 1.0, aim, aim_point);
    absorb->set_reference_frame_geometry(abs_origin, aim_point, 0.0);
    absorb->set_name("Absorber");
    absorb->enable();
    ret = st2->add_element(absorb);
    EXPECT_TRUE(SolTrace::Data::Element::is_success(ret));

    my_sim.add_stage(st1);
    my_sim.add_stage(st2);

    auto sun = SolTrace::Data::make_ray_source<Sun>();
    sun->set_position(sun_pos);
    sun->set_shape(SolTrace::Data::SunShape::GAUSSIAN, 1.0, 0.0, 0.0);
    my_sim.add_ray_source(sun);

    // // We can go over all the elements added
    // for (auto iter = my_sim.get_iterator();
    //      !my_sim.is_at_end(iter);
    //      ++iter)
    // {
    //     // iter is a iterator over the storing container which is a map
    //     // so that the iterator gives the key value pair
    //     element_id id = iter->first;
    //     // `element_ptr` is a std::shared_pointer to an Element
    //     element_ptr el = iter->second;
    //     if (el->is_stage())
    //     {
    //         continue;
    //     }
    //     std::cout << "------------\n"
    //               << "Element ID: " << id
    //               << "\nElement name: " << el->get_name()
    //               << "\nIs Stage: " << el->is_stage()
    //               << "\nIs Composite: " << el->is_composite()
    //               << "\nIs Single: " << el->is_single()
    //               // Below are all the same in this case
    //               << "\nOrigin (ref): " << el->get_origin_ref()
    //               << "\nOrigin (stage): " << el->get_origin_stage()
    //               << "\nOrigin (global): " << el->get_origin_global()
    //               << "\nAim (ref): " << el->get_aim_vector_ref()
    //               << "\nAim (stage): " << el->get_aim_vector_stage()
    //               << "\nAim (global): " << el->get_aim_vector_global()
    //               << "\n";
    // }

    // std::cout << "Number of elements in sim data: "
    //           << my_sim.get_number_of_elements()
    //           << std::endl;

    RunnerStatus sts = my_runner.initialize();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    // Setup runs but is not complete
    sts = my_runner.setup_simulation(&my_sim);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    // Run simulation runs but returns RunnerStatus::ERROR
    sts = my_runner.run_simulation();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);

    const TSystem *sys = my_runner.get_system();
    // sys->RayData.Print();
    const TRayData *ray_data = &(sys->RayData);
    size_t n = ray_data->Count();
    uint_fast64_t num_absorbed = count_absorbed_native(ray_data);

    std::cout << "Number Absorbed: " << num_absorbed << std::endl;
    std::cout << "Number Interactions: " << n << std::endl;

    // ray_data->Print();

    EXPECT_TRUE(n >= NRAYS);
    EXPECT_TRUE(num_absorbed > N_ABSORBED_THRESH);
}

TEST(Heliostat, TraceOffAxisCanting)
{
    constexpr uint_fast64_t NRAYS = 10000;
    constexpr uint_fast64_t N_ABSORBED_THRESH = NRAYS / 10;
    const Vector3d zero(0.0, 0.0, 0.0);
    const Vector3d khat(0.0, 0.0, 1.0);

    SimulationData my_sim;
    // Set parameters
    SimulationParameters& params = my_sim.get_simulation_parameters();
    params.number_of_rays = NRAYS;
    params.max_number_of_rays = params.number_of_rays * 100;
    params.include_optical_errors = false;
    params.include_sun_shape_errors = false;
    params.seed = 12345;

    NativeRunner my_runner;
    my_runner.disable_power_tower();
    my_runner.disable_point_focus();

    OpticalProperties mirror;
    mirror.set_ideal_reflection();

    stage_ptr st1 = SolTrace::Data::make_stage(1);
    st1->set_reference_frame_geometry(zero, khat, 0.0);
    stage_ptr st2 = SolTrace::Data::make_stage(2);
    st2->set_reference_frame_geometry(zero, khat, 0.0);

    Vector3d hs_origin(50.0, 50.0, 5.0);
    Vector3d abs_origin(0.0, 0.0, 5.0);
    double canting_azimuth = 135.0;
    double canting_zenith = 90.0;

    auto hs = SolTrace::Data::make_element<Heliostat>();
    hs->set_mirror_optics(mirror);
    hs->set_reference_frame_geometry(hs_origin, khat, 0.0);
    hs->set_aperture_size(12.0, 8.0);
    hs->set_number_panels(5, 5);
    hs->set_gaps(0.1, 0.1);
    hs->set_focal_length(0.0);
    hs->set_canting(Heliostat::OFF_AXIS, canting_azimuth, canting_zenith);
    hs->set_target_position(abs_origin);
    hs->create_geometry();
    hs->set_name("Heliostat");
    hs->enable();
    hs->update_geometry(canting_azimuth, 90.0 - canting_zenith); // Set to canting angles

    auto ret = st1->add_element(hs);
    EXPECT_TRUE(SolTrace::Data::Element::is_success(ret));

    auto absorb = SolTrace::Data::make_element<SingleElement>();
    absorb->get_front_optical_properties()->set_ideal_absorption();
    absorb->get_back_optical_properties()->set_ideal_absorption();
    absorb->set_aperture(SolTrace::Data::make_aperture<SolTrace::Data::Rectangle>(10.0, 10.0)); // TODO: Set a tight aperture (2.35, 1.55)
    absorb->set_surface(SolTrace::Data::make_surface<SolTrace::Data::Flat>());
    Vector3d v1 = { 0.0, 1.0, 0.0 };
    Vector3d aim_point;
    vector_add(1.0, abs_origin, 1.0, v1, aim_point);
    absorb->set_reference_frame_geometry(abs_origin, aim_point, 0.0);
    absorb->set_name("Absorber");
    absorb->enable();
    ret = st2->add_element(absorb);
    EXPECT_TRUE(SolTrace::Data::Element::is_success(ret));

    my_sim.add_stage(st1);
    my_sim.add_stage(st2);

    Vector3d sun_pos;
    sun_position_vector_degrees(sun_pos, canting_azimuth, 90.0 - canting_zenith);
    auto sun = SolTrace::Data::make_ray_source<Sun>();
    sun->set_position(sun_pos);
    sun->set_shape(SolTrace::Data::SunShape::PILLBOX, 0.0, 4.65, 0.0);
    my_sim.add_ray_source(sun);

    // // We can go over all the elements added
    // for (auto iter = my_sim.get_iterator();
    //      !my_sim.is_at_end(iter);
    //      ++iter)
    // {
    //     // iter is a iterator over the storing container which is a map
    //     // so that the iterator gives the key value pair
    //     element_id id = iter->first;
    //     // `element_ptr` is a std::shared_pointer to an Element
    //     element_ptr el = iter->second;
    //     if (el->is_stage())
    //     {
    //         continue;
    //     }
    //     std::cout << "------------\n"
    //               << "Element ID: " << id
    //               << "\nElement name: " << el->get_name()
    //               << "\nIs Stage: " << el->is_stage()
    //               << "\nIs Composite: " << el->is_composite()
    //               << "\nIs Single: " << el->is_single()
    //               // Below are all the same in this case
    //               << "\nOrigin (ref): " << el->get_origin_ref()
    //               << "\nOrigin (stage): " << el->get_origin_stage()
    //               << "\nOrigin (global): " << el->get_origin_global()
    //               << "\nAim (ref): " << el->get_aim_vector_ref()
    //               << "\nAim (stage): " << el->get_aim_vector_stage()
    //               << "\nAim (global): " << el->get_aim_vector_global()
    //               << "\n";
    // }

    // std::cout << "Number of elements in sim data: "
    //           << my_sim.get_number_of_elements()
    //           << std::endl;

    RunnerStatus sts = my_runner.initialize();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    // Setup runs but is not complete
    sts = my_runner.setup_simulation(&my_sim);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    sts = my_runner.run_simulation();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    SimulationResult result;
    sts = my_runner.report_simulation(&result, 0);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    //result.write_csv_file("native_runner_result_dump.csv");

    const TSystem* sys = my_runner.get_system();
    const TRayData* ray_data = &(sys->RayData);
    size_t n = ray_data->Count();
    uint_fast64_t num_absorbed = count_absorbed_native(ray_data);

    std::cout << "Number Absorbed: " << num_absorbed << std::endl;
    std::cout << "Number Interactions: " << n << std::endl;

    EXPECT_TRUE(n >= NRAYS);
    EXPECT_TRUE(num_absorbed > N_ABSORBED_THRESH);
}

TEST(Heliostat, ErrorChecking_UpdateGeometry)
{
    OpticalProperties mirror;
    mirror.set_ideal_reflection();

    Vector3d sun_pos(0.0, 0.0, 1000.0);
    Vector3d hs_origin(1.0, 1.0, 0.0);
    Vector3d abs_origin(0.0, 0.0, 10.0);
    Vector3d v1;
    Vector3d v2;
    Vector3d aim;
    Vector3d aim_point;
    vector_add(1.0, sun_pos, -1.0, hs_origin, v1);
    vector_add(1.0, abs_origin, -1.0, hs_origin, v2);
    vector_add(0.5, v1, 0.5, v2, aim);
    vector_add(1.0, hs_origin, 1.0, aim, aim_point);

    auto hs = SolTrace::Data::make_element<Heliostat>();
    hs->set_mirror_optics(mirror);
    // hs->set_origin(hs_origin);
    // hs->set_aim_vector(0.0, 0.0, 2.0);
    // hs->set_zrot(0.0);
    // hs->compute_coordinate_rotations();
    // hs->convert_global_to_local(aim, abs_origin);
    hs->set_reference_frame_geometry(hs_origin, aim_point, 0.0);
    hs->set_aperture_size(12.0, 12.0);
    hs->set_number_panels(3, 4);
    hs->set_gaps(0.1, 0.1);
    hs->set_focal_length(0.0);
    hs->set_canting(Heliostat::NONE, 0.0, 0.0);
    hs->set_target_position(abs_origin);

    EXPECT_THROW(hs->update_geometry(10.0, -10.0), std::invalid_argument);
    EXPECT_THROW(hs->update_geometry(10.0, 100.0), std::invalid_argument);
    EXPECT_THROW(hs->update_geometry(-200.0, 40.0), std::invalid_argument);
    EXPECT_THROW(hs->update_geometry(200.0, 40.0), std::invalid_argument);
    EXPECT_THROW(hs->update_geometry(0.0, 40.0), std::runtime_error);

    hs->create_geometry();
    hs->set_name("Heliostat");
    hs->enable();

    EXPECT_NO_THROW(hs->update_geometry(0.0, 40.0));
}

TEST(Heliostat, UpdateGeometry)
{
    constexpr uint_fast64_t NRAYS = 10000;
    constexpr uint_fast64_t N_ABSORBED_THRESH = NRAYS / 10;

    const double sun_az = 180.0;
    const double sun_el = 45.0;
    const double TOL = 1e-12;

    SimulationData my_sim;
    // Set parameters
    SimulationParameters &params = my_sim.get_simulation_parameters();
    params.number_of_rays = NRAYS;
    params.max_number_of_rays = params.number_of_rays * 100;
    params.include_optical_errors = false;
    params.include_sun_shape_errors = false;
    params.seed = 12345;

    NativeRunner my_runner;
    my_runner.disable_power_tower();
    my_runner.enable_point_focus();

    OpticalProperties mirror;
    mirror.set_ideal_reflection();

    Vector3d sun_pos;
    sun_position_vector_degrees(sun_pos, sun_az, sun_el);
    // Vector3d hs_origin(1.0, 1.0, 0.0);
    Vector3d abs_origin(0.0, 0.0, 2.0);

    auto hs = SolTrace::Data::make_element<Heliostat>();
    hs->set_mirror_optics(mirror);
    hs->set_origin(1.0, 1.0, 0.0);
    hs->set_aperture_size(12.0, 12.0);
    hs->set_number_panels(3, 4);
    hs->set_gaps(0.1, 0.1);
    hs->set_focal_length(0.0);
    hs->set_canting(Heliostat::NONE, 0.0, 0.0);
    hs->set_target_position(abs_origin);
    hs->create_geometry();
    hs->set_name("Heliostat");
    hs->enable();

    auto ret = my_sim.add_element(hs);
    EXPECT_TRUE(SolTrace::Data::Element::is_success(ret));

    hs->update_geometry(sun_az, sun_el);
    Vector3d result;
    vector_add(-1.0, hs->get_origin_global(),
               1.0, hs->get_aim_vector_global(),
               result);
    result.make_unit();
    Vector3d temp;
    vector_add(1.0, abs_origin, -1.0, hs->get_origin_global(), temp);
    temp.make_unit();
    double phi1 = acos(dot_product(result, sun_pos)) * SolTrace::Data::R2D;
    double phi2 = acos(dot_product(result, temp)) * SolTrace::Data::R2D;
    double phi3 = acos(dot_product(sun_pos, temp)) * SolTrace::Data::R2D;

    EXPECT_NEAR(phi1, phi2, TOL);
    EXPECT_NEAR(phi1 + phi2, phi3, TOL);

    // TODO: Test for correct z-rotation...

    auto absorb = SolTrace::Data::make_element<SingleElement>();
    absorb->get_front_optical_properties()->set_ideal_absorption();
    absorb->get_back_optical_properties()->set_ideal_absorption();
    absorb->set_aperture(SolTrace::Data::make_aperture<SolTrace::Data::Rectangle>(5.0, 5.0));
    absorb->set_surface(SolTrace::Data::make_surface<SolTrace::Data::Flat>());
    // absorb->set_origin(abs_origin);
    // absorb->set_aim_vector(0.0, 0.0, -1.0);
    // absorb->set_zrot(0.0);
    // absorb->compute_coordinate_rotations();
    // aim.scalar_mult(-1.0);
    // vector_add(1.0, hs_origin, -1.0, abs_origin, aim);
    // vector_add(1.0, abs_origin, 1.0, aim, aim_point);
    Vector3d aim_point(0.0, 0.0, 1.0);
    vector_add(1.0, abs_origin, 1.0, aim_point);
    absorb->set_reference_frame_geometry(abs_origin, aim_point, 0.0);
    absorb->set_name("Absorber");
    absorb->enable();
    ret = my_sim.add_element(absorb);
    EXPECT_TRUE(SolTrace::Data::Element::is_success(ret));

    auto sun = SolTrace::Data::make_ray_source<Sun>();
    sun->set_position(sun_pos);
    sun->set_shape(SolTrace::Data::SunShape::GAUSSIAN, 1.0, 0.0, 0.0);
    my_sim.add_ray_source(sun);

    // std::cout << "Sun Position: " << sun_pos
    //           << "\nHeliostat Aim Point: " << hs->get_aim_vector_global()
    //           << "\nHeliostat ZRot: " << hs->get_zrot()
    //           << "\nHeliostat Elevation Axis: " << hs->get_elevation_axis()
    //           << std::endl;

    // // We can go over all the elements added
    // for (auto iter = my_sim.get_iterator();
    //      !my_sim.is_at_end(iter);
    //      ++iter)
    // {
    //     // iter is a iterator over the storing container which is a map
    //     // so that the iterator gives the key value pair
    //     element_id id = iter->first;
    //     // `element_ptr` is a std::shared_pointer to an Element
    //     element_ptr el = iter->second;
    //     if (el->is_stage())
    //     {
    //         continue;
    //     }
    //     std::cout << "------------\n"
    //               << "Element ID: " << id
    //               << "\nElement name: " << el->get_name()
    //               << "\nIs Stage: " << el->is_stage()
    //               << "\nIs Composite: " << el->is_composite()
    //               << "\nIs Single: " << el->is_single()
    //               // Below are all the same in this case
    //               << "\nOrigin (ref): " << el->get_origin_ref()
    //               << "\nOrigin (stage): " << el->get_origin_stage()
    //               << "\nOrigin (global): " << el->get_origin_global()
    //               << "\nAim (ref): " << el->get_aim_vector_ref()
    //               << "\nAim (stage): " << el->get_aim_vector_stage()
    //               << "\nAim (global): " << el->get_aim_vector_global()
    //               << "\n";
    // }

    // std::cout << "Number of elements in sim data: "
    //           << my_sim.get_number_of_elements()
    //           << std::endl;

    RunnerStatus sts = my_runner.initialize();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    // Setup runs but is not complete
    sts = my_runner.setup_simulation(&my_sim);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    // Run simulation runs but returns RunnerStatus::ERROR
    sts = my_runner.run_simulation();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);

    const TSystem *sys = my_runner.get_system();
    // sys->RayData.Print();
    const TRayData *ray_data = &(sys->RayData);
    uint_fast64_t num_absorbed = count_absorbed_native(ray_data);
    size_t n = ray_data->Count();

    std::cout << "Number Absorbed: " << num_absorbed << std::endl;
    std::cout << "Number Interactions: " << n << std::endl;

    // ray_data->Print();

    EXPECT_TRUE(n >= NRAYS);
    EXPECT_TRUE(num_absorbed > N_ABSORBED_THRESH);
}

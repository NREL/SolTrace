#include <gtest/gtest.h>

#include <native_runner.hpp>
#include <native_runner_types.hpp>
#include <simulation_data.hpp>
#include <simulation_result_export.hpp>
#include <sun.hpp>
#include <utilities.hpp>

#include <cst_templates/arclength.hpp>
#include <cst_templates/linear_fresnel.hpp>

#include "common.hpp"
#include "count_absorbed_native.h"

using SolTrace::Data::LinearFresnel;

using SolTrace::Runner::RunnerStatus;
using SolTrace::NativeRunner::NativeRunner;
using SolTrace::NativeRunner::TRayData;
using SolTrace::NativeRunner::TSystem;

// Error Checking Tests for LinearFresnel
TEST(LinearFresnel, ErrorChecking_SetApertureSize)
{
    auto lf = SolTrace::Data::make_element<LinearFresnel>();

    // Test negative aperture size
    EXPECT_THROW(lf->set_aperture_size(-1.0, 5.0), std::invalid_argument);
    EXPECT_THROW(lf->set_aperture_size(5.0, -1.0), std::invalid_argument);
    EXPECT_THROW(lf->set_aperture_size(-1.0, -1.0), std::invalid_argument);

    // Test zero aperture size
    EXPECT_THROW(lf->set_aperture_size(0.0, 5.0), std::invalid_argument);
    EXPECT_THROW(lf->set_aperture_size(5.0, 0.0), std::invalid_argument);

    // Test valid aperture sizes
    EXPECT_NO_THROW(lf->set_aperture_size(1.0, 2.0));
    EXPECT_NO_THROW(lf->set_aperture_size(10.5, 15.3));
}

TEST(LinearFresnel, ErrorChecking_SetAngles)
{
    auto lf = SolTrace::Data::make_element<LinearFresnel>();

    // Test invalid azimuth angles
    EXPECT_THROW(lf->set_angles(-181.0, 45.0), std::invalid_argument);
    EXPECT_THROW(lf->set_angles(181.0, 45.0), std::invalid_argument);

    // Test invalid tilt angles
    EXPECT_THROW(lf->set_angles(0.0, -91.0), std::invalid_argument);
    EXPECT_THROW(lf->set_angles(0.0, 91.0), std::invalid_argument);

    // Test valid angles
    EXPECT_NO_THROW(lf->set_angles(-180.0, 0.0));
    EXPECT_NO_THROW(lf->set_angles(180.0, 90.0));
    EXPECT_NO_THROW(lf->set_angles(0.0, 45.0));
}

TEST(LinearFresnel, ErrorChecking_SetGaps)
{
    auto lf = SolTrace::Data::make_element<LinearFresnel>();

    // Test negative gap values
    EXPECT_THROW(lf->set_gaps(-0.1, 0.1, 0.1), std::invalid_argument);
    EXPECT_THROW(lf->set_gaps(0.1, -0.1, 0.1), std::invalid_argument);
    EXPECT_THROW(lf->set_gaps(0.1, 0.1, -0.1), std::invalid_argument);

    // Test valid gap values (including zero)
    EXPECT_NO_THROW(lf->set_gaps(0.0, 0.0, 0.0));
    EXPECT_NO_THROW(lf->set_gaps(0.1, 0.2, 0.3));
}

TEST(LinearFresnel, ErrorChecking_SetNumberPanels)
{
    auto lf = SolTrace::Data::make_element<LinearFresnel>();

    // Test invalid panel counts
    EXPECT_THROW(lf->set_number_panels(0, 5), std::invalid_argument);
    EXPECT_THROW(lf->set_number_panels(5, 0), std::invalid_argument);
    EXPECT_THROW(lf->set_number_panels(-1, 5), std::invalid_argument);
    EXPECT_THROW(lf->set_number_panels(5, -1), std::invalid_argument);

    // Test valid panel counts
    EXPECT_NO_THROW(lf->set_number_panels(1, 1));
    EXPECT_NO_THROW(lf->set_number_panels(10, 20));
}

TEST(LinearFresnel, ErrorChecking_SetReceiverHeight)
{
    auto lf = SolTrace::Data::make_element<LinearFresnel>();

    // Test invalid receiver heights
    EXPECT_THROW(lf->set_receiver_height(0.0), std::invalid_argument);
    EXPECT_THROW(lf->set_receiver_height(-1.0), std::invalid_argument);

    // Test valid receiver heights
    EXPECT_NO_THROW(lf->set_receiver_height(0.1));
    EXPECT_NO_THROW(lf->set_receiver_height(10.0));
}

TEST(LinearFresnel, ErrorChecking_SetReceiverDimensions)
{
    auto lf = SolTrace::Data::make_element<LinearFresnel>();

    // Test invalid absorber diameter
    EXPECT_THROW(lf->set_receiver_dimensions(0.0, 1.0, 0.1), std::invalid_argument);
    EXPECT_THROW(lf->set_receiver_dimensions(-0.5, 1.0, 0.1), std::invalid_argument);

    // Test invalid envelope diameter
    EXPECT_THROW(lf->set_receiver_dimensions(0.5, 0.0, 0.1), std::invalid_argument);
    EXPECT_THROW(lf->set_receiver_dimensions(0.5, -1.0, 0.1), std::invalid_argument);

    // Test invalid envelope thickness
    EXPECT_THROW(lf->set_receiver_dimensions(0.5, 1.0, -0.1), std::invalid_argument);

    // Test absorber diameter >= envelope diameter
    EXPECT_THROW(lf->set_receiver_dimensions(1.0, 1.0, 0.1), std::invalid_argument);
    EXPECT_THROW(lf->set_receiver_dimensions(1.5, 1.0, 0.1), std::invalid_argument);

    // Test valid receiver dimensions
    EXPECT_NO_THROW(lf->set_receiver_dimensions(0.5, 1.0, 0.0));
    EXPECT_NO_THROW(lf->set_receiver_dimensions(0.07, 0.115, 0.003));
}

TEST(LinearFresnel, ErrorChecking_CreateGeometryWithoutParameters)
{
    auto lf = SolTrace::Data::make_element<LinearFresnel>();

    // Test create_geometry without setting required parameters
    EXPECT_THROW(lf->create_geometry(), std::invalid_argument);

    // Set aperture size and test again
    lf->set_aperture_size(5.0, 10.0);
    EXPECT_THROW(lf->create_geometry(), std::invalid_argument);

    // Set number of panels and test again
    lf->set_number_panels(5, 2);
    EXPECT_THROW(lf->create_geometry(), std::invalid_argument);

    // Set receiver height and test again
    lf->set_receiver_height(3.0);
    EXPECT_THROW(lf->create_geometry(), std::invalid_argument);

    // Set receiver dimensions and it should work
    lf->set_receiver_dimensions(0.07, 0.115, 0.003);
    EXPECT_NO_THROW(lf->create_geometry());
}

TEST(LinearFresnel, Build)
{
    OpticalProperties mirror;
    mirror.set_ideal_reflection();

    OpticalProperties absorber;
    absorber.set_ideal_absorption();

    OpticalProperties envelop_out;
    envelop_out.set_ideal_transmission();
    envelop_out.refraction_index_front = 1.46;
    envelop_out.refraction_index_back = 1.0;

    OpticalProperties envelop_in;
    envelop_in.set_ideal_transmission();
    envelop_in.refraction_index_front = 1.0;
    envelop_in.refraction_index_back = 1.46;

    auto lf = SolTrace::Data::make_element<LinearFresnel>();
    lf->set_optics(mirror, absorber, envelop_out, envelop_in);
    lf->set_origin(10.0, -10.0, 0.0);
    lf->set_aperture_size(6.0, 12.0);
    lf->set_number_panels(10, 4);
    lf->set_gaps(0.15, 0.02, 0.15);
    lf->set_focused_panels(true);
    lf->set_receiver_height(2.0);
    lf->set_receiver_dimensions(0.07, 0.115, 0.003);
    lf->set_angles(0.0, 0.0);
    lf->create_geometry();

    lf = SolTrace::Data::make_element<LinearFresnel>();
    lf->set_optics(mirror, absorber, envelop_out, envelop_in);
    lf->set_origin(10.0, -10.0, 0.0);
    lf->set_aperture_size(6.0, 12.0);
    lf->set_number_panels(1, 1);
    lf->set_gaps(0.0, 0.01, 0.0);
    lf->set_focused_panels(false);
    lf->set_receiver_height(2.0);
    lf->set_receiver_dimensions(0.07, 0.115, 0.003);
    lf->set_angles(0.0, 0.0);
    lf->create_geometry();

    // TODO: Check that everything ends up in the proper position
}

TEST(LinearFresnel, Tracing)
{
    constexpr uint_fast64_t NRAYS = 10000;
    constexpr uint_fast64_t N_ABSORBED_THRESH = NRAYS / 10;

    SimulationData my_sim;
    // Set parameters
    SimulationParameters &params = my_sim.get_simulation_parameters();
    params.number_of_rays = NRAYS;
    params.max_number_of_rays = params.number_of_rays * 1000;
    params.include_optical_errors = false;
    params.include_sun_shape_errors = false;
    params.seed = 123;

    NativeRunner my_runner;
    my_runner.disable_power_tower();
    my_runner.disable_point_focus();

    OpticalProperties mirror;
    mirror.set_ideal_reflection();
    mirror.slope_error = 1.5;
    mirror.specularity_error = 0.5;

    OpticalProperties absorber;
    absorber.set_ideal_absorption();
    absorber.slope_error = 1e-5;
    absorber.specularity_error = 1e-5;

    OpticalProperties envelop_out;
    envelop_out.set_ideal_transmission();
    // envelop_out.refraction_index_front = 1.46;
    envelop_out.refraction_index_front = 1.0;
    envelop_out.refraction_index_back = 1.0;
    envelop_out.slope_error = 1e-4;
    envelop_out.specularity_error = 1e-4;

    OpticalProperties envelop_in;
    envelop_in.set_ideal_transmission();
    envelop_in.refraction_index_front = 1.0;
    // envelop_in.refraction_index_back = 1.46;
    envelop_in.refraction_index_back = 1.0;
    envelop_in.slope_error = 1e-4;
    envelop_in.specularity_error = 1e-4;

    auto lf = SolTrace::Data::make_element<LinearFresnel>();
    lf->set_optics(mirror, absorber, envelop_out, envelop_in);
    lf->set_origin(0.0, 0.0, 0.0);
    lf->set_aperture_size(6.0, 12.0);
    lf->set_number_panels(2, 2);
    lf->set_gaps(0.05, 0.02, 0.15);
    lf->set_focused_panels(true);
    lf->set_receiver_height(2.0);
    lf->set_receiver_dimensions(0.07, 0.115, 0.003);
    lf->set_angles(0.0, 0.0);
    lf->create_geometry();
    lf->set_name("LinearFresnel");

    auto sun = SolTrace::Data::make_ray_source<Sun>();
    sun->set_position(0.0, 0.0, 1000.0);
    // double NaN = std::numeric_limits<double>::quiet_NaN();
    sun->set_shape(SolTrace::Data::SunShape::GAUSSIAN, 1.0, 0.0, 0.0);
    my_sim.add_ray_source(sun);

    // Assumes that reference and global coordinates are the same
    // Vector3d pt_aim_point;
    // vector_add(1.0, sun->get_position(), -1.0, lf->get_origin_ref(), pt_aim_point);
    // lf->set_aim_vector(pt_aim_point);
    lf->set_aim_vector(sun->get_position());
    lf->set_zrot(0.0);
    lf->compute_coordinate_rotations();
    lf->enable();
    my_sim.add_element(lf);

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
    // for (size_t i = 0; i < n; i++)
    // {
    //     double pos[3], cos[3];
    //     int elm, stage;
    //     unsigned int ray;
    //     RayEvent rev;
    //     if (ray_data->Query(i, pos, cos, &elm, &stage, &ray, &rev))
    //     {
    //         if (elm < 0)
    //             ++num_absorbed;
    //     }
    // }

    std::cout << "Number Absorbed: " << num_absorbed << std::endl;
    std::cout << "Number Interactions: " << n << std::endl;

    // ray_data->Print();

    EXPECT_TRUE(n >= NRAYS);
    EXPECT_TRUE(num_absorbed > N_ABSORBED_THRESH);

    // SimulationResult my_res;
    // my_runner.report_simulation(&my_res, 0);

    // std::cout << my_res << std::endl;

    // my_res.write_csv_file("linear_fresnel_test.csv", 12);
}

TEST(LinearFresnel, UpdateGeometry)
{
    constexpr uint_fast64_t NRAYS = 10000;
    constexpr uint_fast64_t N_ABSORBED_THRESH = NRAYS / 10;

    const double sun_az = 180.0;
    const double sun_el = 45.0;
    const double TOL = 1e-12;
    // const double sun_az = 90.0;
    // const double sun_el = 00.0;

    SimulationData my_sim;
    // Set parameters
    SimulationParameters &params = my_sim.get_simulation_parameters();
    params.number_of_rays = NRAYS;
    params.max_number_of_rays = params.number_of_rays * 1000;
    params.include_optical_errors = true;
    params.include_sun_shape_errors = true;
    params.seed = 123;

    NativeRunner my_runner;
    my_runner.disable_power_tower();
    my_runner.disable_point_focus();

    OpticalProperties mirror;
    mirror.set_ideal_reflection();
    mirror.slope_error = 1.5;
    mirror.specularity_error = 0.5;

    OpticalProperties absorber;
    absorber.set_ideal_absorption();
    absorber.slope_error = 1e-5;
    absorber.specularity_error = 1e-5;

    OpticalProperties envelop_out;
    envelop_out.set_ideal_transmission();
    envelop_out.refraction_index_front = 1.46;
    envelop_out.refraction_index_back = 1.0;
    envelop_out.slope_error = 1e-4;
    envelop_out.specularity_error = 1e-4;

    OpticalProperties envelop_in;
    envelop_in.set_ideal_transmission();
    envelop_in.refraction_index_front = 1.0;
    envelop_in.refraction_index_back = 1.46;
    envelop_in.slope_error = 1e-4;
    envelop_in.specularity_error = 1e-4;

    auto sun = SolTrace::Data::make_ray_source<Sun>();
    Vector3d sun_pos;
    sun_position_vector_degrees(sun_pos, sun_az, sun_el);
    sun_pos.scalar_mult(1000.0);
    sun->set_position(sun_pos);
    sun->set_shape(SolTrace::Data::SunShape::GAUSSIAN, 1.0, 0.0, 0.0);
    my_sim.add_ray_source(sun);

    auto lf = SolTrace::Data::make_element<LinearFresnel>();
    lf->set_optics(mirror, absorber, envelop_out, envelop_in);
    lf->set_origin(10.0, 0.0, 0.0);
    // lf->set_origin(0.0, 0.0, 0.0);
    lf->set_aperture_size(6.0, 12.0);
    lf->set_number_panels(2, 2);
    lf->set_gaps(0.05, 0.02, 0.15);
    lf->set_focused_panels(true);
    lf->set_receiver_height(2.0);
    lf->set_receiver_dimensions(0.07, 0.115, 0.003);
    lf->set_angles(0.0, 40.0);
    // lf->set_angles(0.0, 0.0);
    lf->set_tracking_limits(-90.0, 90.0);
    lf->create_geometry();
    lf->set_name("LinearFresnel");
    lf->update_geometry(sun_az, sun_el);
    lf->enable();
    my_sim.add_element(lf);

    // EXPECT_NEAR(dot_product(lf->get_rotation_vector(),
    //                         lf->get_neutral_normal()),
    //             0.0, 1e-12);

    EXPECT_NEAR(dot_product(lf->get_tracking_origin(),
                            lf->get_rotation_vector()),
                0.0, TOL);
    EXPECT_NEAR(dot_product(lf->get_tracking_origin(),
                            lf->get_neutral_normal()),
                0.0, TOL);
    EXPECT_NEAR(dot_product(lf->get_rotation_vector(),
                            lf->get_neutral_normal()),
                0.0, TOL);

    // std::cout << "Rotation Axis: " << lf->get_rotation_vector()
    //           << "\nNeutral Normal: " << lf->get_neutral_normal()
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

    // SimulationResult my_res;
    // my_runner.report_simulation(&my_res, 0);
    // std::cout << my_res << std::endl;
}

TEST(LinearFresnel, UpdateGeometry_TrackingLimits)
{
    using SolTrace::Data::D2R;

    const double sun_az = 90.0;
    const double sun_el = 0.0;
    // const double sun_az = 90.0;
    // const double sun_el = 00.0;
    const double LOWER = -30.0;
    const double UPPER = 25.0;
    const double TOL = 1e-12;

    auto lf = SolTrace::Data::make_element<LinearFresnel>();
    lf->set_origin(0.0, 0.0, 0.0);
    lf->set_aperture_size(4.0, 8.0);
    lf->set_number_panels(2, 1);
    lf->set_gaps(0.0, 0.0, 0.0);
    lf->set_focused_panels(false);
    lf->set_receiver_height(2.0);
    lf->set_receiver_dimensions(0.07, 0.115, 0.003);
    lf->set_angles(0.0, 0.0);
    lf->set_tracking_limits(LOWER, UPPER);
    lf->create_geometry();
    lf->set_name("LinearFresnel");
    lf->enable();

    Vector3d normal;
    double theta;

    lf->update_geometry(-sun_az, sun_el);

    for (auto citer : lf->get_mirrors())
    {
        vector_add(1.0, citer->get_aim_vector_global(),
                   -1.0, citer->get_origin_global(),
                   normal);
        normal.make_unit();
        // Dot product with [0, 0, 1]
        theta = acos(normal[2]);
        if (normal[0] < 0.0)
            theta = -theta;

        EXPECT_NEAR(theta, LOWER * D2R, TOL);
        EXPECT_NEAR(normal[0], sin(LOWER * D2R), TOL);
        EXPECT_NEAR(normal[1], 0.0, TOL);
        EXPECT_NEAR(normal[2], cos(LOWER * D2R), TOL);
    }

    lf->update_geometry(sun_az, sun_el);
    for (auto citer : lf->get_mirrors())
    {
        vector_add(1.0, citer->get_aim_vector_global(),
                   -1.0, citer->get_origin_global(),
                   normal);
        normal.make_unit();
        // Dot product with [1, 0, 0]
        theta = acos(normal[2]);
        if (normal[0] < 0.0)
            theta = -theta;

        EXPECT_NEAR(theta, UPPER * D2R, TOL);
        EXPECT_NEAR(normal[0], sin(UPPER * D2R), TOL);
        EXPECT_NEAR(normal[1], 0.0, TOL);
        EXPECT_NEAR(normal[2], cos(UPPER * D2R), TOL);
    }
}

// Additional error condition tests for LinearFresnel
TEST(LinearFresnel, ErrorChecking_SetTrackingLimits)
{
    auto lf = SolTrace::Data::make_element<LinearFresnel>();
    // Lower limit > upper limit
    EXPECT_THROW(lf->set_tracking_limits(90.0, 0.0), std::invalid_argument);
    // Lower limit == upper limit (should be allowed, so no throw)
    EXPECT_NO_THROW(lf->set_tracking_limits(45.0, 45.0));
}

TEST(LinearFresnel, ErrorChecking_UpdateGeometryInvalidArgs)
{
    auto lf = SolTrace::Data::make_element<LinearFresnel>();
    lf->set_aperture_size(5.0, 10.0);
    lf->set_number_panels(2, 2);
    lf->set_receiver_height(3.0);
    lf->set_receiver_dimensions(0.07, 0.115, 0.003);
    lf->create_geometry();
    // Elevation out of range
    EXPECT_THROW(lf->update_geometry(0.0, -1.0), std::invalid_argument);
    EXPECT_THROW(lf->update_geometry(0.0, 91.0), std::invalid_argument);
    // Azimuth out of range
    EXPECT_THROW(lf->update_geometry(-181.0, 45.0), std::invalid_argument);
    EXPECT_THROW(lf->update_geometry(181.0, 45.0), std::invalid_argument);
}

TEST(LinearFresnel, ErrorChecking_UpdateGeometryBeforeCreate)
{
    auto lf = SolTrace::Data::make_element<LinearFresnel>();
    lf->set_aperture_size(5.0, 10.0);
    lf->set_number_panels(2, 2);
    lf->set_receiver_height(3.0);
    lf->set_receiver_dimensions(0.07, 0.115, 0.003);
    // Do not call create_geometry
    EXPECT_THROW(lf->update_geometry(0.0, 45.0), std::invalid_argument);
}

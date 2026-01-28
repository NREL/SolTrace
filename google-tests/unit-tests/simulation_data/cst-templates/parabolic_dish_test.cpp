#include <gtest/gtest.h>

#include <native_runner.hpp>
#include <native_runner_types.hpp>
#include <simulation_data.hpp>
#include <sun.hpp>
#include <utilities.hpp>

#include <cst_templates/arclength.hpp>
#include <cst_templates/parabolic_dish.hpp>

#include "common.hpp"
#include "count_absorbed_native.h"

using ParabolicDish = SolTrace::Data::ParabolicDish;

using SolTrace::Runner::RunnerStatus;
using SolTrace::NativeRunner::NativeRunner;
using SolTrace::NativeRunner::TRayData;
using SolTrace::NativeRunner::TSystem;

TEST(ParabolicDish, ArcLength)
{
    const double TOL = 1e-6;
    const double ARC_LENGTH = 4.105679289785514;

    double x0 = -1.0;
    double x1 = 2.0;
    double cx = 1.0;
    // double xstart = -hw;
    // double arc_length = hw * sqrt(hw * this->cx * hw * this->cx + 1) +
    //                     asinh(hw * this->cx) / this->cx;
    double val = SolTrace::Data::parabolic_arc_length(cx, x0, x1, 0.5 * TOL);
    EXPECT_NEAR(val, ARC_LENGTH, TOL);

    double xtest = SolTrace::Data::parabolic_determine_x_coordinate(
        cx, x0, ARC_LENGTH, 0.5 * TOL);
    EXPECT_NEAR(xtest, x1, TOL);
}

TEST(ParabolicDish, Build)
{
    OpticalProperties mirror;
    mirror.set_ideal_reflection();

    OpticalProperties absorber;
    absorber.set_ideal_absorption();

    auto dish = SolTrace::Data::make_element<ParabolicDish>();
    dish->set_optics(mirror, absorber);
    dish->set_origin(20.0, -20.0, 30.0);
    dish->set_aperture_size(10.0);
    dish->set_number_of_panels(2, 2);
    dish->set_gaps(0.02, 0.01, 0.5);
    dish->set_focal_length(7.5);
    dish->set_receiver_dimensions(0.25, 7.25);
    dish->create_geometry();

    dish = SolTrace::Data::make_element<ParabolicDish>();
    dish->set_optics(mirror, absorber);
    dish->set_origin(20.0, -20.0, 30.0);
    dish->set_aperture_size(10.0);
    dish->set_number_of_panels(1, 1);
    dish->set_gaps(0.02, 0.01, 0.5);
    dish->set_focal_length(7.5);
    dish->set_receiver_dimensions(0.25, 7.25);
    dish->create_geometry();

    // TODO: Check that everything ends up in the proper position
}

TEST(ParabolicDish, Tracing)
{
    const uint_fast64_t NRAYS = 10000;

    SimulationData my_sim;
    // Set parameters
    SimulationParameters &params = my_sim.get_simulation_parameters();
    params.number_of_rays = NRAYS;
    params.max_number_of_rays = params.number_of_rays * 100;
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

    auto dish = SolTrace::Data::make_element<ParabolicDish>();
    dish->set_optics(mirror, absorber);
    dish->set_origin(2.0, -2.0, 0.0);
    dish->set_aperture_size(10.0);
    dish->set_number_of_panels(2, 2);
    dish->set_gaps(0.02, 0.01, 0.5);
    dish->set_focal_length(7.5);
    dish->set_receiver_dimensions(0.5, 7.25);
    dish->set_name("ParabolicDish");
    dish->create_geometry();

    auto sun = SolTrace::Data::make_ray_source<Sun>();
    sun->set_position(0.0, 0.0, 1000.0);
    // double NaN = std::numeric_limits<double>::quiet_NaN();
    sun->set_shape(SolTrace::Data::SunShape::PILLBOX, 0.0, 1.0, 0.0);
    my_sim.add_ray_source(sun);

    // Assumes that reference and global coordinates are the same
    dish->set_aim_vector(sun->get_position());
    dish->set_zrot(0.0);
    dish->compute_coordinate_rotations();
    my_sim.add_element(dish);

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
    sts = my_runner.setup_simulation(&my_sim);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
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
    EXPECT_TRUE(num_absorbed > 0);
}

TEST(ParabolicDish, UpdateGeometry)
{
    constexpr uint_fast64_t NRAYS = 10000;
    constexpr uint_fast64_t N_ABSORBED_THRESH = NRAYS / 10;
    const double sun_az = 180.0;
    const double sun_el = 45.0;

    SimulationData my_sim;
    // Set parameters
    SimulationParameters &params = my_sim.get_simulation_parameters();
    params.number_of_rays = NRAYS;
    params.max_number_of_rays = params.number_of_rays * 100;
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

    auto dish = SolTrace::Data::make_element<ParabolicDish>();
    dish->set_optics(mirror, absorber);
    dish->set_origin(10.0, 2.0, 0.0);
    dish->set_aperture_size(10.0);
    dish->set_number_of_panels(2, 2);
    // dish->set_number_of_panels(1, 1);
    dish->set_gaps(0.02, 0.01, 0.5);
    dish->set_focal_length(7.5);
    dish->set_receiver_dimensions(0.5, 7.25);
    dish->set_name("ParabolicDish");
    dish->create_geometry();

    auto sun = SolTrace::Data::make_ray_source<Sun>();
    // sun->set_position(0.0, 0.0, 1000.0);
    // double NaN = std::numeric_limits<double>::quiet_NaN();
    sun->set_shape(SolTrace::Data::SunShape::PILLBOX, 0.0, 1.0, 0.0);
    sun_position_vector_degrees(sun->get_position(), sun_az, sun_el);
    // sun->get_position().scalar_mult(1000.0);
    // std::cout << "Sun Position: " << sun->get_position() << std::endl;
    my_sim.add_ray_source(sun);

    my_sim.add_element(dish);
    dish->update_geometry(sun_az, sun_el);

    // std::cout << std::setprecision(14);
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
    //               << "\nOrigin (ref): " << el->get_origin_ref()
    //               << "\nOrigin (stage): " << el->get_origin_stage()
    //               << "\nOrigin (global): " << el->get_origin_global()
    //               << "\nAim (ref): " << el->get_aim_vector_ref()
    //               << "\nAim (stage): " << el->get_aim_vector_stage()
    //               << "\nAim (global): " << el->get_aim_vector_global()
    //               << "\n";
    //     if (el->is_single())
    //     {
    //         std::cout << "\nAperture Type: " << el->get_aperture()->get_type()
    //                   << "\nAperture Diameter: " << el->get_aperture()->diameter_circumscribed_circle()
    //                   << "\nAperture Area: " << el->get_aperture()->aperture_area()
    //                   << "\n";
    //     }
    //     else if (el->is_composite())
    //     {
    //         dish = std::dynamic_pointer_cast<ParabolicDish>(el);
    //         Vector3d aim_loc;
    //         dish->convert_reference_to_local(aim_loc, dish->get_aim_vector_ref());
    //         std::cout << "\nElevation Axis: " << dish->get_elevation_axis()
    //                   << "\nAim (local): " << aim_loc
    //                   << "\nz-rotation: " << dish->get_zrot()
    //                   << "\nLocal to Ref: " << dish->get_local_to_reference()
    //                   << "\nLocal to Stage: " << dish->get_local_to_stage()
    //                   << "\nLocal to Global: " << dish->get_local_to_global()
    //                   << "\n";
    //     }
    // }

    // std::cout << "Aim vector: " << dish->get_aim_vector_global() << std::endl;

    RunnerStatus sts = my_runner.initialize();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    sts = my_runner.setup_simulation(&my_sim);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
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
    EXPECT_TRUE(num_absorbed >= N_ABSORBED_THRESH);
}

// Error Checking Tests for ParabolicDish
TEST(ParabolicDish, ErrorChecking_SetApertureSize)
{
    auto dish = SolTrace::Data::make_element<ParabolicDish>();

    // Test negative aperture size
    EXPECT_THROW(dish->set_aperture_size(-10.0), std::invalid_argument);
    EXPECT_THROW(dish->set_aperture_size(-1.0), std::invalid_argument);

    // Test zero aperture size
    EXPECT_THROW(dish->set_aperture_size(0.0), std::invalid_argument);

    // Test valid aperture sizes
    EXPECT_NO_THROW(dish->set_aperture_size(10.0));
    EXPECT_NO_THROW(dish->set_aperture_size(15.5));
}

TEST(ParabolicDish, ErrorChecking_SetFocalLength)
{
    auto dish = SolTrace::Data::make_element<ParabolicDish>();

    // Test negative focal length
    EXPECT_THROW(dish->set_focal_length(-7.5), std::invalid_argument);
    EXPECT_THROW(dish->set_focal_length(-0.1), std::invalid_argument);

    // Test zero focal length
    EXPECT_THROW(dish->set_focal_length(0.0), std::invalid_argument);

    // Test valid focal lengths
    EXPECT_NO_THROW(dish->set_focal_length(7.5));
    EXPECT_NO_THROW(dish->set_focal_length(15.0));
}

TEST(ParabolicDish, ErrorChecking_SetNumberPanels)
{
    auto dish = SolTrace::Data::make_element<ParabolicDish>();

    // Test invalid panel counts
    EXPECT_THROW(dish->set_number_of_panels(0, 2), std::invalid_argument);
    EXPECT_THROW(dish->set_number_of_panels(2, 0), std::invalid_argument);
    EXPECT_THROW(dish->set_number_of_panels(-1, 2), std::invalid_argument);
    EXPECT_THROW(dish->set_number_of_panels(2, -1), std::invalid_argument);

    // Test valid panel counts
    EXPECT_NO_THROW(dish->set_number_of_panels(1, 1));
    EXPECT_NO_THROW(dish->set_number_of_panels(2, 2));
    EXPECT_NO_THROW(dish->set_number_of_panels(10, 20));
}

TEST(ParabolicDish, ErrorChecking_SetGaps)
{
    auto dish = SolTrace::Data::make_element<ParabolicDish>();

    // Test negative gap values
    EXPECT_THROW(dish->set_gaps(-0.02, 0.01, 0.5), std::invalid_argument);
    EXPECT_THROW(dish->set_gaps(0.02, -0.01, 0.5), std::invalid_argument);
    // Note: center_radius can be negative (meaning no center gap)

    // Test valid gap values (including zero)
    EXPECT_NO_THROW(dish->set_gaps(0.0, 0.0, 0.0));
    EXPECT_NO_THROW(dish->set_gaps(0.02, 0.01, 0.5));
    EXPECT_NO_THROW(dish->set_gaps(0.02, 0.01, -1.0)); // Negative center gap is valid
}

TEST(ParabolicDish, ErrorChecking_SetReceiverDimensions)
{
    auto dish = SolTrace::Data::make_element<ParabolicDish>();

    // Test invalid receiver diameter
    EXPECT_THROW(dish->set_receiver_dimensions(0.0, 7.25), std::invalid_argument);
    EXPECT_THROW(dish->set_receiver_dimensions(-0.25, 7.25), std::invalid_argument);

    // Test invalid receiver distance
    EXPECT_THROW(dish->set_receiver_dimensions(0.25, -7.25), std::invalid_argument);

    // Test valid receiver dimensions
    EXPECT_NO_THROW(dish->set_receiver_dimensions(0.25, 0.0)); // Zero distance is valid
    EXPECT_NO_THROW(dish->set_receiver_dimensions(0.25, 7.25));
}

TEST(ParabolicDish, ErrorChecking_CreateGeometryWithoutParameters)
{
    auto dish = SolTrace::Data::make_element<ParabolicDish>();

    // Test create_geometry without setting required parameters
    EXPECT_THROW(dish->create_geometry(), std::invalid_argument);

    // Set aperture size and test again
    dish->set_aperture_size(10.0);
    EXPECT_THROW(dish->create_geometry(), std::invalid_argument);

    // Set focal length and test again
    dish->set_focal_length(7.5);
    EXPECT_THROW(dish->create_geometry(), std::invalid_argument);

    // Set number of panels and test again
    dish->set_number_of_panels(2, 2);
    EXPECT_THROW(dish->create_geometry(), std::invalid_argument);

    // Set receiver dimensions and it should work
    dish->set_receiver_dimensions(0.25, 7.25);
    EXPECT_NO_THROW(dish->create_geometry());
}

TEST(ParabolicDish, ErrorChecking_UpdateGeometry)
{
    auto dish = SolTrace::Data::make_element<ParabolicDish>();
    dish->set_aperture_size(10.0);
    dish->set_number_of_panels(2, 2);
    dish->set_gaps(0.02, 0.01, 0.5);
    dish->set_focal_length(7.5);
    dish->set_receiver_dimensions(0.25, 7.25);

    EXPECT_THROW(dish->update_geometry(20.0, 30.0), std::invalid_argument);

    dish->create_geometry();
    EXPECT_THROW(dish->update_geometry(20.0, -30.0), std::invalid_argument);

    EXPECT_NO_THROW(dish->update_geometry(20.0, 30.0));
}

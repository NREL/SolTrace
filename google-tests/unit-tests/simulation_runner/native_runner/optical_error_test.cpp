#include <gtest/gtest.h>

#include <optical_properties.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result_export.hpp>

#include <native_runner.hpp>

TEST(OpticalErrors, Gaussian)
{
    const uint_fast64_t NRAYS = 10000;

    using SolTrace::Runner::RunnerStatus;
    // Setup Runner
    SolTrace::NativeRunner::NativeRunner runner;
    RunnerStatus sts = runner.initialize();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    // Make default simulation data
    SimulationData sd;
    // Sun
    auto sun = make_ray_source<Sun>();
    sun->set_position(0, 0, 100);
    sd.add_ray_source(sun);

    // Make stage
    auto stage = make_stage(0);
    stage->set_origin(0, 0, 0);
    stage->set_aim_vector(0, 0, 1);
    stage->set_name("stage");

    // Make reflective flat plate
    auto plate = make_element<SingleElement>();
    plate->set_origin(0, 0, 0);
    plate->set_aim_vector(0, 0, 100); // Face up towards sun
    plate->set_surface(make_surface<Flat>());
    plate->set_aperture(make_aperture<Rectangle>(5, 5));
    plate->set_name("plate");

    SolTrace::Data::OpticalPropertySet plate_optics(SolTrace::Data::InteractionType::REFLECTION, "plate_optics_gaussian");
    plate_optics.set_ideal_reflection(SolTrace::Data::OpticalSide::Both);
    plate_optics.set_errors(SolTrace::Data::OpticalSide::Front, SolTrace::Data::DistributionType::GAUSSIAN, 1.0, 1e-3);
    plate_optics.set_errors(SolTrace::Data::OpticalSide::Back, SolTrace::Data::DistributionType::NONE, 0.0, 0.0);
    auto plate_optics_ref = sd.add_optical_property_set(plate_optics);
    plate->set_optical_property_set(plate_optics_ref);

    // Add element to stage
    stage->add_element(plate);

    // Add stage to sd
    sd.add_stage(stage);

    // Set parameters
    SimulationParameters &params = sd.get_simulation_parameters();
    params.number_of_rays = NRAYS;
    params.max_number_of_rays = NRAYS * 100;
    params.include_optical_errors = true;
    params.include_sun_shape_errors = false;
    params.seed = 123;

    // Run simulation with errors
    sts = runner.setup_simulation(&sd);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.run_simulation();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    // Collect results
    RayHistoryResult result_error;
    sts = runner.report_simulation(&result_error, 0);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    ASSERT_EQ(result_error.get_number_of_records(), NRAYS);

    element_id plate_id = plate->get_id();
    glm::dvec3 nhat(0.0, 0.0, 1.0);
    glm::dvec3 u;

    // auto it_ideal = result_ideal.get_ray_record_iterator();
    auto it_error = result_error.get_ray_record_iterator();

    while (!result_error.is_at_end(it_error))
    {
        auto err = *it_error;
        EXPECT_EQ(err->get_number_of_interactions(), 3);

        // The way this test is setup, all rays without errors come in
        // parallel (but opposite direction) to the normal of the plane and so
        // should bounce straight back--all departure directions, without
        // errors are khat = (0, 0, 1).

        EXPECT_EQ(err->get_element(1), plate_id);

        err->get_direction(1, u);
        // Extend u so that dot_product(nhat, u - nhat) == 0
        double alpha = 1.0 / glm::dot(nhat, u);
        u = alpha * u;
        // u = u - nhat
        u = u - nhat;

        // u is now the original perturbation vector. Do tests
        // on it.
        EXPECT_GT(glm::length(u), 0.0);

        // TODO: Devise some better statistical tests.

        ++it_error;
    }

    EXPECT_TRUE(result_error.is_at_end(it_error));
}

TEST(OpticalErrors, Uniform)
{
    const uint_fast64_t NRAYS = 10000;

    using SolTrace::Runner::RunnerStatus;
    // Setup Runner
    SolTrace::NativeRunner::NativeRunner runner;
    RunnerStatus sts = runner.initialize();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    // Make default simulation data
    SimulationData sd;
    // Sun
    auto sun = make_ray_source<Sun>();
    sun->set_position(0, 0, 100);
    sd.add_ray_source(sun);

    // Make stage
    auto stage = make_stage(0);
    stage->set_origin(0, 0, 0);
    stage->set_aim_vector(0, 0, 1);
    stage->set_name("stage");

    // Make reflective flat plate
    auto plate = make_element<SingleElement>();
    plate->set_origin(0, 0, 0);
    plate->set_aim_vector(0, 0, 100); // Face up towards sun
    plate->set_surface(make_surface<Flat>());
    plate->set_aperture(make_aperture<Rectangle>(5, 5));
    plate->set_name("plate");

    SolTrace::Data::OpticalPropertySet plate_optics(SolTrace::Data::InteractionType::REFLECTION, "plate_optics_uniform");
    plate_optics.set_ideal_reflection(SolTrace::Data::OpticalSide::Both);
    plate_optics.set_errors(SolTrace::Data::OpticalSide::Front, SolTrace::Data::DistributionType::PILLBOX, 1.0, 1e-3);
    plate_optics.set_errors(SolTrace::Data::OpticalSide::Back, SolTrace::Data::DistributionType::NONE, 0.0, 0.0);
    auto plate_optics_ref = sd.add_optical_property_set(plate_optics);
    plate->set_optical_property_set(plate_optics_ref);

    // Add element to stage
    stage->add_element(plate);

    // Add stage to sd
    sd.add_stage(stage);

    // Set parameters
    SimulationParameters &params = sd.get_simulation_parameters();
    params.number_of_rays = NRAYS;
    params.max_number_of_rays = NRAYS * 100;
    params.include_optical_errors = true;
    params.include_sun_shape_errors = false;
    params.seed = 123;

    // Run simulation with errors
    sts = runner.setup_simulation(&sd);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.run_simulation();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    // Collect results
    RayHistoryResult result_error;
    sts = runner.report_simulation(&result_error, 0);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    ASSERT_EQ(result_error.get_number_of_records(), NRAYS);

    element_id plate_id = plate->get_id();
    glm::dvec3 nhat(0.0, 0.0, 1.0);
    glm::dvec3 u;

    // auto it_ideal = result_ideal.get_ray_record_iterator();
    auto it_error = result_error.get_ray_record_iterator();

    while (!result_error.is_at_end(it_error))
    {
        auto err = *it_error;
        EXPECT_EQ(err->get_number_of_interactions(), 3);
        EXPECT_GT(err->get_number_of_interactions(), 2);

        // The way this test is setup, all rays without errors come in
        // parallel (but opposite direction) to the normal of the plane and so
        // should bounce straight back--all departure directions, without
        // errors are khat = (0, 0, 1).

        EXPECT_EQ(err->get_element(1), plate_id);

        err->get_direction(1, u);
        // Extend u so that dot_product(nhat, u - nhat) == 0
        double alpha = 1.0 / glm::dot(nhat, u);
        u = alpha * u;
        // u = u - nhat
        u = u - nhat;

        // u is now the original perturbation vector. Do tests
        // on it.
        EXPECT_GT(glm::length(u), 0.0);

        // TODO: Devise some better statistical tests.

        ++it_error;
    }

    EXPECT_TRUE(result_error.is_at_end(it_error));
}

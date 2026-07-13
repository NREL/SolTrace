#include <gtest/gtest.h>

#include <optical_properties.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result_export.hpp>

#include <optix_runner.hpp>

static OpticalPropertySetReference add_plate_optics(SimulationData& sd,
                                                  SolTrace::Data::DistributionType distribution)
{
    SolTrace::Data::OpticalPropertySet plate_optics(SolTrace::Data::InteractionType::REFLECTION,
        "plate_optics");
    plate_optics.set_ideal_reflection(OpticalSide::Both);
    plate_optics.set_errors(OpticalSide::Front, distribution, 1, 1e-3);

    return sd.add_optical_property_set(plate_optics);
}

TEST(OpticalErrors, Disabled)
{
    const uint_fast64_t NRAYS = 10000;

    using SolTrace::Runner::RunnerStatus;
    // Setup Runner
    OptixRunner runner;
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
    plate->set_optical_property_set(add_plate_optics(sd, SolTrace::Data::DistributionType::GAUSSIAN));

    // Add element to stage
    stage->add_element(plate);

    // Add stage to sd
    sd.add_stage(stage);

    // Set parameters
    SimulationParameters &params = sd.get_simulation_parameters();
    params.number_of_rays = NRAYS;
    params.max_number_of_rays = NRAYS * 100;
    params.include_optical_errors = false;
    params.include_sun_shape_errors = false;
    params.seed = 123;

    // Run simulation
    sts = runner.setup_simulation(&sd);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.run_simulation();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    // Collect results
    RayHistoryResult result;
    sts = runner.report_simulation(&result, 0);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    ASSERT_EQ(result.get_number_of_records(), NRAYS);

    element_id plate_id = plate->get_id();
    glm::dvec3 nhat(0.0, 0.0, 1.0);
    glm::dvec3 u;

    auto it = result.get_ray_record_iterator();

    while (!result.is_at_end(it))
    {
        auto rr = *it;
        EXPECT_GE(rr->get_number_of_interactions(), 2);

        // The way this test is setup, all rays without errors come in
        // parallel (but opposite direction) to the normal of the plane and so
        // should bounce straight back--all departure directions, without
        // errors are khat = (0, 0, 1).

        EXPECT_EQ(rr->get_element(1), plate_id);

        // // TODO: Need to get direction information in results before
        // // the below can be implemented.

        // err->get_direction(1, u);
        // // Extend u so that dot_product(nhat, u - nhat) == 0
        // double alpha = 1.0 / dot_product(nhat, u);
        // u.scalar_mult(alpha);
        // // u = u - nhat
        // vector_add(-1.0, nhat, 1.0, u);

        // // u is now the original perturbation vector. Do tests
        // // on it. It should be 0 since we errors are off.
	// EXPECT_NEAR(u.norm(), 1e-12);

        // // TODO: Devise some better statistical tests.

        ++it;
    }

    EXPECT_TRUE(result.is_at_end(it));
}


TEST(OpticalErrors, None)
{
    const uint_fast64_t NRAYS = 10000;

    using SolTrace::Runner::RunnerStatus;
    // Setup Runner
    OptixRunner runner;
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
    plate->set_optical_property_set(add_plate_optics(sd, SolTrace::Data::DistributionType::NONE));

    // Add element to stage
    stage->add_element(plate);

    // Add stage to sd
    sd.add_stage(stage);

    // Set parameters
    SimulationParameters &params = sd.get_simulation_parameters();
    params.number_of_rays = NRAYS;
    params.max_number_of_rays = NRAYS * 100;
    params.include_optical_errors = false;
    params.include_sun_shape_errors = false;
    params.seed = 123;

    // Run simulation
    sts = runner.setup_simulation(&sd);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.run_simulation();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    // Collect results
    RayHistoryResult result;
    sts = runner.report_simulation(&result, 0);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    ASSERT_EQ(result.get_number_of_records(), NRAYS);

    element_id plate_id = plate->get_id();
    glm::dvec3 nhat(0.0, 0.0, 1.0);
    glm::dvec3 u;

    auto it = result.get_ray_record_iterator();

    while (!result.is_at_end(it))
    {
        auto rr = *it;
        EXPECT_GE(rr->get_number_of_interactions(), 2);

        // The way this test is setup, all rays without errors come in
        // parallel (but opposite direction) to the normal of the plane and so
        // should bounce straight back--all departure directions, without
        // errors are khat = (0, 0, 1).

        EXPECT_EQ(rr->get_element(1), plate_id);

        // // TODO: Need to get direction information in results before
        // // the below can be implemented.

        // err->get_direction(1, u);
        // // Extend u so that dot_product(nhat, u - nhat) == 0
        // double alpha = 1.0 / dot_product(nhat, u);
        // u.scalar_mult(alpha);
        // // u = u - nhat
        // vector_add(-1.0, nhat, 1.0, u);

        // // u is now the original perturbation vector. Do tests
        // // on it. It should be 0 since we errors are off.
	// EXPECT_NEAR(u.norm(), 1e-12);

        // // TODO: Devise some better statistical tests.

        ++it;
    }

    EXPECT_TRUE(result.is_at_end(it));
}


TEST(OpticalErrors, Gaussian)
{
    const uint_fast64_t NRAYS = 10000;

    using SolTrace::Runner::RunnerStatus;
    // Setup Runner
    OptixRunner runner;
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
    plate->set_optical_property_set(add_plate_optics(sd, SolTrace::Data::DistributionType::GAUSSIAN));

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

    auto it_error = result_error.get_ray_record_iterator();

    while (!result_error.is_at_end(it_error))
    {
        auto err = *it_error;
        EXPECT_GE(err->get_number_of_interactions(), 2);

        // The way this test is setup, all rays without errors come in
        // parallel (but opposite direction) to the normal of the plane and so
        // should bounce straight back--all departure directions, without
        // errors are khat = (0, 0, 1).

        EXPECT_EQ(err->get_element(1), plate_id);

        // TODO: Need to get direction information in results before
        // the below can be implemented.

        // err->get_direction(1, u);
        // // Extend u so that dot_product(nhat, u - nhat) == 0
        // double alpha = 1.0 / dot_product(nhat, u);
        // u.scalar_mult(alpha);
        // // u = u - nhat
        // vector_add(-1.0, nhat, 1.0, u);

        // // u is now the original perturbation vector. Do tests
        // // on it.
        // EXPECT_GT(u.norm(), 0.0);

        // // TODO: Devise some better statistical tests.

        ++it_error;
    }

    EXPECT_TRUE(result_error.is_at_end(it_error));
}

TEST(OpticalErrors, PILLBOX)
{
    const uint_fast64_t NRAYS = 10000;

    using SolTrace::Runner::RunnerStatus;
    // Setup Runner
    OptixRunner runner;
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
    plate->set_optical_property_set(add_plate_optics(sd, SolTrace::Data::DistributionType::PILLBOX));

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

    auto it_error = result_error.get_ray_record_iterator();

    while (!result_error.is_at_end(it_error))
    {
        auto err = *it_error;
        EXPECT_GE(err->get_number_of_interactions(), 2);

        // The way this test is setup, all rays without errors come in
        // parallel (but opposite direction) to the normal of the plane and so
        // should bounce straight back--all departure directions, without
        // errors are khat = (0, 0, 1).

        EXPECT_EQ(err->get_element(1), plate_id);

        // TODO: Need to get direction information in results before
        // the below can be implemented.

        // err->get_direction(1, u);
        // // Extend u so that dot_product(nhat, u - nhat) == 0
        // double alpha = 1.0 / dot_product(nhat, u);
        // u.scalar_mult(alpha);
        // // u = u - nhat
        // vector_add(-1.0, nhat, 1.0, u);

        // // u is now the original perturbation vector. Do tests
        // // on it.
        // EXPECT_GT(u.norm(), 0.0);

        // // TODO: Devise some better statistical tests.

        ++it_error;
    }

    EXPECT_TRUE(result_error.is_at_end(it_error));
}

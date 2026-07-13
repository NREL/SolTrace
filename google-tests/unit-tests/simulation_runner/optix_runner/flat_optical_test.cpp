#include <gtest/gtest.h>

#include <optix_runner.hpp>
#include <optical_properties.hpp>
#include <simulation_data.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result_export.hpp>

using SolTrace::Runner::RunnerStatus;

using SolTrace::Data::OpticalPropertySet;

void make_default_sd(SimulationData& sd, element_ptr& plate)
{
	sd.clear();

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
	plate = make_element<SingleElement>();
	plate->set_origin(0, 0, 50);
	plate->set_aim_vector(0, 0, 100);	// Face up towards sun
	plate->set_surface(make_surface<Flat>());
	plate->set_aperture(make_aperture<Rectangle>(5, 5));
	InteractionType itype = InteractionType::REFLECTION;
	DistributionType dtype = DistributionType::NONE;	// No errors
	double transmissivity = 0;
	double reflectivity = 1;
	double slope_err = 0;	// Error not supported
	double spec_err = 0;
	double ri_front = 0;	// Refraction not supported
	double ri_back = 0;

	OpticalPropertySet plate_opt(itype, ri_front, ri_back, "PlateOptics");
	plate_opt.set_properties(OpticalSide::Both, dtype, transmissivity, reflectivity,
		slope_err, spec_err);
	auto plate_opt_ref = sd.add_optical_property_set(plate_opt);
	plate->set_optical_property_set(plate_opt_ref);

	plate->set_name("plate");

	// Add element to stage
	stage->add_element(plate);

	// Add stage to sd
	sd.add_stage(stage);

	// Set parameters
	SimulationParameters& params = sd.get_simulation_parameters();
	params.number_of_rays = 10000;
	params.max_number_of_rays = params.number_of_rays * 100;
	params.include_optical_errors = false;
	params.include_sun_shape_errors = false;
	params.seed = 123;
}

void count_hits(const RayHistoryResult& result,
	int& absorbed_count, int& transmitted_count,
	int& reflected_count)
{
	absorbed_count = 0;
	transmitted_count = 0;
	reflected_count = 0;
	int n_records = result.get_number_of_records();

	for (int i = 0; i < n_records; i++)
	{
		ray_record_ptr rec = result[i];

		int n_interactions = rec->get_number_of_interactions();
		for (int j = 0; j < n_interactions; j++)
		{
			RayEvent rev = rec->get_event(j);

			if (rev == RayEvent::ABSORB)
				absorbed_count++;
			else if (rev == RayEvent::TRANSMIT)
				transmitted_count++;
			else if (rev == RayEvent::REFLECT)
				reflected_count++;
		}
	}

	return;
}

TEST(FlatOptixOptical, Transmissivity)
{
	// Make default simulation data
	SimulationData sd;
	element_ptr plate;
	make_default_sd(sd, plate);

	// Set plate properties (transmissive)
	double transmissivity = 0.8;
	OpticalPropertySet* plate_opt_set = sd.get_mutable_optical_property_set(*plate);
	ASSERT_NE(plate_opt_set, nullptr);
	plate_opt_set->set_interaction_type(InteractionType::REFRACTION);
	plate_opt_set->set_transmissivity(OpticalSide::Both, transmissivity);

	// Run simulation
	OptixRunner runner;
	RunnerStatus sts = runner.initialize();
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);
	sts = runner.setup_simulation(&sd);
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);
	sts = runner.run_simulation();
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);

	// Collect results
	RayHistoryResult result;
	sts = runner.report_simulation(&result, 0);
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);

	// Calculate transmissivity
	int absorbed_count, transmitted_count, reflected_count;
	count_hits(result, absorbed_count, transmitted_count, reflected_count);
	int total_hits = absorbed_count + transmitted_count + reflected_count;
	double trans_calc = (double)transmitted_count / (double)total_hits;
	ASSERT_NEAR(trans_calc, transmissivity, 0.01);

}

TEST(FlatOptixOptical, Reflectivity)
{
	// Make default simulation data
	SimulationData sd;
	element_ptr plate;
	make_default_sd(sd, plate);
	
	// Set plate properties (reflective)
	double reflectivity = 0.6;
	OpticalPropertySet* plate_opt_set = sd.get_mutable_optical_property_set(*plate);
	ASSERT_NE(plate_opt_set, nullptr);
	plate_opt_set->set_interaction_type(InteractionType::REFLECTION);
	plate_opt_set->set_reflectivity(OpticalSide::Both, reflectivity);

	// Run simulation
	OptixRunner runner;
	RunnerStatus sts = runner.initialize();
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);
	sts = runner.setup_simulation(&sd);
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);
	sts = runner.run_simulation();
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);

	// Collect results
	RayHistoryResult result;
	sts = runner.report_simulation(&result, 0);
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);

	// Calculate reflectivity
	int absorbed_count, transmitted_count, reflected_count;
	count_hits(result, absorbed_count, transmitted_count, reflected_count);
	int total_hits = absorbed_count + transmitted_count + reflected_count;
	double refl_calc = (double)reflected_count / (double)total_hits;
	ASSERT_NEAR(refl_calc, reflectivity, 0.01);
}

TEST(FlatOptixOptical, SimResults)
{
	// Make default simulation data
	SimulationData sd;
	element_ptr plate;
	make_default_sd(sd, plate);

	// Set plate properties (reflective)
	double reflectivity = 0.6;
	OpticalPropertySet* plate_opt_set = sd.get_mutable_optical_property_set(*plate);
	ASSERT_NE(plate_opt_set, nullptr);
	plate_opt_set->set_interaction_type(InteractionType::REFLECTION);
	plate_opt_set->set_reflectivity(OpticalSide::Both, reflectivity);

	// Run simulation
	OptixRunner runner;
	RunnerStatus sts = runner.initialize();
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);
	sts = runner.setup_simulation(&sd);
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);
	sts = runner.run_simulation();
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);

	// Collect results
	RayHistoryResult result;
	sts = runner.report_simulation(&result, 0);
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);

	int n_records = result.get_number_of_records();

	// Check sim result order
	for (int i = 0; i < n_records; i++)
	{
		ray_record_ptr rec = result[i];

		int n_interactions = rec->get_number_of_interactions();
		RayEvent prev_rev = RayEvent::UNKNOWN;
		for (int j = 0; j < n_interactions; j++)
		{
			RayEvent rev = rec->get_event(j);

			// Check that ray event is assigned
			ASSERT_FALSE(rev == RayEvent::UNKNOWN);

			// Check that first interaction is create
			if (j == 0)
				ASSERT_TRUE(rev == RayEvent::CREATE);
			if (rev == RayEvent::CREATE)
				ASSERT_TRUE(j == 0);

			// Check that nothing happens after absorb or exit
			if (j > 0)
			{
				ASSERT_FALSE(prev_rev == RayEvent::ABSORB);
				ASSERT_FALSE(prev_rev == RayEvent::EXIT);
			}

			prev_rev = rev;
		}
	}
	
}

TEST(FlatOptixOptical, Absorption)
{
	// Make default simulation data
	SimulationData sd;
	element_ptr plate;
	make_default_sd(sd, plate);

	// Set plate properties (absorptive - no reflection, no transmission)
	OpticalPropertySet* plate_opt_set = sd.get_mutable_optical_property_set(*plate);
	ASSERT_NE(plate_opt_set, nullptr);
	plate_opt_set->set_ideal_absorption(OpticalSide::Both);

	// Run simulation
	OptixRunner runner;
	RunnerStatus sts = runner.initialize();
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);
	sts = runner.setup_simulation(&sd);
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);
	sts = runner.run_simulation();
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);

	// Collect results
	RayHistoryResult result;
	sts = runner.report_simulation(&result, 0);
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);

	// Calculate absorption - should be nearly 100%
	int absorbed_count, transmitted_count, reflected_count;
	count_hits(result, absorbed_count, transmitted_count, reflected_count);
	int total_hits = absorbed_count + transmitted_count + reflected_count;
	double abs_calc = (double)absorbed_count / (double)total_hits;
	ASSERT_DOUBLE_EQ(abs_calc, 1);
}

TEST(FlatOptixOptical, IdealReflection)
{
	// Make default simulation data
	SimulationData sd;
	element_ptr plate;
	make_default_sd(sd, plate);

	// Set plate properties (ideal reflection)
	OpticalPropertySet* plate_opt_set = sd.get_mutable_optical_property_set(*plate);
	ASSERT_NE(plate_opt_set, nullptr);
	plate_opt_set->set_ideal_reflection(OpticalSide::Both);

	// Run simulation
	OptixRunner runner;
	RunnerStatus sts = runner.initialize();
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);
	sts = runner.setup_simulation(&sd);
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);
	sts = runner.run_simulation();
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);

	// Collect results
	RayHistoryResult result;
	sts = runner.report_simulation(&result, 0);
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);

	// Calculate reflectivity - should be nearly 100%
	int absorbed_count, transmitted_count, reflected_count;
	count_hits(result, absorbed_count, transmitted_count, reflected_count);
	int total_hits = absorbed_count + transmitted_count + reflected_count;
	double refl_calc = (double)reflected_count / (double)total_hits;
	ASSERT_DOUBLE_EQ(refl_calc, 1.0);
}

TEST(FlatOptixOptical, SeedReproducibility)
{
	// Run simulation twice with same seed
	auto run_simulation = [](int seed) -> int {
		SimulationData sd;
		element_ptr plate;
		make_default_sd(sd, plate);

		// Set specific seed
		SimulationParameters& params = sd.get_simulation_parameters();
		params.seed = seed;

		// Set plate properties
		double reflectivity = 0.5;
		OpticalPropertySet* plate_opt_set = sd.get_mutable_optical_property_set(*plate);
		EXPECT_NE(plate_opt_set, nullptr);
		plate_opt_set->set_interaction_type(InteractionType::REFLECTION);
		plate_opt_set->set_reflectivity(OpticalSide::Both, reflectivity);

		OptixRunner runner;
		runner.initialize();
		runner.setup_simulation(&sd);
		runner.run_simulation();

		RayHistoryResult result;
		runner.report_simulation(&result, 0);

		int absorbed_count, transmitted_count, reflected_count;
		count_hits(result, absorbed_count, transmitted_count, reflected_count);
		return reflected_count;
	};

	int result1 = run_simulation(42);
	int result2 = run_simulation(42);
	int result3 = run_simulation(99);

	EXPECT_EQ(result1, result2);  // Same seed should give same result
	EXPECT_NE(result1, result3);  // Different seed should give different result
}

TEST(FlatOptixOptical, RayCountScaling)
{
	double reflectivity = 0.7;

	auto run_with_rays = [reflectivity](uint_fast64_t num_rays) -> double {
		SimulationData sd;
		element_ptr plate;
		make_default_sd(sd, plate);

		SimulationParameters& params = sd.get_simulation_parameters();
		params.number_of_rays = num_rays;
		params.max_number_of_rays = num_rays * 100;

		OpticalPropertySet* plate_opt_set = sd.get_mutable_optical_property_set(*plate);
		EXPECT_NE(plate_opt_set, nullptr);
		plate_opt_set->set_interaction_type(InteractionType::REFLECTION);
		plate_opt_set->set_reflectivity(OpticalSide::Both, reflectivity);

		OptixRunner runner;
		runner.initialize();
		runner.setup_simulation(&sd);
		runner.run_simulation();

		RayHistoryResult result;
		runner.report_simulation(&result, 0);

		int absorbed_count, transmitted_count, reflected_count;
		count_hits(result, absorbed_count, transmitted_count, reflected_count);
		int total_hits = absorbed_count + transmitted_count + reflected_count;
		return (double)reflected_count / (double)total_hits;
	};

	double result_1k = run_with_rays(1000);
	double result_10k = run_with_rays(10000);
	double result_50k = run_with_rays(50000);

	// All should be close to target reflectivity
	EXPECT_NEAR(result_1k, reflectivity, 0.05);   // Larger tolerance for fewer rays
	EXPECT_NEAR(result_10k, reflectivity, 0.02);
	EXPECT_NEAR(result_50k, reflectivity, 0.01);  // Tighter tolerance for more rays
}

TEST(FlatOptixOptical, MixedReflectionAbsorption)
{
	// Make default simulation data
	SimulationData sd;
	element_ptr plate;
	make_default_sd(sd, plate);

	// Set plate properties with partial reflection and absorption
	double reflectivity = 0.4;
	OpticalPropertySet* plate_opt_set = sd.get_mutable_optical_property_set(*plate);
	ASSERT_NE(plate_opt_set, nullptr);
	plate_opt_set->set_interaction_type(InteractionType::REFLECTION);
	plate_opt_set->set_reflectivity(OpticalSide::Both, reflectivity);

	// Run simulation
	OptixRunner runner;
	RunnerStatus sts = runner.initialize();
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);
	sts = runner.setup_simulation(&sd);
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);
	sts = runner.run_simulation();
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);

	// Collect results
	RayHistoryResult result;
	sts = runner.report_simulation(&result, 0);
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);

	// Verify both reflection and absorption occur
	int absorbed_count, transmitted_count, reflected_count;
	count_hits(result, absorbed_count, transmitted_count, reflected_count);
	int total_hits = absorbed_count + transmitted_count + reflected_count;

	double refl_calc = (double)reflected_count / (double)total_hits;
	double abs_calc = (double)absorbed_count / (double)total_hits;

	ASSERT_NEAR(refl_calc, reflectivity, 0.02);
	ASSERT_NEAR(abs_calc, 1.0 - reflectivity, 0.02);
}
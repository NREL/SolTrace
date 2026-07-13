#include <gtest/gtest.h>

#include <optical_properties.hpp>
#include <optix_runner.hpp>
#include <simulation_data.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result_export.hpp>

using SolTrace::Runner::RunnerStatus;

void make_two_plate_sd(SimulationData& sd, element_ptr& plate1, element_ptr& plate2)
{
	sd.clear();

	// Sun directly overhead
	auto sun = make_ray_source<Sun>();
	sun->set_position(0, 0, 100);
	sd.add_ray_source(sun);

	// Make stage
	auto stage = make_stage(0);
	stage->set_origin(0, 0, 0);
	stage->set_aim_vector(0, 0, 1);
	stage->set_name("stage");

	// Plate 1: Angled 45 degrees, receives rays from sun and reflects toward plate 2
	// Position at origin, tilted to reflect rays in +Y direction
	plate1 = make_element<SingleElement>();
	plate1->set_origin(0, 0, 50);
	plate1->set_aim_vector(0, 50, 100);  // Tilted 45 degrees toward +Y
	plate1->set_surface(make_surface<Flat>());
	plate1->set_aperture(make_aperture<Rectangle>(5, 5));
	SolTrace::Data::OpticalPropertySet reflective_optics(
		SolTrace::Data::InteractionType::REFLECTION,
		0.0,
		0.0,
		"two_plate_reflector_optics");
	reflective_optics.set_ideal_reflection(OpticalSide::Both);
	auto reflective_optics_ref = sd.add_optical_property_set(reflective_optics);
	plate1->set_optical_property_set(reflective_optics_ref);
	plate1->set_name("plate1");

	// Plate 2: Positioned to receive reflected rays from plate 1
	// At Y=50, tilted 45 degrees to face plate 1
	plate2 = make_element<SingleElement>();
	plate2->set_origin(0, 50, 50);
	plate2->set_aim_vector(0, 0, 100);  // Tilted 45 degrees toward -Y (facing plate 1)
	plate2->set_surface(make_surface<Flat>());
	plate2->set_aperture(make_aperture<Rectangle>(10, 10));  // Larger to catch reflected rays
	SolTrace::Data::OpticalPropertySet absorbing_optics(
		SolTrace::Data::InteractionType::REFLECTION,
		0.0,
		0.0,
		"two_plate_absorber_optics");
	absorbing_optics.set_ideal_absorption(OpticalSide::Both);
	auto absorbing_optics_ref = sd.add_optical_property_set(absorbing_optics);
	plate2->set_optical_property_set(absorbing_optics_ref);
	plate2->set_name("plate2");

	// Add elements to stage
	stage->add_element(plate1);
	stage->add_element(plate2);

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

void count_hits_two_plate(const RayHistoryResult& result,
	int plate1_id, int plate2_id,
	int& absorbed_plate1, int& transmitted_plate1, int& reflected_plate1,
	int& absorbed_plate2, int& transmitted_plate2, int& reflected_plate2,
	int& absorbed_total, int& transmitted_total, int& reflected_total)
{
	absorbed_plate1 = transmitted_plate1 = reflected_plate1 = 0;
	absorbed_plate2 = transmitted_plate2 = reflected_plate2 = 0;
	absorbed_total = transmitted_total = reflected_total = 0;

	int n_records = result.get_number_of_records();
	for (int i = 0; i < n_records; i++)
	{
		ray_record_ptr rec = result[i];
		int n_interactions = rec->get_number_of_interactions();
		for (int j = 0; j < n_interactions; j++)
		{
			RayEvent rev = rec->get_event(j);
			int elid = static_cast<int>(rec->get_element(j));

			// Totals
			if (rev == RayEvent::ABSORB) absorbed_total++;
			else if (rev == RayEvent::TRANSMIT) transmitted_total++;
			else if (rev == RayEvent::REFLECT) reflected_total++;

			// Per-plate
			if (elid == plate1_id)
			{
				if (rev == RayEvent::ABSORB) absorbed_plate1++;
				else if (rev == RayEvent::TRANSMIT) transmitted_plate1++;
				else if (rev == RayEvent::REFLECT) reflected_plate1++;
			}
			else if (elid == plate2_id)
			{
				if (rev == RayEvent::ABSORB) absorbed_plate2++;
				else if (rev == RayEvent::TRANSMIT) transmitted_plate2++;
				else if (rev == RayEvent::REFLECT) reflected_plate2++;
			}
		}
	}
}

TEST(TwoPlateOptix, ReflectionToAbsorber)
{
	// Make two-plate simulation data
	SimulationData sd;
	element_ptr plate1, plate2;
	make_two_plate_sd(sd, plate1, plate2);

	// Run simulation
	OptixRunner runner;
	RunnerStatus sts = runner.initialize();
	ASSERT_EQ(sts, RunnerStatus::SUCCESS);
	sts = runner.setup_simulation(&sd);
	ASSERT_EQ(sts, RunnerStatus::SUCCESS);
	sts = runner.run_simulation();
	ASSERT_EQ(sts, RunnerStatus::SUCCESS);

	// Collect results
	RayHistoryResult result;
	sts = runner.report_simulation(&result, 0);
	ASSERT_EQ(sts, RunnerStatus::SUCCESS);

	// Count events per-plate and totals
	int a1, t1, r1, a2, t2, r2, at, tt, rt;
	count_hits_two_plate(result, plate1->get_id(), plate2->get_id(), a1, t1, r1, a2, t2, r2, at, tt, rt);

	// Check number of hits
	EXPECT_EQ(sd.get_simulation_parameters().number_of_rays, result.get_number_of_records());

	// Rays should reflect off plate1 and be absorbed by plate2
	EXPECT_EQ(a1, 0);	// plate1 is ideal reflective
	EXPECT_EQ(r2, 0);	// plate2 is ideal absorptive
	EXPECT_GT(r1, 0);
	EXPECT_GT(a2, 0);

	// All the rays hit something
	EXPECT_GT(rt + at, sd.get_simulation_parameters().number_of_rays);

	// More absorptions than reflections
	EXPECT_GT(at, rt);

	// Sun Ray Checks
	int N_sun_rays = runner.get_N_sun_rays();
	OptixCSP::SolTraceSystem* sys = runner.get_optix_system();
	
	std::vector<float4> hp_vec;
	std::vector<uint_fast64_t> raynumber_vec;
	std::vector<int32_t> element_id_vec;
	std::vector<uint8_t> hit_type_vec;
	sys->get_hp_output(hp_vec, raynumber_vec, element_id_vec, hit_type_vec);
	// std::vector<uint_fast64_t> sunraynumber_vec = sys->get_sunraynumber_vec();

	EXPECT_EQ(hp_vec.size(), raynumber_vec.size());						// Hit results are same size
	EXPECT_EQ(raynumber_vec.size(), element_id_vec.size());
	EXPECT_EQ(element_id_vec.size(), hit_type_vec.size());
	// EXPECT_EQ(N_sun_rays, sunraynumber_vec.back());						// Reported sun rays is the sun ray id of last hit
	EXPECT_TRUE(N_sun_rays <= sd.get_simulation_parameters().max_number_of_rays);	// Only generated max number of rays or fewer
}

TEST(TwoPlateOptix, BatchMaxRayLimit)
{
	// Make two-plate simulation data
	SimulationData sd;
	element_ptr plate1, plate2;
	make_two_plate_sd(sd, plate1, plate2);

	// Set max rays equal to desired
	int NRays = 1000;
	int NMaxRays = NRays;
	auto& sim_par = sd.get_simulation_parameters();
	sim_par.number_of_rays = NRays;
	sim_par.max_number_of_rays = NMaxRays;

	// Run simulation
	OptixRunner runner;
	RunnerStatus sts = runner.initialize();
	ASSERT_EQ(sts, RunnerStatus::SUCCESS);
	sts = runner.setup_simulation(&sd);
	ASSERT_EQ(sts, RunnerStatus::SUCCESS);
	sts = runner.run_simulation();
	ASSERT_EQ(sts, RunnerStatus::SUCCESS);

	// Collect results
	RayHistoryResult result;
	sts = runner.report_simulation(&result, 0);
	EXPECT_EQ(sts, RunnerStatus::SUCCESS);

	// Should NOT reach desired number of rays
	ASSERT_TRUE(result.get_number_of_records() < NRays);
}

TEST(TwoPlateOptix, SimResults)
{
	// Make two-plate simulation data
	SimulationData sd;
	element_ptr plate1, plate2;
	make_two_plate_sd(sd, plate1, plate2);

	// Run simulation
	OptixRunner runner;
	RunnerStatus sts = runner.initialize();
	ASSERT_EQ(sts, RunnerStatus::SUCCESS);
	sts = runner.setup_simulation(&sd);
	ASSERT_EQ(sts, RunnerStatus::SUCCESS);
	sts = runner.run_simulation();
	ASSERT_EQ(sts, RunnerStatus::SUCCESS);

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

TEST(TwoPlateOptix, TrimExcessRaysOption)
{
	SimulationData sd;
	element_ptr plate1, plate2;
	make_two_plate_sd(sd, plate1, plate2);
	const int n_rays = sd.get_simulation_parameters().number_of_rays;

	// Default: trim enabled — result has exactly n_rays records
	{
		OptixRunner runner;
		EXPECT_TRUE(runner.get_trim_excess_rays());  // default is true

		ASSERT_EQ(runner.initialize(), RunnerStatus::SUCCESS);
		ASSERT_EQ(runner.setup_simulation(&sd), RunnerStatus::SUCCESS);
		ASSERT_EQ(runner.run_simulation(), RunnerStatus::SUCCESS);

		RayHistoryResult result;
		ASSERT_EQ(runner.report_simulation(&result, 0), RunnerStatus::SUCCESS);
		EXPECT_EQ(result.get_number_of_records(), n_rays);
	}

	// Trim disabled — result has at least n_rays records (batch overshoot is not removed)
	{
		OptixRunner runner;
		runner.set_trim_excess_rays(false);
		EXPECT_FALSE(runner.get_trim_excess_rays());

		ASSERT_EQ(runner.initialize(), RunnerStatus::SUCCESS);
		ASSERT_EQ(runner.setup_simulation(&sd), RunnerStatus::SUCCESS);
		ASSERT_EQ(runner.run_simulation(), RunnerStatus::SUCCESS);

		RayHistoryResult result;
		ASSERT_EQ(runner.report_simulation(&result, 0), RunnerStatus::SUCCESS);
		EXPECT_GE(result.get_number_of_records(), n_rays);
	}
}

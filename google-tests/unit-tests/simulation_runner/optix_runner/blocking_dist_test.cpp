#include <gtest/gtest.h>

#include <composite_element.hpp>
#include <constants.hpp>
#include <element.hpp>
#include <optical_properties.hpp>
#include <optix_runner.hpp>
#include <sun.hpp>
#include <simulation_data.hpp>
#include <simulation_data_export.hpp>
#include <single_element.hpp>
#include <stage_element.hpp>
#include <simulation_result_export.hpp>

#include <cmath>
#include <fstream>
#include <map>
#include <iostream>
#include <sstream>

using SolTrace::Runner::RunnerStatus;


static void setup_tower_sd(SimulationData& sd, element_ptr& reflector, element_ptr& absorber, element_ptr& receiver, uint_fast64_t nrays, DistributionType dist)
{

    auto sun = make_ray_source<Sun>();
    sun->set_position(0.0, 0.0, 1000.0);
    sun->set_shape(SolTrace::Data::SunShape::PILLBOX, 0, 4.65, 0);
    sun->set_gen_type(SolTrace::Data::GenType::HALTON);
    sd.add_ray_source(sun);

    auto stage = make_stage(0);
	stage->set_origin(0, 0, 0);
	stage->set_aim_vector(0, 0, 1);
	stage->set_name("stage");

    absorber = make_element<SingleElement>();
    absorber->set_origin(0.0, 500.0, 5.65);
    absorber->set_aim_vector(0.0, 0.0, 169.0);
    absorber->set_surface(make_surface<Flat>());
    absorber->set_aperture(make_aperture<Rectangle>(11.415, 10.42));
    SolTrace::Data::OpticalPropertySet absorber_optics(InteractionType::REFLECTION, "Absorber");
    absorber_optics.set_ideal_absorption(OpticalSide::Both);
    absorber_optics.set_errors(SolTrace::Data::OpticalSide::Both, dist, 2.0, 0.0);
    auto absorber_optics_ref = sd.add_optical_property_set(absorber_optics);
    absorber->set_optical_property_set(absorber_optics_ref);
    absorber->enable();

    reflector = make_element<SingleElement>();
    reflector->set_origin(0.0, 700.0, 5.65);
    reflector->set_aim_vector(0.0, 0.0, 169.0);
    reflector->set_surface(make_surface<Flat>());
    reflector->set_aperture(make_aperture<Rectangle>(11.415, 10.42));
    SolTrace::Data::OpticalPropertySet reflector_optics(InteractionType::REFLECTION, "Reflector");
    reflector_optics.set_ideal_absorption(OpticalSide::Back);
    reflector_optics.set_ideal_reflection(OpticalSide::Front);
    reflector_optics.set_errors(SolTrace::Data::OpticalSide::Back, dist, 2.0, 0.0);
    reflector_optics.set_errors(SolTrace::Data::OpticalSide::Front, SolTrace::Data::DistributionType::GAUSSIAN, 2.0, 0.0);
    auto reflector_optics_ref = sd.add_optical_property_set(reflector_optics);
    reflector->set_optical_property_set(reflector_optics_ref);
    reflector->enable();

    SolTrace::Data::OpticalPropertySet rec_opt_set(SolTrace::Data::InteractionType::REFLECTION, "Receiver");
    rec_opt_set.set_ideal_absorption(SolTrace::Data::OpticalSide::Both);
    auto rec_ref = sd.add_optical_property_set(rec_opt_set);
    receiver = SolTrace::Data::make_element<SingleElement>();
    receiver->set_optical_property_set(rec_ref);
    receiver->set_aperture(SolTrace::Data::make_aperture<SolTrace::Data::Rectangle>(12.0, 18.0));
    receiver->set_surface(SolTrace::Data::make_surface<SolTrace::Data::Flat>());
    glm::dvec3 v1 = {1.0, 0.0, 0.0}; // Pointing North TODO: change to point towards heliostat
    glm::dvec3 rec_origin(0.0, 0.0, 169.0);
    glm::dvec3 aim_point = rec_origin + v1;
    receiver->set_reference_frame_geometry(rec_origin, aim_point, 0.0);
    receiver->set_name("Receiver");
    receiver->enable();

    // Add elements to stage
	stage->add_element(absorber);
    stage->add_element(reflector);
	stage->add_element(receiver);

	// Add stage to sd
	sd.add_stage(stage);

    SimulationParameters& params = sd.get_simulation_parameters();
	params.number_of_rays = nrays;
	params.max_number_of_rays = params.number_of_rays * 100;
	params.include_optical_errors = false;
	params.include_sun_shape_errors = false;
	params.seed = 123;

}


void calculate_ray_counts(SimulationResult result, int reflector_id, int absorber_id, int receiver_id,
    uint_fast64_t& helio_hit_count, uint_fast64_t& reflect_count, uint_fast64_t& helio_absorb_count,
    uint_fast64_t& rec_absorb_count, uint_fast64_t& miss_count, uint_fast64_t& rec_hit_count,
    uint_fast64_t& rec_direct_hit_count, uint_fast64_t& rec_via_helio_hit_count) 
{
    // Reset counts
    helio_hit_count = 0;
    reflect_count = 0;
    helio_absorb_count = 0;
    rec_absorb_count = 0;
    miss_count = 0;
      
    rec_hit_count = 0;
    rec_direct_hit_count = 0;
    rec_via_helio_hit_count = 0;

    for (size_t i = 0; i < result.get_number_of_records(); i++) {
        const ray_record_ptr rr = result[i];

        for (size_t j = 0; j < rr->interactions.size(); j++) {
            auto hit_element = rr->get_element(j);
            SolTrace::Result::RayEvent rev = rr->get_event(j);

            if (rev == RayEvent::EXIT) miss_count++;
            if ((int)rev <= (int)RayEvent::CREATE || (int)rev >= (int)RayEvent::EXIT) continue;  // create or exit

            // Check heliostat elements
            if (hit_element == absorber_id || hit_element == reflector_id) {
                helio_hit_count++;
                if (rev == RayEvent::REFLECT) reflect_count++;
                if (rev == RayEvent::ABSORB) helio_absorb_count++;

            }

            // Check receiver element
            if (hit_element == receiver_id) {
                rec_hit_count++;

                // Check order of hit
                if (j == 1)
                    rec_direct_hit_count++;
                else
                    rec_via_helio_hit_count++;

                // Check absorbed
                if (rev == RayEvent::ABSORB) rec_absorb_count++;
            }
        }
    }

}


void save_outputs(std::string filename, const std::map<std::string, double>& results)
{
    std::ofstream outputFile(filename, std::ios::out | std::ios::trunc);
    if (!outputFile.is_open())
    {
        std::cerr << "Error: Could not open the file " << filename << std::endl;
        return;
    }

    for (const auto& kv : results)
    {
        outputFile << kv.first << "," << kv.second << std::endl;
    }

    outputFile.close();
}

static void write_to_dict(std::string key_name, double val_a,
	double val_b, std::map<std::string, double>& dict_a,
	std::map<std::string, double>& dict_b)
{
	dict_a[key_name] = val_a;
	dict_b[key_name] = val_b;
}

TEST(FlatOptixOptical, SlopeErrorDistAbsorption)
{
    const uint_fast64_t NRAYS = 2e6;

    element_ptr g_absorber, g_receiver, g_reflector, n_reflector, n_absorber, n_receiver;
    
    SimulationData gauss;
    setup_tower_sd(gauss, g_reflector, g_absorber, g_receiver, NRAYS, DistributionType::GAUSSIAN);
    OptixRunner g_runner;
    RunnerStatus g_sts = g_runner.initialize();
    ASSERT_EQ(g_sts, RunnerStatus::SUCCESS);
    g_sts = g_runner.setup_simulation(&gauss);
    ASSERT_EQ(g_sts, RunnerStatus::SUCCESS);
    g_sts = g_runner.run_simulation();
	ASSERT_EQ(g_sts, RunnerStatus::SUCCESS);
    SimulationResult g_result;
	g_sts = g_runner.report_simulation(&g_result, 0);
    ASSERT_EQ(g_sts, RunnerStatus::SUCCESS);

    SimulationData none;
    setup_tower_sd(none, n_reflector, n_absorber, n_receiver, NRAYS, DistributionType::NONE);
    OptixRunner n_runner;
    RunnerStatus n_sts = n_runner.initialize();
    ASSERT_EQ(n_sts, RunnerStatus::SUCCESS);
    n_sts = n_runner.setup_simulation(&none);
    ASSERT_EQ(n_sts, RunnerStatus::SUCCESS);
    n_sts = n_runner.run_simulation();
	ASSERT_EQ(n_sts, RunnerStatus::SUCCESS);
    SimulationResult n_result;
	n_sts = n_runner.report_simulation(&n_result, 0);
    ASSERT_EQ(n_sts, RunnerStatus::SUCCESS);

    uint_fast64_t g_helio_hit_count, g_reflect_count, g_helio_absorb_count, g_rec_absorb_count, g_miss_count;
    uint_fast64_t g_rec_hit_count, g_rec_direct_hit_count, g_rec_via_helio_hit_count;
    calculate_ray_counts(g_result, g_reflector->get_id(), g_absorber->get_id(), g_receiver->get_id(),
        g_helio_hit_count, g_reflect_count, g_helio_absorb_count,
        g_rec_absorb_count, g_miss_count, g_rec_hit_count,
        g_rec_direct_hit_count, g_rec_via_helio_hit_count);

    uint_fast64_t n_helio_hit_count, n_reflect_count, n_helio_absorb_count, n_rec_absorb_count, n_miss_count;
    uint_fast64_t n_rec_hit_count, n_rec_direct_hit_count, n_rec_via_helio_hit_count;
    calculate_ray_counts(n_result, n_reflector->get_id(), n_absorber->get_id(), n_receiver->get_id(),
        n_helio_hit_count, n_reflect_count, n_helio_absorb_count,
        n_rec_absorb_count, n_miss_count, n_rec_hit_count,
        n_rec_direct_hit_count, n_rec_via_helio_hit_count);

    std::map<std::string, double> dict_none;
	std::map<std::string, double> dict_gauss;

    write_to_dict("00_tot_helio_hits", n_helio_hit_count, g_helio_hit_count, dict_none, dict_gauss);
	write_to_dict("01_tot_helio_absorb_count", n_helio_absorb_count, g_helio_absorb_count, dict_none, dict_gauss);
	write_to_dict("02_tot_reflect_count", n_reflect_count, g_reflect_count, dict_none, dict_gauss);
	write_to_dict("03_rec_absorb_count", n_rec_absorb_count, g_rec_absorb_count, dict_none, dict_gauss);
	write_to_dict("05_miss_count", n_miss_count, g_miss_count, dict_none, dict_gauss);
	write_to_dict("06_tot_rec_hits", n_rec_hit_count, g_rec_hit_count, dict_none, dict_gauss);
	write_to_dict("07_rec_direct_count", n_rec_direct_hit_count, g_rec_direct_hit_count, dict_none, dict_gauss);
	write_to_dict("08_rec_via_helio_count", n_rec_via_helio_hit_count, g_rec_via_helio_hit_count, dict_none, dict_gauss);

    std::string file_outputs_gauss = "gauss_outputs.csv";
	save_outputs(file_outputs_gauss, dict_gauss);

    std::string file_outputs_none = "none_outputs.csv";
	save_outputs(file_outputs_none, dict_none);

    EXPECT_EQ(g_helio_hit_count, n_helio_hit_count);
    EXPECT_EQ(g_reflect_count, n_reflect_count);
    EXPECT_EQ(g_helio_absorb_count, n_helio_absorb_count);
    EXPECT_EQ(g_rec_absorb_count, n_rec_absorb_count);
    EXPECT_EQ(g_miss_count, n_miss_count);
    EXPECT_EQ(g_rec_hit_count, n_rec_hit_count);
    EXPECT_EQ(g_rec_direct_hit_count, n_rec_direct_hit_count);
    EXPECT_EQ(g_rec_via_helio_hit_count, n_rec_via_helio_hit_count);

}
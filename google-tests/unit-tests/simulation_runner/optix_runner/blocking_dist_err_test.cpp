#pragma once

#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include <simulation_data.hpp>
#include <simulation_result_export.hpp>
#include <stage_element.hpp>
#include <sun.hpp>
#include <utilities.hpp>

#include <../../hpvm.h>

#include <cst_templates/heliostat.hpp>

#include "common.hpp"

#include <optix_runner.hpp>

using Heliostat = SolTrace::Data::Heliostat;
using SolTrace::Runner::RunnerStatus;

void setup_scene(SimulationData &simData, DistributionType blocking_dist, element_ptr& receiver){

    simData.clear();

    auto sun = make_ray_source<Sun>();
    glm::dvec3 sun_pos = {0.0, 0.0, 1000.0};
    sun = SolTrace::Data::make_ray_source<Sun>();
    sun->set_position(sun_pos);
    sun->set_shape(SolTrace::Data::SunShape::PILLBOX, 0.0, 4.65, 0.0);
    sun->set_gen_type(SolTrace::Data::GenType::HALTON);
    simData.add_ray_source(sun);

    const glm::dvec3 zero = {0.0, 0.0, 0.0}; // Global origin
    const glm::dvec3 khat = {0.0, 0.0, 1.0}; // Global z-axis

    glm::dvec3 rec_origin = {0.0, 0.0, 180.}; // This is the receiver center
    receiver = SolTrace::Data::make_element<SingleElement>();
    SolTrace::Data::OpticalPropertySet receiver_opt_set(SolTrace::Data::InteractionType::REFLECTION, "ReceiverOptics");
    receiver_opt_set.set_ideal_absorption(SolTrace::Data::OpticalSide::Both);
    auto receiver_ref = simData.add_optical_property_set(receiver_opt_set);
    receiver->set_optical_property_set(receiver_ref);
    receiver->set_aperture(SolTrace::Data::make_aperture<SolTrace::Data::Rectangle>(15, 25));
    receiver->set_surface(SolTrace::Data::make_surface<SolTrace::Data::Cylinder>(7.5));
    glm::dvec3 offset = {0.0, 0.0, 0.0};
    glm::dvec3 rec_origin_offset = rec_origin + offset;
    glm::dvec3 v1 = {0.0, -1.0, 0.0};
    glm::dvec3 aim_point = rec_origin_offset + v1;
    receiver->set_reference_frame_geometry(rec_origin_offset, aim_point, 180.0);
    receiver->set_name("Receiver");
    receiver->enable();

    SolTrace::Data::OpticalPropertySet mirror_opt_set(SolTrace::Data::InteractionType::REFLECTION, "BlockingHeliostatMirrorOptics");
    mirror_opt_set.set_ideal_absorption(SolTrace::Data::OpticalSide::Both);
    mirror_opt_set.set_errors(SolTrace::Data::OpticalSide::Both, blocking_dist, 0.0, 0.0);
    auto blocking_mirror_ref = simData.add_optical_property_set(mirror_opt_set);
    auto heliostat = SolTrace::Data::make_element<Heliostat>();
    heliostat = SolTrace::Data::make_element<Heliostat>();
    heliostat->set_optics(blocking_mirror_ref);
    glm::dvec3 heliostat_origin = {0, 303.45, 4.49}; //add heliostat
    heliostat->set_reference_frame_geometry(heliostat_origin, khat, 0.0);
    heliostat->set_aperture_size(10.0, 10.0);   // Width, Height
    heliostat->set_number_panels(1, 1);
    heliostat->set_gaps(0, 0);
    heliostat->set_canting(Heliostat::NONE, 0.0, 0.0);
    heliostat->set_target_position(rec_origin);
    heliostat->set_focal_length(0.0);
    heliostat->create_geometry();
    heliostat->set_name("Heliostat");
    heliostat->enable();

    SolTrace::Data::sun_position_vector_degrees(sun_pos, 74.95, 26.26);
    sun->set_position(sun_pos);
    heliostat->update_geometry(74.95, 26.26);




    stage_ptr stage = SolTrace::Data::make_stage(0);
    stage->set_reference_frame_geometry(zero, khat, 0.0);

    auto ret = stage->add_element(heliostat);
    EXPECT_TRUE(SolTrace::Data::Element::is_success(ret));

    ret = stage->add_element(receiver);
    EXPECT_TRUE(SolTrace::Data::Element::is_success(ret));

    simData.add_stage(stage);

    SimulationParameters& params = simData.get_simulation_parameters();
	params.number_of_rays = 2e6;
	params.max_number_of_rays = params.number_of_rays * 100;
	params.include_optical_errors = true;
	params.include_sun_shape_errors = true;
	params.seed = 123;

}

void count_hits(const SimulationResult& result, int rec_id,
	uint_fast64_t& helio_hit_count, uint_fast64_t& rec_hit_count)
{
	helio_hit_count = 0;
	rec_hit_count = 0;
	int n_records = result.get_number_of_records();

	for (int i = 0; i < n_records; i++)
	{
		ray_record_ptr rec = result[i];

		int n_interactions = rec->get_number_of_interactions();
		for (int j = 0; j < n_interactions; j++)
		{
            auto hit_element = rec->get_element(j);
			RayEvent rev = rec->get_event(j);

			if (hit_element == rec_id)
				rec_hit_count++;
			else if(rev == RayEvent::ABSORB || rev == RayEvent::REFLECT)
				helio_hit_count++;
		}
	}

	return;
}


TEST(FlatOptixOptical, BlockingErrorDistribution){

    SimulationData g_sd;
    SimulationData n_sd;

    element_ptr g_receiver, n_receiver;

    setup_scene(g_sd, SolTrace::Data::DistributionType::GAUSSIAN, g_receiver);
    setup_scene(n_sd, SolTrace::Data::DistributionType::NONE, n_receiver);

    OptixRunner g_runner;
	RunnerStatus g_sts = g_runner.initialize();
	EXPECT_EQ(g_sts, RunnerStatus::SUCCESS);
	g_sts = g_runner.setup_simulation(&g_sd);
	EXPECT_EQ(g_sts, RunnerStatus::SUCCESS);
	g_sts = g_runner.run_simulation();
	EXPECT_EQ(g_sts, RunnerStatus::SUCCESS);

	// Collect results
	SimulationResult g_result;
	g_sts = g_runner.report_simulation(&g_result, 0);
	EXPECT_EQ(g_sts, RunnerStatus::SUCCESS);

    OptixRunner n_runner;
	RunnerStatus n_sts = n_runner.initialize();
	EXPECT_EQ(n_sts, RunnerStatus::SUCCESS);
	n_sts = n_runner.setup_simulation(&n_sd);
	EXPECT_EQ(n_sts, RunnerStatus::SUCCESS);
	n_sts = n_runner.run_simulation();
	EXPECT_EQ(n_sts, RunnerStatus::SUCCESS);

	// Collect results
	SimulationResult n_result;
	n_sts = n_runner.report_simulation(&n_result, 0);
	EXPECT_EQ(n_sts, RunnerStatus::SUCCESS);

    uint_fast64_t gauss_helio_hit_count, gauss_rec_hit_count, none_helio_hit_count, none_rec_hit_count;

    count_hits(g_result, g_receiver->get_id(), gauss_helio_hit_count, gauss_rec_hit_count);
    count_hits(n_result, n_receiver->get_id(), none_helio_hit_count, none_rec_hit_count);

    EXPECT_EQ(gauss_helio_hit_count,none_helio_hit_count);
    EXPECT_EQ(gauss_rec_hit_count, none_rec_hit_count);

}
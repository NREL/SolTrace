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
#include "count_absorbed_native.h"

using Heliostat = SolTrace::Data::Heliostat;
using SolTrace::Runner::RunnerStatus;

static void save_hit_pos_to_file(const SimulationResult& result, std::string filename)
{
    std::ofstream outputFile(filename, std::ios::out | std::ios::trunc);
    if (!outputFile.is_open()) {
        std::cerr << "Error: Could not open the file " << filename << std::endl;
    }
    for (size_t i = 0; i < result.get_number_of_records(); i++) {
        const ray_record_ptr rr = result[i];
        for (size_t j = 0; j < rr->interactions.size(); j++) {
            SolTrace::Result::RayEvent rev = rr->get_event(j);
            if (j == 1)
            {
                glm::dvec3 pos;
                rr->get_position(j, pos);
                outputFile << pos[0] << "," << pos[1] << "," << pos[2];
                outputFile << std::endl;
            }
        }
    }
}

template <typename RunnerT>
class HeliostatFieldSimulationHelper {
public:
    bool high_accuracy = false;     // Runs 20 Million rays and tighter tolerance on checks
    bool print_info = false;        // Prints information on from simulation results (sun calculations, ray counts, flux calculations)
    bool save_results = false;      // Saves flux map results to CSV files
    bool save_raydata = false;      // Saves ray data to CSV file

    int seed = 123;
    SolTrace::Data::GenType sun_gen_type = SolTrace::Data::GenType::HALTON;
    bool use_optical_errors = true;
    bool use_sunshape_errors = true;

    const glm::dvec3 zero = {0.0, 0.0, 0.0}; // Global origin
    const glm::dvec3 khat = {0.0, 0.0, 1.0}; // Global z-axis

    // Receiver parameters
    glm::dvec3 rec_origin = {0.0, 0.0, 180.33}; // This is the receiver center
    double rec_width = 0.1; // Dummy value for flux map bin size
    double rec_radius = 15.45 / 2.0;
    double rec_height = 18.59;
    double rec_heat_shield_height = 3.2331;

    std::vector<double> scatter_aim_elevation;

    std::shared_ptr<SolTrace::Data::Sun> sun;
    std::vector<std::shared_ptr<Heliostat>> heliostat_field;
    std::shared_ptr<SolTrace::Data::SingleElement> top_heat_shield;
    std::shared_ptr<SolTrace::Data::SingleElement> receiver;
    std::shared_ptr<SolTrace::Data::SingleElement> bottom_heat_shield;

    SimulationData simData;
    RunnerT runner;

    // Sun outputs
    double sun_width;
    double sun_height;
    double A_sun_box;
    int nsun_rays;
    double power_per_ray;

    // Ray counts
    uint_fast64_t tot_helio_hits;
    uint_fast64_t tot_reflect_count;
    uint_fast64_t tot_helio_absorb_count;
    uint_fast64_t tot_helio_block_count;
    uint_fast64_t rec_absorb_count;
    uint_fast64_t heat_shield_absorb_count;
    uint_fast64_t miss_count;
    uint_fast64_t tot_rec_hits;
    uint_fast64_t rec_direct_count;
    uint_fast64_t rec_via_helio_count;
    uint_fast64_t rec_via_rec_count;        // These are rays that reflect inside the cylinder

    // Flux map
    HPM2D fluxGrid;
    std::vector<double> xValues, yValues;
    double binszx, binszy;
    double PeakFlux, PeakFluxUncertainty;
    double AveFlux, AveFluxUncertainty;
    double MinFlux, SigmaFlux, Uniformity;
    glm::dvec3 Centroid;
    double zScale;
    size_t NumberOfRays;
    int NotBinned;
    double max_neg_x_flux_err = 0;
    double max_pos_x_flux_err = 0;

    // Results
    double absorption_efficiency;
    double blocking_efficiency;
    double spillage_efficiency;
    double cosine_efficiency;
    double total_power;

    // Expected results
    double expected_power;
    double expected_peak_flux;
    double expected_flux_RMS;
    double expected_max_flux_RMS;

    double expected_cosine_efficiency;
    double expected_absorption_efficiency;
    double expected_blocking_efficiency;
    double expected_spillage_efficiency;

    HPM2D expected_fluxGrid;

    // Helper initialization (what used to be in SetUp)
    void initialize() {
        // Set parameters
        set_default_params();

        // Initialize runner
        RunnerStatus sts = runner.initialize();
        EXPECT_EQ(sts, RunnerStatus::SUCCESS);

        // Initial setup of receiver
        receiver = SolTrace::Data::make_element<SingleElement>();
        SolTrace::Data::OpticalPropertySet receiver_opt_set(SolTrace::Data::InteractionType::REFLECTION, "ReceiverOptics");
        receiver_opt_set.set_ideal_absorption(SolTrace::Data::OpticalSide::Front);
        receiver_opt_set.set_ideal_reflection(SolTrace::Data::OpticalSide::Back);
        auto receiver_ref = simData.add_optical_property_set(receiver_opt_set);
        receiver->set_optical_property_set(receiver_ref);
        receiver->set_aperture(SolTrace::Data::make_aperture<SolTrace::Data::Rectangle>(rec_radius * 2.0, rec_height));
        receiver->set_surface(SolTrace::Data::make_surface<SolTrace::Data::Cylinder>(rec_radius));
        glm::dvec3 offset = {0.0, 0.0, 0.0};
        glm::dvec3 rec_origin_offset = rec_origin + offset;
        glm::dvec3 v1 = {0.0, -1.0, 0.0};
        glm::dvec3 aim_point = rec_origin_offset + v1;
        receiver->set_reference_frame_geometry(rec_origin_offset, aim_point, 180.0);
        receiver->set_name("Receiver");
        receiver->enable();

        top_heat_shield = SolTrace::Data::make_element<SingleElement>();
        SolTrace::Data::OpticalPropertySet top_heat_shield_opt_set(SolTrace::Data::InteractionType::REFLECTION, "TopHeatShieldOptics");
        top_heat_shield_opt_set.set_ideal_absorption(SolTrace::Data::OpticalSide::Front);
        top_heat_shield_opt_set.set_ideal_reflection(SolTrace::Data::OpticalSide::Back);
        auto top_heat_shield_ref = simData.add_optical_property_set(top_heat_shield_opt_set);
        top_heat_shield->set_optical_property_set(top_heat_shield_ref);
        top_heat_shield->set_aperture(SolTrace::Data::make_aperture<SolTrace::Data::Rectangle>(rec_radius * 2.0, rec_heat_shield_height));
        top_heat_shield->set_surface(SolTrace::Data::make_surface<SolTrace::Data::Cylinder>(rec_radius));
        offset = { 0.0, 0.0, (rec_height + rec_heat_shield_height)/2.};
        rec_origin_offset = rec_origin + offset;
        aim_point = rec_origin_offset + v1;
        top_heat_shield->set_reference_frame_geometry(rec_origin_offset, aim_point, 0.0);
        top_heat_shield->set_name("Top Heat Shield");
        top_heat_shield->enable();

        bottom_heat_shield = SolTrace::Data::make_element<SingleElement>();
        SolTrace::Data::OpticalPropertySet bottom_heat_shield_opt_set(SolTrace::Data::InteractionType::REFLECTION, "BottomHeatShieldOptics");
        bottom_heat_shield_opt_set.set_ideal_absorption(SolTrace::Data::OpticalSide::Front);
        bottom_heat_shield_opt_set.set_ideal_reflection(SolTrace::Data::OpticalSide::Back);
        auto bottom_heat_shield_ref = simData.add_optical_property_set(bottom_heat_shield_opt_set);
        bottom_heat_shield->set_optical_property_set(bottom_heat_shield_ref);
        bottom_heat_shield->set_aperture(SolTrace::Data::make_aperture<SolTrace::Data::Rectangle>(rec_radius * 2.0, rec_heat_shield_height));
        bottom_heat_shield->set_surface(SolTrace::Data::make_surface<SolTrace::Data::Cylinder>(rec_radius));
        offset = { 0.0, 0.0, -(rec_height + rec_heat_shield_height) / 2. };
        rec_origin_offset = rec_origin + offset;
        aim_point = rec_origin_offset + v1;
        bottom_heat_shield->set_reference_frame_geometry(rec_origin_offset, aim_point, 0.0);
        bottom_heat_shield->set_name("Bottom Heat Shield");
        bottom_heat_shield->enable();
    }

    void set_high_accuracy_params() {
        SimulationParameters& params = simData.get_simulation_parameters();
        params.number_of_rays = 20.e6;
        params.max_number_of_rays = params.number_of_rays * 100;
        params.include_optical_errors = use_optical_errors;
        params.include_sun_shape_errors = use_sunshape_errors;
    }

    void set_default_params() {
        SimulationParameters& params = simData.get_simulation_parameters();
        params.number_of_rays = 5.e5;
        params.max_number_of_rays = params.number_of_rays * 100;
        params.include_optical_errors = use_optical_errors;
        params.include_sun_shape_errors = use_sunshape_errors;
        params.seed = seed;
    }

    void create_heliostat_field() {
        // Define mirror optical properties
        SolTrace::Data::OpticalPropertySet mirror_opt_set(SolTrace::Data::InteractionType::REFLECTION, "HeliostatMirrorOptics");
        mirror_opt_set.set_reflectivity(SolTrace::Data::OpticalSide::Front, 0.9);
        mirror_opt_set.set_ideal_absorption(SolTrace::Data::OpticalSide::Back);
        mirror_opt_set.set_errors(SolTrace::Data::OpticalSide::Front, SolTrace::Data::DistributionType::GAUSSIAN, 2.0, 0.0);
        mirror_opt_set.set_errors(SolTrace::Data::OpticalSide::Back, SolTrace::Data::DistributionType::NONE, 0.0, 0.0);
        auto mirror_ref = simData.add_optical_property_set(mirror_opt_set);

        // Reading in field layout and aimpoints
        std::vector<double> x_coords;
        std::vector<double> y_coords;
        scatter_aim_elevation.clear();
        std::string field_file_path = std::string(PROJECT_DIR) + std::string("/round_robin_study/field_positions_aimpoints.csv");
        std::ifstream file(field_file_path);
        if (!file.is_open()) {
            std::cout << "Could not open field_positions_aimpoints.csv" << std::endl;
            return;
        }
        std::string line;
        int skip_lines = 1;
        int lines_skipped = 0;
        while (std::getline(file, line)) {
            if (lines_skipped < skip_lines) {
                lines_skipped++;
                continue;
            }
            std::stringstream ss(line);
            std::string item;
            std::vector<std::string> tokens;
            while (std::getline(ss, item, ',')) {
                tokens.push_back(item);
            }
            if (tokens.size() >= 4) {
                x_coords.push_back(std::stod(tokens[1]));
                y_coords.push_back(std::stod(tokens[2]));
                scatter_aim_elevation.push_back(std::stod(tokens[3]));
            }
        }

        // Generate heliostat field
        heliostat_field.clear();
        for (size_t i = 0; i < x_coords.size(); i++) {
            glm::dvec3 heliostat_origin(x_coords[i], y_coords[i], 4.49);
            auto heliostat = SolTrace::Data::make_element<Heliostat>();
            heliostat->set_optics(mirror_ref);
            heliostat->set_reference_frame_geometry(heliostat_origin, khat, 0.0);
            heliostat->set_aperture_size(10.38, 9.73);   // Width, Height
            heliostat->set_number_panels(1, 1);
            heliostat->set_gaps(0, 0);
            heliostat->set_canting(Heliostat::NONE, 0.0, 0.0);
            // Set aim point on the receiver surface
            double distance = sqrt(pow(x_coords[i], 2) + pow(y_coords[i], 2));
            glm::dvec3 aim_point = {rec_radius * (x_coords[i] / distance),
                                    rec_radius * (y_coords[i] / distance),
                                    rec_origin.z};
            heliostat->set_target_position(aim_point);
            // At-slant focal length - to center of receiver
            glm::dvec3 slant = -rec_origin + heliostat_origin; //aim_point
            double focal_length = glm::length(slant);
            heliostat->set_focal_length(focal_length);

            heliostat->create_geometry();
            heliostat->set_name("Heliostat_" + std::to_string(i + 1));
            heliostat->enable();
            heliostat_field.push_back(heliostat);
        }
    }

    void set_scatter_aimpoints() {
        int helio_idx = 0;
        for (const auto& heliostat : heliostat_field) {
            glm::dvec3 heliostat_origin = heliostat->get_origin_global();
            double distance = sqrt(pow(heliostat_origin.x, 2) + pow(heliostat_origin.y, 2));
            glm::dvec3 aim_point = {rec_radius * (heliostat_origin.x / distance),
                                    rec_radius * (heliostat_origin.y / distance),
                                    rec_origin.z + scatter_aim_elevation[helio_idx]};
            heliostat->set_target_position(aim_point);
            heliostat->create_geometry();
            helio_idx++;
        }
    }

    void assign_focal_lengths_banded() {
        for (const auto& heliostat : heliostat_field) {
            glm::dvec3 heliostat_origin = heliostat->get_origin_global();
            double distance = sqrt(pow(heliostat_origin.x, 2) + pow(heliostat_origin.y, 2));
            double focal_length = 0.0;
            if (distance <= 502.5)
                focal_length = 353.8;
            else if (distance > 502.5 && distance <= 878.0)
                focal_length = 704.8;
            else if (distance > 878.0 && distance <= 1253.5)
                focal_length = 1072.5;
            else if (distance > 1253.5 && distance <= 1650.0)
                focal_length = 1444.3;
            else
                throw std::runtime_error("Heliostat distance out of range for banded focal lengths.");

            heliostat->set_focal_length(focal_length);
            heliostat->create_geometry();
        }
    }

    void assign_focal_lengths_canting_banded() {
        // NOTE: This was created for task 1b where the focal lengths were set to the canting distances.
        for (const auto& heliostat : heliostat_field) {
            glm::dvec3 heliostat_origin = heliostat->get_origin_global();
            double distance = sqrt(pow(heliostat_origin.x, 2) + pow(heliostat_origin.y, 2));
            double focal_length = 0.0;
            if (distance <= 502.0)
                focal_length = 516;
            else if (distance > 502.0 && distance <= 885.0)
                focal_length = 668.0;
            else if (distance > 885.0 && distance <= 1267.0)
                focal_length = 959.0;
            else if (distance > 1267.0 && distance <= 1650.0)
                focal_length = 1500.0;
            else
                throw std::runtime_error("Heliostat distance out of range for banded focal lengths.");

            heliostat->set_focal_length(focal_length);
            heliostat->create_geometry();
        }
    }

    void assign_canted_slant(bool flat_facets) {
        for (const auto& heliostat : heliostat_field) {
            glm::dvec3 slant = -rec_origin + heliostat->get_origin_global();
            double slant_distance = glm::length(slant);
            heliostat->set_number_panels(6, 5);
            heliostat->set_gaps(0.02, 0.02);
            heliostat->set_canting(Heliostat::CantingType::ON_AXIS, slant_distance, 0.0);

            if (flat_facets) heliostat->set_focal_length(0.0);   // Flat facets

            heliostat->create_geometry();
        }
    }

    void assign_canted_banded(bool flat_facets) {
        // Modify heliostat canting to bands - to center of receiver
        for (const auto& heliostat : heliostat_field) {
            glm::dvec3 heliostat_origin = heliostat->get_origin_global();
            double distance = sqrt(pow(heliostat_origin.x, 2) + pow(heliostat_origin.y, 2));
            double cant_distance = 0.0;
            if (distance <= 502.0)
                cant_distance = 516;
            else if (distance > 502.0 && distance <= 885.0)
                cant_distance = 668.0;
            else if (distance > 885.0 && distance <= 1267.0)
                cant_distance = 959.0;
            else if (distance > 1267.0 && distance <= 1650.0)
                cant_distance = 1500.0;
            else
                throw std::runtime_error("Heliostat distance out of range for banded canting lengths.");

            heliostat->set_number_panels(6, 5);
            heliostat->set_gaps(0.02, 0.02);
            heliostat->set_canting(Heliostat::CantingType::ON_AXIS, cant_distance, 0.0);
            if (flat_facets) heliostat->set_focal_length(0.0);   // Flat facets
            heliostat->create_geometry();
        }
    }

    void setup_simData() {
        // Set up sun
        glm::dvec3 sun_pos = {0.0, 0.0, 1000.0};
        sun = SolTrace::Data::make_ray_source<Sun>();
        sun->set_position(sun_pos);
        sun->set_shape(SolTrace::Data::SunShape::PILLBOX, 0.0, 4.65, 0.0);
        sun->set_gen_type(sun_gen_type);
        simData.add_ray_source(sun);

        // Set up stages
        stage_ptr st1 = SolTrace::Data::make_stage(1);
        st1->set_reference_frame_geometry(zero, khat, 0.0);

        for (const auto& heliostat : heliostat_field) {
            auto ret = st1->add_element(heliostat);
            EXPECT_TRUE(SolTrace::Data::Element::is_success(ret));
        }

        stage_ptr st2 = SolTrace::Data::make_stage(2);
        st2->set_reference_frame_geometry(zero, khat, 0.0);
        auto ret = st2->add_element(receiver);
        EXPECT_TRUE(SolTrace::Data::Element::is_success(ret));
        ret = st2->add_element(top_heat_shield);
        EXPECT_TRUE(SolTrace::Data::Element::is_success(ret));
        ret = st2->add_element(bottom_heat_shield);
        EXPECT_TRUE(SolTrace::Data::Element::is_success(ret));

        simData.add_stage(st1);
        simData.add_stage(st2);
    }

    void update_simulation_geometry(double azimuth, double elevation) {
        // Set sun position
        glm::dvec3 sun_pos;
        SolTrace::Data::sun_position_vector_degrees(sun_pos, azimuth, elevation);
        sun->set_position(sun_pos);
        // Update heliostat positions
        for (const auto& heliostat : heliostat_field) {
            heliostat->update_geometry(azimuth, elevation);
        }
    }

    void simulate(SimulationResult* result, int N_rays = -1) {
        if (high_accuracy) set_high_accuracy_params();
        else set_default_params();

        if (N_rays != -1)
        {
            auto& params = simData.get_simulation_parameters();
            params.number_of_rays = N_rays;
            params.max_number_of_rays = static_cast<std::uint_fast64_t>(N_rays) * 10000ULL;
        }

        RunnerStatus sts = runner.setup_simulation(&simData);
        EXPECT_EQ(sts, RunnerStatus::SUCCESS);
        sts = runner.run_simulation();
        EXPECT_EQ(sts, RunnerStatus::SUCCESS);
        sts = runner.report_simulation(result, 0);
        EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    }

    void calculate_sun_size(SimulationResult& result) {
        double dni = 1000.0; // W/m2 (constant for all tests)

        result.get_sun_dimensions(this->sun_width, this->sun_height);
        nsun_rays = result.get_sun_ray_count();
        A_sun_box = result.get_sun_A_box();
        power_per_ray = A_sun_box / nsun_rays * dni;

        if (print_info) {
            std::cout << "Power per ray: " << power_per_ray << std::endl;
            std::cout << "Sun box: " << sun_width << " x " << sun_height << std::endl;
            std::cout << "Sun ray count: " << nsun_rays << std::endl;
        }
    }

    void calculate_ray_counts(const SimulationResult& result) {
        tot_helio_hits = 0;
        tot_reflect_count = 0;
        tot_helio_absorb_count = 0;
        tot_helio_block_count = 0;

        rec_absorb_count = 0;
        heat_shield_absorb_count = 0;
        miss_count = 0;

        tot_rec_hits = 0;
        rec_direct_count = 0;
        rec_via_helio_count = 0;
        rec_via_rec_count = 0;

        for (size_t i = 0; i < result.get_number_of_records(); i++) {
            const ray_record_ptr rr = result[i];
            bool first_helio_hit = false;
            for (size_t j = 0; j < rr->interactions.size(); j++) {
                auto hit_element = rr->get_element(j);
                SolTrace::Result::RayEvent rev = rr->get_event(j);

                // Create or exit of rays
                if (hit_element < 0) {
                    if (rev == RayEvent::CREATE) first_helio_hit = true;
                    if (rev == RayEvent::EXIT) miss_count++;
                    continue;
                }

                // Check receiver element
                if (hit_element == receiver->get_id()) {
                    tot_rec_hits++;
                    if (rev == RayEvent::ABSORB) rec_absorb_count++;
                    if (j == 1) rec_direct_count++;
                    else if (j == 2)
                    {
                        rec_via_helio_count++;
                    }
                    else
                    {
                        rec_via_rec_count++;
                    }
                        
                    continue;
                }

                // Check heat shield elements
                if (hit_element == top_heat_shield->get_id() ||
                    hit_element == bottom_heat_shield->get_id()) {
                    if (rev == RayEvent::ABSORB) heat_shield_absorb_count++;
                    continue;
                }

                // Track hits of elements
                if (first_helio_hit) {
                    first_helio_hit = false;
                    tot_helio_hits++;
                    if (rev == RayEvent::REFLECT) tot_reflect_count++;
                    if (rev == RayEvent::ABSORB) tot_helio_absorb_count++;
                    continue;
                }
                else {
                    // Subsequent hits on a heliostat count as blocking
                    if (rev == RayEvent::ABSORB) tot_helio_block_count++;
                    continue;
                }
            }
        }

        if (print_info) {
            std::cout << "Heliostat Hit Count: " << tot_helio_hits << std::endl;
            std::cout << "Reflect Rays: " << tot_reflect_count << std::endl;
            std::cout << "Heliostat Absorbed Rays: " << tot_helio_absorb_count << std::endl;
            std::cout << "Heliostat Blocked Rays: " << tot_helio_block_count << std::endl;
            std::cout << "Receiver Absorbed Rays: " << rec_absorb_count << std::endl;
            std::cout << "Heat Shield Absorbed Rays: " << heat_shield_absorb_count << std::endl;
            std::cout << "Miss Rays: " << miss_count << std::endl;
        }
    }

    void calculate_outputs(const SimulationResult& result, bool ignore_direct = false) {
        
        absorption_efficiency = (double)tot_reflect_count / (double)tot_helio_hits;
        blocking_efficiency = 1.0 - (double)tot_helio_block_count / (double)tot_reflect_count;
        spillage_efficiency = (double)rec_absorb_count / (double)(tot_reflect_count - tot_helio_block_count);

        double field_area = 0.0;
        for (const auto& heliostat : heliostat_field) {
            field_area += heliostat->get_area();
        }
        cosine_efficiency = ((double)tot_helio_hits / (double)nsun_rays) * (A_sun_box / field_area);

        total_power = rec_absorb_count * power_per_ray / 1.e3;   // W to kW

        calculate_receiver_flux_map(result, 60, 23, true, ignore_direct);
    }

    void read_expected_summary_results(std::string filepath) {
        // Reading in fluxData summary
        std::ifstream file(filepath);
        if (!file.is_open()) {
            std::cout << "Could not open expected results file." << std::endl;
            std::cout << "Filepath:" << filepath << std::endl;
            return;
        }

        int line_number = 1;
        std::string line;
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string item;
            std::vector<std::string> tokens;
            while (std::getline(ss, item, ',')) {
                tokens.push_back(item);
            }
            // TODO: we could average values across the different models
            if (line_number == 2)
                expected_power = std::stod(tokens[1]);
            if (line_number == 3)
                expected_peak_flux = std::stod(tokens[1]);
            if (line_number == 6)
                expected_flux_RMS = std::stod(tokens[1]); // TODO: Maybe not required
            if (line_number == 8)
                expected_max_flux_RMS = std::stod(tokens[1]);
            line_number++;
        }
    }

    void read_expected_efficiency_results(std::string filepath) {
        // Reading in fluxData summary
        std::ifstream file(filepath);
        if (!file.is_open()) {
            std::cout << "Could not open expected results file." << std::endl;
            std::cout << "Filepath:" << filepath << std::endl;
            return;
        }

        int line_number = 1;
        std::string line;
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string item;
            std::vector<std::string> tokens;
            while (std::getline(ss, item, ',')) {
                tokens.push_back(item);
            }
            if (line_number > 2) {
                double avg_eff = (std::stod(tokens[2]) + std::stod(tokens[4]) + std::stod(tokens[6]))
                                 / 3.0;
                if (line_number == 3) expected_cosine_efficiency = avg_eff;
                if (line_number == 4) expected_absorption_efficiency = avg_eff;
                if (line_number == 5) expected_blocking_efficiency = avg_eff;
                if (line_number == 6) expected_spillage_efficiency = avg_eff;
            }
            line_number++;
        }
    }

    void read_expected_flux_map(std::string filepath) {
        // Reading in fluxData summary
        std::ifstream file(filepath);
        if (!file.is_open()) {
            std::cout << "Could not open expected flux map file." << std::endl;
            std::cout << "Filepath:" << filepath << std::endl;
            return;
        }

        expected_fluxGrid.clear();
        expected_fluxGrid.resize(23, 60);
        expected_fluxGrid.fill(0.0);

        int row_number = 0;
        std::string line;
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::string item;
            std::vector<std::string> tokens;
            while (std::getline(ss, item, ',')) {
                tokens.push_back(item);
            }

            for (size_t col_number = 0; col_number < tokens.size(); col_number++) {
                double val = std::stod(tokens[col_number]);
                expected_fluxGrid.at(row_number, col_number) = val;
            }
            row_number++;
        }
    }

    void read_expected_all_results(std::string task_number, std::string aim_strategy, std::string hour) {
        std::string round_robin_base = "/round_robin_study/results/phase_II/";
        std::string summary_file = "fluxDataSummary_Task_" + task_number + "_AimStrat_" + aim_strategy + "_Hour_" + hour + ".csv";
        read_expected_summary_results(std::string(PROJECT_DIR) + round_robin_base + summary_file);
        std::string efficiency_file = "efficiency_Task_" + task_number + "_AimStrat_" + aim_strategy + "_Hour_" + hour + ".csv";
        read_expected_efficiency_results(std::string(PROJECT_DIR) + round_robin_base + efficiency_file);
        std::string fluxmap_file = "soltrace_Task_" + task_number + "_AimStrat_" + aim_strategy + "_" + hour + "_fluxmap_further_aimpoint.csv";
        read_expected_flux_map(std::string(PROJECT_DIR) + round_robin_base + fluxmap_file);
    }

    void reset_flux_map() {
        fluxGrid.clear();
        xValues.clear();
        yValues.clear();
        binszx = binszy = 0.0;
        PeakFlux = PeakFluxUncertainty = 0.0;
        AveFlux = AveFluxUncertainty = 0.0;
        MinFlux = SigmaFlux = Uniformity = 0.0;
        Centroid = glm::dvec3{0.0};
        zScale = 0.0;
        NumberOfRays = 0;
    }

    bool calculate_receiver_flux_map(const SimulationResult& result, int nbinsx, int nbinsy, 
        bool is_cylinder, bool ignore_direct) 
    {
        reset_flux_map();

        double minx = -rec_width / 2.0;
        double maxx = rec_width / 2.0;

        double maxy = rec_height / 2.0;
        double miny = -rec_height / 2.0;

        glm::dvec3 rec_origin = receiver->get_origin_global();

        // Autoscale
        if (false) {
            minx = miny = 1e199;
            maxx = maxy = -1e199;
            glm::dvec3 local_position;
            glm::dvec3 global_position;

            // automatically size the min/max x and y
            for (size_t i = 0; i < result.get_number_of_records(); i++) {
                const ray_record_ptr rr = result[i];
                for (size_t j = 0; j < rr->interactions.size(); j++) {
                    if (receiver->get_id() == rr->get_element(j)) {
                        rr->get_position(j, global_position);
                        receiver->convert_global_to_local(local_position, global_position);

                        double x = local_position.x;
                        double y = local_position.y;

                        if (x < minx) minx = x;
                        if (x > maxx) maxx = x;
                        if (y < miny) miny = y;
                        if (y > maxy) maxy = y;
                    }
                }
            }
        }

        if (nbinsx <= 1 || nbinsy <= 1
            || maxx <= minx || maxy <= miny)
        {
            return false;
        }

        double gridszx = maxx - minx;
        double gridszy = maxy - miny;

        if (is_cylinder) {
            minx = -PI * rec_radius;
            maxx = PI * rec_radius;
            gridszx = maxx - minx;
        }

        binszx = gridszx / (double)nbinsx;
        binszy = gridszy / (double)nbinsy;

        xValues.resize(nbinsx);
        yValues.resize(nbinsy);
        fluxGrid.resize(nbinsy, nbinsx);

        // Bin rays here -> BinRaysXY()
        double x, y, z;
        int GridIncrementX, GridIncrementY;

        size_t RayCount = 0;
        NotBinned = 0;
        max_neg_x_flux_err = 0;
        max_pos_x_flux_err = 0;
        size_t npoints = 0;
        Centroid.x = Centroid.y = Centroid.z = 0.0;
        fluxGrid.fill(0.0);

        for (size_t i = 0; i < result.get_number_of_records(); i++) {
            glm::dvec3 local_position;
            glm::dvec3 global_position;

            const ray_record_ptr rr = result[i];
            for (size_t j = 0; j < rr->interactions.size(); j++) {
                if (receiver->get_id() == rr->get_element(j)) {
                    if (rr->get_event(j) == RayEvent::ABSORB) {
                        if (j > 1 || !ignore_direct)
                        {
                            rr->get_position(j, global_position);

                            receiver->convert_global_to_local(local_position, global_position);

                            x = local_position.x;
                            y = local_position.y;
                            z = local_position.z;

                            Centroid += local_position;

                            npoints++;

                            Centroid[0] += x;
                            Centroid[1] += y;
                            Centroid[2] += z;
                            npoints++;

                            if (is_cylinder) {
                                //if (z <= 0.0)
                                //    x = rec_radius * asin(x / rec_radius);
                                //else if (z > 0.0) {
                                //    if (x < 0) x = -(PI * rec_radius / 2.0 + rec_radius * acos(fabs(x) / rec_radius));
                                //    if (x >= 0) x = PI * rec_radius / 2.0 + rec_radius * acos(x / rec_radius);
                                //}

                                double rho = std::hypot(x, z);
                                if (rho <= 0.0)
                                {
                                    NotBinned++;
                                    continue;
                                }

                                // Project to cylinder surface
                                double xp = rec_radius * x / rho;
                                double zp = rec_radius * z / rho;

                                // Unwrap angle in [-pi, pi], with 0 at the front center (z < 0)
                                double theta = std::atan2(xp, -zp);

                                x = rec_radius * theta;
                            }

                            GridIncrementX = -1; //initialize grid increment counters
                            GridIncrementY = -1;

                            //determine which bin the ray falls into
                            while ((minx + (GridIncrementX + 1) * binszx) < x)
                                GridIncrementX++;

                            while ((miny + (GridIncrementY + 1) * binszy) < y)
                                GridIncrementY++;

                            GridIncrementX = nbinsx - GridIncrementX - 1;
                            GridIncrementY = nbinsy - GridIncrementY - 1;

                            if (GridIncrementX >= 0 && GridIncrementX < (int)fluxGrid.ncols()
                                && GridIncrementY >= 0 && GridIncrementY < (int)fluxGrid.nrows())
                            {
                                fluxGrid.at(GridIncrementY, GridIncrementX) += 1;//if ray falls inside a bin, increment count for that bin
                                RayCount++;  //increment ray intersection counter
                            }
                            else
                            {
                                NotBinned++;
                                //	qDebug("Not binned: [%d %d],  x=%lg, y=%lg", GridIncrementX, GridIncrementY, x, y);
                            }
                        }
                    }
                }
            }
        }

        //calculate midpoints of bins
        for (size_t i = 0; i < xValues.size(); i++)
            xValues[i] = minx + binszx / 2.0 + i * binszx;

        for (size_t i = 0; i < yValues.size(); i++)
            yValues[yValues.size() - i - 1] = miny + binszy / 2.0 + i * binszy;

        if (npoints > 0)
        {
            Centroid /= npoints;
        }

        double SumFlux, SumFlux2;
        PeakFlux = SumFlux = SumFlux2 = 0;
        MinFlux = 1e99;

        int NRaysInMinFluxBin = -1, NRaysInPeakFluxBin = -1;

        zScale = power_per_ray / (binszx * binszy);
        for (size_t r = 0; r < fluxGrid.nrows(); r++)
        {
            for (size_t c = 0; c < fluxGrid.ncols(); c++)
            {
                double z = fluxGrid.at(r, c) * zScale;
                SumFlux += z;
                SumFlux2 += z * z;

                if (z > PeakFlux) {
                    PeakFlux = z;
                    NRaysInPeakFluxBin = (int)fluxGrid.at(r, c);
                }

                if (z < MinFlux) {
                    MinFlux = z;
                    NRaysInMinFluxBin = (int)fluxGrid.at(r, c);
                }
            }
        }

        AveFlux = SumFlux / (nbinsx * nbinsy);
        SigmaFlux = sqrt((nbinsx * nbinsy * SumFlux2 - SumFlux * SumFlux) / (nbinsx * nbinsy * nbinsx * nbinsy));
        Uniformity = SigmaFlux / AveFlux;
        PeakFluxUncertainty = 100 / sqrt((double)NRaysInPeakFluxBin);
        AveFluxUncertainty = 100 / sqrt((double)result.get_number_of_records());
        // TODO: Should the be number of rays hitting the surface, not total rays traced?

        if (print_info) {
            std::cout << "Receiver auto bounds:" << std::endl;
            std::cout << "X min: " << minx << " X max: " << maxx << std::endl;
            std::cout << "Y min: " << miny << " Y max: " << maxy << std::endl;

            std::cout << "Peak flux: " << PeakFlux << std::endl;
            std::cout << "Peak flux uncertainty: +/- " << PeakFluxUncertainty << " %" << std::endl;
            std::cout << "Min flux: " << MinFlux << std::endl;
            std::cout << "Sigma flux: " << SigmaFlux << std::endl;
            std::cout << "Avg. flux: " << AveFlux << std::endl;
            std::cout << "Avg. flux uncertainty: +/- " << AveFluxUncertainty << " %" << std::endl;
            std::cout << "Uniformity: " << Uniformity << std::endl;
            std::cout << "Centroid: (" << Centroid.x << ", " << Centroid.y << ", " << Centroid.z
                      << ")" << std::endl;
        }

        return true;

    }

    void check_outputs(const SimulationResult& result) {
        SimulationParameters& params = simData.get_simulation_parameters();
        EXPECT_EQ(tot_helio_hits, params.number_of_rays);
        EXPECT_EQ(tot_helio_absorb_count + tot_reflect_count, tot_helio_hits);
        EXPECT_EQ(rec_absorb_count + heat_shield_absorb_count + miss_count + tot_helio_block_count, tot_reflect_count);

        double tol = high_accuracy ? 1.e-3 : 5.e-3;
        // Check efficiencies
        EXPECT_NEAR(absorption_efficiency, expected_absorption_efficiency, tol);
        EXPECT_NEAR(blocking_efficiency, expected_blocking_efficiency, tol * 2.0);
        EXPECT_NEAR(spillage_efficiency, expected_spillage_efficiency, tol);
        EXPECT_NEAR(cosine_efficiency, expected_cosine_efficiency, tol);

        // Total power absorbed
        EXPECT_NEAR(total_power, expected_power, tol * expected_power);

        // Peak flux value
        double peak_tol = high_accuracy ? 2.e-2 : 0.25;
        EXPECT_NEAR(PeakFlux / 1.e3, expected_peak_flux, peak_tol * expected_peak_flux);

        double rmse = 0.0;
        for (size_t r = 0; r < fluxGrid.nrows(); r++) {
            for (size_t c = 0; c < fluxGrid.ncols(); c++) {
                double sim_flux = fluxGrid.at(r, c) * zScale / 1.e3;
                double exp_flux = expected_fluxGrid.at(r, c);
                rmse += pow(sim_flux - exp_flux, 2);
            }
        }
        rmse = sqrt(rmse / (fluxGrid.nrows() * fluxGrid.ncols()));

        // RMS of flux values
        EXPECT_EQ(fluxGrid.nrows(), expected_fluxGrid.nrows());
        EXPECT_EQ(fluxGrid.ncols(), expected_fluxGrid.ncols());

        //EXPECT_LE(rmse, 25.0); // expected_flux_RMS
        double rmse_tol = high_accuracy ? 0.02 : 0.08;  // 2% of peak flux // 0.04
        EXPECT_LE(rmse / (PeakFlux / 1.e3), rmse_tol);

        if (print_info) {
            std::cout << "Cosine efficiency: " << cosine_efficiency * 100.0 << " %" << std::endl;
            std::cout << "Absorption efficiency: " << absorption_efficiency * 100.0 << " %" << std::endl;
            std::cout << "Blocking efficiency: " << blocking_efficiency * 100.0 << " %" << std::endl;
            std::cout << "Spillage efficiency: " << spillage_efficiency * 100.0 << " %" << std::endl;
            std::cout << "Total power: " << total_power << " kW" << std::endl;
            std::cout << "Flux map RMS error: " << rmse << " kW/m2" << std::endl;
            std::cout << "Peak flux: " << PeakFlux / 1.e3 << " kW/m2" << std::endl;
        }
    }

    void save_flux_map_to_file(std::string filename) {
        // Print out flux comparison to file
        std::ofstream outputFile(filename, std::ios::out | std::ios::trunc);
        if (!outputFile.is_open()) {
            std::cerr << "Error: Could not open the file " << filename << std::endl;
        }
        for (size_t r = 0; r < fluxGrid.nrows(); r++) {
            for (size_t c = 0; c < fluxGrid.ncols(); c++) {
                outputFile << fluxGrid.at(r, c) * zScale / 1.e3 << ",";
            }
            outputFile << std::endl;
        }
    }

    void save_flux_comparison_to_file(std::string filename) {
        // Print out flux comparison to file
        std::ofstream outputFile(filename, std::ios::out | std::ios::trunc);
        if (!outputFile.is_open()) {
            std::cerr << "Error: Could not open the file " << filename << std::endl;
        }
        for (size_t r = 0; r < fluxGrid.nrows(); r++) {
            for (size_t c = 0; c < fluxGrid.ncols(); c++) {
                double sim_flux = fluxGrid.at(r, c) * zScale / 1.e3;
                double exp_flux = expected_fluxGrid.at(r, c);
                outputFile << sim_flux << "," << exp_flux << "," << (sim_flux - exp_flux) << "," << (sim_flux - exp_flux) / exp_flux << std::endl;
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

    void update_from_hour(std::string hour)
    {
        // Update simulation geometry based on hour
        if (hour == "8")
            update_simulation_geometry(74.95, 26.26);  // Solar position at 8 AM
        else if (hour == "12")
            update_simulation_geometry(0.0, 61.97); // Solar Noon
        else
            throw std::invalid_argument("Hour not supported for simulate_check_outputs.");
    }

    void simulate_check_outputs(std::string task_number, std::string aim_strategy, std::string hour) {

        if (print_info) {
            std::cout << "\n\nTask: " << task_number << ", Aim: " << aim_strategy << ", Hour: " << hour << std::endl;
        }

        // Update simulation geometry based on hour
        update_from_hour(hour);

        SimulationResult result;
        simulate(&result);
        calculate_sun_size(result);
        calculate_ray_counts(result);
        read_expected_all_results(task_number, aim_strategy, hour);

        calculate_outputs(result);
        check_outputs(result);

        // Save results to file
        if (!save_results) return;
        std::string filename_postfix = "_Task_" + task_number + "_AimStrat_" + aim_strategy + "_Hour_" + hour;
        std::string flux_result_filename = "field_fluxmap" + filename_postfix + ".csv";
        save_flux_map_to_file(flux_result_filename);
        std::string flux_comparison_filename = "field_fluxmap_comparison" + filename_postfix + ".csv";
        save_flux_comparison_to_file(flux_comparison_filename);

        if (!save_raydata) return;
        std::string full_field_filename = "full_field_raydata" + filename_postfix + ".csv";
        result.write_csv_file(full_field_filename);
    }

};

// GTest fixture that reuses the helper logic
// This keeps existing TEST_F-based tests working

template <typename RunnerT>
class HeliostatFieldSimulation : public ::testing::Test, public HeliostatFieldSimulationHelper<RunnerT> {
protected:
    void SetUp() override {
        this->initialize();
    }

    void TearDown() override { }
};

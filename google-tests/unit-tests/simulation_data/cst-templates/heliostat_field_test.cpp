#include <iomanip>
#include <fstream>
#include <string>
#include <sstream>
#include <vector> 
#include <numeric> 

#include <gtest/gtest.h>

#include <native_runner.hpp>
#include <native_runner_types.hpp>
#include <simulation_data.hpp>
#include <simulation_result_export.hpp>
#include <stage_element.hpp>
#include <sun.hpp>

#include <../../hpvm.h>

#include <cst_templates/heliostat.hpp>
#include <cst_templates/utilities.hpp>

#include "common.hpp"
#include "count_absorbed_native.h"

using Heliostat = SolTrace::Data::Heliostat;

using SolTrace::Runner::RunnerStatus;
using SolTrace::NativeRunner::NativeRunner;
using SolTrace::NativeRunner::TRayData;
using SolTrace::NativeRunner::TSystem;
using SolTrace::NativeRunner::TSun;

class HeliostatFieldSimulation : public ::testing::Test {
public:
    bool high_accuracy = false;     // Runs 20 Million rays and tighter tolerance on checks
    bool print_info = false;        // Prints information on from simulation results (sun calculations, ray counts, flux calculations)
    bool save_results = false;      // Saves flux map results to CSV files
    bool save_raydata = false;      // Saves ray data to CSV file

    const Vector3d zero = { 0.0, 0.0, 0.0 }; // Global origin
    const Vector3d khat = { 0.0, 0.0, 1.0 }; // Global z-axis

    // Receiver parameters
    Vector3d rec_origin = { 0.0, 0.0, 180.33 };   // This is the receiver center
    double rec_radius = 15.45 / 2.0;
    double rec_height = 18.59;
    double rec_heat_shield_height = 3.2331;

    std::vector<double> scatter_aim_elevation;

    std::shared_ptr<SolTrace::Data::Sun> sun;
    std::vector<std::shared_ptr<Heliostat>> heliostat_field;
    std::shared_ptr<SolTrace::Data::SingleElement> top_heat_shield;
    std::shared_ptr<SolTrace::Data::SingleElement> receiver;
    std::shared_ptr<SolTrace::Data::SingleElement> bottom_heat_shield;

    // Tolerances
    // TODO: we could bring these into the class provide defaults and modify per test (if needed)

protected:
    SimulationData simData;
    NativeRunner runner;

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

    // Flux map
    HPM2D fluxGrid;
    std::vector<double> xValues, yValues;
    double binszx, binszy;
    double PeakFlux, PeakFluxUncertainty;
    double AveFlux, AveFluxUncertainty;
    double MinFlux, SigmaFlux, Uniformity;
    double Centroid[3];
    double zScale;
    size_t NumberOfRays;

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

    void SetUp() override {
        // Set parameters
        SimulationParameters& params = simData.get_simulation_parameters();
        params.number_of_rays = 5.e5;
        params.max_number_of_rays = params.number_of_rays * 100;
        params.include_optical_errors = true;
        params.include_sun_shape_errors = true;
        params.seed = 123;

        // Initialize runner
        RunnerStatus sts = runner.initialize();
        EXPECT_EQ(sts, RunnerStatus::SUCCESS);
        runner.enable_power_tower();
        runner.enable_point_focus();
        runner.set_number_of_threads(14);

        // Initial setup of receiver
        receiver = SolTrace::Data::make_element<SingleElement>();
        receiver->get_front_optical_properties()->set_ideal_absorption();
        receiver->get_back_optical_properties()->set_ideal_reflection();
        receiver->set_aperture(SolTrace::Data::make_aperture<SolTrace::Data::Rectangle>(rec_radius * 2.0, rec_height));
        receiver->set_surface(SolTrace::Data::make_surface<SolTrace::Data::Cylinder>(rec_radius));
        Vector3d offset = { 0.0, rec_radius, 0.0 };    // Cylinder origin is on the edge
        Vector3d rec_origin_offset;
        vector_add(1.0, rec_origin, 1.0, offset, rec_origin_offset);
        Vector3d v1 = { 0.0, -1.0, 0.0 };
        Vector3d aim_point;
        vector_add(1.0, rec_origin_offset, 1.0, v1, aim_point);
        receiver->set_reference_frame_geometry(rec_origin_offset, aim_point, 180.0);
        receiver->set_name("Receiver");
        receiver->enable();

        top_heat_shield = SolTrace::Data::make_element<SingleElement>();
        top_heat_shield->get_front_optical_properties()->set_ideal_absorption();
        top_heat_shield->get_back_optical_properties()->set_ideal_reflection();
        top_heat_shield->set_aperture(SolTrace::Data::make_aperture<SolTrace::Data::Rectangle>(rec_radius * 2.0, rec_heat_shield_height));
        top_heat_shield->set_surface(SolTrace::Data::make_surface<SolTrace::Data::Cylinder>(rec_radius));
        offset = { 0.0, rec_radius, (rec_height + rec_heat_shield_height)/2.};    // Cylinder origin is on the edge
        vector_add(1.0, rec_origin, 1.0, offset, rec_origin_offset);
        vector_add(1.0, rec_origin_offset, 1.0, v1, aim_point);
        top_heat_shield->set_reference_frame_geometry(rec_origin_offset, aim_point, 0.0);
        top_heat_shield->set_name("Top Heat Shield");
        top_heat_shield->enable();

        bottom_heat_shield = SolTrace::Data::make_element<SingleElement>();
        bottom_heat_shield->get_front_optical_properties()->set_ideal_absorption();
        bottom_heat_shield->get_back_optical_properties()->set_ideal_reflection();
        bottom_heat_shield->set_aperture(SolTrace::Data::make_aperture<SolTrace::Data::Rectangle>(rec_radius * 2.0, rec_heat_shield_height));
        bottom_heat_shield->set_surface(SolTrace::Data::make_surface<SolTrace::Data::Cylinder>(rec_radius));
        offset = { 0.0, rec_radius, -(rec_height + rec_heat_shield_height) / 2. };    // Cylinder origin is on the edge
        vector_add(1.0, rec_origin, 1.0, offset, rec_origin_offset);
        vector_add(1.0, rec_origin_offset, 1.0, v1, aim_point);
        bottom_heat_shield->set_reference_frame_geometry(rec_origin_offset, aim_point, 0.0);
        bottom_heat_shield->set_name("Bottom Heat Shield");
        bottom_heat_shield->enable();
    }

    void set_high_accuracy_params() {
        SimulationParameters& params = simData.get_simulation_parameters();
        params.number_of_rays = 20.e6;
        params.max_number_of_rays = params.number_of_rays * 100;
    }

    void create_heliostat_field() {
        // Define mirror optical properties
        OpticalProperties mirror;
        mirror.set_ideal_reflection();
        mirror.reflectivity = 0.9;
        mirror.slope_error = 2.0;       // mrad
        mirror.specularity_error = 0.0; // mrad
        mirror.error_distribution_type = DistributionType::GAUSSIAN;
        // Backside is absorbing
        OpticalProperties mirror_back;
        mirror_back.set_ideal_absorption();

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
            Vector3d heliostat_origin(x_coords[i], y_coords[i], 4.49);
            auto heliostat = SolTrace::Data::make_element<Heliostat>();
            heliostat->set_optics(mirror, mirror_back);
            heliostat->set_reference_frame_geometry(heliostat_origin, khat, 0.0);
            heliostat->set_aperture_size(10.38, 9.73);   // Width, Height
            heliostat->set_number_panels(1, 1);
            heliostat->set_gaps(0, 0);
            heliostat->set_canting(Heliostat::NONE, 0.0, 0.0);
            // Set aim point on the receiver surface
            double distance = sqrt(pow(x_coords[i], 2) + pow(y_coords[i], 2));
            Vector3d aim_point = { rec_radius * (x_coords[i] / distance),
                                   rec_radius * (y_coords[i] / distance),
                                   rec_origin[2]};
            heliostat->set_target_position(aim_point);
            // At-slant focal length - to center of receiver
            Vector3d slant;
            vector_add(-1.0, rec_origin, 1.0, heliostat_origin, slant);  //aim_point
            double focal_length = vector_norm(slant);
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
            Vector3d heliostat_origin = heliostat->get_origin_global();
            double distance = sqrt(pow(heliostat_origin[0], 2) + pow(heliostat_origin[1], 2));
            Vector3d aim_point = { rec_radius * (heliostat_origin[0] / distance),
                                   rec_radius * (heliostat_origin[1] / distance),
                                   rec_origin[2] + scatter_aim_elevation[helio_idx] };
            heliostat->set_target_position(aim_point);
            heliostat->create_geometry();
            helio_idx++;
        }
    }

    void assign_focal_lengths_banded() {
        for (const auto& heliostat : heliostat_field) {
            Vector3d heliostat_origin = heliostat->get_origin_global();
            double distance = sqrt(pow(heliostat_origin[0], 2) + pow(heliostat_origin[1], 2));
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
            Vector3d heliostat_origin = heliostat->get_origin_global();
            double distance = sqrt(pow(heliostat_origin[0], 2) + pow(heliostat_origin[1], 2));
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
            Vector3d slant;
            vector_add(-1.0, rec_origin, 1.0, heliostat->get_origin_global(), slant);
            double slant_distance = vector_norm(slant);
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
            Vector3d heliostat_origin = heliostat->get_origin_global();
            double distance = sqrt(pow(heliostat_origin[0], 2) + pow(heliostat_origin[1], 2));
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
        Vector3d sun_pos = { 0.0, 0.0, 1000.0 };
        sun = SolTrace::Data::make_ray_source<Sun>();
        sun->set_position(sun_pos);
        sun->set_shape(SolTrace::Data::SunShape::PILLBOX, 0.0, 4.65, 0.0);
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
        Vector3d sun_pos;
        sun_position_vector_degrees(sun_pos, azimuth, elevation);
        sun->set_position(sun_pos);
        // Update heliostat positions
        for (const auto& heliostat : heliostat_field) {
            heliostat->update_geometry(azimuth, elevation);
        }
    }

    void simulate(SimulationResult* result) {
        if (high_accuracy) set_high_accuracy_params();
        else { // Default parameters
            SimulationParameters& params = simData.get_simulation_parameters();
            params.number_of_rays = 5.e5;
            params.max_number_of_rays = params.number_of_rays * 100;
        }

        RunnerStatus sts = runner.setup_simulation(&simData);
        EXPECT_EQ(sts, RunnerStatus::SUCCESS);
        sts = runner.run_simulation();
        EXPECT_EQ(sts, RunnerStatus::SUCCESS);
        sts = runner.report_simulation(result, 0);
        EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    }

    void calculate_sun_size() {
        double dni = 1000.0; // W/m2 (constant for all tests)
        const TSystem* sys = runner.get_system();
        const TSun* sun = &(sys->Sun);
        sun_width = (sun->MaxXSun - sun->MinXSun);
        sun_height = (sun->MaxYSun - sun->MinYSun);
        A_sun_box = sun_width * sun_height;
        nsun_rays = sys->SunRayCount;
        power_per_ray = A_sun_box / nsun_rays * dni;

        if (print_info) {
            std::cout << "Power per ray: " << power_per_ray << std::endl;
            std::cout << "Sun box: " << sun_width << " x " << sun_height << std::endl;
            std::cout << "Sun ray count: " << nsun_rays << std::endl;
        }
    }

    void calculate_ray_counts(const SimulationResult &result) {
        tot_helio_hits = 0;
        tot_reflect_count = 0;
        tot_helio_absorb_count = 0;
        tot_helio_block_count = 0;

        rec_absorb_count = 0;
        heat_shield_absorb_count = 0;
        miss_count = 0;

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
                    if (rev == RayEvent::ABSORB) rec_absorb_count++;
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
            if (line_number == 2) expected_power = std::stod(tokens[1]);
            if (line_number == 3) expected_peak_flux = std::stod(tokens[1]);
            if (line_number == 6) expected_flux_RMS = std::stod(tokens[1]);     // TODO: Maybe not required
            if (line_number == 8) expected_max_flux_RMS = std::stod(tokens[1]);
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
                double avg_eff = (std::stod(tokens[2]) + std::stod(tokens[4]) + std::stod(tokens[6])) / 3.0;
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
        Centroid[0] = Centroid[1] = Centroid[2] = 0.0;
        zScale = 0.0;
        NumberOfRays = 0;
    }

    bool calculate_receiver_flux_map(const SimulationResult &result, int nbinsx, int nbinsy, bool is_cylinder) {
        reset_flux_map();

        double minx, maxx, miny, maxy;
        minx = maxx = miny = maxy = 0.0;
        Vector3d rec_origin = receiver->get_origin_global();

        // Autoscale
        if (true) {
            minx = miny = 1e199;
            maxx = maxy = -1e199;
            Vector3d local_position;
            Vector3d global_position;

            // automatically size the min/max x and y
            for (size_t i = 0; i < result.get_number_of_records(); i++) {
                const ray_record_ptr rr = result[i];
                for (size_t j = 0; j < rr->interactions.size(); j++) {
                    if (receiver->get_id() == rr->get_element(j)) {
                        rr->get_position(j, global_position);
                        receiver->convert_global_to_local(local_position, global_position);

                        double x = local_position[0];
                        double y = local_position[1];

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
            gridszx = 2.0 * PI * rec_radius;
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
        size_t NotBinned = 0;
        size_t npoints = 0;
        Centroid[0] = Centroid[1] = Centroid[2] = 0.0;
        fluxGrid.fill(0.0);

        for (size_t i = 0; i < result.get_number_of_records(); i++) {
            Vector3d local_position;
            Vector3d global_position;

            const ray_record_ptr rr = result[i];
            for (size_t j = 0; j < rr->interactions.size(); j++) {
                if (receiver->get_id() == rr->get_element(j)) {
                    if (rr->get_event(j) == RayEvent::ABSORB) {
                        rr->get_position(j, global_position);

                        receiver->convert_global_to_local(local_position, global_position);

                        x = local_position[0];
                        y = local_position[1];
                        z = local_position[2];

                        Centroid[0] += x;
                        Centroid[1] += y;
                        Centroid[2] += z;
                        npoints++;

                        if (is_cylinder) {
                            if (z <= rec_radius)
                                x = rec_radius * asin(x / rec_radius);
                            else if (z > rec_radius) {
                                if (x < 0) x = -(PI * rec_radius / 2.0 + rec_radius * acos(fabs(x) / rec_radius));
                                if (x >= 0) x = PI * rec_radius / 2.0 + rec_radius * acos(x / rec_radius);
                            }
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

        //calculate midpoints of bins
        for (size_t i = 0; i < xValues.size(); i++)
            xValues[i] = minx + binszx / 2.0 + i * binszx;

        for (size_t i = 0; i < yValues.size(); i++)
            yValues[yValues.size() - i - 1] = miny + binszy / 2.0 + i * binszy;

        if (npoints > 0)
        {
            Centroid[0] /= npoints;
            Centroid[1] /= npoints;
            Centroid[2] /= npoints;
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
            std::cout << "Centroid: (" << Centroid[0] << ", " << Centroid[1] << ", " << Centroid[2] << ")" << std::endl;
        }

        return true;

    }

    void check_outputs(const SimulationResult &result) {
        SimulationParameters& params = simData.get_simulation_parameters();
        EXPECT_EQ(tot_helio_hits, params.number_of_rays);
        EXPECT_EQ(tot_helio_absorb_count + tot_reflect_count, tot_helio_hits);
        EXPECT_EQ(rec_absorb_count + heat_shield_absorb_count + miss_count + tot_helio_block_count, tot_reflect_count);

        double tol = high_accuracy ? 1.e-3 : 5.e-3;
        // Check efficiencies
        double absorption_efficiency = (double)tot_reflect_count / (double)tot_helio_hits;
        EXPECT_NEAR(absorption_efficiency, expected_absorption_efficiency, tol);
        double blocking_efficiency = 1.0 - (double)tot_helio_block_count / (double)tot_reflect_count;
        EXPECT_NEAR(blocking_efficiency, expected_blocking_efficiency, tol * 2.0);  
        double spillage_efficiency = (double) rec_absorb_count / (double)(tot_reflect_count - tot_helio_block_count);
        EXPECT_NEAR(spillage_efficiency, expected_spillage_efficiency, tol);

        double field_area = 0.0;
        for (const auto& heliostat : heliostat_field) {
            field_area += heliostat->get_area();
        }
        double cosine_efficiency = ((double)tot_helio_hits / (double)nsun_rays) * (A_sun_box / field_area);
        EXPECT_NEAR(cosine_efficiency, expected_cosine_efficiency, tol);

        // Total power absorbed
        double total_power = rec_absorb_count * power_per_ray / 1.e3;   // W to kW
        EXPECT_NEAR(total_power, expected_power, tol * expected_power);

        // Peak flux value
        double peak_tol = high_accuracy ? 2.e-2 : 0.25;
        calculate_receiver_flux_map(result, 60, 23, true);
        EXPECT_NEAR(PeakFlux / 1.e3, expected_peak_flux, peak_tol * expected_peak_flux);

        // RMS of flux values
        EXPECT_EQ(fluxGrid.nrows(), expected_fluxGrid.nrows());
        EXPECT_EQ(fluxGrid.ncols(), expected_fluxGrid.ncols());
        double rmse = 0.0;
        for (size_t r = 0; r < fluxGrid.nrows(); r++) {
            for (size_t c = 0; c < fluxGrid.ncols(); c++) {
                double sim_flux = fluxGrid.at(r, c) * zScale / 1.e3;
                double exp_flux = expected_fluxGrid.at(r, c);
                rmse += pow(sim_flux - exp_flux, 2);
            }
        }
        rmse = sqrt(rmse / (fluxGrid.nrows() * fluxGrid.ncols()));

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

    void simulate_check_outputs(std::string task_number, std::string aim_strategy, std::string hour) {
        
        if (print_info) {
            std::cout << "\n\nTask: " << task_number << ", Aim: " << aim_strategy << ", Hour: " << hour << std::endl;
        }

        // Update simulation geometry based on hour
        if (hour == "8") 
            update_simulation_geometry(74.95, 26.26);  // Solar position at 8 AM
        else if (hour == "12") 
            update_simulation_geometry(0.0, 61.97); // Solar Noon
        else
            throw std::invalid_argument("Hour not supported for simulate_check_outputs.");        
        
        SimulationResult result;
        simulate(&result);
        calculate_sun_size();
        calculate_ray_counts(result);
        read_expected_all_results(task_number, aim_strategy, hour);
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

    void TearDown() override {
        // Code here will be called immediately after each test (right
        // before the destructor).
    }

};


TEST_F(HeliostatFieldSimulation, singleFacet_SlantFocused)
{
    // Centerline aimpoints
    create_heliostat_field();
    setup_simData();
    simulate_check_outputs("1a", "1", "8");
    simulate_check_outputs("1a", "1", "12");

    // Scatter aimpoints
    set_scatter_aimpoints();
    simulate_check_outputs("1a", "2", "8");
    simulate_check_outputs("1a", "2", "12");
}

TEST_F(HeliostatFieldSimulation, singleFacet_BandFocused)
{
    // Centerline aimpoints
    create_heliostat_field();
    assign_focal_lengths_canting_banded();
    setup_simData();
    simulate_check_outputs("1b", "1", "8");
    simulate_check_outputs("1b", "1", "12");

    // Scatter aimpoints
    set_scatter_aimpoints();
    simulate_check_outputs("1b", "2", "8");
    simulate_check_outputs("1b", "2", "12");
}

TEST_F(HeliostatFieldSimulation, multiFacet_SlantCanted)
{
    // Centerline aimpoints
    create_heliostat_field();
    assign_canted_slant(true);      // Flat facets

    setup_simData();
    simulate_check_outputs("2a", "1", "8");
    simulate_check_outputs("2a", "1", "12");

    // Scatter aimpoints
    set_scatter_aimpoints();
    simulate_check_outputs("2a", "2", "8");
    simulate_check_outputs("2a", "2", "12");
}

TEST_F(HeliostatFieldSimulation, multiFacet_BandCanted)
{
    // Centerline aimpoints
    create_heliostat_field();
    assign_canted_banded(true);     // Flat facets

    setup_simData();
    simulate_check_outputs("2b", "1", "8");
    simulate_check_outputs("2b", "1", "12");

    // Scatter aimpoints
    set_scatter_aimpoints();
    simulate_check_outputs("2b", "2", "8");
    simulate_check_outputs("2b", "2", "12");
}

TEST_F(HeliostatFieldSimulation, multiFacet_SlantFocused_SlantCanted)
{
    // Center aimpoints;
    create_heliostat_field();
    assign_canted_slant(false);     // Slant focused (default)

    setup_simData();
    simulate_check_outputs("3a", "1", "8");
    simulate_check_outputs("3a", "1", "12");

    // Scatter aimpoints
    set_scatter_aimpoints();
    simulate_check_outputs("3a", "2", "8");
    simulate_check_outputs("3a", "2", "12");
}

TEST_F(HeliostatFieldSimulation, multiFacet_BandFocused_BandCanted)
{
    // Center aimpoints
    create_heliostat_field();
    assign_canted_banded(false);    // Canted by band
    assign_focal_lengths_banded();  // Focused by band

    setup_simData();
    simulate_check_outputs("3b", "1", "8");
    simulate_check_outputs("3b", "1", "12");

    // Scatter aimpoints
    set_scatter_aimpoints();
    simulate_check_outputs("3b", "2", "8");
    simulate_check_outputs("3b", "2", "12");
}
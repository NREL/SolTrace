#pragma once

#include <gtest/gtest.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>

#include <native_runner.hpp>
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

// Helper class for single heliostat simulation logic reusable inside tests
// RunnerT must implement the same interface as NativeRunner/EmbreeRunner/OptixRunner
template <typename RunnerT>
class SingleHeliostatSimulationHelper {
public:

    bool high_accuracy = false;
    bool print_info = false;
    bool save_results = false;
    bool save_raydata = false;

    const glm::dvec3 zero = {0.0, 0.0, 0.0}; // Global origin
    const glm::dvec3 khat = {0.0, 0.0, 1.0}; // Global z-axis

    double solar_azimuth = 180.0;
    double solar_elevation = 59.96377;

    double rec_radius = 0.0;
    double rec_height = 18.0;
    double rec_width = 12.0;

    std::shared_ptr<SolTrace::Data::Sun> sun;
    std::shared_ptr<Heliostat> heliostat;
    std::shared_ptr<SolTrace::Data::SingleElement> receiver;

    bool use_optical_errors = true;
    bool use_sunshape_errors = true;
    DistributionType error_dist = DistributionType::GAUSSIAN;
    double spec_error = 0.0;
    double slope_error = 2.0;
    int seed = 123;
    SolTrace::Data::GenType sun_gen_type = SolTrace::Data::GenType::RANDOM;
    SunShape sun_shape = SolTrace::Data::SunShape::PILLBOX;
    double half_width = 4.65;
    double gauss_sigma = 0;
    double csr = 0;
    std::vector<double> user_angle = {};
    std::vector<double> user_intensity = {};

    SimulationData simData;
    RunnerT runner;
    SimulationResult result;

    double sun_width;
    double sun_height;
    double A_sun_box;
    double power_per_ray;

    uint_fast64_t sun_ray_count = 0;

    uint_fast64_t helio_hit_count;
    uint_fast64_t reflect_count;
    uint_fast64_t helio_absorb_count;
    uint_fast64_t rec_absorb_count;
    uint_fast64_t miss_count;
    uint_fast64_t rec_hit_count = 0;
    uint_fast64_t rec_direct_hit_count = 0;
    uint_fast64_t rec_via_helio_hit_count = 0;

    HPM2D fluxGrid;
    std::vector<double> xValues, yValues;
    double binszx, binszy;
    double PeakFlux, PeakFluxUncertainty;
    double AveFlux, AveFluxUncertainty;
    double MinFlux, SigmaFlux, Uniformity;
    double Centroid[3];
    double zScale;
    size_t NumberOfRays;

    double expected_power;
    double expected_peak_flux;
    double expected_flux_RMS;

    HPM2D expected_fluxGrid;

    // Helper initialization (what used to be in SetUp)
    void initialize() {
        // Set parameters
        set_default_params();

        // Initialize runner
        RunnerStatus sts = runner.initialize();
        EXPECT_EQ(sts, RunnerStatus::SUCCESS);
        // Native runner speedups
        //runner.disable_power_tower();
        //runner.disable_point_focus();

        // Define mirror optical properties
        OpticalProperties mirror;
        mirror.set_ideal_reflection();
        mirror.reflectivity = 0.9;
        mirror.slope_error = this->slope_error;       // mrad
        mirror.specularity_error = this->spec_error; // mrad
        mirror.error_distribution_type = this->error_dist;
        // Backside is absorbing
        OpticalProperties mirror_back;
        mirror_back.set_ideal_absorption();

        // Initial setup of heliostat
        glm::dvec3 heliostat_origin(0.0, 500.0, 5.65);
        glm::dvec3 rec_origin(0.0, 0.0, 169.0);
        heliostat = SolTrace::Data::make_element<Heliostat>();
        heliostat->set_optics(mirror, mirror_back);
        heliostat->set_reference_frame_geometry(heliostat_origin, khat, 0.0);
        heliostat->set_aperture_size(11.415, 10.42);   // Width, Height
        heliostat->set_number_panels(1, 1);
        heliostat->set_gaps(0, 0);
        heliostat->set_focal_length(0.0);
        heliostat->set_canting(Heliostat::NONE, 0.0, 0.0);
        heliostat->set_target_position(rec_origin);
        heliostat->create_geometry();
        heliostat->set_name("Heliostat");
        heliostat->enable();

        // Initial setup of receiver
        receiver = SolTrace::Data::make_element<SingleElement>();
        receiver->get_front_optical_properties()->set_ideal_absorption();
        receiver->get_back_optical_properties()->set_ideal_reflection();
        receiver->set_aperture(SolTrace::Data::make_aperture<SolTrace::Data::Rectangle>(rec_width, rec_height));
        receiver->set_surface(SolTrace::Data::make_surface<SolTrace::Data::Flat>());
        glm::dvec3 v1 = {0.0, 1.0, 0.0}; // Pointing North TODO: change to point towards heliostat
        glm::dvec3 aim_point = rec_origin + v1;
        receiver->set_reference_frame_geometry(rec_origin, aim_point, 0.0);
        receiver->set_name("Receiver");
        receiver->enable();
    }

    void setup_simData() {
        // Set up sun
        glm::dvec3 sun_pos = {0.0, 0.0, 1000.0};
        sun = SolTrace::Data::make_ray_source<Sun>();
        sun->set_position(sun_pos);
        sun->set_shape(sun_shape, gauss_sigma, half_width, csr, user_angle, user_intensity);
        sun->set_gen_type(sun_gen_type);
        simData.add_ray_source(sun);

        // Set up stages
        stage_ptr st1 = SolTrace::Data::make_stage(1);
        st1->set_reference_frame_geometry(zero, khat, 0.0);
        auto ret = st1->add_element(heliostat);
        EXPECT_TRUE(SolTrace::Data::Element::is_success(ret));

        stage_ptr st2 = SolTrace::Data::make_stage(2);
        st2->set_reference_frame_geometry(zero, khat, 0.0);
        ret = st2->add_element(receiver);
        EXPECT_TRUE(SolTrace::Data::Element::is_success(ret));

        simData.add_stage(st1);
        simData.add_stage(st2);
    }

    void set_high_accuracy_params() {
        SimulationParameters& params = simData.get_simulation_parameters();
        params.number_of_rays = 20.e6;
        params.max_number_of_rays = params.number_of_rays * 100;
        params.include_optical_errors = this->use_optical_errors;
        params.include_sun_shape_errors = this->use_sunshape_errors;
    }

    void set_default_params() {
        SimulationParameters& params = simData.get_simulation_parameters();
        params.number_of_rays = 5.e5;
        params.max_number_of_rays = params.number_of_rays * 100;
        params.include_optical_errors = this->use_optical_errors;
        params.include_sun_shape_errors = this->use_sunshape_errors;
        params.seed = this->seed;
    }

    void set_heliostat_to_southeast() {
        // Update heliostat position to southeast of tower
        glm::dvec3 helio_origin(200.0, -200.0, 5.65); // Southeast of tower
        heliostat->set_reference_frame_geometry(helio_origin, khat, 0.0);

        // Point receiver to heliostat without tilting down
        helio_origin[2] = 0.0; // Project to ground plane
        glm::dvec3 rec_origin = receiver->get_origin_ref();
        glm::dvec3 aim_point = rec_origin + helio_origin;
        receiver->set_reference_frame_geometry(rec_origin, aim_point, 90.0);
    }

    void set_slant_focal_length() {
        // Set focal length to slant range
        glm::dvec3 distance = -receiver->get_origin_global() + heliostat->get_origin_global();
        double focal_length = glm::length(distance);
        heliostat->set_focal_length(focal_length);
        heliostat->create_geometry();
    }

    void set_flat_multi_facet() {
        heliostat->set_number_panels(7, 5);
        heliostat->set_gaps(0.03, 0.03);
        heliostat->create_geometry();
    }

    void set_onaxis_slant_canting() {
        heliostat->set_number_panels(7, 5);
        heliostat->set_gaps(0.03, 0.03);
        // Set on-axis canting to slant range
        glm::dvec3 distance = -receiver->get_origin_global() + heliostat->get_origin_global();
        double focal_length = glm::length(distance);
        heliostat->set_canting(Heliostat::CantingType::ON_AXIS, focal_length, 0.0);
        heliostat->create_geometry();
    }

    void update_simulation_geometry(double azimuth, double elevation) {
        // Set sun position
        glm::dvec3 sun_pos;
        SolTrace::Data::sun_position_vector_degrees(sun_pos, azimuth, elevation);
        sun->set_position(sun_pos);
        // Update heliostat position
        heliostat->update_geometry(azimuth, elevation);
    }

    void simulate(SimulationResult* result, int N_rays = -1) {
        if (high_accuracy) set_high_accuracy_params();
        else set_default_params();

        if (N_rays != -1)
        {
            SimulationParameters& params = simData.get_simulation_parameters();
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
        sun_ray_count = result.get_sun_ray_count();
        A_sun_box = result.get_sun_A_box();
        power_per_ray = A_sun_box / sun_ray_count * dni;

        if (print_info) {
            std::cout << "Power per ray: " << power_per_ray << std::endl;
            std::cout << "Sun box: " << sun_width << " x " << sun_height << std::endl;
            std::cout << "Sun ray count: " << this->sun_ray_count << std::endl;
        }
    }

    void calculate_ray_counts(SimulationResult result) {
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
                for (auto iter = heliostat->get_const_iterator(); !heliostat->is_at_end(iter); ++iter) {
                    element_id facet_id = iter->second->get_id();
                    if (hit_element == facet_id) {
                        helio_hit_count++;
                        if (rev == RayEvent::REFLECT) reflect_count++;
                        if (rev == RayEvent::ABSORB) helio_absorb_count++;
                    }
                }

                // Check receiver element
                if (hit_element == receiver->get_id()) {
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

        if (print_info) {
            std::cout << "Heliostat Hit Count: " << helio_hit_count << std::endl;
            std::cout << "Reflect Rays: " << reflect_count << std::endl;
            std::cout << "Heliostat Absorbed Rays: " << helio_absorb_count << std::endl;
            std::cout << "Receiver Absorbed Rays: " << rec_absorb_count << std::endl;
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
            // calculate the average value for lines with multiple entries
            if (line_number == 2 || line_number == 3 || line_number == 6) {
                double value = 0.0;
                for (size_t i = 1; i < tokens.size(); i++) {
                    value += std::stod(tokens[i]);
                }
                value /= (double)(tokens.size() - 1);

                if (line_number == 2) expected_power = value;
                if (line_number == 3) expected_peak_flux = value;
                if (line_number == 6) expected_flux_RMS = value;
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
        expected_fluxGrid.resize(150, 100);
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

    void read_expected_all_results(std::string task_number, std::string heliostat_position) {
        std::string round_robin_base = "/round_robin_study/results/phase_I/";
        std::string summary_file = "fluxDataSummary_Task_" + task_number + "_" + heliostat_position + ".csv";
        read_expected_summary_results(std::string(PROJECT_DIR) + round_robin_base + summary_file);
        std::string fluxmap_file = "soltrace_Task_" + task_number + "_" + heliostat_position + "_fluxmap.csv";
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

    bool calculate_receiver_flux_map(SimulationResult result, int nbinsx, int nbinsy, bool is_cylinder,
        bool ignore_direct = false) {
        reset_flux_map();

        double minx = -rec_width / 2.0;
        double maxx = rec_width / 2.0;
        double miny = -rec_height / 2.0;
        double maxy = rec_height / 2.0;

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
            glm::dvec3 local_position;
            glm::dvec3 global_position;

            const ray_record_ptr rr = result[i];
            for (size_t j = 0; j < rr->interactions.size(); j++) {
                if (receiver->get_id() == rr->get_element(j)) {
                    if (rr->get_event(j) == RayEvent::ABSORB) {
                        if (j > 1 || !ignore_direct) // Skip if direct hit
                        {
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
                                if (z <= 0.0)
                                    x = rec_radius * asin(x / rec_radius);
                                else if (z > 0.0) {
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
                                //  qDebug("Not binned: [%d %d],  x=%lg, y=%lg", GridIncrementX, GridIncrementY, x, y);
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

    void check_outputs(SimulationResult result, std::string position) {
        // Check heliostat aim vector and z-rotation
        if (position == "N") {
            EXPECT_NEAR(heliostat->get_aim_vector_ref().x, 0.0, 1.e-3);
            EXPECT_NEAR(heliostat->get_aim_vector_ref().y, -276.838, 1.e-3);
            EXPECT_NEAR(heliostat->get_aim_vector_ref().z, 635.35, 1.e-3);
            EXPECT_NEAR(heliostat->get_zrot(), 180.0, 1.e-4); // TODO: This should be zero
        }
        else if (position == "SE") {
            EXPECT_NEAR(heliostat->get_aim_vector_ref().x, -207.952, 1.e-3);
            EXPECT_NEAR(heliostat->get_aim_vector_ref().y, -125.53, 1.e-3);
            EXPECT_NEAR(heliostat->get_aim_vector_ref().z, 915.611, 1.e-3);
            EXPECT_NEAR(heliostat->get_zrot(), -80.5688, 1.e-4);
        }


        SimulationParameters& params = simData.get_simulation_parameters();
        EXPECT_EQ(helio_hit_count, params.number_of_rays);  // Only if using stages
        EXPECT_EQ(helio_absorb_count + reflect_count, helio_hit_count);
        EXPECT_EQ(rec_absorb_count + miss_count, reflect_count);

        double tol = high_accuracy ? 3.5e-3 : 8.e-3;

        EXPECT_NEAR((double)reflect_count / (double)helio_hit_count, 0.9, tol);

        // Total power absorbed
        double total_power = rec_absorb_count * power_per_ray / 1.e3;   // W to kW
        EXPECT_NEAR(total_power, expected_power, tol * expected_power);

        double peak_tol = high_accuracy ? 1.5e-2 : 0.25;
        if (!high_accuracy) {
            // Low accuracy runs use low-resolution flux map for peak flux comparison
            calculate_receiver_flux_map(result, 30, 30, false);
        }
        else {
            calculate_receiver_flux_map(result, 100, 150, false);
        }
        EXPECT_NEAR(PeakFlux / 1.e3, expected_peak_flux, peak_tol * expected_peak_flux);

        // RMS of flux values
        if (!high_accuracy) calculate_receiver_flux_map(result, 100, 150, false);  // Re-calculate for low-accuracy runs
        EXPECT_EQ(fluxGrid.nrows(), expected_fluxGrid.nrows());
        EXPECT_EQ(fluxGrid.ncols(), expected_fluxGrid.ncols());
        double rmse = 0.0;
        double average_flux = 0.0;
        for (size_t r = 0; r < fluxGrid.nrows(); r++) {
            for (size_t c = 0; c < fluxGrid.ncols(); c++) {
                double sim_flux = fluxGrid.at(r, c) * zScale / 1.e3;
                double exp_flux = expected_fluxGrid.at(r, c);
                rmse += pow(sim_flux - exp_flux, 2);
                average_flux += sim_flux;
            }
        }
        rmse = sqrt(rmse / (fluxGrid.nrows() * fluxGrid.ncols()));
        average_flux /= (fluxGrid.nrows() * fluxGrid.ncols());

        double rmse_tol = high_accuracy ? 0.025 : 0.11;  // 2% of peak flux // 0.04
        EXPECT_LE(rmse / (PeakFlux / 1.e3), rmse_tol);

        if (print_info) {
            std::cout << "Power on receiver: " << total_power << " kW" << std::endl;
            std::cout << "Flux map RMS error: " << rmse << " kW/m2" << std::endl;
            std::cout << "Peak flux: " << PeakFlux / 1.e3 << " kW/m2" << std::endl;
            std::cout << "Average flux: " << average_flux << " kW/m2" << std::endl;
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
        outputFile.close();
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
        outputFile.close();
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

    void simulate_check_outputs(std::string task_number, std::string position) {
        if (print_info) {
            std::cout << "\n\nTask: " << task_number << ", Heliostat Position: " << position << std::endl;
        }

        update_simulation_geometry(solar_azimuth, solar_elevation);

        SimulationResult result;
        simulate(&result);

        calculate_sun_size(result);
        calculate_ray_counts(result);
        read_expected_all_results(task_number, position);
        check_outputs(result, position);

        this->result = result;
        if (!save_results) return;
        std::string filename_postfix = "_Task_" + task_number + "_Position_" + position;
        std::string flux_result_filename = "single_heliostat_fluxmap" + filename_postfix + ".csv";
        save_flux_map_to_file(flux_result_filename);
        std::string flux_comparison_filename = "single_heliostat_fluxmap_comparison" + filename_postfix + ".csv";
        save_flux_comparison_to_file(flux_comparison_filename);
        if (print_info) {
            std::cout << "Flux map saved to: " << flux_result_filename << std::endl;
            std::cout << "Flux comparison results saved to: " << flux_comparison_filename << std::endl;
        }

        if (!save_raydata) return;
        std::string heliostat_filename = "single_heliostat_raydata" + filename_postfix + ".csv";
        result.write_csv_file(heliostat_filename);

        if (print_info) {
            std::cout << "Raydata saved to: " << heliostat_filename << std::endl;
        }

        return;
    }

};

// GTest fixture that reuses the helper logic
// This keeps existing TEST_F-based tests working

template <typename RunnerT>
class SingleHeliostatSimulation : public ::testing::Test, public SingleHeliostatSimulationHelper<RunnerT> {
protected:
    void SetUp() override {
        this->initialize();
    }

    void TearDown() override { }
};

static void write_to_dict(std::string key_name, double val_native,
    double val_optix, std::map<std::string, double>& dict_native,
    std::map<std::string, double>& dict_optix)
{
    dict_native[key_name] = val_native;
    dict_optix[key_name] = val_optix;
}

struct CompareRunnerOptions
{
    int seed = 123;
    bool ignore_direct = true;
    bool save = false;
    bool save_flux = false;

    CompareRunnerOptions() {};

    CompareRunnerOptions(int seed_arg, bool ignore_direct_arg, 
        bool save_arg, bool save_flux_arg)
    {
        seed = seed_arg;
        ignore_direct = ignore_direct_arg;
        save = save_arg;
        save_flux = save_flux_arg;
    }
};

template <typename RunnerA, typename RunnerB>
static void CompareRunners(SingleHeliostatSimulationHelper<RunnerA>& sim_a,
    SingleHeliostatSimulationHelper<RunnerB>& sim_b, int N_rays, const std::string& file_label = "",
    bool skip_a = false, const CompareRunnerOptions& options = {})
{
    double err_frac = 0.01;
    double err_abs = err_frac * (double)N_rays;

    // Run cases
    sim_a.seed = options.seed;
    sim_a.sun_gen_type = SolTrace::Data::GenType::HALTON;
    sim_a.setup_simData();
    sim_a.update_simulation_geometry(sim_a.solar_azimuth, sim_a.solar_elevation);
    SimulationResult result_a;
    if (!skip_a)
        sim_a.simulate(&result_a, N_rays);
    sim_a.calculate_ray_counts(result_a);
    sim_a.calculate_sun_size(result_a);
    sim_a.read_expected_all_results("1a", "N");
    sim_a.calculate_receiver_flux_map(result_a, 30, 30, false, options.ignore_direct);

    sim_b.seed = options.seed;
    sim_b.sun_gen_type = SolTrace::Data::GenType::HALTON;
    sim_b.setup_simData();
    sim_b.update_simulation_geometry(sim_b.solar_azimuth, sim_b.solar_elevation);
    SimulationResult result_b;
    sim_b.simulate(&result_b, N_rays);
    sim_b.calculate_ray_counts(result_b);
    sim_b.calculate_sun_size(result_b);
    sim_b.read_expected_all_results("1a", "N");
    sim_b.calculate_receiver_flux_map(result_b, 30, 30, false, options.ignore_direct);

    std::map<std::string, double> dict_a;
    std::map<std::string, double> dict_b;
    if (options.save_flux)
    {
        std::string file_fluxmap_native = "native_flux_" + file_label + SolTrace::Data::GenTypeMap.at(sim_a.sun_gen_type)
            + "_" + std::to_string(int(N_rays / 1e3)) + "k.csv";
        sim_a.save_flux_map_to_file(file_fluxmap_native);

        std::string file_fluxmap_optix = "optix_flux_" + file_label + SolTrace::Data::GenTypeMap.at(sim_b.sun_gen_type)
            + "_" + std::to_string(int(N_rays / 1e3)) + "k.csv";
        sim_b.save_flux_map_to_file(file_fluxmap_optix);
    }

    // Compare
    EXPECT_EQ(result_a.get_number_of_records(), result_b.get_number_of_records());
    EXPECT_NEAR(sim_a.helio_hit_count, sim_b.helio_hit_count, err_abs);
    EXPECT_NEAR(sim_a.reflect_count, sim_b.reflect_count, err_abs);
    EXPECT_NEAR(sim_a.helio_absorb_count, sim_b.helio_absorb_count, err_abs);
    EXPECT_NEAR(sim_a.rec_absorb_count, sim_b.rec_absorb_count, err_abs);
    EXPECT_NEAR(sim_a.rec_hit_count, sim_b.rec_hit_count, err_abs);
    EXPECT_NEAR(sim_a.rec_direct_hit_count, sim_b.rec_direct_hit_count, err_abs);
    EXPECT_NEAR(sim_a.rec_via_helio_hit_count, sim_b.rec_via_helio_hit_count, err_abs);

    write_to_dict("00_helio_hit_count", sim_a.helio_hit_count, sim_b.helio_hit_count, dict_a, dict_b);
    write_to_dict("01_reflect_count", sim_a.reflect_count, sim_b.reflect_count, dict_a, dict_b);
    write_to_dict("02_helio_absorb_count", sim_a.helio_absorb_count, sim_b.helio_absorb_count, dict_a, dict_b);
    write_to_dict("03_rec_absorb_count", sim_a.rec_absorb_count, sim_b.rec_absorb_count, dict_a, dict_b);
    write_to_dict("04_rec_hit_count", sim_a.rec_hit_count, sim_b.rec_hit_count, dict_a, dict_b);
    write_to_dict("05_rec_direct_hit_count", sim_a.rec_direct_hit_count, sim_b.rec_direct_hit_count, dict_a, dict_b);
    write_to_dict("06_rec_via_helio_hit_count", sim_a.rec_via_helio_hit_count, sim_b.rec_via_helio_hit_count, dict_a, dict_b);

    // Helio rays add up
    EXPECT_EQ(sim_a.reflect_count + sim_a.helio_absorb_count, sim_a.helio_hit_count);
    EXPECT_EQ(sim_b.reflect_count + sim_b.helio_absorb_count, sim_b.helio_hit_count);

    // Receiver rays add up
    EXPECT_EQ(sim_a.rec_direct_hit_count + sim_a.rec_via_helio_hit_count, sim_a.rec_hit_count);
    EXPECT_EQ(sim_b.rec_direct_hit_count + sim_b.rec_via_helio_hit_count, sim_b.rec_hit_count);

    // Reflectivity
    double reflectivity_a = (double)sim_a.reflect_count / (double)sim_a.helio_hit_count;
    double reflectivity_b = (double)sim_b.reflect_count / (double)sim_b.helio_hit_count;
    EXPECT_NEAR(reflectivity_a, 0.9, 0.02);
    EXPECT_NEAR(reflectivity_b, 0.9, 0.02);

    write_to_dict("07_reflectivity", reflectivity_a, reflectivity_b, dict_a, dict_b);

    write_to_dict("08_sun_count", sim_a.sun_ray_count, sim_b.sun_ray_count, dict_a, dict_b);

    // Fraction of hits after helio
    double frac_rec_via_helio_a = (double)sim_a.rec_via_helio_hit_count / (double)sim_a.reflect_count;
    double frac_rec_via_helio_b = (double)sim_b.rec_via_helio_hit_count / (double)sim_b.reflect_count;
    EXPECT_NEAR(frac_rec_via_helio_a, frac_rec_via_helio_b, err_frac);

    write_to_dict("09_frac_rec_via_helio", frac_rec_via_helio_a, frac_rec_via_helio_b, dict_a, dict_b);

    // Compare power per ray
    double power_tol = (5. / (double)N_rays) * 1e3;
    EXPECT_NEAR(sim_a.power_per_ray, sim_b.power_per_ray, power_tol);

    write_to_dict("10_power_per_ray", sim_a.power_per_ray, sim_b.power_per_ray, dict_a, dict_b);

    std::cerr << "Ray Count: " << sim_a.reflect_count << std::endl;

    // Total power absorbed
    double tol = 8.e-3;
    double total_power_a = (double)sim_a.rec_absorb_count * sim_a.power_per_ray * 1.e-3; // [kW]
    double total_power_b = (double)sim_b.rec_absorb_count * sim_b.power_per_ray * 1.e-3; // [kW]
    EXPECT_NEAR(total_power_a, total_power_b, tol * sim_a.expected_power);
    double total_power_diff = total_power_b - total_power_a;
    double total_power_diff_frac = total_power_diff / sim_a.expected_power;

    write_to_dict("11_total_power", total_power_a, total_power_b, dict_a, dict_b);

    // Peak flux
    double peak_tol = 0.25;
    double peak_flux_a = sim_a.PeakFlux / 1.e3;
    double peak_flux_b = sim_b.PeakFlux / 1.e3;
    EXPECT_NEAR(peak_flux_a, peak_flux_b, peak_tol * sim_a.expected_peak_flux);
    double peak_flux_diff = peak_flux_b - peak_flux_a;
    double peak_flux_diff_frac = peak_flux_diff / sim_a.expected_peak_flux;

    write_to_dict("12_peak_flux", peak_flux_a, peak_flux_b, dict_a, dict_b);

    double centroid0_a = sim_a.Centroid[0];
    double centroid0_b = sim_b.Centroid[0];
    double centroid1_a = sim_a.Centroid[1];
    double centroid1_b = sim_b.Centroid[1];

    write_to_dict("13_centroid0", centroid0_a, centroid0_b, dict_a, dict_b);
    write_to_dict("14_centroid1", centroid1_a, centroid1_b, dict_a, dict_b);

    write_to_dict("15_sigmaflux", sim_a.SigmaFlux, sim_b.SigmaFlux, dict_a, dict_b);

    // RMS
    sim_a.calculate_receiver_flux_map(result_a, 100, 150, false, options.ignore_direct);  // Re-calculate for low-accuracy runs
    sim_b.calculate_receiver_flux_map(result_b, 100, 150, false, options.ignore_direct);
    EXPECT_EQ(sim_a.fluxGrid.nrows(), sim_b.fluxGrid.nrows());
    EXPECT_EQ(sim_a.fluxGrid.ncols(), sim_b.fluxGrid.ncols());
    double rmse = 0.0;
    for (size_t r = 0; r < sim_a.fluxGrid.nrows(); r++) {
        for (size_t c = 0; c < sim_a.fluxGrid.ncols(); c++) {
            double flux_a = sim_a.fluxGrid.at(r, c) * sim_a.zScale / 1.e3;
            double flux_b = sim_b.fluxGrid.at(r, c) * sim_b.zScale / 1.e3;
            rmse += pow(flux_a - flux_b, 2);
        }
    }

    // RMS
    rmse = sqrt(rmse / (sim_a.fluxGrid.nrows() * sim_a.fluxGrid.ncols()));
    double rmse_tol = 0.11; // Should be 0.11
    EXPECT_LE(rmse / (peak_flux_a), rmse_tol);


    // Average flux
    EXPECT_NEAR(sim_a.AveFlux / 1000.0, sim_b.AveFlux / 1000.0, rmse_tol);

    write_to_dict("16_average_flux", sim_a.AveFlux / 1000.0, sim_b.AveFlux / 1000.0, dict_a, dict_b);

    // Uniformity
    write_to_dict("17_uniformity", sim_a.Uniformity, sim_b.Uniformity, dict_a, dict_b);

    write_to_dict("18_rmse", rmse, rmse, dict_a, dict_b);
    double rmse_over_peak = rmse / (peak_flux_a);
    write_to_dict("19_rmse_over_peak", rmse_over_peak, rmse_over_peak, dict_a, dict_b);

    if (options.save)
    {
        std::string file_outputs_native = "native_outputs_" + file_label + SolTrace::Data::GenTypeMap.at(sim_a.sun_gen_type)
            + "_" + std::to_string(int(N_rays / 1e3)) + "k.csv";
        sim_a.save_outputs(file_outputs_native, dict_a);

        std::string file_outputs_optix = "optix_outputs_" + file_label + SolTrace::Data::GenTypeMap.at(sim_b.sun_gen_type)
            + "_" + std::to_string(int(N_rays / 1e3)) + "k.csv";
        sim_b.save_outputs(file_outputs_optix, dict_b);
    }


    int x = 0;
}

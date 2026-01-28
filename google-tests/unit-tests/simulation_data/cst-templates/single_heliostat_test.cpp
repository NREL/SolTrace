#include <iomanip>

#include <gtest/gtest.h>

#include <embree_runner.hpp>
#include <native_runner.hpp>
#include <native_runner_types.hpp>
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
using SolTrace::NativeRunner::NativeRunner;
using SolTrace::EmbreeRunner::EmbreeRunner;

using SolTrace::NativeRunner::TRayData;
using SolTrace::NativeRunner::TSystem;
using SolTrace::NativeRunner::TSun;

class SingleHeliostatSimulation : public ::testing::Test {
public:
    bool high_accuracy = false;     // Runs 20 Million rays and tighter tolerance on checks
    bool print_info = false;        // Prints information on from simulation results (sun calculations, ray counts, flux calculations)
    bool save_results = false;      // Saves flux map results to CSV files
    bool save_raydata = false;      // Saves ray data to CSV file

    const Vector3d zero = { 0.0, 0.0, 0.0 }; // Global origin
    const Vector3d khat = { 0.0, 0.0, 1.0 }; // Global z-axis

    double solar_azimuth = 180.0;
    double solar_elevation = 59.96377;

    double rec_radius = 0.0;
    double rec_height = 18.0;
    double rec_width = 12.0;

    std::shared_ptr<SolTrace::Data::Sun> sun;
    std::shared_ptr<Heliostat> heliostat;
    std::shared_ptr<SolTrace::Data::SingleElement> receiver;
protected:

    SimulationData simData;
    //NativeRunner runner;
    EmbreeRunner runner;

    double sun_width;
    double sun_height;
    double A_sun_box;
    double power_per_ray;

    // Ray counts
    uint_fast64_t helio_hit_count;
    uint_fast64_t reflect_count;
    uint_fast64_t helio_absorb_count;
    uint_fast64_t rec_absorb_count;
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

    HPM2D expected_fluxGrid;

    void SetUp() override {
        // Set parameters
        set_default_params();

        // Initialize runner
        RunnerStatus sts = runner.initialize();
        EXPECT_EQ(sts, RunnerStatus::SUCCESS);
        // Native runner speedups
        //runner.disable_power_tower();
        //runner.disable_point_focus();
        runner.set_number_of_threads(10);

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

        // Initial setup of heliostat
        Vector3d heliostat_origin(0.0, 500.0, 5.65);
        Vector3d rec_origin(0.0, 0.0, 169.0);
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
        Vector3d v1 = { 0.0, 1.0, 0.0 }; // Pointing North TODO: change to point towards heliostat
        Vector3d aim_point;
        vector_add(1.0, rec_origin, 1.0, v1, aim_point);
        receiver->set_reference_frame_geometry(rec_origin, aim_point, 0.0);
        receiver->set_name("Receiver");
        receiver->enable();
    }

    void set_high_accuracy_params() {
        SimulationParameters& params = simData.get_simulation_parameters();
        params.number_of_rays = 20.e6;
        params.max_number_of_rays = params.number_of_rays * 100;
    }

    void set_default_params() {
        SimulationParameters& params = simData.get_simulation_parameters();
        params.number_of_rays = 5.e5;
        params.max_number_of_rays = params.number_of_rays * 100;
        params.include_optical_errors = true;
        params.include_sun_shape_errors = true;
        params.seed = 123;
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
        auto ret = st1->add_element(heliostat);
        EXPECT_TRUE(SolTrace::Data::Element::is_success(ret));

        stage_ptr st2 = SolTrace::Data::make_stage(2);
        st2->set_reference_frame_geometry(zero, khat, 0.0);
        ret = st2->add_element(receiver);
        EXPECT_TRUE(SolTrace::Data::Element::is_success(ret));

        simData.add_stage(st1);
        simData.add_stage(st2);
    }

    void set_heliostat_to_southeast() {
        // Update heliostat position to southeast of tower
        Vector3d helio_origin(200.0, -200.0, 5.65);     // Southeast of tower
        heliostat->set_reference_frame_geometry(helio_origin, khat, 0.0);

        // Point receiver to heliostat without tilting down
        helio_origin[2] = 0.0; // Project to ground plane
        Vector3d rec_origin = receiver->get_origin_ref();
        Vector3d aim_point;
        vector_add(1.0, rec_origin, 1.0, helio_origin, aim_point);
        receiver->set_reference_frame_geometry(rec_origin, aim_point, 90.0);
    }

    void set_slant_focal_length() {
        // Set focal length to slant range
        Vector3d distance;
        vector_add(-1.0, receiver->get_origin_global(), 1.0, heliostat->get_origin_global(), distance);
        double focal_length = vector_norm(distance);
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
        Vector3d distance;
        vector_add(-1.0, receiver->get_origin_global(), 1.0, heliostat->get_origin_global(), distance);
        double focal_length = vector_norm(distance);
        heliostat->set_canting(Heliostat::CantingType::ON_AXIS, focal_length, 0.0);
        heliostat->create_geometry();
    }

    void update_simulation_geometry(double azimuth, double elevation) {
        // Set sun position
        Vector3d sun_pos;
        sun_position_vector_degrees(sun_pos, azimuth, elevation);
        sun->set_position(sun_pos);
        // Update heliostat position
        heliostat->update_geometry(azimuth, elevation);
    }

    void simulate(SimulationResult* result) {
        if (high_accuracy) set_high_accuracy_params();
        else set_default_params();

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
        power_per_ray = A_sun_box / sys->SunRayCount * dni;

        if (print_info) {
            std::cout << "Power per ray: " << power_per_ray << std::endl;
            std::cout << "Sun box: " << sun_width << " x " << sun_height << std::endl;
            std::cout << "Sun ray count: " << sys->SunRayCount << std::endl;
        }
    }

    void calculate_ray_counts(SimulationResult result) {
        // Reset counts
        helio_hit_count = 0;
        reflect_count = 0;
        helio_absorb_count = 0;
        rec_absorb_count = 0;
        miss_count = 0;

        for (size_t i = 0; i < result.get_number_of_records(); i++) {
            const ray_record_ptr rr = result[i];

            for (size_t j = 0; j < rr->interactions.size(); j++) {
                auto hit_element = rr->get_element(j);
                SolTrace::Result::RayEvent rev = rr->get_event(j);

                if (rev == RayEvent::EXIT) miss_count++;
                if (hit_element < 0) continue;  // create or exit
                
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

    bool calculate_receiver_flux_map(SimulationResult result, int nbinsx, int nbinsy, bool is_cylinder) {
        reset_flux_map();

        double minx, maxx, miny, maxy;
        minx = -rec_width / 2.0;
        maxx = rec_width / 2.0;
        miny = -rec_height / 2.0;
        maxy = rec_height / 2.0;

        Vector3d rec_origin = receiver->get_origin_global();

        // Autoscale
        if (false) {
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

    void check_outputs(SimulationResult result, std::string position) {
        // Check heliostat aim vector and z-rotation
        if (position == "N") {
            EXPECT_NEAR(heliostat->get_aim_vector_ref().data[0], 0.0, 1.e-3);
            EXPECT_NEAR(heliostat->get_aim_vector_ref().data[1], -276.838, 1.e-3);
            EXPECT_NEAR(heliostat->get_aim_vector_ref().data[2], 635.35, 1.e-3);
            EXPECT_NEAR(heliostat->get_zrot(), 180.0, 1.e-4); // TODO: This should be zero
        }
        else if (position == "SE") {
            EXPECT_NEAR(heliostat->get_aim_vector_ref().data[0], -207.952, 1.e-3);
            EXPECT_NEAR(heliostat->get_aim_vector_ref().data[1], -125.53, 1.e-3);
            EXPECT_NEAR(heliostat->get_aim_vector_ref().data[2], 915.611, 1.e-3);
            EXPECT_NEAR(heliostat->get_zrot(), -80.5688, 1.e-4);
        }


        SimulationParameters& params = simData.get_simulation_parameters();
        EXPECT_EQ(helio_hit_count, params.number_of_rays);
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

    void simulate_check_outputs(std::string task_number, std::string position) {

        if (print_info) {
            std::cout << "\n\nTask: " << task_number << ", Heliostat Position: " << position << std::endl;
        }

        update_simulation_geometry(solar_azimuth, solar_elevation);

        SimulationResult result;
        simulate(&result);

        calculate_sun_size();
        calculate_ray_counts(result);
        read_expected_all_results(task_number, position);
        check_outputs(result, position);

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
    }

    void TearDown() override {
        // Code here will be called immediately after each test (right
        // before the destructor).
    }

};

TEST_F(SingleHeliostatSimulation, SingleFacetFlat_North) 
{
    setup_simData();
    simulate_check_outputs("1a", "N");
    EXPECT_NEAR(sun_width, 15.4557, 1.e-4);
    EXPECT_NEAR(sun_height, 15.4557, 1.e-4);
}

TEST_F(SingleHeliostatSimulation, SingleFacetFlat_Southeast)
{
    set_heliostat_to_southeast();
    setup_simData();
    simulate_check_outputs("1a", "SE");
    EXPECT_NEAR(sun_width, 15.4557, 1.e-4);
    EXPECT_NEAR(sun_height, 15.4557, 1.e-4);
}

TEST_F(SingleHeliostatSimulation, SingleFacetFocused_North)
{
    set_slant_focal_length();
    setup_simData();
    simulate_check_outputs("1b", "N");
    EXPECT_NEAR(sun_width, 15.4557, 1.e-4);
    EXPECT_NEAR(sun_height, 15.4557, 1.e-4);
}

TEST_F(SingleHeliostatSimulation, SingleFacetFocused_Southeast)
{
    set_heliostat_to_southeast();
    set_slant_focal_length();       // reset focal length after moving heliostat
    setup_simData();
    simulate_check_outputs("1b", "SE");
    EXPECT_NEAR(sun_width, 15.4557, 1.e-4);
    EXPECT_NEAR(sun_height, 15.4557, 1.e-4);
}

TEST_F(SingleHeliostatSimulation, MultiFacetFlat_NoCanting_North)
{
    set_flat_multi_facet();
    setup_simData();
    simulate_check_outputs("2", "N");
    EXPECT_NEAR(sun_width, 12.4214, 1.e-4);     // TODO: Why different from single facet?
    EXPECT_NEAR(sun_height, 10.4195, 1.e-4);
}

TEST_F(SingleHeliostatSimulation, MultiFacetFlat_NoCanting_Southeast)
{
    set_flat_multi_facet();
    set_heliostat_to_southeast();
    setup_simData();
    simulate_check_outputs("2", "SE");
    EXPECT_NEAR(sun_width, 11.8574, 1.e-3);
    EXPECT_NEAR(sun_height, 11.5183, 1.e-3);
}

TEST_F(SingleHeliostatSimulation, MultiFacetFlat_SlantCanting_North)
{
    set_onaxis_slant_canting();
    setup_simData();
    simulate_check_outputs("3", "N");
    EXPECT_NEAR(sun_width, 12.4214, 1.e-4);     
    EXPECT_NEAR(sun_height, 10.4236, 1.e-4);  // TODO: Why different than flat facet case?
}

TEST_F(SingleHeliostatSimulation, MultiFacetFlat_SlantCanting_Southeast)
{
    set_heliostat_to_southeast();
    set_onaxis_slant_canting();
    setup_simData();
    simulate_check_outputs("3", "SE");
    EXPECT_NEAR(sun_width, 11.8574, 1.e-3);
    EXPECT_NEAR(sun_height, 11.5183, 1.e-3);
}

TEST_F(SingleHeliostatSimulation, MultiFacetFocused_SlantCanting_North)
{
    set_slant_focal_length();
    set_onaxis_slant_canting();
    setup_simData();
    simulate_check_outputs("4", "N");
    EXPECT_NEAR(sun_width, 12.4214, 1.e-4);     // TODO: Why different from single facet?
    EXPECT_NEAR(sun_height, 10.4236, 1.e-4);
}

TEST_F(SingleHeliostatSimulation, MultiFacetFocused_SlantCanting_Southeast)
{
    set_heliostat_to_southeast();
    set_slant_focal_length();
    set_onaxis_slant_canting();
    setup_simData();
    simulate_check_outputs("4", "SE");
    EXPECT_NEAR(sun_width, 11.8574, 1.e-3);
    EXPECT_NEAR(sun_height, 11.5183, 1.e-3);
}

// TODO: add off-axis cases

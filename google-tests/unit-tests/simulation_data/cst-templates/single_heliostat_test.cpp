#include <iomanip>

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

// TODO: Refactor to remove duplicate code, read in files, create flux maps and compare.

class SingleHeliostatSimulation : public ::testing::Test {
public:
    const Vector3d zero = { 0.0, 0.0, 0.0 }; // Global origin
    const Vector3d khat = { 0.0, 0.0, 1.0 }; // Global z-axis

    double solar_azimuth = 180.0;
    double solar_elevation = 59.96377;

    std::shared_ptr<Heliostat> heliostat;
    std::shared_ptr<SolTrace::Data::SingleElement> receiver;
protected:

    SimulationData simData;
    NativeRunner runner;
    SimulationResult result;

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


    void SetUp() override {
        // Set parameters
        SimulationParameters& params = simData.get_simulation_parameters();
        params.number_of_rays = 1.e5;
        params.max_number_of_rays = params.number_of_rays * 100;
        params.include_optical_errors = true;
        params.include_sun_shape_errors = true;
        params.seed = 123;

        runner.disable_power_tower();
        runner.disable_point_focus();

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
        receiver->set_aperture(SolTrace::Data::make_aperture<SolTrace::Data::Rectangle>(12.0, 18.0));
        receiver->set_surface(SolTrace::Data::make_surface<SolTrace::Data::Flat>());
        Vector3d v1 = { 0.0, 1.0, 0.0 }; // Pointing North TODO: change to point towards heliostat
        Vector3d aim_point;
        vector_add(1.0, rec_origin, 1.0, v1, aim_point);
        receiver->set_reference_frame_geometry(rec_origin, aim_point, 0.0);
        receiver->set_name("Receiver");
        receiver->enable();
    }

    void setup_simData() {
        // Set up sun
        Vector3d sun_pos;
        sun_position_vector_degrees(sun_pos, solar_azimuth, solar_elevation);
        {
            auto sun = SolTrace::Data::make_ray_source<Sun>();
            sun->set_position(sun_pos);
            sun->set_shape(SolTrace::Data::SunShape::PILLBOX, 0.0, 4.65, 0.0);
            simData.add_ray_source(sun);
        }

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

    void simulate() {
        RunnerStatus sts = runner.initialize();
        EXPECT_EQ(sts, RunnerStatus::SUCCESS);
        // Setup runs but is not complete
        sts = runner.setup_simulation(&simData);
        EXPECT_EQ(sts, RunnerStatus::SUCCESS);
        sts = runner.run_simulation();
        EXPECT_EQ(sts, RunnerStatus::SUCCESS);
        sts = runner.report_simulation(&result, 0);
        EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    }

    void calculate_sun_size(double dni, bool print_sun_info) {
        const TSystem* sys = runner.get_system();
        const TSun* sun = &(sys->Sun);
        sun_width = (sun->MaxXSun - sun->MinXSun);
        sun_height = (sun->MaxYSun - sun->MinYSun);
        A_sun_box = sun_width * sun_height;
        power_per_ray = A_sun_box / sys->SunRayCount * dni;

        if (print_sun_info) {
            std::cout << "Power per ray: " << power_per_ray << std::endl;
            std::cout << "Sun box: " << sun_width << " x " << sun_height << std::endl;
            std::cout << "Sun ray count: " << sys->SunRayCount << std::endl;
        }
    }

    void calculate_ray_counts(bool print_info) {
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

    bool calculate_receiver_flux_map(
        int nbinsx, int nbinsy, bool print_info) {
        reset_flux_map();

        double minx, maxx, miny, maxy;
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
            || maxx <= minx || maxy <= miny) {
            return false;
        }

        double gridszx = maxx - minx;
        double gridszy = maxy - miny;

        binszx = gridszx / (double)nbinsx;
        binszy = gridszy / (double)nbinsy;

        xValues.resize(nbinsx);
        yValues.resize(nbinsy);
        fluxGrid.resize(nbinsx, nbinsy);

        // Bin rays here -> BinRaysXY()
        double zval = 0;
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

                        zval = z;

                        GridIncrementX = -1; //initialize grid increment counters
                        GridIncrementY = -1;

                        //determine which bin the ray falls into
                        while ((minx + (GridIncrementX + 1) * binszx) < x)
                            GridIncrementX++;

                        while ((miny + (GridIncrementY + 1) * binszy) < y)
                            GridIncrementY++;

                        if (GridIncrementX >= 0 && GridIncrementX < (int)fluxGrid.nrows()
                            && GridIncrementY >= 0 && GridIncrementY < (int)fluxGrid.ncols())
                        {
                            fluxGrid.at(GridIncrementX, GridIncrementY) += 1;//if ray falls inside a bin, increment count for that bin
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

                if (z < MinFlux){
                    MinFlux = z;
                    NRaysInMinFluxBin = (int)fluxGrid.at(r, c);
                }
            }
        }

        AveFlux = SumFlux / (nbinsx * nbinsy);
        SigmaFlux = sqrt((nbinsx * nbinsy * SumFlux2 - SumFlux * SumFlux) / (nbinsx * nbinsy * nbinsx * nbinsy));
        Uniformity = SigmaFlux / AveFlux;
        PeakFluxUncertainty = 100 / sqrt((double)NRaysInPeakFluxBin);
        AveFluxUncertainty = 100 / sqrt((double)result.get_number_of_records());     // TODO: Should the be number of rays traced not total interactions?

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

    void TearDown() override {
        // Code here will be called immediately after each test (right
        // before the destructor).
    }

};

TEST_F(SingleHeliostatSimulation, SingleFacetFlat_North) 
{
    bool high_accuracy = false;
    bool print_info = false;

    SimulationParameters& params = simData.get_simulation_parameters();
    if (high_accuracy) {
        params.number_of_rays = 20.e6;
        params.max_number_of_rays = params.number_of_rays * 100;
    }
    heliostat->update_geometry(solar_azimuth, solar_elevation);

    // Check heliostat aim vector and z-rotation
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[0], 0.0, 1.e-3);
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[1], -276.838, 1.e-3);
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[2], 635.35, 1.e-3);
    EXPECT_NEAR(heliostat->get_zrot(), 180.0, 1.e-4); // TODO: This should be zero

    setup_simData();
    simulate();
    //result.write_csv_file("singlefacetflat_north_raydata.csv");

    double dni = 1000.0;
    calculate_sun_size(dni, print_info);
    EXPECT_NEAR(sun_width, 15.4557, 1.e-4);
    EXPECT_NEAR(sun_height, 15.4557, 1.e-4);

    calculate_ray_counts(print_info);
    double total_power = rec_absorb_count * power_per_ray;

    double tol = high_accuracy? 5.e-4 : 5.e-3;
    //double expected_power = 99967.3;    // No sunshape or surfaces errors
    //double expected_power = 93219.6;    // Sunshape only
    //double expected_power = 88033.6;    // Surface errors only
    double expected_power = 84984.8;      // Sunshape + surface errors

    uint_fast64_t NRAYS = params.number_of_rays;
    EXPECT_EQ(helio_hit_count, NRAYS);
    EXPECT_EQ(helio_absorb_count + reflect_count, helio_hit_count);
    EXPECT_EQ(rec_absorb_count + miss_count, reflect_count);

    EXPECT_NEAR((double)reflect_count / (double)helio_hit_count, 0.9, tol);
    EXPECT_NEAR(total_power, expected_power, tol * expected_power);

    if (high_accuracy) {
        calculate_receiver_flux_map(100, 150, print_info);
        double expected_peak = 853.65157013;
        EXPECT_NEAR(PeakFlux, expected_peak, 2.0);
        // TODO: Create a flux map and compare to expected distribution
    }
    else {
        calculate_receiver_flux_map(30, 30, print_info);
        double expected_peak = 906.458;
        EXPECT_NEAR(PeakFlux, expected_peak, 30.0);
    }
}
TEST_F(SingleHeliostatSimulation, SingleFacetFlat_Southeast)
{
    bool high_accuracy = false;
    bool print_info = false;

    SimulationParameters& params = simData.get_simulation_parameters();
    if (high_accuracy) {
        params.number_of_rays = 20.e6;
        params.max_number_of_rays = params.number_of_rays * 100;
    }
    // Update heliostat position to southeast of tower
    Vector3d helio_origin(200.0, -200.0, 5.65);     // Southeast of tower
    heliostat->set_reference_frame_geometry(helio_origin, khat, 0.0);
    
    // Point receiver to heliostat without tilting down
    helio_origin[2] = 0.0; // Project to ground plane
    Vector3d rec_origin = receiver->get_origin_ref();
    Vector3d aim_point;
    vector_add(1.0, rec_origin, 1.0, helio_origin, aim_point);
    receiver->set_reference_frame_geometry(rec_origin, aim_point, 90.0);
    heliostat->update_geometry(solar_azimuth, solar_elevation);

    // Check heliostat aim vector and z-rotation
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[0], -207.952, 1.e-3);
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[1], -125.53, 1.e-3);
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[2], 915.611, 1.e-3);
    EXPECT_NEAR(heliostat->get_zrot(), -80.5688, 1.e-4);

    setup_simData();
    simulate();
    //result.write_csv_file("singlefacetflat_southeast_raydata.csv");

    double dni = 1000.0;
    calculate_sun_size(dni, print_info);
    EXPECT_NEAR(sun_width, 15.4557, 1.e-4);
    EXPECT_NEAR(sun_height, 15.4557, 1.e-4);

    calculate_ray_counts(print_info);
    double total_power = rec_absorb_count * power_per_ray;

    double tol = high_accuracy ? 5.e-4 : 5.e-3;
    //double expected_power = 77251.8;    // No sunshape or surfaces errors
    //double expected_power = 76345.1;    // Sunshape only
    //double expected_power = 75585.1;    // Surface errors only
    double expected_power = 74695.7;      // Sunshape + surface errors

    uint_fast64_t NRAYS = params.number_of_rays;
    EXPECT_EQ(helio_hit_count, NRAYS);
    EXPECT_EQ(helio_absorb_count + reflect_count, helio_hit_count);
    EXPECT_EQ(rec_absorb_count + miss_count, reflect_count);

    EXPECT_NEAR((double)reflect_count / (double)helio_hit_count, 0.9, tol);
    EXPECT_NEAR(total_power, expected_power, tol * expected_power);

    if (high_accuracy) {
        calculate_receiver_flux_map(100, 150, print_info);
        double expected_peak = 824.965626;
        EXPECT_NEAR(PeakFlux, expected_peak, 10.0);
        // TODO: Create a flux map and compare to expected distribution
    }
    else {
        calculate_receiver_flux_map(30, 30, print_info);
        double expected_peak = 921.871;
        EXPECT_NEAR(PeakFlux, expected_peak, 30.0);
    }
}
TEST_F(SingleHeliostatSimulation, SingleFacetFocused_North)
{
    bool high_accuracy = false;
    bool print_info = false;

    SimulationParameters& params = simData.get_simulation_parameters();
    if (high_accuracy) {
        params.number_of_rays = 20.e6;
        params.max_number_of_rays = params.number_of_rays * 100;
    }

    // Set focal length to slant range
    Vector3d distance;
    vector_add(-1.0, receiver->get_origin_global(), 1.0, heliostat->get_origin_global(), distance);
    double focal_length = vector_norm(distance);
    heliostat->set_focal_length(focal_length);
    heliostat->create_geometry();
    heliostat->update_geometry(solar_azimuth, solar_elevation);

    // Check heliostat aim vector and z-rotation
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[0], 0.0, 1.e-3);
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[1], -276.838, 1.e-3);
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[2], 635.35, 1.e-3);
    EXPECT_NEAR(heliostat->get_zrot(), 180.0, 1.e-4); // TODO: This should be zero

    setup_simData();
    simulate();
    //result.write_csv_file("singlefacetfocused_north_raydata.csv");

    double dni = 1000.0;
    calculate_sun_size(dni, print_info);
    EXPECT_NEAR(sun_width, 15.4557, 1.e-4);
    EXPECT_NEAR(sun_height, 15.4557, 1.e-4);

    calculate_ray_counts(print_info);
    double total_power = rec_absorb_count * power_per_ray;

    double tol = high_accuracy ? 5.e-4 : 5.e-3;
    double expected_power = 99057.1;      // Sunshape + surface errors

    uint_fast64_t NRAYS = params.number_of_rays;
    EXPECT_EQ(helio_hit_count, NRAYS);
    EXPECT_EQ(helio_absorb_count + reflect_count, helio_hit_count);
    EXPECT_EQ(rec_absorb_count + miss_count, reflect_count);

    EXPECT_NEAR((double)reflect_count / (double)helio_hit_count, 0.9, tol);
    EXPECT_NEAR(total_power, expected_power, tol * expected_power);

    if (high_accuracy) {
        calculate_receiver_flux_map(100, 150, print_info);
        double expected_peak = 2634.595862;

        EXPECT_NEAR(PeakFlux, expected_peak, 10.0);
        // TODO: Create a flux map and compare to expected distribution
    }
    else {
        calculate_receiver_flux_map(30, 30, print_info);
        double expected_peak = 2676.09;
        EXPECT_NEAR(PeakFlux, expected_peak, 30.0);
    }
}
TEST_F(SingleHeliostatSimulation, SingleFacetFocused_Southeast)
{
    bool high_accuracy = false;
    bool print_info = false;

    SimulationParameters& params = simData.get_simulation_parameters();
    if (high_accuracy) {
        params.number_of_rays = 20.e6;
        params.max_number_of_rays = params.number_of_rays * 100;
    }
    // Update heliostat position to southeast of tower
    Vector3d helio_origin(200.0, -200.0, 5.65);     // Southeast of tower
    heliostat->set_reference_frame_geometry(helio_origin, khat, 0.0);

    // Set focal length to slant range
    Vector3d distance;
    vector_add(-1.0, receiver->get_origin_global(), 1.0, heliostat->get_origin_global(), distance);
    double focal_length = vector_norm(distance);
    heliostat->set_focal_length(focal_length);
    heliostat->create_geometry();       // Re-create geometry with new focal length

    // Point receiver to heliostat without tilting down
    helio_origin[2] = 0.0; // Project to ground plane
    Vector3d rec_origin = receiver->get_origin_ref();
    Vector3d aim_point;
    vector_add(1.0, rec_origin, 1.0, helio_origin, aim_point);
    receiver->set_reference_frame_geometry(rec_origin, aim_point, 90.0);

    heliostat->update_geometry(solar_azimuth, solar_elevation);

    // Check heliostat aim vector and z-rotation
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[0], -207.952, 1.e-3);
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[1], -125.53, 1.e-3);
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[2], 915.611, 1.e-3);
    EXPECT_NEAR(heliostat->get_zrot(), -80.5688, 1.e-4);

    setup_simData();
    simulate();
    //result.write_csv_file("singlefacetfocused_southeast_raydata.csv");

    double dni = 1000.0;
    calculate_sun_size(dni, print_info);
    EXPECT_NEAR(sun_width, 15.4557, 1.e-4);
    EXPECT_NEAR(sun_height, 15.4557, 1.e-4);

    calculate_ray_counts(print_info);
    double total_power = rec_absorb_count * power_per_ray;

    double tol = high_accuracy ? 5.e-4 : 5.e-3;
    double expected_power = 80322.3;      // Sunshape + surface errors

    uint_fast64_t NRAYS = params.number_of_rays;
    EXPECT_EQ(helio_hit_count, NRAYS);
    EXPECT_EQ(helio_absorb_count + reflect_count, helio_hit_count);
    EXPECT_EQ(rec_absorb_count + miss_count, reflect_count);

    EXPECT_NEAR((double)reflect_count / (double)helio_hit_count, 0.9, tol);
    EXPECT_NEAR(total_power, expected_power, tol * expected_power);

    if (high_accuracy) {
        calculate_receiver_flux_map(100, 150, print_info);
        double expected_peak = 4273.82;

        EXPECT_NEAR(PeakFlux, expected_peak, 10.0);
        // TODO: Create a flux map and compare to expected distribution
    }
    else {
        calculate_receiver_flux_map(30, 30, print_info);
        double expected_peak = 4386.455;
        EXPECT_NEAR(PeakFlux, expected_peak, 30.0);
    }
}
TEST_F(SingleHeliostatSimulation, MultiFacetFlat_NoCanting_North)
{
    bool high_accuracy = false;
    bool print_info = false;

    SimulationParameters& params = simData.get_simulation_parameters();
    if (high_accuracy) {
        params.number_of_rays = 20.e6;
        params.max_number_of_rays = params.number_of_rays * 100;
    }
    // Multi-facet heliostat
    heliostat->set_number_panels(7, 5);
    heliostat->set_gaps(0.03, 0.03);
    heliostat->create_geometry();

    heliostat->update_geometry(solar_azimuth, solar_elevation);

    // Check heliostat aim vector and z-rotation
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[0], 0.0, 1.e-3);
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[1], -276.838, 1.e-3);
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[2], 635.35, 1.e-3);
    EXPECT_NEAR(heliostat->get_zrot(), 180.0, 1.e-4); // TODO: This should be zero

    setup_simData();
    simulate();
    //result.write_csv_file("raydata.csv");

    double dni = 1000.0;
    calculate_sun_size(dni, print_info);
    EXPECT_NEAR(sun_width, 12.4214, 1.e-4);
    EXPECT_NEAR(sun_height, 10.4195, 1.e-4);

    calculate_ray_counts(print_info);
    double total_power = rec_absorb_count * power_per_ray;

    double tol = high_accuracy ? 5.e-4 : 5.e-3;
    double expected_power = 82587.50674;

    uint_fast64_t NRAYS = params.number_of_rays;
    EXPECT_EQ(helio_hit_count, NRAYS);
    EXPECT_EQ(helio_absorb_count + reflect_count, helio_hit_count);
    EXPECT_EQ(rec_absorb_count + miss_count, reflect_count);

    EXPECT_NEAR((double)reflect_count / (double)helio_hit_count, 0.9, tol);
    EXPECT_NEAR(total_power, expected_power, tol * expected_power);

    if (high_accuracy) {
        calculate_receiver_flux_map(100, 150, print_info);
        double expected_peak = 827.775965;
        EXPECT_NEAR(PeakFlux, expected_peak, 10.0);
        // TODO: Create a flux map and compare to expected distribution
    }
    else {
        calculate_receiver_flux_map(30, 30, print_info);
        double expected_peak = 906.458;
        EXPECT_NEAR(PeakFlux, expected_peak, 30.0);
    }
}

TEST_F(SingleHeliostatSimulation, MultiFacetFlat_NoCanting_Southeast)
{
    bool high_accuracy = false;
    bool print_info = false;

    SimulationParameters& params = simData.get_simulation_parameters();
    if (high_accuracy) {
        params.number_of_rays = 20.e6;
        params.max_number_of_rays = params.number_of_rays * 100;
    }
    // Update heliostat position to southeast of tower
    Vector3d helio_origin(200.0, -200.0, 5.65);     // Southeast of tower
    heliostat->set_reference_frame_geometry(helio_origin, khat, 0.0);

    // Multi-facet heliostat
    heliostat->set_number_panels(7, 5);
    heliostat->set_gaps(0.03, 0.03);
    heliostat->create_geometry();

    // Point receiver to heliostat without tilting down
    helio_origin[2] = 0.0; // Project to ground plane
    Vector3d rec_origin = receiver->get_origin_ref();
    Vector3d aim_point;
    vector_add(1.0, rec_origin, 1.0, helio_origin, aim_point);
    receiver->set_reference_frame_geometry(rec_origin, aim_point, 90.0);

    heliostat->update_geometry(solar_azimuth, solar_elevation);

    // Check heliostat aim vector and z-rotation
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[0], -207.952, 1.e-3);
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[1], -125.53, 1.e-3);
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[2], 915.611, 1.e-3);
    EXPECT_NEAR(heliostat->get_zrot(), -80.5688, 1.e-4);

    setup_simData();
    simulate();
    //result.write_csv_file("raydata.csv");

    double dni = 1000.0;
    calculate_sun_size(dni, print_info);
    EXPECT_NEAR(sun_width, 11.8574, 1.e-3);
    EXPECT_NEAR(sun_height, 11.5183, 1.e-3);

    calculate_ray_counts(print_info);
    double total_power = rec_absorb_count * power_per_ray;

    double tol = high_accuracy ? 3.e-3 : 5.e-3; // TODO: Understand why this case is less accurate
    double expected_power = 72486.83912;

    uint_fast64_t NRAYS = params.number_of_rays;
    EXPECT_EQ(helio_hit_count, NRAYS);
    EXPECT_EQ(helio_absorb_count + reflect_count, helio_hit_count);
    EXPECT_EQ(rec_absorb_count + miss_count, reflect_count);

    EXPECT_NEAR((double)reflect_count / (double)helio_hit_count, 0.9, tol);
    EXPECT_NEAR(total_power, expected_power, tol * expected_power);

    if (high_accuracy) {
        calculate_receiver_flux_map(100, 150, print_info);
        double expected_peak = 834.794525;
        EXPECT_NEAR(PeakFlux, expected_peak, 35.0); // TODO: Understand why this case is less accurate
        // TODO: Create a flux map and compare to expected distribution
    }
    else {
        calculate_receiver_flux_map(30, 30, print_info);
        double expected_peak = 877.152;
        EXPECT_NEAR(PeakFlux, expected_peak, 30.0);
    }
}

TEST_F(SingleHeliostatSimulation, MultiFacetFlat_OnAxisSlantCanting_North)
{
    bool high_accuracy = false;
    bool print_info = false;

    SimulationParameters& params = simData.get_simulation_parameters();
    if (high_accuracy) {
        params.number_of_rays = 20.e6;
        params.max_number_of_rays = params.number_of_rays * 100;
    }
    // Multi-facet heliostat
    heliostat->set_number_panels(7, 5);
    heliostat->set_gaps(0.03, 0.03);
    // Set on-axis canting to slant range
    Vector3d distance;
    vector_add(-1.0, receiver->get_origin_global(), 1.0, heliostat->get_origin_global(), distance);
    double focal_length = vector_norm(distance);
    heliostat->set_canting(Heliostat::CantingType::ON_AXIS, focal_length, 0.0);
    heliostat->create_geometry();

    heliostat->update_geometry(solar_azimuth, solar_elevation);

    // Check heliostat aim vector and z-rotation
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[0], 0.0, 1.e-3);
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[1], -276.838, 1.e-3);
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[2], 635.35, 1.e-3);
    EXPECT_NEAR(heliostat->get_zrot(), 180.0, 1.e-4); // TODO: This should be zero

    setup_simData();
    simulate();
    //result.write_csv_file("raydata.csv");

    double dni = 1000.0;
    calculate_sun_size(dni, print_info);
    EXPECT_NEAR(sun_width, 12.4214, 1.e-4);
    EXPECT_NEAR(sun_height, 10.4236, 1.e-4);

    calculate_ray_counts(print_info);
    double total_power = rec_absorb_count * power_per_ray;

    double tol = high_accuracy ? 5.e-4 : 5.e-3;
    double expected_power = 96208.97147;

    uint_fast64_t NRAYS = params.number_of_rays;
    EXPECT_EQ(helio_hit_count, NRAYS);
    EXPECT_EQ(helio_absorb_count + reflect_count, helio_hit_count);
    EXPECT_EQ(rec_absorb_count + miss_count, reflect_count);

    EXPECT_NEAR((double)reflect_count / (double)helio_hit_count, 0.9, tol);
    EXPECT_NEAR(total_power, expected_power, tol * expected_power);

    if (high_accuracy) {
        calculate_receiver_flux_map(100, 150, print_info);
        double expected_peak = 2462.377564;
        EXPECT_NEAR(PeakFlux, expected_peak, 25.0);
        // TODO: Create a flux map and compare to expected distribution
    }
    else {
        calculate_receiver_flux_map(30, 30, print_info);
        double expected_peak = 2499.6;
        EXPECT_NEAR(PeakFlux, expected_peak, 30.0);
    }
}

TEST_F(SingleHeliostatSimulation, MultiFacetFlat_OnAxisSlantCanting_Southeast)
{
    bool high_accuracy = false;
    bool print_info = false;

    SimulationParameters& params = simData.get_simulation_parameters();
    if (high_accuracy) {
        params.number_of_rays = 20.e6;
        params.max_number_of_rays = params.number_of_rays * 100;
    }
    // Update heliostat position to southeast of tower
    Vector3d helio_origin(200.0, -200.0, 5.65);     // Southeast of tower
    heliostat->set_reference_frame_geometry(helio_origin, khat, 0.0);

    // Multi-facet heliostat
    heliostat->set_number_panels(7, 5);
    heliostat->set_gaps(0.03, 0.03);
    // Set on-axis canting to slant range
    Vector3d distance;
    vector_add(-1.0, receiver->get_origin_global(), 1.0, heliostat->get_origin_global(), distance);
    double focal_length = vector_norm(distance);
    heliostat->set_canting(Heliostat::CantingType::ON_AXIS, focal_length, 0.0);
    heliostat->create_geometry();

    // Point receiver to heliostat without tilting down
    helio_origin[2] = 0.0; // Project to ground plane
    Vector3d rec_origin = receiver->get_origin_ref();
    Vector3d aim_point;
    vector_add(1.0, rec_origin, 1.0, helio_origin, aim_point);
    receiver->set_reference_frame_geometry(rec_origin, aim_point, 90.0);

    heliostat->update_geometry(solar_azimuth, solar_elevation);

    // Check heliostat aim vector and z-rotation
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[0], -207.952, 1.e-3);
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[1], -125.53, 1.e-3);
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[2], 915.611, 1.e-3);
    EXPECT_NEAR(heliostat->get_zrot(), -80.5688, 1.e-4);

    setup_simData();
    simulate();
    //result.write_csv_file("raydata.csv");

    double dni = 1000.0;
    calculate_sun_size(dni, print_info);
    EXPECT_NEAR(sun_width, 11.8574, 1.e-3);
    EXPECT_NEAR(sun_height, 11.5183, 1.e-3);

    calculate_ray_counts(print_info);
    double total_power = rec_absorb_count * power_per_ray;

    double tol = high_accuracy ? 3.e-3 : 5.e-3;
    double expected_power = 78140.42109;

    uint_fast64_t NRAYS = params.number_of_rays;
    EXPECT_EQ(helio_hit_count, NRAYS);
    EXPECT_EQ(helio_absorb_count + reflect_count, helio_hit_count);
    EXPECT_EQ(rec_absorb_count + miss_count, reflect_count);

    EXPECT_NEAR((double)reflect_count / (double)helio_hit_count, 0.9, tol);
    EXPECT_NEAR(total_power, expected_power, tol * expected_power);

    if (high_accuracy) {
        calculate_receiver_flux_map(100, 150, print_info);
        double expected_peak = 3835.962992;
        EXPECT_NEAR(PeakFlux, expected_peak, 25.0);
        // TODO: Create a flux map and compare to expected distribution
    }
    else {
        calculate_receiver_flux_map(30, 30, print_info);
        double expected_peak = 3835.962992;
        EXPECT_NEAR(PeakFlux, expected_peak, 150.0);
    }
}

TEST_F(SingleHeliostatSimulation, MultiFacetFocused_OnAxisSlantCanting_North)
{
    bool high_accuracy = false;
    bool print_info = false;

    SimulationParameters& params = simData.get_simulation_parameters();
    if (high_accuracy) {
        params.number_of_rays = 20.e6;
        params.max_number_of_rays = params.number_of_rays * 100;
    }
    // Multi-facet heliostat
    heliostat->set_number_panels(7, 5);
    heliostat->set_gaps(0.03, 0.03);
    // Set on-axis canting to slant range
    Vector3d distance;
    vector_add(-1.0, receiver->get_origin_global(), 1.0, heliostat->get_origin_global(), distance);
    double focal_length = vector_norm(distance);
    heliostat->set_canting(Heliostat::CantingType::ON_AXIS, focal_length, 0.0);
    heliostat->set_focal_length(focal_length);
    heliostat->create_geometry();

    heliostat->update_geometry(solar_azimuth, solar_elevation);

    // Check heliostat aim vector and z-rotation
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[0], 0.0, 1.e-3);
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[1], -276.838, 1.e-3);
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[2], 635.35, 1.e-3);
    EXPECT_NEAR(heliostat->get_zrot(), 180.0, 1.e-4); // TODO: This should be zero

    setup_simData();
    simulate();
    //result.write_csv_file("raydata.csv");

    double dni = 1000.0;
    calculate_sun_size(dni, print_info);
    EXPECT_NEAR(sun_width, 12.4214, 1.e-4);
    EXPECT_NEAR(sun_height, 10.4236, 1.e-4);

    calculate_ray_counts(print_info);
    double total_power = rec_absorb_count * power_per_ray;

    double tol = high_accuracy ? 5.e-4 : 5.e-3;
    double expected_power = 96370.82753;

    uint_fast64_t NRAYS = params.number_of_rays;
    EXPECT_EQ(helio_hit_count, NRAYS);
    EXPECT_EQ(helio_absorb_count + reflect_count, helio_hit_count);
    EXPECT_EQ(rec_absorb_count + miss_count, reflect_count);

    EXPECT_NEAR((double)reflect_count / (double)helio_hit_count, 0.9, tol);
    EXPECT_NEAR(total_power, expected_power, tol * expected_power);

    if (high_accuracy) {
        calculate_receiver_flux_map(100, 150, print_info);
        double expected_peak = 2583.58203;
        EXPECT_NEAR(PeakFlux, expected_peak, 25.0);
        // TODO: Create a flux map and compare to expected distribution
    }
    else {
        calculate_receiver_flux_map(30, 30, print_info);
        double expected_peak = 2499.6;
        EXPECT_NEAR(PeakFlux, expected_peak, 30.0);
    }
}

TEST_F(SingleHeliostatSimulation, MultiFacetFocused_OnAxisSlantCanting_Southeast)
{
    bool high_accuracy = false;
    bool print_info = false;

    SimulationParameters& params = simData.get_simulation_parameters();
    if (high_accuracy) {
        params.number_of_rays = 20.e6;
        params.max_number_of_rays = params.number_of_rays * 100;
    }
    // Update heliostat position to southeast of tower
    Vector3d helio_origin(200.0, -200.0, 5.65);     // Southeast of tower
    heliostat->set_reference_frame_geometry(helio_origin, khat, 0.0);

    // Multi-facet heliostat
    heliostat->set_number_panels(7, 5);
    heliostat->set_gaps(0.03, 0.03);
    // Set on-axis canting to slant range
    Vector3d distance;
    vector_add(-1.0, receiver->get_origin_global(), 1.0, heliostat->get_origin_global(), distance);
    double focal_length = vector_norm(distance);
    heliostat->set_canting(Heliostat::CantingType::ON_AXIS, focal_length, 0.0);
    heliostat->set_focal_length(focal_length);
    heliostat->create_geometry();

    // Point receiver to heliostat without tilting down
    helio_origin[2] = 0.0; // Project to ground plane
    Vector3d rec_origin = receiver->get_origin_ref();
    Vector3d aim_point;
    vector_add(1.0, rec_origin, 1.0, helio_origin, aim_point);
    receiver->set_reference_frame_geometry(rec_origin, aim_point, 90.0);

    heliostat->update_geometry(solar_azimuth, solar_elevation);

    // Check heliostat aim vector and z-rotation
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[0], -207.952, 1.e-3);
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[1], -125.53, 1.e-3);
    EXPECT_NEAR(heliostat->get_aim_vector_ref().data[2], 915.611, 1.e-3);
    EXPECT_NEAR(heliostat->get_zrot(), -80.5688, 1.e-4);

    setup_simData();
    simulate();
    //result.write_csv_file("raydata.csv");

    double dni = 1000.0;
    calculate_sun_size(dni, print_info);
    EXPECT_NEAR(sun_width, 11.8574, 1.e-3);
    EXPECT_NEAR(sun_height, 11.5183, 1.e-3);

    calculate_ray_counts(print_info);
    double total_power = rec_absorb_count * power_per_ray;

    double tol = high_accuracy ? 3.e-3 : 5.e-3;
    double expected_power = 78143.54596;

    uint_fast64_t NRAYS = params.number_of_rays;
    EXPECT_EQ(helio_hit_count, NRAYS);
    EXPECT_EQ(helio_absorb_count + reflect_count, helio_hit_count);
    EXPECT_EQ(rec_absorb_count + miss_count, reflect_count);

    EXPECT_NEAR((double)reflect_count / (double)helio_hit_count, 0.9, tol);
    EXPECT_NEAR(total_power, expected_power, tol * expected_power);

    if (high_accuracy) {
        calculate_receiver_flux_map(100, 150, print_info);
        double expected_peak = 4132.093047;
        EXPECT_NEAR(PeakFlux, expected_peak, 25.0);
        // TODO: Create a flux map and compare to expected distribution
    }
    else {
        calculate_receiver_flux_map(30, 30, print_info);
        double expected_peak = 4132.093047;
        EXPECT_NEAR(PeakFlux, expected_peak, 150.0);
    }
}

// TODO: add off-axis cases
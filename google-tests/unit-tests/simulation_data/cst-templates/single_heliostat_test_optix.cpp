#include "single_heliostat_test_template.hpp"

#include <optix_runner.hpp>
#include <native_runner.hpp>

using SolTrace::NativeRunner::NativeRunner;
using OptixRunnerType = OptixRunner;
using NativeRunnerType = NativeRunner;

constexpr int N_rays_glob = 2e6;
constexpr int seed = 123;
constexpr bool save = false;
constexpr bool save_flux = false;
constexpr bool ignore_direct = true;

int N_threads = static_cast<int>(std::max(1u, std::min(std::thread::hardware_concurrency(), 10u)));

static void write_to_dict(std::string key_name, double val_native,
    double val_optix, std::map<std::string, double>& dict_native,
    std::map<std::string, double>& dict_optix)
{
    dict_native[key_name] = val_native;
    dict_optix[key_name] = val_optix;
}

static void CompareRunners(SingleHeliostatSimulationHelper<NativeRunner>& sim_native,
    SingleHeliostatSimulationHelper<OptixRunner>& sim_optix, int N_rays, const std::string& file_label = "",
    bool skip_native = false)
{
    double err_frac = 0.01;
    double err_abs = err_frac * (double)N_rays;

    // Run cases
    sim_native.seed = seed;
    sim_native.sun_gen_type = SolTrace::Data::GenType::HALTON;
    sim_native.setup_simData();
    sim_native.update_simulation_geometry(sim_native.solar_azimuth, sim_native.solar_elevation);
    SimulationResult result_native;
    if(!skip_native)
        sim_native.simulate(&result_native, N_rays);
    sim_native.calculate_ray_counts(result_native);
    sim_native.calculate_sun_size(result_native);
    sim_native.read_expected_all_results("1a", "N");
    sim_native.calculate_receiver_flux_map(result_native, 30, 30, false, ignore_direct);

    sim_optix.seed = seed;
    sim_optix.sun_gen_type = SolTrace::Data::GenType::HALTON;
    sim_optix.setup_simData();
    sim_optix.update_simulation_geometry(sim_optix.solar_azimuth, sim_optix.solar_elevation);
    SimulationResult result_optix;
    sim_optix.simulate(&result_optix, N_rays);
    sim_optix.calculate_ray_counts(result_optix);
    sim_optix.calculate_sun_size(result_optix);
    sim_optix.read_expected_all_results("1a", "N");
    sim_optix.calculate_receiver_flux_map(result_optix, 30, 30, false, ignore_direct);

    std::map<std::string, double> dict_native;
    std::map<std::string, double> dict_optix;
    if (save_flux)
    {
        std::string file_fluxmap_native = "native_flux_" + file_label + SolTrace::Data::GenTypeMap.at(sim_native.sun_gen_type)
            + "_" + std::to_string(int(N_rays / 1e3)) + "k.csv";
        sim_native.save_flux_map_to_file(file_fluxmap_native);

        std::string file_fluxmap_optix = "optix_flux_" + file_label + SolTrace::Data::GenTypeMap.at(sim_optix.sun_gen_type)
            + "_" + std::to_string(int(N_rays / 1e3)) + "k.csv";
        sim_optix.save_flux_map_to_file(file_fluxmap_optix);
    }

    // Compare
    EXPECT_EQ(result_native.get_number_of_records(), result_optix.get_number_of_records());
    EXPECT_NEAR(sim_native.helio_hit_count, sim_optix.helio_hit_count, err_abs);
    EXPECT_NEAR(sim_native.reflect_count, sim_optix.reflect_count, err_abs);
    EXPECT_NEAR(sim_native.helio_absorb_count, sim_optix.helio_absorb_count, err_abs);
    EXPECT_NEAR(sim_native.rec_absorb_count, sim_optix.rec_absorb_count, err_abs);
    EXPECT_NEAR(sim_native.rec_hit_count, sim_optix.rec_hit_count, err_abs);
    EXPECT_NEAR(sim_native.rec_direct_hit_count, sim_optix.rec_direct_hit_count, err_abs);
    EXPECT_NEAR(sim_native.rec_via_helio_hit_count, sim_optix.rec_via_helio_hit_count, err_abs);

    write_to_dict("00_helio_hit_count", sim_native.helio_hit_count, sim_optix.helio_hit_count, dict_native, dict_optix);
    write_to_dict("01_reflect_count", sim_native.reflect_count, sim_optix.reflect_count, dict_native, dict_optix);
    write_to_dict("02_helio_absorb_count", sim_native.helio_absorb_count, sim_optix.helio_absorb_count, dict_native, dict_optix);
    write_to_dict("03_rec_absorb_count", sim_native.rec_absorb_count, sim_optix.rec_absorb_count, dict_native, dict_optix);
    write_to_dict("04_rec_hit_count", sim_native.rec_hit_count, sim_optix.rec_hit_count, dict_native, dict_optix);
    write_to_dict("05_rec_direct_hit_count", sim_native.rec_direct_hit_count, sim_optix.rec_direct_hit_count, dict_native, dict_optix);
    write_to_dict("06_rec_via_helio_hit_count", sim_native.rec_via_helio_hit_count, sim_optix.rec_via_helio_hit_count, dict_native, dict_optix);

    // Helio rays add up
    EXPECT_EQ(sim_native.reflect_count + sim_native.helio_absorb_count, sim_native.helio_hit_count);
    EXPECT_EQ(sim_optix.reflect_count + sim_optix.helio_absorb_count, sim_optix.helio_hit_count);

    // Receiver rays add up
    EXPECT_EQ(sim_native.rec_direct_hit_count + sim_native.rec_via_helio_hit_count, sim_native.rec_hit_count);
    EXPECT_EQ(sim_optix.rec_direct_hit_count + sim_optix.rec_via_helio_hit_count, sim_optix.rec_hit_count);

    // Reflectivity
    double reflectivity_native = (double)sim_native.reflect_count / (double)sim_native.helio_hit_count;
    double reflectivity_optix = (double)sim_optix.reflect_count / (double)sim_optix.helio_hit_count;
    EXPECT_NEAR(reflectivity_native, 0.9, 0.02);
    EXPECT_NEAR(reflectivity_optix, 0.9, 0.02);

    write_to_dict("07_reflectivity", reflectivity_native, reflectivity_optix, dict_native, dict_optix);

    write_to_dict("08_sun_count", sim_native.sun_ray_count, sim_optix.sun_ray_count, dict_native, dict_optix);

    // Fraction of hits after helio
    double frac_rec_via_helio_native = (double)sim_native.rec_via_helio_hit_count / (double)sim_native.reflect_count;
    double frac_rec_via_helio_optix = (double)sim_optix.rec_via_helio_hit_count / (double)sim_optix.reflect_count;
    EXPECT_NEAR(frac_rec_via_helio_native, frac_rec_via_helio_optix, err_frac);

    write_to_dict("09_frac_rec_via_helio", frac_rec_via_helio_native, frac_rec_via_helio_optix, dict_native, dict_optix);

    // Compare power per ray
    double power_tol = (5. / (double)N_rays) * 1e3;
    EXPECT_NEAR(sim_native.power_per_ray, sim_optix.power_per_ray, power_tol);

    write_to_dict("10_power_per_ray", sim_native.power_per_ray, sim_optix.power_per_ray, dict_native, dict_optix);

    std::cerr << "Ray Count: " << sim_native.reflect_count << std::endl;

    // Total power absorbed
    double tol = 8.e-3;
    double total_power_native = (double)sim_native.rec_absorb_count * sim_native.power_per_ray * 1.e-3; // [kW]
    double total_power_optix = (double)sim_optix.rec_absorb_count * sim_optix.power_per_ray * 1.e-3; // [kW]
    EXPECT_NEAR(total_power_native, total_power_optix, tol * sim_native.expected_power);
    double total_power_diff = total_power_optix - total_power_native;
    double total_power_diff_frac = total_power_diff / sim_native.expected_power;

    write_to_dict("11_total_power", total_power_native, total_power_optix, dict_native, dict_optix);

    // Peak flux
    double peak_tol = 0.25;
    double peak_flux_native = sim_native.PeakFlux / 1.e3;
    double peak_flux_optix = sim_optix.PeakFlux / 1.e3;
    EXPECT_NEAR(peak_flux_native, peak_flux_optix, peak_tol * sim_native.expected_peak_flux);
    double peak_flux_diff = peak_flux_optix - peak_flux_native;
    double peak_flux_diff_frac = peak_flux_diff / sim_native.expected_peak_flux;

    write_to_dict("12_peak_flux", peak_flux_native, peak_flux_optix, dict_native, dict_optix);

    double centroid0_native = sim_native.Centroid[0];
    double centroid0_optix = sim_optix.Centroid[0];
    double centroid1_native = sim_native.Centroid[1];
    double centroid1_optix = sim_optix.Centroid[1];

    write_to_dict("13_centroid0", centroid0_native, centroid0_optix, dict_native, dict_optix);
    write_to_dict("14_centroid1", centroid1_native, centroid1_optix, dict_native, dict_optix);

    write_to_dict("15_sigmaflux", sim_native.SigmaFlux, sim_optix.SigmaFlux, dict_native, dict_optix);

    // RMS
    sim_native.calculate_receiver_flux_map(result_native, 100, 150, false, ignore_direct);  // Re-calculate for low-accuracy runs
    sim_optix.calculate_receiver_flux_map(result_optix, 100, 150, false, ignore_direct);
    EXPECT_EQ(sim_native.fluxGrid.nrows(), sim_optix.fluxGrid.nrows());
    EXPECT_EQ(sim_native.fluxGrid.ncols(), sim_optix.fluxGrid.ncols());
    double rmse = 0.0;
    for (size_t r = 0; r < sim_native.fluxGrid.nrows(); r++) {
        for (size_t c = 0; c < sim_native.fluxGrid.ncols(); c++) {
            double flux_native = sim_native.fluxGrid.at(r, c) * sim_native.zScale / 1.e3;
            double flux_optix = sim_optix.fluxGrid.at(r, c) * sim_optix.zScale / 1.e3;
            rmse += pow(flux_native - flux_optix, 2);
        }
    }

    // RMS
    rmse = sqrt(rmse / (sim_native.fluxGrid.nrows() * sim_native.fluxGrid.ncols()));
    double rmse_tol = 0.11; // Should be 0.11
    EXPECT_LE(rmse / (peak_flux_native), rmse_tol);


    // Average flux
    EXPECT_NEAR(sim_native.AveFlux/1000.0, sim_optix.AveFlux/1000.0, rmse_tol);

    write_to_dict("16_average_flux", sim_native.AveFlux/1000.0, sim_optix.AveFlux/1000.0, dict_native, dict_optix);

    // Uniformity
    write_to_dict("17_uniformity", sim_native.Uniformity, sim_optix.Uniformity, dict_native, dict_optix);

    write_to_dict("18_rmse", rmse, rmse, dict_native, dict_optix);
    double rmse_over_peak = rmse / (peak_flux_native);
    write_to_dict("19_rmse_over_peak", rmse_over_peak, rmse_over_peak, dict_native, dict_optix);

    if (save)
    {
        std::string file_outputs_native = "native_outputs_" + file_label + SolTrace::Data::GenTypeMap.at(sim_native.sun_gen_type)
            + "_" + std::to_string(int(N_rays / 1e3)) + "k.csv";
        sim_native.save_outputs(file_outputs_native, dict_native);

        std::string file_outputs_optix = "optix_outputs_" + file_label + SolTrace::Data::GenTypeMap.at(sim_optix.sun_gen_type)
            + "_" + std::to_string(int(N_rays / 1e3)) + "k.csv";
        sim_optix.save_outputs(file_outputs_optix, dict_optix);
    }
    

    int x = 0;
}


TEST(SingleHeliostatOptixNative, SingleFacetFlat_North)
{
    // Make native
    SingleHeliostatSimulationHelper<NativeRunner> sim_native;
    sim_native.runner.disable_stages();
    sim_native.runner.set_number_of_threads(N_threads);
    sim_native.initialize();

    // Make optix
    SingleHeliostatSimulationHelper<OptixRunner> sim_optix;
    sim_optix.initialize();

    CompareRunners(sim_native, sim_optix, N_rays_glob, "_1a_N_");
}

TEST(SingleHeliostatOptixNative, SingleFacetFlat_Southeast)
{
    // Make native
    SingleHeliostatSimulationHelper<NativeRunner> sim_native;
    sim_native.runner.disable_stages();
    sim_native.runner.set_number_of_threads(N_threads);
    sim_native.initialize();
    sim_native.set_heliostat_to_southeast();

    // Make optix
    SingleHeliostatSimulationHelper<OptixRunner> sim_optix;
    sim_optix.initialize();
    sim_optix.set_heliostat_to_southeast();

    CompareRunners(sim_native, sim_optix, N_rays_glob, "_1a_SE_");
}

TEST(SingleHeliostatOptixNative, SingleFacetFocused_North)
{
    // Make native
    SingleHeliostatSimulationHelper<NativeRunner> sim_native;
    sim_native.runner.disable_stages();
    sim_native.runner.set_number_of_threads(N_threads);
    sim_native.initialize();
    sim_native.set_slant_focal_length();

    // Make optix
    SingleHeliostatSimulationHelper<OptixRunner> sim_optix;
    sim_optix.initialize();
    sim_optix.set_slant_focal_length();

    CompareRunners(sim_native, sim_optix, N_rays_glob, "_1b_N_");
}

TEST(SingleHeliostatOptixNative, SingleFacetFocused_Southeast)
{
    // Make native
    SingleHeliostatSimulationHelper<NativeRunner> sim_native;
    sim_native.runner.disable_stages();
    sim_native.runner.set_number_of_threads(N_threads);
    sim_native.initialize();
    sim_native.set_heliostat_to_southeast();
    sim_native.set_slant_focal_length();

    // Make optix
    SingleHeliostatSimulationHelper<OptixRunner> sim_optix;
    sim_optix.initialize();
    sim_optix.set_heliostat_to_southeast();
    sim_optix.set_slant_focal_length();

    CompareRunners(sim_native, sim_optix, N_rays_glob, "_1b_SE_");
}

TEST(SingleHeliostatOptixNative, MultiFacetFlat_NoCanting_North)
{
    // Make native
    SingleHeliostatSimulationHelper<NativeRunner> sim_native;
    sim_native.runner.disable_stages();
    sim_native.runner.set_number_of_threads(N_threads);
    sim_native.initialize();
    sim_native.set_flat_multi_facet();

    // Make optix
    SingleHeliostatSimulationHelper<OptixRunner> sim_optix;
    sim_optix.initialize();
    sim_optix.set_flat_multi_facet();

    CompareRunners(sim_native, sim_optix, N_rays_glob, "_2_N_");
}

TEST(SingleHeliostatOptixNative, MultiFacetFlat_NoCanting_Southeast)
{
    // Make native
    SingleHeliostatSimulationHelper<NativeRunner> sim_native;
    sim_native.runner.disable_stages();
    sim_native.runner.set_number_of_threads(N_threads);
    sim_native.use_optical_errors = true;
    sim_native.use_sunshape_errors = false;
    sim_native.initialize();
    sim_native.set_flat_multi_facet();
    sim_native.set_heliostat_to_southeast();

    // Make optix
    SingleHeliostatSimulationHelper<OptixRunner> sim_optix;
    sim_optix.use_optical_errors = true;
    sim_optix.use_sunshape_errors = false;
    sim_optix.initialize();
    sim_optix.set_flat_multi_facet();
    sim_optix.set_heliostat_to_southeast();

    CompareRunners(sim_native, sim_optix, N_rays_glob, "_2_SE_");
}

TEST(SingleHeliostatOptixNative, MultiFacetFlat_SlantCanting_North)
{
    // Make native
    SingleHeliostatSimulationHelper<NativeRunner> sim_native;
    sim_native.runner.disable_stages();
    sim_native.runner.set_number_of_threads(N_threads);
    sim_native.initialize();
    sim_native.set_onaxis_slant_canting();

    // Make optix
    SingleHeliostatSimulationHelper<OptixRunner> sim_optix;
    sim_optix.initialize();
    sim_optix.set_onaxis_slant_canting();

    CompareRunners(sim_native, sim_optix, N_rays_glob, "_3_N_");
}

TEST(SingleHeliostatOptixNative, MultiFacetFlat_SlantCanting_Southeast)
{
    // Make native
    SingleHeliostatSimulationHelper<NativeRunner> sim_native;
    sim_native.runner.disable_stages();
    sim_native.runner.set_number_of_threads(N_threads);
    sim_native.initialize();
    sim_native.set_heliostat_to_southeast();
    sim_native.set_onaxis_slant_canting();

    // Make optix
    SingleHeliostatSimulationHelper<OptixRunner> sim_optix;
    sim_optix.initialize();
    sim_optix.set_heliostat_to_southeast();
    sim_optix.set_onaxis_slant_canting();

    CompareRunners(sim_native, sim_optix, N_rays_glob, "_3_SE_");
}

TEST(SingleHeliostatOptixNative, MultiFacetFocused_SlantCanting_North)
{
    // Make native
    SingleHeliostatSimulationHelper<NativeRunner> sim_native;
    sim_native.runner.disable_stages();
    sim_native.runner.set_number_of_threads(N_threads);
    sim_native.initialize();
    sim_native.set_slant_focal_length();
    sim_native.set_onaxis_slant_canting();

    // Make optix
    SingleHeliostatSimulationHelper<OptixRunner> sim_optix;
    sim_optix.initialize();
    sim_optix.set_slant_focal_length();
    sim_optix.set_onaxis_slant_canting();

    CompareRunners(sim_native, sim_optix, N_rays_glob, "_4_N_");
}

TEST(SingleHeliostatOptixNative, MultiFacetFocused_SlantCanting_Southeast)
{
    // Make native
    SingleHeliostatSimulationHelper<NativeRunner> sim_native;
    sim_native.runner.disable_stages();
    sim_native.runner.set_number_of_threads(N_threads);
    sim_native.initialize();
    sim_native.set_heliostat_to_southeast();
    sim_native.set_slant_focal_length();
    sim_native.set_onaxis_slant_canting();

    // Make optix
    SingleHeliostatSimulationHelper<OptixRunner> sim_optix;
    sim_optix.initialize();
    sim_optix.set_heliostat_to_southeast();
    sim_optix.set_slant_focal_length();
    sim_optix.set_onaxis_slant_canting();

    CompareRunners(sim_native, sim_optix, N_rays_glob, "_4_SE_");
}

TEST(SingleHeliostatOptixNative_ErrorTesting, NoErrors)
{
    bool use_optical = false;
    bool use_sunshape = false;

    // Make native
    SingleHeliostatSimulationHelper<NativeRunner> sim_native;
    sim_native.runner.disable_stages(); // Disable stages
    sim_native.use_optical_errors = use_optical;
    sim_native.use_sunshape_errors = use_sunshape;
    sim_native.runner.set_number_of_threads(N_threads);
    sim_native.initialize();
    
    // Make optix
    SingleHeliostatSimulationHelper<OptixRunner> sim_optix;
    sim_optix.use_optical_errors = use_optical;
    sim_optix.use_sunshape_errors = use_sunshape;
    sim_optix.initialize();

    CompareRunners(sim_native, sim_optix, N_rays_glob);
}

TEST(SingleHeliostatOptixNative_ErrorTesting, SunShapePillBoxOnly)
{
    bool use_optical = false;
    bool use_sunshape = true;

    // Make native
    SingleHeliostatSimulationHelper<NativeRunner> sim_native;
    sim_native.runner.disable_stages(); // Disable stages
    sim_native.use_optical_errors = use_optical;
    sim_native.use_sunshape_errors = use_sunshape;
    sim_native.runner.set_number_of_threads(N_threads);
    sim_native.initialize();

    // Make optix
    SingleHeliostatSimulationHelper<OptixRunner> sim_optix;
    sim_optix.use_optical_errors = use_optical;
    sim_optix.use_sunshape_errors = use_sunshape;
    sim_optix.initialize();

    CompareRunners(sim_native, sim_optix, N_rays_glob);
}

TEST(SingleHeliostatOptixNative_ErrorTesting, SunShapeGaussOnly)
{
    bool use_optical = false;
    bool use_sunshape = true;

    // Make native
    SingleHeliostatSimulationHelper<NativeRunner> sim_native;
    sim_native.runner.disable_stages(); // Disable stages
    sim_native.use_optical_errors = use_optical;
    sim_native.use_sunshape_errors = use_sunshape;
    sim_native.sun_shape = SunShape::GAUSSIAN;
    sim_native.gauss_sigma = 2;
    sim_native.runner.set_number_of_threads(N_threads);
    sim_native.initialize();

    // Make optix
    SingleHeliostatSimulationHelper<OptixRunner> sim_optix;
    sim_optix.use_optical_errors = use_optical;
    sim_optix.use_sunshape_errors = use_sunshape;
    sim_optix.sun_shape = SunShape::GAUSSIAN;
    sim_optix.gauss_sigma = 2;
    sim_optix.initialize();

    CompareRunners(sim_native, sim_optix, N_rays_glob);
}

TEST(SingleHeliostatOptixNative_ErrorTesting, SlopeGaussOnly)
{
    bool use_optical = true;
    bool use_sunshape = false;

    // Make native
    SingleHeliostatSimulationHelper<NativeRunner> sim_native;
    sim_native.runner.disable_stages(); // Disable stages
    sim_native.use_optical_errors = use_optical;
    sim_native.use_sunshape_errors = use_sunshape;
    sim_native.runner.set_number_of_threads(N_threads);
    sim_native.initialize();

    // Make optix
    SingleHeliostatSimulationHelper<OptixRunner> sim_optix;
    sim_optix.use_optical_errors = use_optical;
    sim_optix.use_sunshape_errors = use_sunshape;
    sim_optix.initialize();

    CompareRunners(sim_native, sim_optix, N_rays_glob);
}

TEST(SingleHeliostatOptixNative_ErrorTesting, SlopePillBoxOnly)
{
    bool use_optical = true;
    bool use_sunshape = false;

    // Make native
    SingleHeliostatSimulationHelper<NativeRunner> sim_native;
    sim_native.runner.disable_stages(); // Disable stages
    sim_native.use_optical_errors = use_optical;
    sim_native.use_sunshape_errors = use_sunshape;
    sim_native.error_dist = DistributionType::PILLBOX;
    sim_native.runner.set_number_of_threads(N_threads);
    sim_native.initialize();

    // Make optix
    SingleHeliostatSimulationHelper<OptixRunner> sim_optix;
    sim_optix.use_optical_errors = use_optical;
    sim_optix.use_sunshape_errors = use_sunshape;
    sim_optix.error_dist = DistributionType::PILLBOX;
    sim_optix.initialize();

    CompareRunners(sim_native, sim_optix, N_rays_glob);
}

TEST(SingleHeliostatOptixNative_ErrorTesting, SunShapeAndSlope)
{
    bool use_optical = true;
    bool use_sunshape = true;

    // Make native
    SingleHeliostatSimulationHelper<NativeRunner> sim_native;
    sim_native.runner.disable_stages(); // Disable stages
    sim_native.use_optical_errors = use_optical;
    sim_native.use_sunshape_errors = use_sunshape;
    sim_native.runner.set_number_of_threads(N_threads);
    sim_native.initialize();

    // Make optix
    SingleHeliostatSimulationHelper<OptixRunner> sim_optix;
    sim_optix.use_optical_errors = use_optical;
    sim_optix.use_sunshape_errors = use_sunshape;
    sim_optix.initialize();

    CompareRunners(sim_native, sim_optix, N_rays_glob);
}

TEST(SingleHeliostatOptixNative_ErrorTesting, SpecGaussOnly)
{
    bool use_optical = true;
    bool use_sunshape = false;

    // Make native
    SingleHeliostatSimulationHelper<NativeRunner> sim_native;
    sim_native.runner.disable_stages(); // Disable stages
    sim_native.use_optical_errors = use_optical;
    sim_native.use_sunshape_errors = use_sunshape;
    sim_native.slope_error = 0;     // Turn off slope error
    sim_native.spec_error = 2;    // Set specularity error
    sim_native.runner.set_number_of_threads(N_threads);
    sim_native.initialize();

    // Make optix
    SingleHeliostatSimulationHelper<OptixRunner> sim_optix;
    sim_optix.use_optical_errors = use_optical;
    sim_optix.use_sunshape_errors = use_sunshape;
    sim_optix.slope_error = 0;     // Turn off slope error
    sim_optix.spec_error = 2;    // Set specularity error
    sim_optix.initialize();

    CompareRunners(sim_native, sim_optix, N_rays_glob);
}

TEST(SingleHeliostatOptixNative_ErrorTesting, SpecPillBoxOnly)
{
    bool use_optical = true;
    bool use_sunshape = false;

    // Make native
    SingleHeliostatSimulationHelper<NativeRunner> sim_native;
    sim_native.runner.disable_stages(); // Disable stages
    sim_native.use_optical_errors = use_optical;
    sim_native.use_sunshape_errors = use_sunshape;
    sim_native.slope_error = 0;     // Turn off slope error
    sim_native.spec_error = 2;      // Set specularity error
    sim_native.error_dist = DistributionType::PILLBOX;
    sim_native.runner.set_number_of_threads(N_threads);
    sim_native.initialize();

    // Make optix
    SingleHeliostatSimulationHelper<OptixRunner> sim_optix;
    sim_optix.use_optical_errors = use_optical;
    sim_optix.use_sunshape_errors = use_sunshape;
    sim_optix.slope_error = 0;      // Turn off slope error
    sim_optix.spec_error = 2;       // Set specularity error
    sim_optix.error_dist = DistributionType::PILLBOX;
    sim_optix.initialize();

    CompareRunners(sim_native, sim_optix, N_rays_glob);
}

TEST(SingleHeliostatOptixNative_ErrorTesting, SunShapeSlopeAndSpec)
{
    bool use_optical = true;
    bool use_sunshape = true;

    // Make native
    SingleHeliostatSimulationHelper<NativeRunner> sim_native;
    sim_native.runner.disable_stages(); // Disable stages
    sim_native.use_optical_errors = use_optical;
    sim_native.use_sunshape_errors = use_sunshape;
    sim_native.slope_error = 2;     // Turn off slope error
    sim_native.spec_error = 2;    // Set specularity error
    sim_native.runner.set_number_of_threads(N_threads);
    sim_native.initialize();

    // Make optix
    SingleHeliostatSimulationHelper<OptixRunner> sim_optix;
    sim_optix.use_optical_errors = use_optical;
    sim_optix.use_sunshape_errors = use_sunshape;
    sim_optix.slope_error = 2;     // Turn off slope error
    sim_optix.spec_error = 2;    // Set specularity error
    sim_optix.initialize();

    CompareRunners(sim_native, sim_optix, N_rays_glob);
}

TEST(SingleHeliostatOptixNative_ErrorTesting, SunShapeLimbDarkenedOnly)
{
    bool use_optical = false;
    bool use_sunshape = true;

    // Make native
    SingleHeliostatSimulationHelper<NativeRunner> sim_native;
    sim_native.runner.disable_stages(); // Disable stages
    sim_native.use_optical_errors = use_optical;
    sim_native.use_sunshape_errors = use_sunshape;
    sim_native.runner.set_number_of_threads(N_threads);
    sim_native.sun_shape = SunShape::LIMBDARKENED;
    sim_native.initialize();

    // Make optix
    SingleHeliostatSimulationHelper<OptixRunner> sim_optix;
    sim_optix.use_optical_errors = use_optical;
    sim_optix.use_sunshape_errors = use_sunshape;
    sim_optix.sun_shape = SunShape::LIMBDARKENED;
    sim_optix.initialize();

    CompareRunners(sim_native, sim_optix, N_rays_glob);
}

TEST(SingleHeliostatOptixNative_ErrorTesting, SunShapeUserDefinedOnly)
{
    const std::vector<double> user_angle = {
        0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.05, 1.2, 1.35,
        1.5, 1.65, 1.8, 1.95, 2.1, 2.25, 2.4, 2.55, 2.7, 2.85,
        3.0, 3.15, 3.3, 3.45, 3.6, 3.75, 3.9, 4.05, 4.2, 4.35,
        4.5, 4.65, 4.8, 4.95, 5.1, 5.25, 5.4, 5.55, 5.7, 5.85,
        6.0, 6.15, 6.3, 6.45, 6.6, 6.75, 6.9, 7.05, 7.2, 7.35,
        7.5, 7.65, 7.8, 7.95
    };

    const std::vector<double> user_intensity = {
        1.0, 0.999872, 0.999485, 0.998837, 0.997923, 0.996734, 0.99526, 0.993487, 0.991399, 0.988976,
        0.986193, 0.983019, 0.979417, 0.975345, 0.970747, 0.965558, 0.959697, 0.953063, 0.945528, 0.936933,
        0.927069, 0.915665, 0.902358, 0.886653, 0.867855, 0.844965, 0.816477, 0.78003, 0.731687, 0.66436,
        0.563875, 0.397159, 5.34e-05, 5.07e-05, 4.82e-05, 4.59e-05, 4.38e-05, 4.18e-05, 3.99e-05, 3.82e-05,
        3.66e-05, 3.51e-05, 3.37e-05, 3.24e-05, 3.11e-05, 3.00e-05, 2.89e-05, 2.78e-05, 2.69e-05, 2.59e-05,
        2.51e-05, 2.42e-05, 2.34e-05, 2.27e-05
    };

    bool use_optical = false;
    bool use_sunshape = true;

    // Make native
    SingleHeliostatSimulationHelper<NativeRunner> sim_native;
    sim_native.runner.disable_stages(); // Disable stages
    sim_native.use_optical_errors = use_optical;
    sim_native.use_sunshape_errors = use_sunshape;
    sim_native.runner.set_number_of_threads(N_threads);
    sim_native.sun_shape = SunShape::USER_DEFINED;
    sim_native.user_angle = user_angle;
    sim_native.user_intensity = user_intensity;
    sim_native.initialize();

    // Make optix
    SingleHeliostatSimulationHelper<OptixRunner> sim_optix;
    sim_optix.use_optical_errors = use_optical;
    sim_optix.use_sunshape_errors = use_sunshape;
    sim_optix.sun_shape = SunShape::USER_DEFINED;
    sim_optix.user_angle = user_angle;
    sim_optix.user_intensity = user_intensity;
    sim_optix.initialize();

    CompareRunners(sim_native, sim_optix, N_rays_glob);
}
#include "isolated_heliostat_test_template.hpp"

#include <embree_runner.hpp>
#include <optix_runner.hpp>
#include <native_runner.hpp>

using SolTrace::NativeRunner::NativeRunner;
using OptixRunnerType = OptixRunner;
using NativeRunnerType = NativeRunner;

using SolTrace::EmbreeRunner::EmbreeRunner;

constexpr uint_fast64_t N_rays_glob = 2e6;
constexpr int seed = 123;
constexpr bool save = true;
constexpr bool save_hits = false;
constexpr bool ignore_direct = true;

static const int N_threads = static_cast<int>(std::max(1u, std::min(std::thread::hardware_concurrency(), 10u)));

static void write_to_dict(std::string key_name, double val_a,
	double val_b, std::map<std::string, double>& dict_a,
	std::map<std::string, double>& dict_b)
{
	dict_a[key_name] = val_a;
	dict_b[key_name] = val_b;
}

static void CompareRunners(IsolatedHeliostatSimulationHelper<EmbreeRunner>& sim_embree,
	IsolatedHeliostatSimulationHelper<OptixRunner>& sim_optix, uint_fast64_t N_rays,
	const std::string& hour, const std::string& file_label = "")
{
	sim_embree.update_from_hour(hour);
	SimulationResult result_embree;
	sim_embree.simulate(&result_embree, N_rays);
	sim_embree.calculate_ray_counts(result_embree);
	sim_embree.calculate_sun_size(result_embree);
	sim_embree.calculate_outputs(result_embree, ignore_direct);

	if (save_hits)
	{
		std::string file_hits_embree = "embree_hits_" + file_label + std::to_string(int(N_rays / 1e3)) + "k.csv";
		save_hit_pos_to_file(result_embree, file_hits_embree);
	}

	sim_optix.update_from_hour(hour);
	SimulationResult result_optix;
	sim_optix.simulate(&result_optix, N_rays);
	sim_optix.calculate_ray_counts(result_optix);
	sim_optix.calculate_sun_size(result_optix);
	sim_optix.calculate_outputs(result_optix, ignore_direct);

	if (save_hits)
	{
		std::string file_hits_optix = "optix_hits_" + file_label + std::to_string(int(N_rays / 1e3)) + "k.csv";
		save_hit_pos_to_file(result_optix, file_hits_optix);
	}
	
	// Error tolerances
	double err_frac = 0.01;  //0.01;
	double err_abs = err_frac * (double)N_rays;

	std::map<std::string, double> dict_embree;
	std::map<std::string, double> dict_optix;

	// Compare hits
	ASSERT_EQ(result_embree.get_number_of_records(), result_optix.get_number_of_records());
	EXPECT_NEAR(sim_embree.tot_helio_hits, sim_optix.tot_helio_hits, err_abs);
	EXPECT_NEAR(sim_embree.tot_helio_absorb_count, sim_optix.tot_helio_absorb_count, err_abs);
	EXPECT_NEAR(sim_embree.tot_reflect_count, sim_optix.tot_reflect_count, err_abs);
	EXPECT_NEAR(sim_embree.rec_absorb_count, sim_optix.rec_absorb_count, err_abs);
	EXPECT_NEAR(sim_embree.tot_helio_block_count, sim_optix.tot_helio_block_count, err_abs);
	//EXPECT_NEAR(sim_embree.heat_shield_absorb_count, sim_optix.heat_shield_absorb_count, err_abs);
	
	EXPECT_NEAR(sim_embree.rec_direct_count, sim_optix.rec_direct_count, err_abs);
	EXPECT_NEAR(sim_embree.rec_via_helio_count, sim_optix.rec_via_helio_count, err_abs);
	
	// Do not check tot_helio_hits or rec_via_rec_count
	// Rays that get 'stuck' inside the cylinder hit the max depth for optix
	// so the receiver registers more hits (tot_rec_hits) but this is because
	// Optix does not fully track the 'trapped' rays inside the receiver
	//EXPECT_NEAR(sim_embree.tot_rec_hits, sim_optix.tot_rec_hits, err_abs);
	//EXPECT_NEAR(sim_embree.rec_via_rec_count, sim_optix.rec_via_rec_count, err_abs);

	write_to_dict("00_tot_helio_hits", sim_embree.tot_helio_hits, sim_optix.tot_helio_hits, dict_embree, dict_optix);
	write_to_dict("01_tot_helio_absorb_count", sim_embree.tot_helio_absorb_count, sim_optix.tot_helio_absorb_count, dict_embree, dict_optix);
	write_to_dict("02_tot_reflect_count", sim_embree.tot_reflect_count, sim_optix.tot_reflect_count, dict_embree, dict_optix);
	write_to_dict("03_rec_absorb_count", sim_embree.rec_absorb_count, sim_optix.rec_absorb_count, dict_embree, dict_optix);
	write_to_dict("04_tot_helio_block_count", sim_embree.tot_helio_block_count, sim_optix.tot_helio_block_count, dict_embree, dict_optix);
	//write_to_dict("05_heat_shield_absorb_count", sim_embree.heat_shield_absorb_count, sim_optix.heat_shield_absorb_count, dict_embree, dict_optix);
	write_to_dict("06_tot_rec_hits", sim_embree.tot_rec_hits, sim_optix.tot_rec_hits, dict_embree, dict_optix);
	write_to_dict("07_rec_direct_count", sim_embree.rec_direct_count, sim_optix.rec_direct_count, dict_embree, dict_optix);
	write_to_dict("08_rec_via_helio_count", sim_embree.rec_via_helio_count, sim_optix.rec_via_helio_count, dict_embree, dict_optix);

	// Helio hits add up
	EXPECT_EQ(sim_embree.tot_helio_hits, sim_embree.tot_helio_absorb_count + sim_embree.tot_reflect_count);
	EXPECT_EQ(sim_optix.tot_helio_hits, sim_optix.tot_helio_absorb_count + sim_optix.tot_reflect_count);

	// Receiver hits add up
	EXPECT_EQ(sim_embree.tot_rec_hits, sim_embree.rec_direct_count + sim_embree.rec_via_helio_count + sim_embree.rec_via_rec_count);
	EXPECT_EQ(sim_optix.tot_rec_hits, sim_optix.rec_direct_count + sim_optix.rec_via_helio_count + sim_optix.rec_via_rec_count);

	double refl_embree = (double)sim_embree.tot_reflect_count / (double)sim_embree.tot_helio_hits;
	double refl_optix = (double)sim_optix.tot_reflect_count / (double)sim_optix.tot_helio_hits;

	write_to_dict("09_reflectivity", refl_embree, refl_optix, dict_embree, dict_optix);

	// Sun Count
	write_to_dict("10_sun_count", sim_embree.nsun_rays, sim_optix.nsun_rays, dict_embree, dict_optix);

	// Fraction reflected hits that hit receiver
	double frac_via_helio_a = ((double)sim_embree.rec_absorb_count - (double)sim_embree.rec_direct_count) / (double)sim_embree.tot_reflect_count;
	double frac_via_helio_b = ((double)sim_optix.rec_absorb_count - (double)sim_optix.rec_direct_count) / (double)sim_optix.tot_reflect_count;
	EXPECT_NEAR(frac_via_helio_a, frac_via_helio_b, err_frac);

	write_to_dict("11_frac_via_helio", frac_via_helio_a, frac_via_helio_b, dict_embree, dict_optix);

	// Compare power per ray
	write_to_dict("12_power_per_ray", sim_embree.power_per_ray, sim_optix.power_per_ray, dict_embree, dict_optix);

	// Total power absorbed
	double tol = 8.e-3;
	EXPECT_NEAR(sim_embree.total_power, sim_optix.total_power, tol * sim_embree.total_power);

	write_to_dict("13_total_power", sim_embree.total_power, sim_optix.total_power, dict_embree, dict_optix);

	// Peak flux
	double peak_tol = 0.25;
	double peak_flux_embree = sim_embree.PeakFlux / 1.e3;
	double peak_flux_optix = sim_optix.PeakFlux / 1.e3;
	EXPECT_NEAR(peak_flux_embree, peak_flux_optix, peak_tol * peak_flux_embree);

	write_to_dict("14_peak_flux", peak_flux_embree, peak_flux_optix, dict_embree, dict_optix);

	// RMS of flux values
	EXPECT_EQ(sim_embree.fluxGrid.nrows(), sim_optix.fluxGrid.nrows());
	EXPECT_EQ(sim_embree.fluxGrid.ncols(), sim_optix.fluxGrid.ncols());
	double rmse = 0.0;
	for (size_t r = 0; r < sim_embree.fluxGrid.nrows(); r++) {
		for (size_t c = 0; c < sim_embree.fluxGrid.ncols(); c++) {
			double flux_embree = sim_embree.fluxGrid.at(r, c) * sim_embree.zScale / 1.e3;
			double flux_optix = sim_optix.fluxGrid.at(r, c) * sim_optix.zScale / 1.e3;
			rmse += pow(flux_embree - flux_optix, 2);
		}
	}

	rmse = sqrt(rmse / (sim_embree.fluxGrid.nrows() * sim_embree.fluxGrid.ncols()));
	double rmse_tol = 0.11;  
	EXPECT_LE(rmse / peak_flux_embree, rmse_tol);

	// Average flux
	//EXPECT_NEAR(sim_embree.AveFlux / 1000.0, sim_optix.AveFlux / 1000.0, rmse_tol);

	write_to_dict("15_average_flux", sim_embree.AveFlux / 1000.0, sim_optix.AveFlux / 1000.0, dict_embree, dict_optix);

	// Binning
	//EXPECT_EQ(sim_embree.NotBinned, sim_optix.NotBinned);

	write_to_dict("16_not_binned", sim_embree.NotBinned, sim_optix.NotBinned, dict_embree, dict_optix);
	write_to_dict("17_neg_x_bin_err", sim_embree.max_neg_x_flux_err, sim_optix.max_neg_x_flux_err, dict_embree, dict_optix);
	write_to_dict("18_pos_x_bin_err", sim_embree.max_pos_x_flux_err, sim_optix.max_pos_x_flux_err, dict_embree, dict_optix);

	write_to_dict("20_sigmaflux", sim_embree.SigmaFlux, sim_optix.SigmaFlux, dict_embree, dict_optix);
	write_to_dict("21_uniformity", sim_embree.Uniformity, sim_optix.Uniformity, dict_embree, dict_optix);
	write_to_dict("22_rmse", rmse, rmse, dict_embree, dict_optix);
	double rmse_over_peak = rmse / (peak_flux_embree);
	write_to_dict("23_rmse_over_peak", rmse_over_peak, rmse_over_peak, dict_embree, dict_optix);

	write_to_dict("24_rec_via_rec_count", sim_embree.rec_via_rec_count, sim_optix.rec_via_rec_count, dict_embree, dict_optix);

	if (save)
	{
		std::string file_fluxmap_native = "embree_field_flux_" + file_label + SolTrace::Data::GenTypeMap.at(sim_embree.sun_gen_type)
			+ "_" + std::to_string(int(N_rays / 1e3)) + "k.csv";
		sim_embree.save_flux_map_to_file(file_fluxmap_native);

		std::string file_fluxmap_optix = "optix_field_flux_" + file_label + SolTrace::Data::GenTypeMap.at(sim_optix.sun_gen_type)
			+ "_" + std::to_string(int(N_rays / 1e3)) + "k.csv";
		sim_optix.save_flux_map_to_file(file_fluxmap_optix);

		std::string file_outputs_native = "embree_outputs_" + file_label + SolTrace::Data::GenTypeMap.at(sim_embree.sun_gen_type)
			+ "_" + std::to_string(int(N_rays / 1e3)) + "k.csv";
		sim_embree.save_outputs(file_outputs_native, dict_embree);

		std::string file_outputs_optix = "optix_outputs_" + file_label + SolTrace::Data::GenTypeMap.at(sim_optix.sun_gen_type)
			+ "_" + std::to_string(int(N_rays / 1e3)) + "k.csv";
		sim_optix.save_outputs(file_outputs_optix, dict_optix);
	}
}

/*
DEFAULT IS FACET FOCUS TO SLANT RANGE
Tests were based off of Phase III of the round robin paper. 
Some stinput files provided had differences from the paper, and some tests needed further adjustments to match the result fluxmaps.
Edits are noted by each test.
*/

//task 4a:
//stinput differences: receiver origin at {0,0,171.035}, flat facets
TEST(IsolatedHeliostatOptixEmbree, multiFacet1332_BlockingShading4a)
{
	// Make embree
	IsolatedHeliostatSimulationHelper<EmbreeRunner> sim_embree;
	sim_embree.runner.disable_stages();
	sim_embree.runner.set_number_of_threads(N_threads);
	sim_embree.initialize();
	sim_embree.seed = seed;
	sim_embree.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Make optix
	IsolatedHeliostatSimulationHelper<OptixRunner> sim_optix;
	sim_optix.initialize();
	sim_optix.seed = seed;
	sim_optix.sun_gen_type = SolTrace::Data::GenType::HALTON;

    sim_embree.set_flat_facets();
    sim_optix.set_flat_facets();

    glm::dvec3 origin = {0,0,171.035};
    sim_embree.set_rec_origin(origin);
    sim_optix.set_rec_origin(origin);

    std::vector<int> active {1332};
    std::vector<int> blocking {1266, 1304, 1306};
    sim_embree.create_active_heliostats(active);
    sim_embree.create_blocking_heliostats(blocking);
    sim_optix.create_active_heliostats(active);
    sim_optix.create_blocking_heliostats(blocking);

    sim_embree.assign_canted_banded(false);
    sim_embree.assign_canted_banded(true);
    sim_optix.assign_canted_banded(false);
    sim_optix.assign_canted_banded(true);

    sim_embree.setup_simData();
    sim_optix.setup_simData();

    CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "4a_1_12_");
}

//task 4b: 
//stinput differences: receiver origin {0,0,171.035}
//edits: no slope error
TEST(IsolatedHeliostatOptixEmbree, multiFacet8993_BlockingShading4b)
{
	// Make embree
	IsolatedHeliostatSimulationHelper<EmbreeRunner> sim_embree;
	sim_embree.runner.disable_stages();
	sim_embree.runner.set_number_of_threads(N_threads);
	sim_embree.initialize();
	sim_embree.seed = seed;
	sim_embree.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Make optix
	IsolatedHeliostatSimulationHelper<OptixRunner> sim_optix;
	sim_optix.initialize();
	sim_optix.seed = seed;
	sim_optix.sun_gen_type = SolTrace::Data::GenType::HALTON;

    sim_embree.set_slope_error(0.0);
    sim_optix.set_slope_error(0.0);

    glm::dvec3 origin = {0,0,171.035};
    sim_embree.set_rec_origin(origin);
    sim_optix.set_rec_origin(origin);

    std::vector<int> active {8993};
    std::vector<int> blocking {9100, 9102, 9208};
    sim_embree.create_active_heliostats(active);
    sim_embree.create_blocking_heliostats(blocking);
    sim_optix.create_active_heliostats(active);
    sim_optix.create_blocking_heliostats(blocking);

    sim_embree.assign_canted_banded(false);
    sim_embree.assign_canted_banded(true);
    sim_optix.assign_canted_banded(false);
    sim_optix.assign_canted_banded(true);

    sim_embree.assign_focal_lengths_banded(false);
    sim_embree.assign_focal_lengths_banded(true);
    sim_optix.assign_focal_lengths_banded(false);
    sim_optix.assign_focal_lengths_banded(true);

    sim_embree.setup_simData();
    sim_optix.setup_simData();

    CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "4b_1_12_");
	
}

//task 4d_8: 
//edits: no slope error
TEST(IsolatedHeliostatOptixEmbree, singleFacet8993_LongerAimpoint4d_8)
{
	// Make embree
	IsolatedHeliostatSimulationHelper<EmbreeRunner> sim_embree;
	sim_embree.runner.disable_stages();
	sim_embree.runner.set_number_of_threads(N_threads);
	sim_embree.initialize();
	sim_embree.seed = seed;
	sim_embree.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Make optix
	IsolatedHeliostatSimulationHelper<OptixRunner> sim_optix;
	sim_optix.initialize();
	sim_optix.seed = seed;
	sim_optix.sun_gen_type = SolTrace::Data::GenType::HALTON;

    sim_embree.set_slope_error(0.0);
    sim_optix.set_slope_error(0.0);

    std::vector<int> active {8993};
    sim_embree.create_active_heliostats(active);
    sim_optix.create_active_heliostats(active);

    sim_embree.setup_simData();
    sim_optix.setup_simData();

    CompareRunners(sim_embree, sim_optix, N_rays_glob, "8", "4d_1_8_");

}

//task 4d_12
//edits: no slope error, no sunshape
TEST(IsolatedHeliostatOptixEmbree, singleFacet8993_LongerAimpoint4d_12)
{
	// Make embree
	IsolatedHeliostatSimulationHelper<EmbreeRunner> sim_embree;
	sim_embree.runner.disable_stages();
	sim_embree.runner.set_number_of_threads(N_threads);
	sim_embree.initialize();
	sim_embree.seed = seed;
	sim_embree.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Make optix
	IsolatedHeliostatSimulationHelper<OptixRunner> sim_optix;
	sim_optix.initialize();
	sim_optix.seed = seed;
	sim_optix.sun_gen_type = SolTrace::Data::GenType::HALTON;

    sim_embree.set_slope_error(0.0);
    sim_optix.set_slope_error(0.0);

    std::vector<int> active {8993};
    sim_embree.create_active_heliostats(active);
    sim_optix.create_active_heliostats(active);

    sim_embree.set_no_sunShape();
    sim_optix.set_no_sunShape();

    sim_embree.setup_simData();
    sim_optix.setup_simData();

    CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "4d_1_12_");
	
}

//task 4e:
TEST(IsolatedHeliostatOptixEmbree,  singleFacet8993_TargetCoordSystemE4e)
{
	// Make embree
	IsolatedHeliostatSimulationHelper<EmbreeRunner> sim_embree;
	sim_embree.runner.disable_stages();
	sim_embree.runner.set_number_of_threads(N_threads);
	sim_embree.initialize();
	sim_embree.seed = seed;
	sim_embree.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Make optix
	IsolatedHeliostatSimulationHelper<OptixRunner> sim_optix;
	sim_optix.initialize();
	sim_optix.seed = seed;
	sim_optix.sun_gen_type = SolTrace::Data::GenType::HALTON;

    glm::dvec3 shift = {-1,0,0};
    sim_embree.shift_aimpoint(shift);
    sim_optix.shift_aimpoint(shift);

    sim_embree.set_slope_error(0.0);
    sim_optix.set_slope_error(0.0);

    std::vector<int> active {8993};
    sim_embree.create_active_heliostats(active);
    sim_optix.create_active_heliostats(active);

    sim_embree.set_no_sunShape();
    sim_optix.set_no_sunShape();

    sim_embree.setup_simData();
    sim_optix.setup_simData();

    CompareRunners(sim_embree, sim_optix, N_rays_glob, "8", "4e_1_8_");
    CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "4e_1_12_");
	
}

//task 4f: 
TEST(IsolatedHeliostatOptixEmbree,  singleFacet8993_TargetCoordSystemW4f)
{
	// Make embree
	IsolatedHeliostatSimulationHelper<EmbreeRunner> sim_embree;
	sim_embree.runner.disable_stages();
	sim_embree.runner.set_number_of_threads(N_threads);
	sim_embree.initialize();
	sim_embree.seed = seed;
	sim_embree.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Make optix
	IsolatedHeliostatSimulationHelper<OptixRunner> sim_optix;
	sim_optix.initialize();
	sim_optix.seed = seed;
	sim_optix.sun_gen_type = SolTrace::Data::GenType::HALTON;

    glm::dvec3 shift = {1,0,0};
    sim_embree.shift_aimpoint(shift);
    sim_optix.shift_aimpoint(shift);

    sim_embree.set_slope_error(0.0);
    sim_optix.set_slope_error(0.0);

    std::vector<int> active {8993};
    sim_embree.create_active_heliostats(active);
    sim_optix.create_active_heliostats(active);

    sim_embree.set_no_sunShape();
    sim_optix.set_no_sunShape();

    sim_embree.setup_simData();
    sim_optix.setup_simData();

    CompareRunners(sim_embree, sim_optix, N_rays_glob, "8", "4f_1_8_");
    CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "4f_1_12_");
	
}

//task 4g:
TEST(IsolatedHeliostatOptixEmbree,  singleFacet8993_TargetCoordSystemU4g)
{
	// Make embree
	IsolatedHeliostatSimulationHelper<EmbreeRunner> sim_embree;
	sim_embree.runner.disable_stages();
	sim_embree.runner.set_number_of_threads(N_threads);
	sim_embree.initialize();
	sim_embree.seed = seed;
	sim_embree.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Make optix
	IsolatedHeliostatSimulationHelper<OptixRunner> sim_optix;
	sim_optix.initialize();
	sim_optix.seed = seed;
	sim_optix.sun_gen_type = SolTrace::Data::GenType::HALTON;

    glm::dvec3 shift = {0,0,1};
    sim_embree.shift_aimpoint(shift);
    sim_optix.shift_aimpoint(shift);

    sim_embree.set_slope_error(0.0);
    sim_optix.set_slope_error(0.0);

    std::vector<int> active {8993};
    sim_embree.create_active_heliostats(active);
    sim_optix.create_active_heliostats(active);

    sim_embree.set_no_sunShape();
    sim_optix.set_no_sunShape();

    sim_embree.setup_simData();
    sim_optix.setup_simData();

    CompareRunners(sim_embree, sim_optix, N_rays_glob, "8", "4g_1_8_");
    CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "4g_1_12_");
	
}

//task 4h:
TEST(IsolatedHeliostatOptixEmbree,  singleFacet8993_TargetCoordSystemD4h)
{
	// Make embree
	IsolatedHeliostatSimulationHelper<EmbreeRunner> sim_embree;
	sim_embree.runner.disable_stages();
	sim_embree.runner.set_number_of_threads(N_threads);
	sim_embree.initialize();
	sim_embree.seed = seed;
	sim_embree.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Make optix
	IsolatedHeliostatSimulationHelper<OptixRunner> sim_optix;
	sim_optix.initialize();
	sim_optix.seed = seed;
	sim_optix.sun_gen_type = SolTrace::Data::GenType::HALTON;

    glm::dvec3 shift = {0,0,-1};
    sim_embree.shift_aimpoint(shift);
    sim_optix.shift_aimpoint(shift);

    sim_embree.set_slope_error(0.0);
    sim_optix.set_slope_error(0.0);

    std::vector<int> active {8993};
    sim_embree.create_active_heliostats(active);
    sim_optix.create_active_heliostats(active);

    sim_embree.set_no_sunShape();
    sim_optix.set_no_sunShape();

    sim_embree.setup_simData();
    sim_optix.setup_simData();

    CompareRunners(sim_embree, sim_optix, N_rays_glob, "8", "4h_1_8_");
    CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "4h_1_12_");

}

//task 4i:
//stinput differences: flat facets 
TEST(IsolatedHeliostatOptixEmbree,  singleFacet8993_4i_8)
{
	// Make embree
	IsolatedHeliostatSimulationHelper<EmbreeRunner> sim_embree;
	sim_embree.runner.disable_stages();
	sim_embree.runner.set_number_of_threads(N_threads);
	sim_embree.initialize();
	sim_embree.seed = seed;
	sim_embree.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Make optix
	IsolatedHeliostatSimulationHelper<OptixRunner> sim_optix;
	sim_optix.initialize();
	sim_optix.seed = seed;
	sim_optix.sun_gen_type = SolTrace::Data::GenType::HALTON;

    sim_embree.set_slope_error(0.0);
    sim_optix.set_slope_error(0.0);

    sim_embree.set_flat_facets();
    sim_optix.set_flat_facets();

    std::vector<int> active {8993};
    sim_embree.create_active_heliostats(active);
    sim_optix.create_active_heliostats(active);

    sim_embree.set_no_sunShape();
    sim_optix.set_no_sunShape();

    sim_embree.setup_simData();
    sim_optix.setup_simData();

    CompareRunners(sim_embree, sim_optix, N_rays_glob, "8", "4i_1_8_");
    CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "4i_1_12_");

}

//task 5a: 
//stinput differences: additional blocking heliostat 5321, receiver origin {0,0,171.035}
TEST(IsolatedHeliostatOptixEmbree, multiFacet1332_BlockingShading5a)
{
	// Make embree
	IsolatedHeliostatSimulationHelper<EmbreeRunner> sim_embree;
	sim_embree.runner.disable_stages();
	sim_embree.runner.set_number_of_threads(N_threads);
	sim_embree.initialize();
	sim_embree.seed = seed;
	sim_embree.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Make optix
	IsolatedHeliostatSimulationHelper<OptixRunner> sim_optix;
	sim_optix.initialize();
	sim_optix.seed = seed;
	sim_optix.sun_gen_type = SolTrace::Data::GenType::HALTON;

    glm::dvec3 origin = {0,0,171.035};
    sim_embree.set_rec_origin(origin);
    sim_optix.set_rec_origin(origin);

    std::vector<int> active {1332};
    std::vector<int> blocking {1266, 1304, 1306, 1334, 5321};
    sim_embree.create_active_heliostats(active);
    sim_embree.create_blocking_heliostats(blocking);
    sim_optix.create_active_heliostats(active);
    sim_optix.create_blocking_heliostats(blocking);

    sim_embree.assign_canted_banded(false);
    sim_embree.assign_canted_banded(true);
    sim_optix.assign_canted_banded(false);
    sim_optix.assign_canted_banded(true);

    sim_embree.assign_focal_lengths_banded(false);
    sim_embree.assign_focal_lengths_banded(true);
    sim_optix.assign_focal_lengths_banded(false);
    sim_optix.assign_focal_lengths_banded(true);

    sim_embree.setup_simData();
    sim_optix.setup_simData();

    CompareRunners(sim_embree, sim_optix, N_rays_glob, "8", "5a_1_8_");
}

//task 6a:
//stinput differences: receiver origin {0,0,171.035}
TEST(IsolatedHeliostatOptixEmbree, multiFacet1332_CantingAccuracy6a)
{
	// Make embree
	IsolatedHeliostatSimulationHelper<EmbreeRunner> sim_embree;
	sim_embree.runner.disable_stages();
	sim_embree.runner.set_number_of_threads(N_threads);
	sim_embree.initialize();
	sim_embree.seed = seed;
	sim_embree.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Make optix
	IsolatedHeliostatSimulationHelper<OptixRunner> sim_optix;
	sim_optix.initialize();
	sim_optix.seed = seed;
	sim_optix.sun_gen_type = SolTrace::Data::GenType::HALTON;

    sim_embree.set_flat_facets();
    sim_optix.set_flat_facets();

    glm::dvec3 origin = {0,0,171.035};
    sim_embree.set_rec_origin(origin);
    sim_optix.set_rec_origin(origin);

    std::vector<int> active {1332};
    sim_embree.create_active_heliostats(active);
    sim_optix.create_active_heliostats(active);

    sim_embree.assign_canted_banded(true);
    sim_optix.assign_canted_banded(true);

    sim_embree.setup_simData();
    sim_optix.setup_simData();

    CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "6a_1_12_");
	
}

//task 7a:
//stinput differences: reveiver origin {0,0,171.035}
//edits: shift in aimpoint due to truncation error in legacy when round robin study generated fluxmaps
TEST(IsolatedHeliostatOptixEmbree, multiFacet5473_CantingFocusingAccuracy7a)
{
	// Make embree
	IsolatedHeliostatSimulationHelper<EmbreeRunner> sim_embree;
	sim_embree.runner.disable_stages();
	sim_embree.runner.set_number_of_threads(N_threads);
	sim_embree.initialize();
	sim_embree.seed = seed;
	sim_embree.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Make optix
	IsolatedHeliostatSimulationHelper<OptixRunner> sim_optix;
	sim_optix.initialize();
	sim_optix.seed = seed;
	sim_optix.sun_gen_type = SolTrace::Data::GenType::HALTON;

	glm::dvec3 shift = {0,0,0.825};
    sim_embree.shift_aimpoint(shift);
    sim_optix.shift_aimpoint(shift);

    sim_embree.set_slope_error(0.0);
    sim_optix.set_slope_error(0.0);

    glm::dvec3 origin = {0,0,171.035};
    sim_embree.set_rec_origin(origin);
    sim_optix.set_rec_origin(origin);

    std::vector<int> active {5473};
    sim_embree.create_active_heliostats(active);
    sim_optix.create_active_heliostats(active);

    sim_embree.assign_canted_slant(true);
    sim_optix.assign_canted_slant(true);

    sim_embree.setup_simData();
    sim_optix.setup_simData();

    CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "7a_1_12_");
	
}

//task 7b:
//stinput differences: reveiver origin {0,0,171.035}
//edits: shift in aimpoint due to truncation error in legacy when round robin study generated fluxmaps
TEST(IsolatedHeliostatOptixEmbree, singleFacet5473_7b)
{
	// Make embree
	IsolatedHeliostatSimulationHelper<EmbreeRunner> sim_embree;
	sim_embree.runner.disable_stages();
	sim_embree.runner.set_number_of_threads(N_threads);
	sim_embree.initialize();
	sim_embree.seed = seed;
	sim_embree.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Make optix
	IsolatedHeliostatSimulationHelper<OptixRunner> sim_optix;
	sim_optix.initialize();
	sim_optix.seed = seed;
	sim_optix.sun_gen_type = SolTrace::Data::GenType::HALTON;

	glm::dvec3 shift = {0,0,0.825};
    sim_embree.shift_aimpoint(shift);
    sim_optix.shift_aimpoint(shift);

    sim_embree.set_slope_error(0.0);
    sim_optix.set_slope_error(0.0);

    glm::dvec3 origin = {0,0,171.035};
    sim_embree.set_rec_origin(origin);
    sim_optix.set_rec_origin(origin);

    std::vector<int> active {5473};
    sim_embree.create_active_heliostats(active);
    sim_optix.create_active_heliostats(active);

    sim_embree.setup_simData();
    sim_optix.setup_simData();

    CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "7b_1_12_");
	
}

//task 7c:
//stinput differences: reveiver origin {0,0,171.035}
//edits: shift in aimpoint due to truncation error in legacy when round robin study generated fluxmaps
TEST(IsolatedHeliostatOptixEmbree, singleFacet5473_BlockingShading7c)
{
	// Make embree
	IsolatedHeliostatSimulationHelper<EmbreeRunner> sim_embree;
	sim_embree.runner.disable_stages();
	sim_embree.runner.set_number_of_threads(N_threads);
	sim_embree.initialize();
	sim_embree.seed = seed;
	sim_embree.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Make optix
	IsolatedHeliostatSimulationHelper<OptixRunner> sim_optix;
	sim_optix.initialize();
	sim_optix.seed = seed;
	sim_optix.sun_gen_type = SolTrace::Data::GenType::HALTON;

	glm::dvec3 shift = {0,0,0.825};
    sim_embree.shift_aimpoint(shift);
    sim_optix.shift_aimpoint(shift);

    sim_embree.set_slope_error(0.0);
    sim_optix.set_slope_error(0.0);

    glm::dvec3 origin = {0,0,171.035};
    sim_embree.set_rec_origin(origin);
    sim_optix.set_rec_origin(origin);

    std::vector<int> active {5473};
    std::vector<int> blocking {5573};
    sim_embree.create_active_heliostats(active);
    sim_embree.create_blocking_heliostats(blocking);
    sim_optix.create_active_heliostats(active);
    sim_optix.create_blocking_heliostats(blocking);

    sim_embree.setup_simData();
    sim_optix.setup_simData();

    CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "7c_1_12_");
	
}

//7d: named 7b_8 in round robin files but matches 7d case
//stinput differences: reveiver origin {0,0,171.035}
TEST(IsolatedHeliostatOptixEmbree, singleFacet5473_7d) //named this in files, but refers to the same test as 7d
{
	// Make embree
	IsolatedHeliostatSimulationHelper<EmbreeRunner> sim_embree;
	sim_embree.runner.disable_stages();
	sim_embree.runner.set_number_of_threads(N_threads);
	sim_embree.initialize();
	sim_embree.seed = seed;
	sim_embree.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Make optix
	IsolatedHeliostatSimulationHelper<OptixRunner> sim_optix;
	sim_optix.initialize();
	sim_optix.seed = seed;
	sim_optix.sun_gen_type = SolTrace::Data::GenType::HALTON;

    sim_embree.set_slope_error(0.0);
    sim_optix.set_slope_error(0.0);

    glm::dvec3 origin = {0,0,171.035};
    sim_embree.set_rec_origin(origin);
    sim_optix.set_rec_origin(origin);

    std::vector<int> active {5473};
    sim_embree.create_active_heliostats(active);
    sim_optix.create_active_heliostats(active);

    sim_embree.setup_simData();
    sim_optix.setup_simData();

    CompareRunners(sim_embree, sim_optix, N_rays_glob, "8", "7b_1_8_");
	
}

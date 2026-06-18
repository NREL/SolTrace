#include "heliostat_field_test_template.hpp"

#include <embree_runner.hpp>
#include <optix_runner.hpp>

using SolTrace::EmbreeRunner::EmbreeRunner;

constexpr int N_rays_glob = 2e6;
constexpr int seed = 123;
constexpr bool save = false;
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

static void CompareRunners(HeliostatFieldSimulationHelper<EmbreeRunner>& sim_embree,
	HeliostatFieldSimulationHelper<OptixRunner>& sim_optix, int N_rays,
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
	double err_frac = 0.005;
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
	EXPECT_NEAR(sim_embree.heat_shield_absorb_count, sim_optix.heat_shield_absorb_count, err_abs);
	
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
	write_to_dict("05_heat_shield_absorb_count", sim_embree.heat_shield_absorb_count, sim_optix.heat_shield_absorb_count, dict_embree, dict_optix);
	write_to_dict("06_tot_rec_hits", sim_embree.tot_rec_hits, sim_optix.tot_rec_hits, dict_embree, dict_optix);
	write_to_dict("07_rec_direct_count", sim_embree.rec_direct_count, sim_optix.rec_direct_count, dict_embree, dict_optix);
	write_to_dict("08_rec_via_helio_count", sim_embree.rec_via_helio_count, sim_optix.rec_via_helio_count, dict_embree, dict_optix);

	// Helio hits add up
	EXPECT_EQ(sim_embree.tot_helio_hits, sim_embree.tot_helio_absorb_count + sim_embree.tot_reflect_count);
	EXPECT_EQ(sim_optix.tot_helio_hits, sim_optix.tot_helio_absorb_count + sim_optix.tot_reflect_count);

	// Receiver hits add up
	EXPECT_EQ(sim_embree.tot_rec_hits, sim_embree.rec_direct_count + sim_embree.rec_via_helio_count + sim_embree.rec_via_rec_count);
	EXPECT_EQ(sim_optix.tot_rec_hits, sim_optix.rec_direct_count + sim_optix.rec_via_helio_count + sim_optix.rec_via_rec_count);

	// Reflectivity
	double refl_embree = (double)sim_embree.tot_reflect_count / (double)sim_embree.tot_helio_hits;
	double refl_optix = (double)sim_optix.tot_reflect_count / (double)sim_optix.tot_helio_hits;
	EXPECT_NEAR(refl_embree, 0.9, 0.001);
	EXPECT_NEAR(refl_optix, 0.9, 0.001);

	write_to_dict("09_reflectivity", refl_embree, refl_optix, dict_embree, dict_optix);

	// Sun Count
	write_to_dict("10_sun_count", sim_embree.nsun_rays, sim_optix.nsun_rays, dict_embree, dict_optix);

	// Fraction reflected hits that hit receiver
	double frac_via_helio_a = (double)sim_embree.rec_absorb_count / (double)sim_embree.tot_reflect_count;
	double frac_via_helio_b = (double)sim_optix.rec_absorb_count / (double)sim_optix.tot_reflect_count;
	EXPECT_NEAR(frac_via_helio_a, frac_via_helio_b, err_frac);

	write_to_dict("11_frac_via_helio", frac_via_helio_a, frac_via_helio_b, dict_embree, dict_optix);

	// Compare power per ray
	write_to_dict("12_power_per_ray", sim_embree.power_per_ray, sim_optix.power_per_ray, dict_embree, dict_optix);

	// Total power absorbed
	double tol = 5.e-3;
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
	double rmse_tol = 0.08;  // 2% of peak flux // 0.04
	EXPECT_LE(rmse / peak_flux_embree, rmse_tol);

	// Average flux
	//EXPECT_NEAR(sim_embree.AveFlux / 1000.0, sim_optix.AveFlux / 1000.0, rmse_tol);

	write_to_dict("15_average_flux", sim_embree.AveFlux / 1000.0, sim_optix.AveFlux / 1000.0, dict_embree, dict_optix);

	// Binning
	//EXPECT_EQ(sim_embree.NotBinned, sim_optix.NotBinned);

	write_to_dict("16_not_binned", sim_embree.NotBinned, sim_optix.NotBinned, dict_embree, dict_optix);
	write_to_dict("17_neg_x_bin_err", sim_embree.max_neg_x_flux_err, sim_optix.max_neg_x_flux_err, dict_embree, dict_optix);
	write_to_dict("18_pos_x_bin_err", sim_embree.max_pos_x_flux_err, sim_optix.max_pos_x_flux_err, dict_embree, dict_optix);

	// Efficiencies
	EXPECT_NEAR(sim_embree.absorption_efficiency, sim_optix.absorption_efficiency, tol);
	EXPECT_NEAR(sim_embree.blocking_efficiency, sim_optix.blocking_efficiency, tol * 2.0);
	EXPECT_NEAR(sim_embree.spillage_efficiency, sim_optix.spillage_efficiency, tol);
	EXPECT_NEAR(sim_embree.cosine_efficiency, sim_optix.cosine_efficiency, tol);

	write_to_dict("19_absorption_efficiency", sim_embree.absorption_efficiency, sim_optix.absorption_efficiency, dict_embree, dict_optix);
	write_to_dict("20_blocking_efficiency", sim_embree.blocking_efficiency, sim_optix.blocking_efficiency, dict_embree, dict_optix);
	write_to_dict("21_spillage_efficiency", sim_embree.spillage_efficiency, sim_optix.spillage_efficiency, dict_embree, dict_optix);
	write_to_dict("22_cosine_efficiency", sim_embree.cosine_efficiency, sim_optix.cosine_efficiency, dict_embree, dict_optix);

	write_to_dict("23_sigmaflux", sim_embree.SigmaFlux, sim_optix.SigmaFlux, dict_embree, dict_optix);
	write_to_dict("24_uniformity", sim_embree.Uniformity, sim_optix.Uniformity, dict_embree, dict_optix);
	write_to_dict("25_rmse", rmse, rmse, dict_embree, dict_optix);
	double rmse_over_peak = rmse / (peak_flux_embree);
	write_to_dict("26_rmse_over_peak", rmse_over_peak, rmse_over_peak, dict_embree, dict_optix);

	write_to_dict("27_rec_via_rec_count", sim_embree.rec_via_rec_count, sim_optix.rec_via_rec_count, dict_embree, dict_optix);

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

TEST(HeliostatFieldOptixEmbree, singleFacet_SlantFocused)
{
	// Make embree
	HeliostatFieldSimulationHelper<EmbreeRunner> sim_embree;
	sim_embree.runner.disable_stages();
	sim_embree.runner.set_number_of_threads(N_threads);
	sim_embree.initialize();
	sim_embree.seed = seed;
	sim_embree.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Make optix
	HeliostatFieldSimulationHelper<OptixRunner> sim_optix;
	sim_optix.initialize();
	sim_optix.seed = seed;
	sim_optix.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Centerline aimpoints
	sim_embree.create_heliostat_field();
	sim_optix.create_heliostat_field();

	sim_embree.setup_simData();
	sim_optix.setup_simData();

	CompareRunners(sim_embree, sim_optix, N_rays_glob, "8", "1a_1_8_");
	CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "1a_1_12_");

	// Scatter aimpoints
	sim_embree.set_scatter_aimpoints();
	sim_optix.set_scatter_aimpoints();

	CompareRunners(sim_embree, sim_optix, N_rays_glob, "8", "1a_2_8_");
	CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "1a_2_12_");
}

TEST(HeliostatFieldOptixEmbree, singleFacet_BandFocused)
{
	// Make embree
	HeliostatFieldSimulationHelper<EmbreeRunner> sim_embree;
	sim_embree.runner.disable_stages();
	sim_embree.runner.set_number_of_threads(N_threads);
	sim_embree.initialize();
	sim_embree.seed = seed;
	sim_embree.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Make optix
	HeliostatFieldSimulationHelper<OptixRunner> sim_optix;
	sim_optix.initialize();
	sim_optix.seed = seed;
	sim_optix.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Centerline aimpoints
	sim_embree.create_heliostat_field();
	sim_optix.create_heliostat_field();

	sim_embree.assign_focal_lengths_canting_banded();
	sim_optix.assign_focal_lengths_canting_banded();

	sim_embree.setup_simData();
	sim_optix.setup_simData();

	CompareRunners(sim_embree, sim_optix, N_rays_glob, "8", "1b_1_8_");
	CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "1b_1_12_");

	// Scatter aimpoints
	sim_embree.set_scatter_aimpoints();
	sim_optix.set_scatter_aimpoints();

	CompareRunners(sim_embree, sim_optix, N_rays_glob, "8", "1b_2_8_");
	CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "1b_2_12_");
}

TEST(HeliostatFieldOptixEmbree, multiFacet_SlantCanted)
{
	// Make embree
	HeliostatFieldSimulationHelper<EmbreeRunner> sim_embree;
	sim_embree.runner.disable_stages();
	sim_embree.runner.set_number_of_threads(N_threads);
	sim_embree.initialize();
	sim_embree.seed = seed;
	sim_embree.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Make optix
	HeliostatFieldSimulationHelper<OptixRunner> sim_optix;
	sim_optix.initialize();
	sim_optix.seed = seed;
	sim_optix.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Centerline aimpoints
	sim_embree.create_heliostat_field();
	sim_optix.create_heliostat_field();

	sim_embree.assign_canted_slant(true);
	sim_optix.assign_canted_slant(true);

	sim_embree.setup_simData();
	sim_optix.setup_simData();

	CompareRunners(sim_embree, sim_optix, N_rays_glob, "8", "2a_1_8_");
	CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "2a_1_12_");

	// Scatter aimpoints
	sim_embree.set_scatter_aimpoints();
	sim_optix.set_scatter_aimpoints();

	CompareRunners(sim_embree, sim_optix, N_rays_glob, "8", "2a_2_8_");
	CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "2a_2_12_");
}

TEST(HeliostatFieldOptixEmbree, multiFacet_BandCanted)
{
	// Make embree
	HeliostatFieldSimulationHelper<EmbreeRunner> sim_embree;
	sim_embree.runner.disable_stages();
	sim_embree.runner.set_number_of_threads(N_threads);
	sim_embree.initialize();
	sim_embree.seed = seed;
	sim_embree.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Make optix
	HeliostatFieldSimulationHelper<OptixRunner> sim_optix;
	sim_optix.initialize();
	sim_optix.seed = seed;
	sim_optix.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Centerline aimpoints
	sim_embree.create_heliostat_field();
	sim_optix.create_heliostat_field();

	sim_embree.assign_canted_banded(true);
	sim_optix.assign_canted_banded(true);

	sim_embree.setup_simData();
	sim_optix.setup_simData();

	CompareRunners(sim_embree, sim_optix, N_rays_glob, "8", "2b_1_8_");
	CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "2b_1_12_");

	// Scatter aimpoints
	sim_embree.set_scatter_aimpoints();
	sim_optix.set_scatter_aimpoints();

	CompareRunners(sim_embree, sim_optix, N_rays_glob, "8", "2b_2_8_");
	CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "2b_2_12_");
}

TEST(HeliostatFieldOptixEmbree, multiFacet_SlantFocused_SlantCanted)
{
	// Make embree
	HeliostatFieldSimulationHelper<EmbreeRunner> sim_embree;
	sim_embree.runner.disable_stages();
	sim_embree.runner.set_number_of_threads(N_threads);
	sim_embree.initialize();
	sim_embree.seed = seed;
	sim_embree.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Make optix
	HeliostatFieldSimulationHelper<OptixRunner> sim_optix;
	sim_optix.initialize();
	sim_optix.seed = seed;
	sim_optix.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Centerline aimpoints
	sim_embree.create_heliostat_field();
	sim_optix.create_heliostat_field();

	sim_embree.assign_canted_slant(false);
	sim_optix.assign_canted_slant(false);

	sim_embree.setup_simData();
	sim_optix.setup_simData();

	CompareRunners(sim_embree, sim_optix, N_rays_glob, "8", "3a_1_8_");
	CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "3a_1_12_");

	// Scatter aimpoints
	sim_embree.set_scatter_aimpoints();
	sim_optix.set_scatter_aimpoints();

	CompareRunners(sim_embree, sim_optix, N_rays_glob, "8", "3a_2_8_");
	CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "3a_2_12_");
}

TEST(HeliostatFieldOptixEmbree, multiFacet_BandFocused_BandCanted)
{
	// Make embree
	HeliostatFieldSimulationHelper<EmbreeRunner> sim_embree;
	sim_embree.runner.disable_stages();
	sim_embree.runner.set_number_of_threads(N_threads);
	sim_embree.initialize();
	sim_embree.seed = seed;
	sim_embree.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Make optix
	HeliostatFieldSimulationHelper<OptixRunner> sim_optix;
	sim_optix.initialize();
	sim_optix.seed = seed;
	sim_optix.sun_gen_type = SolTrace::Data::GenType::HALTON;

	// Centerline aimpoints
	sim_embree.create_heliostat_field();
	sim_optix.create_heliostat_field();

	sim_embree.assign_canted_banded(false);    // Canted by band
	sim_embree.assign_focal_lengths_banded();  // Focused by band

	sim_optix.assign_canted_banded(false);     // Canted by band
	sim_optix.assign_focal_lengths_banded();   // Focused by band

	sim_embree.setup_simData();
	sim_optix.setup_simData();

	CompareRunners(sim_embree, sim_optix, N_rays_glob, "8", "3b_1_8_");
	CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "3b_1_12_");

	// Scatter aimpoints
	sim_embree.set_scatter_aimpoints();
	sim_optix.set_scatter_aimpoints();

	CompareRunners(sim_embree, sim_optix, N_rays_glob, "8", "3b_2_8_");
	CompareRunners(sim_embree, sim_optix, N_rays_glob, "12", "3b_2_12_");
}
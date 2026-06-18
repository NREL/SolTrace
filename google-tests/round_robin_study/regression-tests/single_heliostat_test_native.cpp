#include <native_runner.hpp>
#include "single_heliostat_test_template.hpp"

using SolTrace::NativeRunner::NativeRunner;
using NativeRunnerType = SolTrace::NativeRunner::NativeRunner;

using SingleHeliostatSimulationNative = SingleHeliostatSimulation<NativeRunnerType>;

static const int N_threads = static_cast<int>(std::max(1u, std::min(std::thread::hardware_concurrency(), 10u)));

const int N_rays_glob = 2e6;

namespace {
    const CompareRunnerOptions options = CompareRunnerOptions(
        123,    // seed
        true,   // ignore_direct
        false,  // save
        false   // save_flux
    );
}

TEST_F(SingleHeliostatSimulationNative, SingleFacetFlat_North)
{
    this->runner.set_number_of_threads(N_threads);
    setup_simData();
    simulate_check_outputs("1a", "N");
    EXPECT_NEAR(sun_width, 15.4557, 1.e-4);
    EXPECT_NEAR(sun_height, 15.4557, 1.e-4);
}

TEST_F(SingleHeliostatSimulationNative, SingleFacetFlat_Southeast)
{
    this->runner.set_number_of_threads(N_threads);
    set_heliostat_to_southeast();
    setup_simData();
    simulate_check_outputs("1a", "SE");
    EXPECT_NEAR(sun_width, 15.4557, 1.e-4);
    EXPECT_NEAR(sun_height, 15.4557, 1.e-4);
}

TEST_F(SingleHeliostatSimulationNative, SingleFacetFocused_North)
{
    this->runner.set_number_of_threads(N_threads);
    set_slant_focal_length();
    setup_simData();
    simulate_check_outputs("1b", "N");
    EXPECT_NEAR(sun_width, 15.4557, 1.e-4);
    EXPECT_NEAR(sun_height, 15.4557, 1.e-4);
}

TEST_F(SingleHeliostatSimulationNative, SingleFacetFocused_Southeast)
{
    this->runner.set_number_of_threads(N_threads);
    set_heliostat_to_southeast();
    set_slant_focal_length();
    setup_simData();
    simulate_check_outputs("1b", "SE");
    EXPECT_NEAR(sun_width, 15.4557, 1.e-4);
    EXPECT_NEAR(sun_height, 15.4557, 1.e-4);
}

TEST_F(SingleHeliostatSimulationNative, MultiFacetFlat_NoCanting_North)
{
    this->runner.set_number_of_threads(N_threads);
    set_flat_multi_facet();
    setup_simData();
    simulate_check_outputs("2", "N");
    EXPECT_NEAR(sun_width, 12.4214, 1.e-4);
    EXPECT_NEAR(sun_height, 10.4195, 1.e-4);
}

TEST_F(SingleHeliostatSimulationNative, MultiFacetFlat_NoCanting_Southeast)
{
    this->runner.set_number_of_threads(N_threads);
    set_flat_multi_facet();
    set_heliostat_to_southeast();
    setup_simData();
    simulate_check_outputs("2", "SE");
    EXPECT_NEAR(sun_width, 11.8574, 1.e-3);
    EXPECT_NEAR(sun_height, 11.5183, 1.e-3);
}

TEST_F(SingleHeliostatSimulationNative, MultiFacetFlat_SlantCanting_North)
{
    this->runner.set_number_of_threads(N_threads);
    set_onaxis_slant_canting();
    setup_simData();
    simulate_check_outputs("3", "N");
    EXPECT_NEAR(sun_width, 12.4214, 1.e-4);
    EXPECT_NEAR(sun_height, 10.4236, 1.e-4);
}

TEST_F(SingleHeliostatSimulationNative, MultiFacetFlat_SlantCanting_Southeast)
{
    this->runner.set_number_of_threads(N_threads);
    set_heliostat_to_southeast();
    set_onaxis_slant_canting();
    setup_simData();
    simulate_check_outputs("3", "SE");
    EXPECT_NEAR(sun_width, 11.8574, 1.e-3);
    EXPECT_NEAR(sun_height, 11.5183, 1.e-3);
}

TEST_F(SingleHeliostatSimulationNative, MultiFacetFocused_SlantCanting_North)
{
    this->runner.set_number_of_threads(N_threads);
    set_slant_focal_length();
    set_onaxis_slant_canting();
    setup_simData();
    simulate_check_outputs("4", "N");
    EXPECT_NEAR(sun_width, 12.4214, 1.e-4);
    EXPECT_NEAR(sun_height, 10.4236, 1.e-4);
}

TEST_F(SingleHeliostatSimulationNative, MultiFacetFocused_SlantCanting_Southeast)
{
    this->runner.set_number_of_threads(N_threads);
    set_heliostat_to_southeast();
    set_slant_focal_length();
    set_onaxis_slant_canting();
    setup_simData();
    simulate_check_outputs("4", "SE");
    EXPECT_NEAR(sun_width, 11.8574, 1.e-3);
    EXPECT_NEAR(sun_height, 11.5183, 1.e-3);
}

TEST(SingleHeliostatNative_ErrorTesting, SunShapePillboxUser)
{
    constexpr double pillbox_half_width = 4.65;

    const std::vector<double> user_angle = {
        0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.05, 1.2, 1.35,
        1.5, 1.65, 1.8, 1.95, 2.1, 2.25, 2.4, 2.55, 2.7, 2.85,
        3.0, 3.15, 3.3, 3.45, 3.6, 3.75, 3.9, 4.05, 4.2, 4.35,
        4.5, 4.65
    };

    const std::vector<double> user_intensity(user_angle.size(), 1.0);

    bool use_optical = false;
    bool use_sunshape = true;

    // Make user defined
    SingleHeliostatSimulationHelper<NativeRunner> sim_user;
    sim_user.runner.disable_stages(); // Disable stages
    sim_user.use_optical_errors = use_optical;
    sim_user.use_sunshape_errors = use_sunshape;
    sim_user.runner.set_number_of_threads(N_threads);
    sim_user.sun_shape = SunShape::USER_DEFINED;
    sim_user.user_angle = user_angle;
    sim_user.user_intensity = user_intensity;
    sim_user.initialize();

    // Make pillbox
    SingleHeliostatSimulationHelper<NativeRunner> sim_pillbox;
    sim_pillbox.runner.disable_stages();
    sim_pillbox.use_optical_errors = use_optical;
    sim_pillbox.use_sunshape_errors = use_sunshape;
    sim_pillbox.runner.set_number_of_threads(N_threads);
    sim_pillbox.sun_shape = SunShape::PILLBOX;
    sim_pillbox.half_width = pillbox_half_width;
    sim_pillbox.initialize();

    CompareRunners(sim_user, sim_pillbox, N_rays_glob, "", false, options);
}

TEST(SingleHeliostatNative_ErrorTesting, SunShapeGaussUser)
{
    constexpr double gauss_sigma = 2.0;

    const std::vector<double> user_angle = {
    0, 0.14999999999999999, 0.29999999999999999, 0.44999999999999996, 0.59999999999999998,
    0.75, 0.90000000000000002, 1.05, 1.2, 1.3499999999999999,
    1.4999999999999998, 1.6499999999999997, 1.7999999999999996, 1.9499999999999995, 2.0999999999999996,
    2.2499999999999996, 2.3999999999999995, 2.5499999999999994, 2.6999999999999993, 2.8499999999999992,
    2.9999999999999991, 3.149999999999999, 3.2999999999999989, 3.4499999999999988, 3.5999999999999988,
    3.7499999999999987, 3.8999999999999986, 4.0499999999999989, 4.1999999999999993, 4.3499999999999996,
    4.5, 4.6500000000000004, 4.8000000000000007, 4.9500000000000011, 5.1000000000000014,
    5.2500000000000018, 5.4000000000000021, 5.5500000000000025, 5.7000000000000028, 5.8500000000000032,
    6.0000000000000036, 6.1500000000000039, 6.3000000000000043, 6.4500000000000046, 6.600000000000005,
    6.7500000000000053, 6.9000000000000057, 7.050000000000006, 7.2000000000000064, 7.3500000000000068,
    7.4338443776996765
    };

    const std::vector<double> user_intensity = {
        1, 0.99719145137284493, 0.98881304461123309, 0.97500517529841779, 0.95599748183309996,
        0.93210249235952758, 0.90370707787319604, 0.87126203806373947, 0.835270211411272, 0.79627354712941933,
        0.75483960198900735, 0.71154792928753685, 0.66697681085847449, 0.62169074774719335, 0.5762290736718001,
        0.53109599103534533, 0.48675225595997185, 0.44360866071177202, 0.40202138309465507, 0.36228919676675592,
        0.32465246735834996, 0.28929379949956191, 0.25634015141507382, 0.22586619783851467, 0.19789869908361493,
        0.17242162389375304, 0.14938177525041826, 0.12869468021566238, 0.11025052530448531, 0.093919945799004464,
        0.079559508718227687, 0.067016762789207263, 0.056134762834133677, 0.046756008847947853, 0.038725770351664281,
        0.031894793362157031, 0.026121409853918164, 0.021273087559681509, 0.017227471311635049, 0.013872976053607069,
        0.011108996538242247, 0.0088458000785545891, 0.0070041671493623727, 0.0055148407637146818, 0.0043178400076330434,
        0.0033616864879322263, 0.0026025852527253772, 0.0020035944229486961, 0.0015338106793244468, 0.0011675911484396227,
        0.0010000000000000002
    };

    bool use_optical = false;
    bool use_sunshape = true;

    // Make user defined
    SingleHeliostatSimulationHelper<NativeRunner> sim_user;
    sim_user.runner.disable_stages(); // Disable stages
    sim_user.use_optical_errors = use_optical;
    sim_user.use_sunshape_errors = use_sunshape;
    sim_user.runner.set_number_of_threads(N_threads);
    sim_user.sun_shape = SunShape::USER_DEFINED;
    sim_user.user_angle = user_angle;
    sim_user.user_intensity = user_intensity;
    sim_user.initialize();

    // Make gaussian
    SingleHeliostatSimulationHelper<NativeRunner> sim_gauss;
    sim_gauss.runner.disable_stages();
    sim_gauss.use_optical_errors = use_optical;
    sim_gauss.use_sunshape_errors = use_sunshape;
    sim_gauss.runner.set_number_of_threads(N_threads);
    sim_gauss.sun_shape = SunShape::GAUSSIAN;
    sim_gauss.gauss_sigma = gauss_sigma;
    sim_gauss.initialize();

    CompareRunners(sim_user, sim_gauss, N_rays_glob, "", false, options);
}

TEST(SingleHeliostatNative_ErrorTesting, SunShapeBuieCsrUser)
{
    constexpr double csr = 0.1;
    const std::vector<double> user_angle = {
        0.0, 0.889795918, 1.779591837, 2.669387755, 3.559183673,
        4.448979592, 5.33877551, 6.228571429, 7.118367347, 8.008163265,
        8.897959184, 9.787755102, 10.67755102, 11.56734694, 12.45714286,
        13.34693878, 14.23673469, 15.12653061, 16.01632653, 16.90612245,
        17.79591837, 18.68571429, 19.5755102, 20.46530612, 21.35510204,
        22.24489796, 23.13469388, 24.0244898, 24.91428571, 25.80408163,
        26.69387755, 27.58367347, 28.47346939, 29.36326531, 30.25306122,
        31.14285714, 32.03265306, 32.92244898, 33.8122449, 34.70204082,
        35.59183673, 36.48163265, 37.37142857, 38.26122449, 39.15102041,
        40.04081633, 40.93061224, 41.82040816, 42.71020408, 43.6
    };

    const std::vector<double> user_intensity = {
        1.0, 0.995369248, 0.979933971, 0.947146323, 0.873326579,
        0.6031945, 0.025245752, 0.017129438, 0.012241349, 0.009101734,
        0.006982219, 0.005493449, 0.004413316, 0.003608278, 0.002994469,
        0.002517267, 0.002139962, 0.001837213, 0.001591112, 0.001388737,
        0.001220591, 0.001079582, 0.000960333, 0.000858713, 0.00077151,
        0.000696201, 0.000630778, 0.000573636, 0.000523476, 0.000479238,
        0.000440054, 0.000405206, 0.000374096, 0.000346224, 0.000321171,
        0.00029858, 0.000278149, 0.00025962, 0.000242771, 0.000227412,
        0.000213377, 0.000200523, 0.000188726, 0.000177877, 0.00016788,
        0.000158651, 0.000150115, 0.000142208, 0.000134871, 0.000128052
    };

    bool use_optical = false;
    bool use_sunshape = true;

    // Make user defined
    SingleHeliostatSimulationHelper<NativeRunner> sim_user;
    sim_user.runner.disable_stages(); // Disable stages
    sim_user.use_optical_errors = use_optical;
    sim_user.use_sunshape_errors = use_sunshape;
    sim_user.runner.set_number_of_threads(N_threads);
    sim_user.sun_shape = SunShape::USER_DEFINED;
    sim_user.user_angle = user_angle;
    sim_user.user_intensity = user_intensity;
    sim_user.initialize();

    // Make buie csr
    SingleHeliostatSimulationHelper<NativeRunner> sim_buie;
    sim_buie.runner.disable_stages();
    sim_buie.use_optical_errors = use_optical;
    sim_buie.use_sunshape_errors = use_sunshape;
    sim_buie.runner.set_number_of_threads(N_threads);
    sim_buie.sun_shape = SunShape::BUIE_CSR;
    sim_buie.csr = csr;
    sim_buie.initialize();
    CompareRunners(sim_user, sim_buie, N_rays_glob, "", false, options);
}

TEST(SingleHeliostatNative_ErrorTesting, SunShapeLimbUser)
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

    // Make user defined
    SingleHeliostatSimulationHelper<NativeRunner> sim_user;
    sim_user.runner.disable_stages(); // Disable stages
    sim_user.use_optical_errors = use_optical;
    sim_user.use_sunshape_errors = use_sunshape;
    sim_user.runner.set_number_of_threads(N_threads);
    sim_user.sun_shape = SunShape::USER_DEFINED;
    sim_user.user_angle = user_angle;
    sim_user.user_intensity = user_intensity;
    sim_user.initialize();

    // Make limb darkened
    SingleHeliostatSimulationHelper<NativeRunner> sim_limb;
    sim_limb.runner.disable_stages();
    sim_limb.use_optical_errors = use_optical;
    sim_limb.use_sunshape_errors = use_sunshape;
    sim_user.runner.set_number_of_threads(N_threads);
    sim_limb.sun_shape = SunShape::LIMBDARKENED;
    sim_limb.initialize();

    CompareRunners(sim_user, sim_limb, N_rays_glob, "", false, options);
}
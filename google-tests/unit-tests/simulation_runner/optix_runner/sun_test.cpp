#include <gtest/gtest.h>

#include <optix_runner.hpp>
#include <simulation_data.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result_export.hpp>

#include <algorithm>
#include <functional>

using SolTrace::Runner::RunnerStatus;

namespace
{
    constexpr double kSolarDiscAngleMrad = 4.65;
    constexpr double kLimbDarkeningCoeff = 0.5138;

    static double angle_mrad(const float3& a, const float3& b);
    static bool is_valid_dir(const float3& v);

    static double clamp01(double x)
    {
        if (x <= 0.0)
            return 0.0;
        if (x >= 1.0)
            return 1.0;
        return x;
    }

    // Source == 1, SunShape::GAUSSIAN:
    // thetax, thetay ~ N(0, sigma), theta = sqrt(thetax^2 + thetay^2)
    // => theta is Rayleigh-distributed with parameter sigma
    static double cdf_gaussian(double theta_mrad, double sigma_mrad)
    {
        if (theta_mrad <= 0.0)
            return 0.0;
        if (sigma_mrad <= 0.0)
            return 1.0;

        const double x = theta_mrad / sigma_mrad;
        return clamp01(1.0 - std::exp(-0.5 * x * x));
    }

    // Source == 1, SunShape::PILLBOX:
    // uniform over a disk of radius delop in (thetax, thetay)
    // => radial CDF is (theta / R)^2
    static double cdf_pillbox(double theta_mrad, double half_width_mrad)
    {
        if (theta_mrad <= 0.0)
            return 0.0;
        if (half_width_mrad <= 0.0)
            return 1.0;
        if (theta_mrad >= half_width_mrad)
            return 1.0;

        const double x = theta_mrad / half_width_mrad;
        return clamp01(x * x);
    }

    // Source == 1, SunShape::LIMBDARKENED:
    // accepted with stest = 1 - 0.5138 * (theta / MaxAngle)^4
    // sampled uniformly in area, so radial density is proportional to:
    // theta * stest(theta)
    static double cdf_limbdarkened(double theta_mrad, double max_angle_mrad)
    {
        if (theta_mrad <= 0.0)
            return 0.0;
        if (max_angle_mrad <= 0.0)
            return 1.0;
        if (theta_mrad >= max_angle_mrad)
            return 1.0;

        const double x = theta_mrad / max_angle_mrad;
        const double x2 = x * x;
        const double x6 = x2 * x2 * x2;

        // Integral of x * (1 - a x^4) dx from 0..x, normalized by 0..1
        const double numerator = 3.0 * x2 - kLimbDarkeningCoeff * x6;
        const double denominator = 3.0 - kLimbDarkeningCoeff;

        return clamp01(numerator / denominator);
    }

    // Source == 1, SunShape::BUIE_CSR:
    // stest(theta) matches tracing_errors.cpp exactly
    static double buie_intensity(double theta_mrad, double buie_kappa, double buie_gamma)
    {
        if (theta_mrad <= kSolarDiscAngleMrad)
            return std::cos(0.326 * theta_mrad) / std::cos(0.308 * theta_mrad);

        return std::exp(buie_kappa) * std::pow(std::abs(theta_mrad), buie_gamma);
    }

    // Radial density is proportional to theta * stest(theta)
    static double buie_radial_integral(double theta_mrad,
        double max_angle_mrad,
        double buie_kappa,
        double buie_gamma,
        int n_steps = 4096)
    {
        if (theta_mrad <= 0.0 || max_angle_mrad <= 0.0)
            return 0.0;

        if (theta_mrad > max_angle_mrad)
            theta_mrad = max_angle_mrad;

        const double h = theta_mrad / static_cast<double>(n_steps);
        double sum = 0.0;

        for (int i = 0; i <= n_steps; ++i)
        {
            const double t = h * static_cast<double>(i);
            const double f = t * buie_intensity(t, buie_kappa, buie_gamma);

            if (i == 0 || i == n_steps)
                sum += 0.5 * f;
            else
                sum += f;
        }

        return sum * h;
    }

    static double cdf_buie_csr(double theta_mrad,
        double max_angle_mrad,
        double buie_kappa,
        double buie_gamma)
    {
        if (theta_mrad <= 0.0)
            return 0.0;
        if (max_angle_mrad <= 0.0)
            return 1.0;
        if (theta_mrad >= max_angle_mrad)
            return 1.0;

        const double norm = buie_radial_integral(max_angle_mrad, max_angle_mrad, buie_kappa, buie_gamma);
        if (norm <= 0.0)
            return 0.0;

        const double value = buie_radial_integral(theta_mrad, max_angle_mrad, buie_kappa, buie_gamma) / norm;
        return clamp01(value);
    }

    static std::vector<double> collect_theta_mrad(const std::vector<float3>& dirs, const float3& sun_dir_nominal)
    {
        std::vector<double> thetas;
        thetas.reserve(dirs.size());

        for (const auto& d : dirs)
        {
            if (!is_valid_dir(d))
                continue;

            thetas.push_back(angle_mrad(d, sun_dir_nominal));
        }

        return thetas;
    }

    static double ks_statistic(const std::vector<double>& samples,
        const std::function<double(double)>& cdf)
    {
        if (samples.empty())
            return 1.0;

        std::vector<double> sorted = samples;
        std::sort(sorted.begin(), sorted.end());

        const double n = static_cast<double>(sorted.size());
        double d = 0.0;

        for (size_t i = 0; i < sorted.size(); ++i)
        {
            const double f = clamp01(cdf(sorted[i]));
            const double emp_lo = static_cast<double>(i) / n;
            const double emp_hi = static_cast<double>(i + 1) / n;
            d = std::max(d, std::abs(f - emp_lo));
            d = std::max(d, std::abs(emp_hi - f));
        }

        return d;
    }

    static double ks_pvalue_asymptotic(double d, size_t n)
    {
        if (n == 0)
            return 0.0;
        if (d <= 0.0)
            return 1.0;

        const double sqrtn = std::sqrt(static_cast<double>(n));
        const double x = (sqrtn + 0.12 + 0.11 / sqrtn) * d;

        double sum = 0.0;
        for (int k = 1; k <= 100; ++k)
        {
            const double term = std::exp(-2.0 * k * k * x * x);
            if (k % 2 == 1)
                sum += term;
            else
                sum -= term;

            if (term < 1.0e-12)
                break;
        }

        return clamp01(2.0 * sum);
    }

    static double ks_pvalue(const std::vector<double>& samples,
        const std::function<double(double)>& cdf)
    {
        return ks_pvalue_asymptotic(ks_statistic(samples, cdf), samples.size());
    }

    static float3 normalize_float3(const glm::dvec3& v)
    {
        float3 out = make_float3(
            static_cast<float>(v[0]),
            static_cast<float>(v[1]),
            static_cast<float>(v[2]));

        double nx = out.x;
        double ny = out.y;
        double nz = out.z;
        const double n = std::sqrt(nx * nx + ny * ny + nz * nz);
        if (n > 0.0)
        {
            out.x = static_cast<float>(nx / n);
            out.y = static_cast<float>(ny / n);
            out.z = static_cast<float>(nz / n);
        }

        out.x = -out.x;
        out.y = -out.y;
        out.z = -out.z;
        return out;
    }

    // Helper to compute angle (in mrad) between two vectors
    static double angle_mrad(const float3& a, const float3& b)
    {
        const double ax = static_cast<double>(a.x);
        const double ay = static_cast<double>(a.y);
        const double az = static_cast<double>(a.z);
        const double bx = static_cast<double>(b.x);
        const double by = static_cast<double>(b.y);
        const double bz = static_cast<double>(b.z);
        double dot = ax * bx + ay * by + az * bz;
        double na = std::sqrt(ax * ax + ay * ay + az * az);
        double nb = std::sqrt(bx * bx + by * by + bz * bz);
        if (na == 0.0 || nb == 0.0)
            return 0.0;
        dot /= (na * nb);
        if (dot > 1.0) dot = 1.0;
        if (dot < -1.0) dot = -1.0;
        double theta = std::acos(dot); // radians
        return theta * 1000.0; // mrad
    }

    // Simple validity check for direction vectors
    static bool is_valid_dir(const float3& v)
    {
        const double ax = v.x;
        const double ay = v.y;
        const double az = v.z;
        const double n2 = ax * ax + ay * ay + az * az;
        // Accept unit-ish vectors; reject zeros/huge
        return (n2 > 0.5 && n2 < 2.0);
    }

    static std::vector<float3> estimate_dirs_from_result(const SimulationResult& result)
    {
        std::vector<float3> dirs;
        dirs.reserve(result.get_number_of_records());

        for (int i = 0; i < result.get_number_of_records(); ++i)
        {
            ray_record_ptr rec = result[i];
            if (!rec)
                continue;

            const int n_interactions = rec->get_number_of_interactions();
            if (n_interactions < 2)
                continue;

            glm::dvec3 p0, p1;
            rec->get_position(0, p0);
            rec->get_position(1, p1);

            const double dx = p1[0] - p0[0];
            const double dy = p1[1] - p0[1];
            const double dz = p1[2] - p0[2];
            const double n = std::sqrt(dx * dx + dy * dy + dz * dz);
            if (n <= 0.0)
                continue;

            dirs.push_back(make_float3(
                static_cast<float>(dx / n),
                static_cast<float>(dy / n),
                static_cast<float>(dz / n)));
        }

        return dirs;
    }

    // Build a simple scene with a single flat plate receiver
    void make_default_sd_sun(SimulationData& sd, element_ptr& plate)
    {
        sd.clear();

        // Make stage
        auto stage = make_stage(0);
        stage->set_origin(0, 0, 0);
        stage->set_aim_vector(0, 0, 1);
        stage->set_name("stage");

        // Make reflective flat plate
        plate = make_element<SingleElement>();
        plate->set_origin(0, 0, 50);
        plate->set_aim_vector(0, 0, 100);  // Face up towards sun
        plate->set_surface(make_surface<Flat>());
        plate->set_aperture(make_aperture<Rectangle>(5, 5));
        InteractionType itype = InteractionType::REFLECTION;
        DistributionType dtype = DistributionType::NONE; // No errors
        double transmissivity = 0;
        double reflectivity = 1;
        double slope_err = 0;  // Error not supported
        double spec_err = 0;
        double ri_front = 0;   // Refraction not supported
        double ri_back = 0;
        OpticalProperties plate_optics(itype, dtype, transmissivity,
            reflectivity, slope_err, spec_err, ri_front, ri_back);
        plate->set_front_optical_properties(plate_optics);
        plate->set_back_optical_properties(plate_optics);
        plate->set_name("plate");

        // Add element to stage
        stage->add_element(plate);

        // Add stage to sd
        sd.add_stage(stage);

        // Set parameters
        SimulationParameters& params = sd.get_simulation_parameters();
        params.number_of_rays = 10000;
        params.max_number_of_rays = params.number_of_rays * 100;
        params.include_optical_errors = false;
        params.include_sun_shape_errors = false;
        params.seed = 123;
    }

    void count_hits_sun(const SimulationResult& result,
        int& absorbed_count, int& transmitted_count,
        int& reflected_count)
    {
        absorbed_count = 0;
        transmitted_count = 0;
        reflected_count = 0;
        int n_records = result.get_number_of_records();

        for (int i = 0; i < n_records; i++)
        {
            ray_record_ptr rec = result[i];

            int n_interactions = rec->get_number_of_interactions();
            for (int j = 0; j < n_interactions; j++)
            {
                RayEvent rev = rec->get_event(j);

                if (rev == RayEvent::ABSORB)
                    absorbed_count++;
                else if (rev == RayEvent::TRANSMIT)
                    transmitted_count++;
                else if (rev == RayEvent::REFLECT)
                    reflected_count++;
            }
        }
    }
}

// Basic smoke test: PILLBOX sun, no sunshape errors
TEST(Sun, SmokeTest)
{
    SimulationData sd;
    element_ptr plate;
    make_default_sd_sun(sd, plate);

    auto sun = make_ray_source<Sun>();
    sun->set_position(0, 0, 100);
    sun->set_shape(SolTrace::Data::SunShape::PILLBOX, 0.0, 4.65, 0.0);
    sd.add_ray_source(sun);

    // Run simulation
    OptixRunner runner;
    RunnerStatus sts = runner.initialize();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.setup_simulation(&sd);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.run_simulation();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    // Collect results
    SimulationResult result;
    sts = runner.report_simulation(&result, 0);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    int absorbed_count, transmitted_count, reflected_count;
    count_hits_sun(result, absorbed_count, transmitted_count, reflected_count);
    int total_hits = absorbed_count + transmitted_count + reflected_count;

    EXPECT_GT(result.get_number_of_records(), 0);
    EXPECT_GT(total_hits, 0);
}

// Helper to build SimulationData with a configurable sun
static void make_sun_sd(SimulationData& sd,
                        SolTrace::Data::SunShape shape,
                        double sigma_mrad,
                        double half_width_mrad,
                        bool include_sun_shape_errors,
                        const glm::dvec3& sun_pos,
                        double csr = 0)
{
    element_ptr plate;
    make_default_sd_sun(sd, plate);

    SimulationParameters& params = sd.get_simulation_parameters();
    params.include_sun_shape_errors = include_sun_shape_errors;

    auto sun = make_ray_source<Sun>();
    sun->set_position(sun_pos[0], sun_pos[1], sun_pos[2]);
    sun->set_shape(shape, sigma_mrad, half_width_mrad, csr);
    sd.add_ray_source(sun);
}

// Assign all sun shapes and verify which are accepted or rejected by OptixRunner::setup_simulation
TEST(Sun, SunShapeSupportMatrix)
{
    // Prepare a minimal simulation data with one element
    SimulationData sd;
    element_ptr plate;
    make_default_sd_sun(sd, plate);
    auto sun = make_ray_source<Sun>();
    sun->set_position(0, 0, 100);
    sd.add_ray_source(sun);

    // Enable sun shape errors so setup_sun enforces shape support
    SimulationParameters& params = sd.get_simulation_parameters();
    params.include_sun_shape_errors = true;

    for (int i = 0; i < static_cast<int>(SunShape::UNKNOWN); ++i)
    {
        SunShape shape = static_cast<SunShape>(i);
        bool is_supported = false;
        for (auto supported_shape : OptixCSP::kSupportedSunshapes)
        {
            if (shape == supported_shape)
            {
                is_supported = true;
                break;
            }
        }

        try
        {
            sun->set_shape(shape, 1, 1, 0.5, { 0 }, { 0 });
        }
        catch (...)
        {
            EXPECT_EQ(is_supported, false);
            continue;
        }

        OptixRunner runner;
        ASSERT_EQ(runner.initialize(), RunnerStatus::SUCCESS);
        RunnerStatus sts = runner.setup_simulation(&sd);

        if (is_supported)
            EXPECT_EQ(sts, RunnerStatus::SUCCESS);
        else
            EXPECT_EQ(sts, RunnerStatus::ERROR);
    }
}

TEST(Sun, GaussianSunAngleDistribution)
{
    const int N_RAYS = 200e3;
    const double SIGMA_MRAD = 4.65;

    SimulationData sd_gaussian;
    auto sun_pos = glm::dvec3(0.0, 0.0, 100.0);

    make_sun_sd(sd_gaussian, SolTrace::Data::SunShape::GAUSSIAN,
                SIGMA_MRAD, 0.0, true, sun_pos);

    sd_gaussian.get_simulation_parameters().number_of_rays = N_RAYS;
    sd_gaussian.get_simulation_parameters().max_number_of_rays = N_RAYS * 10;

    OptixRunner runner_gaussian;
    ASSERT_EQ(runner_gaussian.initialize(), RunnerStatus::SUCCESS);
    ASSERT_EQ(runner_gaussian.setup_simulation(&sd_gaussian), RunnerStatus::SUCCESS);
    ASSERT_EQ(runner_gaussian.run_simulation(), RunnerStatus::SUCCESS);

    SimulationResult result;
    ASSERT_EQ(runner_gaussian.report_simulation(&result, 0), RunnerStatus::SUCCESS);

    const std::vector<float3> dirs = estimate_dirs_from_result(result);
    EXPECT_FALSE(dirs.empty());

    const float3 sun_dir_nominal = normalize_float3(sun_pos);
    const std::vector<double> thetas = collect_theta_mrad(dirs, sun_dir_nominal);
    ASSERT_FALSE(thetas.empty());

    int count_valid = 0;
    int count_within_3sigma = 0;
    int count_beyond_sigma = 0;

    for (const double theta : thetas)
    {
        ++count_valid;
        if (theta <= 3.0 * SIGMA_MRAD)
            ++count_within_3sigma;
        if (theta > SIGMA_MRAD)
            ++count_beyond_sigma;
    }

    EXPECT_GT(count_valid, 0);

    const double frac_within_3sigma = static_cast<double>(count_within_3sigma) /
                                      static_cast<double>(count_valid);
    const double frac_beyond_sigma = static_cast<double>(count_beyond_sigma) /
                                     static_cast<double>(count_valid);

    EXPECT_GT(frac_within_3sigma, 0.98);
    EXPECT_GT(frac_beyond_sigma, 0.1);

    const double p_value = ks_pvalue(thetas, [SIGMA_MRAD](double theta)
        {
            return cdf_gaussian(theta, SIGMA_MRAD);
        });
    EXPECT_GT(p_value, 1.0e-6);
}

TEST(Sun, PillboxSunAngleDistribution)
{
    const int N_RAYS = 200e3;
    const double HALF_WIDTH_MRAD = 4.65;

    SimulationData sd_pillbox;
    auto sun_pos = glm::dvec3(0.0, 0.0, 100.0);

    make_sun_sd(sd_pillbox, SolTrace::Data::SunShape::PILLBOX,
                0.0, HALF_WIDTH_MRAD, true, sun_pos);

    sd_pillbox.get_simulation_parameters().number_of_rays = N_RAYS;
    sd_pillbox.get_simulation_parameters().max_number_of_rays = N_RAYS * 10;

    OptixRunner runner_pillbox;
    ASSERT_EQ(runner_pillbox.initialize(), RunnerStatus::SUCCESS);
    ASSERT_EQ(runner_pillbox.setup_simulation(&sd_pillbox), RunnerStatus::SUCCESS);
    ASSERT_EQ(runner_pillbox.run_simulation(), RunnerStatus::SUCCESS);

    SimulationResult result;
    ASSERT_EQ(runner_pillbox.report_simulation(&result, 0), RunnerStatus::SUCCESS);

    const std::vector<float3> dirs = estimate_dirs_from_result(result);
    EXPECT_FALSE(dirs.empty());

    const float3 sun_dir_nominal = normalize_float3(sun_pos);
    const std::vector<double> thetas = collect_theta_mrad(dirs, sun_dir_nominal);
    ASSERT_FALSE(thetas.empty());

    double max_theta = 0.0;
    int count_valid = 0;

    for (const double theta : thetas)
    {
        ++count_valid;
        if (theta > max_theta)
            max_theta = theta;
    }

    EXPECT_GT(count_valid, 0);
    EXPECT_LE(max_theta, HALF_WIDTH_MRAD + 0.1);

    const double p_value = ks_pvalue(thetas, [HALF_WIDTH_MRAD](double theta)
        {
            return cdf_pillbox(theta, HALF_WIDTH_MRAD);
        });
    EXPECT_GT(p_value, 1.0e-6);
}

TEST(Sun, LimbDarkenedSunAngleDistribution)
{
    const int N_RAYS = 200e3;
    const double MAX_ANGLE_MRAD = 4.65;
    const double UNIFORM_DISK_MEAN_RADIUS_FRAC = 2.0 / 3.0;

    SimulationData sd_limbdarkened;
    auto sun_pos = glm::dvec3(0.0, 0.0, 100.0);

    make_sun_sd(sd_limbdarkened, SolTrace::Data::SunShape::LIMBDARKENED,
                0.0, 0.0, true, sun_pos);

    sd_limbdarkened.get_simulation_parameters().number_of_rays = N_RAYS;
    sd_limbdarkened.get_simulation_parameters().max_number_of_rays = N_RAYS * 10;

    OptixRunner runner_limbdarkened;
    ASSERT_EQ(runner_limbdarkened.initialize(), RunnerStatus::SUCCESS);
    ASSERT_EQ(runner_limbdarkened.setup_simulation(&sd_limbdarkened), RunnerStatus::SUCCESS);
    ASSERT_EQ(runner_limbdarkened.run_simulation(), RunnerStatus::SUCCESS);

    SimulationResult result;
    ASSERT_EQ(runner_limbdarkened.report_simulation(&result, 0), RunnerStatus::SUCCESS);

    const std::vector<float3> dirs = estimate_dirs_from_result(result);
    EXPECT_FALSE(dirs.empty());

    const float3 sun_dir_nominal = normalize_float3(sun_pos);
    const std::vector<double> thetas = collect_theta_mrad(dirs, sun_dir_nominal);
    ASSERT_FALSE(thetas.empty());

    double max_theta = 0.0;
    int count_valid = 0;
    double theta_sum = 0.0;

    for (const double theta : thetas)
    {
        ++count_valid;
        theta_sum += theta;
        if (theta > max_theta)
            max_theta = theta;
    }

    EXPECT_GT(count_valid, 0);
    EXPECT_LE(max_theta, MAX_ANGLE_MRAD + 0.1);

    const double mean_theta = theta_sum / static_cast<double>(count_valid);
    const double mean_theta_frac = mean_theta / MAX_ANGLE_MRAD;

    EXPECT_LT(mean_theta_frac, UNIFORM_DISK_MEAN_RADIUS_FRAC);
    EXPECT_GT(mean_theta_frac, 0.55);

    const double p_value = ks_pvalue(thetas, [MAX_ANGLE_MRAD](double theta)
        {
            return cdf_limbdarkened(theta, MAX_ANGLE_MRAD);
        });
    EXPECT_GT(p_value, 1.0e-6);
}

TEST(Sun, BuieCSRSunAngleDistribution)
{
    const int N_RAYS = 200e3;
    const double MAX_ANGLE_MRAD = 43.6;
    const double CSR = 0.1;

    SimulationData sd_buie;
    auto sun_pos = glm::dvec3(0.0, 0.0, 100.0);

    make_sun_sd(sd_buie, SolTrace::Data::SunShape::BUIE_CSR,
                0.0, 0.0, true, sun_pos, CSR);

    sd_buie.get_simulation_parameters().number_of_rays = N_RAYS;
    sd_buie.get_simulation_parameters().max_number_of_rays = N_RAYS * 10;

    OptixRunner runner_buie;
    ASSERT_EQ(runner_buie.initialize(), RunnerStatus::SUCCESS);
    ASSERT_EQ(runner_buie.setup_simulation(&sd_buie), RunnerStatus::SUCCESS);
    ASSERT_EQ(runner_buie.run_simulation(), RunnerStatus::SUCCESS);

    SimulationResult result;
    ASSERT_EQ(runner_buie.report_simulation(&result, 0), RunnerStatus::SUCCESS);

    const std::vector<float3> dirs = estimate_dirs_from_result(result);
    EXPECT_FALSE(dirs.empty());

    const float3 sun_dir_nominal = normalize_float3(sun_pos);
    const std::vector<double> thetas = collect_theta_mrad(dirs, sun_dir_nominal);
    ASSERT_FALSE(thetas.empty());

    double max_theta = 0.0;
    int count_valid = 0;
    int count_beyond_disc = 0;

    for (const double theta : thetas)
    {
        ++count_valid;
        if (theta > max_theta)
            max_theta = theta;
        if (theta > 4.65)
            ++count_beyond_disc;
    }

    EXPECT_GT(count_valid, 0);
    EXPECT_LE(max_theta, MAX_ANGLE_MRAD + 0.5);

    const double frac_beyond_disc = static_cast<double>(count_beyond_disc) /
                                     static_cast<double>(count_valid);
    EXPECT_GT(frac_beyond_disc, 0.05);

    const double buie_gamma = 2.2 * std::log(0.52 * CSR) * std::pow(CSR, 0.43) - 0.1;
    const double buie_kappa = 0.9 * std::log(13.5 * CSR) * std::pow(CSR, -0.3);
    const double p_value = ks_pvalue(thetas, [MAX_ANGLE_MRAD, buie_kappa, buie_gamma](double theta)
        {
            return cdf_buie_csr(theta, MAX_ANGLE_MRAD, buie_kappa, buie_gamma);
        });
    EXPECT_GT(p_value, 1.0e-6);
}

TEST(Sun, UserDefinedSunAngleDistribution)
{
    const int N_RAYS = 200e3;
    const double MAX_ANGLE_MRAD = 7.95;

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

    SimulationData sd_user_defined;
    element_ptr plate;
    make_default_sd_sun(sd_user_defined, plate);

    SimulationParameters& params = sd_user_defined.get_simulation_parameters();
    params.include_sun_shape_errors = true;
    params.number_of_rays = N_RAYS;
    params.max_number_of_rays = N_RAYS * 10;

    auto sun_pos = glm::dvec3(0.0, 0.0, 100.0);
    auto sun = make_ray_source<Sun>();
    sun->set_position(sun_pos[0], sun_pos[1], sun_pos[2]);
    sun->set_shape(SolTrace::Data::SunShape::USER_DEFINED, 0.0, 0.0, 0.0, user_angle, user_intensity);
    sd_user_defined.add_ray_source(sun);

    OptixRunner runner_user_defined;
    ASSERT_EQ(runner_user_defined.initialize(), RunnerStatus::SUCCESS);
    ASSERT_EQ(runner_user_defined.setup_simulation(&sd_user_defined), RunnerStatus::SUCCESS);
    ASSERT_EQ(runner_user_defined.run_simulation(), RunnerStatus::SUCCESS);

    SimulationResult result;
    ASSERT_EQ(runner_user_defined.report_simulation(&result, 0), RunnerStatus::SUCCESS);

    const std::vector<float3> dirs = estimate_dirs_from_result(result);
    EXPECT_FALSE(dirs.empty());

    float3 sun_dir_nominal = make_float3(
        static_cast<float>(sun_pos[0]),
        static_cast<float>(sun_pos[1]),
        static_cast<float>(sun_pos[2]));
    {
        double nx = sun_dir_nominal.x;
        double ny = sun_dir_nominal.y;
        double nz = sun_dir_nominal.z;
        double n = std::sqrt(nx * nx + ny * ny + nz * nz);
        ASSERT_GT(n, 0.0);
        sun_dir_nominal.x = static_cast<float>(nx / n);
        sun_dir_nominal.y = static_cast<float>(ny / n);
        sun_dir_nominal.z = static_cast<float>(nz / n);
    }
    sun_dir_nominal.x = -sun_dir_nominal.x;
    sun_dir_nominal.y = -sun_dir_nominal.y;
    sun_dir_nominal.z = -sun_dir_nominal.z;

    double max_theta = 0.0;
    int count_valid = 0;
    int count_beyond_disc = 0;
    int count_beyond_table = 0;

    for (const auto& d : dirs)
    {
        if (!is_valid_dir(d))
            continue;

        double theta = angle_mrad(d, sun_dir_nominal);
        ++count_valid;
        if (theta > max_theta)
            max_theta = theta;
        if (theta > 4.65)
            ++count_beyond_disc;
        if (theta > MAX_ANGLE_MRAD)
            ++count_beyond_table;
    }

    EXPECT_GT(count_valid, 0);
    EXPECT_LE(max_theta, MAX_ANGLE_MRAD + 0.1);
    EXPECT_EQ(count_beyond_table, 0);

    const double frac_beyond_disc = static_cast<double>(count_beyond_disc) /
                                    static_cast<double>(count_valid);
    EXPECT_GT(frac_beyond_disc, 0.0);
    EXPECT_LT(frac_beyond_disc, 0.05);
}

#include <gtest/gtest.h>

#include <optix_runner.hpp>
#include <simulation_data.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result_export.hpp>

using SolTrace::Runner::RunnerStatus;

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

    int count_valid = 0;
    int count_within_3sigma = 0;
    int count_beyond_sigma = 0;

    for (const auto& d : dirs)
    {
        if (!is_valid_dir(d))
            continue;

        double theta = angle_mrad(d, sun_dir_nominal);
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

    for (const auto& d : dirs)
    {
        if (!is_valid_dir(d))
            continue;

        double theta = angle_mrad(d, sun_dir_nominal);
        ++count_valid;
        if (theta > max_theta)
            max_theta = theta;
    }

    EXPECT_GT(count_valid, 0);
    EXPECT_LE(max_theta, HALF_WIDTH_MRAD + 0.1);
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
    double theta_sum = 0.0;

    for (const auto& d : dirs)
    {
        if (!is_valid_dir(d))
            continue;

        double theta = angle_mrad(d, sun_dir_nominal);
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
}

TEST(Sun, BuieCSRSunAngleDistribution)
{
    const int N_RAYS = 200e3;
    const double MAX_ANGLE_MRAD = 43.6;

    SimulationData sd_buie;
    auto sun_pos = glm::dvec3(0.0, 0.0, 100.0);

    make_sun_sd(sd_buie, SolTrace::Data::SunShape::BUIE_CSR,
                0.0, 0.0, true, sun_pos, 0.1);

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
    }

    EXPECT_GT(count_valid, 0);
    EXPECT_LE(max_theta, MAX_ANGLE_MRAD + 0.5);

    const double frac_beyond_disc = static_cast<double>(count_beyond_disc) /
                                     static_cast<double>(count_valid);
    EXPECT_GT(frac_beyond_disc, 0.05);
}



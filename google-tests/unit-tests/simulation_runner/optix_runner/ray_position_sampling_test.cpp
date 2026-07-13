// ray_position_sampling_test.cpp
//
// Tests for sun parallelogram position sampling (GenType::HALTON and
// GenType::RANDOM).  Each test fires rays straight down onto a large flat plate
// with no sun-shape or optical errors, then inspects position[0] of every ray
// record — the ray-generation point on the sun parallelogram — to verify the
// spatial distribution.

#include <gtest/gtest.h>

#include <optical_properties.hpp>
#include <optix_runner.hpp>
#include <simulation_data.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result_export.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>

using SolTrace::Runner::RunnerStatus;

namespace
{
    // ---- Minimal KS-test infrastructure ------------------------------------
    // (duplicated from sun_test.cpp — kept local to this TU)

    static double clamp01(double x)
    {
        return x <= 0.0 ? 0.0 : (x >= 1.0 ? 1.0 : x);
    }

    static double ks_statistic(const std::vector<double>& samples,
                               const std::function<double(double)>& cdf)
    {
        if (samples.empty()) return 1.0;
        std::vector<double> sorted = samples;
        std::sort(sorted.begin(), sorted.end());
        const double n = static_cast<double>(sorted.size());
        double d = 0.0;
        for (size_t i = 0; i < sorted.size(); ++i)
        {
            const double f      = clamp01(cdf(sorted[i]));
            const double emp_lo = static_cast<double>(i) / n;
            const double emp_hi = static_cast<double>(i + 1) / n;
            d = std::max(d, std::abs(f - emp_lo));
            d = std::max(d, std::abs(emp_hi - f));
        }
        return d;
    }

    static double ks_pvalue_asymptotic(double d, size_t n)
    {
        if (n == 0) return 0.0;
        if (d <= 0.0) return 1.0;
        const double sqrtn = std::sqrt(static_cast<double>(n));
        const double x = (sqrtn + 0.12 + 0.11 / sqrtn) * d;
        double sum = 0.0;
        for (int k = 1; k <= 100; ++k)
        {
            const double term = std::exp(-2.0 * k * k * x * x);
            sum += (k % 2 == 1) ? term : -term;
            if (term < 1.0e-12) break;
        }
        return clamp01(2.0 * sum);
    }

    static double ks_pvalue(const std::vector<double>& samples,
                            const std::function<double(double)>& cdf)
    {
        return ks_pvalue_asymptotic(ks_statistic(samples, cdf), samples.size());
    }

    // Test that `coords` (a 1-D projection of source positions) is consistent
    // with a uniform distribution.  Coordinates are normalised to [0,1] using
    // the empirical range before the KS test.  For N >= 1 000 from a true
    // uniform distribution the normalisation bias is O(1/N) — negligible
    // against the ~1% critical value.
    static double ks_pvalue_uniform1d(const std::vector<double>& coords)
    {
        if (coords.size() < 2) return 0.0;
        const double lo = *std::min_element(coords.begin(), coords.end());
        const double hi = *std::max_element(coords.begin(), coords.end());
        if (hi - lo < 1.0e-9) return 0.0;

        std::vector<double> u;
        u.reserve(coords.size());
        for (double c : coords)
            u.push_back((c - lo) / (hi - lo));

        return ks_pvalue(u, [](double x) { return x; });
    }

    // ---- Scene helpers -----------------------------------------------------

    // Build a SimulationData with a large flat plate at z=50 (200 × 200 world
    // units), perfectly reflective, no optical or sun-shape errors.  The plate
    // is intentionally much larger than any realistic sun parallelogram so that
    // all generated rays hit it and every source position is recorded.
    static void make_large_plate_scene(SimulationData& sd, int seed = 42)
    {
        sd.clear();

        auto stage = make_stage(0);
        stage->set_origin(0, 0, 0);
        stage->set_aim_vector(0, 0, 1);
        stage->set_name("stage");

        auto plate = make_element<SingleElement>();
        plate->set_origin(0, 0, 50);
        plate->set_aim_vector(0, 0, 100);
        plate->set_surface(make_surface<Flat>());
        plate->set_aperture(make_aperture<Rectangle>(200, 200));
        SolTrace::Data::OpticalPropertySet plate_optics(
            SolTrace::Data::InteractionType::REFLECTION,
            0.0,
            0.0,
            "ray_position_sampling_plate_optics");
        plate_optics.set_ideal_reflection(OpticalSide::Both);
        auto plate_optics_ref = sd.add_optical_property_set(plate_optics);
        plate->set_optical_property_set(plate_optics_ref);
        plate->set_name("plate");

        stage->add_element(plate);
        sd.add_stage(stage);

        SimulationParameters& params       = sd.get_simulation_parameters();
        params.number_of_rays              = 20000;
        params.max_number_of_rays          = params.number_of_rays * 10;
        params.include_optical_errors      = false;
        params.include_sun_shape_errors    = false;
        params.seed                        = seed;
    }

    // Add a sun pointing straight down (position overhead at z=100) to an
    // already-configured SimulationData, using the requested gen_type.
    static void add_sun(SimulationData& sd, SolTrace::Data::GenType gen_type)
    {
        auto sun = make_ray_source<Sun>();
        sun->set_position(0, 0, 100);
        // PILLBOX shape with no sun-shape errors keeps rays exactly parallel,
        // so hit positions on the plate are the source positions projected down.
        sun->set_shape(SolTrace::Data::SunShape::PILLBOX, 0.0, 4.65, 0.0);
        sun->set_gen_type(gen_type);
        sd.add_ray_source(sun);
    }

    // Run the simulation and populate `result`.  Returns false on any failure.
    static bool run_sim(SimulationData& sd, RayHistoryResult& result)
    {
        OptixRunner runner;
        if (runner.initialize()           != RunnerStatus::SUCCESS) return false;
        if (runner.setup_simulation(&sd)  != RunnerStatus::SUCCESS) return false;
        if (runner.run_simulation()       != RunnerStatus::SUCCESS) return false;
        if (runner.report_simulation(&result, 0) != RunnerStatus::SUCCESS) return false;
        return true;
    }

    // Extract the X and Y components of position[0] (the ray-generation point)
    // for every ray that has at least one surface interaction recorded.
    static void collect_source_xy(const RayHistoryResult& result,
                                  std::vector<double>& xs,
                                  std::vector<double>& ys)
    {
        xs.clear();
        ys.clear();
        for (int i = 0; i < result.get_number_of_records(); ++i)
        {
            auto rec = result[i];
            if (!rec || rec->get_number_of_interactions() < 1) continue;
            glm::dvec3 p;
            rec->get_position(0, p);
            xs.push_back(p[0]);
            ys.push_back(p[1]);
        }
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

// GenType::RANDOM: the marginal X and Y distributions of source positions
// should be statistically consistent with a uniform distribution over the
// sun parallelogram.
TEST(RayPositionSampling, Random_UniformMarginals)
{
    SimulationData sd;
    make_large_plate_scene(sd);
    add_sun(sd, SolTrace::Data::GenType::RANDOM);

    RayHistoryResult result;
    ASSERT_TRUE(run_sim(sd, result));
    ASSERT_GT(result.get_number_of_records(), 0);

    std::vector<double> xs, ys;
    collect_source_xy(result, xs, ys);
    ASSERT_GE(xs.size(), 1000u);

    const double p_x = ks_pvalue_uniform1d(xs);
    const double p_y = ks_pvalue_uniform1d(ys);

    EXPECT_GT(p_x, 1.0e-6) << "X marginal deviates significantly from uniform";
    EXPECT_GT(p_y, 1.0e-6) << "Y marginal deviates significantly from uniform";
}

// GenType::HALTON: same uniformity requirement, via the Halton low-discrepancy
// sequence (bases 2 and 3 for the two parallelogram axes).
TEST(RayPositionSampling, Halton_UniformMarginals)
{
    SimulationData sd;
    make_large_plate_scene(sd);
    add_sun(sd, SolTrace::Data::GenType::HALTON);

    RayHistoryResult result;
    ASSERT_TRUE(run_sim(sd, result));
    ASSERT_GT(result.get_number_of_records(), 0);

    std::vector<double> xs, ys;
    collect_source_xy(result, xs, ys);
    ASSERT_GE(xs.size(), 1000u);

    const double p_x = ks_pvalue_uniform1d(xs);
    const double p_y = ks_pvalue_uniform1d(ys);

    EXPECT_GT(p_x, 1.0e-6) << "X marginal deviates significantly from uniform";
    EXPECT_GT(p_y, 1.0e-6) << "Y marginal deviates significantly from uniform";
}

// GenType::HALTON is deterministic: the seed field has no effect because the
// sequence depends only on the ray index.  Two runs with different seeds must
// produce the same set of source positions.
TEST(RayPositionSampling, Halton_Deterministic)
{
    // Run 1: seed = 1
    SimulationData sd1;
    make_large_plate_scene(sd1, /*seed=*/1);
    add_sun(sd1, SolTrace::Data::GenType::HALTON);

    RayHistoryResult r1;
    ASSERT_TRUE(run_sim(sd1, r1));

    // Run 2: seed = 99999 (different, should be ignored by Halton)
    SimulationData sd2;
    make_large_plate_scene(sd2, /*seed=*/99999);
    add_sun(sd2, SolTrace::Data::GenType::HALTON);

    RayHistoryResult r2;
    ASSERT_TRUE(run_sim(sd2, r2));

    ASSERT_EQ(r1.get_number_of_records(), r2.get_number_of_records());

    std::vector<double> xs1, ys1, xs2, ys2;
    collect_source_xy(r1, xs1, ys1);
    collect_source_xy(r2, xs2, ys2);
    ASSERT_EQ(xs1.size(), xs2.size());

    // Sort both and compare element-by-element.
    std::sort(xs1.begin(), xs1.end());
    std::sort(xs2.begin(), xs2.end());
    std::sort(ys1.begin(), ys1.end());
    std::sort(ys2.begin(), ys2.end());

    for (size_t i = 0; i < xs1.size(); ++i)
    {
        EXPECT_NEAR(xs1[i], xs2[i], 1.0e-4)
            << "Halton X differs at sorted index " << i
            << " — sequence is not deterministic across seeds";
        EXPECT_NEAR(ys1[i], ys2[i], 1.0e-4)
            << "Halton Y differs at sorted index " << i
            << " — sequence is not deterministic across seeds";
    }
}

// GenType::RANDOM IS seed-dependent: two runs with different seeds must
// produce distinguishably different position sets.
TEST(RayPositionSampling, Random_SeedDependent)
{
    SimulationData sd1;
    make_large_plate_scene(sd1, /*seed=*/1);
    add_sun(sd1, SolTrace::Data::GenType::RANDOM);

    SimulationData sd2;
    make_large_plate_scene(sd2, /*seed=*/2);
    add_sun(sd2, SolTrace::Data::GenType::RANDOM);

    RayHistoryResult r1, r2;
    ASSERT_TRUE(run_sim(sd1, r1));
    ASSERT_TRUE(run_sim(sd2, r2));
    ASSERT_GT(r1.get_number_of_records(), 0);
    ASSERT_GT(r2.get_number_of_records(), 0);

    std::vector<double> xs1, ys1, xs2, ys2;
    collect_source_xy(r1, xs1, ys1);
    collect_source_xy(r2, xs2, ys2);
    ASSERT_FALSE(xs1.empty());
    ASSERT_FALSE(xs2.empty());

    const double mean_x1 = std::accumulate(xs1.begin(), xs1.end(), 0.0) / xs1.size();
    const double mean_x2 = std::accumulate(xs2.begin(), xs2.end(), 0.0) / xs2.size();
    const double mean_y1 = std::accumulate(ys1.begin(), ys1.end(), 0.0) / ys1.size();
    const double mean_y2 = std::accumulate(ys2.begin(), ys2.end(), 0.0) / ys2.size();

    // With 20 000 rays and different seeds the probability that both means
    // coincide to within 1e-4 is negligible.
    const bool x_differs = std::abs(mean_x1 - mean_x2) > 1.0e-4;
    const bool y_differs = std::abs(mean_y1 - mean_y2) > 1.0e-4;

    EXPECT_TRUE(x_differs || y_differs)
        << "Seeds 1 and 2 produced indistinguishable mean source positions; "
           "seed variation may not be wired up for RANDOM gen_type.";
}

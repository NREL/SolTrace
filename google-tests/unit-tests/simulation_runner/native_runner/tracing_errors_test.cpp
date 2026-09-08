#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <vector>

#include <constants.hpp>
#include <error_distributions.hpp>
#include <optical_properties.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result_export.hpp>

#include <mtrand.hpp>
#include <native_runner.hpp>
#include <native_runner_types.hpp>
#include <tracing_errors.hpp>

using SolTrace::Data::DistributionType;
using SolTrace::Data::InteractionType;
using SolTrace::Data::OpticalPropertySet;
using SolTrace::Data::OpticalSide;
using SolTrace::Data::SunShape;
using SolTrace::NativeRunner::ApplySlopeError;
using SolTrace::NativeRunner::ApplySpecularityError;
using SolTrace::NativeRunner::ApplySunShape;
using SolTrace::NativeRunner::MTRand;
using SolTrace::NativeRunner::TSun;

namespace
{

constexpr int      kSamples  = 200000;
constexpr uint32_t kSeed     = 20260902;
constexpr double   kFrameTol = 1e-12;

// Axes chosen to cover every branch of the frame construction in
// tracing_errors.cpp, including the two axis.z == 0 degenerate cases.
const std::vector<glm::dvec3> kTestAxes = {
    glm::normalize(glm::dvec3(0.0, 0.0, 1.0)),   // +z
    glm::normalize(glm::dvec3(0.0, 0.0, -1.0)),  // -z
    glm::normalize(glm::dvec3(0.0, 1.0, 0.0)),   // z == 0, x == 0
    glm::normalize(glm::dvec3(1.0, 0.0, 0.0)),   // z == 0, x != 0
    glm::normalize(glm::dvec3(1.0, 1.0, 0.0)),   // z == 0, x != 0
    glm::normalize(glm::dvec3(1.0, -2.0, 3.0)),  // generic
    glm::normalize(glm::dvec3(-0.3, 0.7, -0.5)), // generic
};

double AngleBetween(const glm::dvec3& a, const glm::dvec3& b)
{
    const double c = glm::dot(glm::normalize(a), glm::normalize(b));
    return std::acos(std::clamp(c, -1.0, 1.0));
}

OpticalPropertySet MakeReflector(DistributionType   dist,
                                 double             slope_mrad,
                                 double             spec_mrad,
                                 const std::string& name = "test_optics")
{
    OpticalPropertySet optics(InteractionType::REFLECTION, name);
    // set_ideal_reflection() clears the error terms, so it must come first.
    optics.set_ideal_reflection(OpticalSide::Both);
    optics.set_errors(OpticalSide::Both, dist, slope_mrad, spec_mrad);
    return optics;
}

TSun MakeSun(SunShape shape, double sigma_mrad, double max_angle_mrad)
{
    TSun sun;
    sun.ShapeIndex   = shape;
    sun.Sigma        = sigma_mrad;
    sun.MaxAngle     = max_angle_mrad;
    sun.MaxIntensity = 1.0;
    sun.buie_kappa   = 0.0;
    sun.buie_gamma   = 0.0;
    return sun;
}

// Buie (2003) shape parameters, mirroring Sun::calculate_buie_parameters().
void BuieParameters(double csr, double& kappa, double& gamma)
{
    double chi;
    if (csr > 0.145)
        chi = -0.04419909985804843 +
              csr * (1.401323894233574 +
                     csr * (-0.3639746714505299 +
                            csr * (-0.9579768560161194 +
                                   1.1550475450828657 * csr)));
    else if (csr > 0.035)
        chi = 0.022652077593662934 +
              csr * (0.5252380349996234 +
                     (2.5484334534423887 - 0.8763755326550412 * csr) * csr);
    else
        chi = 0.004733749294807862 +
              csr * (4.716738065192151 +
                     csr * (-463.506669149804 +
                            csr * (24745.88727411664 +
                                   csr * (-606122.7511711778 +
                                          5521693.445014727 * csr))));

    kappa = 0.9 * std::log(13.5 * chi) * std::pow(chi, -0.3);
    gamma = 2.2 * std::log(0.52 * chi) * std::pow(chi, 0.43) - 0.1;
}

// Numerically integrated CDF of the polar angle for a rejection-sampled shape.
// The sampler draws uniformly in a (theta_x, theta_y) box and accepts with
// probability s(theta), so the radial density is proportional to
// theta*s(theta).
class RadialCdf
{
public:
    RadialCdf(std::function<double(double)> intensity,
              double                        max_angle,
              int                           bins = 20000)
        : m_max(max_angle), m_cdf(bins + 1, 0.0)
    {
        const double h = max_angle / bins;
        for (int i = 1; i <= bins; ++i)
        {
            const double t0 = (i - 1) * h;
            const double t1 = i * h;
            const double w0 = t0 * intensity(t0);
            const double w1 = t1 * intensity(t1);
            m_cdf[i]        = m_cdf[i - 1] + 0.5 * (w0 + w1) * h;
        }
        const double total = m_cdf.back();
        for (double& v : m_cdf)
            v /= total;
    }

    double operator()(double theta) const
    {
        if (theta <= 0.0) return 0.0;
        if (theta >= m_max) return 1.0;

        const double x  = theta / m_max * (m_cdf.size() - 1);
        const int    i  = static_cast<int>(x);
        const double fr = x - i;
        return m_cdf[i] * (1.0 - fr) + m_cdf[i + 1] * fr;
    }

private:
    double              m_max;
    std::vector<double> m_cdf;
};

// Two-sided Kolmogorov-Smirnov statistic against an analytic CDF.
double KsStatistic(std::vector<double>                  samples,
                   const std::function<double(double)>& cdf)
{
    std::sort(samples.begin(), samples.end());
    const double n = static_cast<double>(samples.size());

    double d = 0.0;
    for (size_t i = 0; i < samples.size(); ++i)
    {
        const double f = cdf(samples[i]);
        d              = std::max(d, std::max((i + 1) / n - f, f - i / n));
    }
    return d;
}

// Generous enough to avoid flakes, tight enough to catch a wrong distribution.
double KsThreshold(int n)
{ return 2.5 / std::sqrt(static_cast<double>(n)); }

// Sun-shape sample angles, in mrad (ApplySunShape() returns radians).
std::vector<double> SampleSunAngles(TSun& sun, int n, uint32_t seed = kSeed)
{
    MTRand           rng(seed);
    const glm::dvec3 axis(0.0, 0.0, 1.0);

    std::vector<double> angles;
    angles.reserve(n);
    for (int i = 0; i < n; ++i)
    {
        const glm::dvec3 out = ApplySunShape(rng, axis, sun);
        angles.push_back(AngleBetween(axis, out) * 1.0e3);
    }
    return angles;
}

// Surface-error sample angles, in mrad.
std::vector<double> SampleSurfaceAngles(const OpticalPropertySet& optics,
                                        int                       n,
                                        uint32_t                  seed = kSeed)
{
    MTRand           rng(seed);
    const glm::dvec3 axis(0.0, 0.0, 1.0);

    std::vector<double> angles;
    angles.reserve(n);
    for (int i = 0; i < n; ++i)
    {
        const glm::dvec3 out =
            ApplySpecularityError(rng, axis, optics, false, axis);
        angles.push_back(AngleBetween(axis, out) * 1.0e3);
    }
    return angles;
}

// Slope-error sample angles, in mrad.
std::vector<double> SampleSlopeAngles(const OpticalPropertySet& optics,
                                      int                       n,
                                      uint32_t                  seed = kSeed)
{
    MTRand           rng(seed);
    const glm::dvec3 axis(0.0, 0.0, 1.0);

    std::vector<double> angles;
    angles.reserve(n);
    for (int i = 0; i < n; ++i)
    {
        const glm::dvec3 out = ApplySlopeError(rng, axis, optics, false);
        angles.push_back(AngleBetween(axis, out) * 1.0e3);
    }
    return angles;
}

double RootMeanSquare(const std::vector<double>& v)
{
    double acc = 0.0;
    for (double x : v)
        acc += x * x;
    return std::sqrt(acc / v.size());
}

double Mean(const std::vector<double>& v)
{
    double acc = 0.0;
    for (double x : v)
        acc += x;
    return acc / v.size();
}

} // namespace

// ---------------------------------------------------------------------------
// Frame construction
// ---------------------------------------------------------------------------

// A zero perturbation must return the axis itself. This exercises the
// round trip through the ray frame for every branch of its construction.
TEST(TracingErrors, ZeroPerturbationReturnsAxisForAllAxes)
{
    MTRand                   rng(kSeed);
    const OpticalPropertySet optics =
        MakeReflector(DistributionType::NONE, 0.0, 0.0);

    for (const glm::dvec3& axis : kTestAxes)
    {
        const glm::dvec3 out = ApplySlopeError(rng, axis, optics, false);

        SCOPED_TRACE(testing::Message() << "axis (" << axis.x << ", " << axis.y
                                        << ", " << axis.z << ")");
        EXPECT_NEAR(out.x, axis.x, kFrameTol);
        EXPECT_NEAR(out.y, axis.y, kFrameTol);
        EXPECT_NEAR(out.z, axis.z, kFrameTol);
    }
}

TEST(TracingErrors, PerturbedDirectionIsUnitLengthForAllAxes)
{
    MTRand                   rng(kSeed);
    const OpticalPropertySet optics =
        MakeReflector(DistributionType::GAUSSIAN, 5.0, 5.0);
    TSun sun = MakeSun(SunShape::GAUSSIAN, 2.73, 4.65);

    for (const glm::dvec3& axis : kTestAxes)
    {
        for (int i = 0; i < 1000; ++i)
        {
            EXPECT_NEAR(glm::length(ApplySlopeError(rng, axis, optics, false)),
                        1.0,
                        1e-12);
            EXPECT_NEAR(glm::length(ApplySpecularityError(
                            rng, axis, optics, false, axis)),
                        1.0,
                        1e-12);
            EXPECT_NEAR(glm::length(ApplySunShape(rng, axis, sun)), 1.0, 1e-12);
        }
    }
}

// The pillbox sampler is bounded, so the perturbation can never exceed the
// configured half-width regardless of which frame branch was taken.
TEST(TracingErrors, PillboxRespectsHalfWidthForAllAxes)
{
    const double             kSlopeMrad = 12.0;
    MTRand                   rng(kSeed);
    const OpticalPropertySet optics =
        MakeReflector(DistributionType::PILLBOX, kSlopeMrad, kSlopeMrad);

    for (const glm::dvec3& axis : kTestAxes)
    {
        for (int i = 0; i < 5000; ++i)
        {
            const glm::dvec3 out = ApplySlopeError(rng, axis, optics, false);
            EXPECT_LE(AngleBetween(axis, out) * 1.0e3, kSlopeMrad + 1e-9);
        }
    }
}

// ---------------------------------------------------------------------------
// Independence invariants that make splitting Errors() safe
// ---------------------------------------------------------------------------

// ---------------------------------------------------------------------------
// Rejection of perturbations that pass through an opaque surface
// ---------------------------------------------------------------------------

TEST(TracingErrors, RejectionKeepsRayAboveSurface)
{
    MTRand rng(kSeed);

    // A large specularity error at grazing incidence drives a substantial
    // fraction of the raw perturbations below the surface.
    const OpticalPropertySet optics =
        MakeReflector(DistributionType::GAUSSIAN, 0.0, 400.0);

    const glm::dvec3 normal = glm::dvec3(0.0, 0.0, 1.0);
    const glm::dvec3 axis   = glm::normalize(glm::dvec3(1.0, 0.0, 1.0));

    for (int i = 0; i < 20000; ++i)
    {
        const glm::dvec3 out =
            ApplySpecularityError(rng, axis, optics, false, normal);

        ASSERT_GE(glm::dot(out, normal), 0.0)
            << "perturbed ray passed through the surface";
    }
}

// With the surface normal opposed to the incoming direction no perturbation
// can ever be accepted, so the retry cap is the only thing that terminates.
TEST(TracingErrors, RejectionCapTerminates)
{
    MTRand rng(kSeed);

    const OpticalPropertySet optics =
        MakeReflector(DistributionType::GAUSSIAN, 0.0, 1.0);

    const glm::dvec3 axis   = glm::dvec3(0.0, 0.0, 1.0);
    const glm::dvec3 normal = -axis;

    glm::dvec3 out = ApplySpecularityError(rng, axis, optics, false, normal);

    // The call returned, and it did so having exhausted the cap rather than
    // having found an acceptable direction.
    EXPECT_LT(glm::dot(out, normal), 0.0);
    EXPECT_NEAR(glm::length(out), 1.0, 1e-12);
}

// ---------------------------------------------------------------------------
// Error magnitudes and units
// ---------------------------------------------------------------------------

// Slope error is specified in mrad. Sampling two independent normal
// components makes the polar angle Rayleigh distributed, so its RMS is
// sigma*sqrt(2).
TEST(TracingErrors, SlopeErrorGaussianMagnitudeInMrad)
{
    const double             kSlopeMrad = 5.0;
    const OpticalPropertySet optics =
        MakeReflector(DistributionType::GAUSSIAN, kSlopeMrad, 0.0);

    const std::vector<double> angles = SampleSlopeAngles(optics, kSamples);

    EXPECT_NEAR(
        RootMeanSquare(angles), kSlopeMrad * std::sqrt(2.0), 0.02 * kSlopeMrad);
}

// Pillbox samples uniformly over a disc of radius sigma, so the RMS polar
// angle is sigma/sqrt(2).
TEST(TracingErrors, SlopeErrorPillboxMagnitudeInMrad)
{
    const double             kSlopeMrad = 5.0;
    const OpticalPropertySet optics =
        MakeReflector(DistributionType::PILLBOX, kSlopeMrad, 0.0);

    const std::vector<double> angles = SampleSlopeAngles(optics, kSamples);

    EXPECT_NEAR(
        RootMeanSquare(angles), kSlopeMrad / std::sqrt(2.0), 0.02 * kSlopeMrad);
}

TEST(TracingErrors, SpecularityErrorGaussianMagnitudeInMrad)
{
    const double             kSpecMrad = 5.0;
    const OpticalPropertySet optics =
        MakeReflector(DistributionType::GAUSSIAN, 0.0, kSpecMrad);

    const std::vector<double> angles = SampleSurfaceAngles(optics, kSamples);

    EXPECT_NEAR(
        RootMeanSquare(angles), kSpecMrad * std::sqrt(2.0), 0.02 * kSpecMrad);
}

TEST(TracingErrors, SpecularityErrorPillboxMagnitudeInMrad)
{
    const double             kSpecMrad = 5.0;
    const OpticalPropertySet optics =
        MakeReflector(DistributionType::PILLBOX, 0.0, kSpecMrad);

    const std::vector<double> angles = SampleSurfaceAngles(optics, kSamples);

    EXPECT_NEAR(
        RootMeanSquare(angles), kSpecMrad / std::sqrt(2.0), 0.02 * kSpecMrad);
}

TEST(TracingErrors, NoneDistributionIsIdentity)
{
    MTRand                   rng(kSeed);
    const OpticalPropertySet optics =
        MakeReflector(DistributionType::NONE, 0.0, 0.0);

    const glm::dvec3 axis = glm::normalize(glm::dvec3(0.2, -0.4, 0.9));

    // Compared component-wise: acos() is ill-conditioned near zero angle.
    for (int i = 0; i < 1000; ++i)
    {
        glm::dvec3 out = ApplySpecularityError(rng, axis, optics, false, axis);
        EXPECT_NEAR(out.x, axis.x, kFrameTol);
        EXPECT_NEAR(out.y, axis.y, kFrameTol);
        EXPECT_NEAR(out.z, axis.z, kFrameTol);

        out = ApplySlopeError(rng, axis, optics, false);
        EXPECT_NEAR(out.x, axis.x, kFrameTol);
        EXPECT_NEAR(out.y, axis.y, kFrameTol);
        EXPECT_NEAR(out.z, axis.z, kFrameTol);
    }
}

// ---------------------------------------------------------------------------
// Diffuse scattering
// ---------------------------------------------------------------------------

// A gray diffuse surface is Lambertian: directions are cosine-weighted over
// the hemisphere about the surface normal, giving p(theta) = 2*sin*cos,
// E[cos(theta)] = 2/3 and CDF F(theta) = sin(theta)^2.
TEST(TracingErrors, DiffuseIsLambertian)
{
    MTRand                   rng(kSeed);
    const OpticalPropertySet optics =
        MakeReflector(DistributionType::DIFFUSE, 0.0, 0.0);

    const glm::dvec3 normal = glm::normalize(glm::dvec3(0.2, -0.4, 0.9));

    std::vector<double> thetas;
    std::vector<double> cosines;
    thetas.reserve(kSamples);
    cosines.reserve(kSamples);

    double max_theta = 0.0;
    for (int i = 0; i < kSamples; ++i)
    {
        const glm::dvec3 out =
            ApplySpecularityError(rng, normal, optics, false, normal);

        const double theta = AngleBetween(normal, out);
        thetas.push_back(theta);
        cosines.push_back(std::cos(theta));
        max_theta = std::max(max_theta, theta);
    }

    // Scattering covers the full hemisphere, not a narrow lobe.
    EXPECT_GT(max_theta, 1.5);

    EXPECT_NEAR(Mean(cosines), 2.0 / 3.0, 0.01);

    const double d = KsStatistic(
        thetas, [](double t) { return std::pow(std::sin(t), 2.0); });
    EXPECT_LT(d, KsThreshold(kSamples));
}

TEST(TracingErrors, DiffuseAzimuthIsUniform)
{
    MTRand                   rng(kSeed);
    const OpticalPropertySet optics =
        MakeReflector(DistributionType::DIFFUSE, 0.0, 0.0);

    const glm::dvec3 normal(0.0, 0.0, 1.0);

    std::vector<double> cos_phi;
    std::vector<double> sin_phi;
    cos_phi.reserve(kSamples);
    sin_phi.reserve(kSamples);

    for (int i = 0; i < kSamples; ++i)
    {
        const glm::dvec3 out =
            ApplySpecularityError(rng, normal, optics, false, normal);

        const double r = std::hypot(out.x, out.y);
        if (r < 1e-15) continue;
        cos_phi.push_back(out.x / r);
        sin_phi.push_back(out.y / r);
    }

    EXPECT_NEAR(Mean(cos_phi), 0.0, 0.01);
    EXPECT_NEAR(Mean(sin_phi), 0.0, 0.01);
}

// Lambertian scattering is independent of the incoming direction, so the mean
// cosine about the plate normal must be 2/3 no matter how the plate is tilted.
TEST(TracingErrors, DiffuseIsIndependentOfIncidenceDirection)
{
    const uint_fast64_t NRAYS = 20000;

    auto mean_cosine_for_tilt = [&](const glm::dvec3& aim)
    {
        using SolTrace::Runner::RunnerStatus;

        SolTrace::NativeRunner::NativeRunner runner;
        EXPECT_EQ(runner.initialize(), RunnerStatus::SUCCESS);

        SimulationData sd;

        auto sun = make_ray_source<Sun>();
        sun->set_position(0, 0, 100);
        sd.add_ray_source(sun);

        auto stage = make_stage(0);
        stage->set_origin(0, 0, 0);
        stage->set_aim_vector(0, 0, 1);
        stage->set_name("stage");

        auto plate = make_element<SingleElement>();
        plate->set_origin(0, 0, 0);
        plate->set_aim_vector(aim.x, aim.y, aim.z);
        plate->set_surface(make_surface<Flat>());
        plate->set_aperture(make_aperture<Rectangle>(20, 20));
        plate->set_name("plate");

        OpticalPropertySet plate_optics =
            MakeReflector(DistributionType::DIFFUSE, 0.0, 0.0, "diffuse_plate");
        plate->set_optical_property_set(
            sd.add_optical_property_set(plate_optics));

        stage->add_element(plate);
        sd.add_stage(stage);

        SimulationParameters& params    = sd.get_simulation_parameters();
        params.number_of_rays           = NRAYS;
        params.max_number_of_rays       = NRAYS * 100;
        params.include_optical_errors   = true;
        params.include_sun_shape_errors = false;
        params.seed                     = kSeed;

        EXPECT_EQ(runner.setup_simulation(&sd), RunnerStatus::SUCCESS);
        EXPECT_EQ(runner.run_simulation(), RunnerStatus::SUCCESS);

        SimulationResult result;
        EXPECT_EQ(runner.report_simulation(&result, 0), RunnerStatus::SUCCESS);

        const glm::dvec3 normal   = glm::normalize(aim);
        const element_id plate_id = plate->get_id();

        double acc   = 0.0;
        int    count = 0;

        auto it = result.get_ray_record_iterator();
        while (!result.is_at_end(it))
        {
            auto rec = *it;
            if (rec->get_number_of_interactions() > 1 &&
                rec->get_element(1) == plate_id)
            {
                glm::dvec3 u(0.0);
                rec->get_direction(1, u);
                acc += glm::dot(glm::normalize(u), normal);
                ++count;
            }
            ++it;
        }

        EXPECT_GT(count, 0);
        return acc / count;
    };

    const double flat   = mean_cosine_for_tilt(glm::dvec3(0.0, 0.0, 100.0));
    const double tilted = mean_cosine_for_tilt(glm::dvec3(0.0, 40.0, 100.0));

    EXPECT_NEAR(flat, 2.0 / 3.0, 0.02);
    EXPECT_NEAR(tilted, 2.0 / 3.0, 0.02);
    EXPECT_NEAR(flat, tilted, 0.02);
}

// ---------------------------------------------------------------------------
// Sun shapes
// ---------------------------------------------------------------------------

TEST(TracingErrors, SunShapeGaussianMatchesRayleigh)
{
    const double kSigma = 2.73;
    TSun         sun    = MakeSun(SunShape::GAUSSIAN, kSigma, 4.65);

    const std::vector<double> angles = SampleSunAngles(sun, kSamples);

    const double d = KsStatistic(
        angles,
        [&](double t)
        { return 1.0 - std::exp(-t * t / (2.0 * kSigma * kSigma)); });
    EXPECT_LT(d, KsThreshold(kSamples));
}

TEST(TracingErrors, SunShapePillboxIsUniformOnDisc)
{
    const double kHalfWidth = 4.65;
    TSun         sun = MakeSun(SunShape::PILLBOX, kHalfWidth, kHalfWidth);

    const std::vector<double> angles = SampleSunAngles(sun, kSamples);

    const double d = KsStatistic(
        angles,
        [&](double t)
        { return std::clamp(t * t / (kHalfWidth * kHalfWidth), 0.0, 1.0); });
    EXPECT_LT(d, KsThreshold(kSamples));
}

TEST(TracingErrors, SunShapeLimbDarkenedMatchesProfile)
{
    const double kMaxAngle = 4.65;
    TSun         sun       = MakeSun(SunShape::LIMBDARKENED, 0.0, kMaxAngle);

    const std::vector<double> angles = SampleSunAngles(sun, kSamples);

    const RadialCdf cdf([&](double t)
                        { return 1.0 - 0.5138 * std::pow(t / kMaxAngle, 4.0); },
                        kMaxAngle);

    const double d = KsStatistic(angles, [&](double t) { return cdf(t); });
    EXPECT_LT(d, KsThreshold(kSamples));
}

TEST(TracingErrors, SunShapeBuieCsrMatchesProfile)
{
    const double kMaxAngle = 43.6;
    const double kCsr      = 0.05;

    double kappa = 0.0, gamma = 0.0;
    BuieParameters(kCsr, kappa, gamma);

    TSun sun       = MakeSun(SunShape::BUIE_CSR, 0.0, kMaxAngle);
    sun.buie_kappa = kappa;
    sun.buie_gamma = gamma;

    const std::vector<double> angles = SampleSunAngles(sun, kSamples);

    const RadialCdf cdf(
        [&](double t)
        {
            if (std::abs(t) <= 4.65)
                return std::cos(0.326 * t) / std::cos(0.308 * t);
            return std::exp(kappa) * std::pow(std::abs(t), gamma);
        },
        kMaxAngle);

    const double d = KsStatistic(angles, [&](double t) { return cdf(t); });
    EXPECT_LT(d, KsThreshold(kSamples));
}

TEST(TracingErrors, SunShapeUserDefinedMatchesProfile)
{
    const double kMaxAngle = 10.0;

    TSun sun          = MakeSun(SunShape::USER_DEFINED, 0.0, kMaxAngle);
    sun.SunShapeAngle = { 0.0, 2.5, 5.0, 7.5, 10.0 };
    // A simple ramp-down profile, linearly interpolated by the sampler.
    sun.SunShapeIntensity = { 1.0, 0.8, 0.5, 0.2, 0.0 };
    sun.MaxIntensity      = 1.0;

    const std::vector<double> angles = SampleSunAngles(sun, kSamples);

    const RadialCdf cdf(
        [&](double t)
        {
            const auto& a = sun.SunShapeAngle;
            const auto& s = sun.SunShapeIntensity;
            size_t      i = 0;
            while (i < a.size() - 1 && a[i] < t)
                ++i;
            if (i == 0) return s[0];
            return s[i - 1] +
                   (s[i] - s[i - 1]) * (t - a[i - 1]) / (a[i] - a[i - 1]);
        },
        kMaxAngle);

    const double d = KsStatistic(angles, [&](double t) { return cdf(t); });
    EXPECT_LT(d, KsThreshold(kSamples));
}

TEST(TracingErrors, NativeRunnerRejectsUnimplementedDistribution)
{
    using SolTrace::Runner::RunnerStatus;

    SolTrace::NativeRunner::NativeRunner runner;
    SimulationData                       simulation;

    auto sun = make_ray_source<Sun>();
    sun->set_position(0.0, 0.0, 100.0);
    simulation.add_ray_source(sun);

    auto mirror = make_element<SingleElement>();
    mirror->set_aperture(make_aperture<Rectangle>(10.0, 10.0));
    mirror->set_surface(make_surface<Flat>());

    OpticalPropertySet optics(InteractionType::REFLECTION,
                              "user_defined_errors");
    optics.set_ideal_reflection(OpticalSide::Both);
    optics.set_errors(
        OpticalSide::Both, DistributionType::USER_DEFINED, 0.0, 0.0);

    mirror->set_optical_property_set(
        simulation.add_optical_property_set(optics));
    simulation.add_element(mirror);

    ASSERT_EQ(runner.initialize(), RunnerStatus::SUCCESS);
    EXPECT_THROW(runner.setup_simulation(&simulation), std::invalid_argument);
}

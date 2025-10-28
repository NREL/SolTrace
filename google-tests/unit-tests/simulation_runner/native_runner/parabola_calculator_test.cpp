#include <gtest/gtest.h>

#include <cmath>

#include <common.hpp>
#include <parabola_calculator.hpp>
#include <simulation_data_export.hpp>
#include <surface.hpp>
#include <vector3d.hpp>

using SolTrace::NativeRunner::ParabolaCalculator;

// NOTES: Equation for a parabola is z = (cx*x^2 + cy*y^2) / 2
// Computing the intersection point comes down to solving
// quadratic equation for a parameter t given by
//      at^2 + bt + c = 0
// The cases here are broken down by cases with regard to
// this equation. The case is noted at the top.  The coefficients
// a, b, and c can be determined by subtituting x(t) = x0 + mx*t,
// y(t) = y0 + my*t, and z(t) = z0 + mz*t into the given equation
// for a parabola. In the comments and naming, we use the
// following terms:
//    Delta = b^2 - 4ac  -- discrimanant of quadratic equation
//    t1 = (-b - sqrt(Delta)) / (2a) -- negative root
//    t2 = (-b + sqrt(Delta)) / (2a) -- positive root
// From these last two, we always have t1 < t2.

double focal_length(double c)
{
    return 0.5 / c;
}

// Constructor validation tests
TEST(ParabolaCalculator, ConstructorNullSurfaceThrows)
{
    auto circ = create_circle_aperture();
    EXPECT_THROW({ ParabolaCalculator calc(nullptr, circ); }, std::invalid_argument);
}

TEST(ParabolaCalculator, ConstructorNullApertureThrows)
{
    auto para = create_parabola_surface(1.0);
    EXPECT_THROW({ ParabolaCalculator calc(para, nullptr); }, std::invalid_argument);
}

TEST(ParabolaCalculator, ConstructorWrongSurfaceTypeThrows)
{
    auto flat = std::make_shared<Flat>();
    auto circ = create_circle_aperture();
    EXPECT_THROW({ ParabolaCalculator calc(flat, circ); }, std::invalid_argument);
}

TEST(ParabolaCalculator, ConstructorZeroFocalLengthAllowed)
{
    auto zero_focal_surface = create_parabola_surface(0.0);
    auto circ = create_circle_aperture();
    EXPECT_NO_THROW({
        ParabolaCalculator calc(zero_focal_surface, circ);
    });
}

TEST(ParabolaCalculator, ConstructorNegativeFocalLengthAllowed)
{
    auto negative_focal_surface = create_parabola_surface(-1.0);
    auto circ = create_circle_aperture();
    EXPECT_NO_THROW({
        ParabolaCalculator calc(negative_focal_surface, circ);
    });
}

TEST(ParabolaCalculator, ConstructorNaNFocalLengthThrows)
{
    auto nan_focal_surface = create_parabola_surface(std::nan(""));
    auto circ = create_circle_aperture();
    EXPECT_THROW({ ParabolaCalculator calc(nan_focal_surface, circ); }, std::invalid_argument);
}

TEST(ParabolaCalculator, ConstructorInfiniteFocalLengthThrows)
{
    auto inf_focal_surface = create_parabola_surface(std::numeric_limits<double>::infinity());
    auto circ = create_circle_aperture();
    EXPECT_THROW({ ParabolaCalculator calc(inf_focal_surface, circ); }, std::invalid_argument);
}

TEST(ParabolaCalculator, ConstructorValidFocalLength)
{
    auto valid_surface = create_parabola_surface(1.0);
    auto circ = create_circle_aperture();
    EXPECT_NO_THROW({
        ParabolaCalculator calc(valid_surface, circ);
    });
}

TEST(ParabolaCalculator, Case1)
{
    // Case: a == 0, t <= 0 -- returns no solution
    Vector3d zero;
    zero.zero();
    // Ray location
    Vector3d x0(0.0, 0.0, 1.0);
    // Ray direction
    Vector3d m(0.0, 0.0, 1.0);

    // Solution values
    double t;
    Vector3d xt;
    Vector3d mt;
    Vector3d gradf;

    auto parabola = SolTrace::Data::make_surface<Parabola>(0.5, 0.25);
    auto circ = create_circle_aperture(10.0);
    ParabolaCalculator pcalc(parabola, circ);
    int sts = pcalc.intersect(x0.data, m.data,
                              xt.data, mt.data, gradf.data, &t);
    EXPECT_EQ(sts, 1);
    EXPECT_EQ(t, 0.0);
    EXPECT_TRUE(is_identical(xt, zero));
    EXPECT_TRUE(is_identical(mt, zero));
    EXPECT_TRUE(is_identical(gradf, zero));
}

TEST(ParabolaCalculator, Case2)
{
    // Case: a == 0, t > 0 -- returns t
    // NOTE: Here the quadratic equation reduces to a linear equation
    // Ray location
    Vector3d x0(0.0, 0.0, -1.0);
    // Ray direction
    Vector3d m(0.0, 0.0, 1.0);
    // Intersection point
    const double T = 1.0;
    const double TOL = 1e-12;
    // Parabola constants
    const double cx = 1.0;
    const double cy = 2.0;

    // Solution values
    double t;
    Vector3d xt;
    Vector3d mt;
    Vector3d gradf;

    auto parabola = SolTrace::Data::make_surface<Parabola>(0.5, 0.25);
    auto circ = create_circle_aperture(10.0);
    ParabolaCalculator pcalc(parabola, circ);
    int sts = pcalc.intersect(x0.data, m.data,
                              xt.data, mt.data, gradf.data, &t);
    EXPECT_EQ(sts, 0);
    EXPECT_NEAR(t, T, TOL);
    EXPECT_NEAR(xt[0], 0.0, TOL);
    EXPECT_NEAR(xt[1], 0.0, TOL);
    EXPECT_NEAR(xt[2], 0.0, TOL);
    EXPECT_TRUE(is_identical(mt, m));
    EXPECT_NEAR(gradf[0], 0.0, TOL);
    EXPECT_NEAR(gradf[1], 0.0, TOL);
    EXPECT_NEAR(gradf[2], 1.0, TOL);
    EXPECT_NEAR(xt[2] - 0.5 * (cx * xt[0] * xt[0] + cy * xt[1] * xt[1]), 0.0, TOL);
}

TEST(ParabolaCalculator, Case3)
{
    // Case: a != 0, Delta = b^2 - 4ac < 0 -- returns no solution
    Vector3d zero;
    zero.zero();
    // Ray location
    Vector3d x0(-2.0, 1.0, 0.0);
    // Ray direction
    Vector3d m(1.0, 1.0, 1.0);

    // Solution values
    double t;
    Vector3d xt;
    Vector3d mt;
    Vector3d gradf;

    auto parabola = SolTrace::Data::make_surface<Parabola>(0.5, 0.25);
    auto circ = create_circle_aperture(10.0);
    ParabolaCalculator pcalc(parabola, circ);
    int sts = pcalc.intersect(x0.data, m.data,
                              xt.data, mt.data, gradf.data, &t);
    EXPECT_EQ(sts, 1);
    EXPECT_EQ(t, 0.0);
    EXPECT_TRUE(is_identical(xt, zero));
    EXPECT_TRUE(is_identical(mt, zero));
    EXPECT_TRUE(is_identical(gradf, zero));
}

TEST(ParabolaCalculator, Case4)
{
    // Case: a != 0, Delta >= 0, t1 > 0 -- returns t1
    // Ray location
    Vector3d x0(-2.0, 1.0, 0.0);
    // Ray direction
    Vector3d m(2.0, 0.5, 3.0);
    // Intersection point
    const double T = 2.0 / 3.0;
    const double TOL = 1e-12;
    // Parabola constants
    const double cx = 1.0;
    const double cy = 2.0;

    // Solution values
    double t;
    Vector3d xt;
    Vector3d mt;
    Vector3d gradf;

    auto parabola = SolTrace::Data::make_surface<Parabola>(focal_length(cx),
                                                           focal_length(cy));
    auto circ = create_circle_aperture(10.0);
    ParabolaCalculator pcalc(parabola, circ);
    int sts = pcalc.intersect(x0.data, m.data,
                              xt.data, mt.data, gradf.data, &t);
    EXPECT_EQ(sts, 0);
    EXPECT_NEAR(t, T, TOL);
    EXPECT_NEAR(xt[0], m[0] * T + x0[0], TOL);
    EXPECT_NEAR(xt[1], m[1] * T + x0[1], TOL);
    EXPECT_NEAR(xt[2], m[2] * T + x0[2], TOL);
    EXPECT_TRUE(is_identical(mt, m));
    EXPECT_NEAR(gradf[0], -cx * (m[0] * T + x0[0]), TOL);
    EXPECT_NEAR(gradf[1], -cy * (m[1] * T + x0[1]), TOL);
    EXPECT_NEAR(gradf[2], 1.0, TOL);
    EXPECT_NEAR(xt[2] - 0.5 * (cx * xt[0] * xt[0] + cy * xt[1] * xt[1]), 0.0, TOL);
}

TEST(ParabolaCalculator, Case5)
{
    // Case: a != 0, Delta >= 0, t1 < 0, t2 > 0 -- returns t2
    // Ray location
    Vector3d x0(1.0, -1.0, 3.0);
    // Ray direction
    Vector3d m(1.0, 0.0, 2.0);
    // Intersection point
    const double T = 3.0;
    const double TOL = 1e-12;
    // Parabola constants
    const double cx = 1.0;
    const double cy = 2.0;

    // Solution values
    double t;
    Vector3d xt;
    Vector3d mt;
    Vector3d gradf;

    auto parabola = SolTrace::Data::make_surface<Parabola>(focal_length(cx),
                                                           focal_length(cy));
    auto circ = create_circle_aperture(10.0);
    ParabolaCalculator pcalc(parabola, circ);
    int sts = pcalc.intersect(x0.data, m.data,
                              xt.data, mt.data, gradf.data, &t);
    EXPECT_EQ(sts, 0);
    EXPECT_NEAR(t, T, TOL);
    EXPECT_NEAR(xt[0], m[0] * T + x0[0], TOL);
    EXPECT_NEAR(xt[1], m[1] * T + x0[1], TOL);
    EXPECT_NEAR(xt[2], m[2] * T + x0[2], TOL);
    EXPECT_TRUE(is_identical(mt, m));
    EXPECT_NEAR(gradf[0], -cx * (m[0] * T + x0[0]), TOL);
    EXPECT_NEAR(gradf[1], -cy * (m[1] * T + x0[1]), TOL);
    EXPECT_NEAR(gradf[2], 1.0, TOL);
    EXPECT_NEAR(xt[2] - 0.5 * (cx * xt[0] * xt[0] + cy * xt[1] * xt[1]), 0.0, TOL);
}

TEST(ParabolaCalculator, Case6)
{
    // Case: a != 0, Delta >= 0, t1 < 0, t2 < 0 -- returns no solution
    Vector3d zero;
    zero.zero();
    // Ray location
    Vector3d x0(3.0, 1.0, 2.0);
    // Ray direction
    Vector3d m(1.0, -1.0, -4.0);
    // Intersection point
    const double T = 1.0 + sqrt(8.0);
    const double TOL = 1e-12;
    // Parabola constants
    const double cx = 1.0;
    const double cy = 2.0;

    // Solution values
    double t;
    Vector3d xt;
    Vector3d mt;
    Vector3d gradf;

    auto parabola = SolTrace::Data::make_surface<Parabola>(focal_length(cx),
                                                           focal_length(cy));
    auto circ = create_circle_aperture(10.0);
    ParabolaCalculator pcalc(parabola, circ);
    int sts = pcalc.intersect(x0.data, m.data,
                              xt.data, mt.data, gradf.data, &t);
    EXPECT_EQ(sts, 1);
    EXPECT_EQ(t, 0.0);
    EXPECT_TRUE(is_identical(xt, zero));
    EXPECT_TRUE(is_identical(mt, zero));
    EXPECT_TRUE(is_identical(gradf, zero));
}

TEST(ParabolaCalculator, ZAperture)
{
    const double TOL = 1e-12;

    double cx = 1.0;
    double cy = 2.0;
    double r = 2.5;
    auto ap = SolTrace::Data::make_aperture<Circle>(2.0 * r);
    auto parabola = SolTrace::Data::make_surface<Parabola>(focal_length(cx),
                                                           focal_length(cy));
    auto circ = create_circle_aperture(10.0);
    ParabolaCalculator pcalc(parabola, circ);
    double zap = pcalc.compute_z_aperture(ap);
    double zmax = 0.5 * cy * r * r;
    EXPECT_NEAR(zap, zmax, TOL);

    // Swap x and y coefficients
    parabola = SolTrace::Data::make_surface<Parabola>(focal_length(cy),
                                                      focal_length(cx));
    ParabolaCalculator pcalc_swap(parabola, circ);
    zap = pcalc.compute_z_aperture(ap);
    EXPECT_NEAR(zap, zmax, TOL);
}

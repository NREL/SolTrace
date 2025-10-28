#include <gtest/gtest.h>

#include <cmath>

#include <common.hpp>
#include <sphere_calculator.hpp>
#include <simulation_data_export.hpp>

using SolTrace::NativeRunner::SphereCalculator;

// NOTES: Equation for a sphere is x^2 + y^2 + (z - r)^2 = r^2
// where r is the radius. Computing the intersection point 
// comes down to solving quadratic equation for a parameter t 
// given by
//      at^2 + bt + c = 0
// The cases here are broken down by cases with regard to
// this equation. The case is noted at the top.  The coefficients
// a, b, and c can be determined by subtituting x(t) = x0 + mx*t,
// y(t) = y0 + my*t, and z(t) = z0 + mz*t into the given equation
// for a sphere. In the comments and naming, we use the
// following terms:
//    Delta = b^2 - 4ac  -- discrimanant of quadratic equation
//    t1 = (-b - sqrt(Delta)) / (2a) -- negative root
//    t2 = (-b + sqrt(Delta)) / (2a) -- positive root
//    z1 = z(t1)
//    z2 = z(t2)
// From these last two, we always have t1 < t2.

// Constructor validation tests
TEST(SphereCalculator, ConstructorNullSurfaceThrows)
{
    EXPECT_THROW({
        SphereCalculator calc(nullptr, nullptr);
    }, std::invalid_argument);
}

TEST(SphereCalculator, ConstructorWrongSurfaceTypeThrows)
{
    auto flat = std::make_shared<Flat>();
    auto circ = create_circle_aperture();
    EXPECT_THROW({
        SphereCalculator calc(flat, circ);
    }, std::invalid_argument);
}

TEST(SphereCalculator, ConstructorZeroVertexCurvatureThrows)
{
    auto zero_curvature_surface = create_sphere_surface(0.0);
    auto circ = create_circle_aperture();
    EXPECT_THROW({
        SphereCalculator calc(zero_curvature_surface, circ);
    }, std::invalid_argument);
}

TEST(SphereCalculator, ConstructorNegativeVertexCurvatureThrows)
{
    auto negative_curvature_surface = create_sphere_surface(-0.1);
    auto circ = create_circle_aperture();
    EXPECT_THROW({
        SphereCalculator calc(negative_curvature_surface, circ);
    }, std::invalid_argument);
}

TEST(SphereCalculator, ConstructorNaNVertexCurvatureThrows)
{
    auto nan_curvature_surface = create_sphere_surface(std::nan(""));
    auto circ = create_circle_aperture();
    EXPECT_THROW({
        SphereCalculator calc(nan_curvature_surface, circ);
    }, std::invalid_argument);
}

TEST(SphereCalculator, ConstructorInfiniteVertexCurvatureThrows)
{
    auto inf_curvature_surface = create_sphere_surface(std::numeric_limits<double>::infinity());
    auto circ = create_circle_aperture();
    EXPECT_THROW({
        SphereCalculator calc(inf_curvature_surface, circ);
    }, std::invalid_argument);
}

TEST(SphereCalculator, ConstructorNullApertureThrows)
{
    auto valid_surface = create_sphere_surface(0.1);
    EXPECT_THROW({
        SphereCalculator calc(valid_surface, nullptr);
    }, std::invalid_argument);
}

TEST(SphereCalculator, ConstructorValid)
{
    auto valid_surface = create_sphere_surface(0.1);
    auto circ = create_circle_aperture();
    EXPECT_NO_THROW({
        SphereCalculator calc(valid_surface, circ);
    });
}

TEST(SphereCalculator, Case1)
{
    // Delta < 0 -- returns no solution
    Vector3d r0(3.0, -1.0, 1.0);
    Vector3d rd(-3.0, 1.0, -2.0);
    double radius = 1.0;

    Vector3d zero;
    zero.zero();

    // Solution values
    int sts;
    double t;
    Vector3d xt;
    Vector3d mt;
    Vector3d gradf;
    
    surface_ptr sph = make_surface<Sphere>(1.0 / radius);
    auto circ = create_circle_aperture(radius);
    SphereCalculator scalc(sph, circ);
    sts = scalc.intersect(r0.data, rd.data, xt.data, mt.data, gradf.data, &t);
    EXPECT_EQ(sts, 1);
    EXPECT_EQ(t, 0.0);
    EXPECT_TRUE(is_identical(xt, zero));
    EXPECT_TRUE(is_identical(mt, zero));
    EXPECT_TRUE(is_identical(gradf, zero));
}

TEST(SphereCalculator, Case2)
{
    // t1 > 0 -- returns t1
    const double TOL = 1e-12;
    
    Vector3d r0(3.0, -1.0, 1.0);
    Vector3d rd(-3.0, 1.0, -0.5);
    double radius = 1.0;
    double T = 0.5 * (20 - sqrt(31)) / 10.25;
    Vector3d s0(0.0, 0.0, radius);

    // Solution values
    int sts;
    double t;
    Vector3d xt;
    Vector3d mt;
    Vector3d gradf;
    Vector3d r;
    
    surface_ptr sph = make_surface<Sphere>(1.0 / radius);
    auto circ = create_circle_aperture(radius);
    SphereCalculator scalc(sph, circ);
    sts = scalc.intersect(r0.data, rd.data, xt.data, mt.data, gradf.data, &t);

    EXPECT_EQ(sts, 0);
    EXPECT_NEAR(t, T, TOL);
    EXPECT_NEAR(xt[0], rd[0] * T + r0[0], TOL);
    EXPECT_NEAR(xt[1], rd[1] * T + r0[1], TOL);
    EXPECT_NEAR(xt[2], rd[2] * T + r0[2], TOL);
    EXPECT_TRUE(is_identical(mt, rd));
    EXPECT_NEAR(gradf[0], -2.0 * (rd[0] * T + r0[0]), TOL);
    EXPECT_NEAR(gradf[1], -2.0 * (rd[1] * T + r0[1]), TOL);
    EXPECT_NEAR(gradf[2], -2.0 * (rd[2] * T + r0[2] - radius), TOL);

    vector_add(1.0, xt, -1.0, s0, r);
    EXPECT_NEAR(vector_norm(r), radius, TOL);
}

TEST(SphereCalculator, Case3)
{
    // t1 > 0, z1 > r, t2 > 0, z2 < r -- returns t2
    const double TOL = 1e-12;
    
    Vector3d r0(3.0, -1.0, 2.0);
    Vector3d rd(-3.0, 1.0, -1.0);
    double radius = 1.0;
    double T = 0.5 * (22 + sqrt(44)) / 11;
    Vector3d s0(0.0, 0.0, radius);

    // Solution values
    int sts;
    double t;
    Vector3d xt;
    Vector3d mt;
    Vector3d gradf;
    Vector3d r;
    
    surface_ptr sph = make_surface<Sphere>(1.0 / radius);
    auto circ = create_circle_aperture(radius);
    SphereCalculator scalc(sph, circ);
    sts = scalc.intersect(r0.data, rd.data, xt.data, mt.data, gradf.data, &t);

    EXPECT_EQ(sts, 0);
    EXPECT_NEAR(t, T, TOL);
    EXPECT_NEAR(xt[0], rd[0] * T + r0[0], TOL);
    EXPECT_NEAR(xt[1], rd[1] * T + r0[1], TOL);
    EXPECT_NEAR(xt[2], rd[2] * T + r0[2], TOL);
    EXPECT_TRUE(is_identical(mt, rd));
    EXPECT_NEAR(gradf[0], -2.0 * (rd[0] * T + r0[0]), TOL);
    EXPECT_NEAR(gradf[1], -2.0 * (rd[1] * T + r0[1]), TOL);
    EXPECT_NEAR(gradf[2], -2.0 * (rd[2] * T + r0[2] - radius), TOL);

    vector_add(1.0, xt, -1.0, s0, r);
    EXPECT_NEAR(vector_norm(r), radius, TOL);
}

TEST(SphereCalculator, Case4)
{
    // t1, t2 > 0, z1, z2 > r -- returns no solution
    Vector3d r0(3.0, -1.0, 2.0);
    Vector3d rd(-3.0, 1.0, -0.5);
    double radius = 1.0;

    Vector3d zero;
    zero.zero();

    // Solution values
    int sts;
    double t;
    Vector3d xt;
    Vector3d mt;
    Vector3d gradf;
    
    surface_ptr sph = make_surface<Sphere>(1.0 / radius);
    auto circ = create_circle_aperture(radius);
    SphereCalculator scalc(sph, circ);
    sts = scalc.intersect(r0.data, rd.data, xt.data, mt.data, gradf.data, &t);
    EXPECT_EQ(sts, 1);
    EXPECT_EQ(t, 0.0);
    EXPECT_TRUE(is_identical(xt, zero));
    EXPECT_TRUE(is_identical(mt, zero));
    EXPECT_TRUE(is_identical(gradf, zero));
}

TEST(SphereCalculator, Case5)
{
    // t1 < 0, t2 > 0, z2 < r -- returns t2
    const double TOL = 1e-12;
    
    Vector3d r0(-0.5, 0.5, 0.5);
    Vector3d rd(-3.0, 1.0, -0.5);
    double radius = 1.0;
    double T = 0.5 * (-4.5 + sqrt(30.5)) / 10.25;

    // Solution values
    int sts;
    double t;
    Vector3d xt;
    Vector3d mt;
    Vector3d gradf;
    
    surface_ptr sph = make_surface<Sphere>(1.0 / radius);
    auto circ = create_circle_aperture(radius);
    SphereCalculator scalc(sph, circ);
    sts = scalc.intersect(r0.data, rd.data, xt.data, mt.data, gradf.data, &t);

    EXPECT_EQ(sts, 0);
    EXPECT_NEAR(t, T, TOL);
    EXPECT_NEAR(xt[0], rd[0] * T + r0[0], TOL);
    EXPECT_NEAR(xt[1], rd[1] * T + r0[1], TOL);
    EXPECT_NEAR(xt[2], rd[2] * T + r0[2], TOL);
    EXPECT_TRUE(is_identical(mt, rd));
    EXPECT_NEAR(gradf[0], -2.0 * (rd[0] * T + r0[0]), TOL);
    EXPECT_NEAR(gradf[1], -2.0 * (rd[1] * T + r0[1]), TOL);
    EXPECT_NEAR(gradf[2], -2.0 * (rd[2] * T + r0[2] - radius), TOL);
}

TEST(SphereCalculator, Case6)
{
    // t1 < 0, t2 > 0, z2 > r -- returns no intersection
    Vector3d r0(-0.5, 0.5, 0.5);
    Vector3d rd(-1.0, 1.0, 5.0);
    double radius = 1.0;

    Vector3d zero;
    zero.zero();

    // Solution values
    int sts;
    double t;
    Vector3d xt;
    Vector3d mt;
    Vector3d gradf;
    
    surface_ptr sph = make_surface<Sphere>(1.0 / radius);
    auto circ = create_circle_aperture(radius);
    SphereCalculator scalc(sph, circ);
    sts = scalc.intersect(r0.data, rd.data, xt.data, mt.data, gradf.data, &t);
    EXPECT_EQ(sts, 1);
    EXPECT_EQ(t, 0.0);
    EXPECT_TRUE(is_identical(xt, zero));
    EXPECT_TRUE(is_identical(mt, zero));
    EXPECT_TRUE(is_identical(gradf, zero));
}

TEST(SphereCalculator, Case7)
{
    // t1 < 0, t2 < 0 -- returns no intersection
    Vector3d r0(-0.5, 0.5, 3.0);
    Vector3d rd(-1.0, 1.0, 5.0);
    double radius = 1.0;

    Vector3d zero;
    zero.zero();

    // Solution values
    int sts;
    double t;
    Vector3d xt;
    Vector3d mt;
    Vector3d gradf;
    
    surface_ptr sph = make_surface<Sphere>(1.0 / radius);
    auto circ = create_circle_aperture(radius);
    SphereCalculator scalc(sph, circ);
    sts = scalc.intersect(r0.data, rd.data, xt.data, mt.data, gradf.data, &t);
    EXPECT_EQ(sts, 1);
    EXPECT_EQ(t, 0.0);
    EXPECT_TRUE(is_identical(xt, zero));
    EXPECT_TRUE(is_identical(mt, zero));
    EXPECT_TRUE(is_identical(gradf, zero));
}

TEST(SphereCalculator, Case8)
{
    // t1 == t2 > 0 -- returns t1
    const double TOL = 1e-12;
    
    Vector3d r0(-1.0, -1.0, 0.0);
    Vector3d rd(1.0, 1.0, 0.0);
    double radius = 1.0;
    double T = 0.5 * (4.0 - sqrt(0)) / 2;

    // Solution values
    int sts;
    double t;
    Vector3d xt;
    Vector3d mt;
    Vector3d gradf;
    
    surface_ptr sph = make_surface<Sphere>(1.0 / radius);
    auto circ = create_circle_aperture(radius);
    SphereCalculator scalc(sph, circ);
    sts = scalc.intersect(r0.data, rd.data, xt.data, mt.data, gradf.data, &t);

    EXPECT_EQ(sts, 0);
    EXPECT_NEAR(t, T, TOL);
    EXPECT_NEAR(xt[0], rd[0] * T + r0[0], TOL);
    EXPECT_NEAR(xt[1], rd[1] * T + r0[1], TOL);
    EXPECT_NEAR(xt[2], rd[2] * T + r0[2], TOL);
    EXPECT_TRUE(is_identical(mt, rd));
    EXPECT_NEAR(gradf[0], -2.0 * (rd[0] * T + r0[0]), TOL);
    EXPECT_NEAR(gradf[1], -2.0 * (rd[1] * T + r0[1]), TOL);
    EXPECT_NEAR(gradf[2], -2.0 * (rd[2] * T + r0[2] - radius), TOL);
}

TEST(SphereCalculator, ZAperture)
{
    const double TOL = 1e-12;

    double r = 1.0;
    surface_ptr sph = make_surface<Sphere>(1.0 / r);
    auto circ = create_circle_aperture(r);
    SphereCalculator scalc(sph, circ);

    double XLEN = 1.0;
    double YLEN = 1.0;
    aperture_ptr ap1 = make_aperture<Rectangle>(XLEN, YLEN);
    double RMAXSQ = 0.25 * (XLEN * XLEN + YLEN * YLEN);
    double zmax1 = scalc.compute_z_aperture(ap1);
    EXPECT_NEAR(zmax1, r - sqrt(r * r - RMAXSQ), TOL);

    aperture_ptr ap2 = make_aperture<Rectangle>(5.0, 5.0);
    double zmax2 = scalc.compute_z_aperture(ap2);
    EXPECT_NEAR(zmax2, r, TOL);
    
}

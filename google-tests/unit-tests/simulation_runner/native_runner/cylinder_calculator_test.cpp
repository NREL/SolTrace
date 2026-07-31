#include <gtest/gtest.h>
#include <memory>
#include <cmath>

#include <aperture.hpp>
#include <cylinder_calculator.hpp>
#include <surface.hpp>

#include "common.hpp"

using SolTrace::NativeRunner::CylinderCalculator;

// NOTES: Equation for a cylinder is x^2 + (z - r)^2 = r^2
// Computing the intersection point comes down to solving
// quadratic equation for a parameter t given by
//      at^2 + bt + c = 0
// The cases here are broken down by cases with regard to
// this equation. The case is noted at the top.  The coefficients
// a, b, and c can be determined by subtituting x(t) = x0 + mx*t,
// y(t) = y0 + my*t, and z(t) = z0 + mz*t into the given equation
// for a cylinder. In the comments and naming, we use the
// following terms:
//    Delta = b^2 - 4ac  -- discrimanant of quadratic equation
//    t1 = (-b - sqrt(Delta)) / (2a) -- negative root
//    t2 = (-b + sqrt(Delta)) / (2a) -- positive root
// From these last two, we always have t1 < t2.

// Intersection case tests
TEST(CylinderCalculator, Case1)
{
    // Case a == 0, no solution
    const double r = 1.0;
    const double ymax = 4.0;
    auto surface = create_cylinder_surface(r);
    auto aperture = create_rectangle_aperture(2.0 * r, 2.0 * ymax);
    CylinderCalculator calc(surface, aperture);

    // Ray position and direction
    glm::dvec3 x0(-1.0, -1.0, -1.0);
    glm::dvec3 m(0.0, 1.0, 0.0);

    // Solution values
    double t;
    glm::dvec3 xt;
    glm::dvec3 mt;
    glm::dvec3 gradf;

    int result = calc.intersect(x0, m, xt, mt, gradf, &t);

    EXPECT_EQ(result, 1);
    EXPECT_EQ(t, 0.0);
}

TEST(CylinderCalculator, Case2)
{
    // Case a != 0, Delta < 0 -- no solution
    const double r = 1.0;
    const double ymax = 4.0;
    auto surface = create_cylinder_surface(r);
    auto aperture = create_rectangle_aperture(2.0 * r, 2.0 * ymax);
    CylinderCalculator calc(surface, aperture);

    // Ray position and direction
    glm::dvec3 x0(5.0, 0.0, 1.0);
    glm::dvec3 m(-1.0, 1.0, 1.0);

    // Solution values
    double t;
    glm::dvec3 xt;
    glm::dvec3 mt;
    glm::dvec3 gradf;

    int result = calc.intersect(x0, m, xt, mt, gradf, &t);

    EXPECT_EQ(result, 1);
    EXPECT_EQ(t, 0.0);
}

TEST(CylinderCalculator, Case3)
{
    // Case a != 0, Delta > 0, t1 > 0, t2 > 0 -- returns t1
    const double TOL = 1e-12;
    const double r = 1.0;
    const double ymax = 4.0;
    auto surface = create_cylinder_surface(r);
    auto aperture = create_rectangle_aperture(2.0 * r, 2.0 * ymax);
    CylinderCalculator calc(surface, aperture);

    // Ray position and direction
    glm::dvec3 x0(5.0, -3.0, 0.0);
    glm::dvec3 m(-1.0, 1.0, 0.0);

    // Solution values
    double t;
    glm::dvec3 xt;
    glm::dvec3 mt;
    glm::dvec3 gradf;

    int result = calc.intersect(x0, m, xt, mt, gradf, &t);

    EXPECT_EQ(result, 0);
    EXPECT_NEAR(t, 4.0, TOL);
    EXPECT_NEAR(xt[0], 1.0, TOL);
    EXPECT_NEAR(xt[1], 1.0, TOL);
    EXPECT_NEAR(xt[2], 0.0, TOL);
    EXPECT_NEAR(mt[0], m[0], TOL);
    EXPECT_NEAR(mt[1], m[1], TOL);
    EXPECT_NEAR(mt[2], m[2], TOL);
}

TEST(CylinderCalculator, Case4)
{
    // Case a != 0, Delta > 0, t1 > 0, t2 > 0 -- returns t1
    const double TOL = 1e-12;
    const double r = 1.0;
    const double ymax = 4.0;
    auto surface = create_cylinder_surface(r);
    auto aperture = create_rectangle_aperture(2.0 * r, 2.0 * ymax);
    CylinderCalculator calc(surface, aperture);

    // Ray position and direction
    glm::dvec3 x0(0.0, -1.0, 0.0);
    glm::dvec3 m(-1.0, 1.0, 0.0);

    // Solution values
    double t;
    glm::dvec3 xt;
    glm::dvec3 mt;
    glm::dvec3 gradf;

    int result = calc.intersect(x0, m, xt, mt, gradf, &t);

    EXPECT_EQ(result, 0);
    EXPECT_NEAR(t, 1.0, TOL);
    EXPECT_NEAR(xt[0], -1.0, TOL);
    EXPECT_NEAR(xt[1], 0.0, TOL);
    EXPECT_NEAR(xt[2], 0.0, TOL);
    EXPECT_NEAR(mt[0], m[0], TOL);
    EXPECT_NEAR(mt[1], m[1], TOL);
    EXPECT_NEAR(mt[2], m[2], TOL);
}

TEST(CylinderCalculator, Case5)
{
    // Case a != 0, Delta > 0, t1 > 0, t2 > 0, yt1 > ymax, yt2 > ymax -- no solution
    const double r = 1.0;
    const double ymax = 4.0;
    auto surface = create_cylinder_surface(r);
    auto aperture = create_rectangle_aperture(2.0 * r, 2.0 * ymax);
    CylinderCalculator calc(surface, aperture);

    // Ray position and direction
    glm::dvec3 x0(5.0, ymax, 1.0);
    glm::dvec3 m(-1.0, 1.0, 0.0);

    // Solution values
    double t;
    glm::dvec3 xt;
    glm::dvec3 mt;
    glm::dvec3 gradf;

    int result = calc.intersect(x0, m, xt, mt, gradf, &t);

    EXPECT_EQ(result, 1);
    EXPECT_EQ(t, 0.0);
}

TEST(CylinderCalculator, Case6)
{
    // Case a != 0, Delta > 0, t1 < 0, t2 < 0
    const double r = 1.0;
    const double ymax = 4.0;
    auto surface = create_cylinder_surface(r);
    auto aperture = create_rectangle_aperture(2.0 * r, 2.0 * ymax);
    CylinderCalculator calc(surface, aperture);

    // Ray position and direction
    glm::dvec3 x0(5.0, ymax, 1.0);
    glm::dvec3 m(-1.0, 1.0, 0.0);

    // Solution values
    double t;
    glm::dvec3 xt;
    glm::dvec3 mt;
    glm::dvec3 gradf;

    int result = calc.intersect(x0, m, xt, mt, gradf, &t);

    EXPECT_EQ(result, 1);
    EXPECT_EQ(t, 0.0);
}

// Constructor validation tests -- includes error tests
TEST(CylinderCalculator, ConstructorValidConstruction)
{
    auto surface = create_cylinder_surface();
    auto aperture = create_rectangle_aperture();
    EXPECT_NO_THROW({
        CylinderCalculator calc(surface, aperture);
    });
}

TEST(CylinderCalculator, ConstructorNullSurfaceThrows)
{
    auto aperture = create_rectangle_aperture();
    EXPECT_THROW({ CylinderCalculator calc(nullptr, aperture); }, std::invalid_argument);
}

TEST(CylinderCalculator, ConstructorNullApertureThrows)
{
    auto surface = create_cylinder_surface();
    EXPECT_THROW({ CylinderCalculator calc(surface, nullptr); }, std::invalid_argument);
}

TEST(CylinderCalculator, ConstructorWrongSurfaceTypeThrows)
{
    auto flat = std::make_shared<Flat>();
    auto aperture = create_rectangle_aperture();
    EXPECT_THROW({ CylinderCalculator calc(flat, aperture); }, std::invalid_argument);
}

TEST(CylinderCalculator, ConstructorWrongApertureTypeThrows)
{
    auto surface = create_cylinder_surface();
    auto circular_ap = create_circle_aperture();
    EXPECT_THROW({ CylinderCalculator calc(surface, circular_ap); }, std::invalid_argument);
}

TEST(CylinderCalculator, ConstructorZeroRadiusThrows)
{
    auto surf = create_cylinder_surface();
    surf->radius = 0.0;
    auto aperture = create_rectangle_aperture();
    EXPECT_THROW({ CylinderCalculator calc(surf, aperture); }, std::invalid_argument);
}

TEST(CylinderCalculator, ConstructorNegativeRadiusThrows)
{
    auto surf = create_cylinder_surface();
    surf->radius = -1.0;
    auto aperture = create_rectangle_aperture();
    EXPECT_THROW({ CylinderCalculator calc(surf, aperture); }, std::invalid_argument);
}

TEST(CylinderCalculator, ConstructorNaNRadiusThrows)
{
    auto surf = create_cylinder_surface();
    surf->radius = std::nan("");
    auto aperture = create_rectangle_aperture();
    EXPECT_THROW({ CylinderCalculator calc(surf, aperture); }, std::invalid_argument);
}

TEST(CylinderCalculator, ConstructorInfiniteRadiusThrows)
{
    auto surf = create_cylinder_surface();
    surf->radius = std::numeric_limits<double>::infinity();
    auto aperture = create_rectangle_aperture();
    EXPECT_THROW({ CylinderCalculator calc(surf, aperture); }, std::invalid_argument);
}

TEST(CylinderCalculator, ConstructorMismatchedApertureDimensionsThrows)
{
    auto surface = create_cylinder_surface(1.0); // radius = 1.0, diameter = 2.0
    // Create aperture with x_length != 2 * radius
    auto mismatched_aperture = create_rectangle_aperture(3.0, 2.0); // x_length = 3.0, but cylinder diameter = 2.0
    EXPECT_THROW({ CylinderCalculator calc(surface, mismatched_aperture); }, std::invalid_argument);
}

TEST(CylinderCalculator, ConstructorZeroApertureDimensionsThrows)
{
    auto surface = create_cylinder_surface();
    auto aperture = create_rectangle_aperture();
    aperture->set_y_length(0.0);
    EXPECT_THROW({ CylinderCalculator calc(surface, aperture); }, std::invalid_argument);
}

// // Basic intersection test
// TEST(CylinderCalculator, ValidIntersection)
// {
//     auto surface = create_cylinder_surface();
//     auto aperture = create_rectangle_aperture();
//     CylinderCalculator calc(surface, aperture);

//     // Use the array-based intersect method
//     // Ray starting outside cylinder and hitting it
//     double pos_loc[3] = {2.0, 0.0, 0.0}; // Start outside cylinder (radius=1) in x direction
//     double cos_loc[3] = {-1.0, 0.0, 0.0}; // Moving toward center in -x direction
//     double pos_xyz[3], cos_klm[3], df_xyz[3];
//     double path_length;

//     int result = calc.intersect(pos_loc, cos_loc, pos_xyz, cos_klm, df_xyz, &path_length);

//     // Should find intersection (result == 0 means success)
//     EXPECT_EQ(result, 0);
//     EXPECT_GT(path_length, 0.0);
// }

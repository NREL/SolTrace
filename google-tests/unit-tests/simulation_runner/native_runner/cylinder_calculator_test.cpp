#include <gtest/gtest.h>
#include <memory>
#include <cmath>

#include <aperture.hpp>
#include <cylinder_calculator.hpp>
#include <surface.hpp>

#include "common.hpp"

using SolTrace::NativeRunner::CylinderCalculator;

// Helper functions to create surfaces and apertures
std::shared_ptr<Cylinder> create_cylinder_surface(double radius = 1.0)
{
    auto cylinder = std::make_shared<Cylinder>(radius);
    return cylinder;
}

std::shared_ptr<Rectangle> create_rectangle_aperture(double x_length = 2.0, double y_length = 2.0)
{
    auto rect = std::make_shared<Rectangle>(x_length, y_length);
    return rect;
}

std::shared_ptr<Circle> create_circular_aperture(double diameter = 2.0)
{
    auto circ = std::make_shared<Circle>(diameter);
    return circ;
}

// Constructor validation tests
TEST(CylinderCalculator, ConstructorValidConstruction)
{
    auto surface = create_cylinder_surface();
    auto aperture = create_rectangle_aperture(); // x_length = 2.0 matches cylinder diameter
    EXPECT_NO_THROW({
        CylinderCalculator calc(surface, aperture);
    });
}

TEST(CylinderCalculator, ConstructorNullSurfaceThrows)
{
    auto aperture = create_rectangle_aperture();
    EXPECT_THROW({
        CylinderCalculator calc(nullptr, aperture);
    }, std::invalid_argument);
}

TEST(CylinderCalculator, ConstructorNullApertureThrows)
{
    auto surface = create_cylinder_surface();
    EXPECT_THROW({
        CylinderCalculator calc(surface, nullptr);
    }, std::invalid_argument);
}

TEST(CylinderCalculator, ConstructorWrongSurfaceTypeThrows)
{
    auto flat = std::make_shared<Flat>();
    auto aperture = create_rectangle_aperture();
    EXPECT_THROW({
        CylinderCalculator calc(flat, aperture);
    }, std::invalid_argument);
}

TEST(CylinderCalculator, ConstructorWrongApertureTypeThrows)
{
    auto surface = create_cylinder_surface();
    auto circular_ap = create_circular_aperture();
    EXPECT_THROW({
        CylinderCalculator calc(surface, circular_ap);
    }, std::invalid_argument);
}

TEST(CylinderCalculator, ConstructorZeroRadiusThrows)
{
    auto zero_radius_surface = create_cylinder_surface(0.0);
    auto aperture = create_rectangle_aperture();
    EXPECT_THROW({
        CylinderCalculator calc(zero_radius_surface, aperture);
    }, std::invalid_argument);
}

TEST(CylinderCalculator, ConstructorNegativeRadiusThrows)
{
    auto negative_radius_surface = create_cylinder_surface(-1.0);
    auto aperture = create_rectangle_aperture();
    EXPECT_THROW({
        CylinderCalculator calc(negative_radius_surface, aperture);
    }, std::invalid_argument);
}

TEST(CylinderCalculator, ConstructorNaNRadiusThrows)
{
    auto nan_radius_surface = create_cylinder_surface(std::nan(""));
    auto aperture = create_rectangle_aperture();
    EXPECT_THROW({
        CylinderCalculator calc(nan_radius_surface, aperture);
    }, std::invalid_argument);
}

TEST(CylinderCalculator, ConstructorInfiniteRadiusThrows)
{
    auto inf_radius_surface = create_cylinder_surface(std::numeric_limits<double>::infinity());
    auto aperture = create_rectangle_aperture();
    EXPECT_THROW({
        CylinderCalculator calc(inf_radius_surface, aperture);
    }, std::invalid_argument);
}

TEST(CylinderCalculator, ConstructorMismatchedApertureDimensionsThrows)
{
    auto surface = create_cylinder_surface(1.0); // radius = 1.0, diameter = 2.0
    // Create aperture with x_length != 2 * radius
    auto mismatched_aperture = create_rectangle_aperture(3.0, 2.0); // x_length = 3.0, but cylinder diameter = 2.0
    EXPECT_THROW({
        CylinderCalculator calc(surface, mismatched_aperture);
    }, std::invalid_argument);
}

TEST(CylinderCalculator, ConstructorZeroApertureDimensionsThrows)
{
    auto surface = create_cylinder_surface();
    auto zero_aperture = create_rectangle_aperture(2.0, 0.0); // y_length = 0
    EXPECT_THROW({
        CylinderCalculator calc(surface, zero_aperture);
    }, std::invalid_argument);
}

// Basic intersection test
TEST(CylinderCalculator, ValidIntersection)
{
    auto surface = create_cylinder_surface();
    auto aperture = create_rectangle_aperture();
    CylinderCalculator calc(surface, aperture);
    
    // Use the array-based intersect method
    // Ray starting outside cylinder and hitting it
    double pos_loc[3] = {2.0, 0.0, 0.0}; // Start outside cylinder (radius=1) in x direction
    double cos_loc[3] = {-1.0, 0.0, 0.0}; // Moving toward center in -x direction
    double pos_xyz[3], cos_klm[3], df_xyz[3];
    double path_length;
    
    int result = calc.intersect(pos_loc, cos_loc, pos_xyz, cos_klm, df_xyz, 0, &path_length);
    
    // Should find intersection (result == 0 means success)
    EXPECT_EQ(result, 0);
    EXPECT_GT(path_length, 0.0);
}

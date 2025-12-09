#include <gtest/gtest.h>

#include <calculator_factory.hpp>
#include <native_runner_types.hpp>
#include <surface.hpp>
#include <aperture.hpp>

#include <cylinder_calculator.hpp>
#include <flat_calculator.hpp>
#include <parabola_calculator.hpp>
#include <sphere_calculator.hpp>

#include "common.hpp"

using SolTrace::NativeRunner::CalculatorFactory;
using SolTrace::NativeRunner::CylinderCalculator;
using SolTrace::NativeRunner::ElementParameters;
using SolTrace::NativeRunner::FlatCalculator;
using SolTrace::NativeRunner::ParabolaCalculator;
using SolTrace::NativeRunner::SphereCalculator;

class CalculatorFactoryTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        factory = CalculatorFactory::get();
        eparams = ElementParameters();
    }

    void TearDown() override
    {
        // Note: CalculatorFactory is a singleton, so we don't delete it
    }

    CalculatorFactory *factory;
    ElementParameters eparams;
};

struct UnknownSurface : public Surface
{
    UnknownSurface() : Surface(SurfaceType::SURFACE_UNKNOWN){}
    ~UnknownSurface() {}
    virtual void write_json(nlohmann::ordered_json& json) const override {}
};

TEST_F(CalculatorFactoryTest, SingletonBehavior)
{
    // Test that multiple calls to get() return the same instance
    auto factory1 = CalculatorFactory::get();
    auto factory2 = CalculatorFactory::get();

    EXPECT_EQ(factory1, factory2);
    EXPECT_NE(factory1, nullptr);
}

TEST_F(CalculatorFactoryTest, CreateParabolaCalculator)
{
    // Create a parabola surface
    auto surf = SolTrace::Data::make_surface<SolTrace::Data::Parabola>(10.0, 10.0);
    auto ap = SolTrace::Data::make_aperture<SolTrace::Data::Circle>(5.0);

    // Create calculator
    auto calc = factory->make_calculator(ap, surf, eparams);

    // Verify the calculator was created and is of correct type
    EXPECT_NE(calc, nullptr);

    // Try to cast to specific type to verify it's correct
    auto parabola_calc = std::dynamic_pointer_cast<ParabolaCalculator>(calc);
    EXPECT_NE(parabola_calc, nullptr);
}

TEST_F(CalculatorFactoryTest, CreateFlatCalculator)
{
    // Create a flat surface
    auto surf = SolTrace::Data::make_surface<SolTrace::Data::Flat>();
    auto ap = SolTrace::Data::make_aperture<SolTrace::Data::Rectangle>(2.0, 3.0);

    // Create calculator
    auto calc = factory->make_calculator(ap, surf, eparams);

    // Verify the calculator was created and is of correct type
    EXPECT_NE(calc, nullptr);

    // Try to cast to specific type to verify it's correct
    auto flat_calc = std::dynamic_pointer_cast<FlatCalculator>(calc);
    EXPECT_NE(flat_calc, nullptr);
}

TEST_F(CalculatorFactoryTest, CreateCylinderCalculator)
{
    // Create a cylinder surface
    auto surf = SolTrace::Data::make_surface<SolTrace::Data::Cylinder>(5.0);
    // CylinderCalculator specifically requires a Rectangle aperture
    auto ap = SolTrace::Data::make_aperture<SolTrace::Data::Rectangle>(10.0, 8.0); // x_length should be 2*radius = 10.0

    // Create calculator
    auto calc = factory->make_calculator(ap, surf, eparams);

    // Verify the calculator was created and is of correct type
    EXPECT_NE(calc, nullptr);

    // Try to cast to specific type to verify it's correct
    auto cylinder_calc = std::dynamic_pointer_cast<CylinderCalculator>(calc);
    EXPECT_NE(cylinder_calc, nullptr);
}

TEST_F(CalculatorFactoryTest, CreateSphereCalculator)
{
    // Create a sphere surface
    auto surf = SolTrace::Data::make_surface<SolTrace::Data::Sphere>(7.5);
    auto ap = SolTrace::Data::make_aperture<SolTrace::Data::Hexagon>(4.0);

    // Create calculator
    auto calc = factory->make_calculator(ap, surf, eparams);

    // Verify the calculator was created and is of correct type
    EXPECT_NE(calc, nullptr);

    // Try to cast to specific type to verify it's correct
    auto sphere_calc = std::dynamic_pointer_cast<SphereCalculator>(calc);
    EXPECT_NE(sphere_calc, nullptr);
}

TEST_F(CalculatorFactoryTest, ErrorHandling_NullSurface)
{
    // Test with null surface pointer
    auto ap = SolTrace::Data::make_aperture<SolTrace::Data::Circle>(5.0);
    surface_ptr null_surf = nullptr;

    // Should throw std::invalid_argument
    EXPECT_THROW(factory->make_calculator(ap, null_surf, eparams),
                 std::invalid_argument);
}

TEST_F(CalculatorFactoryTest, ErrorHandling_NullAperture)
{
    // Test with null aperture pointer
    aperture_ptr null_ap = nullptr;
    surface_ptr surf = SolTrace::Data::make_surface<SolTrace::Data::Flat>();

    // Should throw std::invalid_argument
    EXPECT_THROW(factory->make_calculator(null_ap, surf, eparams),
                 std::invalid_argument);
}

TEST_F(CalculatorFactoryTest, ErrorHandling_UnsupportedSurfaceType)
{
    // Create a surface with an unsupported type (CONE is not handled)
    auto cone_surf = SolTrace::Data::make_surface<SolTrace::Data::Cone>(0.5); // half angle
    auto ap = SolTrace::Data::make_aperture<SolTrace::Data::Circle>(5.0);

    // Should throw std::invalid_argument for unsupported surface type
    EXPECT_THROW({ factory->make_calculator(ap, cone_surf, eparams); }, std::invalid_argument);
}

TEST_F(CalculatorFactoryTest, MultipleCalculatorCreation)
{
    // Test creating multiple calculators in sequence
    auto ap = SolTrace::Data::make_aperture<SolTrace::Data::Rectangle>(2.0, 3.0);

    // Create multiple flat calculators
    auto calc1 = factory->make_calculator(ap, SolTrace::Data::make_surface<SolTrace::Data::Flat>(), eparams);
    auto calc2 = factory->make_calculator(ap, SolTrace::Data::make_surface<SolTrace::Data::Flat>(), eparams);
    auto calc3 = factory->make_calculator(ap, SolTrace::Data::make_surface<SolTrace::Data::Flat>(), eparams);

    // All should be valid but different instances
    EXPECT_NE(calc1, nullptr);
    EXPECT_NE(calc2, nullptr);
    EXPECT_NE(calc3, nullptr);
    EXPECT_NE(calc1.get(), calc2.get());
    EXPECT_NE(calc2.get(), calc3.get());
    EXPECT_NE(calc1.get(), calc3.get());
}

TEST_F(CalculatorFactoryTest, ExceptionMessageContent)
{
    // Test that exception messages contain useful information

    // Test null surface exception message
    auto ap = SolTrace::Data::make_aperture<SolTrace::Data::Circle>(5.0);
    try
    {
        factory->make_calculator(ap, nullptr, eparams);
        FAIL() << "Expected std::invalid_argument to be thrown";
    }
    catch (const std::invalid_argument &e)
    {
        std::string msg = e.what();
        EXPECT_TRUE(msg.find("Surface pointer cannot be null") != std::string::npos);
        EXPECT_TRUE(msg.find("CalculatorFactory::make_calculator") != std::string::npos);
    }

    // Test cylinder with null aperture exception message
    auto cylinder_surf = SolTrace::Data::make_surface<SolTrace::Data::Cylinder>(5.0);
    try
    {
        factory->make_calculator(nullptr, cylinder_surf, eparams);
        FAIL() << "Expected std::invalid_argument to be thrown";
    }
    catch (const std::invalid_argument &e)
    {
        std::string msg = e.what();
        EXPECT_TRUE(msg.find("Aperture pointer cannot be null") != std::string::npos);
        EXPECT_TRUE(msg.find("CalculatorFactory::make_calculator") != std::string::npos);
    }

    // Test unsupported surface type exception message
    auto cone_surf = SolTrace::Data::make_surface<UnknownSurface>();
    try
    {
        factory->make_calculator(ap, cone_surf, eparams);
        FAIL() << "Expected std::invalid_argument to be thrown";
    }
    catch (const std::invalid_argument &e)
    {
        std::string msg = e.what();
        EXPECT_TRUE(msg.find("Unsupported surface type") != std::string::npos);
        EXPECT_TRUE(msg.find("CalculatorFactory::make_calculator") != std::string::npos);
    }
}

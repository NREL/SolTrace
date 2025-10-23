#include <gtest/gtest.h>

#include <sun.hpp>
#include <limits>
#include <cmath>

#include "common.hpp"

// Test Gaussian distribution error handling
TEST(Sun, GaussianDistributionErrors)
{
    Sun sun;
    
    // Test negative sigma
    EXPECT_THROW(
        sun.set_shape(SolTrace::Data::SunShape::GAUSSIAN, -1.0, 0.0, 0.0, {}, {}),
        std::invalid_argument
    );
    
    // Test zero sigma
    EXPECT_THROW(
        sun.set_shape(SolTrace::Data::SunShape::GAUSSIAN, 0.0, 0.0, 0.0, {}, {}),
        std::invalid_argument
    );
    
    // Test NaN sigma
    EXPECT_THROW(
        sun.set_shape(SolTrace::Data::SunShape::GAUSSIAN, std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0, {}, {}),
        std::invalid_argument
    );
    
    // Test infinite sigma
    EXPECT_THROW(
        sun.set_shape(SolTrace::Data::SunShape::GAUSSIAN, std::numeric_limits<double>::infinity(), 0.0, 0.0, {}, {}),
        std::invalid_argument
    );
}

// Test Pillbox distribution error handling
TEST(Sun, PillboxDistributionErrors)
{
    Sun sun;
    
    // Test negative half_width
    EXPECT_THROW(
        sun.set_shape(SolTrace::Data::SunShape::PILLBOX, 0.0, -1.0, 0.0, {}, {}),
        std::invalid_argument
    );
    
    // Test zero half_width
    EXPECT_THROW(
        sun.set_shape(SolTrace::Data::SunShape::PILLBOX, 0.0, 0.0, 0.0, {}, {}),
        std::invalid_argument
    );
    
    // Test NaN half_width
    EXPECT_THROW(
        sun.set_shape(SolTrace::Data::SunShape::PILLBOX, 0.0, std::numeric_limits<double>::quiet_NaN(), 0.0, {}, {}),
        std::invalid_argument
    );
    
    // Test infinite half_width
    EXPECT_THROW(
        sun.set_shape(SolTrace::Data::SunShape::PILLBOX, 0.0, std::numeric_limits<double>::infinity(), 0.0, {}, {}),
        std::invalid_argument
    );
}

// Test Buie CSR distribution error handling
TEST(Sun, BuieCsrDistributionErrors)
{
    Sun sun;

    // Test negative CSR
    EXPECT_THROW(
        sun.set_shape(SolTrace::Data::SunShape::BUIE_CSR, 0.0, 0.0, -1.0, {}, {}),
        std::invalid_argument
    );

    // Test high CSR
    EXPECT_THROW(
        sun.set_shape(SolTrace::Data::SunShape::BUIE_CSR, 0.0, 0.0, 0.81, {}, {}),
        std::invalid_argument
    );

    // Test NaN CSR
    EXPECT_THROW(
        sun.set_shape(SolTrace::Data::SunShape::BUIE_CSR, 0.0, 0.0, std::numeric_limits<double>::quiet_NaN(), {}, {}),
        std::invalid_argument
    );

    // Test infinite CSR
    EXPECT_THROW(
        sun.set_shape(SolTrace::Data::SunShape::BUIE_CSR, 0.0, 0.0, std::numeric_limits<double>::infinity(), {}, {}),
        std::invalid_argument
    );
}

// Test User-defined distribution error handling
TEST(Sun, UserDefinedDistributionErrors)
{
    Sun sun;
    
    // Test empty angle vector
    EXPECT_THROW(
        sun.set_shape(SolTrace::Data::SunShape::USER_DEFINED, 0.0, 0.0, 0.0, {}, {1.0}),
        std::invalid_argument
    );
    
    // Test empty intensity vector
    EXPECT_THROW(
        sun.set_shape(SolTrace::Data::SunShape::USER_DEFINED, 0.0, 0.0, 0.0, {1.0}, {}),
        std::invalid_argument
    );
    
    // Test both vectors empty
    EXPECT_THROW(
        sun.set_shape(SolTrace::Data::SunShape::USER_DEFINED, 0.0, 0.0, 0.0, {}, {}),
        std::invalid_argument
    );
    
    // Test mismatched vector sizes
    EXPECT_THROW(
        sun.set_shape(SolTrace::Data::SunShape::USER_DEFINED, 0.0, 0.0, 0.0, {1.0, 2.0}, {1.0}),
        std::invalid_argument
    );
    
    // Test negative angle
    EXPECT_THROW(
        sun.set_shape(SolTrace::Data::SunShape::USER_DEFINED, 0.0, 0.0, 0.0, {-1.0, 2.0}, {1.0, 2.0}),
        std::invalid_argument
    );
    
    // Test negative intensity
    EXPECT_THROW(
        sun.set_shape(SolTrace::Data::SunShape::USER_DEFINED, 0.0, 0.0, 0.0, {1.0, 2.0}, {-1.0, 2.0}),
        std::invalid_argument
    );
    
    // Test NaN in angles
    EXPECT_THROW(
        sun.set_shape(SolTrace::Data::SunShape::USER_DEFINED, 0.0, 0.0, 0.0,
                     {std::numeric_limits<double>::quiet_NaN(), 2.0}, {1.0, 2.0}),
        std::invalid_argument
    );
    
    // Test NaN in intensities
    EXPECT_THROW(
        sun.set_shape(SolTrace::Data::SunShape::USER_DEFINED, 0.0, 0.0, 0.0,
                     {1.0, 2.0}, {1.0, std::numeric_limits<double>::quiet_NaN()}),
        std::invalid_argument
    );
    
    // Test unsorted angles
    EXPECT_THROW(
        sun.set_shape(SolTrace::Data::SunShape::USER_DEFINED, 0.0, 0.0, 0.0, {2.0, 1.0}, {1.0, 2.0}),
        std::invalid_argument
    );
}

// Test valid cases
TEST(Sun, ValidDistributions)
{
    Sun sun;
    
    // Valid Gaussian distribution
    EXPECT_NO_THROW(
        sun.set_shape(SolTrace::Data::SunShape::GAUSSIAN, 1.0, 0.0, 0.0, {}, {})
    );
    EXPECT_EQ(sun.get_shape(), SolTrace::Data::SunShape::GAUSSIAN);
    
    // Valid Pillbox distribution
    EXPECT_NO_THROW(
        sun.set_shape(SolTrace::Data::SunShape::PILLBOX, 0.0, 2.0, 0.0, {}, {})
    );
    EXPECT_EQ(sun.get_shape(), SolTrace::Data::SunShape::PILLBOX);

	// Valid Limb Darkened distribution
    EXPECT_NO_THROW(
        sun.set_shape(SolTrace::Data::SunShape::LIMBDARKENED, 0.0, 0.0, 0.0, {}, {})
	);
	EXPECT_EQ(sun.get_shape(), SolTrace::Data::SunShape::LIMBDARKENED);

	// Valid Buie CSR distribution
    EXPECT_NO_THROW(
        sun.set_shape(SolTrace::Data::SunShape::BUIE_CSR, 0.0, 0.0, 0.1, {}, {})
	);
	EXPECT_EQ(sun.get_shape(), SolTrace::Data::SunShape::BUIE_CSR);
    
    // Valid User-defined distribution
    EXPECT_NO_THROW(
        sun.set_shape(SolTrace::Data::SunShape::USER_DEFINED, 0.0, 0.0, 0.0, {0.0, 1.0, 2.0}, {1.0, 0.8, 0.6})
    );
    EXPECT_EQ(sun.get_shape(), SolTrace::Data::SunShape::USER_DEFINED);
    
    // Valid User-defined distribution with equal angles (edge case)
    EXPECT_NO_THROW(
        sun.set_shape(SolTrace::Data::SunShape::USER_DEFINED, 0.0, 0.0, 0.0, {1.0, 1.0, 2.0}, {1.0, 0.8, 0.6})
    );
}

// Test unknown distribution type
TEST(Sun, UnknownDistributionType)
{
    Sun sun;
    
    // Cast an invalid enum value to test default case
    EXPECT_THROW(
        sun.set_shape(static_cast<SolTrace::Data::SunShape>(999), 1.0, 1.0, 0.0, {}, {}),
        std::invalid_argument
    );
}

// Test basic Sun functionality
TEST(Sun, BasicFunctionality)
{
    Sun sun;
    
    // Test position setting and getting
    sun.set_position(1.0, 2.0, 3.0);
    Vector3d pos = sun.get_position();
    EXPECT_DOUBLE_EQ(pos[0], 1.0);
    EXPECT_DOUBLE_EQ(pos[1], 2.0);
    EXPECT_DOUBLE_EQ(pos[2], 3.0);
    
    // Test position setting with Vector3d
    Vector3d new_pos(4.0, 5.0, 6.0);
    sun.set_position(new_pos);
    pos = sun.get_position();
    EXPECT_DOUBLE_EQ(pos[0], 4.0);
    EXPECT_DOUBLE_EQ(pos[1], 5.0);
    EXPECT_DOUBLE_EQ(pos[2], 6.0);
}

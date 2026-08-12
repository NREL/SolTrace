#include "utilities/euler_angles.h"

#include <gtest/gtest.h>

#include <glm/geometric.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/trigonometric.hpp>

#include <cmath>
#include <vector>

namespace
{

constexpr double kAngleTolerance = 1.0e-10;
constexpr double kQuatTolerance  = 1.0e-12;

bool same_rotation(glm::dquat const& a,
                   glm::dquat const& b,
                   double            tolerance = kQuatTolerance)
{ return glm::length(a - b) <= tolerance || glm::length(a + b) <= tolerance; }

void expect_vec_near(glm::dvec3 actual,
                     glm::dvec3 expected,
                     double     tolerance = kAngleTolerance)
{
    EXPECT_NEAR(actual.x, expected.x, tolerance);
    EXPECT_NEAR(actual.y, expected.y, tolerance);
    EXPECT_NEAR(actual.z, expected.z, tolerance);
}

glm::dvec3 degrees(glm::dvec3 value)
{
    return { glm::degrees(value.x),
             glm::degrees(value.y),
             glm::degrees(value.z) };
}

glm::dvec3 radians(glm::dvec3 value)
{
    return { glm::radians(value.x),
             glm::radians(value.y),
             glm::radians(value.z) };
}

} // namespace

TEST(EulerAnglesXYZ, RoundTripsAwayFromSingularity)
{
    std::vector<glm::dvec3> const cases {
        { 0.0, 0.0, 0.0 },    { 12.5, -30.0, 44.0 },    { -80.0, 35.0, 170.0 },
        { 90.0, 25.0, 40.0 }, { -135.0, -40.0, 210.0 },
    };

    for (auto const& degrees_case : cases)
    {
        SCOPED_TRACE("euler=(" + std::to_string(degrees_case.x) + ", " +
                     std::to_string(degrees_case.y) + ", " +
                     std::to_string(degrees_case.z) + ")");

        auto expected = radians(degrees_case);
        auto quat     = db::euler_xyz_to_quat(expected);
        auto actual   = db::compatible_euler_xyz_from_quat(quat, expected);

        expect_vec_near(actual, expected);
        EXPECT_TRUE(same_rotation(db::euler_xyz_to_quat(actual), quat));
    }
}

TEST(EulerAnglesXYZ, PreservesCompatibleValuesPastOneRevolution)
{
    glm::dvec3 previous = radians({ 370.0, -20.0, 725.0 });
    auto       quat     = db::euler_xyz_to_quat(previous);

    auto actual = db::compatible_euler_xyz_from_quat(quat, previous);

    expect_vec_near(actual, previous);
    EXPECT_TRUE(same_rotation(db::euler_xyz_to_quat(actual), quat));
}

TEST(EulerAnglesXYZ, KeepsXChannelStableAtNinetyDegrees)
{
    glm::dvec3 previous = radians({ 90.0, 0.0, 0.0 });
    auto       quat     = db::euler_xyz_to_quat(previous);

    auto refreshed = db::compatible_euler_xyz_from_quat(quat, previous);
    expect_vec_near(degrees(refreshed), { 90.0, 0.0, 0.0 });

    glm::dvec3 edited_y = radians({ 90.0, 22.5, 0.0 });
    quat                = db::euler_xyz_to_quat(edited_y);
    refreshed           = db::compatible_euler_xyz_from_quat(quat, previous);
    expect_vec_near(degrees(refreshed), { 90.0, 22.5, 0.0 });

    glm::dvec3 edited_z = radians({ 90.0, 22.5, -37.0 });
    quat                = db::euler_xyz_to_quat(edited_z);
    refreshed           = db::compatible_euler_xyz_from_quat(quat, edited_y);
    expect_vec_near(degrees(refreshed), { 90.0, 22.5, -37.0 });
}

TEST(EulerAnglesXYZ, ChoosesEquivalentSolutionNearestPrevious)
{
    glm::dvec3 previous  = radians({ 350.0, 10.0, -355.0 });
    glm::dvec3 canonical = radians({ -10.0, 10.0, 5.0 });
    auto       quat      = db::euler_xyz_to_quat(canonical);

    auto actual = db::compatible_euler_xyz_from_quat(quat, previous);

    expect_vec_near(degrees(actual), { 350.0, 10.0, -355.0 });
    EXPECT_TRUE(same_rotation(db::euler_xyz_to_quat(actual), quat));
}

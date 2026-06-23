#include "database/conversion.h"

#include "single_element.hpp"

#include <gtest/gtest.h>

#define GLM_ENABLE_EXPERIMENTAL 1

#include <glm/geometric.hpp>
#include <glm/gtc/constants.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/vec3.hpp>

#include <cmath>
#include <vector>

namespace {

constexpr double kTolerance      = 1.0e-12;
constexpr double kAngleTolerance = 1.0e-10;

void expect_vec_near(glm::dvec3 const& actual,
                     glm::dvec3 const& expected,
                     double            tolerance = kTolerance) {
    EXPECT_NEAR(actual.x, expected.x, tolerance);
    EXPECT_NEAR(actual.y, expected.y, tolerance);
    EXPECT_NEAR(actual.z, expected.z, tolerance);
}

bool same_rotation(glm::dquat const& a,
                   glm::dquat const& b,
                   double            tolerance = kTolerance) {
    return glm::length(a - b) <= tolerance || glm::length(a + b) <= tolerance;
}

struct Basis {
    glm::dvec3 x;
    glm::dvec3 y;
    glm::dvec3 z;
};

Basis materialize_in_simdata(glm::dvec3 const& direction, double roll) {
    SolTrace::Data::SingleElement element;
    element.set_origin(glm::dvec3(0.0));
    element.set_aim_vector(glm::normalize(direction) * 100.0);
    element.set_zrot_radians(roll);

    EXPECT_EQ(element.compute_coordinate_rotations(), 0);

    auto const materialized = element.get_local_to_reference();

    return Basis {
        .x = glm::normalize(materialized * glm::dvec3(1.0, 0.0, 0.0)),
        .y = glm::normalize(materialized * glm::dvec3(0.0, 1.0, 0.0)),
        .z = glm::normalize(materialized * glm::dvec3(0.0, 0.0, 1.0)),
    };
}

void expect_basis_near(Basis const& actual,
                       Basis const& expected,
                       double       tolerance = kAngleTolerance) {
    expect_vec_near(actual.x, expected.x, tolerance);
    expect_vec_near(actual.y, expected.y, tolerance);
    expect_vec_near(actual.z, expected.z, tolerance);
}

} // namespace

TEST(DirRollToQuat, IdentityForZeroDirection) {
    auto q = db::dir_roll_to_quat({ 0.0, 0.0, 0.0 }, 1.25);

    EXPECT_TRUE(same_rotation(q, glm::identity<glm::dquat>()));
}

TEST(DirRollToQuat, MapsLocalForwardToDirection) {
    glm::dvec3 direction = glm::normalize(glm::dvec3(1.0, 2.0, 3.0));

    auto q = db::dir_roll_to_quat(direction, 0.0);

    expect_vec_near(q * glm::dvec3(0.0, 0.0, 1.0), direction);
}

TEST(DirRollToQuat, AppliesRollAroundDirection) {
    glm::dvec3 direction = glm::dvec3(0.0, 0.0, 1.0);
    double     roll      = glm::half_pi<double>();

    auto q = db::dir_roll_to_quat(direction, roll);

    expect_vec_near(q * glm::dvec3(0.0, 0.0, 1.0), direction);
    expect_vec_near(q * glm::dvec3(0.0, 1.0, 0.0), glm::dvec3(-1.0, 0.0, 0.0));
}

TEST(DirRollToQuat, RoundTripsThroughQuatToDirRollAndSimData) {
    struct Case {
        glm::dvec3 direction;
        double     roll;
    };

    std::vector<Case> const cases {
        { glm::normalize(glm::dvec3(0.0, 0.0, 1.0)), 0.0 },
        { glm::normalize(glm::dvec3(0.0, 0.0, 1.0)), glm::half_pi<double>() },
        { glm::normalize(glm::dvec3(1.0, 2.0, 3.0)), 0.7 },
        { glm::normalize(glm::dvec3(-2.0, 5.0, 1.0)), -1.2 },
        { glm::normalize(glm::dvec3(0.3, -0.7, 0.64)), glm::pi<double>() },
    };

    for (auto const& item : cases) {
        SCOPED_TRACE("dir=(" + std::to_string(item.direction.x) + ", " +
                     std::to_string(item.direction.y) + ", " +
                     std::to_string(item.direction.z) +
                     "), roll=" + std::to_string(item.roll));

        auto basis0 = materialize_in_simdata(item.direction, item.roll);
        auto q0     = glm::quat_cast(glm::dmat3(basis0.x, basis0.y, basis0.z));

        glm::dvec3 decoded_direction;
        double     decoded_roll = 0.0;
        db::quat_to_dir_roll(q0, decoded_direction, decoded_roll);

        auto basis1 = materialize_in_simdata(decoded_direction, decoded_roll);

        expect_vec_near(decoded_direction, item.direction, kAngleTolerance);
        expect_basis_near(basis1, basis0);
    }
}

TEST(DirRollToQuat, PreservesPowerTowerSurroundElement10028Roll) {
    glm::dvec3 const expected_x_axis = glm::normalize(
        glm::dvec3(0.694036915425, 0.719938293569, -0.00127022849067));
    glm::dvec3 const expected_y_axis = glm::normalize(
        glm::dvec3(-0.525964115736, 0.50824458421, 0.681945152912));
    glm::dvec3 const expected_z_axis = glm::normalize(
        glm::dvec3(0.491604016446, -0.472627015811, 0.731402211468));

    auto q0 = glm::quat_cast(
        glm::dmat3(expected_x_axis, expected_y_axis, expected_z_axis));

    glm::dvec3 decoded_direction;
    double     decoded_roll = 0.0;
    db::quat_to_dir_roll(q0, decoded_direction, decoded_roll);

    auto basis = materialize_in_simdata(decoded_direction, decoded_roll);

    expect_vec_near(decoded_direction, expected_z_axis, kAngleTolerance);
    expect_vec_near(basis.x, expected_x_axis, kAngleTolerance);
    expect_vec_near(basis.y, expected_y_axis, kAngleTolerance);
    expect_vec_near(basis.z, expected_z_axis, kAngleTolerance);
}

TEST(DirRollToQuat, MaterializesPowerTowerSurroundElement10028InSimData) {
    glm::dvec3 const expected_x_axis = glm::normalize(
        glm::dvec3(0.694036915425, 0.719938293569, -0.00127022849067));
    glm::dvec3 const expected_y_axis = glm::normalize(
        glm::dvec3(-0.525964115736, 0.50824458421, 0.681945152912));
    glm::dvec3 const expected_z_axis = glm::normalize(
        glm::dvec3(0.491604016446, -0.472627015811, 0.731402211468));

    auto q0 = glm::quat_cast(
        glm::dmat3(expected_x_axis, expected_y_axis, expected_z_axis));

    glm::dvec3 exported_direction;
    double     exported_roll = 0.0;
    db::quat_to_dir_roll(q0, exported_direction, exported_roll);

    glm::dvec3 const origin(-859.977, 573.788, -0.520187);

    SolTrace::Data::SingleElement element;
    element.set_origin(origin);
    element.set_aim_vector(origin + exported_direction * 100.0);
    element.set_zrot_radians(exported_roll);

    ASSERT_EQ(element.compute_coordinate_rotations(), 0);

    auto const materialized = element.get_local_to_reference();

    expect_vec_near(glm::normalize(materialized * glm::dvec3(1.0, 0.0, 0.0)),
                    expected_x_axis,
                    kAngleTolerance);
    expect_vec_near(glm::normalize(materialized * glm::dvec3(0.0, 1.0, 0.0)),
                    expected_y_axis,
                    kAngleTolerance);
    expect_vec_near(glm::normalize(materialized * glm::dvec3(0.0, 0.0, 1.0)),
                    expected_z_axis,
                    kAngleTolerance);
}

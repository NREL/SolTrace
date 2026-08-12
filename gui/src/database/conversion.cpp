#include "conversion.h"
#include "vector_utility.hpp"

#include <glm/gtx/quaternion.hpp>

#include <QDebug>
#include <QtMath>

namespace db {

static constexpr double PI     = M_PI;
static constexpr double TWO_PI = 2.0 * PI;

static inline double wrap_pi(double a) {
    a = std::fmod(a + PI, TWO_PI);
    if (a < 0) a += TWO_PI;
    return a - PI;
}

/// Compute a 'safe' rotation from one unit vector to another
static inline glm::dquat align_unit_vector(glm::dvec3 const& fromUnit,
                                           glm::dvec3 const& toUnit) {
    double d = glm::dot(fromUnit, toUnit);

    if (d > 1.0 - 1e-7) { return glm::identity<glm::dquat>(); }

    if (d < -1.0 + 1e-7) {
        glm::dvec3 ortho = (std::abs(fromUnit.z) < 0.9) ? glm::dvec3(0, 0, 1)
                                                        : glm::dvec3(0, 1, 0);
        glm::dvec3 axis  = glm::normalize(glm::cross(fromUnit, ortho));
        return glm::angleAxis(glm::pi<double>(), axis);
    }

    return glm::rotation(fromUnit, toUnit);
}

glm::dquat dir_roll_to_quat(glm::dvec3 const& directionWorld,
                            double            zRollRadians) {
    glm::dvec3 dir = directionWorld;
    double     len = glm::length(dir);
    if (len < 1e-8) return glm::identity<glm::dquat>();

    dir /= len;

    const glm::dvec3 localForward(0, 0, 1);
    auto            qAlign = align_unit_vector(localForward, dir);

    auto qRoll = glm::angleAxis(zRollRadians, dir);

    return glm::normalize(qRoll * qAlign);
}

void quat_to_dir_roll(glm::dquat const& qIn,
                      glm::dvec3&       outDirectionWorld,
                      double&           outZRollRadians) {
    auto q = glm::normalize(qIn);

    const glm::dvec3 localForward(0, 0, 1);
    auto             dir = q * localForward;

    double len = glm::length(dir);
    if (len < 1e-8) {
        outDirectionWorld = glm::dvec3(0, 0, 1);
        outZRollRadians   = 0.0;
        return;
    }

    dir /= len;
    outDirectionWorld = dir;

    // SimulationData does not treat zrot as an arbitrary quaternion twist
    // around the aim vector. It stores zrot as gamma in the Spencer/Murty
    // Euler convention used by CalculateTransformMatrices().
    auto local_to_reference = glm::mat3_cast(q);
    auto reference_to_local = glm::transpose(local_to_reference);

    // CalculateTransformMatrices fills:
    // RRefToLoc[1][0] = -cos(beta) * sin(gamma)
    // RRefToLoc[1][1] =  cos(beta) * cos(gamma)
    // GLM indexes matrices as m[column][row].
    outZRollRadians =
        wrap_pi(std::atan2(-reference_to_local[1][0],
                           reference_to_local[1][1]));
}


} // namespace db

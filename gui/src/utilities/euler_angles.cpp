#include "utilities/euler_angles.h"

#include <glm/geometric.hpp>
#include <glm/gtc/constants.hpp>

#include <cmath>

namespace db {

namespace {

constexpr double euler_hypot_epsilon = 0.0000375;

glm::dquat normalized_or_identity(glm::dquat quat) {
    double length = glm::length(quat);
    if (length == 0.0 || !std::isfinite(length)) {
        return { 1.0, 0.0, 0.0, 0.0 };
    }
    return quat / length;
}

struct EulerMatrix {
    glm::dmat3 value;

    // Note that this is backwards from GLM convention, but matches the blender
    // math we are using here. WATCH OUT
    double operator()(int row, int column) const { return value[row][column]; }
};

EulerMatrix quat_to_euler_matrix(glm::dquat quat) {
    auto q = normalized_or_identity(quat);

    constexpr double sqrt2 = 1.4142135623730950488;

    double q0 = sqrt2 * q.w;
    double q1 = sqrt2 * q.x;
    double q2 = sqrt2 * q.y;
    double q3 = sqrt2 * q.z;

    double qda = q0 * q1;
    double qdb = q0 * q2;
    double qdc = q0 * q3;
    double qaa = q1 * q1;
    double qab = q1 * q2;
    double qac = q1 * q3;
    double qbb = q2 * q2;
    double qbc = q2 * q3;
    double qcc = q3 * q3;

    // do NOT change the row major ordering here
    return { {
        { 1.0 - qbb - qcc, qdc + qab, -qdb + qac },
        { -qdc + qab, 1.0 - qaa - qcc, qda + qbc },
        { qdb + qac, -qda + qbc, 1.0 - qaa - qbb },
    } };
}

void compatible_euler(glm::dvec3& euler, glm::dvec3 previous) {
    constexpr double pi     = glm::pi<double>();
    constexpr double two_pi = glm::two_pi<double>();

    glm::dvec3 delta = euler - previous;

    for (int i = 0; i < 3; ++i) {
        if (delta[i] > pi) {
            euler[i] -= std::floor((delta[i] / two_pi) + 0.5) * two_pi;
            delta[i] = euler[i] - previous[i];
        } else if (delta[i] < -pi) {
            euler[i] += std::floor((-delta[i] / two_pi) + 0.5) * two_pi;
            delta[i] = euler[i] - previous[i];
        }
    }

    int j = 1;
    int k = 2;
    for (int i = 0; i < 3; j = k, k = i++) {
        if (std::abs(delta[i]) > pi &&
            std::abs(delta[j]) < glm::half_pi<double>() &&
            std::abs(delta[k]) < glm::half_pi<double>()) {
            euler[i] += delta[i] > 0.0 ? -two_pi : two_pi;
        }
    }
}

} // namespace

glm::dquat euler_xyz_to_quat(glm::dvec3 euler) {

    glm::dvec3 ts = euler * 0.5;

    auto coss = glm::cos(ts);
    auto sins = glm::sin(ts);

    double ci = coss.x;
    double cj = coss.y;
    double ch = coss.z;
    double si = sins.x;
    double sj = sins.y;
    double sh = sins.z;

    double cc = ci * ch;
    double cs = ci * sh;
    double sc = si * ch;
    double ss = si * sh;

    return normalized_or_identity(glm::dquat(cj * cc + sj * ss,
                                             cj * sc - sj * cs,
                                             cj * ss + sj * cc,
                                             cj * cs - sj * sc));
}

glm::dvec3 compatible_euler_xyz_from_quat(glm::dquat quat,
                                          glm::dvec3 previous) {
    auto matrix = quat_to_euler_matrix(quat);

    double cy = std::hypot(matrix(0, 0), matrix(0, 1));

    glm::dvec3 euler1 { 0.0 };
    glm::dvec3 euler2 { 0.0 };

    if (cy > euler_hypot_epsilon) {
        euler1.x = std::atan2(matrix(1, 2), matrix(2, 2));
        euler1.y = std::atan2(-matrix(0, 2), cy);
        euler1.z = std::atan2(matrix(0, 1), matrix(0, 0));

        euler2.x = std::atan2(-matrix(1, 2), -matrix(2, 2));
        euler2.y = std::atan2(-matrix(0, 2), -cy);
        euler2.z = std::atan2(-matrix(0, 1), -matrix(0, 0));
    } else {
        euler1.x = std::atan2(-matrix(2, 1), matrix(1, 1));
        euler1.y = std::atan2(-matrix(0, 2), cy);
        euler1.z = 0.0;
        euler2   = euler1;
    }

    compatible_euler(euler1, previous);
    compatible_euler(euler2, previous);

    auto d1 = std::abs(euler1.x - previous.x) +
              std::abs(euler1.y - previous.y) + std::abs(euler1.z - previous.z);
    auto d2 = std::abs(euler2.x - previous.x) +
              std::abs(euler2.y - previous.y) + std::abs(euler2.z - previous.z);

    return d1 > d2 ? euler2 : euler1;
}

} // namespace db

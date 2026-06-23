#include "conversion.h"
#include "vector_utility.hpp"

#include <glm/gtx/quaternion.hpp>

#include <QDebug>
#include <QtMath>

#include <random>

namespace db {

static constexpr double PI     = M_PI;
static constexpr double TWO_PI = 2.0 * PI;

static inline double wrap_pi(double a) {
    a = std::fmod(a + PI, TWO_PI);
    if (a < 0) a += TWO_PI;
    return a - PI;
}

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


// static inline glm::dvec3 random_unit_vec(std::mt19937& rng) {
//     std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
//     glm::dvec3                            v;
//     do {
//         v = glm::dvec3(dist(rng), dist(rng), dist(rng));
//     } while (glm::length2(v) < 1e-8f);
//     return glm::normalize(v);
// }

// static inline void runRoundTripTests(std::uint32_t seed  = 12345,
//                                      int           iters = 10000) {
//     std::mt19937                          rng(seed);
//     std::uniform_real_distribution<float> rollDist(-PI, PI);

//     float maxDirErr    = 0.0f;
//     float maxRollErr   = 0.0f;
//     int   quatFailures = 0;

//     for (int i = 0; i < iters; ++i) {
//         glm::dvec3 dir0  = random_unit_vec(rng);
//         float     roll0 = rollDist(rng);

//         auto q0 = dir_roll_to_quat(dir0, roll0);

//         glm::dvec3 dir1;
//         double    roll1;
//         quat_to_dir_roll(q0, dir1, roll1);

//         auto q1 = dir_roll_to_quat(dir1, roll1);

//         float dirErr  = glm::length(dir0 - dir1);
//         float rollErr = angle_diff(roll0, roll1);

//         maxDirErr  = std::max(maxDirErr, dirErr);
//         maxRollErr = std::max(maxRollErr, rollErr);

//         if (!approx_equal_quat(q0, q1, 1e-4f)) { ++quatFailures; }
//     }

//     qDebug() << "Round-trip tests: " << iters;
//     qDebug() << "Max direction error: " << maxDirErr;
//     qDebug() << "Max roll error (rad): " << maxRollErr;
//     qDebug() << "Quaternion re-encode mismatches (tolerant): " <<
//     quatFailures;
// }

// int test() {
//     auto  dir  = glm::normalize(glm::dvec3(1, 2, 3));
//     float roll = 0.7f;

//     auto q = dir_roll_to_quat(dir, roll);

//     glm::dvec3 dir2;
//     double    roll2;
//     quat_to_dir_roll(q, dir2, roll2);

//     qDebug() << "Original dir: (" << dir.x << "," << dir.y << "," << dir.z
//              << ")";
//     qDebug() << "Decoded  dir: (" << dir2.x << "," << dir2.y << "," << dir2.z
//              << ")";
//     qDebug() << "Original roll: " << roll;
//     qDebug() << "Decoded  roll: " << roll2;

//     runRoundTripTests(12345, 20000);
//     return 0;
// }

// int test2() {

//     std::vector<glm::dvec3> positions = {
//         { 300.4760, 670.9580, 3.6927 },    { 300.0650, 669.9370, 2.6376 },
//         { 299.6530, 668.9150, 1.5826 },    { 293.9950, 671.1970, 1.5826 },
//         { -1214.6600, 163.1450, -0.5042 }, { -1082.3500, 360.9840, -3.5815 },
//         { 209.9780, -847.9350, 2.9950 },   { -1003.4800, -476.8940, 3.2507 },
//     };
//     std::vector<glm::dvec3> aims = {
//         { 38.8543, 27.7015, 723.2590 },    { 38.7178, 27.3632, 722.9140 },
//         { 38.5817, 27.0256, 722.5680 },    { 36.7186, 27.7769, 722.5680 },
//         { -591.0470, -55.8965, 749.9110 }, { -501.2880, 28.8448, 739.4200 },
//         { 14.6201, -319.2480, 829.0250 },  { -413.4470, -318.4600, 794.9350
//         },
//     };
//     std::vector<float> rots =
//     { 29.184700,  29.184700,  29.184700,  29.184700,
//                                 -75.322200, -67.126400, -23.855000, 77.872600
//                                 };

//     for (int i = 0; i < positions.size(); i++) {
//         auto p = positions.at(i);
//         auto a = aims.at(i);
//         auto r = qDegreesToRadians(rots.at(i));

//         a = glm::normalize(a - p);

//         auto q = dir_roll_to_quat(a, r);

//         glm::dvec3 dir2;
//         double    roll2;
//         quat_to_dir_roll(q, dir2, roll2);

//         qDebug() << "-----";
//         qDebug() << "Original dir: (" << a.x << "," << a.y << "," << a.z <<
//         ")"; qDebug() << "Decoded  dir: (" << dir2.x << "," << dir2.y << ","
//                  << dir2.z << ")";
//         qDebug() << "Original roll: " << r;
//         qDebug() << "Decoded  roll: " << roll2;
//         qDebug() << glm::length(a - p);
//     }
//     return 0;
// }


// static int VT = test();
// static int VT2 = test2();

} // namespace db

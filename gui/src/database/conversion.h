#pragma once

#include <glm/gtc/quaternion.hpp>
#include <glm/vec3.hpp>

namespace db {

/// From a direction + roll, compute a quaternion
glm::dquat dir_roll_to_quat(glm::dvec3 const& directionWorld,
                            double            zRollRadians);

/// From a quaternion, compute a direction and roll
void quat_to_dir_roll(glm::dquat const& qIn,
                      glm::dvec3&       outDirectionWorld,
                      double&           outZRollRadians);

} // namespace db

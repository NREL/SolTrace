#pragma once

#include <glm/gtc/quaternion.hpp>
#include <glm/vec3.hpp>

namespace db {

/// Compose XYZ Euler angles in radians using Blender's XYZ convention.
glm::dquat euler_xyz_to_quat(glm::dvec3 euler);

/// Extract XYZ Euler angles in radians, choosing the solution nearest previous.
glm::dvec3 compatible_euler_xyz_from_quat(glm::dquat quat,
                                          glm::dvec3 previous);

} // namespace db

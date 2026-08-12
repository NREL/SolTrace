#pragma once

#include "database/mesh.h"
#include "utilities/asynctask.h"
#include "utilities/grid3d.h"

#include <glm/vec3.hpp>

namespace analysis {

/// Generate an isosurface mesh from a sparse scalar volume.
///
/// This is currently a work-in-progress marching-cubes style utility.
/// TODO: Move to mc33
Result<db::Mesh, QString>
volume_to_mesh(TaskControl& output, SparseGrid3D<float> volume, float isoval);
} // namespace analysis

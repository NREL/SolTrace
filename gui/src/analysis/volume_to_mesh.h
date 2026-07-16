#pragma once

#include "database/mesh.h"
#include "utilities/asynctask.h"
#include "utilities/grid3d.h"

#include <glm/vec3.hpp>

namespace analysis {

Result<db::Mesh, QString>
volume_to_mesh(TaskControl& output, SparseGrid3D<float> volume, float isoval);
}

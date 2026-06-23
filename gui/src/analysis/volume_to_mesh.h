#pragma once

#include "utilities/grid3d.h"
#include "database/mesh.h"

#include <QPromise>
#include <glm/vec3.hpp>

namespace analysis {

void volume_to_mesh(QPromise<db::Mesh>& output,
                    SparseGrid3D<float> volume,
                    float               isoval);
}

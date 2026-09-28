#pragma once

#include "data/mesh.h"
#include "support/asynctask.h"
#include "support/grid3d.h"

#include <glm/vec3.hpp>

namespace SolTrace::GUI::Analysis {

/// Generate an isosurface mesh from a sparse scalar volume.
///
/// This is currently a work-in-progress marching-cubes style utility.
/// TODO: Move to mc33
Support::Result<SolTrace::GUI::Data::Mesh, QString>
volume_to_mesh(Support::TaskControl&        output,
               Support::SparseGrid3D<float> volume,
               float                        isoval);
} // namespace SolTrace::GUI::Analysis

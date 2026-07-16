#pragma once

#include "database/simulationresult.h"
#include "utilities/asynctask.h"
#include "utilities/grid3d.h"

#include <glm/glm.hpp>
#include <glm/vec3.hpp>

namespace analysis {

Result<analysis::SparseGrid3D<float>, QString>
compute_ray_volume_raster(TaskControl&            promise,
                          unsigned                resolution,
                          db::SimulationResultPtr results);

} // namespace analysis

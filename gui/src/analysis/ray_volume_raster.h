#pragma once

#include "database/simulationresult.h"
#include "utilities/asynctask.h"
#include "utilities/grid3d.h"

#include <glm/glm.hpp>
#include <glm/vec3.hpp>

namespace analysis {

/// Rasterize ray paths from a result set into a sparse 3D volume.
///
/// This is currently a work-in-progress analysis primitive used by the flux
/// module for volume visualization.
Result<analysis::SparseGrid3D<float>, QString>
compute_ray_volume_raster(TaskControl&            promise,
                          unsigned                resolution,
                          db::SimulationResultPtr results);

} // namespace analysis

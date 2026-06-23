#pragma once

#include "database/simulationresult.h"
#include "utilities/grid3d.h"

#include <QPromise>

#include <glm/glm.hpp>
#include <glm/vec3.hpp>

namespace analysis {

void compute_ray_volume_raster(QPromise<analysis::SparseGrid3D<float>>& promise,
                               unsigned                resolution,
                               db::SimulationResultPtr results);

} // namespace analysis

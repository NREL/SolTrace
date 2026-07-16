#include "ray_volume_raster.h"
#include "vector_utility.hpp"

#include <QtConcurrent/qtconcurrentrun.h>
#include <glm/gtx/component_wise.hpp>

#include <QDebug>

namespace analysis {


inline glm::ivec3
world_to_voxel(glm::vec3 const& p, glm::vec3 const& origin, float cellSize) {
    auto rel = (p - origin) / cellSize;
    return glm::floor(rel);
}


void raster_segment(QVector<glm::ivec3>& grid,
                    glm::vec3 const&     p0_grid,
                    glm::vec3 const&     p1_grid) {
    // Direction in grid space
    glm::vec3 d = p1_grid - p0_grid;

    // Start and end voxels (integer indices)
    auto start_voxel = glm::ivec3(glm::floor(p0_grid));
    auto end_voxel   = glm::ivec3(glm::floor(p1_grid));

    // Step direction per axis (-1, 0, or +1)
    int stepX = (d.x > 0.f) ? 1 : (d.x < 0.f ? -1 : 0);
    int stepY = (d.y > 0.f) ? 1 : (d.y < 0.f ? -1 : 0);
    int stepZ = (d.z > 0.f) ? 1 : (d.z < 0.f ? -1 : 0);

    float tMaxX, tMaxY, tMaxZ;
    float tDeltaX, tDeltaY, tDeltaZ;

    auto INF = std::numeric_limits<float>::infinity();

    // X axis
    if (stepX != 0) {
        float nextVoxelBoundary =
            (stepX > 0) ? (static_cast<float>(start_voxel.x) + 1.0f)
                        : static_cast<float>(start_voxel.x);
        tMaxX   = (nextVoxelBoundary - p0_grid.x) / d.x;
        tDeltaX = 1.0f / std::fabs(d.x);
    } else {
        tMaxX   = INF;
        tDeltaX = INF;
    }

    // Y axis
    if (stepY != 0) {
        float nextVoxelBoundary =
            (stepY > 0) ? (static_cast<float>(start_voxel.y) + 1.0f)
                        : static_cast<float>(start_voxel.y);

        tMaxY   = (nextVoxelBoundary - p0_grid.y) / d.y;
        tDeltaY = 1.0f / std::fabs(d.y);
    } else {
        tMaxY   = INF;
        tDeltaY = INF;
    }

    // Z axis
    if (stepZ != 0) {
        float nextVoxelBoundary =
            (stepZ > 0) ? (static_cast<float>(start_voxel.z) + 1.0f)
                        : static_cast<float>(start_voxel.z);
        tMaxZ   = (nextVoxelBoundary - p0_grid.z) / d.z;
        tDeltaZ = 1.0f / std::fabs(d.z);
    } else {
        tMaxZ   = INF;
        tDeltaZ = INF;
    }

    // Make sure start voxel is inside before touching the grid
    grid << start_voxel;

    // Traverse until we reach the voxel containing p1
    while (start_voxel.x != end_voxel.x || start_voxel.y != end_voxel.y ||
           start_voxel.z != end_voxel.z) {
        // Advance to next voxel boundary in the dimension with smallest tMax
        if (tMaxX < tMaxY) {
            if (tMaxX < tMaxZ) {
                start_voxel.x += stepX;
                tMaxX += tDeltaX;
            } else {
                start_voxel.z += stepZ;
                tMaxZ += tDeltaZ;
            }
        } else {
            if (tMaxY < tMaxZ) {
                start_voxel.y += stepY;
                tMaxY += tDeltaY;
            } else {
                start_voxel.z += stepZ;
                tMaxZ += tDeltaZ;
            }
        }

        grid << start_voxel;
    }
}

void compute_raster_chunk(QPromise<QVector<glm::ivec3>>& promise,
                          glm::vec3                      grid_scale,
                          db::SimulationResultPtr        results,
                          std::span<db::RayRecord>       records) {

    qDebug() << Q_FUNC_INFO << records.size();

    auto const extent = results->bounds_max - results->bounds_min;

    auto to_grid_coords = [&](glm::dvec3 const& p) {
        glm::vec3 rel(0.0f);

        if (extent.x > 0.0) {
            rel.x = static_cast<float>((p.x - results->bounds_min.x) /
                                       extent.x * grid_scale.x);
        }
        if (extent.y > 0.0) {
            rel.y = static_cast<float>((p.y - results->bounds_min.y) /
                                       extent.y * grid_scale.y);
        }
        if (extent.z > 0.0) {
            rel.z = static_cast<float>((p.z - results->bounds_min.z) /
                                       extent.z * grid_scale.z);
        }

        return rel;
    };

    QVector<glm::ivec3> grid;

    // size_t complete = 1;

    for (auto const& ray : records) {

        if (ray.events.empty()) continue;

        auto last_p = to_grid_coords(ray.events[0].location);

        for (auto i = 1; i < ray.events.size(); ++i) {
            if (promise.isCanceled()) {
                promise.emplaceResult(QVector<glm::ivec3>());
                return;
            }
            auto const& this_interaction = ray.events[i];
            auto this_interaction_p = to_grid_coords(this_interaction.location);

            if (last_p == this_interaction_p) continue;

            analysis::raster_segment(grid, last_p, this_interaction_p);

            last_p = this_interaction_p;
        }

        // if (complete % 10 == 0) {
        //     qDebug() << Q_FUNC_INFO << complete << "/" << records.size();
        // }
        // complete += 1;
    }

    promise.emplaceResult(grid);
}

Result<analysis::SparseGrid3D<float>, QString>
compute_ray_volume_raster(TaskControl&            promise,
                          unsigned                resolution,
                          db::SimulationResultPtr results) {
    auto const extent = results->bounds_max - results->bounds_min;

    if (glm::any(glm::equal(extent, glm::dvec3(0)))) {
        return analysis::SparseGrid3D<float>();
    }

    // Compute volume
    auto grid_size = ceil(extent / glm::compMax(extent) * double(resolution));

    auto const grid_scale = glm::vec3(grid_size);


    qDebug() << Q_FUNC_INFO << "Computing ray grid" << grid_size[0]
             << grid_size[1] << grid_size[2];

    // split into chunks
    size_t num_chunks = std::thread::hardware_concurrency();

    if (num_chunks < 1) num_chunks = 10;

    size_t chunk_size = results->records.size() / num_chunks;

    if (results->records.size() < num_chunks) {
        num_chunks = 1;
        chunk_size = results->records.size();
    }

    std::vector<QFuture<QVector<glm::ivec3>>> chunks;

    for (int i = 0; i < num_chunks; i++) {
        auto chunk_start = i * chunk_size;
        auto chunk_ext   = chunk_size;

        if (i == num_chunks - 1) {
            chunk_ext = (results->records.size()) - chunk_start;
        }

        auto sp = std::span(results->records).subspan(chunk_start, chunk_ext);

        chunks.push_back(
            QtConcurrent::run(compute_raster_chunk, grid_scale, results, sp));
    }

    // wait for all to be done

    analysis::SparseGrid3D<float> grid;

    glm::vec3 grid_transform_scale = grid_scale / glm::vec3(extent);

    auto grid_transform_translate =
        -glm::vec3(results->bounds_min) * grid_transform_scale;
    grid.set_transform(grid_transform_translate, grid_transform_scale);

    while (chunks.size()) {

        if (promise.cancelRequested()) {

            for (auto& chunk : chunks) {
                chunk.cancel();
            }

            return return_failure(QStringLiteral("Cancelled"));
        }

        if (!chunks.back().isFinished()) {
            std::this_thread::yield();
            continue;
        }

        qDebug() << "Chunk done," << chunks.size() << "left";

        auto other = chunks.back().result();

        // std::sort(other.begin(), other.end(), [](glm::ivec3 a, glm::ivec3 b)
        // {
        //     return glm::length2(a) < glm::length2(b);
        // });

        qDebug() << "Merging partial...";
        for (auto p : std::as_const(other)) {
            grid(p) += 1;
        }

        // grid.merge_in(other);
        qDebug() << "Merging partial...done";

        chunks.pop_back();
    }

    grid.update_active_bounds();

    // normalize
    float largest = 0.0;
    for (auto& block : grid) {
        for (auto& v : block.values) {
            largest = std::max(v, largest);
        }
    }

    if (largest != 0.0) {
        for (auto& block : grid) {
            for (auto& v : block.values) {
                v /= largest;
            }
        }
    }

    qDebug() << Q_FUNC_INFO << "Grid largest value" << largest;

    return grid;
}


} // namespace analysis

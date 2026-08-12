#include "flux_map.h"

#include "utilities/grid2d.h"
#include "utilities/asynctask.h"
#include "vector_utility.hpp"

#include <QFile>
#include <QImage>
#include <QPainter>
#include <QTextStream>

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>

#define GLM_ENABLE_EXPERIMENTAL 1
#include <glm/gtx/intersect.hpp>

namespace analysis {

/// A point projected to a specific triangle of a mesh.
struct TriangleProjection {
    size_t    triangle_index = 0;
    glm::vec3 barycentric    = glm::vec3(0.0);
    glm::vec2 uv             = glm::vec2(0.0);
};

/// A per-triangle bin of ray data.
struct TriangleFluxBin {
    /// The area of the triangle in world space.
    float world_area = 0.0;

    /// Accumulated energy for this bin.
    float accumulated_energy = 0.0;

    /// How many rays intersected this triangle.
    size_t hit_count = 0;

    /// The final flux value
    float flux = 0.0f; // W/m^2 or normalized equivalent
};


/// Barycentric interpolation.
template <class T>
static T interpolate(T const& a, T const& b, T const& c, glm::vec3 bary) {
    return a * bary.x + b * bary.y + c * bary.z;
}

/// For a UV, snap to a pixel
static QPointF uv_to_pixel(glm::vec2 uv, QSize const& size) {
    auto width  = std::max(0, size.width() - 1);
    auto height = std::max(0, size.height() - 1);

    return QPointF(uv.x * width, uv.y * height);
}

/// Get the area of a triangle in world space.
static float
triangle_area(glm::vec3 const& a, glm::vec3 const& b, glm::vec3 const& c) {
    return 0.5f * glm::length(glm::cross(b - a, c - a));
}


/// Obtain barycentric coords for a point on a triangle.
static std::optional<glm::vec3> barycentric_for_point(QPointF const& p,
                                                      QPointF const& a,
                                                      QPointF const& b,
                                                      QPointF const& c) {
    auto v0 = b - a;
    auto v1 = c - a;
    auto v2 = p - a;

    auto denom = v0.x() * v1.y() - v1.x() * v0.y();

    if (std::abs(denom) < 1e-8f) { return {}; }

    float u = (v2.x() * v1.y() - v1.x() * v2.y()) / denom;
    float v = (v0.x() * v2.y() - v2.x() * v0.y()) / denom;
    float w = 1.0f - u - v;

    return glm::vec3(w, u, v);
}

/// Results for the point-to-closest triangle algorithm
struct ClosestPointResult {
    glm::vec3 point;
    glm::vec3 bary;
    float     sq_dist;
};

/// Classic find-the-closest-point-on-triangle algorithm
ClosestPointResult
closest_point_on_triangle(glm::vec3 p, glm::vec3 a, glm::vec3 b, glm::vec3 c) {
    const glm::vec3 ab = b - a;
    const glm::vec3 ac = c - a;
    const glm::vec3 ap = p - a;

    const float d1 = glm::dot(ab, ap);
    const float d2 = glm::dot(ac, ap);
    if (d1 <= 0.f && d2 <= 0.f) {
        return {
            .point   = a,
            .bary    = { 1, 0, 0 },
            .sq_dist = glm::distance2(p, a),
        };
    }

    const glm::vec3 bp = p - b;
    const float     d3 = glm::dot(ab, bp);
    const float     d4 = glm::dot(ac, bp);
    if (d3 >= 0.f && d4 <= d3) {
        return {
            .point   = b,
            .bary    = { 0, 1, 0 },
            .sq_dist = glm::distance2(p, b),
        };
    }

    const glm::vec3 cp = p - c;
    const float     d5 = glm::dot(ab, cp);
    const float     d6 = glm::dot(ac, cp);
    if (d6 >= 0.f && d5 <= d6) {
        return {
            .point   = c,
            .bary    = { 0, 0, 1 },
            .sq_dist = glm::distance2(p, c),
        };
    }

    const float vc = d1 * d4 - d3 * d2;
    if (vc <= 0.f && d1 >= 0.f && d3 <= 0.f) {
        const float v = d1 / (d1 - d3);

        auto closest = a + v * ab;
        return {
            .point   = closest,
            .bary    = { 1.0f - v, v, 0.0f },
            .sq_dist = glm::distance2(closest, p),
        };
    }

    const float vb = d5 * d2 - d1 * d6;
    if (vb <= 0.f && d2 >= 0.f && d6 <= 0.f) {
        const float w = d2 / (d2 - d6);

        auto closest = a + w * ac;
        return {
            .point   = closest,
            .bary    = { 1.0f - w, 0.0f, w },
            .sq_dist = glm::distance2(closest, p),
        };
    }

    const float va = d3 * d6 - d5 * d4;
    if (va <= 0.f && (d4 - d3) >= 0.f && (d5 - d6) >= 0.f) {
        const float v = (d4 - d3) / ((d4 - d3) + (d5 - d6));

        auto closest = b + v * (c - b);
        return {
            .point   = closest,
            .bary    = { 0.0f, 1.0f - v, v },
            .sq_dist = glm::distance2(closest, p),
        };
    }

    const float denom = 1.f / (va + vb + vc);
    const float v     = vb * denom;
    const float w     = vc * denom;

    auto closest = a + v * ab + w * ac;
    return {
        .point   = closest,
        .bary    = { 1 - v - w, v, w },
        .sq_dist = glm::distance2(closest, p),
    };
}


/// For a given point, find the closest point on a triangle mesh.
/// TODO: Add accelleration structure for mesh.
static std::optional<TriangleProjection>
project_point_to_triangle(db::Mesh const& mesh, glm::vec3 p) {
    constexpr float INIT = 100000000.0f;

    // Scan all triangles, and find the closest point on each. Best wins.
    // TODO: add cutoff so really far away points don't match.

    ptrdiff_t closest_triangle = -1;

    ClosestPointResult best;

    best.sq_dist = INIT;

    for (size_t i = 0; i < mesh.triangles.size(); i++) {
        auto const& v1 = mesh.vertex[mesh.triangles[i].x];
        auto const& v2 = mesh.vertex[mesh.triangles[i].y];
        auto const& v3 = mesh.vertex[mesh.triangles[i].z];

        auto next =
            closest_point_on_triangle(p, v1.position, v2.position, v3.position);

        if (next.sq_dist < best.sq_dist) {
            best             = next;
            closest_triangle = i;
        }
    }

    if (closest_triangle < 0) { return {}; }

    auto const& v1 = mesh.vertex[mesh.triangles[closest_triangle].x];
    auto const& v2 = mesh.vertex[mesh.triangles[closest_triangle].y];
    auto const& v3 = mesh.vertex[mesh.triangles[closest_triangle].z];

    return TriangleProjection {
        .triangle_index = size_t(closest_triangle),
        .barycentric    = best.bary,
        .uv             = interpolate(v1.uv, v2.uv, v3.uv, best.bary),
    };
}

/// Make triangle bins; the per-triangle stat bin
static std::vector<TriangleFluxBin> make_triangle_bins(db::Mesh const& mesh) {
    std::vector<TriangleFluxBin> triangles;
    triangles.reserve(mesh.triangles.size());

    for (size_t i = 0; i < mesh.triangles.size(); i++) {
        auto const& v1 = mesh.vertex[mesh.triangles[i].x];
        auto const& v2 = mesh.vertex[mesh.triangles[i].y];
        auto const& v3 = mesh.vertex[mesh.triangles[i].z];

        auto world = triangle_area(v1.position, v2.position, v3.position);

        triangles.emplace_back(TriangleFluxBin {
            .world_area = world,
        });
    }

    return triangles;
}

/// For all bins, compute the flux after all energy has been accumulated
static float
compute_triangle_flux_and_max(std::vector<TriangleFluxBin>& triangles,
                              float                         total_energy) {
    if (total_energy <= 0.0f) { return 0.0f; }

    float max_density = 0.0f;

    for (auto& triangle : triangles) {
        if (triangle.world_area <= 0.0f) { continue; }

        triangle.flux = triangle.accumulated_energy / triangle.world_area;

        max_density = std::max(max_density, triangle.flux);
    }

    return max_density;
}

/// Using triangle bins, burn stats to a 2D grid
static Grid2D<float>
raster_triangle_flux(std::vector<TriangleFluxBin> const& triangles,
                     db::Mesh const&                     mesh,
                     QSize const&                        image_size,
                     TaskControl&                        control,
                     int                                 progress_low,
                     int                                 progress_high) {
    Grid2D<float> raster(image_size.width(), image_size.height());
    raster.fill(0.0f);

    auto report_progress = [&](int item, int max_item) {
        auto a = max_item > 0 ? float(item) / float(max_item) : 1.0f;
        int  p = glm::mix(float(progress_low), float(progress_high), a);
        control.setProgressValue(p);
    };

    size_t tri_index = 0;
    for (size_t tri_index = 0; tri_index < mesh.triangles.size(); ++tri_index) {
        auto const& tri = triangles[tri_index];

        if (tri.flux <= 0.0f) {
            report_progress(int(tri_index + 1), int(triangles.size()));
            continue;
        }

        // Fill covered pixels directly before nearby-sample interpolation.

        auto const& v1 = mesh.vertex[mesh.triangles[tri_index].x];
        auto const& v2 = mesh.vertex[mesh.triangles[tri_index].y];
        auto const& v3 = mesh.vertex[mesh.triangles[tri_index].z];

        QPointF a = uv_to_pixel(v1.uv, image_size);
        QPointF b = uv_to_pixel(v2.uv, image_size);
        QPointF c = uv_to_pixel(v3.uv, image_size);

        int min_x =
            std::max(0, int(std::floor(std::min({ a.x(), b.x(), c.x() }))));
        int max_x = std::min(image_size.width() - 1,
                             int(std::ceil(std::max({ a.x(), b.x(), c.x() }))));
        int min_y =
            std::max(0, int(std::floor(std::min({ a.y(), b.y(), c.y() }))));
        int max_y = std::min(image_size.height() - 1,
                             int(std::ceil(std::max({ a.y(), b.y(), c.y() }))));

        for (int y = min_y; y <= max_y; ++y) {
            for (int x = min_x; x <= max_x; ++x) {
                QPointF p(x + 0.5, y + 0.5);

                auto bary = barycentric_for_point(p, a, b, c);
                if (!bary.has_value()) continue;

                constexpr float eps = -1e-5f;
                if (bary->x < eps || bary->y < eps || bary->z < eps) continue;

                // Assuming non-overlapping UVs, but still...
                raster(x, y) += tri.flux;
            }
        }

        if (control.cancelRequested()) return raster;

        report_progress(int(tri_index + 1), int(triangles.size()));
    }

    return raster;
}

/// Take a raster bin and color map, and burn that to a QImage
static void colorize_raster(QImage&              image,
                            Grid2D<float> const& raster,
                            QImage const&        color_map,
                            float                max_density) {
    for (int x = 0; x < image.width(); ++x) {
        for (int y = 0; y < image.height(); ++y) {
            float normalized = 0.0f;

            if (max_density > 0.0f) {
                normalized = std::clamp(raster(x, y) / max_density, 0.0f, 1.0f);
            }

            auto sample = QPoint(normalized * (color_map.width() - 1),
                                 color_map.height() / 2);

            image.setPixelColor(x, y, color_map.pixelColor(sample));
        }
    }
}

/// Draw the UV structure to an Image
static void raster_mesh_overlay(QPainter&       painter,
                                db::Mesh const& mesh,
                                QSize const&    image_size,
                                QColor const&   line_color) {
    if (!line_color.isValid()) { return; }

    painter.save();
    painter.setPen(QPen(line_color, 1.0));

    for (size_t i = 0; i < mesh.triangles.size(); i++) {
        auto const& v1 = mesh.vertex[mesh.triangles[i].x];
        auto const& v2 = mesh.vertex[mesh.triangles[i].y];
        auto const& v3 = mesh.vertex[mesh.triangles[i].z];

        QPolygonF triangle;
        triangle << uv_to_pixel(v1.uv, image_size)
                 << uv_to_pixel(v2.uv, image_size)
                 << uv_to_pixel(v3.uv, image_size);

        painter.drawPolygon(triangle);
    }

    painter.restore();
}

static QImage make_points_debug_image(std::vector<glm::vec2> const& uvs,
                                      QSize const&                  image_size,
                                      db::Mesh const&               mesh) {
    QImage image(image_size, QImage::Format_RGB32);
    image.fill(Qt::white);

    auto painter = QPainter(&image);
    raster_mesh_overlay(painter, mesh, image_size, QColor(0, 0, 0, 128));

    painter.save();
    painter.setPen(QPen(Qt::red, 2.0));

    for (auto const& uv : uvs) {
        QPointF pixel = uv_to_pixel(uv, image_size);
        painter.drawPoint(pixel);
    }

    painter.restore();

    return image;
}

static bool dump_interaction_points_csv(std::vector<glm::vec3> const& points,
                                        QString const&                path) {
    QFile file(path);
    if (!file.open(QIODevice::WriteOnly | QIODevice::Truncate |
                   QIODevice::Text)) {
        return false;
    }

    QTextStream stream(&file);
    stream << "x,y,z\n";

    for (auto const& point : points) {
        stream << point.x << ',' << point.y << ',' << point.z << '\n';
    }

    return true;
}

static BakedFluxMapStats
compute_flux_map_stats(Grid2D<float> const&                    raster,
                       std::vector<TriangleFluxBin> const&     triangles,
                       std::vector<glm::vec3> const&           interaction_points,
                       std::size_t                             source_ray_count,
                       double                                  power_per_ray) {
    BakedFluxMapStats stats;
    stats.source_ray_count  = static_cast<quint64>(source_ray_count);
    stats.plotted_ray_count = static_cast<quint64>(interaction_points.size());
    stats.power_per_ray     = power_per_ray;
    stats.plotted_power     = stats.plotted_ray_count * stats.power_per_ray;

    if (!interaction_points.empty()) {
        glm::dvec3 centroid(0.0);
        for (auto const& point : interaction_points) {
            centroid += glm::dvec3(point);
        }
        centroid /= static_cast<double>(interaction_points.size());
        stats.centroid = QVector3D(centroid.x, centroid.y, centroid.z);
    }

    if (raster.size() == 0) return stats;

    double sum       = 0.0;
    double sum_sq    = 0.0;
    double min_flux  = std::numeric_limits<double>::max();
    double peak_flux = 0.0;

    for (unsigned i = 0; i < raster.size(); ++i) {
        const double value = raster[i];
        sum += value;
        sum_sq += value * value;
        min_flux  = std::min(min_flux, value);
        peak_flux = std::max(peak_flux, value);
    }

    const double bin_count = raster.size();
    stats.peak_flux        = peak_flux;
    stats.min_flux         = min_flux == std::numeric_limits<double>::max()
                                 ? 0.0
                                 : min_flux;
    stats.average_flux     = sum / bin_count;

    const double variance =
        (bin_count * sum_sq - sum * sum) / (bin_count * bin_count);
    stats.sigma_flux = std::sqrt(std::max(0.0, variance));
    stats.uniformity = stats.average_flux > 0.0
                           ? stats.sigma_flux / stats.average_flux
                           : 0.0;

    std::size_t peak_hit_count = 0;
    double      peak_triangle_flux = 0.0;
    for (auto const& triangle : triangles) {
        if (triangle.flux > peak_triangle_flux) {
            peak_triangle_flux = triangle.flux;
            peak_hit_count     = triangle.hit_count;
        }
    }

    if (peak_hit_count > 0) {
        stats.peak_flux_uncertainty =
            100.0 / std::sqrt(static_cast<double>(peak_hit_count));
    }
    if (source_ray_count > 0) {
        stats.average_flux_uncertainty =
            100.0 / std::sqrt(static_cast<double>(source_ray_count));
    }

    return stats;
}

/// Main fluxmap compute function
Result<BakedFluxMapPtr, QString>
execute_map_generation_for(TaskControl&            control,
                           FluxMapBakeOptions      opts,
                           entt::entity            entity,
                           db::SimulationResultPtr results,
                           db::Mesh                mesh) {

    qDebug() << Q_FUNC_INFO << "starting map generation";

    // Image size makes no sense, bail
    if (!glm::all(glm::lessThan(glm::uvec2(1), opts.image_resolution))) {
        qDebug() << Q_FUNC_INFO << "image resolution is too small";
        return QStringLiteral("Image resolution is not sufficient");
    }

    // We need to know what rays have hit this entity
    auto iter = results->entity_to_ray_ids.find(entity);

    if (iter == results->entity_to_ray_ids.end()) {
        qDebug() << Q_FUNC_INFO << "no rays for this element";
        return QStringLiteral("No rays have interacted with selected element");
    }

    // TODO: rescope. might also be good to have a little util for all this
    constexpr int PROGRESS_SETUP      = 10;
    constexpr int PROGRESS_ACCUMULATE = 50;
    constexpr int PROGRESS_RASTER     = 90;
    constexpr int PROGRESS_COMPLETE   = 100;

    // Starting setup
    control.setProgressValue(PROGRESS_SETUP);

    qDebug() << Q_FUNC_INFO << "setup complete";

    // Creating image
    auto img = QImage(
        opts.image_resolution.x, opts.image_resolution.y, QImage::Format_RGB32);

    // Fill triangle bins
    auto triangles = make_triangle_bins(mesh);

    qDebug() << Q_FUNC_INFO << "ready triangle bins";

    constexpr double       power_per_ray    = 1.0;
    float                  total_ray_impact = 0.0f;
    std::vector<glm::vec2> interaction_uvs;
    std::vector<glm::vec3> interaction_points;
    interaction_uvs.reserve(iter->second.size());
    interaction_points.reserve(iter->second.size());

    auto report_progress =
        [&](int item, int max_item, int prog_low, int prog_high) {
            auto a = (float)item / (float)max_item;
            int  p = glm::mix((float)prog_low, (float)prog_high, a);
            control.setProgressValue(p);
        };

    size_t ray_count = 0;
    size_t ray_count_chunk =
        std::clamp<size_t>(iter->second.size() / 10, 1, 500);

    // Burn rays to triangle bins
    for (auto ray_index : iter->second) {
        ray_count++;

        if (ray_count % ray_count_chunk == 0) {
            report_progress(ray_count,
                            iter->second.size(),
                            PROGRESS_SETUP,
                            PROGRESS_ACCUMULATE);
        }

        ASYNC_TASK_SYNC_POINT(control);

        for (auto const& interaction : results->records.at(ray_index).events) {
            // TODO: we need to double check that this is ok
            if (interaction.entity != entity) { continue; }
            switch (interaction.event) {
            case db::RayEventType::CREATE:
            case db::RayEventType::VIRTUAL:
            case db::RayEventType::EXIT:
            case db::RayEventType::UNKNOWN: continue;
            case db::RayEventType::ABSORB:
            case db::RayEventType::REFLECT:
            case db::RayEventType::TRANSMIT: break;
            }

            auto projection =
                project_point_to_triangle(mesh, interaction.location);

            if (!projection.has_value()) { continue; }

            auto& triangle = triangles[projection->triangle_index];

            float hit_energy = static_cast<float>(power_per_ray);
            triangle.hit_count += 1;
            triangle.accumulated_energy += hit_energy;
            total_ray_impact += hit_energy;
            interaction_points.push_back(interaction.location);
            interaction_uvs.push_back(projection->uv);
        }
    }

    qDebug() << Q_FUNC_INFO << "triangle bins filled";

    ASYNC_TASK_SYNC_POINT(control);

    control.setProgressValue(PROGRESS_ACCUMULATE);


    // Burn triangle bins to raster

    auto max_density =
        compute_triangle_flux_and_max(triangles, total_ray_impact);

    auto raster = raster_triangle_flux(triangles,
                                       mesh,
                                       img.size(),
                                       control,
                                       PROGRESS_ACCUMULATE,
                                       PROGRESS_RASTER);

    qDebug() << Q_FUNC_INFO << "rastered bins";

    ASYNC_TASK_SYNC_POINT(control);

    colorize_raster(img, raster, opts.color_map, max_density);

    control.setProgressValue(PROGRESS_RASTER);

    {
        auto painter = QPainter(&img);
        raster_mesh_overlay(painter, mesh, img.size(), opts.grid_line_color);
    }

    auto points_img =
        make_points_debug_image(interaction_uvs, img.size(), mesh);
    auto stats = compute_flux_map_stats(raster,
                                        triangles,
                                        interaction_points,
                                        results->records.size(),
                                        power_per_ray);

    qDebug() << Q_FUNC_INFO << "complete";

    img = img.convertToFormat(QImage::Format_RGBA8888);

    control.setProgressValue(PROGRESS_COMPLETE);

    return std::make_shared<BakedFluxMap>(BakedFluxMap {
        .counts    = std::move(raster),
        .bin_map   = img,
        .point_map = points_img,
        .stats     = stats,
    });
}

FluxMapComputer::FluxMapComputer(QObject* parent) : QObject(parent) { }

FluxMapComputer::~FluxMapComputer() = default;

void FluxMapComputer::set_results(db::SimulationResultPtr p) {
    m_database = p;

    cancel_all();
}


/// Precondition:
/// Mesh must NOT have overlapping UVs
bool FluxMapComputer::start_generate_for(db::Entity         e,
                                         db::Mesh           mesh,
                                         FluxMapBakeOptions options) {

    if (!m_database) {
        qDebug() << Q_FUNC_INFO << "No database, bailing";
        return false;
    }

    if (options.color_map.isNull()) {
        options.color_map = QImage(":/assets/images/b_to_r_wide.png");

        if (options.color_map.isNull()) {
            qCritical() << "Missing colormaps!";
            return false;
        }
    }

    qDebug() << Q_FUNC_INFO << "Loaded colormap" << options.color_map.size()
             << options.color_map.sizeInBytes();

    auto task = launch_async_task<BakedFluxMapPtr, QString>(
        e,
        this,
        &FluxMapComputer::image_ready,
        &FluxMapComputer::image_failed,
        execute_map_generation_for,
        options,
        e,
        m_database,
        mesh);


    // Set up cancelling
    connect(this, &FluxMapComputer::cancel_all, task, &AsyncTaskBase::cancel);

    // Set up targeted cancelling
    connect(this,
            &FluxMapComputer::cancel_specific,
            task,
            [e, task](db::Entity item) {
                if (item == e) { task->cancel(); }
            });

    qDebug() << Q_FUNC_INFO << "Job kicked off to thread.";

    return true;
}

} // namespace analysis

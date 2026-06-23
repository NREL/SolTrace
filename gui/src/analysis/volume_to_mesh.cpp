#include "volume_to_mesh.h"

#include <QDebug>

#include <algorithm>
#include <array>
#include <cmath>

#define GLM_ENABLE_EXPERIMENTAL 1
#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>

namespace analysis {

namespace {

struct InterpolatedVertex {
    glm::vec3 position;
    glm::vec3 normal;
};

glm::vec3 gradient_at(SparseGrid3DLookupCache<float> const& volume,
                      glm::ivec3                            p) {
    auto left  = p - 1;
    auto right = p + 1;

    constexpr auto delta_denom = 2.0;

    glm::vec3 gradient(0.0f);

    gradient.x =
        (volume(right.x, p.y, p.z) - volume(left.x, p.y, p.z)) / delta_denom;
    gradient.y =
        (volume(p.x, right.y, p.z) - volume(p.x, left.y, p.z)) / delta_denom;
    gradient.z =
        (volume(p.x, p.y, right.z) - volume(p.x, p.y, left.z)) / delta_denom;

    return gradient;
}

InterpolatedVertex interpolate_vertex(glm::vec3 const& p0,
                                      glm::vec3 const& p1,
                                      glm::vec3 const& n0,
                                      glm::vec3 const& n1,
                                      float            v0,
                                      float            v1,
                                      float            isoval) {
    float t  = 0.5f;
    auto  dv = v1 - v0;
    if (std::abs(dv) > 1e-6f) { t = (isoval - v0) / dv; }
    t = std::clamp(t, 0.0f, 1.0f);

    auto position = glm::mix(p0, p1, t);
    auto normal   = glm::mix(n0, n1, t);
    if (glm::length2(normal) > 1e-12f) {
        normal = glm::normalize(normal);
    } else {
        auto edge = p1 - p0;
        normal    = glm::length2(edge) > 1e-12f ? glm::normalize(edge)
                                                : glm::vec3(0.0f, 0.0f, 1.0f);
    }

    return { position, normal };
}

void append_polygon(db::Mesh&                              mesh,
                    std::vector<InterpolatedVertex> const& polygon,
                    glm::vec3 const&                       preferred_normal) {
    if (polygon.size() < 3) { return; }

    glm::vec3 center(0.0f);
    for (auto const& vertex : polygon) {
        center += vertex.position;
    }
    center /= static_cast<float>(polygon.size());

    auto plane_normal = preferred_normal;
    if (glm::length2(plane_normal) <= 1e-12f) {
        for (size_t i = 0; i < polygon.size(); ++i) {
            auto const& a = polygon[i].position;
            auto const& b = polygon[(i + 1) % polygon.size()].position;
            plane_normal += glm::cross(a - center, b - center);
        }
    }
    if (glm::length2(plane_normal) <= 1e-12f) {
        plane_normal = glm::vec3(0.0f, 0.0f, 1.0f);
    } else {
        plane_normal = glm::normalize(plane_normal);
    }

    auto tangent = glm::cross(plane_normal, glm::vec3(1.0f, 0.0f, 0.0f));

    if (glm::length2(tangent) <= 1e-12f) {
        tangent = glm::cross(plane_normal, glm::vec3(0.0f, 1.0f, 0.0f));
    }

    tangent        = glm::normalize(tangent);
    auto bitangent = glm::normalize(glm::cross(plane_normal, tangent));

    std::vector<size_t> order(polygon.size());
    for (size_t i = 0; i < polygon.size(); ++i) {
        order[i] = i;
    }

    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) {
        auto ra = polygon[a].position - center;
        auto rb = polygon[b].position - center;
        auto angleA =
            std::atan2(glm::dot(ra, bitangent), glm::dot(ra, tangent));
        auto angleB =
            std::atan2(glm::dot(rb, bitangent), glm::dot(rb, tangent));
        return angleA < angleB;
    });

    glm::vec3 ordered_normal(0.0f);
    for (size_t i = 0; i < order.size(); ++i) {
        auto const& a = polygon[order[i]].position;
        auto const& b = polygon[order[(i + 1) % order.size()]].position;
        ordered_normal += glm::cross(a - center, b - center);
    }
    if (glm::dot(ordered_normal, plane_normal) < 0.0f) {
        std::reverse(order.begin(), order.end());
    }

    auto emit_vertex = [&](InterpolatedVertex const& vertex) -> uint32_t {
        auto normal = vertex.normal;
        if (glm::dot(normal, plane_normal) < 0.0f) {
            normal = -normal;
        }

        mesh.vertex.push_back(db::Vertex {
            .position = vertex.position,
            .normal   = normal,
            .uv       = glm::vec2(0.0f),
        });
        return static_cast<uint32_t>(mesh.vertex.size() - 1);
    };

    auto const base = emit_vertex(polygon[order[0]]);
    for (size_t i = 1; i + 1 < order.size(); ++i) {
        auto const i1 = emit_vertex(polygon[order[i]]);
        auto const i2 = emit_vertex(polygon[order[i + 1]]);
        mesh.triangles.push_back({ base, i1, i2 });
    }
}


/// Given a tetrahedra, extract faces for a given iso val
void polygonize_tetrahedron(db::Mesh&                        mesh,
                            std::vector<InterpolatedVertex>& cached_polygon,
                            std::array<glm::vec3, 4> const&  positions,
                            std::array<glm::vec3, 4> const&  gradients,
                            std::array<float, 4> const&      values,
                            float                            isoval) {
    cached_polygon.clear();

    static constexpr std::array<std::array<int, 2>, 6> tetra_edges = { {
        { 0, 1 },
        { 0, 2 },
        { 0, 3 },
        { 1, 2 },
        { 1, 3 },
        { 2, 3 },
    } };

    for (auto const& edge : tetra_edges) {
        auto const a  = edge[0];
        auto const b  = edge[1];
        auto const va = values[a];
        auto const vb = values[b];
        auto const crosses =
            (va < isoval && vb >= isoval) || (vb < isoval && va >= isoval);
        if (!crosses) { continue; }

        auto vertex = interpolate_vertex(positions[a],
                                         positions[b],
                                         gradients[a],
                                         gradients[b],
                                         va,
                                         vb,
                                         isoval);

        bool duplicate = false;
        for (auto const& existing : cached_polygon) {
            if (glm::length2(existing.position - vertex.position) < 1e-10f) {
                duplicate = true;
                break;
            }
        }
        if (!duplicate) { cached_polygon.push_back(vertex); }
    }

    if (cached_polygon.size() < 3) { return; }

    glm::vec3 preferred_normal(0.0f);
    for (auto const& vertex : cached_polygon) {
        preferred_normal += vertex.normal;
    }
    if (glm::length2(preferred_normal) > 1e-12f) {
        preferred_normal = glm::normalize(preferred_normal);
    }

    append_polygon(mesh, cached_polygon, preferred_normal);
}

} // namespace

void volume_to_mesh(QPromise<db::Mesh>& output,
                    SparseGrid3D<float> input_volume,
                    float               isoval) {


    qDebug() << Q_FUNC_INFO << "generating isosurf @" << isoval;

    auto active_span = input_volume.active_span();
    auto grid_scale  = input_volume.scale();

    auto reader = SparseGrid3DLookupCache(input_volume);

    if (active_span.x < 2 || active_span.y < 2 || active_span.z < 2 ||
        grid_scale.x == 0.0f || grid_scale.y == 0.0f ||
        grid_scale.z == 0.0f) {
        qDebug() << Q_FUNC_INFO << "invalid volume dimensions";
        output.emplaceResult(db::Mesh {});
        return;
    }

    auto to_world_position = [&](glm::vec3 const& grid_position) {
        return input_volume.grid_to_world(grid_position + 0.5f);
    };

    // Cube index decomposition to tetrahedra
    static constexpr std::array<std::array<int, 4>, 6> tetrahedra = {
        std::array<int, 4> { 0, 5, 1, 6 }, std::array<int, 4> { 0, 1, 2, 6 },
        std::array<int, 4> { 0, 2, 3, 6 }, std::array<int, 4> { 0, 3, 7, 6 },
        std::array<int, 4> { 0, 7, 4, 6 }, std::array<int, 4> { 0, 4, 5, 6 },
    };

    db::Mesh mesh;

    std::vector<InterpolatedVertex> cached_polygon;
    cached_polygon.reserve(10);


    // TODO: multithread

    auto active_min = input_volume.active_min();
    auto active_max = input_volume.active_max();

    for (int z = active_min.z; z + 1 < active_max.z; ++z) {
        for (int y = active_min.y; y + 1 < active_max.y; ++y) {
            // dont want to check too often
            if (output.isCanceled()) { return; }

            for (int x = active_min.x; x + 1 < active_max.x; ++x) {
                std::array<glm::ivec3, 8> const corner_coords = { {
                    { x, y, z },
                    { x + 1, y, z },
                    { x + 1, y + 1, z },
                    { x, y + 1, z },
                    { x, y, z + 1 },
                    { x + 1, y, z + 1 },
                    { x + 1, y + 1, z + 1 },
                    { x, y + 1, z + 1 },
                } };

                std::array<glm::vec3, 8> positions;
                std::array<glm::vec3, 8> gradients;
                std::array<float, 8>     values;

                bool has_below = false;
                bool has_above = false;

                for (size_t i = 0; i < corner_coords.size(); ++i) {
                    auto const& c = corner_coords[i];
                    positions[i]  = to_world_position(c);
                    values[i]     = reader(c[0], c[1], c[2]);
                    gradients[i]  = gradient_at(reader, c);
                    has_below |= values[i] < isoval;
                    has_above |= values[i] >= isoval;
                }

                if (!(has_below && has_above)) { continue; }

                for (auto const& tetra : tetrahedra) {
                    std::array<glm::vec3, 4> tetra_positions;
                    std::array<glm::vec3, 4> tetra_gradients;
                    std::array<float, 4>     tetra_values;

                    for (size_t i = 0; i < tetra.size(); ++i) {
                        tetra_positions[i] = positions[tetra[i]];
                        tetra_gradients[i] = gradients[tetra[i]];
                        tetra_values[i]    = values[tetra[i]];
                    }

                    polygonize_tetrahedron(mesh,
                                           cached_polygon,
                                           tetra_positions,
                                           tetra_gradients,
                                           tetra_values,
                                           isoval);
                }
            }
        }
    }

    if (mesh.triangles.empty()) {
        qDebug() << Q_FUNC_INFO << "no triangles generated";
    }

    output.emplaceResult(std::move(mesh));
}

} // namespace analysis

#include "triangle_bvh.h"

#define GLM_ENABLE_EXPERIMENTAL 1

#include <glm/common.hpp>
#include <glm/geometric.hpp>
#include <glm/gtx/norm.hpp>

#include <algorithm>
#include <limits>

namespace analysis {

namespace {

constexpr size_t LEAF_TRIANGLE_COUNT = 12;

float component(glm::vec3 const& value, int axis) {
    return value[axis];
}

} // namespace

void TriangleBvh::Aabb::expand(glm::vec3 point) {
    if (!valid) {
        min   = point;
        max   = point;
        valid = true;
        return;
    }

    min = glm::min(min, point);
    max = glm::max(max, point);
}

void TriangleBvh::Aabb::expand(Aabb const& box) {
    if (!box.valid) return;

    expand(box.min);
    expand(box.max);
}

float TriangleBvh::Aabb::distance_squared(glm::vec3 point) const {
    glm::vec3 clamped = glm::clamp(point, min, max);
    return glm::distance2(point, clamped);
}

int TriangleBvh::Aabb::widest_axis() const {
    glm::vec3 extent = max - min;

    if (extent.y > extent.x && extent.y >= extent.z) return 1;
    if (extent.z > extent.x && extent.z > extent.y) return 2;
    return 0;
}

bool TriangleBvh::Node::leaf() const {
    return left < 0 && right < 0;
}

TriangleBvh::TriangleBvh(db::Mesh const& mesh) : m_mesh(&mesh) {
    m_indices.reserve(mesh.triangles.size());

    for (size_t i = 0; i < mesh.triangles.size(); ++i) {
        auto const& tri = mesh.triangles[i];

        if (tri.x >= mesh.vertex.size() || tri.y >= mesh.vertex.size() ||
            tri.z >= mesh.vertex.size()) {
            continue;
        }

        m_indices.push_back(i);
    }

    if (!m_indices.empty()) {
        m_nodes.reserve(m_indices.size() * 2);
        build_node(0, m_indices.size());
    }
}

bool TriangleBvh::empty() const {
    return m_nodes.empty();
}

std::optional<TriangleBvh::ClosestResult>
TriangleBvh::closest_point(glm::vec3 point) const {
    if (m_nodes.empty()) return { };

    float best_sq = std::numeric_limits<float>::max();
    return closest_in_node(0, point, best_sq);
}

TriangleBvh::Aabb TriangleBvh::triangle_bounds(size_t triangle_index) const {
    auto const& tri = m_mesh->triangles[triangle_index];

    Aabb bounds;
    bounds.expand(m_mesh->vertex[tri.x].position);
    bounds.expand(m_mesh->vertex[tri.y].position);
    bounds.expand(m_mesh->vertex[tri.z].position);
    return bounds;
}

glm::vec3 TriangleBvh::triangle_centroid(size_t triangle_index) const {
    auto const& tri = m_mesh->triangles[triangle_index];

    return (m_mesh->vertex[tri.x].position + m_mesh->vertex[tri.y].position +
            m_mesh->vertex[tri.z].position) /
           3.0f;
}

int TriangleBvh::build_node(size_t begin, size_t count) {
    Node node;
    node.begin = begin;
    node.count = count;

    Aabb centroid_bounds;
    for (size_t i = begin; i < begin + count; ++i) {
        auto triangle_index = m_indices[i];
        node.bounds.expand(triangle_bounds(triangle_index));
        centroid_bounds.expand(triangle_centroid(triangle_index));
    }

    const int node_index = static_cast<int>(m_nodes.size());
    m_nodes.push_back(node);

    if (count <= LEAF_TRIANGLE_COUNT) return node_index;

    const int   axis         = centroid_bounds.widest_axis();
    const float split_extent = component(centroid_bounds.max, axis) -
                               component(centroid_bounds.min, axis);

    if (split_extent <= 0.0f) return node_index;

    const size_t mid = begin + count / 2;
    std::nth_element(m_indices.begin() + begin,
                     m_indices.begin() + mid,
                     m_indices.begin() + begin + count,
                     [this, axis](size_t a, size_t b) {
                         return component(triangle_centroid(a), axis) <
                                component(triangle_centroid(b), axis);
                     });

    m_nodes[node_index].left  = build_node(begin, mid - begin);
    m_nodes[node_index].right = build_node(mid, begin + count - mid);

    return node_index;
}

std::optional<TriangleBvh::ClosestResult>
TriangleBvh::closest_in_node(int       node_index,
                             glm::vec3 point,
                             float&    best_sq) const {
    auto const& node = m_nodes[node_index];

    if (node.bounds.distance_squared(point) > best_sq) return { };

    if (node.leaf()) {
        std::optional<ClosestResult> best;

        for (size_t i = node.begin; i < node.begin + node.count; ++i) {
            auto next = closest_on_triangle(m_indices[i], point);

            if (next.sq_distance < best_sq) {
                best_sq = next.sq_distance;
                best    = next;
            }
        }

        return best;
    }

    auto const& left  = m_nodes[node.left];
    auto const& right = m_nodes[node.right];

    const float left_dist  = left.bounds.distance_squared(point);
    const float right_dist = right.bounds.distance_squared(point);

    std::optional<ClosestResult> best;

    auto visit = [&](int child) {
        auto next = closest_in_node(child, point, best_sq);
        if (next) best = next;
    };

    if (left_dist < right_dist) {
        visit(node.left);
        visit(node.right);
    } else {
        visit(node.right);
        visit(node.left);
    }

    return best;
}

TriangleBvh::ClosestResult
TriangleBvh::closest_on_triangle(size_t triangle_index, glm::vec3 point) const {
    auto const& tri = m_mesh->triangles[triangle_index];

    const glm::vec3 a = m_mesh->vertex[tri.x].position;
    const glm::vec3 b = m_mesh->vertex[tri.y].position;
    const glm::vec3 c = m_mesh->vertex[tri.z].position;

    const glm::vec3 ab = b - a;
    const glm::vec3 ac = c - a;
    const glm::vec3 ap = point - a;

    const float d1 = glm::dot(ab, ap);
    const float d2 = glm::dot(ac, ap);
    if (d1 <= 0.0f && d2 <= 0.0f) {
        return {
            .triangle_index = triangle_index,
            .point          = a,
            .barycentric    = { 1.0f, 0.0f, 0.0f },
            .sq_distance    = glm::distance2(point, a),
        };
    }

    const glm::vec3 bp = point - b;
    const float     d3 = glm::dot(ab, bp);
    const float     d4 = glm::dot(ac, bp);
    if (d3 >= 0.0f && d4 <= d3) {
        return {
            .triangle_index = triangle_index,
            .point          = b,
            .barycentric    = { 0.0f, 1.0f, 0.0f },
            .sq_distance    = glm::distance2(point, b),
        };
    }

    const glm::vec3 cp = point - c;
    const float     d5 = glm::dot(ab, cp);
    const float     d6 = glm::dot(ac, cp);
    if (d6 >= 0.0f && d5 <= d6) {
        return {
            .triangle_index = triangle_index,
            .point          = c,
            .barycentric    = { 0.0f, 0.0f, 1.0f },
            .sq_distance    = glm::distance2(point, c),
        };
    }

    const float vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0f && d1 >= 0.0f && d3 <= 0.0f) {
        const float v       = d1 / (d1 - d3);
        const auto  closest = a + v * ab;

        return {
            .triangle_index = triangle_index,
            .point          = closest,
            .barycentric    = { 1.0f - v, v, 0.0f },
            .sq_distance    = glm::distance2(closest, point),
        };
    }

    const float vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0f && d2 >= 0.0f && d6 <= 0.0f) {
        const float w       = d2 / (d2 - d6);
        const auto  closest = a + w * ac;

        return {
            .triangle_index = triangle_index,
            .point          = closest,
            .barycentric    = { 1.0f - w, 0.0f, w },
            .sq_distance    = glm::distance2(closest, point),
        };
    }

    const float va = d3 * d6 - d5 * d4;
    if (va <= 0.0f && (d4 - d3) >= 0.0f && (d5 - d6) >= 0.0f) {
        const float v       = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        const auto  closest = b + v * (c - b);

        return {
            .triangle_index = triangle_index,
            .point          = closest,
            .barycentric    = { 0.0f, 1.0f - v, v },
            .sq_distance    = glm::distance2(closest, point),
        };
    }

    const float denom   = 1.0f / (va + vb + vc);
    const float v       = vb * denom;
    const float w       = vc * denom;
    const auto  closest = a + v * ab + w * ac;

    return {
        .triangle_index = triangle_index,
        .point          = closest,
        .barycentric    = { 1.0f - v - w, v, w },
        .sq_distance    = glm::distance2(closest, point),
    };
}

} // namespace analysis

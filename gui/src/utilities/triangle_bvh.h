#pragma once

#include "database/mesh.h"

#include <glm/vec3.hpp>

#include <cstddef>
#include <optional>
#include <vector>

namespace analysis {

/// Exact closest-point queries against a triangle mesh, accelerated by an AABB
/// tree over triangle bounds.
///
/// Used explicitly for our flux map comp.
class TriangleBvh {
public:
    struct ClosestResult {
        size_t    triangle_index = 0;
        glm::vec3 point          = glm::vec3(0.0f);
        glm::vec3 barycentric    = glm::vec3(0.0f);
        float     sq_distance    = 0.0f;
    };

    explicit TriangleBvh(db::Mesh const& mesh);

    std::optional<ClosestResult> closest_point(glm::vec3 point) const;

    bool empty() const;

private:
    struct Aabb {
        glm::vec3 min   = glm::vec3(0.0f);
        glm::vec3 max   = glm::vec3(0.0f);
        bool      valid = false;

        void  expand(glm::vec3 point);
        void  expand(Aabb const& box);
        float distance_squared(glm::vec3 point) const;
        int   widest_axis() const;
    };

    struct Node {
        Aabb   bounds;
        size_t begin = 0;
        size_t count = 0;
        int    left  = -1;
        int    right = -1;

        bool leaf() const;
    };

    db::Mesh const*     m_mesh = nullptr;
    std::vector<size_t> m_indices;
    std::vector<Node>   m_nodes;

    Aabb      triangle_bounds(size_t triangle_index) const;
    glm::vec3 triangle_centroid(size_t triangle_index) const;
    int       build_node(size_t begin, size_t count);

    std::optional<ClosestResult>
    closest_in_node(int node_index, glm::vec3 point, float& best_sq) const;

    ClosestResult closest_on_triangle(size_t    triangle_index,
                                      glm::vec3 point) const;
};

} // namespace analysis

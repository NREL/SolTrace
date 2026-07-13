#ifndef SOLTRACE_NATIVE_RUNNER_BBOX_CALCULATOR_H
#define SOLTRACE_NATIVE_RUNNER_BBOX_CALCULATOR_H

#include "native_runner_types.hpp"

namespace SolTrace::NativeRunner
{

    /**
     * Compute the global axis-aligned bounding box for an element.
     *
     * The element's aperture and surface bounding boxes are transformed from
     * element-local coordinates to global coordinates via the element and stage
     * rotation matrices and origins.
     *
     * @param el         Element to compute bounds for (must have parent_stage set)
     * @param min_global Output: lower bound of the global AABB
     * @param max_global Output: upper bound of the global AABB
     * @return true on success
     */
    bool compute_element_bounds(const TElement *el,
                                glm::dvec3 &min_global,
                                glm::dvec3 &max_global);

} // namespace SolTrace::NativeRunner

#endif

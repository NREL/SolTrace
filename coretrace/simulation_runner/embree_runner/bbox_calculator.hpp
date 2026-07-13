#ifndef SOLTRACE_BBOX_CALCULATOR_H
#define SOLTRACE_BBOX_CALCULATOR_H

#include <native_runner_types.hpp>

namespace SolTrace::EmbreeRunner
{
    /**
     * Read the precomputed global AABB from a TElement.
     * Bounds are computed at element construction time by the native runner.
     */
    bool get_bounds(const SolTrace::NativeRunner::TElement *st_element,
                    glm::vec3& min_coord_global,
                    glm::vec3& max_coord_global);

} // namespace SolTrace::EmbreeRunner

#endif

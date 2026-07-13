#include "bbox_calculator.hpp"

#include <native_runner_types.hpp>

namespace SolTrace::EmbreeRunner
{
    using SolTrace::NativeRunner::TElement;

    bool get_bounds(const TElement *st_element,
                    glm::vec3 &min_coord_global,
                    glm::vec3 &max_coord_global)
    {
        // Bounds are precomputed in make_telement() via compute_element_bounds().
        min_coord_global = glm::vec3(st_element->BBoxMin);
        max_coord_global = glm::vec3(st_element->BBoxMax);
        return true;
    }

} // namespace SolTrace::EmbreeRunner

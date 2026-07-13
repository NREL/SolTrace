#include "bbox_calculator.hpp"

#include <glm/gtc/type_ptr.hpp>

namespace SolTrace::NativeRunner
{

static void transform_to_global(const glm::dvec3 &coord_element,
                                 const tstage_ptr &st_stage,
                                 const TElement *el,
                                 glm::dvec3 &coord_global)
{
    glm::dvec3 pos_stage = el->RLocToRef * coord_element + el->Origin;
    coord_global = st_stage->RLocToRef * pos_stage + st_stage->Origin;
}

static void transform_bounds(const glm::dvec3 &min_elem,
                              const glm::dvec3 &max_elem,
                              const tstage_ptr &st_stage,
                              const TElement *el,
                              glm::dvec3 &min_global,
                              glm::dvec3 &max_global)
{
    glm::dvec3 corners[8] = {
        {min_elem.x, min_elem.y, min_elem.z},
        {min_elem.x, min_elem.y, max_elem.z},
        {min_elem.x, max_elem.y, min_elem.z},
        {min_elem.x, max_elem.y, max_elem.z},
        {max_elem.x, min_elem.y, min_elem.z},
        {max_elem.x, min_elem.y, max_elem.z},
        {max_elem.x, max_elem.y, min_elem.z},
        {max_elem.x, max_elem.y, max_elem.z},
    };

    glm::dvec3 corners_global[8];
    for (int i = 0; i < 8; ++i)
        transform_to_global(corners[i], st_stage, el, corners_global[i]);

    min_global = corners_global[0];
    max_global = corners_global[0];
    for (int i = 1; i < 8; ++i)
    {
        min_global = glm::min(corners_global[i], min_global);
        max_global = glm::max(corners_global[i], max_global);
    }
}

bool compute_element_bounds(const TElement *el,
                             glm::dvec3 &min_global,
                             glm::dvec3 &max_global)
{
    const tstage_ptr &st_stage = el->parent_stage;

    glm::dvec2 x_minmax, y_minmax, z_minmax;
    el->aperture->bounding_box(x_minmax.x, x_minmax.y, y_minmax.x, y_minmax.y);
    el->surface->bounding_box(glm::value_ptr(x_minmax), glm::value_ptr(y_minmax),
                               z_minmax.x, z_minmax.y);

    // Expand slightly to account for floating-point precision
    constexpr double eps = 1e-3;
    glm::dvec2 expand = {-eps, eps};
    x_minmax += expand;
    y_minmax += expand;
    z_minmax += expand;

    glm::dvec3 min_elem = {x_minmax.x, y_minmax.x, z_minmax.x};
    glm::dvec3 max_elem = {x_minmax.y, y_minmax.y, z_minmax.y};

    transform_bounds(min_elem, max_elem, st_stage, el, min_global, max_global);
    return true;
}

} // namespace SolTrace::NativeRunner

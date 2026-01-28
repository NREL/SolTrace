#ifndef SOLTRACE_BBOX_CALCULATOR_H
#define SOLTRACE_BBOX_CALCULATOR_H

// #include "types.h"

#include <native_runner_types.hpp>

namespace SolTrace::EmbreeRunner
{
    bool get_bounds(const SolTrace::NativeRunner::TElement *st_element,
                    float (&min_coord_global)[3],
                    float (&max_coord_global)[3]);

} // namespace SolTrace::EmbreeRunner

#endif

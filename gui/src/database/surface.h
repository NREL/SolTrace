#pragma once

#include <array>
#include <cstdint>
#include <optional>

#include <glm/vec3.hpp>

#include "aperture.hpp"
#include "surface.hpp"

#include "mesh.h"

namespace SD = SolTrace::Data;


namespace db {

struct SurfaceGenerationOptions {
    glm::uvec2 height_field_resolution       = { 24, 24 };
    uint32_t   radial_subdivisions           = 24;
    uint32_t   perimeter_subdivisions        = 64;
    uint32_t   cylinder_angular_subdivisions = 64;
    uint32_t   cylinder_length_subdivisions  = 24;
    bool       add_thickness                 = false;
    double     thickness                     = 0.01;

    // Fidelity goes from 1 to 10
    static SurfaceGenerationOptions from_resolution_and_thickness(unsigned,
                                                                  float);
};

std::optional<Mesh>
generate_surface(SD::surface_ptr const&          surface,
                 SD::aperture_ptr const&         aperture,
                 SurfaceGenerationOptions const& options = {});

} // namespace db

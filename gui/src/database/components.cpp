#include "components.h"
#include <QtGui/qmatrix4x4.h>
#include <stdexcept>

namespace db {

glm::dmat4 TransformComponent::as_matrix() const {
    glm::dmat4 m { 1.0 };

    m = glm::translate(m, position);
    m = m * glm::mat4_cast(rotation);

    return m;
}

glm::dmat4 GlobalTransformComponent::as_matrix() const {
    glm::dmat4 m { 1.0 };

    m = glm::translate(m, position);
    m = m * glm::mat4_cast(rotation);

    return m;
}

GlobalTransformComponent
GlobalTransformComponent::compute_for(entt::registry const& reg,
                                      entt::entity          entity) {

    GlobalTransformComponent out;
    out.position = glm::dvec3 { 0.0 };
    out.rotation = glm::dquat { 1.0, 0.0, 0.0, 0.0 };

    entt::entity current = entity;

    while (current != entt::null) {
        if (auto* t = reg.try_get<TransformComponent>(current)) {
            // apply current local, then whatever accumulated so
            // far.
            out.position = t->position + t->rotation * out.position;
            out.rotation = t->rotation * out.rotation;
        }

        if (auto* child_of = reg.try_get<ChildOfComponent>(current)) {
            current = child_of->parent;

            if (auto* global = reg.try_get<GlobalTransformComponent>(current)) {
                out.position = global->position + global->rotation * out.position;
                out.rotation = global->rotation * out.rotation;
                break;
            }
        } else {
            break;
        }
    }

    return out;
}

// --- Aperture value equality -----------------------------------------------

static bool operator==(SD::Annulus const& a, SD::Annulus const& b) {
    return std::tie(a.inner_radius, a.outer_radius, a.arc_angle) ==
           std::tie(b.inner_radius, b.outer_radius, b.arc_angle);
}
static bool operator==(SD::Circle const& a, SD::Circle const& b) {
    return a.diameter == b.diameter;
}
static bool operator==(SD::Hexagon const& a, SD::Hexagon const& b) {
    return a.circumscribe_diameter == b.circumscribe_diameter;
}
static bool operator==(SD::Rectangle const& a, SD::Rectangle const& b) {
    return a.x_length() == b.x_length() && a.y_length() == b.y_length() &&
           a.x_coord() == b.x_coord() && a.y_coord() == b.y_coord();
}
static bool operator==(SD::EquilateralTriangle const& a,
                       SD::EquilateralTriangle const& b) {
    return a.circumscribe_diameter == b.circumscribe_diameter;
}
static bool operator==(SD::IrregularTriangle const& a,
                       SD::IrregularTriangle const& b) {
    return std::tie(a.x1, a.y1, a.x2, a.y2, a.x3, a.y3) ==
           std::tie(b.x1, b.y1, b.x2, b.y2, b.x3, b.y3);
}
static bool operator==(SD::IrregularQuadrilateral const& a,
                       SD::IrregularQuadrilateral const& b) {
    return std::tie(a.x1, a.y1, a.x2, a.y2, a.x3, a.y3, a.x4, a.y4) ==
           std::tie(b.x1, b.y1, b.x2, b.y2, b.x3, b.y3, b.x4, b.y4);
}

static bool is_equal(SD::Aperture const& a, SD::Aperture const& b) {
    if (a.my_type != b.my_type) return false;

    switch (a.my_type) {

    case SolTrace::Data::ANNULUS:
        return *dynamic_cast<SD::Annulus const*>(&a) ==
               *dynamic_cast<SD::Annulus const*>(&b);
    case SolTrace::Data::CIRCLE:
        return *dynamic_cast<SD::Circle const*>(&a) ==
               *dynamic_cast<SD::Circle const*>(&b);
    case SolTrace::Data::HEXAGON:
        return *dynamic_cast<SD::Hexagon const*>(&a) ==
               *dynamic_cast<SD::Hexagon const*>(&b);
    case SolTrace::Data::RECTANGLE:
        return *dynamic_cast<SD::Rectangle const*>(&a) ==
               *dynamic_cast<SD::Rectangle const*>(&b);
    case SolTrace::Data::EQUILATERAL_TRIANGLE:
        return *dynamic_cast<SD::EquilateralTriangle const*>(&a) ==
               *dynamic_cast<SD::EquilateralTriangle const*>(&b);
    case SolTrace::Data::SINGLE_AXIS_CURVATURE_SECTION:
        // TODO: When that class is implemented, we can implement this
        return false;
    case SolTrace::Data::IRREGULAR_TRIANGLE:
        return *dynamic_cast<SD::IrregularTriangle const*>(&a) ==
               *dynamic_cast<SD::IrregularTriangle const*>(&b);
    case SolTrace::Data::IRREGULAR_QUADRILATERAL:
        return *dynamic_cast<SD::IrregularQuadrilateral const*>(&a) ==
               *dynamic_cast<SD::IrregularQuadrilateral const*>(&b);
    case SolTrace::Data::APERTURE_UNKNOWN: return false;
    }

    return false;
}

static bool is_equal(SD::aperture_ptr const& a, SD::aperture_ptr const& b) {
    if (!a || !b) return a == b;
    return is_equal(*a, *b);
}

// --- Surface value equality -------------------------------------------------

static bool operator==(SD::Cone const& a, SD::Cone const& b) {
    return a.half_angle == b.half_angle;
}
static bool operator==(SD::Cylinder const& a, SD::Cylinder const& b) {
    return a.radius == b.radius;
}
static bool operator==(SD::Flat const&, SD::Flat const&) {
    return true;
}
static bool operator==(SD::Parabola const& a, SD::Parabola const& b) {
    return std::tie(a.focal_length_x, a.focal_length_y) ==
           std::tie(b.focal_length_x, b.focal_length_y);
}
static bool operator==(SD::Sphere const& a, SD::Sphere const& b) {
    return a.vertex_curv == b.vertex_curv;
}

static bool is_equal(SD::Surface const& a, SD::Surface const& b) {
    if (a.my_type != b.my_type) return false;

    switch (a.my_type) {
    case SolTrace::Data::CONE:
        return *dynamic_cast<SD::Cone const*>(&a) ==
               *dynamic_cast<SD::Cone const*>(&b);
    case SolTrace::Data::CYLINDER:
        return *dynamic_cast<SD::Cylinder const*>(&a) ==
               *dynamic_cast<SD::Cylinder const*>(&b);
    case SolTrace::Data::FLAT:
        return *dynamic_cast<SD::Flat const*>(&a) ==
               *dynamic_cast<SD::Flat const*>(&b);
    case SolTrace::Data::PARABOLA:
        return *dynamic_cast<SD::Parabola const*>(&a) ==
               *dynamic_cast<SD::Parabola const*>(&b);
    case SolTrace::Data::SPHERE:
        return *dynamic_cast<SD::Sphere const*>(&a) ==
               *dynamic_cast<SD::Sphere const*>(&b);
    case SolTrace::Data::HYPER:
    case SolTrace::Data::GENERAL_SPENCER_MURTY:
    case SolTrace::Data::TORUS:
    case SolTrace::Data::SURFACE_UNKNOWN: return false;
    }

    return false;
}

static bool is_equal(SD::surface_ptr const& a, SD::surface_ptr const& b) {
    if (!a || !b) return a == b;
    return is_equal(*a, *b);
}

bool MaterialComponent::operator==(db::MaterialComponent const& b) const {
    return optics == b.optics;
}

RaySourceResource RaySourceResource::clone() const {
    if (!source) {
        return { .source = {}, .type = type };
    }

    if (auto sun = std::dynamic_pointer_cast<SD::Sun>(source)) {
        std::vector<double> user_angle;
        std::vector<double> user_intensity;
        sun->get_user_data(user_angle, user_intensity);

        auto copy = SD::make_ray_source<SD::Sun>();
        copy->set_position(sun->get_position());
        copy->set_shape(sun->get_shape(),
                        sun->get_sigma(),
                        sun->get_half_width(),
                        sun->get_circumsolar_ratio(),
                        user_angle,
                        user_intensity);

        return { .source = copy, .type = type };
    }

    throw std::runtime_error("Unsupported ray source type in clone()");
}

DatabaseNameResource DatabaseNameResource::clone() const {
    return DatabaseNameResource {
        .name = this->name,
    };
}

GeometryComponent GeometryComponent::clone() const {
    GeometryComponent copy;

    if (aperture) {
        copy.aperture = aperture->make_copy();
    }

    if (surface) {
        nlohmann::ordered_json node;
        surface->write_json(node);
        copy.surface = SD::make_surface_from_json(node);
    }

    return copy;
}

bool GeometryComponent::operator==(db::GeometryComponent const& b) const {
    return is_equal(aperture, b.aperture) and is_equal(surface, b.surface);
}


} // namespace db

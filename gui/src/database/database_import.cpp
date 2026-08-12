#include "database/database.h"
#include "conversion.h"
#include "database/components.h"
#include "database/database_notification.h"
#include "utilities/math_utility.h"

#include <algorithm>
#include <unordered_map>

#define GLM_ENABLE_EXPERIMENTAL 1
#include <glm/gtx/io.hpp>

#include <QDebug>

/// Classic workable hash combiner
template <class T>
void hash_combine(size_t& seed, T const& v) {
    seed ^= std::hash<T>()(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}

template <>
struct std::hash<db::GeometryComponent> {
    std::size_t operator()(db::GeometryComponent const& a) const {
        size_t seed = 0;

        if (a.aperture) hash_combine(seed, *a.aperture);
        if (a.surface) hash_combine(seed, *a.surface);

        return seed;
    }
};

static size_t hash_aperture(SD::Aperture const& a);
static size_t hash_surface(SD::Surface const& a);

template <>
struct std::hash<SD::Aperture> {
    std::size_t operator()(SD::Aperture const& a) const {
        return hash_aperture(a);
    }
};

template <>
struct std::hash<SD::Surface> {
    std::size_t operator()(SD::Surface const& a) const {
        return hash_surface(a);
    }
};

// --- Hashing ----------------------------------------------------------------

static size_t hash_aperture(SD::Aperture const& a) {
    size_t seed = 0;
    hash_combine(seed, a.my_type);

    switch (a.my_type) {
    case SolTrace::Data::ANNULUS: {
        auto* aa = dynamic_cast<SD::Annulus const*>(&a);
        if (!aa) break;
        hash_combine(seed, aa->inner_radius);
        hash_combine(seed, aa->outer_radius);
        hash_combine(seed, aa->arc_angle);
        break;
    }
    case SolTrace::Data::CIRCLE: {
        auto* aa = dynamic_cast<SD::Circle const*>(&a);
        if (!aa) break;
        hash_combine(seed, aa->diameter);
        break;
    }
    case SolTrace::Data::HEXAGON: {
        auto* aa = dynamic_cast<SD::Hexagon const*>(&a);
        if (!aa) break;
        hash_combine(seed, aa->circumscribe_diameter);
        break;
    }
    case SolTrace::Data::RECTANGLE: {
        auto* aa = dynamic_cast<SD::Rectangle const*>(&a);
        if (!aa) break;
        hash_combine(seed, aa->x_length());
        hash_combine(seed, aa->y_length());
        hash_combine(seed, aa->x_coord());
        hash_combine(seed, aa->y_coord());
        break;
    }
    case SolTrace::Data::EQUILATERAL_TRIANGLE: {
        auto* aa = dynamic_cast<SD::EquilateralTriangle const*>(&a);
        if (!aa) break;
        hash_combine(seed, aa->circumscribe_diameter);
        break;
    }
    case SolTrace::Data::SINGLE_AXIS_CURVATURE_SECTION: break;
    case SolTrace::Data::IRREGULAR_TRIANGLE: {
        auto* aa = dynamic_cast<SD::IrregularTriangle const*>(&a);
        if (!aa) break;
        hash_combine(seed, aa->x1);
        hash_combine(seed, aa->y1);
        hash_combine(seed, aa->x2);
        hash_combine(seed, aa->y2);
        hash_combine(seed, aa->x3);
        hash_combine(seed, aa->y3);
        break;
    }
    case SolTrace::Data::IRREGULAR_QUADRILATERAL: {
        auto* aa = dynamic_cast<SD::IrregularQuadrilateral const*>(&a);
        if (!aa) break;
        hash_combine(seed, aa->x1);
        hash_combine(seed, aa->y1);
        hash_combine(seed, aa->x2);
        hash_combine(seed, aa->y2);
        hash_combine(seed, aa->x3);
        hash_combine(seed, aa->y3);
        hash_combine(seed, aa->x4);
        hash_combine(seed, aa->y4);
        break;
    }
    case SolTrace::Data::APERTURE_UNKNOWN: break;
    }

    return seed;
}

static size_t hash_surface(SD::Surface const& a) {
    size_t seed = 0;
    hash_combine(seed, a.my_type);

    switch (a.my_type) {
    case SolTrace::Data::CONE: {
        auto* aa = dynamic_cast<SD::Cone const*>(&a);
        if (!aa) break;
        hash_combine(seed, aa->half_angle);
        break;
    }
    case SolTrace::Data::CYLINDER: {
        auto* aa = dynamic_cast<SD::Cylinder const*>(&a);
        if (!aa) break;
        hash_combine(seed, aa->radius);
        break;
    }
    case SolTrace::Data::FLAT: break;
    case SolTrace::Data::PARABOLA: {
        auto* aa = dynamic_cast<SD::Parabola const*>(&a);
        if (!aa) break;
        hash_combine(seed, aa->focal_length_x);
        hash_combine(seed, aa->focal_length_y);
        break;
    }
    case SolTrace::Data::SPHERE: {
        auto* aa = dynamic_cast<SD::Sphere const*>(&a);
        if (!aa) break;
        hash_combine(seed, aa->vertex_curv);
        break;
    }
    case SolTrace::Data::HYPER:
    case SolTrace::Data::GENERAL_SPENCER_MURTY:
    case SolTrace::Data::TORUS:
    case SolTrace::Data::SURFACE_UNKNOWN: break;
    }

    return seed;
}

namespace db {

TransformComponent extract_tf(SD::Element const& e) {
    auto pos = e.get_origin_ref();

    auto aim = e.get_aim_vector_ref();

    aim = glm::normalize(aim - pos);

    auto quat = dir_roll_to_quat(aim, e.get_zrot_radians());
    (void)quat;

    // TODO: Verify transform conversion against SolTrace coordinate conventions.
    auto quat2 = glm::quat_cast(e.get_local_to_reference());

    return TransformComponent { .position = pos, .rotation = quat2 };
}

TransformComponent extract_tf_stage(SD::Element const& e) {
    auto pos = e.get_origin_stage();

    auto aim = e.get_aim_vector_stage();

    aim = glm::normalize(aim - pos);

    auto quat = dir_roll_to_quat(aim, e.get_zrot_radians());
    (void)quat;

    auto quat2 = glm::quat_cast(e.get_local_to_reference());

    return TransformComponent { .position = pos, .rotation = quat2 };
}

// =============================================================================

struct MaterialMap {
    std::unordered_map<std::shared_ptr<const SD::OpticalPropertySet>,
                       entt::entity>
        map;
};

static void import_optics(
    Database&                                            reg,
    entt::entity                                         entity,
    SD::Element const&                                   item,
    MaterialMap&                                         material_groups,
    std::unordered_map<GeometryComponent, entt::entity>& geometry_groups,
    size_t&                                              group_counter) {

    auto& registry = reg.as_registry();

    auto property_sptr = item.get_optical_property_set();

    // If it has no optical data, skip importing optics.
    if (!property_sptr and !item.get_aperture() and !item.get_surface()) {
        // does not have geometry.
        qDebug() << "Skipping optics on" << entt::to_integral(entity);
        return;
    }

    // It should have all or nothing
    if (!property_sptr or !item.get_aperture() or !item.get_surface()) {

        emplace_patch<ImportErrorComponent>(reg, entity, [](auto& c) {
            c.reason += "Missing aperture and surface";
        });
        qWarning() << "Entity" << entt::to_integral(entity)
                   << "missing surface properties";
        return;
    }

    // auto local_mat = MaterialComponent {
    //     .optics_front = *item.get_front_optical_properties(),
    //     .optics_back  = *item.get_back_optical_properties(),
    // };

    auto local_geom = GeometryComponent {
        .aperture = item.get_aperture(),
        .surface  = item.get_surface(),
    };

    if (auto group = material_groups.map.find(property_sptr);
        group != material_groups.map.end()) {
        // we have such a group

        reg.assign_material(entity, group->second);

    } else {
        // no such group

        auto new_group_params = MaterialComponent { *property_sptr };

        auto new_group = MaterialGroupComponent { };

        auto group_entity = reg.create();

        registry.emplace<MaterialGroupComponent>(group_entity, new_group);
        registry.emplace<MaterialComponent>(group_entity, new_group_params);
        auto material_name = QString::fromStdString(property_sptr->get_name());
        if (material_name.isEmpty()) {
            material_name = QString("Material %1").arg(group_counter);
        }

        registry.emplace<IdentityComponent>(
            group_entity, IdentityComponent { .name = material_name });

        group_counter++;

        material_groups.map.try_emplace(property_sptr, group_entity);

        reg.assign_material(entity, group_entity);
    }

    if (auto group = geometry_groups.find(local_geom);
        group != geometry_groups.end()) {
        // we have such a group

        reg.assign_geometry(entity, group->second);

    } else {
        // no such group

        auto new_group_params = local_geom;

        auto new_group = GeometryGroupComponent { };

        auto group_entity = reg.create();

        registry.emplace<GeometryGroupComponent>(group_entity, new_group);
        registry.emplace<GeometryComponent>(group_entity, new_group_params);
        registry.emplace<IdentityComponent>(
            group_entity,
            IdentityComponent {
                .name = QString("Geometry %1").arg(group_counter),
            });

        group_counter++;

        geometry_groups.try_emplace(local_geom, group_entity);

        reg.assign_geometry(entity, group_entity);
    }
}

struct StageComponent {
    TransformComponent stage_tf;
    TransformComponent this_in_stage;
};

void Database::import(SD::SimulationData& data, bool legacy_import) {

    // Assuming we are not re-using registries, which we are not for the moment
    auto imported_tag = create_tag("imported");

    // we use pointers as element IDs are scoped to global AND composite
    std::unordered_map<SD::Element*, entt::entity> element_to_entity;

    auto get_or_create_entity = [&](SD::Element* ptr) -> entt::entity {
        if (!ptr) return entt::null;

        if (auto eiter = element_to_entity.find(ptr);
            eiter != element_to_entity.end()) {
            return eiter->second;
        } else {
            auto ent               = m_registry.create();
            element_to_entity[ptr] = ent;

            m_registry.emplace<ElementComponent>(ent);

            assign_tag(ent, imported_tag);

            return ent;
        }
    };

    // MAP MAY NOT BE RE-USED!
    MaterialMap                                         material_groups;
    std::unordered_map<GeometryComponent, entt::entity> geometry_groups;

    size_t group_counter = 0;

    for (auto iter = data.get_iterator(); !data.is_at_end(iter); ++iter) {
        auto const& element = *(iter->second);

        if (element.is_stage()) {
            // we want to avoid materializing stage elements
            auto c = std::dynamic_pointer_cast<SD::StageElement>(iter->second);

            TransformComponent stage_tf = extract_tf(element);

            for (auto iter = c->get_const_iterator(); !c->is_at_end(iter);
                 ++iter) {

                auto child_ent = get_or_create_entity(iter->second.get());

                m_registry.emplace_or_replace<StageComponent>(
                    child_ent,
                    StageComponent {
                        .stage_tf      = stage_tf,
                        .this_in_stage = extract_tf_stage(*iter->second),
                    });
            }

            // we do NOT add these as children, we instead make them globals
            continue;
        }

        entt::entity ent = get_or_create_entity(iter->second.get());

        if (!m_registry.valid(ent)) {
            qWarning() << "Unable to mirror element" << &element;
        }

        // auto position = element.get_origin_ref();

        // auto rotation = convert(element.get_local_to_reference());

        m_registry.emplace<TransformComponent>(ent, extract_tf(element));

        if (!element.get_name().empty()) {

            // for legacy, sometimes the name is just a number

            auto name = QString::fromStdString(element.get_name());

            bool ok = false;

            name.toLong(&ok);

            if (legacy_import && ok) {
                // Legacy numeric names are not descriptive; add context when possible.

                auto id = element.get_id();

                name = QStringLiteral("Element %1").arg(id);
            }

            m_registry.emplace<IdentityComponent>(ent, name);
        }

        if (!element.is_enabled()) {
            m_registry.emplace<DisabledComponent>(ent);
        }

        if (element.is_virtual()) {
            qDebug() << "Adding virtual tag to " << ent;
            m_registry.emplace<VirtualTagComponent>(ent);
        }

        if (element.is_composite()) {
            auto& c = *static_cast<SD::CompositeElement const*>(&element);

            for (auto iter = c.get_const_iterator(); !c.is_at_end(iter);
                 ++iter) {

                auto child_ent = get_or_create_entity(iter->second.get());

                if (m_registry.any_of<ChildOfComponent>(child_ent)) {
                    // Multiple parents indicate invalid hierarchy input.
                    throw std::runtime_error("Multiple parents for element");
                }

                set_parent(child_ent, ent);
            }
        }

        import_optics(*this,
                      ent,
                      element,
                      material_groups,
                      geometry_groups,
                      group_counter);
    }

    {
        // burn in globals from a stage

        auto view = m_registry.view<StageComponent const>();

        for (auto const& [e, stage] : view.each()) {
            auto new_pos =
                stage.stage_tf.position +
                stage.stage_tf.rotation * stage.this_in_stage.position;
            auto new_rot =
                stage.stage_tf.rotation * stage.this_in_stage.rotation;

            m_registry.emplace_or_replace<TransformComponent>(
                e,
                TransformComponent {
                    .position = new_pos,
                    .rotation = new_rot,
                });
        }

        m_registry.clear<StageComponent>();
    }

    ray_source_resource.set(RaySourceResource {
        .source = data.get_ray_source(),
    });

    simulation_parameters_resource.set(
        data.get_simulation_parameters());

    qInfo() << "Imported" << this->m_registry.view<ElementComponent>()->size()
            << "elements";

    qInfo() << "Imported"
            << this->m_registry.view<MaterialGroupComponent>()->size()
            << "materials";
    qInfo() << "Imported"
            << this->m_registry.view<GeometryGroupComponent>()->size()
            << "geometries";
}

} // namespace db

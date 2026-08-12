#include "database/database.h"
#include "conversion.h"
#include "database/components.h"
#include "database/database_clone.h"

#include "simulation_data_api.hpp"

#include <unordered_map>

#include <QDateTime>
#include <QDebug>

namespace db {

// TODO: Remove once vector aim export semantics are finalized.
static void install_transform(SD::element_ptr           ptr,
                              TransformComponent const& tf_comp) {

    glm::dvec3 aim  = { };
    double     roll = { };
    quat_to_dir_roll(tf_comp.rotation, aim, roll);

    auto origin = tf_comp.position;
    ptr->set_origin(origin.x, origin.y, origin.z);

    aim = (origin + aim * 100.0);

    ptr->set_aim_vector(aim.x, aim.y, aim.z);

    ptr->set_zrot_radians(roll);
}

static void install_group(SD::element_ptr                 ptr,
                          SD::OpticalPropertySetReference param,
                          GeometryComponent const&        geom_param) {
    ptr->set_aperture(geom_param.aperture);
    ptr->set_surface(geom_param.surface);

    ptr->set_optical_property_set(param);
}

static SD::OpticalPropertySet
with_optical_name(SD::OpticalPropertySet const& optics, QString const& name) {
    nlohmann::ordered_json node;
    optics.write_json(node);
    node["my_name"] = name.toStdString();
    return SD::OpticalPropertySet(node);
}

Result<std::shared_ptr<DatabaseExport>, QString> Database::export_to_simdata() {
    SD::SimulationData ret;

    auto param_ptr = simulation_parameters_resource.get();
    if (!param_ptr) {
        return QStringLiteral("Internal error: simulation is missing "
                              "simulation parameters resource.");
    }

    ret.get_simulation_parameters() = *param_ptr;

    auto ray_source_ptr = ray_source_resource.get();
    if (!ray_source_ptr) {
        return QStringLiteral("Ray source is not defined for this simulation.");
    }

    ret.add_ray_source(ray_source_ptr->source);

    {
        auto view = m_registry.view<const ElementComponent>();
        if (view->size() == 0) {
            return QStringLiteral("There are no elements in this simulation.");
        }
    }

    std::unordered_map<entt::entity, SD::element_ptr> entity_element_map;

    // Mirror all elements

    {
        auto view = m_registry.view<const ElementComponent>();
        for (auto const& [e] : view.each()) {
            SD::element_ptr ptr;

            if (children_of(e).size()) {
                auto n = std::make_shared<SD::CompositeElement>();
                ptr    = n;
            } else if (m_registry.all_of<VirtualTagComponent>(e)) {
                auto n = std::make_shared<SD::VirtualElement>();
                ptr    = n;
            } else {
                auto n = std::make_shared<SD::SingleElement>();
                ptr    = n;
            }

            ptr->set_name(name_of(e).toStdString());

            entity_element_map[e] = ptr;
        }
    }

    {
        auto view = m_registry.view<const TransformComponent>();
        for (auto const& [e, tf] : view.each()) {
            install_transform(entity_element_map.at(e), tf);
        }
    }

    {
        auto view =
            m_registry.view<const ElementComponent, const IdentityComponent>();
        for (auto const& [e, ident] : view.each()) {
            entity_element_map.at(e)->set_name(ident.name.toStdString());
        }
    }

    {
        auto view = m_registry.view<const DisabledComponent>();
        for (auto const& [e] : view.each()) {
            entity_element_map.at(e)->disable();
        }
    }

    std::unordered_map<entt::entity, SD::OpticalPropertySetReference> prop_map;

    {
        auto view = m_registry.view<const MaterialGroupComponent,
                                    const MaterialComponent>();
        for (auto const& [e, mat_members, mat_data] : view.each()) {
            auto named_optics =
                with_optical_name(mat_data.optics, name_of(e));
            auto ref = ret.find_or_add_optical_property_set(named_optics);

            prop_map.try_emplace(e, ref);
        }
    }

    {
        auto view = m_registry.view<const MaterialGroupMemberComponent,
                                    const GeometryGroupMemberComponent>();
        for (auto const& [e, mat, geom] : view.each()) {

            // get group, we assume this is valid
            auto const& mat_group = prop_map.at(mat.group);
            auto const& geom_group =
                m_registry.get<GeometryComponent>(geom.group);

            auto element = entity_element_map.at(e);

            if (auto ptr = m_registry.try_get<ChildrenComponent>(e); ptr) {
                if (!ptr->children.empty()) {
                    // it has children AND a group component. add a proxy that
                    // will only hold the geometry

                    auto composite =
                        dynamic_cast<SD::CompositeElement*>(element.get());
                    Q_ASSERT(composite);

                    auto n = std::make_shared<SD::SingleElement>();
                    composite->add_element(n);
                    install_group(n, mat_group, geom_group);

                    continue;
                }
            }

            // has no children

            install_group(element, mat_group, geom_group);
        }
    }

    {
        auto view = m_registry.view<const ChildrenComponent>();
        for (auto const& [e, children] : view.each()) {

            auto composite = dynamic_cast<SD::CompositeElement*>(
                entity_element_map.at(e).get());
            Q_ASSERT(composite);

            for (auto child : children.children) {

                auto iter = entity_element_map.find(child);

                if (iter == entity_element_map.end()) {
                    qCritical() << "Missing child element";
                    continue;
                }

                if (iter->second->is_composite()) {
                    qCritical() << "Composite child under non-stage composite "
                                   "is not supported";
                    continue;
                }

                auto ret = composite->add_element(iter->second);

                if (ret < 0) { qCritical() << "Add failed"; }
            }
        }
    }

    // Install all elements into the sim engine
    for (auto const& iter : entity_element_map) {
        if (m_registry.any_of<ChildOfComponent>(iter.first)) {
            // Only add top-level elements; composites will add subelements.
            continue;
        }
        auto const& ptr = iter.second;
        try {
            ret.add_element(iter.second);
        } catch (std::exception const& e) {
            qCritical() << "Unable to export element"
                        << entt::to_integral(iter.first)
                        << iter.second->get_name() << e.what();

            return QStringLiteral(
                "Internal error: an element could not be exported.");
        }
    }

    std::unordered_map<SD::element_id, entt::entity> entity_rev_map;

    for (auto iter = entity_element_map.begin();
         iter != entity_element_map.end();
         ++iter) {
        entity_rev_map[iter->second->get_id()] = iter->first;
    }

    auto new_name = QString("%1 v%2")
                        .arg(this->name())
                        .arg(QDateTime::currentSecsSinceEpoch());

    auto clone_result =
        clone_database_with_entity_map(*this, new_name, nullptr);

    for (auto& [element_id, entity] : entity_rev_map) {
        auto iter = clone_result.old_to_new_map.find(entity);
        if (iter != clone_result.old_to_new_map.end()) {
            entity = iter->second;
        } else {
            entity = entt::null;
        }
    }

    DatabaseExport export_ret;
    export_ret.data = std::make_shared<SD::SimulationData>(std::move(ret));
    export_ret.element_map = std::move(entity_rev_map);
    export_ret.source_database = std::move(clone_result.database);

    return std::make_shared<DatabaseExport>(std::move(export_ret));
}

} // namespace db

#pragma once

#include <entt/entt.hpp>

#include "data/components.h"
#include "data_io/database_export.h"
#include "data_io/database_notification.h"
#include "entity.h"
#include "magic_enum/magic_enum.hpp"
#include "support/result.h"

#include <QDebug>
#include <QPointer>
#include <QStringListModel>
#include <QtTypes>
#include <qqmlintegration.h>

namespace SolTrace::GUI::Data {


/// Helper function, find a corresponding key to a value in a map
/// Slow, but OK for UI work
template <class K, class V>
std::optional<K> reverse_lookup(std::map<K, V> const& map, V const& value) {
    for (auto const& [k, v] : map) {
        if (v == value) return k;
    }
    return std::nullopt;
}

/// Write all enum options to a string list model
template <class K>
void build_options(QStringListModel& dest) {
    QStringList items;

    for (auto const& iter : magic_enum::enum_entries<K>()) {
        auto val = QString(iter.second.data());

        items.push_back(val);
    }

    dest.setStringList(items);
}


/// A simulation scene
class Database : public QObject {
    Q_OBJECT
    entt::registry m_registry;

    Q_PROPERTY(QString name READ name WRITE set_name NOTIFY name_changed FINAL)

public:
    /// Create a new simulation database
    explicit Database(QString database_name, QObject* p = nullptr);

    virtual ~Database();

    /// Duplicate this database with an optional QObject parent
    Database* clone(QString new_database_name, QObject* p = nullptr) const;

    /// Merge in simulation data. Note, this should be called closely after
    /// the database constructor. We have this split here so that we can
    /// allocate a database on one thread and fill it in another.
    /// Thus: DO NOT DO QObject THINGS IN THIS FUNCTION, only fill the reg.
    void import(SD::SimulationData&, bool legacy_import = false);

    /// Convert a database back into a Soltrace dataset
    Support::Result<std::shared_ptr<DatabaseExport>, QString>
    export_to_simdata();

    QString name() const;
    void    set_name(QString);

public:
    /// Get underlying registry of this database

    operator entt::registry&();
    operator entt::registry const&() const;

    entt::registry&       as_registry();
    entt::registry const& as_registry() const;

public:
    ComponentAPIUpdate<IdentityComponent>   identity;
    ComponentAPIUpdate<TransformComponent>  transform;
    ComponentAPI<GlobalTransformComponent>  global_transform;
    ComponentAPIUpdate<InvisibleComponent>  invisible;
    ComponentAPIUpdate<DisabledComponent>   disabled;
    ComponentAPIUpdate<VirtualTagComponent> virtual_tag;
    ComponentAPI<ChildOfComponent>          parent;
    ComponentAPI<TagComponent>              tag_root;

    ComponentAPI<ElementComponent> element_tag;

    ComponentAPI<MaterialGroupComponent>       material_root;
    ComponentAPIUpdate<MaterialComponent>      material_parameters;
    ComponentAPI<MaterialGroupMemberComponent> material_group_membership;

    ComponentAPI<GeometryGroupComponent>       geometry_root;
    ComponentAPIUpdate<GeometryComponent>      geometry_parameters;
    ComponentAPI<GeometryGroupMemberComponent> geometry_group_membership;

    ComponentAPI<ChildrenComponent> children;

    ComponentAPI<TagMembershipComponent> tag_membership;

    ComponentAPIUpdate<SelectedComponent> selected;

    ComponentAPIUpdate<ColorComponent> color;

    ComponentAPI<HasFluxMapComponent> flux_map;

    SingletonComponentAPI<DatabaseNameResource> database_name_resource;
    SingletonComponentAPI<RaySourceResource>    ray_source_resource;
    SingletonComponentAPI<SD::SimulationParameters>
        simulation_parameters_resource;

public:
    /// Helper function: patch a component, creating it if it does not exist.
    template <class Component, class Function>
    void emplace_patch(entt::entity entity, Function&& f) {
        emplace_patch(m_registry, entity, f);
    }

public:
    entt::entity create();

    bool valid(entt::entity) const;

    /// Clear the active parent of an entity
    void unset_parent(entt::entity child);

    /// Set the parent of an entity
    void set_parent(entt::entity child, entt::entity parent);

    /// Get the list of children of this entity. Returns an empty list if there
    /// are none.
    std::span<Entity const> children_of(entt::entity parent) const;

    /// Get the parent of this entity. Returns entt::null if there is none.
    entt::entity parent_of(entt::entity child) const;

    /// Assign an entity to a material group
    void assign_material(entt::entity child, entt::entity group);

    /// Remove an entity from a material group
    void remove_material(entt::entity child);

    /// Assign an entity to a geometry group
    void assign_geometry(entt::entity child, entt::entity group);

    /// Remove an entity from a geometry group
    void remove_geometry(entt::entity child);

    /// Create a new tag. Note that tag names should be unique.

    QString sanitize_tag_name(QString);

    entt::entity create_tag(QString name);

    /// Ask if an entity has been given a tag
    bool is_tagged(Entity item, Entity tag) const;

    /// Assign an entity to a specific tag
    void assign_tag(Entity item, Entity tag);

    /// Remove a tag from an entity
    void unassign_tag(entt::entity item, entt::entity tag);

    /// Clear and destroy a specific tag
    void delete_tag(entt::entity tag);

    /// Get all the tags for an entity
    std::span<Entity const> tags_for(entt::entity item) const;

    /// Get the global ray source of the database
    SD::ray_source_ptr get_ray_source() const;

    /// Get the global simulation parameters
    SD::SimulationParameters const& get_sim_params() const;

public slots:
    /// Get the name of an entity, either using the Identity component, or by
    /// using the entity ID.
    QString name_of(SolTrace::GUI::Data::Entity item) const;

    void set_name_of(SolTrace::GUI::Data::Entity item, QString new_name);

    QString sanitize_element_name(QString);
    QString sanitize_entity_name(QString);

    SolTrace::GUI::Data::Entity add_element(QString new_name, SolTrace::GUI::Data::Entity parent = { });

    void delete_element(SolTrace::GUI::Data::Entity to_delete);

    bool is_virtual_element(SolTrace::GUI::Data::Entity element) const;
    void set_virtual_element(SolTrace::GUI::Data::Entity element, bool is_virtual);

    /// Materials
    QString sanitize_material_name(QString);

    SolTrace::GUI::Data::Entity add_material_group(QString             new_name,
                                  QVector<SolTrace::GUI::Data::Entity> members    = { },
                                  SolTrace::GUI::Data::Entity          clone_from = { });

    size_t material_use_count(SolTrace::GUI::Data::Entity material);

    void delete_material_group(SolTrace::GUI::Data::Entity to_delete, SolTrace::GUI::Data::Entity move_to = { });

    SolTrace::GUI::Data::Entity material_of(SolTrace::GUI::Data::Entity element) const;

    /// Geometry

    QString sanitize_geometry_name(QString);

    SolTrace::GUI::Data::Entity add_geometry_group(QString             new_name,
                                  QVector<SolTrace::GUI::Data::Entity> members    = { },
                                  SolTrace::GUI::Data::Entity          clone_from = { });

    size_t geometry_use_count(SolTrace::GUI::Data::Entity geometry);

    void delete_geometry_group(SolTrace::GUI::Data::Entity to_delete, SolTrace::GUI::Data::Entity move_to = { });

    SolTrace::GUI::Data::Entity geometry_of(SolTrace::GUI::Data::Entity element) const;

    /// Selection methods
    void select(SolTrace::GUI::Data::Entity to_select);
    void add_to_selection(SolTrace::GUI::Data::Entity to_select);

    void deselect(SolTrace::GUI::Data::Entity to_deselect);

    void toggle_selection(SolTrace::GUI::Data::Entity to_toggle_selection);

    void clear_selection();

    bool is_selected(SolTrace::GUI::Data::Entity e) const;

    void select_all_with_material(SolTrace::GUI::Data::Entity);
    void select_all_with_geometry(SolTrace::GUI::Data::Entity);
    void deselect_all_with_material(SolTrace::GUI::Data::Entity);
    void deselect_all_with_geometry(SolTrace::GUI::Data::Entity);

    /// Color
    void set_color(SolTrace::GUI::Data::Entity to_color, QColor new_color);

signals:
    void bulk_selection_changed();
    void name_changed();
};

} // namespace SolTrace::GUI::Data


Q_DECLARE_METATYPE(SolTrace::GUI::Data::Entity);

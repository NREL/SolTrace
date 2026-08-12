#include "database.h"
#include "conversion.h"
#include "database/components.h"
#include "database/database_notification.h"
#include "database/database_clone.h"

#include "simulation_data_api.hpp"

#include <algorithm>
#include <optional>
#include <unordered_map>
#include <unordered_set>

#include <QDebug>
#include <QSet>

namespace db {

/// Get a global transform regardless if the entity has a local transform
static GlobalTransformComponent
compute_global_without_local_transform(entt::registry& reg,
                                       entt::entity    entity) {
    if (auto* child_of = reg.try_get<ChildOfComponent>(entity)) {
        return GlobalTransformComponent::compute_for(reg, child_of->parent);
    }

    return GlobalTransformComponent {
        .position = glm::dvec3 { 0.0 },
        .rotation = glm::dquat { 1.0, 0.0, 0.0, 0.0 },
    };
}

/// Compute global transforms for a subtree
static void update_global_transform_subtree(
    entt::registry&                         reg,
    entt::entity                            entity,
    std::optional<GlobalTransformComponent> root_transform = std::nullopt) {
    auto global = root_transform.value_or(
        GlobalTransformComponent::compute_for(reg, entity));

    reg.emplace_or_replace<GlobalTransformComponent>(entity, global);

    auto* ptr = reg.try_get<ChildrenComponent>(entity);
    if (!ptr) return;

    // Copy before recursion in case a callback mutates the child list.
    auto cpy = ptr->children;

    for (auto e : std::as_const(cpy)) {
        update_global_transform_subtree(reg, e);
    }
}

/// On local transform change, update globals
static void tf_change_callback(entt::registry& reg, entt::entity entity) {
    update_global_transform_subtree(reg, entity);
}

/// On local transform change, update globals
static void tf_destroy_callback(entt::registry& reg, entt::entity entity) {
    update_global_transform_subtree(
        reg, entity, compute_global_without_local_transform(reg, entity));
}

Database::Database(QString database_name, QObject* p)
    : QObject(p),
      m_registry(),
      identity(m_registry),
      transform(m_registry),
      global_transform(m_registry),
      invisible(m_registry),
      disabled(m_registry),
      virtual_tag(m_registry),
      parent(m_registry),
      tag_root(m_registry),
      element_tag(m_registry),
      material_root(m_registry),
      material_parameters(m_registry),
      material_group_membership(m_registry),
      geometry_root(m_registry),
      geometry_parameters(m_registry),
      geometry_group_membership(m_registry),
      children(m_registry),
      tag_membership(m_registry),
      selected(m_registry),
      color(m_registry),
      flux_map(m_registry),
      database_name_resource(m_registry),
      ray_source_resource(m_registry),
      simulation_parameters_resource(m_registry) {

    m_registry.on_construct<TransformComponent>()
        .template connect<&tf_change_callback>();
    m_registry.on_update<TransformComponent>()
        .template connect<&tf_change_callback>();
    m_registry.on_destroy<TransformComponent>()
        .template connect<&tf_destroy_callback>();

    database_name_resource.set(DatabaseNameResource { .name = database_name });

    simulation_parameters_resource.set({});

    qDebug() << Q_FUNC_INFO;
}

Database::~Database() {
    qDebug() << Q_FUNC_INFO;
}

Database* Database::clone(QString new_database_name, QObject* p) const {
    return clone_database_with_entity_map(*this, new_database_name, p)
        .database.release();
}

Database::operator entt::registry&() {
    return m_registry;
}

Database::operator const entt::registry&() const {
    return m_registry;
}

entt::registry& Database::as_registry() {
    return m_registry;
}

const entt::registry& Database::as_registry() const {
    return m_registry;
}

QString Database::name() const {
    auto ptr = database_name_resource.get();

    if (ptr) { return ptr->name; }

    return "Untitled";
}

void Database::set_name(QString s) {

    if (s == name()) { return; }

    database_name_resource.set(DatabaseNameResource { .name = s });

    emit name_changed();
}

entt::entity Database::create() {
    return m_registry.create();
}

bool Database::valid(entt::entity e) const {
    return m_registry.valid(e);
}

void Database::unset_parent(entt::entity child) {

    auto child_comp = m_registry.try_get<ChildOfComponent>(child);

    if (!child_comp) return;

    auto const parent = child_comp->parent;

    auto has_parent_comp = m_registry.all_of<ChildrenComponent>(parent);

    if (!has_parent_comp) {
        m_registry.erase<ChildOfComponent>(child);
        update_global_transform_subtree(m_registry, child);
        return;
    }

    m_registry.erase<ChildOfComponent>(child);

    m_registry.patch<ChildrenComponent>(parent, [child](ChildrenComponent& c) {
        erase(c.children, Entity(child));
    });

    update_global_transform_subtree(m_registry, child);
}

void Database::set_parent(entt::entity child, entt::entity parent) {
    if (!m_registry.valid(parent) or !m_registry.valid(child)) {
        qWarning() << "Invalid parent or child in set_parent";
        return;
    }

    // remove existing parent, if any
    unset_parent(child);

    if (!m_registry.all_of<ChildrenComponent>(parent)) {
        // If the parent has no child component, insert this as the first
        // child
        m_registry.emplace<ChildrenComponent>(parent,
                                              ChildrenComponent {
                                                  .children = { child },
                                              });
    } else {
        // it has a child component, add to it
        m_registry.patch<ChildrenComponent>(
            parent,
            [child](ChildrenComponent& a) { a.children.push_back(child); });
    }

    // set the child's parent
    m_registry.emplace<ChildOfComponent>(child,
                                         ChildOfComponent { .parent = parent });

    update_global_transform_subtree(m_registry, child);
}

std::span<Entity const> Database::children_of(entt::entity parent) const {
    auto parent_comp = m_registry.try_get<ChildrenComponent>(parent);

    if (!parent_comp) { return { }; }

    return parent_comp->children;
}

entt::entity Database::parent_of(entt::entity child) const {
    if (!m_registry.valid(child)) return entt::null;
    auto* ptr = m_registry.try_get<ChildOfComponent>(child);

    if (!ptr) return entt::null;

    return ptr->parent;
}

void Database::remove_material(entt::entity child) {
    // is this a member of a group?

    auto child_comp = m_registry.try_get<MaterialGroupMemberComponent>(child);

    if (!child_comp) return;

    // Check if the group parent exists
    auto parent_comp =
        m_registry.try_get<MaterialGroupComponent>(child_comp->group);

    // remove us from a member of the group
    m_registry.erase<MaterialGroupMemberComponent>(child);

    if (!parent_comp) { return; }

    // remove us from the parent

    m_registry.patch<MaterialGroupComponent>(
        child_comp->group, [child](MaterialGroupComponent& comp) {
            erase(comp.members, Entity(child));
        });
}

void Database::assign_material(entt::entity child, entt::entity group) {
    if (!valid(child) or !valid(group)) return;

    if (!m_registry.all_of<MaterialGroupComponent>(group)) return;

    remove_material(child);

    m_registry.emplace_or_replace<MaterialGroupMemberComponent>(
        child, MaterialGroupMemberComponent { .group = group });

    ::db::emplace_patch<MaterialGroupComponent>(
        m_registry, group, [child](MaterialGroupComponent& c) {
            c.members.push_back(child);
        });
}

void Database::remove_geometry(entt::entity child) {
    // is this a member of a group?
    auto child_comp = m_registry.try_get<GeometryGroupMemberComponent>(child);

    if (!child_comp) return;

    // Check if the group parent exists
    auto parent_comp =
        m_registry.try_get<GeometryGroupComponent>(child_comp->group);

    // remove us from a member of the group
    m_registry.erase<GeometryGroupMemberComponent>(child);

    if (!parent_comp) { return; }

    // remove us from the parent
    m_registry.patch<GeometryGroupComponent>(
        child_comp->group, [child](GeometryGroupComponent& comp) {
            erase(comp.members, Entity(child));
        });
}

void Database::assign_geometry(entt::entity child, entt::entity group) {
    if (!valid(child) or !valid(group)) return;

    if (!m_registry.all_of<GeometryGroupComponent>(group)) return;

    remove_geometry(child);

    m_registry.emplace_or_replace<GeometryGroupMemberComponent>(
        child, GeometryGroupMemberComponent { .group = group });

    ::db::emplace_patch<GeometryGroupComponent>(
        m_registry, group, [child](GeometryGroupComponent& c) {
            c.members.push_back(child);
        });
}

SD::ray_source_ptr Database::get_ray_source() const {
    if (auto ptr = ray_source_resource.get(); ptr) {
        return ptr->source;
    }

    return { };
}

SD::SimulationParameters const& Database::get_sim_params() const {
    if (auto ptr = simulation_parameters_resource.get(); ptr) {
        return *ptr;
    }

    throw std::runtime_error("missing simulation parameters");
}

template <class Component>
QString
sanitize_new_name(entt::registry& registry, QString name, QString class_type) {
    // Names are case-sensitive: "Tag" and "tag" are allowed as distinct names.

    name = name.trimmed();

    if (name.isEmpty()) { name = class_type; }

    QSet<QString> current_name_list;

    for (auto const& pack :
         registry.view<Component, IdentityComponent>().each()) {
        auto const& ident = std::get<IdentityComponent&>(pack);
        current_name_list.insert(ident.name);
    }

    auto stem    = name;
    int  counter = 1;

    // If name ends in "_N", continue from the base name.
    const int underscore = name.lastIndexOf('_');

    if (underscore > 0 && underscore < name.size() - 1) {
        bool ok     = false;
        int  number = QStringView(name).mid(underscore + 1).toInt(&ok);

        if (ok && number > 0) {
            stem    = name.left(underscore);
            counter = number + 1;
        }
    }

    while (current_name_list.contains(name)) {
        name = QString("%1_%2").arg(stem).arg(counter);
        counter++;
    }

    return name;
}

QString Database::sanitize_tag_name(QString name) {
    return sanitize_new_name<TagComponent>(m_registry, name, "Tag");
}

entt::entity Database::create_tag(QString name) {
    auto ret = m_registry.create();

    m_registry.emplace<TagComponent>(ret);
    m_registry.emplace<IdentityComponent>(ret,
                                          IdentityComponent { .name = name });

    return ret;
}

bool Database::is_tagged(Entity item, Entity tag) const {
    if (auto* ptr = m_registry.try_get<TagMembershipComponent>(item); ptr) {
        auto& t = ptr->tags;

        return std::find(t.begin(), t.end(), tag) != t.end();
    }
    return false;
}

void Database::assign_tag(Entity item, Entity tag) {
    if (!m_registry.all_of<TagComponent>(tag)) { return; }

    if (is_tagged(item, tag)) return;

    auto& storage =
        m_registry.storage<ATagMemberComponent>(entt::to_integral(tag.value));

    storage.emplace(item);

    ::db::emplace_patch<TagMembershipComponent>(
        m_registry, item, [tag](TagMembershipComponent& tc) {
            tc.tags.push_back(tag);
        });
}

void Database::unassign_tag(entt::entity item, entt::entity tag) {
    if (!m_registry.all_of<TagComponent>(tag)) { return; }

    if (!m_registry.all_of<TagMembershipComponent>(item)) { return; }

    m_registry.patch<TagMembershipComponent>(
        item,
        [tag](TagMembershipComponent& tc) { erase(tc.tags, Entity(tag)); });

    auto& storage =
        m_registry.storage<ATagMemberComponent>(entt::to_integral(tag));

    if (storage.contains(item)) { storage.erase(item); }

    // if (storage.empty()) { reg.reset(entt::to_integral(tag)); }
}

void Database::delete_tag(entt::entity tag) {
    auto& storage =
        m_registry.storage<ATagMemberComponent>(entt::to_integral(tag));

    for (auto x : storage) {

        m_registry.patch<TagMembershipComponent>(
            x, [tag](TagMembershipComponent& comp) {
                erase(comp.tags, Entity(tag));
            });
    }

    m_registry.reset(entt::to_integral(tag));

    m_registry.destroy(tag);
}

std::span<Entity const> Database::tags_for(entt::entity item) const {
    if (auto ptr = m_registry.try_get<TagMembershipComponent>(item); ptr) {
        return ptr->tags;
    }

    return { };
}

QString Database::name_of(Entity item) const {
    if (!m_registry.valid(item)) return { };

    if (auto ptr = m_registry.try_get<IdentityComponent>(item); ptr) {
        return ptr->name;
    }

    return QString("Element %1").arg(entt::to_integral(item.value));
}

void Database::set_name_of(db::Entity item, QString new_name) {
    if (!m_registry.valid(item)) return;

    m_registry.emplace_or_replace<IdentityComponent>(
        item, IdentityComponent { .name = new_name });
}

QString Database::sanitize_element_name(QString name) {
    return sanitize_new_name<ElementComponent>(m_registry, name, "Element");
}

QString Database::sanitize_entity_name(QString name) {
    return sanitize_element_name(name);
}

db::Entity Database::add_element(QString new_name, db::Entity parent) {
    auto entity = create();

    m_registry.emplace<ElementComponent>(entity);
    m_registry.emplace<TransformComponent>(
        entity,
        TransformComponent {
            .position = glm::dvec3 { 0.0 },
            .rotation = glm::dquat { 1.0, 0.0, 0.0, 0.0 },
        });
    m_registry.emplace<IdentityComponent>(
        entity,
        IdentityComponent {
            .name = new_name,
        });

    if (parent.is_valid()) { set_parent(entity, parent); }

    return entity;
}

void Database::delete_element(db::Entity to_delete) {
    if (!valid(to_delete)) return;

    if (!m_registry.all_of<ElementComponent>(to_delete)) return;

    if (auto children = children_of(to_delete); !children.empty()) {
        QVector<Entity> copy;
        copy.reserve(children.size());
        for (auto child : children) {
            copy.push_back(child);
        }
        for (auto child : std::as_const(copy)) {
            unset_parent(child);
        }
    }

    unset_parent(to_delete);
    remove_material(to_delete);
    remove_geometry(to_delete);

    auto tags = tags_for(to_delete);
    QVector<Entity> tag_copy;
    tag_copy.reserve(tags.size());
    for (auto tag : tags) {
        tag_copy.push_back(tag);
    }
    for (auto tag : std::as_const(tag_copy)) {
        unassign_tag(to_delete, tag);
    }

    m_registry.destroy(to_delete);
}

bool Database::is_virtual_element(db::Entity element) const {
    return valid(element) && m_registry.all_of<VirtualTagComponent>(element);
}

void Database::set_virtual_element(db::Entity element, bool is_virtual) {
    if (!valid(element)) return;
    if (!m_registry.all_of<ElementComponent>(element)) return;
    if (is_virtual_element(element) == is_virtual) return;

    if (is_virtual) {
        virtual_tag.set(element, VirtualTagComponent {});
    } else {
        virtual_tag.remove(element);
    }
}

QString Database::sanitize_material_name(QString name) {
    return sanitize_new_name<MaterialComponent>(m_registry, name, "Material");
}

db::Entity Database::add_material_group(QString             new_name,
                                        QVector<db::Entity> members,
                                        db::Entity          clone_from) {
    auto set = std::unordered_set(members.begin(), members.end());

    MaterialComponent params;

    if (m_registry.valid(clone_from) and
        m_registry.all_of<MaterialComponent>(clone_from)) {

        auto& other_p = m_registry.get<MaterialComponent>(clone_from);

        params = other_p;
    } else {
        params.optics.set_ideal_absorption(SD::OpticalSide::Back);
        params.optics.set_ideal_reflection(SD::OpticalSide::Front);
    }

    for (auto mem : std::as_const(members)) {
        this->remove_material(mem);
    }

    auto ent = m_registry.create();
    m_registry.emplace<MaterialGroupComponent>(
        ent,
        MaterialGroupComponent {
            .members = QVector<Entity>(set.begin(), set.end()),
        });
    m_registry.emplace<MaterialComponent>(ent, params);

    m_registry.emplace<IdentityComponent>(
        ent, IdentityComponent { .name = new_name });

    for (auto child : set) {
        m_registry.emplace_or_replace<MaterialGroupMemberComponent>(
            child, MaterialGroupMemberComponent { .group = ent });
    }

    return ent;
}

size_t Database::material_use_count(db::Entity material) {
    auto* ptr = m_registry.try_get<MaterialGroupComponent>(material);

    if (!ptr) { return 0; }

    return ptr->members.size();
}

void Database::delete_material_group(db::Entity to_delete, db::Entity move_to) {
    qDebug() << Q_FUNC_INFO << to_delete << move_to;
    if (to_delete == move_to) {
        qWarning()
            << "Trying to delete a group and move members to the same group!";
        return;
    }

    // if not a group, bail
    if (!m_registry.all_of<MaterialGroupComponent>(to_delete)) {
        qDebug() << Q_FUNC_INFO << "Not a group";
        return;
    }

    // steal current member list
    auto members =
        std::move(m_registry.get<MaterialGroupComponent>(to_delete).members);

    // destroy current group entity
    m_registry.destroy(to_delete);

    if (m_registry.valid(move_to) and
        m_registry.all_of<MaterialGroupComponent>(move_to)) {
        // moving to valid target

        // reset member list membership
        for (auto child : members) {
            m_registry.emplace_or_replace<MaterialGroupMemberComponent>(
                child, MaterialGroupMemberComponent { .group = move_to });
        }

        m_registry.patch<MaterialGroupComponent>(
            move_to, [&](MaterialGroupComponent& a) {
                a.members.append(members.begin(), members.end());
            });
    } else {
        // invalid target. clear

        qDebug() << Q_FUNC_INFO << "Removing component";

        m_registry.remove<MaterialGroupMemberComponent>(members.begin(),
                                                        members.end());
    }
}

db::Entity Database::material_of(db::Entity element) const {
    if (auto* m = m_registry.try_get<MaterialGroupMemberComponent>(element)) {
        return m->group;
    }
    return {};
}

QString Database::sanitize_geometry_name(QString name) {
    return sanitize_new_name<GeometryComponent>(m_registry, name, "Geometry");
}

db::Entity Database::add_geometry_group(QString             new_name,
                                        QVector<db::Entity> members,
                                        db::Entity          clone_from) {
    auto set = std::unordered_set(members.begin(), members.end());

    GeometryComponent params;

    if (m_registry.valid(clone_from) and
        m_registry.all_of<GeometryComponent>(clone_from)) {

        auto& other_p = m_registry.get<GeometryComponent>(clone_from);

        // Library classes do not consistently provide clone(), so copy through
        // their serialization API.
        nlohmann::ordered_json node;

        other_p.surface->write_json(node);

        params.aperture = other_p.aperture->make_copy();
        params.surface  = SD::make_surface_from_json(node);
    } else {
        params.aperture = SD::make_aperture<SD::Circle>(1.0);
        params.surface  = SolTrace::Data::make_surface_from_type(
            SolTrace::Data::SurfaceType::FLAT, { 1.0, 1.0 });
    }

    for (auto mem : members) {
        this->remove_geometry(mem);
    }

    auto ent = m_registry.create();
    m_registry.emplace<GeometryGroupComponent>(
        ent,
        GeometryGroupComponent {
            .members = QVector<Entity>(set.begin(), set.end()),
        });
    m_registry.emplace<GeometryComponent>(ent, params);

    m_registry.emplace<IdentityComponent>(
        ent, IdentityComponent { .name = new_name });

    for (auto child : set) {
        m_registry.emplace_or_replace<GeometryGroupMemberComponent>(
            child, GeometryGroupMemberComponent { .group = ent });
    }

    return ent;
}

size_t Database::geometry_use_count(db::Entity geometry) {
    auto* ptr = m_registry.try_get<GeometryGroupComponent>(geometry);

    if (!ptr) { return 0; }

    return ptr->members.size();
}

void Database::delete_geometry_group(db::Entity to_delete, db::Entity move_to) {
    if (to_delete == move_to) {
        qWarning()
            << "Trying to delete a group and move members to the same group!";
        return;
    }

    // if not a group, bail
    if (!m_registry.all_of<GeometryGroupComponent>(to_delete)) { return; }

    // steal current member list
    auto members =
        std::move(m_registry.get<GeometryGroupComponent>(to_delete).members);

    // destroy current group entity
    m_registry.destroy(to_delete);

    if (m_registry.valid(move_to) and
        m_registry.all_of<GeometryGroupComponent>(move_to)) {
        // moving to valid target

        // reset member list membership
        for (auto child : members) {
            m_registry.emplace_or_replace<GeometryGroupMemberComponent>(
                child, GeometryGroupMemberComponent { .group = move_to });
        }

        m_registry.patch<GeometryGroupComponent>(
            move_to, [&](GeometryGroupComponent& a) {
                a.members.append(members.begin(), members.end());
            });
    } else {
        // invalid target. clear

        m_registry.remove<GeometryGroupMemberComponent>(members.begin(),
                                                        members.end());
    }
}

db::Entity Database::geometry_of(db::Entity element) const {
    if (auto* m = m_registry.try_get<GeometryGroupMemberComponent>(element)) {
        return m->group;
    }
    return {};
}

void Database::select(db::Entity to_select) {
    clear_selection();
    selected.set(to_select, SelectedComponent { });
}

void Database::add_to_selection(db::Entity to_select) {
    selected.set(to_select, SelectedComponent { });
}

void Database::deselect(db::Entity to_deselect) {
    selected.remove(to_deselect);
}

void Database::toggle_selection(db::Entity to_toggle_selection) {
    if (m_registry.all_of<SelectedComponent>(to_toggle_selection)) {
        selected.remove(to_toggle_selection);
    } else {
        selected.set(to_toggle_selection, SelectedComponent { });
    }
}

void Database::clear_selection() {
    auto                      view = selected.view();
    std::vector<entt::entity> to_deselect(view.begin(), view.end());
    for (auto e : to_deselect)
        selected.remove(e);
}

bool Database::is_selected(db::Entity e) const {
    return m_registry.valid(e) && m_registry.all_of<SelectedComponent>(e);
}

void Database::select_all_with_material(db::Entity material) {
    if (material.value == entt::null) return;
    auto view = m_registry.view<GeometryGroupMemberComponent>();
    for (auto e : view) {
        if (material_of(e) == db::Entity(material.value)) {
            add_to_selection(e);
        }
    }
}

void Database::select_all_with_geometry(db::Entity geometry) {
    if (geometry.value == entt::null) return;
    auto view = m_registry.view<GeometryGroupMemberComponent>();
    for (auto e : view) {
        if (geometry_of(e) == db::Entity(geometry.value)) {
            add_to_selection(e);
        }
    }
}

void Database::deselect_all_with_material(db::Entity material) {
    if (material.value == entt::null) return;
    auto view = m_registry.view<GeometryGroupMemberComponent>();
    for (auto e : view) {
        if (material_of(e) == db::Entity(material.value)) { deselect(e); }
    }
}

void Database::deselect_all_with_geometry(db::Entity geometry) {
    if (geometry.value == entt::null) return;
    auto view = m_registry.view<GeometryGroupMemberComponent>();
    for (auto e : view) {
        if (geometry_of(e) == db::Entity(geometry.value)) { deselect(e); }
    }
}

void Database::set_color(db::Entity to_color, QColor new_color) {
    this->color.set(to_color, ColorComponent { .color = new_color });
}

} // namespace db

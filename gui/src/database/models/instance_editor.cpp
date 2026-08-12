#include "database/models/instance_editor.h"

#include "database/components.h"
#include "database/conversion.h"
#include "utilities/euler_angles.h"
#include "utilities/math_utility.h"

#include <glm/geometric.hpp>
#include <glm/trigonometric.hpp>

#include <unordered_set>

namespace db {

namespace {

glm::dvec3 radians(QVector3D degrees) {
    return { glm::radians(static_cast<double>(degrees.x())),
             glm::radians(static_cast<double>(degrees.y())),
             glm::radians(static_cast<double>(degrees.z())) };
}

QVector3D degrees(glm::dvec3 radians) {
    return { static_cast<float>(glm::degrees(radians.x)),
             static_cast<float>(glm::degrees(radians.y)),
             static_cast<float>(glm::degrees(radians.z)) };
}

} // namespace

#define FIND(MEM)                                                              \
    if (!m_host) return;                                                       \
    if (!m_host->valid(m_entity)) return;                                      \
    auto& component = m_host->MEM;

void AnInstanceEditor::recompute() {
    if (m_host && m_host->valid(m_entity)) {
        auto previous = m_euler_angles_xyz_valid
                ? radians(m_euler_angles_xyz)
                : glm::dvec3 { 0.0 };
        auto euler = compatible_euler_xyz_from_quat(convert(orientation()),
                                                    previous);
        m_euler_angles_xyz       = degrees(euler);
        m_euler_angles_xyz_valid = true;
    } else {
        m_euler_angles_xyz_valid = false;
    }

    emit position_changed();
    emit global_position_changed();
    emit orientation_changed();
    emit euler_angles_xyz_changed();
    emit color_changed();
    emit hidden_changed();
    emit disabled_changed();
    emit virtual_element_changed();
    emit material_group_changed();
    emit current_material_changed();
    emit current_material_name_changed();
    emit geometry_group_changed();
    emit current_geometry_changed();
    emit current_geometry_name_changed();
    emit parent_changed();
    emit parent_name_changed();
    emit tags_changed();
}

void AnInstanceEditor::an_entity_changed(db::Entity e) {
    if (m_entity == e) {
        recompute();
        return;
    }

    if (material_group() == e) { emit current_material_name_changed(); }
    if (geometry_group() == e) { emit current_geometry_name_changed(); }
    if (parent() == e) { emit parent_name_changed(); }
}

AnInstanceEditor::AnInstanceEditor(QObject* parent) : QObject(parent) { }

void AnInstanceEditor::set(db::Entity ent) {
    m_euler_angles_xyz_valid = false;
    set_entity(ent);
    recompute();
}

void AnInstanceEditor::reset(Database* database) {
    if (m_host) {
        QObject::disconnect(m_host->identity.self(), nullptr, this, nullptr);
        QObject::disconnect(m_host->transform.self(), nullptr, this, nullptr);
        QObject::disconnect(
            m_host->global_transform.self(), nullptr, this, nullptr);
        QObject::disconnect(m_host->invisible.self(), nullptr, this, nullptr);
        QObject::disconnect(m_host->disabled.self(), nullptr, this, nullptr);
        QObject::disconnect(m_host->virtual_tag.self(), nullptr, this, nullptr);
        QObject::disconnect(m_host->color.self(), nullptr, this, nullptr);
        QObject::disconnect(
            m_host->material_group_membership.self(), nullptr, this, nullptr);
        QObject::disconnect(
            m_host->geometry_group_membership.self(), nullptr, this, nullptr);
        QObject::disconnect(
            m_host->tag_membership.self(), nullptr, this, nullptr);
        QObject::disconnect(m_host->parent.self(), nullptr, this, nullptr);
    }

    m_host = database;
    m_euler_angles_xyz_valid = false;

    if (!database) {
        recompute();
        return;
    }

    connect(database->identity.self(),
            &ComponentAPIBase::changed,
            this,
            &AnInstanceEditor::an_entity_changed);

    connect(database->transform.self(),
            &ComponentAPIBase::changed,
            this,
            &AnInstanceEditor::an_entity_changed);

    connect(database->global_transform.self(),
            &ComponentAPIBase::changed,
            this,
            [this](db::Entity e) {
                if (m_entity == e) { emit global_position_changed(); }
            });

    connect(database->invisible.self(),
            &ComponentAPIBase::changed,
            this,
            &AnInstanceEditor::an_entity_changed);

    connect(database->invisible.self(),
            &ComponentAPIBase::removed,
            this,
            &AnInstanceEditor::an_entity_changed);

    connect(database->disabled.self(),
            &ComponentAPIBase::changed,
            this,
            &AnInstanceEditor::an_entity_changed);

    connect(database->disabled.self(),
            &ComponentAPIBase::removed,
            this,
            &AnInstanceEditor::an_entity_changed);

    connect(database->virtual_tag.self(),
            &ComponentAPIBase::changed,
            this,
            &AnInstanceEditor::an_entity_changed);

    connect(database->virtual_tag.self(),
            &ComponentAPIBase::removed,
            this,
            &AnInstanceEditor::an_entity_changed);

    connect(database->color.self(),
            &ComponentAPIBase::changed,
            this,
            &AnInstanceEditor::an_entity_changed);

    connect(database->color.self(),
            &ComponentAPIBase::removed,
            this,
            &AnInstanceEditor::an_entity_changed);

    connect(database->material_group_membership.self(),
            &ComponentAPIBase::changed,
            this,
            &AnInstanceEditor::an_entity_changed);

    connect(database->material_group_membership.self(),
            &ComponentAPIBase::removed,
            this,
            &AnInstanceEditor::an_entity_changed);

    connect(database->geometry_group_membership.self(),
            &ComponentAPIBase::changed,
            this,
            &AnInstanceEditor::an_entity_changed);

    connect(database->geometry_group_membership.self(),
            &ComponentAPIBase::removed,
            this,
            &AnInstanceEditor::an_entity_changed);

    connect(database->tag_membership.self(),
            &ComponentAPIBase::changed,
            this,
            &AnInstanceEditor::an_entity_changed);

    connect(database->parent.self(),
            &ComponentAPIBase::changed,
            this,
            &AnInstanceEditor::an_entity_changed);

    recompute();
}

QVector3D AnInstanceEditor::position() const {

    if (m_host and m_host->valid(m_entity)) {
        if (auto tf = m_host->transform.get(m_entity); tf) {
            return convert(tf->position);
        }
    }

    return { };
}

void AnInstanceEditor::set_position(const QVector3D& newPosition) {
    if (position() == newPosition) return;

    FIND(transform);

    component.patch(m_entity, [&](TransformComponent& a) {
        a.position = convert(newPosition);
    });

    emit position_changed();
    emit global_position_changed();
}

QVector3D AnInstanceEditor::global_position() const {
    if (m_host and m_host->valid(m_entity)) {
        if (auto tf = m_host->global_transform.get(m_entity); tf) {
            return convert(tf->position);
        }
    }

    return { };
}

void AnInstanceEditor::set_global_position(const QVector3D& newPosition) {
    if (global_position() == newPosition) return;
    if (!m_host) return;
    if (!m_host->valid(m_entity)) return;

    auto local_position = convert(newPosition);
    if (auto parent_component = m_host->parent.get(m_entity);
        parent_component && m_host->valid(parent_component->parent)) {
        if (auto parent_global =
                m_host->global_transform.get(parent_component->parent);
            parent_global) {
            local_position = glm::inverse(parent_global->rotation) *
                             (local_position - parent_global->position);
        }
    }

    set_position(convert(local_position));
}

QQuaternion AnInstanceEditor::orientation() const {
    if (m_host and m_host->valid(m_entity)) {
        if (auto tf = m_host->transform.get(m_entity); tf) {
            return convert(tf->rotation);
        }
    }
    return { };
}

void AnInstanceEditor::set_orientation(const QQuaternion& newOrientation) {
    if (orientation() == newOrientation) return;

    FIND(transform);

    auto previous = m_euler_angles_xyz_valid
            ? radians(m_euler_angles_xyz)
            : glm::dvec3 { 0.0 };
    auto euler = compatible_euler_xyz_from_quat(convert(newOrientation),
                                                previous);

    component.patch(m_entity, [&](TransformComponent& a) {
        a.rotation = convert(newOrientation);
    });

    m_euler_angles_xyz       = degrees(euler);
    m_euler_angles_xyz_valid = true;

    emit orientation_changed();
    emit euler_angles_xyz_changed();
}

QVector3D AnInstanceEditor::euler_angles_xyz() const {
    if (!m_euler_angles_xyz_valid) {
        auto previous = glm::dvec3 { 0.0 };
        auto euler = compatible_euler_xyz_from_quat(convert(orientation()),
                                                    previous);
        m_euler_angles_xyz       = degrees(euler);
        m_euler_angles_xyz_valid = true;
    }

    return m_euler_angles_xyz;
}

void AnInstanceEditor::set_euler_angles_xyz(const QVector3D& angles) {
    bool changed = !m_euler_angles_xyz_valid || m_euler_angles_xyz != angles;

    m_euler_angles_xyz       = angles;
    m_euler_angles_xyz_valid = true;

    auto new_orientation = convert(euler_xyz_to_quat(radians(angles)));
    if (orientation() != new_orientation) {
        FIND(transform);

        component.patch(m_entity, [&](TransformComponent& a) {
            a.rotation = convert(new_orientation);
        });

        emit orientation_changed();
        changed = true;
    }

    if (changed) { emit euler_angles_xyz_changed(); }
}

QColor AnInstanceEditor::color() const {
    if (m_host and m_host->valid(m_entity)) {
        if (auto tf = m_host->color.get(m_entity); tf) { return tf->color; }
    }

    return Qt::white;
}

void AnInstanceEditor::set_color(QColor newColor) {
    if (color() == newColor) return;
    if (!m_host) return;
    if (!m_host->valid(m_entity)) return;

    m_host->set_color(m_entity, newColor);

    emit color_changed();
}

bool AnInstanceEditor::hidden() const {
    if (m_host and m_host->valid(m_entity)) {
        if (m_host->as_registry().any_of<InvisibleComponent>(m_entity)) {
            return true;
        }
    }

    return false;
}

void AnInstanceEditor::set_hidden(bool newHidden) {
    if (hidden() == newHidden) return;
    FIND(invisible);

    if (newHidden) {
        component.set(m_entity, InvisibleComponent { });
    } else {
        component.remove(m_entity);
    }

    emit hidden_changed();
}

bool AnInstanceEditor::disabled() const {
    if (m_host and m_host->valid(m_entity)) {
        return m_host->as_registry().any_of<DisabledComponent>(m_entity);
    }

    return false;
}

void AnInstanceEditor::set_disabled(bool newDisabled) {
    if (disabled() == newDisabled) return;
    FIND(disabled);

    if (newDisabled) {
        component.set(m_entity, DisabledComponent { });
    } else {
        component.remove(m_entity);
    }

    emit disabled_changed();
}

bool AnInstanceEditor::virtual_element() const {
    return m_host && m_host->is_virtual_element(m_entity);
}

void AnInstanceEditor::set_virtual_element(bool newVirtualElement) {
    if (virtual_element() == newVirtualElement) return;
    if (!m_host) return;
    if (!m_host->valid(m_entity)) return;

    m_host->set_virtual_element(m_entity, newVirtualElement);

    emit virtual_element_changed();
}

db::Entity AnInstanceEditor::material_group() const {
    if (m_host and m_host->valid(m_entity)) {
        if (auto tf = m_host->material_group_membership.get(m_entity); tf) {
            return tf->group;
        }
    }

    return { };
}

void AnInstanceEditor::set_material_group(db::Entity newGroup) {
    if (material_group() == newGroup) return;
    if (!m_host) return;
    if (!m_host->valid(m_entity)) return;

    m_host->assign_material(m_entity, newGroup);

    emit material_group_changed();
    emit current_material_changed();
    emit current_material_name_changed();
}

db::Entity AnInstanceEditor::geometry_group() const {
    if (m_host and m_host->valid(m_entity)) {
        if (auto tf = m_host->geometry_group_membership.get(m_entity); tf) {
            return tf->group;
        }
    }

    return { };
}

void AnInstanceEditor::set_geometry_group(db::Entity newGroup) {
    if (geometry_group() == newGroup) return;
    if (!m_host) return;
    if (!m_host->valid(m_entity)) return;

    m_host->assign_geometry(m_entity, newGroup);

    emit geometry_group_changed();
    emit current_geometry_changed();
    emit current_geometry_name_changed();
}

Entity AnInstanceEditor::current_material() const {
    return material_group();
}

void AnInstanceEditor::set_current_material(Entity newGroup) {
    set_material_group(newGroup);
}

QString AnInstanceEditor::current_material_name() const {
    if (m_host and m_host->valid(current_material())) {
        return m_host->name_of(current_material());
    }

    return { };
}

Entity AnInstanceEditor::current_geometry() const {
    return geometry_group();
}

void AnInstanceEditor::set_current_geometry(Entity newGroup) {
    set_geometry_group(newGroup);
}

QString AnInstanceEditor::current_geometry_name() const {
    if (m_host and m_host->valid(current_geometry())) {
        return m_host->name_of(current_geometry());
    }

    return { };
}

db::Entity AnInstanceEditor::parent() const {
    if (m_host and m_host->valid(m_entity)) {
        if (auto tf = m_host->parent.get(m_entity); tf) { return tf->parent; }
    }

    return { };
}

void AnInstanceEditor::set_parent(db::Entity newParent) {
    if (parent() == newParent) return;

    if (newParent == m_entity) {
        emit notify(
            ANotification::error("An element cannot be its own parent."));
        return;
    }

    // Guard cycles
    {
        auto cursor = newParent;

        while (true) {
            auto ptr = m_host->parent.get(cursor);

            if (!ptr) { break; }

            if (ptr->parent == m_entity) {
                // cycle
                emit notify(
                    ANotification::error("That parent would create a loop in "
                                         "the element hierarchy."));
                return;
            }

            cursor = ptr->parent;
        }
    }

    m_host->set_parent(m_entity, newParent);

    emit parent_changed();
    emit parent_name_changed();
    emit global_position_changed();
}

QString AnInstanceEditor::parent_name() const {
    if (m_host and m_host->valid(parent())) {
        return m_host->name_of(parent());
    }

    return { };
}

QVector<db::Entity> AnInstanceEditor::tags() const {
    if (m_host and m_host->valid(m_entity)) {
        if (auto tf = m_host->tag_membership.get(m_entity); tf) {
            return tf->tags;
        }
    }

    return { };
}

void AnInstanceEditor::set_tags(QVector<db::Entity> const& newTags) {
    auto current_tags = tags();

    std::unordered_set<db::Entity> incoming(newTags.begin(), newTags.end());
    std::unordered_set<db::Entity> current(current_tags.begin(),
                                           current_tags.end());

    if (incoming == current) return;
    if (!m_host) return;
    if (!m_host->valid(m_entity)) return;

    for (auto new_tag : incoming) {
        if (!current.contains(new_tag)) {
            m_host->assign_tag(m_entity, new_tag);
        }
    }

    // Remove tags no longer present.
    for (auto old_tag : current) {
        if (!incoming.contains(old_tag)) {
            m_host->unassign_tag(m_entity, old_tag);
        }
    }

    emit tags_changed();
}

QString AnInstanceEditor::entity_name() const {
    if (m_host and m_host->valid(m_entity)) {
        if (auto tf = m_host->identity.get(m_entity); tf) { return tf->name; }
    }

    return { };
}

void AnInstanceEditor::set_entity_name(const QString& newEntity_name) {
    if (entity_name() == newEntity_name) return;

    FIND(identity);

    component.set(m_entity, IdentityComponent { .name = newEntity_name });

    emit entity_name_changed();
}

void AnInstanceEditor::set_from_angles(QVector3D angles) {
    set_euler_angles_xyz(angles);
}

void AnInstanceEditor::look_at_world_position(QVector3D targetPosition) {
    if (!m_host) return;
    if (!m_host->valid(m_entity)) return;

    auto* global = m_host->global_transform.get(m_entity);
    if (!global) return;

    auto target    = convert(targetPosition);
    auto direction = target - global->position;

    if (glm::length(direction) < 1e-8) {
        emit notify(ANotification::error(
            "Choose a different target position for this element."));
        return;
    }

    auto target_global_rotation = dir_roll_to_quat(direction, 0.0);

    glm::dquat parent_global_rotation { 1.0, 0.0, 0.0, 0.0 };
    if (auto* parent = m_host->parent.get(m_entity); parent) {
        if (auto* parent_global = m_host->global_transform.get(parent->parent);
            parent_global) {
            parent_global_rotation = parent_global->rotation;
        }
    }

    auto local_rotation =
        glm::inverse(parent_global_rotation) * target_global_rotation;

    set_orientation(convert(glm::normalize(local_rotation)));
}

void AnInstanceEditor::look_at_entity(Entity target) {
    if (!m_host) return;
    if (!target.is_valid() || !m_host->valid(target)) return;

    auto* global = m_host->global_transform.get(target);
    if (!global) return;

    look_at_world_position(convert(global->position));
}

void AnInstanceEditor::clear_parent() {
    if (m_host and m_host->valid(m_entity)) {
        m_host->unset_parent(m_entity);
        emit parent_changed();
        emit parent_name_changed();
        emit global_position_changed();
    }
}

} // namespace db

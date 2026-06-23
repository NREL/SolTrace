#include "database_models.h"

#include "components.h"
#include "conversion.h"
#include "database/database.h"
#include "database/database_notification.h"
#include "utilities/math_utility.h"

namespace db {

void NameModel::recompute(db::Entity e) {
    if (!m_host or !e.is_valid() or !m_host->valid(e)) return;

    set_name(m_host->name_of(e));
}

NameModel::NameModel(QObject* parent) : QObject(parent) { }

void NameModel::reset(Database* database) {
    if (m_host) {
        QObject::disconnect(m_host->identity.self(), nullptr, this, nullptr);
    }

    m_host = database;

    if (!m_host) return;

    connect(m_host->identity.self(),
            &ComponentAPIBase::changed,
            this,
            &NameModel::recompute);
}

void NameModel::update_name(QString name) {
    if (!m_host or !m_node.is_valid() or !m_host->valid(m_node)) return;

    if (m_host->name_of(m_node) != m_name) {
        m_host->identity.patch(
            m_node, [&](IdentityComponent& ident) { ident.name = m_name; });
    }
}

// =============================================================================

void BreadcrumbModel::recompute() {
    m_path.clear();

    if (!m_host) return;

    if (!m_host->valid(m_node)) return;

    std::unordered_set<entt::entity> seen;

    entt::entity looking_at = m_node;

    while (true) {

        if (!m_host->valid(looking_at)) { break; }

        if (seen.contains(looking_at)) {
            // cycles??
            break;
        }

        seen.insert(looking_at);

        m_path.push_back(looking_at);

        looking_at = m_host->parent_of(looking_at);
    }

    std::reverse(m_path.begin(), m_path.end());

    QVector<EntityNamePair> ret;

    for (auto e : std::as_const(m_path)) {
        ret << EntityNamePair { m_host->name_of(e), e };
    }

    this->store_reset(ret);
}

BreadcrumbModel::BreadcrumbModel(QObject* parent) : StructModelAdapter(parent) {
    connect(this,
            &BreadcrumbModel::node_changed,
            this,
            &BreadcrumbModel::recompute);
}

void BreadcrumbModel::reset(Database* database) {
    if (m_host) {
        QObject::disconnect(m_host->identity.self(), nullptr, this, nullptr);
        QObject::disconnect(m_host->parent.self(), nullptr, this, nullptr);
    }

    m_host = database;
    m_node = {};
    recompute();

    if (database) {
        connect(database->identity.self(),
                &ComponentAPIBase::changed,
                this,
                [this](entt::entity e) {
                    if (m_path.contains(Entity(e))) { this->recompute(); }
                });

        connect(database->parent.self(),
                &ComponentAPIBase::changed,
                this,
                [this](entt::entity e) {
                    if (m_path.contains(Entity(e))) { this->recompute(); }
                });
    }
}

// =============================================================================

EntityNamePair EntityNamePair::record_for_entity(Database&  db,
                                                 db::Entity entity) {
    return EntityNamePair { .name = db.name_of(entity), .entity = entity };
}

QVector<EntityNamePair> ChildModel::rebuild_lists() {
    QVector<EntityNamePair> new_recs;
    m_reverse.clear();

    if (!m_host) return {};

    if (!m_host->valid(m_node)) return {};

    auto children = m_host->children_of(m_node);

    // Copy here, can revise later
    new_recs.reserve(children.size());

    for (auto x : children) {
        new_recs << EntityNamePair::record_for_entity(*m_host, x);
    }

    for (size_t i = 0; i < new_recs.size(); i++) {
        m_reverse[new_recs[i].entity] = i;
    }

    return new_recs;
}

void ChildModel::recompute() {
    auto r = rebuild_lists();

    this->store_reset(r);
}

void ChildModel::ident_changed(db::Entity e) {
    if (!m_host) return;

    if (auto iter = m_reverse.find(e); iter != m_reverse.end()) {
        this->store_push_update(iter->second,
                                EntityNamePair::record_for_entity(*m_host, e));
    }
}

ChildModel::ChildModel(QObject* parent) : StructModelAdapter(parent) {
    connect(this, &ChildModel::node_changed, this, &ChildModel::recompute);
}

void ChildModel::reset(Database* database) {
    if (m_host) {
        QObject::disconnect(m_host->children.self(), nullptr, this, nullptr);
        QObject::disconnect(m_host->identity.self(), nullptr, this, nullptr);
    }

    m_host = database;
    m_node = {};
    recompute();

    if (database) {
        connect(database->children.self(),
                &ComponentAPIBase::changed,
                this,
                [this](entt::entity e) {
                    if ((entt::entity)node() == e) { recompute(); }
                });

        connect(database->children.self(),
                &ComponentAPIBase::removed,
                this,
                [this](entt::entity e) {
                    if ((entt::entity)node() == e) { recompute(); }
                });

        connect(database->identity.self(),
                &ComponentAPIBase::changed,
                this,
                &ChildModel::ident_changed);
    }
}

// =============================================================================

QVector<EntityNamePair> MaterialGroupsModel::rebuild_lists() {
    QVector<EntityNamePair> new_recs;
    m_reverse.clear();

    if (!m_host) return {};

    auto view = m_host->as_registry().view<MaterialGroupComponent>();

    for (auto const& [e, group] : view.each()) {
        new_recs.push_back(EntityNamePair::record_for_entity(*m_host, e));
    }

    for (size_t i = 0; i < new_recs.size(); i++) {
        m_reverse[new_recs[i].entity] = i;
    }

    return new_recs;
}

void MaterialGroupsModel::recompute() {
    auto r = rebuild_lists();

    this->store_reset(r);
}

void MaterialGroupsModel::group_changed(db::Entity e) {
    if (!m_host) return;

    auto iter = m_reverse.find(e);

    if (iter == m_reverse.end()) {
        if (!m_host->as_registry().all_of<MaterialGroupComponent>(e)) return;

        return recompute();
    }

    this->store_push_update(iter->second,
                            EntityNamePair::record_for_entity(*m_host, e));
}
void MaterialGroupsModel::group_removed(db::Entity e) {
    qDebug() << Q_FUNC_INFO << db::Entity(e);
    recompute();
}

MaterialGroupsModel::MaterialGroupsModel(QObject* parent)
    : StructModelAdapter(parent) { }

void MaterialGroupsModel::reset(Database* database) {
    if (m_host) {
        QObject::disconnect(m_host->identity.self(), nullptr, this, nullptr);
        QObject::disconnect(
            m_host->material_root.self(), nullptr, this, nullptr);
    }
    qDebug() << Q_FUNC_INFO << database;
    m_host = database;
    recompute();

    if (!database) { return; }

    connect(database->identity.self(),
            &ComponentAPIBase::changed,
            this,
            &MaterialGroupsModel::group_changed);

    connect(database->material_root.self(),
            &ComponentAPIBase::changed,
            this,
            &MaterialGroupsModel::group_changed);

    connect(database->material_root.self(),
            &ComponentAPIBase::removed,
            this,
            &MaterialGroupsModel::group_removed);
}

QVariant MaterialGroupsModel::get(int index) {
    auto rec = get_at(index);
    if (!rec) return {};
    return QVariant::fromValue(rec->entity);
}

int MaterialGroupsModel::index_of(db::Entity entity) const {
    auto iter = m_reverse.find(entity);
    if (iter == m_reverse.end()) return -1;
    return iter->second;
}

// =============================================================================

QVector<EntityNamePair> GeometryGroupsModel::rebuild_lists() {
    QVector<EntityNamePair> new_recs;
    m_reverse.clear();

    if (!m_host) return {};

    auto view = m_host->as_registry().view<GeometryGroupComponent>();

    for (auto const& [e, group] : view.each()) {
        new_recs.push_back(EntityNamePair::record_for_entity(*m_host, e));
    }

    for (size_t i = 0; i < new_recs.size(); i++) {
        m_reverse[new_recs[i].entity] = i;
    }

    return new_recs;
}

void GeometryGroupsModel::recompute() {
    auto r = rebuild_lists();

    this->store_reset(r);
}

void GeometryGroupsModel::group_changed(db::Entity e) {
    if (!m_host) return;

    auto iter = m_reverse.find(e);

    if (iter == m_reverse.end()) {
        if (!m_host->as_registry().all_of<GeometryGroupComponent>(e)) return;

        return recompute();
    }

    this->store_push_update(iter->second,
                            EntityNamePair::record_for_entity(*m_host, e));
}
void GeometryGroupsModel::group_removed(db::Entity e) {
    recompute();
}

GeometryGroupsModel::GeometryGroupsModel(QObject* parent)
    : StructModelAdapter(parent) { }

int GeometryGroupsModel::index_of(db::Entity entity) const {
    auto iter = m_reverse.find(entity);
    if (iter == m_reverse.end()) return -1;
    return iter->second;
}

void GeometryGroupsModel::reset(Database* database) {
    if (m_host) {
        QObject::disconnect(m_host->identity.self(), nullptr, this, nullptr);
        QObject::disconnect(
            m_host->geometry_root.self(), nullptr, this, nullptr);
    }
    qDebug() << Q_FUNC_INFO << database;
    m_host = database;
    recompute();

    if (!database) { return; }

    connect(database->identity.self(),
            &ComponentAPIBase::changed,
            this,
            &GeometryGroupsModel::group_changed);

    connect(database->geometry_root.self(),
            &ComponentAPIBase::changed,
            this,
            &GeometryGroupsModel::group_changed);

    connect(database->geometry_root.self(),
            &ComponentAPIBase::removed,
            this,
            &GeometryGroupsModel::group_removed);
}

QVariant GeometryGroupsModel::get(int index) {
    auto rec = get_at(index);
    if (!rec) return {};
    return QVariant::fromValue(rec->entity);
}

// =============================================================================

QVector<EntityNamePair> TagsModel::rebuild_lists() {
    QVector<EntityNamePair> new_recs;
    m_reverse.clear();

    if (!m_host) return {};

    auto view = m_host->as_registry().view<TagComponent>();

    for (auto const& e : view.each()) {
        new_recs.push_back(
            EntityNamePair::record_for_entity(*m_host, std::get<0>(e)));
    }

    for (size_t i = 0; i < new_recs.size(); i++) {
        m_reverse[new_recs[i].entity] = i;
    }

    return new_recs;
}

void TagsModel::recompute() {
    auto r = rebuild_lists();

    this->store_reset(r);
}

void TagsModel::tag_changed(db::Entity e) {
    if (!m_host) return;

    auto iter = m_reverse.find(e);

    if (iter == m_reverse.end()) { return recompute(); }

    this->store_push_update(iter->second,
                            EntityNamePair::record_for_entity(*m_host, e));
}
void TagsModel::tag_removed(db::Entity e) {
    recompute();
}

TagsModel::TagsModel(QObject* parent) : StructModelAdapter(parent) { }

void TagsModel::reset(Database* database) {
    if (m_host) {
        QObject::disconnect(m_host->tag_root.self(), nullptr, this, nullptr);
        QObject::disconnect(m_host->identity.self(), nullptr, this, nullptr);
    }

    m_host = database;
    recompute();

    if (!database) return;

    connect(database->tag_root.self(),
            &ComponentAPIBase::changed,
            this,
            &TagsModel::tag_changed);

    connect(database->tag_root.self(),
            &ComponentAPIBase::removed,
            this,
            &TagsModel::tag_removed);

    connect(database->identity.self(),
            &ComponentAPIBase::changed,
            this,
            &TagsModel::tag_changed);
}

// =============================================================================


#define FIND(MEM)                                                              \
    if (!m_host) return;                                                       \
    if (!m_host->valid(m_entity)) return;                                      \
    auto& component = m_host->MEM;

void AnInstanceEditor::recompute() {
    emit position_changed();
    emit global_position_changed();
    emit orientation_changed();
    emit color_changed();
    emit hidden_changed();
    emit disabled_changed();
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

AnInstanceEditor::AnInstanceEditor(QObject* parent) : QObject(parent) {
    /*
     *     auto lock = m_host.lock();

     if (!lock) { return set_empty(); }

      auto note = get_notifier(*lock);

     if (!note) { return set_empty(); }

      if (!lock->valid(m_entity)) { return set_empty(); }

     */
}

void AnInstanceEditor::set(db::Entity ent) {
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

    return {};
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

    return {};
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
    return {};
}

void AnInstanceEditor::set_orientation(const QQuaternion& newOrientation) {
    if (orientation() == newOrientation) return;

    FIND(transform);

    component.patch(m_entity, [&](TransformComponent& a) {
        a.rotation = convert(newOrientation);
    });

    emit orientation_changed();
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
        component.set(m_entity, InvisibleComponent {});
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
        component.set(m_entity, DisabledComponent {});
    } else {
        component.remove(m_entity);
    }

    emit disabled_changed();
}

db::Entity AnInstanceEditor::material_group() const {
    if (m_host and m_host->valid(m_entity)) {
        if (auto tf = m_host->material_group_membership.get(m_entity); tf) {
            return tf->group;
        }
    }

    return {};
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

    return {};
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

    return {};
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

    return {};
}

db::Entity AnInstanceEditor::parent() const {
    if (m_host and m_host->valid(m_entity)) {
        if (auto tf = m_host->parent.get(m_entity); tf) { return tf->parent; }
    }

    return {};
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
                emit notify(ANotification::error(
                    "That parent would create a loop in the element hierarchy."));
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

    return {};
}

QVector<db::Entity> AnInstanceEditor::tags() const {
    if (m_host and m_host->valid(m_entity)) {
        if (auto tf = m_host->tag_membership.get(m_entity); tf) {
            return tf->tags;
        }
    }

    return {};
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

    return {};
}

void AnInstanceEditor::set_entity_name(const QString& newEntity_name) {
    if (entity_name() == newEntity_name) return;

    FIND(identity);

    component.set(m_entity, IdentityComponent { .name = newEntity_name });

    emit entity_name_changed();
}

void AnInstanceEditor::set_from_angles(QVector3D angles) {
    set_orientation(QQuaternion::fromEulerAngles(angles));
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

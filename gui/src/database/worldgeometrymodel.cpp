#include "worldgeometrymodel.h"
#include "utilities/math_utility.h"

namespace db {

void InstancedElements::on_geometry_group_change(entt::entity group) {
    if (group != entt::null && Entity(group) != m_target_group) return;

    m_instance_data.clear();
    m_member_cache.clear();
    m_rev_cache.clear();


    if (!m_database) {
        markDirty();
        return;
    }

    auto* ptr = m_database->geometry_root.get(group);

    if (!ptr) {
        markDirty();
        return;
    }

    for (auto member : std::as_const(ptr->members)) {

        auto global = m_database->global_transform.get(member);
        if (!global) continue;

        if (m_database->as_registry().all_of<HasFluxMapComponent>(member)) {
            continue;
        }

        QColor color = m_database->is_selected(member) ? Qt::yellow
                       : m_database->color.get(member)
                           ? m_database->color.get(member)->color
                           : Qt::white;

        // qDebug() << convert(global->position) << convert(global->rotation);

        auto entry =
            calculateTableEntryFromQuaternion(convert(global->position),
                                              QVector3D(1, 1, 1),
                                              convert(global->rotation),
                                              color);

        m_instance_data.append(reinterpret_cast<const char*>(&entry),
                               sizeof(entry));

        m_rev_cache[member] = m_member_cache.size();
        m_member_cache.push_back(member);
    }

    // qDebug() << Q_FUNC_INFO << "group" << entt::to_integral(group) << "->"
    //          << m_member_cache.size();

    markDirty();
}

void InstancedElements::on_geometry_group_membership_change(
    entt::entity entity) {
    auto* ptr = m_database->geometry_group_membership.get(entity);

    if (!ptr) return;

    if (!m_target_group.is_valid()) return;

    if (ptr->group != m_target_group) return;

    on_geometry_group_change(m_target_group);
}

void InstancedElements::on_instance_changed(entt::entity e) {
    auto iter = m_rev_cache.find(e);

    if (iter == m_rev_cache.end()) { return; }

    // TODO update only a single instance on change

    on_geometry_group_change(m_target_group);
}

void InstancedElements::on_selection_changed(entt::entity e) {
    auto iter = m_rev_cache.find(e);

    if (iter == m_rev_cache.end()) { return; }

    on_geometry_group_change(m_target_group);
}

entt::entity InstancedElements::entity_at(int index) {
    if (index < 0 || index >= m_member_cache.size()) { return entt::null; }
    return m_member_cache[index];
}

void InstancedElements::toggle_selection(int index) {
    entt::entity instance = entity_at(index);
    if (instance == entt::null) return;
    m_database->toggle_selection(instance);
    on_geometry_group_change(m_target_group);
}

void InstancedElements::set_color(int index, QColor color) {
    entt::entity instance = entity_at(index);
    if (instance == entt::null) return;
    m_database->set_color(instance, color);
    on_geometry_group_change(m_target_group);
}

void InstancedElements::set_all_color(QColor color) {
    auto* ptr = m_database->geometry_root.get(m_target_group);
    for (auto member : m_member_cache) {
        m_database->set_color(member, color);
    }
    on_geometry_group_change(m_target_group);
}

Entity InstancedElements::at(int index) {
    if (index < 0 || index >= m_member_cache.size()) { return { }; }
    return m_member_cache[index];
}

Entity InstancedElements::material_of_group() {
    return m_database->material_of(m_target_group);
}

Entity InstancedElements::geometry_of_group() {
    return m_database->geometry_of(m_target_group);
}

Entity InstancedElements::material_of(int index) {
    entt::entity instance = entity_at(index);
    if (instance == entt::null) return { };
    // qDebug() << m_database->material_of(instance);
    return m_database->material_of(instance);
}

Entity InstancedElements::geometry_of(int index) {
    entt::entity instance = entity_at(index);
    if (instance == entt::null) return { };
    return m_database->geometry_of(instance);
}

InstancedElements::InstancedElements(Database*       db,
                                     entt::entity    group,
                                     QQuick3DObject* parent)
    : QQuick3DInstancing(parent), m_database(db), m_target_group(group) {
    if (!db) return;

    connect(db->geometry_group_membership.self(),
            &ComponentAPIBase::changed,
            this,
            &InstancedElements::on_geometry_group_membership_change);

    connect(db->geometry_group_membership.self(),
            &ComponentAPIBase::removed,
            this,
            &InstancedElements::on_geometry_group_membership_change);

    connect(db->global_transform.self(),
            &ComponentAPIBase::changed,
            this,
            &InstancedElements::on_instance_changed);

    connect(db->geometry_root.self(),
            &ComponentAPIBase::changed,
            this,
            &InstancedElements::on_geometry_group_change);

    connect(db->geometry_root.self(),
            &ComponentAPIBase::removed,
            this,
            &InstancedElements::on_geometry_group_change);

    // Selection connections
    connect(db->selected.self(),
            &ComponentAPIBase::changed,
            this,
            &InstancedElements::on_selection_changed);

    connect(db->selected.self(),
            &ComponentAPIBase::removed,
            this,
            &InstancedElements::on_selection_changed);

    on_geometry_group_change(group);
}

QByteArray InstancedElements::getInstanceBuffer(int* instanceCount) {

    if (instanceCount) { *instanceCount = m_member_cache.size(); }

    // qDebug() << Q_FUNC_INFO << m_target_group << m_member_cache.size()
    //          << m_instance_data.size();

    return m_instance_data;
}

// =============================================================================

static VisibleGroup vis_assets_for_entity(Database& db, entt::entity e) {
    auto vg = VisibleGroup {
        .geometry_group_entity = e,
        .group_instances       = std::make_shared<InstancedElements>(&db, e),
        .group_geometry        = std::make_shared<SurfaceGeometry>(),
    };

    auto grp = db.geometry_parameters.get(e);

    if (!grp) {
        qWarning() << "Unable to get group parameters";
        return vg;
    }

    vg.group_geometry->set(&db, e);
    vg.group_geometry->set_add_thickness(true);
    vg.group_geometry->set_thickness(.05);

    return vg;
}

QVector<VisibleGroup> WorldGeometryModel::rebuild_lists() {
    qDebug() << Q_FUNC_INFO;
    QVector<VisibleGroup> new_recs;
    m_reverse.clear();

    if (!m_host) return { };

    auto view = m_host->as_registry().view<GeometryGroupComponent>();

    for (auto const& [e, group] : view.each()) {
        new_recs.push_back(vis_assets_for_entity(*m_host, e));
    }

    for (size_t i = 0; i < new_recs.size(); i++) {
        m_reverse[new_recs[i].geometry_group_entity] = i;
    }

    qDebug() << Q_FUNC_INFO << "Done" << new_recs.size();

    return new_recs;
}

void WorldGeometryModel::recompute() {
    auto r = rebuild_lists();

    this->store_reset(r);
}

void WorldGeometryModel::group_changed(entt::entity e) {
    if (!m_host) return;

    auto iter = m_reverse.find(e);

    if (iter == m_reverse.end()) { return recompute(); }
}
void WorldGeometryModel::group_removed(entt::entity e) {
    recompute();
}

void WorldGeometryModel::set_all_color(QColor color) {
    for (auto const& vg : m_records) {
        vg.group_instances->set_all_color(color);
    }
}

WorldGeometryModel::WorldGeometryModel(QObject* parent)
    : StructModelAdapter(parent) { }

void WorldGeometryModel::reset(Database* database) {
    if (m_host) {
        disconnect(m_host->geometry_root.self(), nullptr, this, nullptr);
        disconnect(m_host->identity.self(), nullptr, this, nullptr);
    }

    qDebug() << Q_FUNC_INFO << database;
    m_host = database;
    recompute();

    if (!database) { return; }

    connect(database->geometry_root.self(),
            &ComponentAPIBase::changed,
            this,
            &WorldGeometryModel::group_changed);

    connect(database->geometry_root.self(),
            &ComponentAPIBase::removed,
            this,
            &WorldGeometryModel::group_removed);

    connect(database->identity.self(),
            &ComponentAPIBase::changed,
            this,
            &WorldGeometryModel::group_changed);
}

} // namespace db

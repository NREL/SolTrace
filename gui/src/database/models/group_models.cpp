#include "database/models/group_models.h"

#include "database/components.h"

#include <QDebug>

namespace db {

QVector<EntityNamePair> MaterialGroupsModel::rebuild_lists() {
    QVector<EntityNamePair> new_recs;
    m_reverse.clear();

    if (!m_host) return { };

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
    if (!rec) return { };
    return QVariant::fromValue(rec->entity);
}

int MaterialGroupsModel::index_of(db::Entity entity) const {
    auto iter = m_reverse.find(entity);
    if (iter == m_reverse.end()) return -1;
    return iter->second;
}

QVector<EntityNamePair> GeometryGroupsModel::rebuild_lists() {
    QVector<EntityNamePair> new_recs;
    m_reverse.clear();

    if (!m_host) return { };

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
    if (!rec) return { };
    return QVariant::fromValue(rec->entity);
}

QVector<EntityNamePair> TagsModel::rebuild_lists() {
    QVector<EntityNamePair> new_recs;
    m_reverse.clear();

    if (!m_host) return { };

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

} // namespace db

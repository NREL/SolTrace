#include "database/models/element_models.h"

#include "database/components.h"
#include "database/database_notification.h"

namespace db {

QVector<EntityNamePair> RootElementsModel::rebuild_lists() {
    QVector<EntityNamePair> new_recs;
    m_reverse.clear();

    if (!m_host) return { };

    auto view = m_host->as_registry().view<ElementComponent>(
        entt::exclude<ChildOfComponent>);

    for (auto [entity] : view.each()) {
        new_recs.push_back(EntityNamePair::record_for_entity(*m_host, entity));
    }

    for (qsizetype i = 0; i < new_recs.size(); ++i) {
        m_reverse[new_recs[i].entity] = static_cast<int>(i);
    }

    return new_recs;
}

void RootElementsModel::recompute() {
    store_reset(rebuild_lists());
}

void RootElementsModel::record_changed(entt::entity entity) {
    if (!m_host) return;

    if (auto iter = m_reverse.find(entity); iter != m_reverse.end()) {
        store_push_update(iter->second,
                          EntityNamePair::record_for_entity(*m_host, entity));
    }
}

RootElementsModel::RootElementsModel(QObject* parent)
    : StructModelAdapter(parent) { }

void RootElementsModel::reset(Database* database) {
    if (m_host) {
        disconnect(m_host->identity.self(), nullptr, this, nullptr);
        disconnect(m_host->parent.self(), nullptr, this, nullptr);
        disconnect(m_host->element_tag.self(), nullptr, this, nullptr);
        disconnect(m_host->children.self(), nullptr, this, nullptr);
    }

    m_host = database;
    recompute();

    if (!database) return;

    connect(database->identity.self(),
            &ComponentAPIBase::changed,
            this,
            &RootElementsModel::record_changed);

    connect(database->children.self(),
            &ComponentAPIBase::changed,
            this,
            &RootElementsModel::record_changed);

    connect(database->children.self(),
            &ComponentAPIBase::removed,
            this,
            &RootElementsModel::record_changed);

    // Root membership is determined by the presence of ChildOfComponent.
    connect(database->parent.self(),
            &ComponentAPIBase::changed,
            this,
            &RootElementsModel::recompute);

    connect(database->parent.self(),
            &ComponentAPIBase::removed,
            this,
            &RootElementsModel::recompute);

    connect(database->element_tag.self(),
            &ComponentAPIBase::changed,
            this,
            &RootElementsModel::recompute);

    connect(database->element_tag.self(),
            &ComponentAPIBase::removed,
            this,
            &RootElementsModel::recompute);
}

QVariant RootElementsModel::get(int index) {
    auto rec = get_at(index);
    if (!rec) return { };
    return QVariant::fromValue(rec->entity);
}

QVector<EntityNamePair> AllElementsModel::rebuild_lists() {
    QVector<EntityNamePair> new_recs;
    m_reverse.clear();

    if (!m_host) return { };

    auto view = m_host->as_registry().view<ElementComponent>();

    for (auto [entity] : view.each()) {
        new_recs.push_back(EntityNamePair::record_for_entity(*m_host, entity));
    }

    for (qsizetype i = 0; i < new_recs.size(); ++i) {
        m_reverse[new_recs[i].entity] = static_cast<int>(i);
    }

    return new_recs;
}

void AllElementsModel::recompute() {
    store_reset(rebuild_lists());
}

void AllElementsModel::record_changed(entt::entity entity) {
    if (!m_host) return;

    if (auto iter = m_reverse.find(entity); iter != m_reverse.end()) {
        store_push_update(iter->second,
                          EntityNamePair::record_for_entity(*m_host, entity));
    }
}

AllElementsModel::AllElementsModel(QObject* parent)
    : StructModelAdapter(parent) { }

void AllElementsModel::reset(Database* database) {
    if (m_host) {
        disconnect(m_host->identity.self(), nullptr, this, nullptr);
        disconnect(m_host->element_tag.self(), nullptr, this, nullptr);
        disconnect(m_host->children.self(), nullptr, this, nullptr);
    }

    m_host = database;
    recompute();

    if (!database) return;

    connect(database->identity.self(),
            &ComponentAPIBase::changed,
            this,
            &AllElementsModel::record_changed);

    connect(database->children.self(),
            &ComponentAPIBase::changed,
            this,
            &AllElementsModel::record_changed);

    connect(database->children.self(),
            &ComponentAPIBase::removed,
            this,
            &AllElementsModel::record_changed);

    connect(database->element_tag.self(),
            &ComponentAPIBase::changed,
            this,
            &AllElementsModel::recompute);

    connect(database->element_tag.self(),
            &ComponentAPIBase::removed,
            this,
            &AllElementsModel::recompute);
}

} // namespace db

#include "database/models/hierarchy_models.h"

#include "database/database_notification.h"

#include <algorithm>
#include <unordered_set>

namespace db {

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
    m_node = { };
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

QVector<EntityNamePair> ChildModel::rebuild_lists() {
    QVector<EntityNamePair> new_recs;
    m_reverse.clear();

    if (!m_host) return { };

    if (!m_host->valid(m_node)) return { };

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

void ChildModel::record_changed(db::Entity e) {
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
    m_node = { };
    recompute();

    if (database) {
        connect(database->children.self(),
                &ComponentAPIBase::changed,
                this,
                [this](entt::entity e) {
                    if ((entt::entity)node() == e) recompute();
                    else
                        record_changed(db::Entity(e));
                });

        connect(database->children.self(),
                &ComponentAPIBase::removed,
                this,
                [this](entt::entity e) {
                    if ((entt::entity)node() == e) recompute();
                    else
                        record_changed(db::Entity(e));
                });

        connect(database->identity.self(),
                &ComponentAPIBase::changed,
                this,
                &ChildModel::record_changed);
    }
}

} // namespace db

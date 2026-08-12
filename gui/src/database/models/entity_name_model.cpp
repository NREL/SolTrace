#include "database/models/entity_name_model.h"

#include "database/components.h"

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

EntityNamePair EntityNamePair::record_for_entity(Database&  db,
                                                 db::Entity entity) {
    return EntityNamePair {
        .name         = db.name_of(entity),
        .entity       = entity,
        .has_children = !db.children_of(entity).empty(),
    };
}

} // namespace db

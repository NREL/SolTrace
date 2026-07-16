#include "rootelementsmodel.h"

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

// =============================================================================


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

// =============================================================================

void InstanceSortFilter::recompute_has_filter() {
    bool do_name = !m_name_filter.isEmpty();
    bool do_mat  = m_material_filter.is_valid();
    bool do_geo  = m_geometry_filter.is_valid();

    set_has_filter(do_name or do_mat or do_geo);
}


void InstanceSortFilter::update_material_name(Entity entity) {
    if (!m_host) return;

    if (entity != m_material_filter) { return; }

    set_material_filter_name(m_host->name_of(entity));
}

void InstanceSortFilter::update_geometry_name(Entity entity) {
    if (!m_host) return;

    if (entity != m_geometry_filter) { return; }

    set_geometry_filter_name(m_host->name_of(entity));
}

InstanceSortFilter::InstanceSortFilter(QObject* parent)
    : QSortFilterProxyModel(parent) {

    // Connect to recompute has filter
    connect(this,
            &InstanceSortFilter::geometry_filter_changed,
            this,
            &InstanceSortFilter::recompute_has_filter);

    connect(this,
            &InstanceSortFilter::material_filter_changed,
            this,
            &InstanceSortFilter::recompute_has_filter);

    connect(this,
            &InstanceSortFilter::name_filter_changed,
            this,
            &InstanceSortFilter::recompute_has_filter);

    // Update names...

    connect(this, &InstanceSortFilter::geometry_filter_changed, this, [this]() {
        update_geometry_name(geometry_filter());
    });

    connect(this, &InstanceSortFilter::material_filter_changed, this, [this]() {
        update_material_name(material_filter());
    });

    // Invalidate filter on change...

    connect(this,
            &InstanceSortFilter::geometry_filter_changed,
            this,
            &InstanceSortFilter::invalidate);

    connect(this,
            &InstanceSortFilter::material_filter_changed,
            this,
            &InstanceSortFilter::invalidate);

    connect(this,
            &InstanceSortFilter::name_filter_changed,
            this,
            &InstanceSortFilter::invalidate);
}


void InstanceSortFilter::reset(Database* database) {
    if (m_host) {
        disconnect(m_host->identity.self(), nullptr, this, nullptr);
        disconnect(
            m_host->geometry_group_membership.self(), nullptr, this, nullptr);
        disconnect(
            m_host->material_group_membership.self(), nullptr, this, nullptr);
    }

    m_host = database;

    clear_all_filters();

    if (!database) return;

    connect(database->identity.self(),
            &ComponentAPIBase::changed,
            this,
            &InstanceSortFilter::invalidate);

    connect(database->identity.self(),
            &ComponentAPIBase::changed,
            this,
            &InstanceSortFilter::update_geometry_name);

    connect(database->identity.self(),
            &ComponentAPIBase::changed,
            this,
            &InstanceSortFilter::update_material_name);


    connect(database->geometry_group_membership.self(),
            &ComponentAPIBase::changed,
            this,
            &InstanceSortFilter::invalidate);

    connect(database->geometry_group_membership.self(),
            &ComponentAPIBase::removed,
            this,
            &InstanceSortFilter::invalidate);

    connect(database->material_group_membership.self(),
            &ComponentAPIBase::changed,
            this,
            &InstanceSortFilter::invalidate);

    connect(database->material_group_membership.self(),
            &ComponentAPIBase::removed,
            this,
            &InstanceSortFilter::invalidate);
}

void InstanceSortFilter::clear_material() {
    set_material_filter({ });
    set_material_filter_name(QString());
}

void InstanceSortFilter::clear_geometry() {
    set_geometry_filter({ });
    set_geometry_filter_name(QString());
}

void InstanceSortFilter::clear_all_filters() {
    clear_material();
    clear_geometry();
    set_name_filter({ });
}

bool db::InstanceSortFilter::filterAcceptsRow(
    int                source_row,
    QModelIndex const& source_parent) const {

    if (!sourceModel()) return true;
    if (!m_host) return true;

    // TODO: clean this up
    bool do_name = !m_name_filter.isEmpty();
    bool do_mat  = m_material_filter.is_valid();
    bool do_geo  = m_geometry_filter.is_valid();

    if (!do_name and !do_mat and !do_geo) return true;

    auto entity = sourceModel()
                      ->data(sourceModel()->index(source_row, 0, source_parent),
                             entity_role)
                      .value<Entity>();

    if (do_name) {
        auto ident = m_host->identity.get(entity);

        if (!ident) return false;

        if (!ident->name.contains(m_name_filter, Qt::CaseInsensitive)) {
            return false;
        }
    }

    if (do_mat) {
        auto ptr = m_host->material_group_membership.get(entity);

        if (!ptr) return false;

        if (Entity(ptr->group) != m_material_filter) return false;
    }

    if (do_geo) {
        auto ptr = m_host->geometry_group_membership.get(entity);

        if (!ptr) return false;

        if (Entity(ptr->group) != m_geometry_filter) return false;
    }

    return true;
}

QVariant RootElementsModel::get(int index) {
    auto rec = get_at(index);
    if (!rec) return { };
    return QVariant::fromValue(rec->entity);
}

} // namespace db

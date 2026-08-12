#include "database/models/instance_sort_filter.h"

#include "database/components.h"
#include "database/database_notification.h"

namespace db {

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

    connect(this, &InstanceSortFilter::geometry_filter_changed, this, [this]() {
        update_geometry_name(geometry_filter());
    });

    connect(this, &InstanceSortFilter::material_filter_changed, this, [this]() {
        update_material_name(material_filter());
    });

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

bool InstanceSortFilter::filterAcceptsRow(
    int                source_row,
    QModelIndex const& source_parent) const {

    if (!sourceModel()) return true;
    if (!m_host) return true;

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

} // namespace db

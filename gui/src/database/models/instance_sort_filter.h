#pragma once

#include "database/models/entity_name_model.h"

#include <QSortFilterProxyModel>

namespace db {

/// Filter proxy for element models by material, geometry, and name.
///
/// Only works with source models that expose EntityNamePair rows backed by
/// Element components.
class InstanceSortFilter : public QSortFilterProxyModel {
    Q_OBJECT

    QPointer<Database> m_host;

    Q_WRITABLE_PROPERTY(Entity, material_filter, { });
    Q_WRITABLE_PROPERTY(Entity, geometry_filter, { });
    Q_WRITABLE_PROPERTY(QString, name_filter, { });

    Q_READONLY_PROPERTY(bool, has_filter);

    Q_READONLY_PROPERTY(QString, material_filter_name);
    Q_READONLY_PROPERTY(QString, geometry_filter_name);

    static constexpr int entity_role = ROLE_FOR_MEMBER(EntityNamePair, entity);

private slots:
    void recompute_has_filter();
    void update_material_name(db::Entity);
    void update_geometry_name(db::Entity);

public:
    explicit InstanceSortFilter(QObject* parent = nullptr);

    /// Observe a database so filters can resolve entity metadata.
    void reset(Database* database);

public slots:
    /// Clear the material-group filter.
    void clear_material();

    /// Clear the geometry-group filter.
    void clear_geometry();

    /// Clear all active filters.
    void clear_all_filters();

protected:
    bool filterAcceptsRow(int                source_row,
                          QModelIndex const& source_parent) const override;
};

} // namespace db

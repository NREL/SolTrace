#pragma once

#include <QtQuick3D/qquick3dinstancing.h>
#include <QVariantMap>

#include "database/database.h"
#include "database/geometryeditor.h"

namespace db {

/// Quick3D instancing adapter for all elements that share one geometry group.
///
/// Instance data is rebuilt on the GUI thread today.
/// TODO: push instance generation to a thread
class InstancedElements : public QQuick3DInstancing {
    Q_OBJECT

    QByteArray                               m_instance_data;
    std::vector<entt::entity>                m_member_cache;
    std::unordered_map<entt::entity, size_t> m_rev_cache;

    QByteArray getInstanceBuffer(int* instanceCount) override;

    QPointer<Database> m_database;
    Entity             m_target_group;
    QColor             m_default_color = Qt::white;

private slots:
    // When the geometry information for this group has changed
    void on_geometry_group_change(entt::entity);

    // When the membership for this group has changed
    void on_geometry_group_membership_change(entt::entity);

    // When other things about an instance (parent, tf) change
    void on_instance_changed(entt::entity);

    // When an instance selection marker changes
    void on_selection_changed(entt::entity);

    // Retrieve entity using instance index
    entt::entity entity_at(int index);

public slots:
    /// Toggle selection for the instance at a rendered instance index.
    void toggle_selection(int index);

    /// Override the display color for the instance at index.
    void set_color(int index, QColor color);

    /// Set the default color for instances without a ColorComponent.
    void set_default_color(QColor color);

    /// Material group assigned to the represented geometry group.
    db::Entity material_of_group();

    /// Geometry group represented by this instancing adapter.
    db::Entity geometry_of_group();

    /// Material group assigned to the instance at index.
    db::Entity material_of(int index);

    /// Geometry group assigned to the instance at index.
    db::Entity geometry_of(int index);

    /// Entity represented by the instance at index.
    db::Entity at(int index);

public:
    explicit InstancedElements(Database*       db,
                               entt::entity    group,
                               QQuick3DObject* parent = nullptr);
};

/// One visible geometry group row consumed by the 3D scene.
struct VisibleGroup {
    entt::entity                       geometry_group_entity;
    std::shared_ptr<InstancedElements> group_instances;
    std::shared_ptr<SurfaceGeometry>   group_geometry;

    RECORD_META(VisibleGroup,
                SM_EXPOSE_RO(group_instances),
                SM_EXPOSE_RO(group_geometry), );
};

/// Model of visible geometry groups and their Quick3D adapters.
class WorldGeometryModel : public StructModelAdapter<VisibleGroup> {
    Q_OBJECT

    QPointer<Database> m_host;

    std::unordered_map<entt::entity, int> m_reverse;
    double                                m_surface_thickness = 0.05;
    unsigned                              m_subdivision_scale = 2;
    QColor                                m_default_color     = Qt::white;

    QVector<VisibleGroup> rebuild_lists();
    void                  apply_surface_options(VisibleGroup const& group);

private slots:
    void recompute();

    void group_changed(entt::entity);
    void group_removed(entt::entity);

public slots:
    /// Compute current visible world-geometry bounds on demand.
    Q_INVOKABLE QVariantMap content_bounds(bool include_flux_mapped) const;

    /// Update default color on all visible instancing adapters.
    void set_default_color(QColor color);

    /// Update generated surface thickness on all visible surface geometries.
    void set_surface_thickness(double thickness);

    /// Update generated surface subdivision on all visible surface geometries.
    void set_subdivision_scale(unsigned scale);

public:
    explicit WorldGeometryModel(QObject* parent = nullptr);
    virtual ~WorldGeometryModel() = default;

    /// Observe a new database and rebuild the visible group list.
    void reset(Database* database);
};

} // namespace db

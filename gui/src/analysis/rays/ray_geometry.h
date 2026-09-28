#pragma once


#include "jobs/job_run_common.h"
#include "support/qt_helpers.h"

#include <QQmlEngine>
#include <QVariantMap>
#include <QtQuick3D/qquick3dgeometry.h>


namespace SolTrace::GUI::Analysis {

/// A not-bitset of ray event types
struct EventTypeContainer {
    std::unordered_set<SolTrace::GUI::Data::RayEventType> events;

    EventTypeContainer() = default;

    EventTypeContainer(std::initializer_list<SolTrace::GUI::Data::RayEventType>);
    EventTypeContainer(QStringList);

    QStringList to_list() const;
};

/// A user-facing record of ray interaction.
struct UIRayRecord {
    /// Did this interaction hit an entity?
    bool has_entity = false;

    /// If it hit an entity, which one?
    qint64 entity = -1;

    /// What is the interaction type?
    QString interaction_type;

    /// Where was the interaction?
    QVector3D location;

    RECORD_META(SolTrace::GUI::Analysis::UIRayRecord,
                SM_EXPOSE_RO(has_entity),
                SM_EXPOSE_RO(entity),
                SM_EXPOSE_RO(interaction_type),
                SM_EXPOSE_RO(location), );
};

/// A table of ray interactions
struct RayInteractionTable : public Support::StructTableModel<UIRayRecord> {
    using StructTableModel::StructTableModel;
};

/// Class that builds/rebuilds ray geometry for QML visualization
///
/// TODO: make all deltas queued up for Concurrent off thread rebuilding of geom
/// TODO: move to tubes and instancing?
class RayGeometry : public QQuick3DGeometry {
    Q_OBJECT
    QML_ELEMENT

    SolTrace::GUI::Data::SimulationResultPtr m_database;

    EventTypeContainer m_include_events;

    bool      m_content_bounds_valid = false;
    QVector3D m_content_bounds_min;
    QVector3D m_content_bounds_max;

public:
    enum class TextureMode { SolidColor, Length, Segment };
    Q_ENUM(TextureMode)

    enum class IntersectionMode { Point, Line };
    Q_ENUM(IntersectionMode)

private:
    /// Which events to include in our geometry?
    Q_WRITABLE_PROPERTY(QStringList, event_include, {});

    /// How many rays to show?
    Q_WRITABLE_PROPERTY(float, show_percent, 100);

    /// How to set the UVs of the geometry?
    Q_WRITABLE_PROPERTY(TextureMode, texture_mode, TextureMode::Length);

    /// How to draw the intersections?
    Q_WRITABLE_PROPERTY(IntersectionMode, isect_mode, IntersectionMode::Line);

    /// How many rays are there?
    Q_READONLY_PROPERTY(quint64, available_rays);

    /// What is/is there a selected ray?
    Q_WRITABLE_PROPERTY(qint64, selected_ray_id, -1);

    /// Interactions for the current ray
    QOBJECT_READONLY_PROPERTY(RayInteractionTable, selected_ray_interactions);

    /// Only show rays that interacted with this entity.
    Q_WRITABLE_PROPERTY(SolTrace::GUI::Data::Entity, entity_filter, { });

    /// Display name for the selected entity filter.
    Q_READONLY_PROPERTY(QString, entity_filter_name);

    /// Rays that exit a sphere of this size are clipped
    Q_WRITABLE_PROPERTY(double, max_ray_distance, 5000);

private slots:
    void inclusion_list_update();
    void entity_filter_update();

    void selected_ray_change_update_table();

public:
    explicit RayGeometry(QQuick3DObject* parent = nullptr);

    /// Reset with new simulation results
    void set_results(SolTrace::GUI::Data::SimulationResultPtr);

public slots:
    /// Rebuild ray geometry from scratch
    void rebuild_geometry();

    /// Compute bounds for the ray geometry currently visible in the viewport.
    Q_INVOKABLE QVariantMap content_bounds() const;

    /// Given a world ray, pick a traced ray.
    /// For performance we will be dropping precision.
    void pick_ray(QVector3D world_position,
                  QVector3D world_direction,
                  float     angle_tolerance_rads);

    /// Set the entity whose interacting rays are visible.
    void select_entity_filter(SolTrace::GUI::Data::Entity entity);

    /// Clear entity filtering.
    void clear_entity_filter();
};


} // namespace SolTrace::GUI::Analysis

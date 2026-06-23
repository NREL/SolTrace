#pragma once

#include <QtQuick3D/qquick3dinstancing.h>

#include "database.h"
#include "database/database_models.h"
#include "database/geometryeditor.h"

namespace db {

// TODO: push instance generation to a thread
class InstancedElements : public QQuick3DInstancing {
    Q_OBJECT

    QByteArray                               m_instance_data;
    std::vector<entt::entity>                m_member_cache;
    std::unordered_map<entt::entity, size_t> m_rev_cache;

    QByteArray getInstanceBuffer(int* instanceCount) override;

    QPointer<Database> m_database;
    Entity             m_target_group;

private slots:
    // When the geometry information for this group has changed
    void on_geometry_group_change(entt::entity);

    // When the membership for this group has changed
    void on_geometry_group_membership_change(entt::entity);

    // When other things about an instance (parent, tf) change
    void on_instance_changed(entt::entity);

    // Retrieve entity using instance index
    entt::entity entity_at(int index);

public slots:
    // Select geometry models
    void toggle_selection(int index);

    void set_color(int index, QColor color);

    db::Entity material_of_group();
    db::Entity geometry_of_group();

    db::Entity material_of(int index);
    db::Entity geometry_of(int index);

    db::Entity at(int index);

    void set_all_color(QColor color);

public:
    explicit InstancedElements(Database*       db,
                               entt::entity    group,
                               QQuick3DObject* parent = nullptr);
};

struct VisibleGroup {
    entt::entity                       geometry_group_entity;
    std::shared_ptr<InstancedElements> group_instances;
    std::shared_ptr<SurfaceGeometry>   group_geometry;

    RECORD_META(VisibleGroup,
                SM_EXPOSE_RO(group_instances),
                SM_EXPOSE_RO(group_geometry), );
};

class WorldGeometryModel : public StructModelAdapter<VisibleGroup> {
    Q_OBJECT

    QPointer<Database> m_host;

    std::unordered_map<entt::entity, int> m_reverse;

    QVector<VisibleGroup> rebuild_lists();

private slots:
    void recompute();

    void group_changed(entt::entity);
    void group_removed(entt::entity);

public slots:
    void set_all_color(QColor color);

public:
    explicit WorldGeometryModel(QObject* parent = nullptr);
    virtual ~WorldGeometryModel() = default;

    void reset(Database* database);
};

} // namespace db

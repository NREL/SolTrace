#pragma once

#include "database/models/entity_name_model.h"
#include "utilities/notification.h"

#include <QColor>
#include <QObject>
#include <QPointer>
#include <QQuaternion>
#include <QVector3D>
#include <QVector>

namespace db {

/// A model that helps edit a material geometry instance.
class AnInstanceEditor : public QObject {
    Q_OBJECT
    QPointer<Database> m_host;

    Q_WRITABLE_PROPERTY(Entity, entity, { });

    mutable bool      m_euler_angles_xyz_valid = false;
    mutable QVector3D m_euler_angles_xyz;

    Q_PROPERTY(QString entity_name READ entity_name WRITE set_entity_name NOTIFY
                   entity_name_changed FINAL)

    Q_PROPERTY(QVector3D position READ position WRITE set_position NOTIFY
                   position_changed FINAL)
    Q_PROPERTY(QVector3D global_position READ global_position WRITE
                   set_global_position NOTIFY global_position_changed FINAL)
    Q_PROPERTY(QQuaternion orientation READ orientation WRITE set_orientation
                   NOTIFY orientation_changed FINAL)
    Q_PROPERTY(QVector3D euler_angles_xyz READ euler_angles_xyz WRITE
                   set_euler_angles_xyz NOTIFY euler_angles_xyz_changed FINAL)
    Q_PROPERTY(
        QColor color READ color WRITE set_color NOTIFY color_changed FINAL)
    Q_PROPERTY(
        bool hidden READ hidden WRITE set_hidden NOTIFY hidden_changed FINAL)
    Q_PROPERTY(bool disabled READ disabled WRITE set_disabled NOTIFY
                   disabled_changed FINAL)
    Q_PROPERTY(bool virtual_element READ virtual_element WRITE
                   set_virtual_element NOTIFY virtual_element_changed FINAL)

    Q_PROPERTY(Entity material_group READ material_group WRITE
                   set_material_group NOTIFY material_group_changed FINAL)

    Q_PROPERTY(Entity geometry_group READ geometry_group WRITE
                   set_geometry_group NOTIFY geometry_group_changed FINAL)

    Q_PROPERTY(Entity current_material READ current_material WRITE
                   set_current_material NOTIFY current_material_changed FINAL)
    Q_PROPERTY(QString current_material_name READ current_material_name NOTIFY
                   current_material_name_changed FINAL)

    Q_PROPERTY(Entity current_geometry READ current_geometry WRITE
                   set_current_geometry NOTIFY current_geometry_changed FINAL)
    Q_PROPERTY(QString current_geometry_name READ current_geometry_name NOTIFY
                   current_geometry_name_changed FINAL)

    Q_PROPERTY(
        Entity parent READ parent WRITE set_parent NOTIFY parent_changed FINAL)
    Q_PROPERTY(
        QString parent_name READ parent_name NOTIFY parent_name_changed FINAL)

    Q_PROPERTY(
        QVector<Entity> tags READ tags WRITE set_tags NOTIFY tags_changed FINAL)

private slots:
    void an_entity_changed(db::Entity);
    void recompute();

public:
    explicit AnInstanceEditor(QObject* parent = nullptr);
    virtual ~AnInstanceEditor() = default;

    /// Observe a database and synchronize the currently selected entity.
    void reset(Database* database);

    /// Select the entity to edit.
    void set(Entity);

public:
    /// Name of the selected entity.
    QString entity_name() const;

    /// Rename the selected entity.
    void set_entity_name(const QString& newEntity_name);

    /// Local transform position.
    QVector3D position() const;

    /// Set local transform position.
    void set_position(const QVector3D& newPosition);

    /// World-space transform position.
    QVector3D global_position() const;

    /// Set position by world-space coordinates.
    void set_global_position(const QVector3D& newPosition);

    /// Local transform orientation.
    QQuaternion orientation() const;

    /// Set local transform orientation.
    void set_orientation(const QQuaternion& newOrientation);

    /// Local transform orientation as compatible XYZ Euler angles, in degrees.
    QVector3D euler_angles_xyz() const;

    /// Set local transform orientation from XYZ Euler angles, in degrees.
    void set_euler_angles_xyz(const QVector3D& angles);

    /// Display color for the selected entity.
    QColor color() const;

    /// Set display color for the selected entity.
    void set_color(QColor newColor);

    /// Whether the selected entity is hidden in the viewport.
    bool hidden() const;

    /// Set whether the selected entity is hidden in the viewport.
    void set_hidden(bool newHidden);

    /// Whether the selected entity is disabled for simulation.
    bool disabled() const;

    /// Set whether the selected entity is disabled for simulation.
    void set_disabled(bool newDisabled);

    /// Whether the selected entity is a virtual/non-interacting element.
    bool virtual_element() const;

    /// Set whether the selected entity is a virtual/non-interacting element.
    void set_virtual_element(bool newVirtualElement);

    /// Material group assigned to the selected entity.
    Entity material_group() const;

    /// Assign a material group to the selected entity.
    void set_material_group(Entity newGroup);

    /// Geometry group assigned to the selected entity.
    Entity geometry_group() const;

    /// Assign a geometry group to the selected entity.
    void set_geometry_group(Entity newGroup);

    /// Currently selected material group in the editor context.
    Entity current_material() const;

    /// Select a material group in the editor context.
    void set_current_material(Entity newGroup);

    /// Display name for current_material().
    QString current_material_name() const;

    /// Currently selected geometry group in the editor context.
    Entity current_geometry() const;

    /// Select a geometry group in the editor context.
    void set_current_geometry(Entity newGroup);

    /// Display name for current_geometry().
    QString current_geometry_name() const;

    /// Parent entity of the selected entity.
    Entity parent() const;

    /// Reparent the selected entity.
    void set_parent(Entity newParent);

    /// Display name for parent().
    QString parent_name() const;

    /// Tags assigned to the selected entity.
    QVector<Entity> tags() const;

    /// Replace tags assigned to the selected entity.
    void set_tags(const QVector<Entity>& newTags);

public slots:
    /// Set orientation from Euler angles in degrees.
    void set_from_angles(QVector3D angles);

    /// Orient the selected entity toward a world-space point.
    void look_at_world_position(QVector3D targetPosition);

    /// Orient the selected entity toward another entity.
    void look_at_entity(db::Entity target);

    /// Remove the selected entity's parent.
    void clear_parent();

signals:
    void position_changed();
    void global_position_changed();
    void orientation_changed();
    void euler_angles_xyz_changed();
    void color_changed();
    void hidden_changed();
    void disabled_changed();
    void virtual_element_changed();
    void material_group_changed();
    void geometry_group_changed();
    void current_material_changed();
    void current_material_name_changed();
    void current_geometry_changed();
    void current_geometry_name_changed();
    void parent_changed();
    void parent_name_changed();
    void tags_changed();
    void entity_name_changed();

    void notify(ANotification);
};

} // namespace db

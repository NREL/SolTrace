#pragma once


#include "database/database.h"
#include "utilities/notification.h"
#include "utilities/qt_helpers.h"
#include "utilities/structmodel.h"

#include <entt/entity/fwd.hpp>

#include <QColor>
#include <QQuaternion>
#include <QStringListModel>
#include <QVector3D>

namespace db {


struct EntityNamePair {
    QString name;
    Entity  entity;
    bool    has_children = false;

    RECORD_META(db::EntityNamePair,
                SM_EXPOSE_RW(name),
                SM_EXPOSE_RO(entity),
                SM_EXPOSE_RO(has_children), );

    static EntityNamePair record_for_entity(Database& db, db::Entity entity);
};

/// Observe an entity's name
class NameModel : public QObject {
    Q_OBJECT

    QPointer<Database> m_host;

    Entity m_target;

    Q_WRITABLE_PROPERTY(Entity, node, {});
    Q_READONLY_PROPERTY(QString, name);

private slots:
    void recompute(db::Entity);

public:
    explicit NameModel(QObject* parent = nullptr);
    virtual ~NameModel() = default;

    void reset(Database* database);

public slots:
    void update_name(QString);
};

// =============================================================================

/// Get the hierarchy of an entity, walks the parent chain and provides a string
/// list, starting from the root on down.
class BreadcrumbModel : public StructModelAdapter<EntityNamePair> {
    Q_OBJECT

    QPointer<Database> m_host;

    QVector<Entity> m_path;

    Q_WRITABLE_PROPERTY(Entity, node, {});

private slots:
    void recompute();

public:
    explicit BreadcrumbModel(QObject* parent = nullptr);
    virtual ~BreadcrumbModel() = default;

    void reset(Database* database);
};

// =============================================================================

/// A model providing all children of a given entity
class ChildModel : public StructModelAdapter<EntityNamePair> {
    Q_OBJECT
    QPointer<Database> m_host;

    Q_WRITABLE_PROPERTY(Entity, node, {});

    // QVector<Entity>                 m_list;
    std::unordered_map<Entity, int> m_reverse;

    QVector<EntityNamePair> rebuild_lists();

private slots:
    void recompute();
    void record_changed(db::Entity);

public:
    explicit ChildModel(QObject* parent = nullptr);
    virtual ~ChildModel() = default;

    void reset(Database* database);
};

// =============================================================================

/// A model providing all children of a given entity
// class AllInstanceModel : public StructModelAdapter<EntityNamePair> {
//     Q_OBJECT
//     QPointer<Database> m_host;

//     Q_WRITABLE_PROPERTY(Entity, node, entt::null);

//     std::unordered_map<Entity, int> m_reverse;

//     QVector<EntityNamePair> rebuild_lists();

// private slots:
//     void recompute();
//     void ident_changed(Entity);

// public:
//     explicit AllEntityModel(QObject* parent = nullptr);
//     virtual ~AllEntityModel() = default;

//     void reset(Database* database);
// };

// =============================================================================

/// A model providing the active material groups in a database
class MaterialGroupsModel : public StructModelAdapter<EntityNamePair> {
    Q_OBJECT

    QPointer<Database> m_host;

    std::unordered_map<Entity, int> m_reverse;

    QVector<EntityNamePair> rebuild_lists();

private slots:
    void recompute();

    void group_changed(db::Entity);
    void group_removed(db::Entity);

public:
    explicit MaterialGroupsModel(QObject* parent = nullptr);
    virtual ~MaterialGroupsModel() = default;

    void reset(Database* database);

public slots:
    QVariant get(int index);
    int      index_of(db::Entity entity) const;
};

/// A model providing the active geometry groups in a database
class GeometryGroupsModel : public StructModelAdapter<EntityNamePair> {
    Q_OBJECT

    QPointer<Database> m_host;

    std::unordered_map<entt::entity, int> m_reverse;

    QVector<EntityNamePair> rebuild_lists();

private slots:
    void recompute();

    void group_changed(db::Entity);
    void group_removed(db::Entity);

public:
    explicit GeometryGroupsModel(QObject* parent = nullptr);
    virtual ~GeometryGroupsModel() = default;

    void reset(Database* database);

public slots:
    QVariant get(int index);
    int      index_of(db::Entity entity) const;
};

// =============================================================================

/// A model providing the available tags in a database
class TagsModel : public StructModelAdapter<EntityNamePair> {
    Q_OBJECT
    QPointer<Database> m_host;

    std::unordered_map<entt::entity, int> m_reverse;

    QVector<EntityNamePair> rebuild_lists();

private slots:
    void recompute();

    void tag_changed(db::Entity);
    void tag_removed(db::Entity);

public:
    explicit TagsModel(QObject* parent = nullptr);
    virtual ~TagsModel() = default;

    void reset(Database* database);
};


// =============================================================================

/// A model that helps edit a geometry instance
class AnInstanceEditor : public QObject {
    Q_OBJECT
    QPointer<Database> m_host;

    Q_WRITABLE_PROPERTY(Entity, entity, {});

    Q_PROPERTY(QString entity_name READ entity_name WRITE set_entity_name NOTIFY
                   entity_name_changed FINAL)

    Q_PROPERTY(QVector3D position READ position WRITE set_position NOTIFY
                   position_changed FINAL)
    Q_PROPERTY(QVector3D global_position READ global_position WRITE
                   set_global_position NOTIFY global_position_changed FINAL)
    Q_PROPERTY(QQuaternion orientation READ orientation WRITE set_orientation
                   NOTIFY orientation_changed FINAL)
    Q_PROPERTY(
        QColor color READ color WRITE set_color NOTIFY color_changed FINAL)
    Q_PROPERTY(
        bool hidden READ hidden WRITE set_hidden NOTIFY hidden_changed FINAL)
    Q_PROPERTY(bool disabled READ disabled WRITE set_disabled NOTIFY
                   disabled_changed FINAL)

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

    void reset(Database* database);

    void set(Entity);

public:
    QString entity_name() const;
    void    set_entity_name(const QString& newEntity_name);

    QVector3D position() const;
    void      set_position(const QVector3D& newPosition);
    QVector3D global_position() const;
    void      set_global_position(const QVector3D& newPosition);

    QQuaternion orientation() const;
    void        set_orientation(const QQuaternion& newOrientation);

    QColor color() const;
    void   set_color(QColor newColor);

    bool hidden() const;
    void set_hidden(bool newHidden);

    bool disabled() const;
    void set_disabled(bool newDisabled);

    Entity material_group() const;
    void   set_material_group(Entity newGroup);
    Entity geometry_group() const;
    void   set_geometry_group(Entity newGroup);

    Entity  current_material() const;
    void    set_current_material(Entity newGroup);
    QString current_material_name() const;

    Entity  current_geometry() const;
    void    set_current_geometry(Entity newGroup);
    QString current_geometry_name() const;

    Entity  parent() const;
    void    set_parent(Entity newParent);
    QString parent_name() const;

    QVector<Entity> tags() const;
    void            set_tags(const QVector<Entity>& newTags);

public slots:
    void set_from_angles(QVector3D angles);
    void look_at_world_position(QVector3D targetPosition);
    void look_at_entity(Entity target);
    // void set_from_dir_up(QVector3D direction, QVector3D up);

    void clear_parent();

signals:
    void position_changed();
    void global_position_changed();
    void orientation_changed();
    void color_changed();
    void hidden_changed();
    void disabled_changed();
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

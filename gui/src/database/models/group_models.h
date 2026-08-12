#pragma once

#include "database/models/entity_name_model.h"

#include <QPointer>
#include <QVariant>
#include <QVector>

#include <entt/entity/fwd.hpp>
#include <unordered_map>

namespace db {

/// A model providing the active material groups in a database.
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

    /// Observe a database and rebuild the material group list.
    void reset(Database* database);

public slots:
    /// Return the row as a QVariant map/object for QML convenience.
    QVariant get(int index);

    /// Return the row index for entity, or -1 if not present.
    int index_of(db::Entity entity) const;
};

/// A model providing the active geometry groups in a database.
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

    /// Observe a database and rebuild the geometry group list.
    void reset(Database* database);

public slots:
    /// Return the row as a QVariant map/object for QML convenience.
    QVariant get(int index);

    /// Return the row index for entity, or -1 if not present.
    int index_of(db::Entity entity) const;
};

/// A model providing the available tags in a database.
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

    /// Observe a database and rebuild the tag list.
    void reset(Database* database);
};

} // namespace db

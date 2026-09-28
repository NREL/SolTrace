#pragma once

#include "scene_models/entity_name_model.h"

#include <QPointer>
#include <QVector>

#include <unordered_map>

namespace SolTrace::GUI::Data {

/// Get the hierarchy of an entity, walking the parent chain from root down.
class BreadcrumbModel : public Support::StructModelAdapter<EntityNamePair> {
    Q_OBJECT

    QPointer<Database> m_host;

    QVector<Entity> m_path;

    Q_WRITABLE_PROPERTY(Entity, node, { });

private slots:
    void recompute();

public:
    explicit BreadcrumbModel(QObject* parent = nullptr);
    virtual ~BreadcrumbModel() = default;

    /// Observe a database and rebuild the breadcrumb path.
    void reset(Database* database);
};

/// A model providing all children of a given entity.
class ChildModel : public Support::StructModelAdapter<EntityNamePair> {
    Q_OBJECT
    QPointer<Database> m_host;

    Q_WRITABLE_PROPERTY(Entity, node, { });

    std::unordered_map<Entity, int> m_reverse;

    QVector<EntityNamePair> rebuild_lists();

private slots:
    void recompute();
    void record_changed(SolTrace::GUI::Data::Entity);

public:
    explicit ChildModel(QObject* parent = nullptr);
    virtual ~ChildModel() = default;

    /// Observe a database and rebuild the child list for node.
    void reset(Database* database);
};

} // namespace SolTrace::GUI::Data

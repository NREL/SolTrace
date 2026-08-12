#pragma once

#include "database/models/entity_name_model.h"

#include <QVariant>

namespace db {

/// A model providing all element entities that do not currently have a parent.
class RootElementsModel : public StructModelAdapter<EntityNamePair> {
    Q_OBJECT

    QPointer<Database> m_host;

    std::unordered_map<entt::entity, int> m_reverse;

    QVector<EntityNamePair> rebuild_lists();

private slots:
    void recompute();
    void record_changed(entt::entity);

public:
    explicit RootElementsModel(QObject* parent = nullptr);
    ~RootElementsModel() override = default;

    void reset(Database* database);

public slots:
    QVariant get(int index);
};

/// A model providing all element entities in a database.
class AllElementsModel : public StructModelAdapter<EntityNamePair> {
    Q_OBJECT

    QPointer<Database> m_host;

    std::unordered_map<entt::entity, int> m_reverse;

    QVector<EntityNamePair> rebuild_lists();

private slots:
    void recompute();
    void record_changed(entt::entity);

public:
    explicit AllElementsModel(QObject* parent = nullptr);
    ~AllElementsModel() override = default;

    /// Observe a database and rebuild the element list.
    void reset(Database* database);
};

} // namespace db

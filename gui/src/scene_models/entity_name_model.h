#pragma once

#include "data/database.h"
#include "support/qt_helpers.h"
#include "support/structmodel.h"

#include <QObject>
#include <QPointer>
#include <QString>

namespace SolTrace::GUI::Data {

/// An Entity + Name struct that can be exposed to QML.
struct EntityNamePair {
    QString name;
    Entity  entity;
    bool    has_children = false;

    RECORD_META(SolTrace::GUI::Data::EntityNamePair,
                SM_EXPOSE_RW(name),
                SM_EXPOSE_RO(entity),
                SM_EXPOSE_RO(has_children), );

    /// Build a display row for an entity using the database name hierarchy.
    static EntityNamePair record_for_entity(Database& db, SolTrace::GUI::Data::Entity entity);
};

/// Observe an entity's name.
class NameModel : public QObject {
    Q_OBJECT

    QPointer<Database> m_host;

    Entity m_target;

    Q_WRITABLE_PROPERTY(Entity, node, { });
    Q_READONLY_PROPERTY(QString, name);

private slots:
    void recompute(SolTrace::GUI::Data::Entity);

public:
    explicit NameModel(QObject* parent = nullptr);
    virtual ~NameModel() = default;

    /// Observe a database and recompute the current entity name.
    void reset(Database* database);

public slots:
    /// Rename the observed entity.
    void update_name(QString);
};

} // namespace SolTrace::GUI::Data

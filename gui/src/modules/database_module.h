#pragma once

#include "database/database.h"
#include "database/simulationresult.h"
#include "utilities/notification.h"
#include "utilities/qt_helpers.h"
#include "utilities/structmodel.h"

#include <QObject>
#include <QUrl>

#include <variant>

namespace SolTrace::GUI::App {

struct DatabaseRecord {
    QPointer<db::Database> database;

    RECORD_META(DatabaseRecord, SM_EXPOSE_RO(database), );
};

struct LoadedFile {
    // TODO: Store this in the DB
    QString       provenance = {};
    db::Database* ptr        = nullptr;
};

struct LoadFileFailed {
    ANotification notification;

    LoadFileFailed(QString message)
        : notification(ANotification::error(message)) { }
};

using LoadResult = std::variant<LoadedFile, LoadFileFailed>;

class DatabaseModule : public StructModelAdapter<DatabaseRecord> {
    Q_OBJECT
    QML_ELEMENT

    QOBJECT_WRITABLE_PROPERTY(db::Database, current_database)

    Q_WRITABLE_PROPERTY(bool, is_loading, false)

    void file_ready(QUrl, LoadResult);
    void file_failed(QUrl, QString);

public:
    DatabaseModule(QObject* parent = nullptr);

public slots:
    void load_url(QUrl, QString name_override = "");
    void load_new();

    bool set_current(int);

    void delete_current();
    void append_new(QString);
    bool append_clone(db::SimulationResultPtr);

signals:
    void notify(ANotification);
    void cancel_current_load(QPrivateSignal);
};

} // namespace SolTrace::GUI::App

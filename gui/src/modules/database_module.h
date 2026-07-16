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
    QString                       provenance = {};
    std::unique_ptr<db::Database> ptr;
};

struct LoadFileFailed {
    ANotification notification;

    LoadFileFailed(QString message)
        : notification(ANotification::error(message)) { }
};

// using LoadResult = std::variant<LoadedFile, LoadFileFailed>;

class DatabaseModule : public StructModelAdapter<DatabaseRecord> {
    Q_OBJECT
    QML_ELEMENT

    QOBJECT_WRITABLE_PROPERTY(db::Database, current_database)
    Q_PROPERTY(QUrl examples_folder READ examples_folder CONSTANT)
    Q_WRITABLE_PROPERTY(QString,
                        default_example_filename,
                        "Power-tower-surround_singlefacet.stinput")
    Q_PROPERTY(QUrl default_example READ default_example)

    QUrl examples_folder() const;
    QUrl default_example() const;

    Q_WRITABLE_PROPERTY(bool, is_loading, false)

    void file_ready(QUrl, LoadedFile);
    void file_failed(QUrl, LoadFileFailed);

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

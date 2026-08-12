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

/// One open database shown in the load/start workflow.
struct DatabaseRecord {
    QPointer<db::Database> database;

    RECORD_META(DatabaseRecord, SM_EXPOSE_RO(database), );
};

/// Successful result from asynchronous database loading.
struct LoadedFile {
    // TODO: Store this in the DB
    QString                       provenance = { };
    std::unique_ptr<db::Database> ptr;
};

/// Failed database load result packaged as a user notification.
struct LoadFileFailed {
    ANotification notification;

    LoadFileFailed(QString message)
        : notification(ANotification::error(message)) { }
};

/// QML-facing controller for opening, saving, creating, and switching
/// databases.
///
/// The model rows represent open databases. The module owns loaded database
/// instances and exposes the selected one through current_database.
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
    /// Open a native file picker and load the selected SolTrace input file.
    bool open_file_dialog();

    /// Open a native file picker and save the current database.
    bool save_current_dialog();

    /// Load a database from a URL. name_override is used for display only.
    void load_url(QUrl, QString name_override = "");

    /// Replace the current selection with a new blank database.
    void load_new();

    /// Select the open database at model row index.
    bool set_current(int);

    /// Save a specific open database to path.
    void save_db_at_index(int, QUrl);

    /// Save current_database to path.
    void save_current(QUrl path);

    /// Remove current_database from the open database list.
    void delete_current();

    /// Add a new blank database to the open database list.
    void append_new(QString);

    /// Clone a simulation result into an editable database.
    bool append_clone(db::SimulationResultPtr);

signals:
    void notify(ANotification);
    void cancel_current_load(QPrivateSignal);
};

} // namespace SolTrace::GUI::App

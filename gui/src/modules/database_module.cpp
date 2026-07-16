#include "database_module.h"
#include "utilities/asynctask.h"
#include "utilities/math_utility.h"
#include "utilities/result.h"

#include <QCoreApplication>
#include <QDir>
#include <QtConcurrent/qtconcurrentrun.h>
#include <QtCore/qfileinfo.h>
#include <QtCore/qfuturewatcher.h>

#include <exception>

namespace SolTrace::GUI::App {


static ::Result<LoadedFile, LoadFileFailed>
load_file(TaskControl& control, QString fname, db::Database* new_db) {

    // Take control of that free pointer...
    std::unique_ptr<db::Database> destination(new_db);
    QString                       stage = "starting";

    try {
        qDebug() << Q_FUNC_INFO << fname;

        stage = "reading file";
        control.setProgressValueAndText(0, "Reading file...");

        auto file = QFileInfo(fname);

        if (!(file.isFile() && file.isReadable())) {
            return ::return_failure(
                QString("Could not open the file for reading: %1").arg(fname));
        }

        auto new_data = std::make_shared<SD::SimulationData>();

        auto str = fname.toStdString();

        stage = "parsing file";
        if (!new_data->import_from_file(str)) {
            return return_failure(
                QString("Could not import the file: %1").arg(fname));
        }

        ASYNC_TASK_SYNC_POINT(control);

        stage = "importing content";
        control.setProgressValueAndText(50, "Importing content...");

        destination->import(*new_data);

        stage = "finalizing";
        control.setProgressValueAndText(100, "Done");

        // We cannot store non-copy types into Qt types, sigh.

        ASYNC_TASK_SYNC_POINT(control);

        return LoadedFile {
            .provenance = fname,
            .ptr        = std::move(destination),
        };
    } catch (std::exception const& e) {
        return return_failure(
            QString("Could not load %1 while %2: %3")
                .arg(fname, stage, QString::fromUtf8(e.what())));
    } catch (...) {
        return return_failure(
            QString("Could not load %1 while %2.").arg(fname, stage));
    }
}

QUrl DatabaseModule::examples_folder() const {
    QDir appDir(QCoreApplication::applicationDirPath());
#ifdef Q_OS_MACOS
    appDir.cdUp(); // Contents/
    appDir.cd("Resources/examples");
#else
    appDir.cd("examples"); // Linux/Windows: alongside binary
#endif
    return QUrl::fromLocalFile(appDir.absolutePath());
}

void DatabaseModule::file_ready(QUrl, LoadedFile result) {
    if (!result.ptr) {
        // Cancelled.
    } else {

        // Set owner for ptr...

        auto* database = result.ptr.release();

        database->setParent(this);
        store_push_append({ .database = database });

        notify(ANotification::info(
            QString("Loaded scene: %1").arg(database->name())));
    }

    set_is_loading(false);
}

void DatabaseModule::file_failed(QUrl, LoadFileFailed reason) {
    emit this->notify(reason.notification);
    set_is_loading(false);
}

void DatabaseModule::load_url(QUrl url, QString name_override) {
    if (is_loading()) {
        emit notify(ANotification::warning(
            "A file is already loading. Please wait for it to finish."));
        return;
        // emit cancel_current_load(QPrivateSignal {});
    }

    set_is_loading(false);

    auto new_source = url;

    if (new_source.isEmpty()) {
        qDebug() << Q_FUNC_INFO << "new database";

        if (name_override.isEmpty()) name_override = "Untitled";

        this->store_push_append(DatabaseRecord {
            .database = new db::Database(name_override, this),
        });

        return;
    }

    qDebug() << Q_FUNC_INFO << new_source;

    set_is_loading(true);

    auto fname = url.fileName();

    if (!name_override.isEmpty()) { fname = name_override; }

    if (fname.isEmpty()) { fname = "Untitled"; }

    // Now we have to do this in this round about way. If we have the thread
    // create the the database, there is no clean way to migrate it to our
    // thread. Setting parent, and setting current thread do not cut it.
    // Therefore we create it here, and send it to the thread for mutation.
    // Since there is nothing listening to signals, and the thread is the only
    // one with control of the database, this is essentially safe.

    // This also gets tricky, as we can't just wrap this in a shared pointer or
    // unique pointer and send it off to the async task, which takes copies.
    // So we just do a raw new, NOT giving it a parent, and immediately send it
    // to the task, which then wraps it.
    auto ptr = new db::Database(fname);

    auto task = launch_async_task<LoadedFile, LoadFileFailed>(
        url,
        this,
        &DatabaseModule::file_ready,
        &DatabaseModule::file_failed,
        load_file,
        new_source.toLocalFile(),
        std::move(ptr));

    connect(this,
            &DatabaseModule::cancel_current_load,
            task,
            &AsyncTaskBase::cancel);
}

void DatabaseModule::load_new() {
    load_url(QUrl());
}

DatabaseModule::DatabaseModule(QObject* parent)
    : StructModelAdapter { parent } {

    connect(this,
            &DatabaseModule::rowsInserted,
            this,
            [this](QModelIndex const& parent, int first, int last) {
                set_current(first);
            });
}

bool DatabaseModule::set_current(int index) {
    auto db = this->get_at(index);

    if (!db) return false;

    set_current_database(db->database);

    return true;
}

void DatabaseModule::delete_current() {
    auto const& v = this->vector();
    auto        iter =
        std::find_if(v.begin(), v.end(), [this](DatabaseRecord const& record) {
            return record.database == m_current_database;
        });

    if (iter == v.end()) return;

    auto index = std::distance(v.begin(), iter);

    db::Database* curr_cache = m_current_database;

    if (rowCount() == 1) {
        // this should mean that index == 0.
        load_new();
        set_current(1);
    } else {
        auto replacement_index = index;
        if (replacement_index > 0) { replacement_index--; }
        set_current(replacement_index);
    }

    this->store_push_remove(index, 1);

    if (curr_cache) curr_cache->deleteLater();
}

void DatabaseModule::append_new(QString new_name) {
    load_url({ }, new_name);
}

bool DatabaseModule::append_clone(db::SimulationResultPtr result) {
    if (!result || !result->database) {
        emit notify(ANotification::warning(
            "Select a simulation result before creating an editable copy."));
        return false;
    }

    auto clone_name = result->database->name() + " Copy";
    auto clone      = result->database->clone(clone_name, this);
    if (!clone) {
        emit notify(ANotification::error(
            "Could not create an editable copy of this result."));
        return false;
    }

    this->store_push_append(DatabaseRecord {
        .database = clone,
    });
    emit notify(ANotification::info(
        QString("Created editable scene: %1").arg(clone_name)));

    return true;
}

QUrl DatabaseModule::default_example() const {
    auto path =
        examples_folder().toLocalFile() + "/" + m_default_example_filename;
    return QUrl::fromLocalFile(path);
}

} // namespace SolTrace::GUI::App

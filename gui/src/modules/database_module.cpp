#include "database_module.h"
#include "utilities/asynctask.h"
#include "utilities/math_utility.h"

#include <QtConcurrent/qtconcurrentrun.h>
#include <QtCore/qfileinfo.h>
#include <QtCore/qfuturewatcher.h>

#include <exception>

namespace SolTrace::GUI::App {


using ResultFuture = QFutureWatcher<LoadResult>;

static void
load_file(QPromise<LoadResult>& result, QString fname, db::Database* new_db) {

    // Take control of that free pointer...
    std::unique_ptr<db::Database> destination(new_db);
    QString                       stage = "starting";

    try {
        result.setProgressRange(0, 100);
        qDebug() << Q_FUNC_INFO << fname;

        stage = "reading file";
        result.setProgressValueAndText(0, "Reading file...");

        auto file = QFileInfo(fname);
        auto suffix = file.suffix().toLower();
        bool legacy_stinput_cylinder_origins =
            suffix == QStringLiteral("stinput") ||
            suffix == QStringLiteral("stimport");

        if (!(file.isFile() && file.isReadable())) {
            result.emplaceResult(LoadFileFailed(
                QString("Could not open the file for reading: %1").arg(fname)));
            return;
        }

        auto new_data = std::make_shared<SD::SimulationData>();

        auto str = fname.toStdString();

        stage = "parsing file";
        if (!new_data->import_from_file(str)) {
            result.emplaceResult(
                LoadFileFailed(QString("Could not import the file: %1")
                                   .arg(fname)));
            return;
        }

        if (result.isCanceled()) {
            result.emplaceResult(LoadedFile {});
            return;
        }

        stage = "importing content";
        result.setProgressValueAndText(50, "Importing content...");

        destination->import(*new_data, legacy_stinput_cylinder_origins);

        stage = "finalizing";
        result.setProgressValueAndText(100, "Done");

        // We cannot store non-copy types into Qt types, sigh.

        if (result.isCanceled()) {
            result.emplaceResult(LoadedFile {});
            return;
        }

        result.emplaceResult(LoadedFile {
            .provenance = fname,
            .ptr        = destination.release(),
        });
    } catch (std::exception const& e) {
        result.emplaceResult(LoadFileFailed(
            QString("Could not load %1 while %2: %3")
                .arg(fname, stage, QString::fromUtf8(e.what()))));
    } catch (...) {
        result.emplaceResult(LoadFileFailed(
            QString("Could not load %1 while %2.")
                .arg(fname, stage)));
    }
}

void DatabaseModule::file_ready(QUrl, LoadResult result) {

    std::visit(overloaded {
                   [this](LoadedFile& arg) {
                       if (!arg.ptr) {
                           // Cancelled.
                       } else {
                           arg.ptr->setParent(this);
                           this->store_push_append(DatabaseRecord {
                               .database = arg.ptr,
                           });
                           this->notify(ANotification::info(
                               QString("Loaded scene: %1")
                                   .arg(arg.ptr->name())));
                       }
                   },
                   [this](LoadFileFailed failure) {
                       emit this->notify(failure.notification);
                   },
               },
               result);

    set_is_loading(false);
}

void DatabaseModule::file_failed(QUrl, QString reason) {
    emit notify(ANotification::error("Could not load the file: " + reason));
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

    auto task = launch_async_task<LoadResult>(url,
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
    load_url({}, new_name);
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

} // namespace SolTrace::GUI::App

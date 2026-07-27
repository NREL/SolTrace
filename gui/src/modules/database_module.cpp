#include "database_module.h"
#include "utilities/asynctask.h"
#include "utilities/math_utility.h"
#include "utilities/result.h"

#include <QCoreApplication>
#include <QDir>
#include <QFile>
#include <QFileInfo>
#include <QTemporaryFile>
#include <QTimer>
#include <QtGlobal>
#include <QtConcurrent/qtconcurrentrun.h>
#include <QtCore/qfileinfo.h>
#include <QtCore/qfuturewatcher.h>

#ifdef Q_OS_WASM
#include <QFileDialog>
#endif

#include <exception>

namespace SolTrace::GUI::App {

namespace {

QString import_name_filter() {
    return QStringLiteral(
        "SolTrace Files (*.stinput *.json);;All Files (*)");
}

QString wasm_temp_import_template(QString const& file_name) {
    auto suffix = QFileInfo(file_name).suffix();
    if (suffix.isEmpty()) { suffix = QStringLiteral("stinput"); }

    return QDir::tempPath() + QStringLiteral("/soltrace-import-XXXXXX.") +
           suffix;
}

#if defined(Q_OS_WASM) && !defined(__EMSCRIPTEN_PTHREADS__)
struct DirectTaskControl : TaskControl {
    void suspendIfRequested() const override { }
    bool cancelRequested() const override { return false; }

    void setProgressValue(int) const override { }
    void setProgressValueAndText(int, QString) const override { }
};
#endif

} // namespace

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

        bool legacy_import = false;

        if (str.ends_with("stinput")) {
            legacy_import = true;
            if (!new_data->import_from_file(str)) {
                return return_failure(
                    QString("Could not import the file: %1").arg(fname));
            }
        } else if (str.ends_with("json")) {
            try {
                new_data->import_json_file(str);
            } catch (std::exception const& e) {
                return return_failure(
                    QString("Could not import the file: %1 %2")
                        .arg(fname)
                        .arg(e.what()));
            }
        } else {
            return return_failure("Unknown file type.");
        }


        ASYNC_TASK_SYNC_POINT(control);

        stage = "importing content";
        control.setProgressValueAndText(50, "Importing content...");

        destination->import(*new_data, legacy_import);

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

#if defined(Q_OS_WASM) && !defined(__EMSCRIPTEN_PTHREADS__)
    auto local_path = new_source.toLocalFile();
    QTimer::singleShot(0, this, [this, url, local_path, ptr]() {
        DirectTaskControl control;
        auto              result = load_file(control, local_path, ptr);

        if (result) {
            file_ready(url, std::move(result.get_success()));
        } else {
            file_failed(url, std::move(result.get_failure()));
        }
    });
#else
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
#endif
}

void DatabaseModule::load_new() {
    load_url(QUrl());
}

bool DatabaseModule::open_file_dialog() {
#ifndef Q_OS_WASM
    return false;
#else
    if (is_loading()) {
        emit notify(ANotification::warning(
            "A file is already loading. Please wait for it to finish."));
        return true;
    }

    QFileDialog::getOpenFileContent(
        import_name_filter(),
        [this](QString const& file_name, QByteArray const& content) {
            if (file_name.isEmpty()) { return; }

            QTemporaryFile file(wasm_temp_import_template(file_name));
            file.setAutoRemove(false);

            if (!file.open()) {
                emit notify(ANotification::error(
                    QStringLiteral(
                        "Could not create a temporary import file: %1")
                        .arg(file.errorString())));
                return;
            }

            auto const bytes_written = file.write(content);
            if (bytes_written != content.size()) {
                emit notify(ANotification::error(
                    QStringLiteral("Could not stage the selected file: %1")
                        .arg(file.errorString())));
                return;
            }

            auto const temp_path = file.fileName();
            file.close();

            load_url(QUrl::fromLocalFile(temp_path),
                     QFileInfo(file_name).fileName());
        });

    return true;
#endif
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

static bool save_common(db::Database&   source,
                        QString         path,
                        DatabaseModule& notification,
                        bool            emit_success = true) {
    auto result = source.export_to_simdata();

    if (!result) {
        emit notification.notify(ANotification::error(
            QStringLiteral("Unable to save database. An error occurred while "
                           "packing content: %1")
                .arg(result.get_failure())));
        return false;
    }

    auto pack = result.get_success();

    try {
        pack->data->export_json_file(path.toStdString());
    } catch (std::exception const& ex) {
        emit notification.notify(ANotification::error(
            QStringLiteral(
                "An exception occurred while trying to save content: %1")
                .arg(ex.what())));

        return false;
    }

    if (emit_success) {
        emit notification.notify(
            ANotification::info(QStringLiteral("File successfully saved.")));
    }

    return true;
}

void DatabaseModule::save_db_at_index(int index, QUrl path) {
    auto db = this->get_at(index);

    if (!db or !db->database) {
        notify(ANotification::error(QStringLiteral(
            "An internal error was encountered trying to save the scene.")));
        return;
    }

    save_common(*(db->database), path.toLocalFile(), *this);
}

void DatabaseModule::save_current(QUrl path) {
    if (!m_current_database) {
        notify(ANotification::error(QStringLiteral(
            "An internal error was encountered trying to save the scene.")));
        return;
    }

    save_common(*m_current_database, path.toLocalFile(), *this);
}

bool DatabaseModule::save_current_dialog() {
#ifndef Q_OS_WASM
    return false;
#else
    if (!m_current_database) {
        notify(ANotification::error(QStringLiteral(
            "An internal error was encountered trying to save the scene.")));
        return true;
    }

    QTemporaryFile file(QDir::tempPath() +
                        QStringLiteral("/soltrace-export-XXXXXX.json"));
    file.setAutoRemove(true);

    if (!file.open()) {
        emit notify(ANotification::error(
            QStringLiteral("Could not create a temporary export file: %1")
                .arg(file.errorString())));
        return true;
    }

    auto const temp_path = file.fileName();
    file.close();

    if (!save_common(*m_current_database,
                     temp_path,
                     *this,
                     false /* emit_success */)) {
        return true;
    }

    QFile exported_file(temp_path);
    if (!exported_file.open(QIODevice::ReadOnly)) {
        emit notify(ANotification::error(
            QStringLiteral("Could not read the exported scene: %1")
                .arg(exported_file.errorString())));
        return true;
    }

    auto const content = exported_file.readAll();
    QFileDialog::saveFileContent(
        content, m_current_database->name() + QStringLiteral(".json"));

    emit notify(
        ANotification::info(QStringLiteral("File successfully saved.")));

    return true;
#endif
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

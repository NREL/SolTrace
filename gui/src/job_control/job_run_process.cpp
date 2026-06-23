#include "job_run_process.h"

#include <QCoreApplication>
#include <QJsonArray>
#include <QJsonDocument>
#include <QJsonObject>
#include <QProcess>
#include <QString>
#include <QTemporaryDir>
#include <QTimer>

#include "job_run.h"

// Stdin polling is not very clean across platforms. We'll use a control file
// for now


constexpr inline static const char* SOLTRACE_WORKER_ENV   = "SOLTRACE_WORKER";
constexpr inline static const char* SOLTRACE_JOB_MANIFEST = "dataset.json";
constexpr inline static const char* SOLTRACE_JOB_RESULT   = "result.json";
constexpr inline static const char* SOLTRACE_JOB_CONTROL  = "control.json";

static void spin_off_messages(QPromise<SimResult>& promise, QByteArray array) {
    auto stream = QTextStream(array);

    QString buffer;

    int last_progress = 0;

    while (stream.readLineInto(&buffer)) {
        auto parts = buffer.split(':');
        buffer.clear();

        auto first  = parts.value(0).trimmed();
        auto second = parts.value(1);

        if (first == "progress") {

            bool ok;
            int  temp_last_progress = second.toInt(&ok);

            if (!ok) {
                qWarning() << "Unable to decode progress" << second;
                continue;
            }

            last_progress = temp_last_progress;


            promise.setProgressValue(last_progress);
        } else if (first == "message") {
            promise.setProgressValueAndText(last_progress, second);
        } else if (first == "done") {
            promise.setProgressValueAndText(100, "Done");
        } else if (first == "error") {
            promise.emplaceResult(second);
        } else {
            qWarning() << "Unknown message from worker" << first;
        }
    }
}

void execute_process_runner(QPromise<SimResult>&      promise,
                            SimDataPtr                data,
                            ThreadRunnerConfig const& config) {
    try {
        promise.setProgressRange(0, 100);

        auto temp_dir = QTemporaryDir();
        temp_dir.setAutoRemove(true);

        auto source = temp_dir.filePath(SOLTRACE_JOB_MANIFEST);

        // throws on write fail
        data->data->export_json_file(source.toStdString());

        auto process = QProcess();

        auto env = QProcessEnvironment::systemEnvironment();

        env.insert(SOLTRACE_WORKER_ENV, temp_dir.path());

        process.setProcessEnvironment(env);
        process.setReadChannel(QProcess::StandardOutput);

        process.start(qApp->applicationFilePath());

        if (!process.waitForStarted()) {
            promise.emplaceResult(
                QStringLiteral("Could not start the simulation worker."));
            return;
        }

        // poll for progess and status

        while (true) {

            if (process.waitForFinished(16)) { break; }

            // not done

            if (process.waitForReadyRead(0)) {
                // get any content

                spin_off_messages(promise, process.readAll());
            }
        }

        // process is finished

        // read in content


    } catch (std::exception& e) {
        promise.emplaceResult(QString(e.what()));
        return;
    }
}


static SimDataPtr worker_start_flow(QByteArray work_dir) {
    // If not a worker, bail, return to normal app flow

    auto source = QDir(work_dir).absoluteFilePath(SOLTRACE_JOB_MANIFEST);

    // From here on out, we use stdout for all messages
    auto write_error = [](QString message) {
        // not clean, but this isnt perf critical
        std::cout << "error:" << message.toStdString() << std::endl;
    };

    if (!QFileInfo::exists(source)) {
        write_error("Manifest missing");
        return nullptr;
    }

    auto dataset = std::make_shared<SD::SimulationData>();

    if (!dataset->import_from_file(source.toStdString())) {
        write_error("Unable to read manifest");
        return nullptr;
    }

    auto database = db::Database("Temporary");

    database.import(*dataset);

    // TODO: Non native runner
    return database.export_to_simdata().get_success();
}

inline QJsonArray to_json(glm::vec3 v) {
    return QJsonArray() << v.x << v.y << v.z;
}

static QJsonArray build_results(db::SimulationResultPtr results) {
    // if (!results) { return {}; }

    // auto& result = results->result;

    // QJsonArray array;

    // for (auto iter = result.get_ray_record_iteratior();
    // !result.is_at_end(iter);
    //      ++iter) {
    //     QJsonArray record;

    //     auto const& collection = *(iter->get());

    //     for (auto const& interaction : collection.interactions) {
    //         record.push_back(QJsonObject {
    //             { QStringLiteral("element"), interaction->element },
    //             { QStringLiteral("event"),
    //               static_cast<int>(interaction->event) },
    //             { QStringLiteral("location"), to_json(interaction->location)
    //             }, { QStringLiteral("direction"),
    //               to_json(interaction->direction) },
    //         });
    //     }

    //     array.push_back(QJsonObject {
    //         { QStringLiteral("id"), collection.id },
    //         { QStringLiteral("collection"), record },
    //     });
    // }

    // return array;

    // TODO: Disabled till we can get process running working
    abort();
}

static bool dump_results(db::SimulationResultPtr results, QByteArray work_dir) {
    auto obj = build_results(results);

    auto content = QJsonDocument(obj).toJson();

    auto dest = QDir(work_dir).absoluteFilePath(SOLTRACE_JOB_RESULT);

    auto dest_file = QFile(dest);

    if (!dest_file.open(QFile::WriteOnly)) { return false; }

    dest_file.write(content);

    return true;
}

void check_if_process_worker(int argc, char* argv[]) {
    QCoreApplication app(argc, argv);

    auto work_dir = qgetenv(SOLTRACE_WORKER_ENV);

    if (work_dir.isEmpty()) return;

    qDebug() << "This is a worker process";

    auto control      = QDir(work_dir).absoluteFilePath(SOLTRACE_JOB_CONTROL);
    auto is_cancelled = [&]() { return QFileInfo::exists(control); };

    auto sim_data = worker_start_flow(work_dir);

    // This will delete itself on done
    auto ptr = new RunningJob(sim_data, RunType::Thread, 0);

    auto timer = new QTimer(ptr);

    timer->setInterval(16);
    timer->start();

    QObject::connect(timer, &QTimer::timeout, [=]() {
        if (is_cancelled()) { ptr->cancel(); }
    });

    QObject::connect(ptr, &RunningJob::progress_update, [](int percent) {
        std::cout << "progress:" << percent << std::endl;
    });

    QObject::connect(
        ptr, &RunningJob::progress_text_update, [](QString message) {
            std::cout << "message:" << message.toStdString() << std::endl;
        });

    QObject::connect(ptr, &RunningJob::finished, [ptr, work_dir]() {
        std::cout << "done" << std::endl;

        // write to disk
        auto p  = ptr->take();
        bool ok = dump_results(p, work_dir);

        qApp->exit(ok ? EXIT_SUCCESS : EXIT_FAILURE);
    });

    QObject::connect(ptr, &RunningJob::error, [](QString message) {
        std::cout << "error:" << message.toStdString() << std::endl;
    });

    app.exec();
}

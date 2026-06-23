#pragma once

#include "database/database.h"
#include "database/simulationresult.h"
#include "database/worldgeometrymodel.h"
#include "job_control/job_run.h"
#include "utilities/notification.h"
#include "utilities/qt_helpers.h"
#include "utilities/structmodel.h"

#include "module_common.h"

#include <QObject>
#include <QQmlEngine>
#include <QString>

namespace SolTrace::GUI::App {

// TODO: Track added simulation results and allow deletion!

class SimulationRunnerModel;

/**
 * @class SimulationModule
 * @brief Simulation execution and progress tracking module.
 *
 * Mediates between QML controls and the job runner.
 * Exposes progress, timing metadata, and execution control to QML.
 *
 * start(), stop(), and pause() delegate to the backend after validating
 * that all required configuration modules are in a Ready or Complete state.
 *
 * QML access pattern: App.simulation.start()
 */
class SimulationModule : public QObject {
    Q_OBJECT

    QPointer<RunningJob> m_running;

    db::SimulationResultPtr m_current_result;

    QVector<std::shared_ptr<db::SimulationResult>> m_completed_sims;

private slots:
    void job_done();
    void job_failed(QString const& message);
    void update_result_world(db::SimulationResultPtr);

public:
    explicit SimulationModule(QObject* parent = nullptr);
    ~SimulationModule();

    QOBJECT_WRITABLE_PROPERTY(db::Database, current_database)
    QOBJECT_READONLY_PROPERTY(StatusComponent, status);
    QOBJECT_READONLY_PROPERTY(SimulationRunnerModel, runners);
    QOBJECT_READONLY_PROPERTY(db::SimulationResultModel, results);
    QOBJECT_READONLY_PROPERTY(db::WorldGeometryModel, world_geometry_model);

    enum Runner { CPU = 0, Embree = 1, GPU = 2 };

    Q_ENUM(Runner)

    Q_WRITABLE_PROPERTY(Runner, runner, Runner::CPU);
    Q_WRITABLE_PROPERTY(uint32_t, ray_count, 10000);
    Q_WRITABLE_PROPERTY(uint32_t, max_ray_count, 100000);
    Q_WRITABLE_PROPERTY(uint32_t, max_threads, 10);
    Q_WRITABLE_PROPERTY(uint32_t, seed_value, 1234)

    Q_WRITABLE_PROPERTY(bool, sun_shape, false)
    Q_WRITABLE_PROPERTY(bool, optical_errors, false)
    Q_WRITABLE_PROPERTY(bool, point_focus_system, false)

    Q_WRITABLE_PROPERTY(QString,
                        current_simulation_result_name,
                        "No Simulation Result")


    /// Is a simulation being run?
    Q_READONLY_PROPERTY(bool, is_running)

    /// Ray tracing progress, 0 to 100
    Q_READONLY_PROPERTY(int, progress)

    /// e.g. "Initializing", "Ray tracing", "Complete"
    Q_READONLY_PROPERTY(QString, current_stage)
    Q_READONLY_PROPERTY(QDateTime, last_run_time)
    Q_READONLY_PROPERTY(double, elapsed_seconds)

public slots:
    void run();
    // void pause(); // no executor support for pause or resume
    // void resume();
    void cancel();
    void select_result(int index);
    void delete_result(int index);
    void rename_result(int index, QString const& name);
    void export_result(int index);
    void duplicate_current_result_for_edit();

signals:
    void new_results(db::SimulationResultPtr);
    void edit_result_copy_requested(db::SimulationResultPtr);
    void notify(ANotification);
};

struct SimulationRunnerRecord {
    QString                  name;
    SimulationModule::Runner runner;

    RECORD_META(SimulationRunnerRecord,
                SM_EXPOSE_RO(name),
                SM_EXPOSE_RO(runner), );
};

class SimulationRunnerModel
    : public StructModelAdapter<SimulationRunnerRecord> {
    Q_OBJECT

public:
    explicit SimulationRunnerModel(QObject* parent = nullptr);

public slots:
    SimulationModule::Runner runner_at(int index) const;
    int                      index_of(SimulationModule::Runner runner) const;
};

} // namespace SolTrace::GUI::App

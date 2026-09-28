#pragma once

#include "data/database.h"
#include "data/simulationresult.h"
#include "jobs/job_run.h"
#include "scene_models/world_geometry_model.h"
#include "support/notification.h"
#include "support/qt_helpers.h"
#include "support/structmodel.h"

#include <QObject>
#include <QQmlEngine>
#include <QString>
#include <QVariantMap>
#include <QVector3D>

namespace SolTrace::GUI::App {

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

    QPointer<Jobs::RunningJob> m_running;

    Data::SimulationResultPtr m_current_result;

    QVector<std::shared_ptr<Data::SimulationResult>> m_completed_sims;

    uint32_t m_running_requested_ray_count     = 0;
    uint32_t m_running_requested_max_ray_count = 0;

private slots:
    void job_done();
    void job_failed(QString const& message);
    void update_result_world(Data::SimulationResultPtr);

private:
    Jobs::ThreadRunnerBackend                  selected_backend() const;
    uint32_t                                   effective_thread_count() const;
    Support::Result<Jobs::SimDataPtr, QString> prepare_simulation_data();
    void connect_running_job(Jobs::RunningJob* job);
    void publish_completed_result(Data::SimulationResultPtr results);

public:
    explicit SimulationModule(QObject* parent = nullptr);
    ~SimulationModule();

    QOBJECT_WRITABLE_PROPERTY(Data::Database, current_database)
    QOBJECT_READONLY_PROPERTY(SimulationRunnerModel, runners);
    QOBJECT_READONLY_PROPERTY(Data::SimulationResultModel, results);
    QOBJECT_READONLY_PROPERTY(Data::WorldGeometryModel, world_geometry_model);
    Q_READONLY_PROPERTY(QVector3D, result_sun_position)
    Q_READONLY_PROPERTY(bool, result_sun_is_point_source)

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
    /// Compute current simulation ray-result bounds on demand.
    Q_INVOKABLE QVariantMap current_result_bounds() const;

    void run();
    // void pause(); // no executor support for pause or resume
    // void resume();
    void cancel();
    void select_result(int index);
    void delete_result(int index);
    void rename_result(int index, QString const& name);
    void export_result(int index);
    void duplicate_current_result_for_edit();

    void update_ray_count(int new_count);
    void update_max_ray_count(int new_max);
signals:
    void new_results(SolTrace::GUI::Data::SimulationResultPtr);
    void edit_result_copy_requested(SolTrace::GUI::Data::SimulationResultPtr);
    void notify(SolTrace::GUI::Support::ANotification);
};

struct SimulationRunnerRecord {
    QString                  name;
    SimulationModule::Runner runner;

    RECORD_META(SimulationRunnerRecord,
                SM_EXPOSE_RO(name),
                SM_EXPOSE_RO(runner), );
};

class SimulationRunnerModel
    : public Support::StructModelAdapter<SimulationRunnerRecord> {
    Q_OBJECT

public:
    explicit SimulationRunnerModel(QObject* parent = nullptr);

public slots:
    SimulationModule::Runner runner_at(int index) const;
    int                      index_of(SimulationModule::Runner runner) const;
};

} // namespace SolTrace::GUI::App

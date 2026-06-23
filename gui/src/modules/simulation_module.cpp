#include "simulation_module.h"

#include <algorithm>
#include <thread>

namespace SolTrace::GUI::App {

SimulationRunnerModel::SimulationRunnerModel(QObject* parent)
    : StructModelAdapter { parent } {
    store_push_append(SimulationRunnerRecord {
        .name   = "Legacy Runner",
        .runner = SimulationModule::CPU,
    });

#ifdef SOLTRACE_HAS_EMBREE_RUNNER
    store_push_append(SimulationRunnerRecord {
        .name   = "Embree Runner",
        .runner = SimulationModule::Embree,
    });
#endif

#ifdef SOLTRACE_HAS_OPTIX_RUNNER
    store_push_append(SimulationRunnerRecord {
        .name   = "GPU Runner",
        .runner = SimulationModule::GPU,
    });
#endif
}

SimulationModule::Runner SimulationRunnerModel::runner_at(int index) const {
    auto record = get_at(index);
    if (!record) return SimulationModule::CPU;

    return record->runner;
}

int SimulationRunnerModel::index_of(SimulationModule::Runner runner) const {
    int index = 0;
    for (auto const& record : *this) {
        if (record.runner == runner) return index;
        ++index;
    }

    return 0;
}

void SimulationModule::update_result_world(db::SimulationResultPtr results) {
    auto* database =
        results ? const_cast<db::Database*>(results->database.get()) : nullptr;
    m_world_geometry_model->reset(database);
}

void SimulationModule::job_done() {
    qDebug() << Q_FUNC_INFO;

    auto* from = qobject_cast<RunningJob*>(sender());
    if (!from) {
        qCritical() << Q_FUNC_INFO << "bad cast";
        return;
    }

    if (m_running != from) {
        qCritical() << Q_FUNC_INFO << m_running << from;
        return;
    }

    set_is_running(false);

    m_running = nullptr;

    auto results = from->take();

    if (!results) {
        emit notify(
            ANotification::error(
                "The simulation completed, but no results were produced."));
        return;
    }

    m_completed_sims.push_back(results);
    m_results->append_result(results);
    m_current_result = results;
    set_current_simulation_result_name(results->database->name());

    qDebug() << Q_FUNC_INFO << "publish";
    emit new_results(results);
}

void SimulationModule::job_failed(QString const& message) {
    qDebug() << Q_FUNC_INFO << message;

    auto* from = qobject_cast<RunningJob*>(sender());
    if (!from) {
        qCritical() << Q_FUNC_INFO << "bad cast";
        return;
    }

    if (m_running != from) {
        qCritical() << Q_FUNC_INFO << m_running << from;
        return;
    }

    set_is_running(false);
    m_running = nullptr;
    set_progress(0);
    set_current_stage("Idle");

    if (message.startsWith(QStringLiteral("Cancelled"))) {
        emit notify(ANotification::info("Simulation cancelled."));
    } else {
        emit notify(ANotification::error(message));
    }
}

SimulationModule::SimulationModule(QObject* parent)
    : QObject { parent },
      m_status(new StatusComponent(this)),
      m_runners(new SimulationRunnerModel(this)),
      m_results(new db::SimulationResultModel(this)),
      m_world_geometry_model(new db::WorldGeometryModel(this)) {

    auto thread_count = std::thread::hardware_concurrency();
    set_max_threads(thread_count <= 0 ? 1 : thread_count);

#ifdef SOLTRACE_HAS_EMBREE_RUNNER
    set_runner(Runner::Embree);
#endif

    connect(this,
            &SimulationModule::new_results,
            this,
            &SimulationModule::update_result_world);

    qDebug() << Q_FUNC_INFO;
}

SimulationModule::~SimulationModule() {
    qDebug() << Q_FUNC_INFO;
}

void SimulationModule::run() {
    qDebug() << Q_FUNC_INFO;
    if (!m_current_database) {
        emit notify(ANotification::warning(
            "Create or load a scene before running a simulation."));
        return;
    }

    if (m_running) {
        emit notify(ANotification::error(
            "A simulation is already running. Wait for it to finish before "
            "starting another one."));
        return;
    }

    qDebug() << Q_FUNC_INFO << "Launch";
    auto exported_result = m_current_database->export_to_simdata();

    if (!exported_result) {
        auto err = exported_result.get_failure();

        emit notify(ANotification::error(
            QString("Could not start the simulation: %1").arg(err)));
        return;
    }

    auto sim_data = exported_result.get_success();

    sim_data->data->set_seed(m_seed_value);
    sim_data->data->set_number_of_rays(m_ray_count);
    sim_data->data->set_max_rays_traced(m_max_ray_count);

    auto& sim_params = sim_data->data->get_simulation_parameters();
    sim_params.include_sun_shape_errors = m_sun_shape;
    sim_params.include_optical_errors   = m_optical_errors;

    auto backend = ThreadRunnerBackend::Native;
#ifdef SOLTRACE_HAS_EMBREE_RUNNER
    if (m_runner == Runner::Embree) { backend = ThreadRunnerBackend::Embree; }
#endif

    qDebug() << Q_FUNC_INFO << magic_enum::enum_name(backend);

    m_running =
        new RunningJob(sim_data, RunType::Thread, m_max_threads, backend, this);

    connect(m_running,
            &RunningJob::progress_update,
            this,
            &SimulationModule::set_progress);
    connect(m_running,
            &RunningJob::progress_text_update,
            this,
            &SimulationModule::set_current_stage);

    connect(
        m_running, &RunningJob::finished, this, &SimulationModule::job_done);
    connect(
        m_running, &RunningJob::error, this, &SimulationModule::job_failed);
    connect(
        m_running, &RunningJob::finished, m_running, &RunningJob::deleteLater);
    connect(m_running, &RunningJob::error, m_running, &RunningJob::deleteLater);
    connect(this, &QObject::destroyed, m_running, &RunningJob::cancel);

    set_is_running(true);
}

void SimulationModule::cancel() {
    if (m_running) { m_running->cancel(); }
}

void SimulationModule::select_result(int index) {
    auto result = m_results->result_at(index);
    if (!result) return;

    m_current_result = result;
    set_current_simulation_result_name(m_results->name_at(index));
    emit new_results(result);
}

void SimulationModule::delete_result(int index) {
    auto result = m_results->result_at(index);
    if (!result) return;

    const bool deleting_current = m_current_result == result;

    m_completed_sims.erase(std::remove(m_completed_sims.begin(),
                                       m_completed_sims.end(),
                                       result),
                           m_completed_sims.end());
    m_results->remove_result(index);

    if (!deleting_current) return;

    db::SimulationResultPtr replacement;
    QString                 replacement_name = "No Simulation Result";
    if (m_results->rowCount() > 0) {
        auto replacement_index = std::min(index, m_results->rowCount() - 1);
        replacement           = m_results->result_at(replacement_index);
        replacement_name      = m_results->name_at(replacement_index);
    }

    m_current_result = replacement;
    set_current_simulation_result_name(replacement_name);
    emit new_results(replacement);
}

void SimulationModule::rename_result(int index, QString const& name) {
    auto result = m_results->result_at(index);
    if (!result) return;

    m_results->rename_result(index, name);

    if (m_current_result == result) {
        set_current_simulation_result_name(name);
    }
}

void SimulationModule::export_result(int index) {
    if (!m_results->result_at(index)) return;

    emit notify(ANotification::info("Result export is not available yet."));
}

void SimulationModule::duplicate_current_result_for_edit() {
    if (!m_current_result) return;

    emit edit_result_copy_requested(m_current_result);
}

} // namespace SolTrace::GUI::App

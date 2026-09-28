#include "simulation_module.h"

#include <algorithm>
#include <thread>
#include <QLoggingCategory>
#include <QVariantMap>

using namespace SolTrace::GUI::Support;

namespace SolTrace::GUI::App {

Q_LOGGING_CATEGORY(simulationLog, "soltrace.gui.simulation")

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

void SimulationModule::update_result_world(SolTrace::GUI::Data::SimulationResultPtr results) {
    auto* database =
        results ? const_cast<SolTrace::GUI::Data::Database*>(results->database.get()) : nullptr;
    m_world_geometry_model->reset(database);

    QVector3D sun_position(0.0f, 0.0f, 1.0f);
    bool      is_point_source = false;

    if (database) {
        auto const* resource = database->ray_source_resource.get();
        if (resource) {
            is_point_source = resource->type == SolTrace::GUI::Data::RaySourceType::PointSource;

            if (resource->source) {
                auto const& position = resource->source->get_position();
                sun_position = QVector3D(position.x, position.y, position.z);
            }
        }
    }

    set_result_sun_position(sun_position);
    set_result_sun_is_point_source(is_point_source);
}

void SimulationModule::job_done() {
    qCDebug(simulationLog) << Q_FUNC_INFO;

    auto* from = qobject_cast<SolTrace::GUI::Jobs::RunningJob*>(sender());
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
        emit notify(ANotification::error(
            "The simulation completed, but no results were produced."));
        return;
    }

    publish_completed_result(results);
}

void SimulationModule::publish_completed_result(
    SolTrace::GUI::Data::SimulationResultPtr results) {
    if (m_running_requested_max_ray_count > 0 &&
        results->sun_ray_count >= m_running_requested_max_ray_count &&
        static_cast<uint64_t>(results->records.size()) <
            m_running_requested_ray_count) {
        emit notify(ANotification::warning(QString(
            "Simulation reached the maximum traced ray limit (%1) before "
            "collecting the requested %2 result rays. Increase Max # Rays "
            "Traced to collect more rays.")
                                               .arg(m_running_requested_max_ray_count)
                                               .arg(m_running_requested_ray_count)));
    } else {
        emit notify(ANotification::info(QString(
            "A simulation has completed. Check analysis for new data.")));
    }

    m_completed_sims.push_back(results);
    m_results->append_result(results);
    m_current_result = results;
    set_current_simulation_result_name(results->database->name());

    qCDebug(simulationLog) << Q_FUNC_INFO << "publish";
    emit new_results(results);
}

void SimulationModule::job_failed(QString const& message) {
    qCDebug(simulationLog) << Q_FUNC_INFO << message;

    auto* from = qobject_cast<SolTrace::GUI::Jobs::RunningJob*>(sender());
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
      m_runners(new SimulationRunnerModel(this)),
      m_results(new SolTrace::GUI::Data::SimulationResultModel(this)),
      m_world_geometry_model(new SolTrace::GUI::Data::WorldGeometryModel(this)) {

    auto thread_count = std::thread::hardware_concurrency();
    set_max_threads(thread_count <= 0 ? 1 : thread_count);

#ifdef Q_OS_WASM
    set_max_threads(1);
    set_ray_count(500);
    set_max_ray_count(500);
#endif

#ifdef SOLTRACE_HAS_EMBREE_RUNNER
    set_runner(Runner::Embree);
#endif

    connect(this,
            &SimulationModule::new_results,
            this,
            &SimulationModule::update_result_world);

    qCDebug(simulationLog) << Q_FUNC_INFO;
}

SimulationModule::~SimulationModule() {
    qCDebug(simulationLog) << Q_FUNC_INFO;
}

QVariantMap SimulationModule::current_result_bounds() const {
    const bool valid = m_current_result &&
                       m_current_result->bounds_min.x <
                           m_current_result->bounds_max.x &&
                       m_current_result->bounds_min.y <
                           m_current_result->bounds_max.y &&
                       m_current_result->bounds_min.z <
                           m_current_result->bounds_max.z;

    return QVariantMap {
        { QStringLiteral("valid"), valid },
        { QStringLiteral("min"),
          valid ? QVector3D(m_current_result->bounds_min.x,
                            m_current_result->bounds_min.y,
                            m_current_result->bounds_min.z)
                : QVector3D() },
        { QStringLiteral("max"),
          valid ? QVector3D(m_current_result->bounds_max.x,
                            m_current_result->bounds_max.y,
                            m_current_result->bounds_max.z)
                : QVector3D() },
    };
}

void SimulationModule::run() {
    qCDebug(simulationLog) << Q_FUNC_INFO;
    if (!m_current_database) {
        emit notify(ANotification::warning(
            "Create or Load Scene before running a simulation."));
        return;
    }

    if (m_running) {
        emit notify(ANotification::error(
            "A simulation is already running. Wait for it to finish before "
            "starting another one."));
        return;
    }

    qCDebug(simulationLog) << Q_FUNC_INFO << "Launch";
    auto exported_result = prepare_simulation_data();

    if (!exported_result) {
        auto err = exported_result.get_failure();

        emit notify(ANotification::error(
            QString("Could not start the simulation: %1").arg(err)));
        return;
    }

    auto sim_data = exported_result.get_success();

    m_running_requested_ray_count     = m_ray_count;
    m_running_requested_max_ray_count = m_max_ray_count;

    auto backend = selected_backend();
    qCDebug(simulationLog) << Q_FUNC_INFO << magic_enum::enum_name(backend);

    m_running =
        new SolTrace::GUI::Jobs::RunningJob(sim_data, effective_thread_count(), backend, this);
    connect_running_job(m_running);
    set_is_running(true);
}

SolTrace::GUI::Jobs::ThreadRunnerBackend SimulationModule::selected_backend() const {
    auto backend = SolTrace::GUI::Jobs::ThreadRunnerBackend::Native;
#ifdef SOLTRACE_HAS_EMBREE_RUNNER
    if (m_runner == Runner::Embree) { backend = SolTrace::GUI::Jobs::ThreadRunnerBackend::Embree; }
#endif
#ifdef SOLTRACE_HAS_OPTIX_RUNNER
    if (m_runner == Runner::GPU) { backend = SolTrace::GUI::Jobs::ThreadRunnerBackend::Optix; }
#endif
    return backend;
}

uint32_t SimulationModule::effective_thread_count() const {
#ifdef Q_OS_WASM
    return 1;
#else
    return m_max_threads;
#endif
}

::Result<SolTrace::GUI::Jobs::SimDataPtr, QString>
SimulationModule::prepare_simulation_data() {
    auto exported_result = m_current_database->export_to_simdata();
    if (!exported_result) { return exported_result; }

    auto sim_data = exported_result.get_success();

    sim_data->data->set_seed(m_seed_value);
    sim_data->data->set_number_of_rays(m_ray_count);
    sim_data->data->set_max_rays_traced(m_max_ray_count);

    auto& sim_params = sim_data->data->get_simulation_parameters();
    sim_params.include_sun_shape_errors = m_sun_shape;
    sim_params.include_optical_errors   = m_optical_errors;

    return sim_data;
}

void SimulationModule::connect_running_job(SolTrace::GUI::Jobs::RunningJob* job) {
    connect(job,
            &SolTrace::GUI::Jobs::RunningJob::progress_update,
            this,
            &SimulationModule::set_progress);
    connect(job,
            &SolTrace::GUI::Jobs::RunningJob::progress_text_update,
            this,
            &SimulationModule::set_current_stage);

    connect(job, &SolTrace::GUI::Jobs::RunningJob::finished, this, &SimulationModule::job_done);
    connect(job, &SolTrace::GUI::Jobs::RunningJob::error, this, &SimulationModule::job_failed);
    connect(job, &SolTrace::GUI::Jobs::RunningJob::finished, job, &SolTrace::GUI::Jobs::RunningJob::deleteLater);
    connect(job, &SolTrace::GUI::Jobs::RunningJob::error, job, &SolTrace::GUI::Jobs::RunningJob::deleteLater);
    connect(this, &QObject::destroyed, job, &SolTrace::GUI::Jobs::RunningJob::cancel);
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

    m_completed_sims.erase(
        std::remove(m_completed_sims.begin(), m_completed_sims.end(), result),
        m_completed_sims.end());
    m_results->remove_result(index);

    if (!deleting_current) return;

    SolTrace::GUI::Data::SimulationResultPtr replacement;
    QString                 replacement_name = "No Simulation Result";
    if (m_results->rowCount() > 0) {
        auto replacement_index = std::min(index, m_results->rowCount() - 1);
        replacement            = m_results->result_at(replacement_index);
        replacement_name       = m_results->name_at(replacement_index);
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

void SimulationModule::update_ray_count(int new_count) {
    if (new_count > max_ray_count()) { set_max_ray_count(new_count); }
    set_ray_count(new_count);
}

void SimulationModule::update_max_ray_count(int new_max) {
    if (new_max < ray_count()) { set_ray_count(new_max); }
    set_max_ray_count(new_max);
}

} // namespace SolTrace::GUI::App

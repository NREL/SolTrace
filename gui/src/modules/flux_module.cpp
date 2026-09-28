#include "flux_module.h"
#include "analysis/rays/ray_volume_raster.h"
#include "analysis/volume/volume_to_mesh.h"
#include "data/components.h"
#include "support/asynctask.h"

#include <QQmlEngine>
#include <QLoggingCategory>
#include <QUuid>

using namespace SolTrace::GUI::Support;

namespace SolTrace::GUI::App {

Q_LOGGING_CATEGORY(fluxLog, "soltrace.gui.flux")

FluxModule::FluxModule(QQmlEngine* engine, QObject* parent)
    : QObject(parent),
      m_entity_model(new SolTrace::GUI::Data::AllElementsModel(this)),
      m_computed_maps_model(new SolTrace::GUI::Data::AllComputedMapsModel(this)),
      m_pending_flux_maps(new SolTrace::GUI::Data::PendingFluxMapModel(this)),
      m_flux_map_world_model(new SolTrace::GUI::Data::FluxMapWorldModel(this)),
      m_ray_iso_volume(new SolTrace::GUI::Data::QMLMesh()) {

    set_ray_volume_flux_in_progress(false);

    m_ray_iso_volume->setParent(this);

    m_image_provider = m_pending_flux_maps->make_new_provider();

    engine->addImageProvider("fluxmap", m_image_provider);

    connect(m_pending_flux_maps,
            &SolTrace::GUI::Data::PendingFluxMapModel::ready,
            m_flux_map_world_model,
            &SolTrace::GUI::Data::FluxMapWorldModel::on_ready);

    connect(m_pending_flux_maps,
            &SolTrace::GUI::Data::PendingFluxMapModel::ready,
            this,
            &FluxModule::flux_map_ready);

    connect(m_pending_flux_maps,
            &SolTrace::GUI::Data::PendingFluxMapModel::cleared,
            m_flux_map_world_model,
            &SolTrace::GUI::Data::FluxMapWorldModel::on_reset);

    connect(m_pending_flux_maps,
            &SolTrace::GUI::Data::PendingFluxMapModel::failed,
            [this](SolTrace::GUI::Data::Entity entity, QString reason) {
                // we record completion, batch or not
                maybe_update_batch(entity);
                this->notify(ANotification::error(reason));
            });

    connect(m_pending_flux_maps, &SolTrace::GUI::Data::PendingFluxMapModel::cleared, this, [this] {
        set_current_flux_stats({});
    });

    connect(
        m_pending_flux_maps, &SolTrace::GUI::Data::PendingFluxMapModel::all_done, this, [this] {
            m_batch.reset(-1);
            batch_progress(-1, -1);
        });
}

void FluxModule::set_results(SolTrace::GUI::Data::SimulationResultPtr p) {
    m_results = p;
    set_current_entity({});
    set_current_entity_name(QString());
    set_current_entity_position({});
    set_current_flux_stats({});
    set_current_image(QString());
    m_entity_model->reset(nullptr);
    m_computed_maps_model->reset(nullptr);
    m_pending_flux_maps->reset(nullptr);
    m_ray_iso_volume->set_current_mesh({});

    if (!p) return;

    auto mptr = const_cast<SolTrace::GUI::Data::Database*>(p->database.get());

    m_entity_model->reset(mptr);
    m_computed_maps_model->reset(mptr);
    m_pending_flux_maps->reset(p);
    m_ray_iso_volume->set_current_mesh({});

    m_batch.reset();
    emit batch_progress(-1, -1);

    // Default to the entity with the most ray hits.
    entt::entity largest = entt::null;
    size_t       best    = 0;

    for (auto& [c, v] : p->entity_to_ray_ids) {
        if (v.size() > best) {
            largest = c;
            best    = v.size();
        }
    }

    select_entity(largest);
}

void FluxModule::select_entity(SolTrace::GUI::Data::Entity entity) {
    set_current_entity(entity);

    if (!m_results || !m_results->database || !entity.is_valid()) {
        set_current_entity_name(QString());
        set_current_entity_position({});
        set_current_flux_stats({});
        return;
    }

    set_current_entity_name(m_results->database->name_of(entity));

    auto* database = m_results->database.get();
    auto* global   = database->global_transform.get(entity);
    auto  transform =
        global ? *global
               : SolTrace::GUI::Data::GlobalTransformComponent::compute_for(database->as_registry(),
                                                           entity);
    set_current_entity_position(QVector3D(transform.position.x,
                                          transform.position.y,
                                          transform.position.z));

    refresh_current_flux_stats();
}

void FluxModule::maybe_update_batch(SolTrace::GUI::Data::Entity entity) {
    if (!m_batch.mark_done(entity)) return;
    emit this->batch_progress(m_batch.completed(), m_batch.max());
}

void FluxModule::refresh_current_flux_stats() {
    for (auto const& item : m_flux_map_world_model->vector()) {
        if (item.flux_entity == current_entity()) {
            set_current_flux_stats(item.flux_stats);
            set_current_image(item.flux_image_path);
            return;
        }
    }

    set_current_flux_stats({});
    set_current_image(QString());
}

void FluxModule::flux_map_ready(SolTrace::GUI::Data::Entity              entity,
                                SolTrace::GUI::Analysis::BakedFluxMapPtr,
                                SolTrace::GUI::Data::Database const*) {
    // we record completion, batch or not
    maybe_update_batch(entity);

    if (entity != current_entity()) return;

    refresh_current_flux_stats();
}

void FluxModule::start_generate() {
    qCDebug(fluxLog) << Q_FUNC_INFO
                     << "Starting fluxmap generation for current entity";

    if (!m_results) {
        emit notify(ANotification::warning(
            "Run a trace before generating a flux map."));
        return;
    }

    if (!current_entity().is_valid()) {
        emit notify(ANotification::warning(
            "Select an element to generate a flux map."));
        return;
    }

    if (m_batch.is_pending()) {
        emit notify(
            ANotification::warning("Please wait for a batch to complete before "
                                   "generating a flux map."));
        return;
    }

    m_pending_flux_maps->set_dni(dni());
    m_pending_flux_maps->start_generate_for(current_entity());
}

void FluxModule::start_generate_volume_flux(unsigned resolution) {
    if (ray_volume_flux_in_progress()) {
        emit notify(ANotification::info(
            "Volume flux generation is already running."));
        return;
    }

    if (!m_results) {
        emit notify(ANotification::warning(
            "Run a trace before generating volume flux."));
        return;
    }

    set_ray_volume_flux_in_progress(true);

    qCDebug(fluxLog) << Q_FUNC_INFO << "Starting volume flux raster";

    launch_async_task<Support::SparseGrid3D<float>, QString>(
        QUuid::createUuid(),
        this,
        &FluxModule::flux_vol_ready,
        &FluxModule::flux_vol_failed,
        SolTrace::GUI::Analysis::compute_ray_volume_raster,
        resolution,
        m_results);
}

void FluxModule::start_generate_isosurface(float value) {
    if (!m_results) {
        emit notify(ANotification::warning(
            "Run a simulation before creating an isosurface."));
        return;
    }
    if (!m_results->ray_volume.size_in_bricks()) {
        emit notify(ANotification::warning(
            "Generate volume flux before creating an isosurface."));
        return;
    }

    qCDebug(fluxLog) << Q_FUNC_INFO << "launching volume generation" << value;
    launch_async_task<SolTrace::GUI::Data::Mesh, QString>(QUuid::createUuid(),
                                         this,
                                         &FluxModule::iso_surf_ready,
                                         &FluxModule::iso_surf_failed,
                                         SolTrace::GUI::Analysis::volume_to_mesh,
                                         m_results->ray_volume,
                                         value);
}

void FluxModule::save_image(QString requested_image, QUrl path) {
    static constexpr char TO_REMOVE[]    = "image://fluxmap/";
    static constexpr auto TO_REMOTE_SIZE = std::size(TO_REMOVE) - 1;

    requested_image = requested_image.mid(TO_REMOTE_SIZE);

    qCDebug(fluxLog) << Q_FUNC_INFO << requested_image << path;
    if (!m_image_provider) {
        emit notify(ANotification::error(
            "Internal error trying to save image: missing image provider"));
        return;
    }

    QSize img_size;

    auto image =
        m_image_provider->requestImage(requested_image, &img_size, QSize());

    if (image.isNull()) {
        emit notify(ANotification::error("Internal error trying to save image: "
                                         "unable to fetch requested image"));
        return;
    }


    if (!image.save(path.toLocalFile())) {
        emit notify(ANotification::error("Internal error trying to save image: "
                                         "unable to save image to given path"));
        return;
    }
}

void FluxModule::start_generate_batch() {
    qCDebug(fluxLog) << Q_FUNC_INFO << "Starting batch fluxmap generation";

    if (!m_results) {
        emit notify(ANotification::warning(
            "Run a trace before generating a flux map."));
        return;
    }

    if (m_pending_flux_maps->rowCount() != 0) {
        emit notify(ANotification::error(
            "Please wait until all fluxmaps have completed generation before "
            "starting a batch."));
        return;
    }

    m_pending_flux_maps->set_dni(dni());

    m_batch.reset();

    // load em all up

    for (auto const& entity_list : m_results->entity_to_ray_ids) {
        auto entity = entity_list.first;

        if (entity_list.second.empty()) continue;

        bool ok = m_pending_flux_maps->start_generate_for(entity);

        if (ok) { m_batch.add_pending(entity); }
    }

    if (m_batch.empty()) {
        m_batch.reset(-1);
        return;
    }

    m_batch.start_from_pending();

    emit started_batch();

    emit batch_progress(0, m_batch.max());
}

void FluxModule::cancel_batch() {
    qCDebug(fluxLog) << Q_FUNC_INFO << "Cancelling batch...";

    m_pending_flux_maps->cancel_all();
}

void FluxModule::flux_vol_ready(QUuid const&                 id,
                                Support::SparseGrid3D<float> grid) {
    if (m_results) m_results->ray_volume = grid;
    set_ray_volume_flux_in_progress(false);
    emit notify(ANotification::info("Volume flux generation complete."));
    qCDebug(fluxLog) << Q_FUNC_INFO << id;
}
void FluxModule::flux_vol_failed(QUuid const& id, QString reason) {
    qCDebug(fluxLog) << Q_FUNC_INFO << id;
    set_ray_volume_flux_in_progress(false);
    qCritical() << "Unable to generate volume flux" << reason;
    emit notify(ANotification::error(
        QString("Could not generate volume flux: %1").arg(reason)));
}

void FluxModule::iso_surf_ready(QUuid const& id, SolTrace::GUI::Data::Mesh mesh) {
    qCDebug(fluxLog) << Q_FUNC_INFO << id;
    m_ray_iso_volume->set_current_mesh(mesh);
    emit notify(ANotification::info("Isosurface generation complete."));
}
void FluxModule::iso_surf_failed(QUuid const& id, QString reason) {
    qCDebug(fluxLog) << Q_FUNC_INFO << id;
    qCritical() << "Unable to generate isosurface" << reason;
    emit notify(ANotification::error(
        QString("Could not generate the isosurface: %1").arg(reason)));
}

} // namespace SolTrace::GUI::App

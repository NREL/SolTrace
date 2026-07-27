#include "flux_module.h"
#include "analysis/ray_volume_raster.h"
#include "analysis/volume_to_mesh.h"
#include "utilities/asynctask.h"

#include <QQmlEngine>
#include <QUuid>

namespace SolTrace::GUI::App {

FluxModule::FluxModule(QQmlEngine* engine, QObject* parent)
    : QObject(parent),
      m_entity_model(new db::AllElementsModel(this)),
      m_computed_maps_model(new db::AllComputedMapsModel(this)),
      m_pending_flux_maps(new db::PendingFluxMapModel(this)),
      m_flux_map_world_model(new db::FluxMapWorldModel(this)),
      m_ray_iso_volume(new db::QMLMesh()) {

    set_ray_volume_flux_in_progress(false);

    m_ray_iso_volume->setParent(this);

    auto provider = m_pending_flux_maps->make_new_provider();

    engine->addImageProvider("fluxmap", provider);

    connect(m_pending_flux_maps,
            &db::PendingFluxMapModel::ready,
            m_flux_map_world_model,
            &db::FluxMapWorldModel::on_ready);

    connect(m_pending_flux_maps,
            &db::PendingFluxMapModel::ready,
            this,
            &FluxModule::flux_map_ready);

    connect(m_pending_flux_maps,
            &db::PendingFluxMapModel::cleared,
            m_flux_map_world_model,
            &db::FluxMapWorldModel::on_reset);

    connect(m_pending_flux_maps, &db::PendingFluxMapModel::cleared, this, [this] {
        set_current_flux_stats({});
    });
}

void FluxModule::set_results(db::SimulationResultPtr p) {
    m_results = p;
    set_current_entity({});
    set_current_entity_name(QString());
    set_current_flux_stats({});
    m_entity_model->reset(nullptr);
    m_pending_flux_maps->reset(nullptr);
    m_ray_iso_volume->set_current_mesh({});

    if (!p) return;

    auto mptr = const_cast<db::Database*>(p->database.get());

    m_entity_model->reset(mptr);
    m_pending_flux_maps->reset(p);
    m_ray_iso_volume->set_current_mesh({});

    // HACK HACK HACK

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

void FluxModule::select_entity(db::Entity entity) {
    set_current_entity(entity);

    if (!m_results || !m_results->database || !entity.is_valid()) {
        set_current_entity_name(QString());
        set_current_flux_stats({});
        return;
    }

    set_current_entity_name(m_results->database->name_of(entity));
    refresh_current_flux_stats();
}

void FluxModule::refresh_current_flux_stats() {
    for (auto const& item : m_flux_map_world_model->vector()) {
        if (item.flux_entity == current_entity()) {
            set_current_flux_stats(item.flux_stats);
            return;
        }
    }

    set_current_flux_stats({});
}

void FluxModule::flux_map_ready(db::Entity              entity,
                                analysis::BakedFluxMapPtr image,
                                db::Database const*) {
    if (entity != current_entity()) return;

    set_current_flux_stats(image ? image->stats
                                 : analysis::BakedFluxMapStats {});
}

void FluxModule::start_generate() {
    qDebug() << Q_FUNC_INFO << "Starting fluxmap generation for current entity";
    if (!m_results) {
        emit notify(ANotification::warning(
            "Run a simulation before generating a flux map."));
        return;
    }

    if (!current_entity().is_valid()) {
        emit notify(ANotification::warning(
            "Select an element before generating a flux map."));
        return;
    }

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
            "Run a simulation before generating volume flux."));
        return;
    }

    set_ray_volume_flux_in_progress(true);

    qDebug() << Q_FUNC_INFO << "Starting volume flux raster";

    launch_async_task<analysis::SparseGrid3D<float>, QString>(
        QUuid::createUuid(),
        this,
        &FluxModule::flux_vol_ready,
        &FluxModule::flux_vol_failed,
        analysis::compute_ray_volume_raster,
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

    qDebug() << Q_FUNC_INFO << "launching volume generation" << value;
    launch_async_task<db::Mesh, QString>(QUuid::createUuid(),
                                         this,
                                         &FluxModule::iso_surf_ready,
                                         &FluxModule::iso_surf_failed,
                                         analysis::volume_to_mesh,
                                         m_results->ray_volume,
                                         value);
}

void FluxModule::flux_vol_ready(QUuid const&                  id,
                                analysis::SparseGrid3D<float> grid) {
    if (m_results) m_results->ray_volume = grid;
    set_ray_volume_flux_in_progress(false);
    emit notify(ANotification::info("Volume flux generation complete."));
    qDebug() << Q_FUNC_INFO << id;
}
void FluxModule::flux_vol_failed(QUuid const& id, QString reason) {
    qDebug() << Q_FUNC_INFO << id;
    set_ray_volume_flux_in_progress(false);
    qCritical() << "Unable to generate volume flux" << reason;
    emit notify(ANotification::error(
        QString("Could not generate volume flux: %1").arg(reason)));
}

void FluxModule::iso_surf_ready(QUuid const& id, db::Mesh mesh) {
    qDebug() << Q_FUNC_INFO << id;
    m_ray_iso_volume->set_current_mesh(mesh);
    emit notify(ANotification::info("Isosurface generation complete."));
}
void FluxModule::iso_surf_failed(QUuid const& id, QString reason) {
    qDebug() << Q_FUNC_INFO << id;
    qCritical() << "Unable to generate isosurface" << reason;
    emit notify(ANotification::error(
        QString("Could not generate the isosurface: %1").arg(reason)));
}

} // namespace SolTrace::GUI::App

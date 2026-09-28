#pragma once

#include "scene_models/fluxmapworldmodel.h"
#include "scene_models/mesh_qml_bridge.h"
#include "scene_models/element_models.h"
#include "modules/flux_batch_state.h"
#include "support/notification.h"
#include "support/qt_helpers.h"
#include <QObject>
#include <QVector3D>

namespace SolTrace::GUI::App {

/**
 * @class FluxModule
 * @brief Flux analysis module.
 *
 * Provides access to flux distribution results from the simulation.
 * Shares ResultsBackend with Intersections — both modules read from
 * the same simulation result data.
 *
 * QML access pattern: App.flux.results
 */
class FluxModule : public QObject {
    Q_OBJECT

    QPointer<SolTrace::GUI::Data::FluxMapProvider> m_image_provider;

    // TODO: add front or back filtering

    SolTrace::GUI::Data::SimulationResultPtr m_results;

    QOBJECT_READONLY_PROPERTY(SolTrace::GUI::Data::AllElementsModel, entity_model);
    QOBJECT_READONLY_PROPERTY(SolTrace::GUI::Data::AllComputedMapsModel, computed_maps_model);
    QOBJECT_READONLY_PROPERTY(SolTrace::GUI::Data::PendingFluxMapModel, pending_flux_maps);
    QOBJECT_READONLY_PROPERTY(SolTrace::GUI::Data::FluxMapWorldModel, flux_map_world_model);

    Q_READONLY_PROPERTY(bool, ray_volume_flux_in_progress);
    QOBJECT_READONLY_PROPERTY(SolTrace::GUI::Data::QMLMesh, ray_iso_volume);

    Q_WRITABLE_PROPERTY(SolTrace::GUI::Data::Entity, current_entity, { });
    Q_READONLY_PROPERTY(QString, current_entity_name);
    Q_READONLY_PROPERTY(QVector3D, current_entity_position);
    Q_READONLY_PROPERTY(SolTrace::GUI::Analysis::BakedFluxMapStats, current_flux_stats);

    Q_WRITABLE_PROPERTY(bool, show_flux_volume, true);
    Q_WRITABLE_PROPERTY(bool, show_other_geometry, false);
    Q_WRITABLE_PROPERTY(double, dni, 1000.0);

    // Hack
    Q_WRITABLE_PROPERTY(QString, current_image, { });

    FluxBatchState m_batch;

private:
    void refresh_current_flux_stats();

    void maybe_update_batch(SolTrace::GUI::Data::Entity);

private slots:
    void flux_map_ready(Data::Entity,
                        Analysis::BakedFluxMapPtr,
                        Data::Database const*);

    void flux_vol_ready(QUuid const&, Support::SparseGrid3D<float>);
    void flux_vol_failed(QUuid const&, QString);

    void iso_surf_ready(QUuid const&, SolTrace::GUI::Data::Mesh);
    void iso_surf_failed(QUuid const&, QString);

public:
    explicit FluxModule(QQmlEngine*, QObject* parent = nullptr);

public slots:
    /// Set the result set used for all flux computations and scene models.
    void set_results(Data::SimulationResultPtr);

    /// Select the entity whose flux map/statistics are shown in the UI.
    void select_entity(Data::Entity);

    /// Generate a surface flux map for current_entity.
    void start_generate();

    /// Generate a volumetric ray-density raster from the current result.
    void start_generate_volume_flux(unsigned resolution);

    /// Generate an isosurface mesh from the current volumetric raster.
    void start_generate_isosurface(float value);

    /// Save a flux image to disk
    void save_image(QString requested_image, QUrl path);

    /// Start a long running process batch generating everything
    void start_generate_batch();

    void cancel_batch();

signals:
    void notify(Support::ANotification);

    void started_batch();
    void batch_progress(qint64, qint64);
};

} // namespace SolTrace::GUI::App

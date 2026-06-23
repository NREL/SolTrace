#pragma once

#include "backend.h"
#include "database/fluxmapworldmodel.h"
#include "database/mesh_qml_bridge.h"
#include "database/rootelementsmodel.h"
#include "module_common.h"
#include "utilities/notification.h"
#include "utilities/qt_helpers.h"
#include <QObject>

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

    // TODO: add front or back filtering

    db::SimulationResultPtr m_results;

    QOBJECT_READONLY_PROPERTY(db::AllElementsModel, entity_model);
    QOBJECT_READONLY_PROPERTY(db::AllComputedMapsModel, computed_maps_model);
    QOBJECT_READONLY_PROPERTY(db::PendingFluxMapModel, pending_flux_maps);
    QOBJECT_READONLY_PROPERTY(db::FluxMapWorldModel, flux_map_world_model);

    Q_READONLY_PROPERTY(bool, ray_volume_flux_in_progress);
    QOBJECT_READONLY_PROPERTY(db::QMLMesh, ray_iso_volume);

    Q_WRITABLE_PROPERTY(db::Entity, current_entity, {});
    Q_READONLY_PROPERTY(QString, current_entity_name);
    Q_READONLY_PROPERTY(analysis::BakedFluxMapStats, current_flux_stats);

    Q_WRITABLE_PROPERTY(bool, show_flux_volume, true);

    // Hack
    Q_WRITABLE_PROPERTY(QString, current_image, {});

private:
    void refresh_current_flux_stats();

private slots:
    void flux_map_ready(db::Entity, analysis::BakedFluxMapPtr, db::Database const*);

    void flux_vol_ready(QUuid const&, analysis::SparseGrid3D<float>);
    void flux_vol_failed(QUuid const&, QString);

    void iso_surf_ready(QUuid const&, db::Mesh);
    void iso_surf_failed(QUuid const&, QString);

public:
    explicit FluxModule(QQmlEngine*, QObject* parent = nullptr);

public slots:
    void set_results(db::SimulationResultPtr);
    void select_entity(db::Entity);

    void start_generate();

    void start_generate_volume_flux(unsigned resolution);
    void start_generate_isosurface(float value);

signals:
    void notify(ANotification);
};

} // namespace SolTrace::GUI::App

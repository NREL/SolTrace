#pragma once

#include "analysis/flux/baked_flux_map.h"
#include "data/entity.h"
#include "data/simulationresult.h"
#include "support/notification.h"
#include "support/qt_helpers.h"

#include <QObject>
#include <QUrl>

#include <unordered_map>

namespace SolTrace::GUI::App {

/// QML-facing export controller for simulation results and derived flux maps.
///
/// Exports operate on the currently selected SimulationResult. Flux maps are
/// cached as they are generated so export_current() can include optional map
/// images without recomputing them.
class ExportModule : public QObject {
    Q_OBJECT

    SolTrace::GUI::Data::SimulationResultPtr                                   m_results;
    std::unordered_map<SolTrace::GUI::Data::Entity, SolTrace::GUI::Analysis::BakedFluxMapPtr> m_flux_maps;

    Q_WRITABLE_PROPERTY(QUrl, export_directory, { });
    Q_WRITABLE_PROPERTY(bool, export_flux_map_images, true);
    Q_WRITABLE_PROPERTY(bool, export_flux_map_meshes, false);
    Q_WRITABLE_PROPERTY(bool, export_rays, true);
    Q_WRITABLE_PROPERTY(bool, export_scene_copy, false);
    Q_WRITABLE_PROPERTY(bool, random_sample_rays, false);
    Q_WRITABLE_PROPERTY(int, random_sample_ray_count, 10000);

    Q_READONLY_PROPERTY(QString, current_result_name);
    Q_READONLY_PROPERTY(int, generated_flux_map_count);
    Q_READONLY_PROPERTY(bool, can_export);

private:
    void    update_can_export();
    QString current_result_file_stem() const;

public:
    explicit ExportModule(QObject* parent = nullptr);

public slots:
    /// Set the simulation result that subsequent export operations use.
    void set_results(SolTrace::GUI::Data::SimulationResultPtr);

    /// Store a generated flux map so it can be included in exports.
    void
    cache_flux_map(SolTrace::GUI::Data::Entity, SolTrace::GUI::Analysis::BakedFluxMapPtr, SolTrace::GUI::Data::Database const*);

    /// Export the current result according to the configured export options.
    void export_current();

signals:
    void notify(Support::ANotification);
};

} // namespace SolTrace::GUI::App

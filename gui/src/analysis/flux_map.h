#pragma once

#include <QImage>
#include <QObject>

#include "baked_flux_map.h"
#include "database/mesh.h"
#include "job_control/job_run_common.h"
#include "utilities/qt_helpers.h"

namespace analysis {

/// Options for rasterizing ray/surface hits into a flux-map image.
struct FluxMapBakeOptions {
    // TODO: change to support UV aspect ratio
    glm::uvec2 image_resolution = { 1024, 1024 };

    /// The color of grid lines (ie, UV grids). If null, does not draw lines.
    QColor grid_line_color = QColor();

    /// The color map to use. If null, uses a classic blue to red map.
    QImage color_map;

    /// Solar direct normal irradiance in W/m^2.
    double dni = 1000.0;
};

/// Asynchronous flux-map generator for simulation results.
///
/// The input mesh is expected to be a UV-unwrapped receiver surface. Generated
/// maps are piecewise-continuous triangle rasters emitted by signal when ready.
class FluxMapComputer : public QObject {
    Q_OBJECT

    db::SimulationResultPtr m_database;

public:
    /// Create a flux-map computer with no attached result set.
    explicit FluxMapComputer(QObject* parent);
    ~FluxMapComputer() override;

    /// Set the current simulation results
    void set_results(db::SimulationResultPtr);

public slots:
    /// Start generating a flux map for the given entity, its mesh, and options.
    /// The mesh MUST NOT have overlapping UVs.
    bool start_generate_for(db::Entity,
                            db::Mesh                     mesh,
                            analysis::FluxMapBakeOptions options);

signals:
    /// Issued when the computation has completed.
    void image_ready(db::Entity, analysis::BakedFluxMapPtr);

    /// Issued when the computation has failed.
    void image_failed(db::Entity, QString reason);

    /// Issued as the map is being generated, 0-100.
    void image_progress(db::Entity, int);

    /// Issued (or can be issued externally) to cancel all map generation.
    void cancel_all();

    /// Issued (or can be issued externally) to cancel a specific map gen.
    void cancel_specific(db::Entity);
};

} // namespace analysis

#pragma once

#include "database/mesh.h"
#include "utilities/grid2d.h"

#include <QImage>
#include <QMetaType>
#include <QObject>
#include <QtGlobal>
#include <QVector3D>
#include <qqmlintegration.h>

namespace analysis {

/// Summary statistics computed with the baked flux map.
struct BakedFluxMapStats {
    Q_GADGET
    QML_VALUE_TYPE(baked_flux_map_stats);
    Q_PROPERTY(quint64 source_ray_count MEMBER source_ray_count)
    Q_PROPERTY(quint64 plotted_ray_count MEMBER plotted_ray_count)
    Q_PROPERTY(double power_per_ray MEMBER power_per_ray)
    Q_PROPERTY(double plotted_power MEMBER plotted_power)
    Q_PROPERTY(double peak_flux MEMBER peak_flux)
    Q_PROPERTY(double min_flux MEMBER min_flux)
    Q_PROPERTY(double average_flux MEMBER average_flux)
    Q_PROPERTY(double sigma_flux MEMBER sigma_flux)
    Q_PROPERTY(double uniformity MEMBER uniformity)
    Q_PROPERTY(double peak_flux_uncertainty MEMBER peak_flux_uncertainty)
    Q_PROPERTY(double average_flux_uncertainty MEMBER average_flux_uncertainty)
    Q_PROPERTY(QVector3D centroid MEMBER centroid)

public:
    /// Total rays in the simulation result used for this map.
    quint64 source_ray_count = 0;

    /// Number of plotted ray/surface interactions that contributed to this map.
    quint64 plotted_ray_count = 0;

    /// Power represented by each plotted ray interaction, in W when scaled.
    double power_per_ray = 1.0;

    /// Integrated power represented by plotted interactions, in W when scaled.
    double plotted_power = 0.0;

    /// Flux statistics over the rasterized map values.
    double peak_flux    = 0.0;
    double min_flux     = 0.0;
    double average_flux = 0.0;
    double sigma_flux   = 0.0;
    double uniformity   = 0.0;

    /// Monte Carlo uncertainty estimates in percent where counts are known.
    double peak_flux_uncertainty    = 0.0;
    double average_flux_uncertainty = 0.0;

    /// Average world-space location of plotted interactions.
    QVector3D centroid;

    bool operator==(BakedFluxMapStats const&) const = default;
};

/// A computed flux map
struct BakedFluxMap {
    /// Per-pixel flux values used with a color map to generate images.
    Grid2D<float> counts;

    /// Mesh used when projecting ray hits into the map.
    db::Mesh mesh;

    /// Per-face world-space area, ordered with mesh.triangles.
    QVector<float> face_area;

    /// Per-face ray interaction count, ordered with mesh.triangles.
    QVector<quint64> face_ray_count;

    /// A generated image using counts and a color map
    QImage bin_map;

    /// A generated image showing UV points of intersecting rays
    QImage point_map;

    /// Summary statistics for the generated map.
    BakedFluxMapStats stats;
};

using BakedFluxMapPtr = std::shared_ptr<BakedFluxMap const>;

} // namespace analysis

Q_DECLARE_METATYPE(analysis::BakedFluxMapStats)

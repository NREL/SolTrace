#pragma once

#include <QMutex>
#include <QObject>
#include <QQuaternion>
#include <QQuick3DTextureData>
#include <QQuickImageProvider>
#include <QVariantMap>

#include "analysis/flux/flux_map.h"
#include "data/database.h"
#include "editing/geometryeditor.h"
#include "scene_models/entity_name_model.h"
#include "data/simulationresult.h"

namespace SolTrace::GUI::Analysis {
class FluxMapComputer;
}

namespace SolTrace::GUI::Data {

class FluxMapProvider;

/// Quick3D texture data backed by a generated flux-map image.
class FluxTextureData : public QQuick3DTextureData {
    Q_OBJECT

public:
    explicit FluxTextureData(QQuick3DObject* parent = nullptr);

    /// Replace texture contents with image data.
    void set_image(QImage const& image);
};

// ============================================================================

/// One pending flux-map computation exposed to QML.
struct FluxMappedPendingItem {
    Entity entity;
    int    progress = 0;

    RECORD_META(FluxMappedPendingItem,
                SM_EXPOSE_RO(entity),
                SM_EXPOSE_RO(progress));
};

/// Tracks in-progress surface flux-map computations for one result set.
///
/// Also owns the FluxMapComputer and creates the image provider used by QML.
class PendingFluxMapModel
    : public Support::StructModelAdapter<FluxMappedPendingItem> {
    Q_OBJECT
    QPointer<Database const>            m_host;
    QPointer<SolTrace::GUI::Analysis::FluxMapComputer> m_compute;

private slots:
    void on_changed();

    void on_ready(SolTrace::GUI::Data::Entity, SolTrace::GUI::Analysis::BakedFluxMapPtr);
    void on_failed(SolTrace::GUI::Data::Entity, QString);
    void on_progress(SolTrace::GUI::Data::Entity, int);

public:
    // TODO: Make sure we have bins counts and bin areas to export
    // power per ray

    Q_WRITABLE_PROPERTY(int, mesh_resolution_multiply, 1);
    Q_WRITABLE_PROPERTY(QSize, image_resolution, (QSize { 1024, 1024 }));
    Q_WRITABLE_PROPERTY(QColor, mesh_line_color, "black");
    Q_WRITABLE_PROPERTY(bool, show_mesh_grid, false);
    Q_WRITABLE_PROPERTY(QString, color_map, "");
    Q_WRITABLE_PROPERTY(double, dni, 1000.0);

    explicit PendingFluxMapModel(QObject* parent = nullptr);
    virtual ~PendingFluxMapModel() = default;

    /// Reset computation state for a new simulation result.
    void reset(SolTrace::GUI::Data::SimulationResultPtr);

    /// Database associated with the current simulation result, if any.
    Database const* database() { return m_host; }

    /// Create a provider for flux map textures. This MUST be done at the start
    /// of the app, and before any databases!
    FluxMapProvider* make_new_provider();

public slots:
    /// Start generating a flux map for entity.
    bool start_generate_for(SolTrace::GUI::Data::Entity);

    /// Cancel the pending map generation for entity.
    void cancel_for(SolTrace::GUI::Data::Entity);

    void cancel_all();

signals:
    /// A new flux map is ready, for an entity, from a given database
    void ready(SolTrace::GUI::Data::Entity, SolTrace::GUI::Analysis::BakedFluxMapPtr, SolTrace::GUI::Data::Database const*);
    void failed(SolTrace::GUI::Data::Entity, QString);
    void cleared();

    void all_done();
};

// ============================================================================

/// One rendered flux map and its transform in the 3D result scene.
struct FluxMappedItem {
    Entity                               flux_entity;
    std::shared_ptr<QQuick3DTextureData> flux_texture_data;
    QString                              flux_image_path;
    QVector3D                            flux_position;
    QQuaternion                          flux_rotation;
    std::shared_ptr<SurfaceGeometry>     flux_geometry;
    SolTrace::GUI::Analysis::BakedFluxMapStats          flux_stats;

    RECORD_META(FluxMappedItem,
                SM_EXPOSE_RO(flux_entity),
                SM_EXPOSE_RO(flux_texture_data),
                SM_EXPOSE_RO(flux_image_path),
                SM_EXPOSE_RO(flux_position),
                SM_EXPOSE_RO(flux_rotation),
                SM_EXPOSE_RO(flux_geometry),
                SM_EXPOSE_RO(flux_stats));
};

/// Model of completed flux maps shown in the 3D result scene.
class FluxMapWorldModel : public Support::StructModelAdapter<FluxMappedItem> {
    Q_OBJECT

public:
    // TODO: power per ray

    explicit FluxMapWorldModel(QObject* parent = nullptr);
    virtual ~FluxMapWorldModel() = default;

public slots:
    /// Compute current completed flux-map geometry bounds on demand.
    Q_INVOKABLE QVariantMap content_bounds() const;

    /// Clear all completed flux map rows.
    void on_reset();

    /// Add or replace the completed map for an entity.
    void on_ready(Entity, SolTrace::GUI::Analysis::BakedFluxMapPtr, Database const*);
};

// ============================================================================

/// Image provider for QML image://fluxmap requests.
class FluxMapProvider : public QQuickImageProvider {
    Q_OBJECT

    QHash<QString, SolTrace::GUI::Analysis::BakedFluxMapPtr> m_store;
    QMutex                                    m_lock;

public:
    FluxMapProvider();

    /// Return a generated flux-map image for id.
    QImage requestImage(QString const& id,
                        QSize*         size,
                        QSize const&   requestedSize) override;

public slots:
    /// Store a completed flux map under its generated image id.
    void on_ready(Entity, SolTrace::GUI::Analysis::BakedFluxMapPtr, Database const*);

    /// Remove all stored images.
    void clear();
};

// =============================================================================

/// Model listing entities that already have computed flux maps.
class AllComputedMapsModel
    : public Support::StructModelAdapter<EntityNamePair> {
    Q_OBJECT

    QPointer<Database> m_host;

    std::unordered_map<entt::entity, int> m_reverse;

    QVector<EntityNamePair> rebuild_lists();

private slots:
    void recompute();
    void ident_changed(entt::entity);

public:
    explicit AllComputedMapsModel(QObject* parent = nullptr);
    ~AllComputedMapsModel() override = default;

    /// Observe a database and rebuild the computed-map list.
    void reset(Database* database);

public slots:
    /// Return the row for a computed flux map entity, or -1 if absent.
    int index_of(SolTrace::GUI::Data::Entity entity) const;

    /// Return the computed flux map entity at row, or an invalid entity.
    SolTrace::GUI::Data::Entity entity_at(int index) const;
};


} // namespace SolTrace::GUI::Data

#include "fluxmapworldmodel.h"

#include "analysis/flux_map.h"
#include "database/components.h"
#include "database/surface.h"

#include <algorithm>
#include <limits>
#include <QVariantMap>

namespace db {

namespace {

// TODO: We can use something like this all over the place.
struct QVectorBoundsAccumulator {
    bool      valid = false;
    QVector3D min;
    QVector3D max;

    QVectorBoundsAccumulator()
        : min(std::numeric_limits<float>::max(),
              std::numeric_limits<float>::max(),
              std::numeric_limits<float>::max()),
          max(-std::numeric_limits<float>::max(),
              -std::numeric_limits<float>::max(),
              -std::numeric_limits<float>::max()) { }

    void include(QVector3D const& point) {
        min.setX(std::min(min.x(), point.x()));
        min.setY(std::min(min.y(), point.y()));
        min.setZ(std::min(min.z(), point.z()));
        max.setX(std::max(max.x(), point.x()));
        max.setY(std::max(max.y(), point.y()));
        max.setZ(std::max(max.z(), point.z()));
        valid = true;
    }
};

QVariantMap bounds_map(QVectorBoundsAccumulator const& bounds) {
    return QVariantMap {
        { QStringLiteral("valid"), bounds.valid },
        { QStringLiteral("min"), bounds.min },
        { QStringLiteral("max"), bounds.max },
    };
}

void include_transformed_box(QVectorBoundsAccumulator& bounds,
                             BoundingBox const&       box,
                             QVector3D const&         position,
                             QQuaternion const&       rotation) {
    for (int ix = 0; ix < 2; ++ix) {
        for (int iy = 0; iy < 2; ++iy) {
            for (int iz = 0; iz < 2; ++iz) {
                QVector3D corner(
                    ix ? box.max.x() : box.min.x(),
                    iy ? box.max.y() : box.min.y(),
                    iz ? box.max.z() : box.min.z());

                bounds.include(position + rotation.rotatedVector(corner));
            }
        }
    }
}

} // namespace

static QString image_name(Entity e) {
    return QString::number(entt::to_integral((entt::entity)e));
}

static QString const POINT_MAP_SUFFIX = QStringLiteral("_point_map");

static QString point_map_name(QString name) {
    return name + POINT_MAP_SUFFIX;
}

FluxTextureData::FluxTextureData(QQuick3DObject* parent)
    : QQuick3DTextureData(parent) { }

void FluxTextureData::set_image(QImage const& image) {
    auto rgba = image.convertToFormat(QImage::Format_RGBA8888);

    setSize(rgba.size());
    setFormat(QQuick3DTextureData::RGBA8);
    setHasTransparency(rgba.hasAlphaChannel());
    setTextureData(QByteArray(reinterpret_cast<const char*>(rgba.constBits()),
                              rgba.sizeInBytes()));

    update();
}

// ============================================================================

PendingFluxMapModel::PendingFluxMapModel(QObject* parent)
    : StructModelAdapter(parent) {

    m_compute = new analysis::FluxMapComputer(this);

    connect(m_compute,
            &analysis::FluxMapComputer::image_progress,
            this,
            &PendingFluxMapModel::on_progress);

    connect(m_compute,
            &analysis::FluxMapComputer::image_ready,
            this,
            &PendingFluxMapModel::on_ready);

    connect(m_compute,
            &analysis::FluxMapComputer::image_failed,
            this,
            &PendingFluxMapModel::on_failed);
}


void PendingFluxMapModel::reset(db::SimulationResultPtr res) {
    emit cleared();

    m_compute->set_results(res);

    on_changed();

    if (res && res->database) {
        m_host = res->database.get();
        for (auto const& [e, c] :
             res->database->as_registry().view<HasFluxMapComponent>().each()) {
            emit ready(e, c.map_info, res->database.get());
        }
    } else {
        m_host = nullptr;
    }
}

void PendingFluxMapModel::on_changed() {
    store_remove_all();
}

void PendingFluxMapModel::on_ready(Entity e, analysis::BakedFluxMapPtr image) {
    if (!m_host) return;

    {
        // super hazardous, but this is how we make sure modification of a
        // 'frozen' database is kept limited

        auto* ptr = const_cast<db::Database*>(m_host.data());

        ptr->as_registry().emplace_or_replace<HasFluxMapComponent>(
            e, HasFluxMapComponent { .map_info = image });
    }


    store_remove_by_predicate([e](auto& record) { return record.entity == e; });

    emit ready(e, image, m_host);
}

void PendingFluxMapModel::on_failed(Entity e, QString reason) {
    emit failed(reason);

    this->store_remove_by_predicate(
        [e](FluxMappedPendingItem const& item) { return item.entity == e; });
}

void PendingFluxMapModel::on_progress(Entity e, int progress) {
    // find and update

    for (int i = 0; i < this->rowCount(); i++) {
        auto* p = get_at(i);

        if (!p) continue;

        if (p->entity == e) {
            auto copy = *p;

            copy.progress = progress;

            store_push_update(i, copy);

            return;
        }
    }
}

static std::optional<Mesh>
find_mesh_for(db::SurfaceGenerationOptions surface_options,
              Database const*              database,
              Entity                       entity) {
    auto* surface_membership = database->geometry_group_membership.get(entity);

    if (!surface_membership) { return {}; }

    auto* surface =
        database->geometry_parameters.get(surface_membership->group);

    if (!surface) { return {}; }

    auto mesh = db::generate_surface(
        surface->surface, surface->aperture, surface_options);
    if (!mesh) { return {}; }

    auto global = database->global_transform.get(entity);
    auto world =
        global ? *global
               : GlobalTransformComponent::compute_for(database->as_registry(),
                                                      entity);

    for (auto& vertex : mesh->vertex) {
        vertex.position =
            glm::vec3(world.position + world.rotation * glm::dvec3(vertex.position));
        vertex.normal =
            glm::normalize(glm::vec3(world.rotation * glm::dvec3(vertex.normal)));
    }

    return mesh;
}

// void PendingFluxMapModel::start_debug() {
//     // find something with name zero
//     if (!m_host) return;
//     for (auto const& [e, c, g] :
//          m_host->as_registry()
//              .view<IdentityComponent, GeometryGroupMemberComponent>()
//              .each()) {
//         if (c.name == "0") {
//             start_generate_for(e);
//             return;
//         }
//     }
// }

bool PendingFluxMapModel::start_generate_for(Entity entity) {
    if (!m_host) return false;

    qDebug() << Q_FUNC_INFO << entity;

    auto mesh_res = std::max(1, mesh_resolution_multiply());

    db::SurfaceGenerationOptions surface_options;
    surface_options.subdivide_polygon_apertures = true;
    surface_options.height_field_resolution *= mesh_res;
    surface_options.radial_subdivisions *= mesh_res;
    surface_options.perimeter_subdivisions *= mesh_res;
    surface_options.cylinder_angular_subdivisions *= mesh_res;
    surface_options.cylinder_length_subdivisions *= mesh_res;

    auto mesh = find_mesh_for(surface_options, m_host, entity);

    if (!mesh) {
        qDebug() << Q_FUNC_INFO << "entity has no mesh, bailing";
        return false;
    }

    qDebug() << Q_FUNC_INFO << mesh->vertex.size() << mesh->triangles.size();

    store_remove_by_predicate([this, entity](auto const& item) {
        if (item.entity == entity) { this->cancel_for(entity); }
        return item.entity == entity;
    });

    auto opts = analysis::FluxMapBakeOptions {
        .image_resolution = { this->image_resolution().width(),
                              this->image_resolution().height(), },
        .grid_line_color =
            this->show_mesh_grid() ? this->mesh_line_color() : QColor(),
        .color_map = QImage(color_map()),
        .dni       = this->dni(),
    };

    store_push_append(FluxMappedPendingItem {
        .entity   = entity,
        .progress = 0,
    });

    if (!m_compute->start_generate_for(entity, *mesh, opts)) {
        store_remove_by_predicate([entity](auto const& item) {
            return item.entity == entity;
        });
        qDebug() << Q_FUNC_INFO << "generation kickoff failed";
        return false;
    }

    qDebug() << Q_FUNC_INFO << "generation kickoff success";

    return true;
}

void PendingFluxMapModel::cancel_for(Entity entity) {
    m_compute->cancel_specific(entity);
}


FluxMapProvider* PendingFluxMapModel::make_new_provider() {
    auto ret = new FluxMapProvider();

    connect(this, &PendingFluxMapModel::ready, ret, &FluxMapProvider::on_ready);
    connect(this, &PendingFluxMapModel::cleared, ret, &FluxMapProvider::clear);

    return ret;
}


// ============================================================================

void FluxMapWorldModel::on_reset() {
    store_remove_all();
}

FluxMapWorldModel::FluxMapWorldModel(QObject* parent)
    : StructModelAdapter(parent) { }

QVariantMap FluxMapWorldModel::content_bounds() const {
    QVectorBoundsAccumulator bounds;

    for (auto const& item : m_records) {
        if (!item.flux_geometry || item.flux_geometry->vertex_count() == 0) {
            continue;
        }

        include_transformed_box(bounds,
                                item.flux_geometry->bounding_box(),
                                item.flux_position,
                                item.flux_rotation);
    }

    return bounds_map(bounds);
}


void FluxMapWorldModel::on_ready(Entity                    e,
                                 analysis::BakedFluxMapPtr img,
                                 Database const*           db) {
    // Avoid duplicating map entries; this list should stay small.

    if (!db) return;

    qDebug() << Q_FUNC_INFO << e << img->bin_map << db;

    for (auto const& item : m_records) {
        if (item.flux_entity == e && item.flux_texture_data) {
            item.flux_texture_data->deleteLater();
        }
    }

    store_remove_by_predicate([e](auto const& item) { return item.flux_entity == e; });

    // get geom params

    auto group = db->geometry_of(e);

    auto geom = std::make_shared<SurfaceGeometry>();

    geom->set(db, group);

    auto global = GlobalTransformComponent::compute_for(db->as_registry(), e);

    auto position =
        QVector3D(global.position.x, global.position.y, global.position.z);

    auto rotation = QQuaternion(global.rotation.w,
                                global.rotation.x,
                                global.rotation.y,
                                global.rotation.z);

    auto texture_data = std::make_shared<FluxTextureData>();
    texture_data->set_image(img->bin_map);

    store_push_append(FluxMappedItem {
        .flux_entity       = e,
        .flux_texture_data = texture_data,
        .flux_image_path   = "image://fluxmap/" + image_name(e),
        .flux_position     = position,
        .flux_rotation     = rotation,
        .flux_geometry     = geom,
        .flux_stats        = img->stats,
    });
}

// ============================================================================

FluxMapProvider::FluxMapProvider()
    : QQuickImageProvider(QQuickImageProvider::ImageType::Image) { }

QImage FluxMapProvider::requestImage(QString const& id,
                                     QSize*         size,
                                     QSize const&   requestedSize) {
    qDebug() << Q_FUNC_INFO << id << requestedSize;

    QString local_id = id;

    QImage ret;

    bool needs_point_map = local_id.endsWith(POINT_MAP_SUFFIX);

    if (needs_point_map) {
        local_id = local_id.left(local_id.size() - POINT_MAP_SUFFIX.size());
    }


    {
        m_lock.lock();
        auto iter = m_store.find(local_id);
        if (iter != m_store.end()) {

            if (needs_point_map) {
                ret = iter.value()->point_map;
            } else {
                ret = iter.value()->bin_map;
            }
        }
        m_lock.unlock();
    }

    if (size) *size = ret.size();

    qDebug() << Q_FUNC_INFO << ret;

    return ret;
}

void FluxMapProvider::on_ready(Entity                    k,
                               analysis::BakedFluxMapPtr v,
                               Database const*) {
    auto name = image_name(k);
    qDebug() << Q_FUNC_INFO << k << v->bin_map << name;
    m_lock.lock();
    m_store[name] = v;
    m_lock.unlock();
}

void FluxMapProvider::clear() {
    m_lock.lock();
    m_store.clear();
    m_lock.unlock();
}

// ============================================================================


QVector<EntityNamePair> AllComputedMapsModel::rebuild_lists() {
    QVector<EntityNamePair> new_recs;
    m_reverse.clear();

    if (!m_host) return {};

    auto view = m_host->as_registry().view<HasFluxMapComponent>();

    for (auto [entity, map] : view.each()) {
        new_recs.push_back(EntityNamePair::record_for_entity(*m_host, entity));
    }

    for (qsizetype i = 0; i < new_recs.size(); ++i) {
        m_reverse[new_recs[i].entity] = static_cast<int>(i);
    }

    return new_recs;
}

void AllComputedMapsModel::recompute() {
    store_reset(rebuild_lists());
}

void AllComputedMapsModel::ident_changed(entt::entity entity) {
    if (!m_host) return;

    if (auto iter = m_reverse.find(entity); iter != m_reverse.end()) {
        store_push_update(iter->second,
                          EntityNamePair::record_for_entity(*m_host, entity));
    }
}


AllComputedMapsModel::AllComputedMapsModel(QObject* parent)
    : StructModelAdapter(parent) { }

void AllComputedMapsModel::reset(Database* database) {
    if (m_host) {
        disconnect(m_host->identity.self(), nullptr, this, nullptr);
        disconnect(m_host->flux_map.self(), nullptr, this, nullptr);
    }

    m_host = database;
    recompute();

    if (!database) return;

    connect(database->identity.self(),
            &ComponentAPIBase::changed,
            this,
            &AllComputedMapsModel::ident_changed);

    connect(database->flux_map.self(),
            &ComponentAPIBase::changed,
            this,
            &AllComputedMapsModel::recompute);

    connect(database->flux_map.self(),
            &ComponentAPIBase::removed,
            this,
            &AllComputedMapsModel::recompute);
}

int AllComputedMapsModel::index_of(db::Entity entity) const {
    if (auto iter = m_reverse.find(entity); iter != m_reverse.end()) {
        return iter->second;
    }

    return -1;
}

db::Entity AllComputedMapsModel::entity_at(int index) const {
    auto const& items = vector();
    if (index < 0 || index >= items.size()) { return {}; }

    return items[index].entity;
}

} // namespace db

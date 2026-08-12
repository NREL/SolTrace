#include "ray_geometry.h"

#include "analysis/ray_volume_raster.h"
#include "utilities/math_utility.h"

#include <QtMath>
#include <algorithm>
#include <cmath>
#include <magic_enum/magic_enum.hpp>

#include <QtConcurrentMap>

#define GLM_ENABLE_EXPERIMENTAL 1

#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>

namespace analysis {

namespace {

/// Number of rays to show by default
constexpr quint64 DEFAULT_VISIBLE_RAY_COUNT = 10000;

/// Percent of a ray count. Used in a few places.
float percent_for_ray_count(quint64 count, quint64 available) {
    if (available == 0) return 0.0f;
    return static_cast<float>(count * 100.0 / available);
}

/// From percent of a count, to the actual count
size_t visible_ray_limit(size_t available, float show_percent) {
    const auto effective_percent = std::clamp(show_percent, 0.0f, 100.0f);
    const auto requested_rays =
        static_cast<double>(available) * effective_percent / 100.0;
    return std::min(available,
                    static_cast<size_t>(std::llround(requested_rays)));
}

/// Ask if our event filter contains an event.
bool includes_event(EventTypeContainer const& filter, db::RayEventType event) {
    return filter.events.contains(event);
}

} // namespace

EventTypeContainer::EventTypeContainer(
    std::initializer_list<db::RayEventType> l)
    : events(l) { }

EventTypeContainer::EventTypeContainer(QStringList l) {
    for (auto const& item : l) {
        auto str = item.toUpper().toStdString();

        auto maybe_enum = magic_enum::enum_cast<db::RayEventType>(str);

        if (!maybe_enum) {
            qDebug() << "Unknown enum name" << item;
            continue;
        }

        events.insert(*maybe_enum);
    }
}

QStringList EventTypeContainer::to_list() const {
    QStringList ret;

    for (auto item : events) {
        ret << QString(magic_enum::enum_name(item).data()).toLower();
    }

    return ret;
}

struct LineVertex {
    QVector3D position;
    QVector2D uv;
};

void RayGeometry::rebuild_geometry() {
    // Kickoff
    qDebug() << Q_FUNC_INFO << "Start";
    clear();

    // No db, no geom for you
    if (!m_database) {
        qDebug() << Q_FUNC_INFO << "No database";
        update(); // note that we have cleared
        return;
    }

    const size_t ray_limit =
        visible_ray_limit(m_database->records.size(), this->show_percent());

    // Collect total vertex count
    size_t vertex_count = 0;
    size_t counted_rays = 0;

    for (auto const& rec : m_database->records) {
        if (counted_rays >= ray_limit) break;
        vertex_count += rec.events.size();
        counted_rays += 1;
    }

    qDebug() << Q_FUNC_INFO << vertex_count;

    std::vector<LineVertex> verts;
    std::vector<uint32_t>   index;
    verts.reserve(vertex_count);
    index.reserve(vertex_count * 2); // close enough

    {
        size_t ray_number     = 0;
        size_t rays_remaining = ray_limit;

        auto target_span = std::span(m_database->records);

        if (m_selected_ray_id >= 0) {
            target_span = target_span.subspan(m_selected_ray_id, 1);
        }

        for (auto const& ray : target_span) {

            if (rays_remaining == 0) { break; }

            rays_remaining -= 1;
            ray_number += 1;


            // Since events are filtered, compute UV ranges from visible events.
            size_t visible_event_count = 0;
            double total_ray_distance  = 0.0;

            QVector3D last_point = { };

            // first compute an idea of the total ray distance
            for (auto const& interaction : ray.events) {

                if (!includes_event(m_include_events, interaction.event))
                    continue;

                auto p = convert(interaction.location);

                // if this is not the first
                if (visible_event_count > 0) {
                    // record delta
                    total_ray_distance += (p - last_point).length();
                }

                last_point = p;
                visible_event_count += 1;
            }

            if (total_ray_distance == 0.0) { total_ray_distance = 1.0; }

            // Reset counter
            size_t visible_point_index  = 0;
            double current_ray_distance = 0.0;
            auto   texture_mode         = this->texture_mode();
            bool   flip_flip            = false;

            for (auto const& interaction : ray.events) {

                if (!includes_event(m_include_events, interaction.event))
                    continue;

                auto p = convert(interaction.location);

                // Segments are flip flopped for simplicity

                if (texture_mode == TextureMode::Segment &&
                    visible_point_index > 0) {

                    float at = float(flip_flip) * .5 + .25;

                    auto prev = static_cast<uint32_t>(verts.size());
                    verts.push_back({
                        .position = last_point,
                        .uv       = { at, 0.0f },
                    });
                    auto cur = static_cast<uint32_t>(verts.size());
                    verts.push_back({
                        .position = p,
                        .uv       = { at, 0.0f },
                    });

                    index.push_back(prev);
                    index.push_back(cur);

                    flip_flip = !flip_flip;

                } else if (texture_mode == TextureMode::Length) {
                    // if this is not the first
                    if (visible_point_index > 0) {
                        // record delta
                        current_ray_distance += (p - last_point).length();
                    }

                    // bool is_active =
                    //     m_selected_ray_id < 0 || ray.id == m_selected_ray_id;

                    QVector2D uv {
                        static_cast<float>(current_ray_distance /
                                           total_ray_distance),
                        0.0,
                    };

                    verts.push_back({
                        .position = p,
                        .uv       = uv,
                    });


                    // install index: connect consecutive vertices within this
                    // ray
                    if (visible_point_index > 0) {
                        auto cur  = static_cast<uint32_t>(verts.size() - 1);
                        auto prev = static_cast<uint32_t>(verts.size() - 2);
                        index.push_back(prev);
                        index.push_back(cur);
                    }
                }

                // Keep this after vertex emission. Segment mode needs the
                // previous visible point to build a flat-colored segment.
                last_point = p;
                visible_point_index += 1;
            }
        }
    }

    qDebug() << Q_FUNC_INFO << "New buffers ready";

    // install all computed geometry

    auto vertex_buffer = QByteArray(reinterpret_cast<const char*>(verts.data()),
                                    verts.size() * sizeof(LineVertex));
    auto index_buffer  = QByteArray(reinterpret_cast<const char*>(index.data()),
                                    index.size() * sizeof(uint32_t));

    addAttribute(QQuick3DGeometry::Attribute::PositionSemantic,
                 0,
                 QQuick3DGeometry::Attribute::ComponentType::F32Type);

    addAttribute(QQuick3DGeometry::Attribute::TexCoord0Semantic,
                 3 * sizeof(float),
                 QQuick3DGeometry::Attribute::ComponentType::F32Type);

    addAttribute(QQuick3DGeometry::Attribute::IndexSemantic,
                 0,
                 QQuick3DGeometry::Attribute::ComponentType::U32Type);

    setStride(sizeof(LineVertex));
    setVertexData(vertex_buffer);
    setIndexData(index_buffer);
    setBounds(QVector3D(m_database->bounds_min.x,
                        m_database->bounds_min.y,
                        m_database->bounds_min.z),
              QVector3D(m_database->bounds_max.x,
                        m_database->bounds_max.y,
                        m_database->bounds_max.z));
    setPrimitiveType(QQuick3DGeometry::PrimitiveType::Lines);

    qDebug() << Q_FUNC_INFO << "Update";
    update();
}

void RayGeometry::inclusion_list_update() {
    qDebug() << Q_FUNC_INFO << "List changed";
    m_include_events = EventTypeContainer(event_include());

    rebuild_geometry();
}

RayGeometry::RayGeometry(QQuick3DObject* parent) : QQuick3DGeometry(parent) {

    m_include_events = EventTypeContainer({
        db::RayEventType::ABSORB,
        db::RayEventType::REFLECT,
        db::RayEventType::TRANSMIT,
    });

    set_event_include(m_include_events.to_list());

    connect(this,
            &RayGeometry::event_include_changed,
            this,
            &RayGeometry::inclusion_list_update);

    connect(this,
            &RayGeometry::show_percent_changed,
            this,
            &RayGeometry::rebuild_geometry);

    connect(this,
            &RayGeometry::selected_ray_id_changed,
            this,
            &RayGeometry::rebuild_geometry);

    connect(this,
            &RayGeometry::texture_mode_changed,
            this,
            &RayGeometry::rebuild_geometry);
}

void RayGeometry::set_results(db::SimulationResultPtr data) {
    qDebug() << Q_FUNC_INFO << "New ray geometry database";
    m_database = std::move(data);
    const auto available =
        m_database ? static_cast<quint64>(m_database->records.size()) : 0;
    set_available_rays(available);

    const auto default_percent = percent_for_ray_count(
        std::min(DEFAULT_VISIBLE_RAY_COUNT, available), available);
    const bool percent_changed = show_percent() != default_percent;
    set_show_percent(default_percent);
    if (!percent_changed) rebuild_geometry();
}

/// Classic; give me the closest point on a line segment to a ray
/// A and B are line points, P is the query ray start, rayDir is the query ray
/// direction.
///
/// Output is the distance, and the closest point is the closestOnSegment out
/// param
static float dist_segment_ray_closest_points(glm::vec3  A,
                                             glm::vec3  B,
                                             glm::vec3  P,
                                             glm::vec3  ray_dir,
                                             glm::vec3& closest_on_segment) {
    glm::vec3 u = B - A;
    glm::vec3 v = ray_dir;
    glm::vec3 w = A - P;

    float a = dot(u, u);
    float b = dot(u, v);
    float c = dot(v, v);
    float d = dot(u, w);
    float e = dot(v, w);

    const float EPS = 1e-8;

    // Degenerate segment: A == B
    if (a < EPS) {
        closest_on_segment = A;

        if (c < EPS) {
            // Degenerate ray too: ray is just point C
            return length(closest_on_segment - P);
        }

        float t            = glm::max(dot(A - P, v) / c, 0.0f);
        auto  closestOnRay = P + t * v;
        return length(closest_on_segment - closestOnRay);
    }

    // Degenerate ray direction: ray is just point C
    if (c < EPS) {
        auto closestOnRay = P;

        float s          = glm::clamp(dot(P - A, u) / a, 0.0f, 1.0f);
        closest_on_segment = A + s * u;
        return length(closest_on_segment - closestOnRay);
    }

    float denom = a * c - b * b;

    float s;
    float t;

    if (denom > EPS) {
        // Closest points on the infinite supporting lines
        s = (b * e - c * d) / denom;
        t = (a * e - b * d) / denom;
    } else {
        // Nearly parallel
        s = 0.0;
        t = e / c;
    }

    // Enforce segment and ray constraints
    s = glm::clamp(s, 0.0f, 1.0f);
    t = glm::max(t, 0.0f);

    // Recompute after clamping to handle endpoint/ray-origin cases
    s = glm::clamp((b * t - d) / a, 0.0f, 1.0f);
    t = glm::max((b * s + e) / c, 0.0f);

    closest_on_segment = A + s * u;
    auto closestOnRay = P + t * v;

    return length(closest_on_segment - closestOnRay);
}

/// Ray cast query result
struct RayCastRayResult {
    /// Ray id, -1 if not found
    int64_t    ray_id = -1;
    /// World position, valid only if valid ray
    glm::dvec3 world_pos;
};

/// For a traced ray, and a query ray, check if this 'intersects' with some
/// tolerance
static RayCastRayResult check_distance(db::RayRecord const& record,
                                       glm::dvec3 const&    world_position,
                                       glm::dvec3 const&    world_direction,
                                       EventTypeContainer const& filter,
                                       float angle_tolerance_rads_cos) {
    if (record.events.empty()) return RayCastRayResult { };


    bool       have_segment_start = false;
    glm::dvec3 segment_a;

    // we assemble ray segments from non-ignored events

    for (auto const& event : record.events) {
        // Skip filtered events
        if (!includes_event(filter, event.event)) continue;

        if (!have_segment_start) {
            segment_a          = event.location;
            have_segment_start = true;
            continue;
        }

        glm::dvec3 segment_b = event.location;

        glm::vec3 closest_segment_point;

        dist_segment_ray_closest_points(segment_a,
                                        segment_b,
                                        world_position,
                                        world_direction,
                                        closest_segment_point);

        auto angle = glm::dot(
            glm::normalize(glm::dvec3(closest_segment_point) - world_position),
            world_direction);

        if (angle > angle_tolerance_rads_cos) {
            return RayCastRayResult {
                .ray_id    = static_cast<int64_t>(record.id),
                .world_pos = closest_segment_point,
            };
        }

        segment_a = segment_b;
    }

    return RayCastRayResult { };
}

void RayGeometry::pick_ray(QVector3D world_position,
                           QVector3D world_direction,
                           float     angle_tolerance_rads) {

    // Sanitize tolerance
    angle_tolerance_rads = std::clamp<float>(angle_tolerance_rads, 0, M_PI);
    float angle_tolerance_rads_cos = std::cos(angle_tolerance_rads);

    // No database, bail
    if (!m_database) return;

    // Sanitize query ray

    auto glm_world_pos =
        glm::dvec3(world_position.x(), world_position.y(), world_position.z());

    auto glm_world_dir = glm::dvec3(
        world_direction.x(), world_direction.y(), world_direction.z());

    if (glm::length2(glm_world_dir) == 0.0) return;

    glm_world_dir = glm::normalize(glm_world_dir);

    // Get filter

    const auto event_filter = m_include_events;

    // Build iteration over all rays, bounded by what we can see

    auto start_iter = m_database->records.begin();
    auto end_iter =
        m_database->records.begin() +
        visible_ray_limit(m_database->records.size(), this->show_percent());

    // Do this in parallel
    RayCastRayResult result = QtConcurrent::blockingMappedReduced(
        start_iter,
        end_iter,
        [=](db::RayRecord const& record) {
            return check_distance(record,
                                  glm_world_pos,
                                  glm_world_dir,
                                  event_filter,
                                  angle_tolerance_rads_cos);
        },
        [=](RayCastRayResult& dest, RayCastRayResult const& next) {
            if (next.ray_id < 0) return;

            if (dest.ray_id < 0) {
                dest = next;
                return;
            }

            auto curr_dist = glm::distance2(glm_world_pos, dest.world_pos);
            auto next_dist = glm::distance2(glm_world_pos, next.world_pos);

            if (next_dist < curr_dist) { dest = next; }
        });

    if (result.ray_id < 0) {
        qDebug() << Q_FUNC_INFO << "No hit";
        return;
    }

    qDebug() << Q_FUNC_INFO << "Hit" << result.ray_id;

    set_selected_ray_id(result.ray_id);
}

} // namespace analysis

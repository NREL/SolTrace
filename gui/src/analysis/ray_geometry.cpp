#include "ray_geometry.h"

#include "analysis/ray_volume_raster.h"
#include "utilities/math_utility.h"

#include <QtMath>
#include <algorithm>
#include <cmath>
#include <magic_enum/magic_enum.hpp>
#include <optional>

#include <QtConcurrentMap>

#define GLM_ENABLE_EXPERIMENTAL 1

#include <glm/glm.hpp>
#include <glm/gtx/intersect.hpp>
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

/// Ray segment clipping
enum class RaySphereTestResult {
    /// Ray segment is fully inside the sphere
    Ok,
    /// Ray segment should be clipped
    Clip,
    /// Ray segment is fully outside the sphere
    Skip
};


/// Snap a segment at a sphere radius
static bool clip_segment(glm::dvec3& p0, glm::dvec3& p1, double radius) {
    // We are just going to use a model of a sphere to make life easy

    auto d = p1 - p0;

    // Quadratic equation coefficients: a*t^2 + b*t + c = 0
    auto a = glm::dot(d, d);

    // Handle edge case where p0 and p1 are the exact same point
    if (a == 0.0) { return glm::dot(p0, p0) <= (radius * radius); }

    auto b = 2.0 * glm::dot(p0, d);

    auto c = glm::dot(p0, p0) - (radius * radius);

    // Calculate discriminant
    auto discriminant = (b * b) - (4.0 * a * c);

    // If negative, the infinite line completely misses the sphere
    if (discriminant < 0.0) { return false; }

    // Find the entry and exit t along the line
    auto disc_root = glm::sqrt(discriminant);
    auto t1        = (-b - disc_root) / (2.0 * a); // Entry point
    auto t2        = (-b + disc_root) / (2.0 * a); // Exit point

    // Find the overlap interval between the sphere [t1, t2] and the segment
    // [0.0, 1.0]
    auto clip_start = std::max(0.0, t1);
    auto clip_end   = std::min(1.0, t2);

    // If the segment doesn't overlap with spheres interior volume
    if (clip_start > clip_end) {
        // Bail
        return false;
    }

    // Modify the original coordinates
    auto original_p0 = p0;

    p0 = original_p0 + clip_start * d;
    p1 = original_p0 + clip_end * d;

    return true;
}

/// Test a line segment against a zero-center sphere, and clip it if it
/// intersects
static RaySphereTestResult
test_ray_sphere(QVector3D& fa, QVector3D& fb, double radius) {
    auto a = convert(fa);
    auto b = convert(fb);

    double rad_squared = radius * radius;

    auto center = glm::dvec3(0);

    auto a_ok = glm::distance2(a, center) < rad_squared;
    auto b_ok = glm::distance2(b, center) < rad_squared;

    if (a_ok and b_ok) {
        // ray is fully contained, bail
        return RaySphereTestResult::Ok;
    }

    if (!a_ok and !b_ok) {
        // ray is fully outside, skip
        return RaySphereTestResult::Skip;
    }

    // straddles, clip

    bool isect = clip_segment(a, b, radius);

    if (!isect) {
        // hmm, ray does not isect but is not fully in or out. skip for now
        return RaySphereTestResult::Skip;
    }

    fa = convert(a);
    fb = convert(b);

    return RaySphereTestResult::Clip;
}

/// Ask if our event filter contains an event.
bool includes_event(EventTypeContainer const& filter, db::RayEventType event) {
    return filter.events.contains(event);
}

bool includes_entity(db::RayRecord const& path, db::Entity entity) {
    if (!entity.is_valid()) return true;

    return std::ranges::any_of(path.events, [entity](auto const& event) {
        return event.entity == entity.value;
    });
}

bool point_inside_sphere(QVector3D const& point, double radius) {
    auto const p = convert(point);
    return glm::distance2(p, glm::dvec3(0)) < radius * radius;
}

struct ClippedSegment {
    QVector3D start;
    QVector3D end;
};

std::optional<ClippedSegment>
visible_segment(QVector3D const& start, QVector3D const& end, double radius) {
    auto clipped_start = start;
    auto clipped_end   = end;

    switch (test_ray_sphere(clipped_start, clipped_end, radius)) {
    case RaySphereTestResult::Ok:
    case RaySphereTestResult::Clip:
        return ClippedSegment {
            .start = clipped_start,
            .end   = clipped_end,
        };
    case RaySphereTestResult::Skip: return std::nullopt;
    }

    return std::nullopt;
}


double extract_total_ray_distance(db::RayRecord const&      path,
                                  EventTypeContainer const& include_events,
                                  double                    filter_sphere) {
    QVector3D last_point;
    bool      have_last_point    = false;
    double    total_ray_distance = 0.0;

    for (auto const& interaction : path.events) {

        if (!includes_event(include_events, interaction.event)) continue;

        auto p = convert(interaction.location);

        if (!have_last_point) {
            last_point      = p;
            have_last_point = true;
            continue;
        }

        auto clipped = visible_segment(last_point, p, filter_sphere);

        if (clipped) {
            total_ray_distance +=
                (clipped->end - clipped->start).length();
        }

        last_point = p;
    }


    return total_ray_distance;
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

QVariantMap RayGeometry::content_bounds() const {
    return QVariantMap {
        { QStringLiteral("valid"), m_content_bounds_valid },
        { QStringLiteral("min"),
          m_content_bounds_valid ? m_content_bounds_min : QVector3D() },
        { QStringLiteral("max"),
          m_content_bounds_valid ? m_content_bounds_max : QVector3D() },
    };
}

void RayGeometry::rebuild_geometry() {
    // Kickoff
    qDebug() << Q_FUNC_INFO << "Start";
    clear();
    m_content_bounds_valid = false;
    m_content_bounds_min   = QVector3D();
    m_content_bounds_max   = QVector3D();

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

    // set up max volume
    double filter_sphere = max_ray_distance();
    auto   texture_mode  = this->texture_mode();
    auto   render_mode   = this->isect_mode();
    bool   is_point_mode = render_mode == IntersectionMode::Point;


    std::vector<LineVertex> verts;
    std::vector<uint32_t>   index;
    verts.reserve(vertex_count);
    index.reserve(vertex_count * 2); // close enough

    auto append_vertex = [this, &verts](LineVertex vertex) {
        if (!m_content_bounds_valid) {
            m_content_bounds_min   = vertex.position;
            m_content_bounds_max   = vertex.position;
            m_content_bounds_valid = true;
        } else {
            m_content_bounds_min.setX(
                std::min(m_content_bounds_min.x(), vertex.position.x()));
            m_content_bounds_min.setY(
                std::min(m_content_bounds_min.y(), vertex.position.y()));
            m_content_bounds_min.setZ(
                std::min(m_content_bounds_min.z(), vertex.position.z()));
            m_content_bounds_max.setX(
                std::max(m_content_bounds_max.x(), vertex.position.x()));
            m_content_bounds_max.setY(
                std::max(m_content_bounds_max.y(), vertex.position.y()));
            m_content_bounds_max.setZ(
                std::max(m_content_bounds_max.z(), vertex.position.z()));
        }

        verts.push_back(vertex);
    };

    {
        size_t ray_number     = 0;
        size_t rays_remaining = ray_limit;

        auto target_span = std::span(m_database->records);

        if (m_selected_ray_id >= 0) {
            target_span = target_span.subspan(m_selected_ray_id, 1);
        }

        for (auto const& path : target_span) {

            if (rays_remaining == 0) { break; }

            rays_remaining -= 1;
            ray_number += 1;

            if (!includes_entity(path, entity_filter())) continue;

            // first compute an idea of the total ray distance
            double total_ray_distance = extract_total_ray_distance(
                path, m_include_events, filter_sphere);

            if (total_ray_distance == 0.0) {
                // No drawable non-zero-length segments for this ray.
                if (!is_point_mode) continue;

                // we set this here to avoid div by zeros
                total_ray_distance = 1.0;
            }


            QVector3D last_point;
            bool      have_last_point = false;
            // Since events are filtered, compute UV ranges from visible events.
            size_t visible_point_index  = 0;
            double current_ray_distance = 0.0;

            bool flip_flip = false;


            for (auto const& interaction : path.events) {

                if (!includes_event(m_include_events, interaction.event))
                    continue;

                auto p = convert(interaction.location);

                if (is_point_mode) {
                    if (point_inside_sphere(p, filter_sphere)) {
                        auto at = 0.0f;

                        switch (texture_mode) {
                        case TextureMode::SolidColor: break;
                        case TextureMode::Length:
                            at = static_cast<float>(visible_point_index) /
                                 static_cast<float>(
                                     std::max<size_t>(1, path.events.size()));
                            break;
                        case TextureMode::Segment:
                            at        = float(flip_flip) * .5f + .25f;
                            flip_flip = !flip_flip;
                            break;
                        }

                        append_vertex({
                            .position = p,
                            .uv       = { at, 0.0f },
                        });

                        visible_point_index += 1;
                    }

                    last_point      = p;
                    have_last_point = true;
                    continue;
                }

                if (!have_last_point) {
                    last_point      = p;
                    have_last_point = true;
                    continue;
                }

                auto clipped = visible_segment(last_point, p, filter_sphere);

                if (!clipped) {
                    last_point = p;
                    continue;
                }

                auto segment_distance =
                    (clipped->end - clipped->start).length();

                if (segment_distance == 0.0) {
                    last_point = p;
                    continue;
                }

                // Segments are flip flopped for simplicity

                switch (texture_mode) {
                case TextureMode::SolidColor: {
                    static QVector2D uv {
                        0.0,
                        0.0,
                    };


                    auto prev = static_cast<uint32_t>(verts.size());
                    append_vertex({
                        .position = clipped->start,
                        .uv       = uv,
                    });

                    auto cur = static_cast<uint32_t>(verts.size());
                    append_vertex({
                        .position = clipped->end,
                        .uv       = uv,
                    });

                    index.push_back(prev);
                    index.push_back(cur);
                    break;
                }
                case TextureMode::Length: {
                    // bool is_active =
                    //     m_selected_ray_id < 0 || ray.id == m_selected_ray_id;

                    QVector2D start_uv {
                        static_cast<float>(current_ray_distance /
                                           total_ray_distance),
                        0.0,
                    };

                    current_ray_distance += segment_distance;

                    QVector2D end_uv {
                        static_cast<float>(current_ray_distance /
                                           total_ray_distance),
                        0.0,
                    };

                    auto prev = static_cast<uint32_t>(verts.size());
                    append_vertex({
                        .position = clipped->start,
                        .uv       = start_uv,
                    });

                    auto cur = static_cast<uint32_t>(verts.size());
                    append_vertex({
                        .position = clipped->end,
                        .uv       = end_uv,
                    });

                    index.push_back(prev);
                    index.push_back(cur);
                    break;
                }
                case TextureMode::Segment: {
                    float at = float(flip_flip) * .5 + .25;

                    auto prev = static_cast<uint32_t>(verts.size());
                    append_vertex({
                        .position = clipped->start,
                        .uv       = { at, 0.0f },
                    });

                    auto cur = static_cast<uint32_t>(verts.size());
                    append_vertex({
                        .position = clipped->end,
                        .uv       = { at, 0.0f },
                    });

                    index.push_back(prev);
                    index.push_back(cur);

                    flip_flip = !flip_flip;
                    break;
                }
                }

                // Keep the raw event point for the next ray segment. Clipped
                // endpoints are render-only.
                last_point = p;
                visible_point_index += 1;
            }
        }
    }

    qDebug() << Q_FUNC_INFO << "New buffers ready";

    // install all computed geometry

    auto vertex_buffer = QByteArray(reinterpret_cast<const char*>(verts.data()),
                                    verts.size() * sizeof(LineVertex));

    QByteArray index_buffer;

    if (index.size()) {
        index_buffer = QByteArray(reinterpret_cast<const char*>(index.data()),
                                  index.size() * sizeof(uint32_t));
    }


    addAttribute(QQuick3DGeometry::Attribute::PositionSemantic,
                 0,
                 QQuick3DGeometry::Attribute::ComponentType::F32Type);

    addAttribute(QQuick3DGeometry::Attribute::TexCoord0Semantic,
                 3 * sizeof(float),
                 QQuick3DGeometry::Attribute::ComponentType::F32Type);

    setStride(sizeof(LineVertex));
    setVertexData(vertex_buffer);

    if (index.size()) {
        addAttribute(QQuick3DGeometry::Attribute::IndexSemantic,
                     0,
                     QQuick3DGeometry::Attribute::ComponentType::U32Type);

        setIndexData(index_buffer);
    }

    setBounds(m_content_bounds_valid ? m_content_bounds_min : QVector3D(),
              m_content_bounds_valid ? m_content_bounds_max : QVector3D());

    if (index.size()) {
        setPrimitiveType(QQuick3DGeometry::PrimitiveType::Lines);
    } else {
        setPrimitiveType(QQuick3DGeometry::PrimitiveType::Points);
    }

    qDebug() << Q_FUNC_INFO << "Update";
    update();
}

void RayGeometry::inclusion_list_update() {
    qDebug() << Q_FUNC_INFO << "List changed";
    m_include_events = EventTypeContainer(event_include());

    rebuild_geometry();
}

void RayGeometry::entity_filter_update() {
    if (!m_database || !m_database->database || !entity_filter().is_valid()) {
        set_entity_filter_name(QString());
    } else {
        set_entity_filter_name(m_database->database->name_of(entity_filter()));
    }

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
            &RayGeometry::entity_filter_changed,
            this,
            &RayGeometry::entity_filter_update);

    connect(this,
            &RayGeometry::texture_mode_changed,
            this,
            &RayGeometry::rebuild_geometry);

    connect(this,
            &RayGeometry::max_ray_distance_changed,
            this,
            &RayGeometry::rebuild_geometry);

    connect(this,
            &RayGeometry::isect_mode_changed,
            this,
            &RayGeometry::rebuild_geometry);
}

void RayGeometry::set_results(db::SimulationResultPtr data) {
    qDebug() << Q_FUNC_INFO << "New ray geometry database";
    m_database = std::move(data);
    set_entity_filter({});
    const auto available =
        m_database ? static_cast<quint64>(m_database->records.size()) : 0;
    set_available_rays(available);

    const auto default_percent = percent_for_ray_count(
        std::min(DEFAULT_VISIBLE_RAY_COUNT, available), available);
    const bool percent_changed = show_percent() != default_percent;
    set_show_percent(default_percent);
    if (!percent_changed) rebuild_geometry();
}

void RayGeometry::select_entity_filter(db::Entity entity) {
    set_entity_filter(entity);
}

void RayGeometry::clear_entity_filter() {
    set_entity_filter({});
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

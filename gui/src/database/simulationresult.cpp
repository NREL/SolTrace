#include "simulationresult.h"
#include "analysis/ray_volume_raster.h"

#include "database.h"

#include <entt/entity/entity.hpp>

#include <QDebug>

namespace db {

/// Switch between SD and GUI event types
static RayEventType convert(SolTrace::Result::RayEvent e) {
    switch (e) {
    case SolTrace::Result::RayEvent::CREATE: return RayEventType::CREATE;
    case SolTrace::Result::RayEvent::ABSORB: return RayEventType::ABSORB;
    case SolTrace::Result::RayEvent::REFLECT: return RayEventType::REFLECT;
    case SolTrace::Result::RayEvent::TRANSMIT: return RayEventType::TRANSMIT;
    case SolTrace::Result::RayEvent::VIRTUAL: return RayEventType::VIRTUAL;
    case SolTrace::Result::RayEvent::EXIT: return RayEventType::EXIT;
    case SolTrace::Result::RayEvent::UNKNOWN: return RayEventType::UNKNOWN;
    }
    return RayEventType::UNKNOWN;
}

/// Convert from an SD ray record and GUI ray record types.
///
/// Rec is a mut ref as the API does not have a const& version
static RayRecord extract(uint64_t                          id,
                         SimulationResultConversion const& opts,
                         SolTrace::Result::RayRecord&      rec) {
    std::vector<RayEvent> pack;

    pack.reserve(rec.interactions.size());

    auto entity_getter = [&](SolTrace::Data::element_id id) -> entt::entity {
        auto iter = opts.map.find(id);
        if (iter != opts.map.end()) { return iter->second; }
        return entt::null;
    };

    for (auto& c : rec.interactions) {
        pack.push_back({
            .location  = c->location,
            .direction = c->direction,
            .entity    = entity_getter(c->element),
            .event     = convert(c->event),
        });
    }

    return RayRecord {
        .id     = id,
        .events = std::move(pack),
    };
}

SimulationResult::SimulationResult() = default;

SimulationResult::~SimulationResult() = default;


std::unique_ptr<SimulationResult>
SimulationResult::convert(SimulationResultConversion const& opts) {
    qDebug() << Q_FUNC_INFO << "Converting results...";

    auto ret = std::make_unique<SimulationResult>();

    opts.result.get_sun_dimensions(ret->sun_width, ret->sun_height);
    ret->sun_area        = opts.result.get_sun_A_box();
    ret->sun_ray_count   = static_cast<quint64>(opts.result.get_sun_ray_count());
    ret->ray_area_weight = opts.result.get_ray_area_weight();

    ret->records.reserve(opts.result.get_number_of_records());

    uint64_t id = 0;

    // Copy over all rays
    for (auto iter = opts.result.get_ray_record_iterator();
         !opts.result.is_at_end(iter);
         ++iter) {

        ret->records.emplace_back(extract(id, opts, **iter));

        id++;
    }

    // Build lookup table
    for (size_t ray_i = 0; ray_i < ret->records.size(); ray_i++) {
        auto const& events = ret->records[ray_i].events;
        for (auto const& event : events) {
            if (event.entity != entt::null) {
                ret->entity_to_ray_ids[event.entity].push_back(ray_i);
            }
        }
    }

    // Clear dupes
    for (auto& [k, v] : ret->entity_to_ray_ids) {
        std::sort(v.begin(), v.end());
        auto last = std::unique(v.begin(), v.end());
        v.erase(last, v.end());
    }

    // Done

    qDebug() << Q_FUNC_INFO << "Converted" << ret->records.size() << "rays";

    // Bound boxes
    {
        constexpr float maxFloat = std::numeric_limits<float>::max();

        glm::dvec3 bounds_min(maxFloat);
        glm::dvec3 bounds_max(-maxFloat);


        for (auto const& r : ret->records) {
            for (auto const& inter : r.events) {
                bounds_min = glm::min(inter.location, bounds_min);
                bounds_max = glm::max(inter.location, bounds_max);
            }
        }

        if (glm::all(glm::lessThan(bounds_min, bounds_max))) {
            ret->bounds_max = bounds_max;
            ret->bounds_min = bounds_min;
        }
    }

    return ret;
}

SimulationResultModel::SimulationResultModel(QObject* parent)
    : StructModelAdapter { parent } { }

void SimulationResultModel::append_result(SimulationResultPtr result) {
    if (!result) return;

    auto name = result->database->name();

    store_push_append(SimulationResultRecord {
        .name      = name,
        .when      = QDateTime::currentDateTime(),
        .ray_count = static_cast<quint64>(result->records.size()),
        .result    = result,
    });
}

SimulationResultPtr SimulationResultModel::result_at(int index) const {
    auto record = get_at(index);
    if (!record) return {};

    return record->result;
}

QString SimulationResultModel::name_at(int index) const {
    auto record = get_at(index);
    if (!record) return {};

    return record->name;
}

void SimulationResultModel::remove_result(int index) {
    if (index < 0 || index >= rowCount()) return;

    store_push_remove(index, 1);
}

void SimulationResultModel::rename_result(int index, QString const& name) {
    if (index < 0 || index >= rowCount()) return;

    auto record = m_records[index];
    if (record.name == name) return;

    record.name = name;
    store_push_update(index, record);
}

void SimulationResultModel::clear() {
    store_remove_all();
}

} // namespace db

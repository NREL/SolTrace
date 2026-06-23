#pragma once

#include "utilities/grid3d.h"
#include "utilities/structmodel.h"

#include <entt/entity/entity.hpp>
#include <simulation_data.hpp>
#include <simulation_result.hpp>

#include <entt/entity/fwd.hpp>
#include <glm/ext/vector_double3.hpp>

#include <vector>

#include <QDateTime>
#include <QObject>
#include <QPointer>
#include <QString>

namespace db {

class Database;

struct SimulationResultConversion {
    SolTrace::Result::SimulationResult&                                 result;
    SolTrace::Data::SimulationData const&                               data;
    std::unordered_map<SolTrace::Data::element_id, entt::entity> const& map;
};

// can probably make this a variant.
// create and exit dont need entities

enum class RayEventType : uint8_t {
    CREATE   = 1,
    ABSORB   = 2,
    REFLECT  = 3,
    TRANSMIT = 4,
    VIRTUAL  = 5,
    EXIT     = 6,
    UNKNOWN  = UINT8_MAX
};

struct RayEvent {
    glm::dvec3   location;
    glm::dvec3   direction;
    entt::entity entity = entt::null;
    RayEventType event;
};

struct RayRecord {
    uint64_t              id;
    std::vector<RayEvent> events;
};

class SimulationResult {
public:
    SimulationResult();
    ~SimulationResult();

    std::vector<RayRecord> records;

    glm::dvec3 bounds_min;
    glm::dvec3 bounds_max;

    analysis::SparseGrid3D<float> ray_volume;

    std::unordered_map<entt::entity, std::vector<uint64_t>> entity_to_ray_ids;

    std::unique_ptr<Database const> database;

    // TODO: Why is this fallible?
    static std::unique_ptr<SimulationResult>
    convert(SimulationResultConversion const&);
};

using SimulationResultPtr = std::shared_ptr<SimulationResult>;

struct SimulationResultRecord {
    QString             name;
    QDateTime           when;
    quint64             ray_count = 0;
    SimulationResultPtr result;

    RECORD_META(SimulationResultRecord,
                SM_EXPOSE_RO(name),
                SM_EXPOSE_RO(when),
                SM_EXPOSE_RO(ray_count), );
};

class SimulationResultModel
    : public StructModelAdapter<SimulationResultRecord> {
    Q_OBJECT

public:
    explicit SimulationResultModel(QObject* parent = nullptr);

    SimulationResultPtr result_at(int index) const;
    QString             name_at(int index) const;

public slots:
    void append_result(SimulationResultPtr result);
    void remove_result(int index);
    void rename_result(int index, QString const& name);
    void clear();
};

} // namespace db

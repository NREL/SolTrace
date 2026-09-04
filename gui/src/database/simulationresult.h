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

/// Inputs needed to translate library simulation results into GUI ray records.
struct SimulationResultConversion {
    SolTrace::Result::SimulationResult&                                 result;
    SolTrace::Data::SimulationData const&                               data;
    std::unordered_map<SolTrace::Data::element_id, entt::entity> const& map;
};

/// Event type for one point along a traced ray path.
///
/// Some events, such as CREATE and EXIT, may not correspond to a database
/// entity.
enum class RayEventType : uint8_t {
    CREATE   = 1,
    ABSORB   = 2,
    REFLECT  = 3,
    TRANSMIT = 4,
    VIRTUAL  = 5,
    EXIT     = 6,
    UNKNOWN  = UINT8_MAX
};

/// One interaction/location along a ray path.
struct RayEvent {
    glm::dvec3   location;
    glm::dvec3   direction;
    entt::entity entity = entt::null;
    RayEventType event;
};

/// Complete event chain for one simulated ray.
struct RayRecord {
    uint64_t              id;
    std::vector<RayEvent> events;
};

/// GUI-owned simulation result plus lookup data used by analysis modules.
class SimulationResult {
public:
    SimulationResult();
    ~SimulationResult();

    std::vector<RayRecord> records;

    glm::dvec3 bounds_min;
    glm::dvec3 bounds_max;

    double  sun_width       = 0.0;
    double  sun_height      = 0.0;
    double  sun_area        = 0.0;
    quint64 sun_ray_count   = 0;
    double  ray_area_weight = 0.0;

    analysis::SparseGrid3D<float> ray_volume;

    std::unordered_map<entt::entity, std::vector<uint64_t>> entity_to_ray_ids;

    std::unique_ptr<Database const> database;

    // TODO: Why is this fallible?
    /// Convert a library result into GUI result records and entity mappings.
    static std::unique_ptr<SimulationResult>
    convert(SimulationResultConversion const&);
};

using SimulationResultPtr = std::shared_ptr<SimulationResult>;

/// One saved result row exposed to QML.
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

/// QML-facing list model of completed simulation results.
class SimulationResultModel
    : public StructModelAdapter<SimulationResultRecord> {
    Q_OBJECT

public:
    explicit SimulationResultModel(QObject* parent = nullptr);

    /// Return the result pointer at index, or null if index is invalid.
    SimulationResultPtr result_at(int index) const;

    /// Return the display name at index, or an empty string if invalid.
    QString name_at(int index) const;

public slots:
    /// Append a completed result with generated display metadata.
    void append_result(db::SimulationResultPtr result);

    /// Remove a result row.
    void remove_result(int index);

    /// Rename a result row.
    void rename_result(int index, QString const& name);

    /// Remove all result rows.
    void clear();
};

} // namespace db

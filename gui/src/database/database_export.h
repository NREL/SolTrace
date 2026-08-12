#pragma once

#include <entt/entity/entity.hpp>

#include "simulation_data.hpp"

#include <memory>
#include <unordered_map>

namespace db {

class Database;

/// A packed export of a database to the SolTrace data format.
struct DatabaseExport {
    /// SolTrace sim data
    std::shared_ptr<SolTrace::Data::SimulationData> data;

    /// A map of soltrace element ids to entities
    std::unordered_map<SolTrace::Data::element_id, entt::entity> element_map;

    /// A cloned database from where this data came from
    std::unique_ptr<Database> source_database;

    DatabaseExport() = default;

    DatabaseExport(DatabaseExport const&)            = delete;
    DatabaseExport& operator=(DatabaseExport const&) = delete;
    DatabaseExport(DatabaseExport&&)                 = default;
    DatabaseExport& operator=(DatabaseExport&&)      = default;
};

} // namespace db

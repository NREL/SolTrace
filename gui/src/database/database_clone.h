#pragma once

#include "database/database.h"

#include <memory>
#include <unordered_map>

namespace db {

/// Cloned database result content
struct DatabaseCloneResult {
    std::unique_ptr<Database> database;

    /// Fixup map from entities that were in the old database to their new
    /// identities in the cloned map
    std::unordered_map<entt::entity, entt::entity> old_to_new_map;
};

/// Clone a given database, with a given name, and a parent
DatabaseCloneResult clone_database_with_entity_map(Database const& from,
                                                   QString new_database_name,
                                                   QObject* p);

} // namespace db

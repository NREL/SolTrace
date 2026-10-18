#ifndef JSONSCHEMA_H
#define JSONSCHEMA_H

#include <nlohmann/json.hpp>

#include <string>

#include <simulation_data.hpp>

namespace SolTrace::Data {

const std::string kSchemaVersion = "2026.07.15";

void write_json_file(SimulationData& sd, std::string filename);

/// @brief Load simulation data from a JSON file, upgrading older schema
///        versions to the current schema as needed.
/// @param sd SimulationData object to populate.
/// @param filename Path to the JSON file to load.
/// @param upgrade_log Optional pointer to a string that will be populated
///        with a human-readable description of any schema upgrades applied
///        during load. Left untouched if no upgrade was needed. Pass
///        nullptr (default) if this information is not needed.
void load_json_file(SimulationData& sd, std::string filename, std::string* upgrade_log = nullptr);

/// @brief Load simulation data from a JSON string, upgrading older schema
///        versions to the current schema as needed.
/// @param sd SimulationData object to populate.
/// @param json_str The JSON string to load.
/// @param upgrade_log Optional pointer to a string that will be populated
///        with a human-readable description of any schema upgrades applied
///        during load. Left untouched if no upgrade was needed. Pass
///        nullptr (default) if this information is not needed.
void load_json_cstr(SimulationData& sd, const char* json_str, std::string* upgrade_log = nullptr);
}

#endif
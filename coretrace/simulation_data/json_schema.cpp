#include "json_schema.hpp"

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <map>
#include <functional>

#include <json_helpers.hpp>
#include <single_element.hpp>
#include <stage_element.hpp>
#include <sun.hpp>



namespace SolTrace::Data {

// Private

bool upgrade_20251112_20260715(nlohmann::ordered_json& jroot,
                                std::string&            change_log) {
    using json = nlohmann::ordered_json;

    // Major change: Optical properties consolidated to shared
    // group in SimulationData. Previously each single element
    // carried its own inline "optics_front"/"optics_back" nodes;
    // now elements reference a shared "opt_id" that indexes into
    // a top level "optical_properties" node.

    // Only handle the specific version this function is named for.
    std::string schema_version = jroot.value("schema_version", "");
    if (schema_version != "2025.11.12") {
        return false;
    }

    try {
        std::map<std::string, optics_id> signature_to_id;
        json joptics_top = json::object();
        optics_id next_id = 0;

        std::function<void(json&)> process_element = [&](json& jelement) {
            if (jelement.contains("elements") &&
                jelement["elements"].is_object()) {
                for (auto& [child_key, jchild] : jelement["elements"].items()) {
                    process_element(jchild);
                }
            }

            if (!jelement.contains("optics_front") ||
                !jelement.contains("optics_back")) {
                return;
            }

            json joptics_front = jelement.at("optics_front");
            json joptics_back  = jelement.at("optics_back");

            json jcombined = json::object();
            jcombined["my_type"] = joptics_front.at("my_type");
            jcombined["refraction_index_front"] =
                joptics_front.at("refraction_index_front");
            jcombined["refraction_index_back"] =
                joptics_front.at("refraction_index_back");

            json jfront = json::object();
            jfront["error_distribution_type"] =
                joptics_front.at("error_distribution_type");
            jfront["transmissivity"]    = joptics_front.at("transmissivity");
            jfront["reflectivity"]      = joptics_front.at("reflectivity");
            jfront["slope_error"]       = joptics_front.at("slope_error");
            jfront["specularity_error"] = joptics_front.at("specularity_error");

            json jback = json::object();
            jback["error_distribution_type"] =
                joptics_back.at("error_distribution_type");
            jback["transmissivity"]    = joptics_back.at("transmissivity");
            jback["reflectivity"]      = joptics_back.at("reflectivity");
            jback["slope_error"]       = joptics_back.at("slope_error");
            jback["specularity_error"] = joptics_back.at("specularity_error");

            jcombined["front"] = jfront;
            jcombined["back"]  = jback;

            std::string signature = jcombined.dump();

            optics_id opt_id;
            auto      found = signature_to_id.find(signature);
            if (found != signature_to_id.end()) {
                opt_id = found->second;
            } else {
                opt_id = next_id++;
                signature_to_id[signature] = opt_id;

                std::stringstream ss;
                ss << "Legacy Optics " << opt_id;
                jcombined["my_name"] = ss.str();

                joptics_top[std::to_string(opt_id)] = jcombined;
            }

            jelement.erase("optics_front");
            jelement.erase("optics_back");
            jelement["opt_id"] = opt_id;
        };

        if (!jroot.contains("elements") || !jroot["elements"].is_object()) {
            return false;
        }

        for (auto& [key, jelement] : jroot["elements"].items()) {
            process_element(jelement);
        }

        jroot["optical_properties"] = joptics_top;
        jroot["schema_version"] = "2026.07.15";

        // Describe what changed, for user-facing reporting.
        std::stringstream log;
        log << "[2025.11.12 -> 2026.07.15] Consolidated per-element "
               "optical properties into a shared \"optical_properties\" "
               "node. Extracted " << joptics_top.size()
            << " unique optical property set(s); each affected element "
               "now references its set via \"opt_id\" instead of "
               "inline \"optics_front\"/\"optics_back\" data.\n";
        change_log += log.str();
    }
    catch (const json::exception& e) {
        std::cerr << "upgrade_20251112_20260715 failed: " << e.what() << '\n';
        return false;
    }

    return true;
}

// Public

void write_json_file(SimulationData& sd, std::string filename) {
    using json = nlohmann::ordered_json;

    // Create empty object
    json root;

    // Write general meta data
    {
        root["schema_version"]     = kSchemaVersion;
        root["number_of_elements"] = sd.get_number_of_elements();
    }

    // Write parameters
    {
        json jpar;
        sd.get_simulation_parameters().write_json(jpar);
        root["simulation_parameters"] = jpar;
    }

    // Write ray sources
    {
        json jsources;
        for (auto it = sd.get_ray_source_iterator();
             !sd.is_ray_source_at_end(it);
             ++it) {
            json jsrc;

            SolTrace::Data::ray_source_id i          = it->first;
            auto                          ray_source = it->second;

            // Check source type
            if (auto sun_ptr = std::dynamic_pointer_cast<SolTrace::Data::Sun>(
                    ray_source)) {
                sun_ptr->write_json(jsrc);
            } else {
                // UNSUPPORTED type
                throw std::runtime_error("Unsupported ray source type");
            }

            jsources[std::to_string(i)] = jsrc;
        }

        root["ray_sources"] = jsources;
    }

    // Write optical properties
    {
        json joptics_top;
        for (auto it = sd.get_optics_iterator(); !sd.is_optics_at_end(it);
             ++it) {
            json joptics;

            SolTrace::Data::optics_id id         = it->first;
            auto                      optics_set = it->second;

            optics_set->write_json(joptics);

            joptics_top[std::to_string(id)] = joptics;
        }
        root["optical_properties"] = joptics_top;
    }

    // Write Elements
    {
        json jelements_top;
        int  i_top = 0;
        for (auto it = sd.get_iterator(); !sd.is_at_end(it); ++it) {
            json jelement;
            auto element = it->second;

            if (element->is_single() && (element->is_top_level() == false)) {
                // Skip single elements that are within other stages/composites
                continue;
            }

            element->write_json(jelement);

            jelements_top[std::to_string(i_top)] = jelement;
            i_top++;
        }

        root["elements"] = jelements_top;
    }

    // Write to disk
    std::ofstream ofs(filename, std::ios::out | std::ios::trunc);
    if (!ofs.is_open()) throw std::runtime_error("Failure writing json");
    ofs << root.dump(kJsonIndentSpaces) << '\n';

    return;
}

void load_json_file(SimulationData& sd, std::string filename, std::string* upgrade_log) {
    using json = nlohmann::ordered_json;

    // Clear simulation data
    sd.clear();

    // Load json file
    std::ifstream ifs(filename);
    if (!ifs.is_open()) throw std::runtime_error("Failure opening json");

    // Load json from file stream
    json root;
    ifs >> root;

    // File meta data
    std::string schema_version     = root.at("schema_version");

    // Upgrade to modern version
    std::string local_upgrade_log;
    if (schema_version == "2025.11.12") {
        if (!upgrade_20251112_20260715(root, local_upgrade_log)) {
            throw std::runtime_error(
                "Failed to upgrade JSON file schema from 2025.11.12 to 2026.07.15");
        }
    }

    if (upgrade_log != nullptr && !local_upgrade_log.empty()) {
        *upgrade_log = local_upgrade_log;
    }

    // Check file is up to date
    if (root.at("schema_version") != kSchemaVersion) {
        std::stringstream ss;
        ss << "Unsupported or unrecognized schema version: "
           << root.at("schema_version").get<std::string>();
        throw std::runtime_error(ss.str());
    }

    // Simulation parameters
    SolTrace::Data::SimulationParameters& sim_par =
        sd.get_simulation_parameters();
    json jpar = root["simulation_parameters"];
    sim_par   = SolTrace::Data::SimulationParameters(jpar);

    // Ray sources
    json jsources = root["ray_sources"];
    for (auto& [key, jsrc] : jsources.items()) {
        std::string source_type = jsrc.at("source_type");
        if (source_type != "Sun") {
            // UNSUPPORTED source type
            throw std::runtime_error("Unsupported ray source type");
        }

        // Make sun for simulation data
        auto sun = make_ray_source<Sun>(jsrc);
        sd.add_ray_source(sun);
    }

    // Optical properties
    json joptics = root.at("optical_properties");
    for (auto& [key, joptic] : joptics.items()) {
        optics_id          opt_id = static_cast<optics_id>(std::stoll(key));
        OpticalPropertySet opt_set(joptic);

        // Check for pre-existing optical property sets
        const OpticalPropertySet* existing =
            sd.get_optical_property_set(opt_id);
        if (existing != nullptr) {
            // This should be a built in optical property set
            if (opt_id >= 0)
                throw std::runtime_error(
                    "Custom optical property set already exists");

            // Ensure loaded built in matches
            if (*existing != opt_set) {
                std::stringstream ss;
                ss << "Built-in optical property set mismatch for id "
                   << opt_id;
                throw std::runtime_error(ss.str());
            }
        } else // Insert new optical property set
        {
            auto ptr  = std::make_shared<OpticalPropertySet>(opt_set);
            bool flag = sd.my_optical_property_sets.insert_item(opt_id, ptr);
            if (!flag) {
                std::stringstream ss;
                ss << "Failed to insert optical property set id from JSON: "
                   << opt_id;
                throw std::runtime_error(ss.str());
            }
        }
    }
    // Set optical property set next id
    sd.my_optical_property_sets.recompute_next_id(0);

    // Elements
    json jelements      = root["elements"];
    auto resolve_optics = [&sd](const optics_id id) {
        auto ptr = sd.my_optical_property_sets.get_item(id);
        return OpticalPropertySetReference { id, ptr };
    };
    for (auto& [key, jelement] : jelements.items()) {
        // Check if stage
        // Note a stage is also a composite, so check stage first
        if (jelement.contains("is_stage") && jelement.at("is_stage") == true) {
            // Make stage
            stage_ptr stage = make_stage(jelement, resolve_optics);
            sd.add_stage(stage);
        }
        // Composite
        else if (jelement.contains("is_composite") &&
                 jelement.at("is_composite") == true) {
            composite_element_ptr comp =
                make_element<CompositeElement>(jelement, resolve_optics);
            sd.add_element(comp);
        }
        // Single Element
        else {
            single_element_ptr single =
                make_element<SingleElement>(jelement, resolve_optics);
            sd.add_element(single);
        }
    }

    return;
}



} // namespace SolTrace::Data
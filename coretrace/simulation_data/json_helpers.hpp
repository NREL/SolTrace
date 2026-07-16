/**
 * @file json_helpers.hpp
 * @brief Lightweight helpers for extracting values from JSON nodes.
 *
 * Defines utility functions used to read specific typed values from
 * nlohmann::ordered_json objects while providing simple fallback behavior.
 * Intended for internal parsing of simulation configuration structures
 * (not full file serialization).
 *
 * Add additional typed accessors here as needed.
 */

#ifndef JSONHELPERS_H
#define JSONHELPERS_H

#include <string>

#include <nlohmann/json.hpp>
#include <glm/vec3.hpp>

namespace SolTrace::Data {

    const int kJsonIndentSpaces = 4;

    inline double json_get_double(const nlohmann::ordered_json& jnode, std::string key)
    {
        auto jval = jnode.at(key);
        if (jval.is_null())
            return std::numeric_limits<double>::quiet_NaN();
        else
            return jval;
    }

    inline
    std::array<double, 3> to_array(glm::dvec3 v)
    {
        return {v.x, v.y, v.z};
    }

    inline
    glm::dvec3 from_array(std::array<double, 3> array)
    {
        return {array[0], array[1], array[2]};
    }


} // namespace SolTrace::Data

#endif

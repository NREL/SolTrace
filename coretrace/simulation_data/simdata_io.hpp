/**
 * @file simdata_io.hpp
 * @brief Simulation data input/output operations
 *
 * Provides functions for reading and writing simulation data
 * to various file formats and data exchange mechanisms.
 * Includes support for SolTrace input files (.stinput) and
 * other data formats.
 */

#ifndef SIMDATA_H
#define SIMDATA_H

#include <string>
#include "simulation_data.hpp"

namespace SolTrace::Data {

/**
* @brief Lookup an enumeration value from its string representation.
*
* Performs a reverse lookup on a forward (enum->string) map to find the
* enumeration value whose associated string exactly matches @p str.
* If no match is found, the provided @p unknown value is returned.
*
* @tparam EnumT Enumeration type.
* @param str String to search for.
* @param forward_map Map from enumeration values to their string names.
* @param unknown Fallback value returned when @p str is not present.
* @return Matching enumeration value, or @p unknown if not found.
*
*/
template<typename EnumT>
EnumT get_enum_from_string(const std::string& str, 
    const std::map<EnumT, std::string>& forward_map, 
    EnumT unknown)
{
    for (const auto& key_pair : forward_map)
    {
        if (key_pair.second == str)
        {
            return key_pair.first;
        }
    }

    return unknown;
}
bool load_stinput_file(SimulationData& sd, std::string filename);

void write_json_file(SimulationData& sd, std::string filename);

void load_json_file(SimulationData& sd, std::string filename);

} // namespace SolTrace::Data

#endif

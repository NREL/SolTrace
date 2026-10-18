/**
 * @file simulation_parameters.hpp
 * @brief Simulation parameter definitions and validation
 *
 * Defines parameters that control ray tracing simulation behavior,
 * including ray counts, error modeling options, output settings,
 * and convergence criteria. Provides validation and default values
 * for simulation configuration.
 */

#ifndef SOLTRACE_SIMULATION_PARAMETERS_H
#define SOLTRACE_SIMULATION_PARAMETERS_H

#include <cstdint>
#include <nlohmann/json.hpp>

namespace SolTrace::Data {

struct SimulationParameters
{
public:
    // TODO: Figure out how to store time...
    DateTime sim_dt;

    std::uint_fast64_t number_of_rays;
    std::uint_fast64_t max_number_of_rays;
    double tolerance;

    double latitude;
    double longitude;

    int seed;

    bool include_sun_shape_errors;
    bool include_optical_errors;

    bool as_power_tower;

    SimulationParameters() : number_of_rays(10000),
                             max_number_of_rays(1000000),
                             tolerance(0.0),
                             latitude(0.0),
                             longitude(0.0),
                             seed(0),
                             include_sun_shape_errors(false),
                             include_optical_errors(false),
                             as_power_tower(false)
    {
    }
    SimulationParameters(const nlohmann::ordered_json& jnode)
    {
        this->include_sun_shape_errors = jnode.at("include_sun_shape_errors");
        this->include_optical_errors = jnode.at("include_optical_errors");
        this->number_of_rays = jnode.at("number_of_rays");
        this->max_number_of_rays = jnode.at("max_number_of_rays");
        this->tolerance = jnode.at("tolerance");
        this->latitude = jnode.at("latitude");
        this->longitude = jnode.at("longitude");
        this->seed = jnode.at("seed");
    }
    ~SimulationParameters() {}

    void write_json(nlohmann::ordered_json& jnode) const
    {
        jnode["include_sun_shape_errors"] = this->include_sun_shape_errors;   // bool
        jnode["include_optical_errors"] = this->include_optical_errors;       // bool
        jnode["number_of_rays"] = this->number_of_rays;                       // int
        jnode["max_number_of_rays"] = this->max_number_of_rays;               // int
        jnode["tolerance"] = this->tolerance;                                 // double
        jnode["latitude"] = this->latitude;                                   // double
        jnode["longitude"] = this->longitude;                                 // double
        jnode["seed"] = this->seed;                                           // int
    }
};

// TODO: Implement the output stream operator.

} // namespace SolTrace::Data

#endif

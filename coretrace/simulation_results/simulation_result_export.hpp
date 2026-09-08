#ifndef SOLTRACE_SIMULATION_RESULT_EXPORT_H
#define SOLTRACE_SIMULATION_RESULT_EXPORT_H

#include "simulation_result_api.hpp"

// Types
using SolTrace::Result::interaction_ptr;
using SolTrace::Result::ElementRecord;
using SolTrace::Result::element_record_ptr;
using SolTrace::Result::InteractionRecord;
using SolTrace::Result::interaction_ptr;
using SolTrace::Result::ray_id;
using SolTrace::Result::RayEvent;
using SolTrace::Result::RayRecord;
using SolTrace::Result::ray_record_ptr;
using SolTrace::Result::SimulationResult;
using SolTrace::Result::GroupResult;

// Functions
using SolTrace::Result::make_element_record;
using SolTrace::Result::make_interaction_record;
using SolTrace::Result::make_ray_record;
using SolTrace::Result::ray_event_string;

#endif

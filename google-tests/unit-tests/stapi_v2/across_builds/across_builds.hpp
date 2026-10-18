/*
testing coverage across builds (x: done, -: skipped)

function \ build                        : native | embree | optix | all

call_stapi_v2_read_input_json              [x]      [x]      [x]    [x]
call_stapi_v2_read_input_json_file         [x]      [x]      [x]    [x]
call_stapi_v2_set_simulation_parameters    [x]      [x]      [x]    [x]
call_stapi_v2_sim_params                   [x]      [x]      [x]    [x]
call_stapi_v2_sim_rays                     [x]      [x]      [x]    [x]
call_stapi_v2_sim_power_tower              [x]      [x]      [x]    [x]
call_stapi_v2_sim_errors                   [x]      [x]      [x]    [x]
call_stapi_v2_sim_location                 [x]      [x]      [x]    [x]
call_stapi_v2_sim_tolerance                [x]      [x]      [x]    [x]
call_stapi_v2_get_simulation_parameters    [x]      [x]      [x]    [x]
call_stapi_v2_add_optics                   [x]      [x]      [x]    [x]
call_stapi_v2_get_optic                    [x]      [x]      [x]    [x]
call_stapi_v2_remove_optics                [x]      [x]      [x]    [x]
call_stapi_v2_add_elements                 [x]      [x]      [x]    [x]
call_stapi_v2_get_element                  [x]      [x]      [x]    [x]
call_stapi_v2_remove_elements              [x]      [x]      [x]    [x]
call_stapi_v2_toggle_element               [x]      [x]      [x]    [x]
call_stapi_v2_element_xyz                  [x]      [x]      [x]    [x]
call_stapi_v2_element_aim                  [x]      [x]      [x]    [x]
call_stapi_v2_element_zrot                 [x]      [x]      [x]    [x]
call_stapi_v2_element_aperture             [x]      [x]      [x]    [x]
call_stapi_v2_element_surface              [x]      [x]      [x]    [x]
call_stapi_v2_element_optic                [x]      [x]      [x]    [x]
call_stapi_v2_element_group                [x]      [x]      [x]    [x]
call_stapi_v2_add_sun                      [x]      [x]      [x]    [x]
call_stapi_v2_get_sun                      [x]      [x]      [x]    [x]
call_stapi_v2_sun_shape                    [x]      [x]      [x]    [x]
call_stapi_v2_sun_xyz                      [x]      [x]      [x]    [x]
call_stapi_v2_sun_userdata                 [x]      [x]      [x]    [x]
call_stapi_v2_solar_calculator             [x]      [x]      [x]    [x]
call_stapi_v2_sim_setup                    [x]      [x]      [x]    [x]
call_stapi_v2_sim_run_v2                   [x]      [x]      [x]    [x]
call_stapi_v2_sim_report                   [x]      [x]      [x]    [x]
call_stapi_v2_write_results_csv            [x]      [x]      [x]    [x]
call_stapi_v2_write_group_results_json     [x]      [x]      [x]    [x]
call_stapi_v2_locations                    [x]      [x]      [x]    [x]
call_stapi_v2_cosines                      [x]      [x]      [x]    [x]
call_stapi_v2_elementmap                   [x]      [x]      [x]    [x]
call_stapi_v2_stagemap                     [x]      [x]      [x]    [x]
call_stapi_v2_raynumbers                   [x]      [x]      [x]    [x]
call_stapi_v2_sun_stats                    [x]      [x]      [x]    [x]
call_stapi_v2_get_results_data             [x]      [x]      [x]    [x]
*/

#ifndef STAPI_V2_ACROSS_BUILDS_H
#define STAPI_V2_ACROSS_BUILDS_H

#include "../../../../api/stapi_v2.h"

using json = nlohmann::ordered_json;

#define SETUP_TEST_CXT()                         \
    st_context_v2_t pcxt;                        \
    st_return_t code = st_create_context(&pcxt); \
    EXPECT_EQ(code, st_return_code::SUCCESS);    \
    EXPECT_NE(pcxt, nullptr);

#define CLEANUP_TEST_CXT()                   \
    code = st_free_context(pcxt);            \
    EXPECT_EQ(code, st_return_code::SUCCESS);

#define LOAD_TEST_JSON()                                  \
    json root = load_json();                              \
    code = st_read_input_json(pcxt, root.dump().c_str()); \
    EXPECT_EQ(code, st_return_code::SUCCESS);

///////////////////////
// Utility Functions //
///////////////////////

json load_json();

////////////////////////////////
// Simlulation Data Functions //
////////////////////////////////

// functions for simulation data management thru json strings
st_return_t call_stapi_v2_read_input_json(st_context_v2_t pcxt);
st_return_t call_stapi_v2_read_input_json_file(st_context_v2_t pcxt);

// functions for simulation data management directly
st_return_t call_stapi_v2_set_simulation_parameters(st_context_v2_t pcxt);
st_return_t call_stapi_v2_sim_params(st_context_v2_t pcxt);
st_return_t call_stapi_v2_sim_rays(st_context_v2_t pcxt);
st_return_t call_stapi_v2_sim_power_tower(st_context_v2_t pcxt);
st_return_t call_stapi_v2_sim_errors(st_context_v2_t pcxt);
st_return_t call_stapi_v2_sim_location(st_context_v2_t pcxt);
st_return_t call_stapi_v2_sim_tolerance(st_context_v2_t pcxt);
st_return_t call_stapi_v2_get_simulation_parameters(st_context_v2_t pcxt);

// function for optical properties
st_return_t call_stapi_v2_add_optics(st_context_v2_t pcxt);
st_return_t call_stapi_v2_get_optic(st_context_v2_t pcxt);
st_return_t call_stapi_v2_remove_optics(st_context_v2_t pcxt);

// functions for elements
st_return_t call_stapi_v2_add_elements(st_context_v2_t pcxt);
st_return_t call_stapi_v2_get_element(st_context_v2_t pcxt);
st_return_t call_stapi_v2_remove_elements(st_context_v2_t pcxt);
st_return_t call_stapi_v2_toggle_element(st_context_v2_t pcxt);
st_return_t call_stapi_v2_element_xyz(st_context_v2_t pcxt);
st_return_t call_stapi_v2_element_aim(st_context_v2_t pcxt);
st_return_t call_stapi_v2_element_zrot(st_context_v2_t pcxt);
st_return_t call_stapi_v2_element_aperture(st_context_v2_t pcxt);
st_return_t call_stapi_v2_element_surface(st_context_v2_t pcxt);
st_return_t call_stapi_v2_element_geometry(st_context_v2_t pcxt);
st_return_t call_stapi_v2_element_optic(st_context_v2_t pcxt);
st_return_t call_stapi_v2_element_group(st_context_v2_t pcxt);

// sun functions
st_return_t call_stapi_v2_add_sun(st_context_v2_t pcxt);
st_return_t call_stapi_v2_get_sun(st_context_v2_t pcxt);
st_return_t call_stapi_v2_sun_shape(st_context_v2_t pcxt);
st_return_t call_stapi_v2_sun_xyz(st_context_v2_t pcxt);
st_return_t call_stapi_v2_sun_userdata(st_context_v2_t pcxt);
st_return_t call_stapi_v2_solar_calculator(st_context_v2_t pcxt);

//////////////////////////////////
// Simlulation Runner Functions //
//////////////////////////////////

st_return_t call_stapi_v2_sim_setup(st_context_v2_t pcxt);
st_return_t call_stapi_v2_sim_run_v2(st_context_v2_t  pcxt,
                                     st_runner_type_t runner_type);
st_return_t call_stapi_v2_sim_report(st_context_v2_t  pcxt, 
                                     st_runner_type_t runner_type);

///////////////////////////////////
// Simlulation Results Functions //
///////////////////////////////////

// TODO: add test for st_write_results_json
st_return_t call_stapi_v2_write_results_csv(st_context_v2_t  pcxt, 
                                            st_runner_type_t runner_type, 
                                            const char       *filename);
st_return_t call_stapi_v2_write_group_results_json(st_context_v2_t  pcxt, 
                                                   st_runner_type_t runner_type, 
                                                   const char       *filename);
// functions to get results directly
st_return_t call_stapi_v2_locations(st_context_v2_t  pcxt,
                                    st_runner_type_t runner_type);
st_return_t call_stapi_v2_cosines(st_context_v2_t  pcxt,
                                  st_runner_type_t runner_type);
st_return_t call_stapi_v2_elementmap(st_context_v2_t  pcxt,
                                     st_runner_type_t runner_type);
st_return_t call_stapi_v2_stagemap(st_context_v2_t  pcxt,
                                   st_runner_type_t runner_type);
st_return_t call_stapi_v2_raynumbers(st_context_v2_t  pcxt,
                                     st_runner_type_t runner_type);
st_return_t call_stapi_v2_sun_stats(st_context_v2_t  pcxt,
                                    st_runner_type_t runner_type);
st_return_t call_stapi_v2_get_results_data(st_context_v2_t  pcxt,
                                           st_runner_type_t runner_type);

#endif

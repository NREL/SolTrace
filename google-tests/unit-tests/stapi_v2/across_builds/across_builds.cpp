/*
functions that are across different api builds (only native, with only embree, with only optix, with both embree and optix).
i.e. shouldn't need google-test.h in here.
intended to be called from tests inside different api build tests, those manage the context that these functions use.
*/

#include <fstream>
#include <vector>

#include "across_builds.hpp"

///////////////////////
// Utility Functions //
///////////////////////

json load_json()
{
    namespace fs = std::filesystem;

    // Build paths
    const fs::path project_root(PROJECT_DIR);
    const fs::path sample_path = project_root / "sample_ver_20251112.json";
    // const fs::path project_root(PROJECT_DIR);
    // const std::string input_str = project_root.string() + "/good_test.json";

    std::ifstream ifs(sample_path);
    if (!ifs.is_open()) throw std::runtime_error("Failure opening json");

    // Load json from file stream
    json root;
    ifs >> root;
    return root;
}

auto check = [](auto val, auto should_be)
{
    // false = 0, true = 1, to make successful test = 0 
    // -> !(a == b) returns 0 if the are equal, 1 otherwise
    return (st_return_t)!(val == (decltype(val))should_be); 
};

auto check_not = [](auto val, auto should_be)
{
    // false = 0, true = 1, to make successful test = 0 
    // -> a == b returns 0 if the are not equal, 1 otherwise
    return (st_return_t)(val == (decltype(val))should_be); 
};

auto check_near = [](auto val, auto should_be_near, double abs_tol)
{
    // false = 0, true = 1, to make successful test = 0 
    // -> (abs(a - b) > abs_tol) returns 0 if the diff 
    //    is less than the tolerance, 1 otherwise
    return (st_return_t)((double)abs(val - (decltype(val))should_be_near) > abs_tol); 
};

st_return_t check_optical_side(optical_set_ptr  set,
                               OpticalSide      side,
                               DistributionType dist,
                               double           tran,
                               double           refl,
                               double           slope,
                               double           spec)
{
    if (set->get_error_distribution(side) != dist
        || set->get_reflectivity(side) != refl
        || set->get_transmissivity(side) != tran
        || set->get_slope_error(side) != slope
        || set->get_specularity_error(side) != spec)
    {
        return 1;
    }
    return 0;
}

////////////////////////////////
// Simlulation Data Functions //
////////////////////////////////

// functions for simulation data management thru json strings
st_return_t call_stapi_v2_read_input_json(st_context_v2_t pcxt)
{
    json root = load_json();
    // add calls to other overloads
    return st_read_input_json(pcxt, root.dump().c_str());
}

st_return_t call_stapi_v2_read_input_json_file(st_context_v2_t pcxt)
{
    namespace fs = std::filesystem;

    // Build paths
    const fs::path project_root(PROJECT_DIR);
    const fs::path sample_path = project_root / "sample_ver_20251112.json";
    return st_read_input_json_file(pcxt, sample_path.string().c_str());
}

// functions for simulation data management directly
st_return_t call_stapi_v2_set_simulation_parameters(st_context_v2_t pcxt)
{
    args_simulation_parameters _params = { 1, 100, .1, 35.962278, -106.5122622, true, true, false };

    // expect += st_return_code::SUCCESS
    st_return_t code = st_set_simulation_parameters(pcxt, &_params);

    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    SimulationParameters &params = cxt->p_data->get_simulation_parameters();
    code += check(params.number_of_rays,           _params.number_of_rays);
    code += check(params.max_number_of_rays,       _params.max_number_of_rays);
    code += check(params.tolerance,                _params.tolerance);
    code += check(params.latitude,                 _params.latitude);
    code += check(params.longitude,                _params.longitude);
    code += check(params.include_sun_shape_errors, _params.include_sun_shape_errors);
    code += check(params.include_optical_errors,   _params.include_optical_errors);
    code += check(params.as_power_tower,           _params.as_power_tower);

    return code;
}

st_return_t call_stapi_v2_sim_params(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    SimulationParameters &params = cxt->p_data->get_simulation_parameters();

    // test that initialized with defaults
    st_return_t code = check(params.number_of_rays, 10000);
    code += check(params.max_number_of_rays, 1000000);
    code += check(params.as_power_tower, false);

    st_sim_params(pcxt, 1, 100, true);

    code += check(params.number_of_rays, 1);
    code += check(params.max_number_of_rays, 100);
    code += check(params.as_power_tower, true);
    return code;
}

st_return_t call_stapi_v2_sim_rays(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    SimulationParameters &params = cxt->p_data->get_simulation_parameters();

    // test that initialized with defaults
    st_return_t code = check(params.number_of_rays, 10000);
    code += check(params.max_number_of_rays, 1000000);
    
    st_sim_rays(pcxt, 1, 100);

    code += check(params.number_of_rays, 1);
    code += check(params.max_number_of_rays, 100);
    return code;
}

st_return_t call_stapi_v2_sim_power_tower(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    SimulationParameters &params = cxt->p_data->get_simulation_parameters();

    // test that initialized with defaults
    st_return_t code = check(params.as_power_tower, false);

    st_sim_power_tower(pcxt, true);

    code += check(params.as_power_tower, true);
    return code;
}

st_return_t call_stapi_v2_sim_errors(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    SimulationParameters &params = cxt->p_data->get_simulation_parameters();

    // test that initialized with defaults
    st_return_t code = check(params.include_sun_shape_errors, false);
    code += check(params.include_optical_errors, false);

    st_sim_errors(pcxt, true, true);

    code += check(params.include_sun_shape_errors, true);
    code += check(params.include_optical_errors, true);
    return code;
}

st_return_t call_stapi_v2_sim_location(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    SimulationParameters &params = cxt->p_data->get_simulation_parameters();

    // test that initialized with defaults
    st_return_t code = check(params.latitude, 0);
    code += check(params.longitude, 0);

    st_sim_location(pcxt, 35.962278, -106.5122622);

    code += check(params.latitude, 35.962278);
    code += check(params.longitude, -106.5122622);
    return code;
}

st_return_t call_stapi_v2_sim_tolerance(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    SimulationParameters &params = cxt->p_data->get_simulation_parameters();

    // test that initialized with defaults
    st_return_t code = check(params.tolerance, 0);

    st_sim_tolerance(pcxt, .1);

    code += check(params.tolerance, .1);
    return code;
}

st_return_t call_stapi_v2_get_simulation_parameters(st_context_v2_t pcxt)
{
    args_simulation_parameters _params = { 1, 100, .1, 35.962278, -106.5122622, true, true, false };

    // expect += st_return_code::SUCCESS
    st_return_t code = st_set_simulation_parameters(pcxt, &_params);
    
    args_simulation_parameters rt_params;
    code += st_get_simulation_parameters(pcxt, &rt_params);
    
    code += check(rt_params.number_of_rays,           _params.number_of_rays);
    code += check(rt_params.max_number_of_rays,       _params.max_number_of_rays);
    code += check(rt_params.tolerance,                _params.tolerance);
    code += check(rt_params.latitude,                 _params.latitude);
    code += check(rt_params.longitude,                _params.longitude);
    code += check(rt_params.include_sun_shape_errors, _params.include_sun_shape_errors);
    code += check(rt_params.include_optical_errors,   _params.include_optical_errors);
    code += check(rt_params.as_power_tower,           _params.as_power_tower);

    return code;
}

// functions to add/remove/set optical properties
st_return_t call_stapi_v2_add_optics(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    SimulationData *data = cxt->p_data;
    auto get_last_opt_set = [data]()
    {
        optical_set_ptr opt_set;
        for (auto it = data->get_optics_iterator(); !data->is_optics_at_end(it); ++it)
            opt_set = it->second;
        return opt_set;
    };

    args_optical_properties_set set = {"test1", 1.1, 1.1, 1};
    args_optical_properties_face f = {.5, .5, .5, .5, 'g'};
    args_optical_properties_face b = {.25, .25, .25, .25, 'g'};

    // check no optics set
    uint_fast64_t num = -1;
    uint_fast64_t id = -1;
    st_num_optics(pcxt, &num);
    st_return_t code = check(num, 0);

    // test add refractive set
    // expect += st_return_code::SUCCESS
    code += st_add_optical_properties_set(pcxt, &set, &f, &b, &id);
    code += check(id, 0);
    st_num_optics(pcxt, &num);
    code += check(num, 1);

    // check values were set
    auto last_opt_set = get_last_opt_set();
    code += check(last_opt_set->get_name(), "test1");
    code += check(last_opt_set->get_interaction_type(), InteractionType::REFRACTION);
    code += check_optical_side(last_opt_set, OpticalSide::Front, DistributionType::GAUSSIAN, .5, .5, .5, .5);
    code += check_optical_side(last_opt_set, OpticalSide::Back, DistributionType::GAUSSIAN, .25, .25, .25, .25);

    // test add reflective set
    set.type = 2;
    set.name = "test2";
    // expect += st_return_code::SUCCESS
    code += st_add_optical_properties_set(pcxt, &set, &f, &b, &id);
    code += check(id, 1);
    st_num_optics(pcxt, &num);
    code += check(num, 2);

    // check values were set
    last_opt_set = get_last_opt_set();
    code += check(last_opt_set->get_name(), "test2");
    code += check(last_opt_set->get_interaction_type(), InteractionType::REFLECTION);
    code += check_optical_side(last_opt_set, OpticalSide::Front, DistributionType::GAUSSIAN, .5, .5, .5, .5);
    code += check_optical_side(last_opt_set, OpticalSide::Back, DistributionType::GAUSSIAN, .25, .25, .25, .25);

    // set bad distribution characters for front
    f.error_distribution_type = 'z';
    // expect += st_return_code::INVALID_ARGUMENTS
    code += st_add_optical_properties_set(pcxt, &set, &f, &b, &id);
    code += check(id, 1);
    st_num_optics(pcxt, &num);
    code += check(num, 2);

    // set bad distribution characters for back
    f.error_distribution_type = 'g';
    b.error_distribution_type = 'z';
    // expect += st_return_code::INVALID_ARGUMENTS
    code += st_add_optical_properties_set(pcxt, &set, &f, &b, &id);
    code += check(id, 1);
    st_num_optics(pcxt, &num);
    code += check(num, 2);

    // expect == 2 * st_return_code::INVALID_ARGUMENTS
    return code;
}

st_return_t call_stapi_v2_get_optic(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    SimulationData *data = cxt->p_data;
    args_optical_properties_set set = {"test1", 1.1, 1.1, 2};
    args_optical_properties_face f = {.5, .5, .5, .5, 'g'};
    args_optical_properties_face b = {.25, .25, .25, .25, 'g'};
    uint_fast64_t id = -1;

    // expect += st_return_code::SUCCESS
    st_return_t code = st_add_optical_properties_set(pcxt, &set, &f, &b, &id);

    args_optical_properties_set rt_set;
    args_optical_properties_face rt_f;
    args_optical_properties_face rt_b;

    // expect += st_return_code::SUCCESS
    code += st_get_optical_properties_set(pcxt, id, &rt_set, &rt_f, &rt_b);

    // check set struct
    std::string rt_str(rt_set.name);
    std::string set_str(set.name);
    code += check(rt_str, set_str);
    code += check(rt_set.refraction_index_front, set.refraction_index_front);
    code += check(rt_set.refraction_index_back, set.refraction_index_back);
    code += check(rt_set.type, set.type);
    
    // check front struct
    code += check(rt_f.error_distribution_type, f.error_distribution_type);
    code += check(rt_f.reflectivity, f.reflectivity);
    code += check(rt_f.transmissivity, f.transmissivity);
    code += check(rt_f.slope_error, f.slope_error);
    code += check(rt_f.specularity_error, f.specularity_error);
    
    // check back struct
    code += check(rt_b.error_distribution_type, b.error_distribution_type);
    code += check(rt_b.reflectivity, b.reflectivity);
    code += check(rt_b.transmissivity, b.transmissivity);
    code += check(rt_b.slope_error, b.slope_error);
    code += check(rt_b.specularity_error, b.specularity_error);

    // expect += st_return_code::DATA_VALUE_NOT_FOUND
    code += st_get_optical_properties_set(pcxt, id + 1, &rt_set, &rt_f, &rt_b);
    
    return code;
}

st_return_t call_stapi_v2_remove_optics(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    SimulationData *data = cxt->p_data;
    args_optical_properties_set set = {"test1", 1.1, 1.1, 0};
    args_optical_properties_face f = {.5, .5, .5, .5, 'g'};
    args_optical_properties_face b = {.25, .25, .25, .25, 'g'};

    // check no optics set
    uint_fast64_t num = -1;
    uint_fast64_t id = -1;
    st_num_optics(pcxt, &num);
    st_return_t code = check(num, 0);

    // expect += st_return_code::SUCCESS
    code += st_add_optical_properties_set(pcxt, &set, &f, &b, &id);
    code += check(id, 0);
    st_num_optics(pcxt, &num);
    code += check(num, 1);
    
    // expect += st_return_code::SUCCESS
    code += st_clear_optics(pcxt);
    st_num_optics(pcxt, &num);
    code += check(num, 0);

    // expect += st_return_code::SUCCESS
    code += st_add_optical_properties_set(pcxt, &set, &f, &b, &id);
    code += check(id, 1);
    st_num_optics(pcxt, &num);
    code += check(num, 1);
    
    // expect += st_return_code::SUCCESS
    code += st_delete_optic(pcxt, 1);
    st_num_optics(pcxt, &num);
    code += check(num, 0);
    
    return code;
}

// functions to add/remove/modify elements
st_return_t call_stapi_v2_add_elements(st_context_v2_t pcxt)
{
    // set up dummy optical properties
    std::string dummy("dummy");
    OpticalPropertySet opt(InteractionType::REFLECTION, dummy);
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    SimulationData *data = cxt->p_data;
    OpticalPropertySetReference res = data->add_optical_property_set(opt);

    // check no elements set
    uint_fast64_t num = -1;
    uint_fast64_t id = -1;
    st_num_elements(pcxt, &num);
    st_return_t code = check(num, 0);

    // set element to non default values
    args_element el_args = { 2, 2, 2, 2, 2, 2, 2, false, true };
    args_geometry geom_args = { { 2 }, { 2, 2 }, 'c', 'p' };

    // test good element
    // expect += st_return_code::SUCCESS
    code += st_add_element(pcxt, &el_args, res.id, &geom_args, &id);
    code += check(id, 1);
    st_num_elements(pcxt, &num);
    code += check(num, 1);

    // test values set
    element_ptr el = data->get_element(id);
    auto sel = std::static_pointer_cast<SingleElement>(el);
    code += check(sel->is_enabled(), false);
    code += check(sel->is_virtual(), true);
    auto origin = sel->get_origin_ref();
    auto aim = sel->get_aim_vector_ref();
    auto zrot = sel->get_zrot();
    code += check(origin[0], 2);
    code += check(origin[1], 2);
    code += check(origin[2], 2);
    code += check(aim[0], 2);
    code += check(aim[1], 2);
    code += check(aim[2], 2);
    code += check(zrot, 2);

    optical_set_ptr opt_set = sel->get_optical_property_set();
    code += check(opt_set->get_name(), dummy);

    aperture_ptr ap = sel->get_aperture();
    code += check(ap->my_type, ApertureType::CIRCLE);
    auto circle = std::static_pointer_cast<Circle>(ap);
    code += check(circle->diameter, 2);

    surface_ptr surf = sel->get_surface();
    code += check(surf->my_type, SurfaceType::PARABOLA);
    auto parabola = std::static_pointer_cast<Parabola>(surf);
    code += check(parabola->focal_length_x, 1. / 4.);
    code += check(parabola->focal_length_y, 1. / 4.);

    // try bad optical id
    // expect code += st_return_code::DATA_VALUE_NOT_FOUND
    code += st_add_element(pcxt, &el_args, res.id + 1, &geom_args, &id);
    code += check(id, 1);
    st_num_elements(pcxt, &num);
    code += check(num, 1);

    // try bad aperture type
    // expect code += st_return_code::INVALID_ARGUMENTS
    args_geometry geom_bad_ap = { { 2 }, { 2, 2 }, 'z', 'p' };
    code += st_add_element(pcxt, &el_args, res.id, &geom_bad_ap, &id);
    code += check(id, 1);
    st_num_elements(pcxt, &num);
    code += check(num, 1);
    
    // try bad surface type
    // expect code += st_return_code::INVALID_ARGUMENTS
    args_geometry geom_bad_surf = { { 2 }, { 2, 2 }, 'c', 'z' };
    code += st_add_element(pcxt, &el_args, res.id, &geom_bad_surf, &id);
    code += check(id, 1);
    st_num_elements(pcxt, &num);
    code += check(num, 1);
    
    // try bad aperture params
    // expect code += st_return_code::INVALID_ARGUMENTS
    args_geometry geom_bad_ap_params = { { -2 }, { 2, 2 }, 'c', 'p' };
    code += st_add_element(pcxt, &el_args, res.id, &geom_bad_ap_params, &id);
    code += check(id, 1);
    st_num_elements(pcxt, &num);
    code += check(num, 1);

    // try bad surface params
    // expect code += st_return_code::INVALID_ARGUMENTS
    args_geometry geom_bad_surf_params = { { 2 }, { std::numeric_limits<double>::quiet_NaN(), 2 }, 'c', 'p'};
    code += st_add_element(pcxt, &el_args, res.id, &geom_bad_surf_params, &id);
    code += check(id, 1);
    st_num_elements(pcxt, &num);
    code += check(num, 1);

    // test another good element
    // expect += st_return_code::SUCCESS
    code += st_add_element(pcxt, &el_args, res.id, &geom_args, &id);
    code += check(id, 2);
    st_num_elements(pcxt, &num);
    code += check(num, 2);

    // expect == st_return_code::DATA_VALUE_NOT_FOUND
    //         + 4 * st_return_code::INVALID_ARGUMENTS 
    return code;
}

st_return_t add_dummy_element(st_context_v2_t             pcxt,
                              SimulationData              *data,
                              args_element                &el_args,
                              OpticalPropertySetReference &res,
                              args_geometry               &geom_args,
                              uint_fast64_t               &id,
                              int32_t                     group = -1)
{
    // set up dummy optical properties
    OpticalPropertySet opt(InteractionType::REFLECTION, std::string("dummy"));
    res = data->add_optical_property_set(opt);
    // set element to non default values
    el_args = { 2, 2, 2, 2, 2, 2, 2, false, true };
    el_args.group = group;
    geom_args = { { 2 }, { 2, 2 }, 'c', 'p' };
    return st_add_element(pcxt, &el_args, res.id, &geom_args, &id);
}

st_return_t call_stapi_v2_get_element(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    SimulationData *data = cxt->p_data;
    args_element  el_args;
    OpticalPropertySetReference res;
    args_geometry geom_args;
    uint_fast64_t id = -1;
    st_return_t code = add_dummy_element(pcxt, data, el_args, res, geom_args, id);
    
    args_element  rt_el_args;
    int_fast64_t  opt_id;
    args_geometry rt_geom_args;
    
    // expect += st_return_code::SUCCESS
    code += st_get_element(pcxt, id, &rt_el_args, &opt_id, &rt_geom_args);

    // test values set
    element_ptr el = data->get_element(1);
    auto sel = std::static_pointer_cast<SingleElement>(el);
    code += check(rt_el_args.enabled_flag, sel->is_enabled());
    code += check(rt_el_args.virtual_flag, sel->is_virtual());
    auto origin = sel->get_origin_ref();
    auto aim = sel->get_aim_vector_ref();
    auto zrot = sel->get_zrot();
    code += check(rt_el_args.x, origin[0]);
    code += check(rt_el_args.y, origin[1]);
    code += check(rt_el_args.z, origin[2]);
    code += check(rt_el_args.ax, aim[0]);
    code += check(rt_el_args.ay, aim[1]);
    code += check(rt_el_args.az, aim[2]);
    code += check(rt_el_args.zrot, zrot);

    code += check(opt_id, res.id);

    aperture_ptr ap = sel->get_aperture();
    code += check(rt_geom_args.ap, aperture_to_char(ap->get_type()));
    auto circle = std::static_pointer_cast<Circle>(ap);
    code += check(rt_geom_args.a_params[0], circle->diameter);

    surface_ptr surf = sel->get_surface();
    code += check(rt_geom_args.surf, surface_to_char(surf->get_type()));
    auto parabola = std::static_pointer_cast<Parabola>(surf);
    code += check(rt_geom_args.s_params[0], parabola->focal_length_x);
    code += check(rt_geom_args.s_params[1], parabola->focal_length_y);

    // expect += st_return_code::DATA_VALUE_NOT_FOUND
    code += st_get_element(pcxt, id + 1, &rt_el_args, &opt_id, &rt_geom_args);
    
    return code;
}

st_return_t call_stapi_v2_remove_elements(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    SimulationData *data = cxt->p_data;
    args_element  el_args;
    OpticalPropertySetReference res;
    args_geometry geom_args;
    uint_fast64_t id = -1;
    st_return_t code = add_dummy_element(pcxt, data, el_args, res, geom_args, id);
    
    uint_fast64_t num = -1;
    st_num_elements(pcxt, &num);
    code += check(num, 1);
    
    // expect += st_return_code::SUCCESS
    code += st_clear_elements(pcxt);
    st_num_elements(pcxt, &num);
    code += check(num, 0);
    
    // expect += st_return_code::SUCCESS
    code += st_add_element(pcxt, &el_args, res.id, &geom_args, &id);
    code += check(id, 2);
    st_num_elements(pcxt, &num);
    code += check(num, 1);

    // expect += st_return_code::WARNING_NOT_FOUND
    code += st_delete_element(pcxt, 1);
    st_num_elements(pcxt, &num);
    code += check(num, 1);
    
    // expect += st_return_code::SUCCESS
    code += st_delete_element(pcxt, 2);
    st_num_elements(pcxt, &num);
    code += check(num, 0);
    
    // expect == st_return_code::WARNING_NOT_FOUND
    return code;
}

st_return_t call_stapi_v2_toggle_element(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    SimulationData *data = cxt->p_data;
    args_element  el_args;
    OpticalPropertySetReference res;
    args_geometry geom_args;
    uint_fast64_t id = -1;
    st_return_t code = add_dummy_element(pcxt, data, el_args, res, geom_args, id);
    
    // expect += st_return_code::SUCCESS
    code += st_element_enabled(pcxt, 1, true);
    
    // expect += st_return_code::SUCCESS
    code += st_element_virtual(pcxt, 1, false);
    
    // expect += st_return_code::WARNING_NOT_FOUND
    code += st_element_enabled(pcxt, 0, false);
    
    // expect += st_return_code::WARNING_NOT_FOUND
    code += st_element_virtual(pcxt, 0, true);

    // test values set
    auto sel = std::static_pointer_cast<SingleElement>(data->get_element(1));
    code += check(sel->is_enabled(), true);
    code += check(sel->is_virtual(), false);

    // expect == 2 * st_return_code::WARNING_NOT_FOUND
    return code; 
}

st_return_t call_stapi_v2_element_xyz(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    SimulationData *data = cxt->p_data;
    args_element  el_args;
    OpticalPropertySetReference res;
    args_geometry geom_args;
    uint_fast64_t id = -1;
    st_return_t code = add_dummy_element(pcxt, data, el_args, res, geom_args, id);
    
    // expect += st_return_code::SUCCESS
    code += st_element_xyz(pcxt, 1, 1, 1, 1);
    
    // expect += st_return_code::WARNING_NOT_FOUND
    code += st_element_xyz(pcxt, 0, 2, 2, 2);

    // test values set
    auto sel = std::static_pointer_cast<SingleElement>(data->get_element(1));
    auto origin = sel->get_origin_ref();
    code += check(origin[0], 1);
    code += check(origin[1], 1);
    code += check(origin[2], 1);

    // expect == st_return_code::WARNING_NOT_FOUND
    return code; 
}

st_return_t call_stapi_v2_element_aim(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    SimulationData *data = cxt->p_data;
    args_element  el_args;
    OpticalPropertySetReference res;
    args_geometry geom_args;
    uint_fast64_t id = -1;
    st_return_t code = add_dummy_element(pcxt, data, el_args, res, geom_args, id);
    
    // expect += st_return_code::SUCCESS
    code += st_element_aim(pcxt, 1, 1, 1, 1);
    
    // expect += st_return_code::WARNING_NOT_FOUND
    code += st_element_aim(pcxt, 0, 2, 2, 2);

    // test values set
    auto sel = std::static_pointer_cast<SingleElement>(data->get_element(1));
    auto aim = sel->get_aim_vector_ref();
    code += check(aim[0], 1);
    code += check(aim[1], 1);
    code += check(aim[2], 1);

    // expect == st_return_code::WARNING_NOT_FOUND
    return code; 
}

st_return_t call_stapi_v2_element_zrot(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    SimulationData *data = cxt->p_data;
    args_element  el_args;
    OpticalPropertySetReference res;
    args_geometry geom_args;
    uint_fast64_t id = -1;
    st_return_t code = add_dummy_element(pcxt, data, el_args, res, geom_args, id);
    
    // expect += st_return_code::SUCCESS
    code += st_element_zrot(pcxt, 1, 1);
    
    // expect += st_return_code::WARNING_NOT_FOUND
    code += st_element_zrot(pcxt, 0, 2);

    // test values set
    auto sel = std::static_pointer_cast<SingleElement>(data->get_element(1));
    auto zrot = sel->get_zrot();
    code += check(zrot, 1);

    // expect == st_return_code::WARNING_NOT_FOUND
    return code; 
}

st_return_t call_stapi_v2_element_aperture(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    SimulationData *data = cxt->p_data;
    args_element  el_args;
    OpticalPropertySetReference res;
    args_geometry geom_args;
    uint_fast64_t id = -1;
    st_return_t code = add_dummy_element(pcxt, data, el_args, res, geom_args, id);

    double params[2] = { 4, 4 };

    // expect += st_return_code::SUCCESS
    code += st_element_aperture(pcxt, 1, 'r', params);

    auto sel = std::static_pointer_cast<SingleElement>(data->get_element(1));
    aperture_ptr ap = sel->get_aperture();
    code += check(ap->my_type, ApertureType::RECTANGLE);
    auto rectangle = std::static_pointer_cast<Rectangle>(ap);
    code += check(rectangle->x_length(), 4);
    code += check(rectangle->y_length(), 4);

    // expect += st_return_code::INVALID_ARGUMENTS
    code += st_element_aperture(pcxt, 1, 'z', params);

    double bad_params[2] = { -4, 4 };
    // expect += st_return_code::INVALID_ARGUMENTS
    code += st_element_aperture(pcxt, 1, 'r', bad_params);
    
    // expect == 2 * st_return_code::INVALID_ARGUMENTS
    return code;
}

st_return_t call_stapi_v2_element_surface(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    SimulationData *data = cxt->p_data;
    args_element  el_args;
    OpticalPropertySetReference res;
    args_geometry geom_args;
    uint_fast64_t id = -1;
    st_return_t code = add_dummy_element(pcxt, data, el_args, res, geom_args, id);
    
    double params[1] = { 4 };
    // expect += st_return_code::SUCCESS
    code += st_element_surface(pcxt, 1, 's', params);
    
    auto sel = std::static_pointer_cast<SingleElement>(data->get_element(1));
    surface_ptr surf = sel->get_surface();
    code += check(surf->my_type, SurfaceType::SPHERE);
    auto sphere = std::static_pointer_cast<Sphere>(surf);
    code += check(sphere->vertex_curv, 4);

    // expect += st_return_code::INVALID_ARGUMENTS
    code += st_element_surface(pcxt, 1, 'z', params);

    double bad_params[1] = { -4 };
    // expect += st_return_code::INVALID_ARGUMENTS
    code += st_element_surface(pcxt, 1, 's', bad_params);

    // expect == 2 * st_return_code::INVALID_ARGUMENTS
    return code;
}

st_return_t call_stapi_v2_element_geometry(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    SimulationData *data = cxt->p_data;
    args_element  el_args;
    OpticalPropertySetReference res;
    args_geometry geom_args;
    uint_fast64_t id = -1;
    st_return_t code = add_dummy_element(pcxt, data, el_args, res, geom_args, id);

    // test other good geometry
    args_geometry geom_other = { { 4, 4 }, { 4 }, 'r', 's'};
    code += st_element_geometry(pcxt, id, &geom_other);

    // try bad aperture type
    // expect code += st_return_code::INVALID_ARGUMENTS
    args_geometry geom_bad_ap = { { 2 }, { 2, 2 }, 'z', 'p' };
    code += st_element_geometry(pcxt, id, &geom_bad_ap);

    // try bad surface type
    // expect code += st_return_code::INVALID_ARGUMENTS
    args_geometry geom_bad_surf = { { 2 }, { 2, 2 }, 'c', 'z' };
    code += st_element_geometry(pcxt, id, &geom_bad_surf);

    // try bad aperture params
    // expect code += st_return_code::INVALID_ARGUMENTS
    args_geometry geom_bad_ap_params = { { -2 }, { 2, 2 }, 'c', 'p' };
    code += st_element_geometry(pcxt, id, &geom_bad_ap_params);

    // try bad surface params
    // expect code += st_return_code::INVALID_ARGUMENTS
    args_geometry geom_bad_surf_params = { { 2 }, { std::numeric_limits<double>::quiet_NaN(), 2 }, 'c', 'p'};
    code += st_element_geometry(pcxt, id, &geom_bad_surf_params);

    // expect == 4 * st_return_code::INVALID_ARGUMENTS 
    return code;
}

st_return_t call_stapi_v2_element_optic(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    SimulationData *data = cxt->p_data;
    args_element  el_args;
    OpticalPropertySetReference res;
    args_geometry geom_args;
    uint_fast64_t id = -1;
    st_return_t code = add_dummy_element(pcxt, data, el_args, res, geom_args, id);

    // expect += st_return_code::DATA_VALUE_NOT_FOUND
    code += st_element_optic(pcxt, 1, res.id + 1);

    // set up other optical properties
    OpticalPropertySet opt2(InteractionType::REFLECTION, std::string("other"));
    OpticalPropertySetReference res2 = data->add_optical_property_set(opt2);
    // expect += st_return_code::SUCCESS
    code += st_element_optic(pcxt, 1, res2.id);
    
    auto sel = std::static_pointer_cast<SingleElement>(data->get_element(1));
    optical_set_ptr opt_set = sel->get_optical_property_set();
    code += check(opt_set->get_name(), std::string("other"));
    
    // expect == st_return_code::DATA_VALUE_NOT_FOUND
    return code;
}

st_return_t call_stapi_v2_element_group(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    SimulationData *data = cxt->p_data;
    args_element  el_args;
    OpticalPropertySetReference res;
    args_geometry geom_args;
    uint_fast64_t id = -1;
    // expect += st_return_code::WARNING_GROUP_IGNORED
    st_return_t code = add_dummy_element(pcxt, data, el_args, res, geom_args, id, -2);

    args_element rt_el_args;
    int_fast64_t opt_id;
    args_geometry rt_geom_args;
    
    // expect += st_return_code::SUCCESS
    code += st_get_element(pcxt, id, &rt_el_args, &opt_id, &rt_geom_args);
    code += check(rt_el_args.group, -1);
    
    // expect += st_return_code::WARNING_GROUP_IGNORED
    code += st_element_group(pcxt, id, -2);
    code += st_get_element(pcxt, id, &rt_el_args, &opt_id, &rt_geom_args);
    code += check(rt_el_args.group, -1);
    
    // expect += st_return_code::SUCCESS
    code += st_element_group(pcxt, id, 1);
    code += st_get_element(pcxt, id, &rt_el_args, &opt_id, &rt_geom_args);
    code += check(rt_el_args.group, 1);
    
    // expect == 2 * st_return_code::WARNING_GROUP_IGNORED
    return code;
}

// sun functions
st_return_t call_stapi_v2_add_sun(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    // expect == 0
    st_return_t code = (st_return_t)cxt->p_data->get_number_of_ray_sources();

    double good_angles[3]      = {0, 1, 2};
    double good_intensities[3] = {0, 1, 2};
    double bad_intensities[2]  = {0, -1};

    args_sun args = {3, 608, 303, 1000, 5, ' '};

    // test bad intensitites
    // expect += st_return_code::EXCEPTION
    code += st_add_sun(pcxt, &args, good_angles, bad_intensities);

    auto sun_0 = std::static_pointer_cast<Sun>(cxt->p_data->get_ray_source(0));
    auto xyz = sun_0->get_position();
    // expect position to be what was set
    if (xyz[0] != 608 || xyz[1] != 303 || xyz[2] != 1000) ++code;
    std::vector<double> angles, intensities;
    sun_0->get_user_data(angles, intensities);
    // expect vectors of no size
    if (angles.size() || intensities.size()) ++code;

    // test good intensities
    // expect += st_return_code::SUCCESS
    code += st_add_sun(pcxt, &args, good_angles, good_intensities);
    
    sun_0->get_user_data(angles, intensities);
    if (angles[0] != good_angles[0]
        || angles[1] != good_angles[1]
        || angles[2] != good_angles[2]
        || intensities[0] != good_intensities[0]
        || intensities[1] != good_intensities[1]
        || intensities[2] != good_intensities[2]) ++code;

    // test built in shapes
    args.npoints = 0;
    
    // expect += st_return_code::WARNING_SUN_SHAPE_IGNORED
    code += st_add_sun(pcxt, &args, good_angles, good_intensities);
    auto sun_1 = std::static_pointer_cast<Sun>(cxt->p_data->get_ray_source(0));
    SunShape shape = sun_1->get_shape();
    double sigma = sun_1->get_sigma();
    // expect both of these to equal defaults
    if (shape != SunShape::GAUSSIAN || sigma != 4.65) ++code;
    // expect sun object to be same
    if (sun_0 != sun_1) ++code;
    // expect one ray source
    code += (st_return_t)cxt->p_data->get_number_of_ray_sources() - 1;

    // expect += st_return_code::SUCCESS
    args.shape = 'g';
    code += st_add_sun(pcxt, &args, good_angles, good_intensities);
    auto sun_2 = std::static_pointer_cast<Sun>(cxt->p_data->get_ray_source(0));
    shape = sun_2->get_shape();
    sigma = sun_2->get_sigma();
    // expect both of these to equal args
    if (shape != SunShape::GAUSSIAN || sigma != 5) ++code;
    // expect sun object to be same
    if (sun_1 != sun_2) ++code;
    // expect one ray source
    code += (st_return_t)cxt->p_data->get_number_of_ray_sources() - 1;

    // expect += st_return_code::SUCCESS
    args.shape = 'p';
    code += st_add_sun(pcxt, &args, good_angles, good_intensities);
    auto sun_3 = std::static_pointer_cast<Sun>(cxt->p_data->get_ray_source(0));
    shape = sun_3->get_shape();
    double hw = sun_3->get_half_width();
    // expect both of these to equal args
    if (shape != SunShape::PILLBOX || hw != 5) ++code;
    // expect sun object to be same
    if (sun_2 != sun_3) ++code;
    // expect one ray source
    code += (st_return_t)cxt->p_data->get_number_of_ray_sources() - 1;

    // expect += st_return_code::WARNING_SUN_SHAPE_IGNORED
    args.shape = 'l';
    code += st_add_sun(pcxt, &args, good_angles, good_intensities);
    auto sun_4 = std::static_pointer_cast<Sun>(cxt->p_data->get_ray_source(0));
    shape = sun_4->get_shape();
    sigma = sun_4->get_sigma();
    // expect both of these to equal defaults
    if (shape != SunShape::GAUSSIAN || sigma != 4.65) ++code;
    // expect sun object to be same
    if (sun_3 != sun_4) ++code;
    // expect one ray source
    code += (st_return_t)cxt->p_data->get_number_of_ray_sources() - 1;

    // expect += 0
    args.shape = 'b';
    args.sigma_halfwidth_csr = .5;
    code += st_add_sun(pcxt, &args, good_angles, good_intensities);
    auto sun_5 = std::static_pointer_cast<Sun>(cxt->p_data->get_ray_source(0));
    shape = sun_5->get_shape();
    double csr = sun_5->get_circumsolar_ratio();
    // expect both of these to equal args
    if (shape != SunShape::BUIE_CSR || csr != .5) ++code;
    // expect sun object to be same
    if (sun_4 != sun_5) ++code;
    // expect one ray source
    code += (st_return_t)cxt->p_data->get_number_of_ray_sources() - 1;

    // expect += st_return_code::WARNING_SUN_SHAPE_IGNORED
    args.shape = 'u';
    code += st_add_sun(pcxt, &args, good_angles, good_intensities);
    auto sun_6 = std::static_pointer_cast<Sun>(cxt->p_data->get_ray_source(0));
    shape = sun_6->get_shape();
    sigma = sun_6->get_sigma();
    // expect both of these to equal defaults
    if (shape != SunShape::GAUSSIAN || sigma != 4.65) ++code;
    // expect sun object to be same
    if (sun_5 != sun_6) ++code;
    // expect one ray source
    code += (st_return_t)cxt->p_data->get_number_of_ray_sources() - 1;

    // expect == 3 * st_return_code::WARNING_SUN_SHAPE_IGNORED
    //           + st_return_code::EXCEPTION
    return code;
}

st_return_t call_stapi_v2_get_sun(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);

    args_sun rt_args;
    double *rt_angle     = new double[3];
    double *rt_intensity = new double[3];

    // expect += st_return_code::DATA_VALUE_NOT_FOUND
    st_return_t code = st_get_sun(pcxt, &rt_args, &rt_angle, &rt_intensity);

    args_sun args = {0, 608, 303, 1000, 5, 'g'};
    // expect += st_return_code::SUCCESS
    code += st_add_sun(pcxt, &args, {}, {});
    
    // expect += st_return_code::SUCCESS
    code += st_get_sun(pcxt, &rt_args, &rt_angle, &rt_intensity);

    auto sun = cxt->p_data->get_ray_source(0);

    code += check(rt_args.npoints, 0);
    code += check(rt_args.shape, sunshape_to_char(sun->get_shape()));
    code += check(rt_args.sigma_halfwidth_csr, sun->get_sigma());
    auto pos = sun->get_position();
    code += check(rt_args.x, pos[0]);
    code += check(rt_args.y, pos[1]);
    code += check(rt_args.z, pos[2]);

    double angles[3]      = {0, 1, 2};
    double intensities[3] = {0, 1, 2};
    args = {3, 608, 303, 1000, 5, 'd'};
    // expect += st_return_code::SUCCESS
    code += st_add_sun(pcxt, &args, angles, intensities);

    // expect += st_return_code::SUCCESS
    code += st_get_sun(pcxt, &rt_args, &rt_angle, &rt_intensity, 3);
    code += check(rt_angle[0], angles[0]);
    code += check(rt_angle[1], angles[1]);
    code += check(rt_angle[2], angles[2]);
    code += check(rt_intensity[0], intensities[0]);
    code += check(rt_intensity[1], intensities[1]);
    code += check(rt_intensity[2], intensities[2]);
    
    // expect += st_return_code::INVALID_ARGUMENTS
    code += st_get_sun(pcxt, &rt_args, &rt_angle, &rt_intensity);
    // expect += st_return_code::INVALID_ARGUMENTS
    code += st_get_sun(pcxt, &rt_args, &rt_angle, &rt_intensity, 4);

    // expect == st_return_code::DATA_VALUE_NOT_FOUND +
    //           2 * st_return_code::INVALID_ARGUMENTS
    return code;
}

st_return_t call_stapi_v2_sun_shape(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    // expect == 0
    st_return_t code = (st_return_t)cxt->p_data->get_number_of_ray_sources();

    double good_angles[3]      = {0, 1, 2};
    double good_intensities[3] = {0, 1, 2};

    args_sun args = {3, 608, 303, 1000, 5, 'g'};

    // expect += st_return_code::SUCCESS
    code += st_add_sun(pcxt, &args, good_angles, good_intensities);

    // expect += st_return_code::WARNING_SUN_SHAPE_IGNORED
    code += st_sun_shape(pcxt, ' ', 0);
    auto sun = std::static_pointer_cast<Sun>(cxt->p_data->get_ray_source(0));
    SunShape shape = sun->get_shape();
    double sigma = sun->get_sigma();
    // expect both of these to equal defaults
    if (shape != SunShape::GAUSSIAN || sigma != 4.65) ++code;

    // tests below follow switch statement order in st_sun_shape

    // expect += 0
    code += st_sun_shape(pcxt, 'g', 5);
    auto sun_2 = std::static_pointer_cast<Sun>(cxt->p_data->get_ray_source(0));
    shape = sun_2->get_shape();
    sigma = sun_2->get_sigma();
    // expect both of these to equal args
    if (shape != SunShape::GAUSSIAN || sigma != 5) ++code;
    
    // expect += 0
    code += st_sun_shape(pcxt, 'p', 5);
    auto sun_3 = std::static_pointer_cast<Sun>(cxt->p_data->get_ray_source(0));
    shape = sun_3->get_shape();
    double hw = sun_3->get_half_width();
    // expect both of these to equal args
    if (shape != SunShape::PILLBOX || hw != 5) ++code;
    
    // expect += st_return_code::WARNING_SUN_SHAPE_IGNORED
    code += st_sun_shape(pcxt, 'l', 5);
    auto sun_4 = std::static_pointer_cast<Sun>(cxt->p_data->get_ray_source(0));
    shape = sun_4->get_shape();
    sigma = sun_4->get_sigma();
    // expect both of these to equal defaults
    if (shape != SunShape::GAUSSIAN || sigma != 4.65) ++code;
    
    // expect += 0
    code += st_sun_shape(pcxt, 'b', .5);
    auto sun_5 = std::static_pointer_cast<Sun>(cxt->p_data->get_ray_source(0));
    shape = sun_5->get_shape();
    double csr = sun_5->get_circumsolar_ratio();
    // expect both of these to equal args
    if (shape != SunShape::BUIE_CSR || csr != .5) ++code;
    
    // expect += st_return_code::WARNING_SUN_SHAPE_IGNORED
    code += st_sun_shape(pcxt, 'u', 0);
    auto sun_6 = std::static_pointer_cast<Sun>(cxt->p_data->get_ray_source(0));
    shape = sun_6->get_shape();
    sigma = sun_6->get_sigma();
    // expect both of these to equal defaults
    if (shape != SunShape::GAUSSIAN || sigma != 4.65) ++code;
    
    // expect code == 3 * st_return_code::WARNING_SUN_SHAPE_IGNORED
    return code;
}

st_return_t call_stapi_v2_sun_xyz(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    // expect == 0
    st_return_t code = (st_return_t)cxt->p_data->get_number_of_ray_sources();

    double good_angles[3]      = {0, 1, 2};
    double good_intensities[3] = {0, 1, 2};

    args_sun args = {3, 608, 303, 1000, 5, 'g'};

    // expect += st_return_code::SUCCESS
    code += st_add_sun(pcxt, &args, good_angles, good_intensities);

    // expect == 0
    code += st_sun_xyz(pcxt, 0, 0, 0);
    auto sun = std::static_pointer_cast<Sun>(cxt->p_data->get_ray_source(0));
    auto xyz = sun->get_position();
    // expect position to be what was set
    if (xyz[0] != 0 || xyz[1] != 0 || xyz[2] != 0) ++code;

    return code;
}

st_return_t call_stapi_v2_sun_userdata(st_context_v2_t pcxt)
{
    st_context *cxt = reinterpret_cast<st_context*>(pcxt);
    // expect == 0
    st_return_t code = (st_return_t)cxt->p_data->get_number_of_ray_sources();

    double good_angles[3]      = {0, 1, 2};
    double good_intensities[3] = {0, 1, 2};

    args_sun args = {3, 608, 303, 1000, 5, 'g'};

    // expect += st_return_code::SUCCESS
    code += st_add_sun(pcxt, &args, good_angles, good_intensities);

    double new_angles[3]      = {0, .1, .2};
    double new_intensities[3] = {0, .1, .2};
    // expect += 0
    code += st_sun_userdata(pcxt, 3, new_angles, new_intensities);
    auto sun = std::static_pointer_cast<Sun>(cxt->p_data->get_ray_source(0));
    
    std::vector<double> angles, intensities;
    sun->get_user_data(angles, intensities);
    if (angles[0] != new_angles[0]
        || angles[1] != new_angles[1]
        || angles[2] != new_angles[2]
        || intensities[0] != new_intensities[0]
        || intensities[1] != new_intensities[1]
        || intensities[2] != new_intensities[2]) ++code;

    double bad_intensities[2] = {0, -1};
    // expect += st_return_code::EXCEPTION
    code += st_sun_userdata(pcxt, 3, good_angles, bad_intensities);
    
    return code;
}

typedef struct expected_solar_values {
    double azimuth;
    double elevation;
    double sun_x;
    double sun_y;
    double sun_z;
} expected_solar_values;

st_return_t test_calculator_method(st_context_v2_t pcxt,
                                   SolTrace::Data::SolarPositionCalculationMethod calc,
                                   expected_solar_values *expected)
{
    args_sun_location loc = { 40.0, -105.0, -7.0 };
    args_sun_datetime dt  = { 2025, 6, 20, 12 };

    double azimuth, elevation, zenith, sun_x, sun_y, sun_z;
    
    // azimuth/elevation
    st_return_t code = st_get_sun_az_el(pcxt, calc, &loc, &dt, &azimuth, &elevation);
    
    // if checking inputs skip the rest
    if (code > st_return_code::SUCCESS) return code;

    code += check_near(azimuth, expected->azimuth, 1e-3);
    code += check_near(elevation, expected->elevation, 1e-3);

    // azimuth/elevation
    st_get_sun_az_zen(pcxt, calc, &loc, &dt, &azimuth, &zenith);
    
    double expected_zenith = 90.0 - expected->elevation;
    code += check_near(azimuth, expected->azimuth, 1e-3);
    code += check_near(zenith, expected_zenith, 1e-3);

    // sun vector
    st_get_sun_vector(pcxt, calc, &loc, &dt, &sun_x, &sun_y, &sun_z);

    // on macOS sun_x diff is 0.000003 -> changed bound from 1e-6 to 5e-6
    code += check_near(sun_x, expected->sun_x, 5e-6);
    code += check_near(sun_y, expected->sun_y, 5e-6);
    code += check_near(sun_z, expected->sun_z, 5e-6);

    return code;
}

st_return_t call_stapi_v2_solar_calculator(st_context_v2_t pcxt)
{
    // st_return_t code = st_return_code::SUCCESS;

    // method: SPA
    // Expected values - confirmed with https://gml.noaa.gov/grad/solcalc/azel.html
    expected_solar_values expected_SPA = {178.61128380, 73.439035265, 0.006908, -0.2849516, 0.9585169};
    st_return_t code = test_calculator_method(pcxt,
                                              SolTrace::Data::SolarPositionCalculationMethod::SPA,
                                              &expected_SPA);
    
    // method: SPA_ORIGINAL
    expected_solar_values expected_SPA_ORIGINAL = {178.614529, 73.443097, 0.006890135, -0.284884148, 0.95853719};
    code += test_calculator_method(pcxt,
                                   SolTrace::Data::SolarPositionCalculationMethod::SPA_ORIGINAL,
                                   &expected_SPA_ORIGINAL);
   
    // method: SOLPOS
    expected_solar_values expected_SOLPOS = {178.6119537, 73.430931, 0.006907870, -0.28508729, 0.95847666};
    code += test_calculator_method(pcxt,
                                   SolTrace::Data::SolarPositionCalculationMethod::SOLPOS,
                                   &expected_SOLPOS);
   
    // method: LEGACY
    expected_solar_values expected_LEGACY = {179.9991897, 73.435378, 4.031894e-06, -0.28509658, 0.958498794};
    code += test_calculator_method(pcxt,
                                   SolTrace::Data::SolarPositionCalculationMethod::LEGACY,
                                   &expected_LEGACY);
   
    // method: DUFFIE
    expected_solar_values expected_DUFFIE = {179.110386, 73.43996784, 0.00442455, -0.284985454, 0.958521629};
    code += test_calculator_method(pcxt,
                                   SolTrace::Data::SolarPositionCalculationMethod::DUFFIE,
                                   &expected_DUFFIE);
   
    // expect += st_return_code::INVALID_ARGUMENTS
    code += test_calculator_method(pcxt,
                                   SolTrace::Data::SolarPositionCalculationMethod::CALCULATOR_COUNT,
                                   &expected_DUFFIE);
    return code;
}

//////////////////////////////////
// Simlulation Runner Functions //
//////////////////////////////////

// functions for SolTrace runner management
st_return_t call_stapi_v2_sim_setup(st_context_v2_t pcxt)
{
    // this fails on first execution after build if optix is built
    // when debugging first execution, it passes 
    json root = load_json();

    st_return_t code = st_read_input_json(pcxt, root.dump().c_str());

    if (code != st_return_code::SUCCESS) return code;

    code =  st_sim_setup(pcxt, st_runner_type_t::NATIVE);
    unsigned int *seeds_test = new unsigned int[2] { 608, 303, };
    code += st_sim_setup(pcxt, st_runner_type_t::NATIVE, 1, seeds_test, 2);
    code += st_sim_setup(pcxt, st_runner_type_t::EMBREE);
    code += st_sim_setup(pcxt, st_runner_type_t::OPTIX);

    return code;
}

st_return_t call_stapi_v2_sim_run_v2(st_context_v2_t pcxt, st_runner_type_t runner_type)
{
    if (runner_type != st_runner_type_t::OPTIX)
    {
        st_context *cxt = reinterpret_cast<st_context*>(pcxt);
        cxt->p_data->set_number_of_rays(1000);
    }

    return st_sim_setup(pcxt, runner_type) + st_sim_run_v2(pcxt);
}

st_return_t call_stapi_v2_sim_report(st_context_v2_t pcxt, st_runner_type_t runner_type)
{
    if (runner_type != st_runner_type_t::OPTIX)
    {
        st_context *cxt = reinterpret_cast<st_context*>(pcxt);
        cxt->p_data->set_number_of_rays(1000);
    }

    st_return_t code = st_sim_setup(pcxt, runner_type);
    // expect += st_return_code::RUNNER_NOT_READY_TO_REPORT
    code += st_sim_report(pcxt, RunnerStatistics::RAY_RECORDS);

    code += st_sim_run_v2(pcxt);

    // expect += 2 * st_return_code::INVALID_ARGUMENTS
    code += st_sim_report(pcxt, RunnerStatistics::RAY_RECORDS - 1);
    code += st_sim_report(pcxt, RunnerStatistics::STATISTICS_COUNT);

    // expect += st_return_code::SUCCESS
    code += st_sim_report(pcxt, RunnerStatistics::GROUPED_COUNTS);

    return code;
}

///////////////////////////////////
// Simlulation Results Functions //
///////////////////////////////////

st_return_t call_stapi_v2_write_results_csv(st_context_v2_t  pcxt, 
                                            st_runner_type_t runner_type, 
                                            const char       *filename)
{
    if (runner_type != st_runner_type_t::OPTIX)
    {
        st_context* cxt = reinterpret_cast<st_context*>(pcxt);
        cxt->p_data->set_number_of_rays(1000);
    }

    st_return_t code = st_sim_setup(pcxt, runner_type); 
    code += st_sim_run_v2(pcxt);
    code += st_sim_report(pcxt, 0);
    code += st_write_results_csv(pcxt, filename);

    std::error_code ec;
    std::filesystem::remove(filename, ec);

    return code;
}

st_return_t call_stapi_v2_write_group_results_json(st_context_v2_t  pcxt, 
                                                   st_runner_type_t runner_type, 
                                                   const char       *filename)
{
    if (runner_type != st_runner_type_t::OPTIX)
    {
        st_context* cxt = reinterpret_cast<st_context*>(pcxt);
        cxt->p_data->set_number_of_rays(1000);
    }

    st_return_t code = st_sim_setup(pcxt, runner_type); 
    code += st_sim_run_v2(pcxt);
    code += st_sim_report(pcxt, 1);
    code += st_write_group_results_json(pcxt, filename);

    std::error_code ec;
    std::filesystem::remove(filename, ec);

    return code;
}

// functions to get results directly

st_return_t call_stapi_v2_locations(st_context_v2_t  pcxt,
                                    st_runner_type_t runner_type)
{
    st_context* cxt = reinterpret_cast<st_context*>(pcxt);
    cxt->p_data->set_number_of_rays(1000);

    st_return_t code = st_sim_setup(pcxt, runner_type); 
    code += st_sim_run_v2(pcxt);
    code += st_sim_report(pcxt, 0);
    
    uint_fast64_t num;
    
    code += st_num_intersections(pcxt, &num);

    double *loc_x = new double[num];
    double *loc_y = new double[num];
    double *loc_z = new double[num];

    code += st_locations(pcxt, loc_x, loc_y, loc_z);

    code += check_not(loc_x, nullptr);
    code += check_not(loc_y, nullptr);
    code += check_not(loc_z, nullptr);

    uint_fast64_t n = 0;
    glm::dvec3 loc;
    for (auto it = cxt->p_results->get_ray_record_iterator(); !cxt->p_results->is_at_end(it); ++it)
    {
        ray_record_ptr rec = *it;
        for (auto jt = rec->get_interaction_record_iterator(); !rec->is_at_end(jt); ++jt)
        {
            loc = (*jt)->location;
            code += check(loc_x[n], loc[0]);
            code += check(loc_y[n], loc[1]);
            code += check(loc_z[n], loc[2]);
            ++n;
        }
    }

    return code;
}

st_return_t call_stapi_v2_cosines(st_context_v2_t  pcxt,
                                  st_runner_type_t runner_type)
{
    st_context* cxt = reinterpret_cast<st_context*>(pcxt);
    cxt->p_data->set_number_of_rays(1000);

    st_return_t code = st_sim_setup(pcxt, runner_type); 
    code += st_sim_run_v2(pcxt);
    code += st_sim_report(pcxt, 0);
    
    uint_fast64_t num;
    
    code += st_num_intersections(pcxt, &num);

    double *cos_x = new double[num];
    double *cos_y = new double[num];
    double *cos_z = new double[num];

    code += st_cosines(pcxt, cos_x, cos_y, cos_z);

    code += check_not(cos_x, nullptr);
    code += check_not(cos_y, nullptr);
    code += check_not(cos_z, nullptr);

    uint_fast64_t n = 0;
    glm::dvec3 cosine;
    for (auto it = cxt->p_results->get_ray_record_iterator(); !cxt->p_results->is_at_end(it); ++it)
    {
        ray_record_ptr rec = *it;
        for (auto jt = rec->get_interaction_record_iterator(); !rec->is_at_end(jt); ++jt)
        {
            cosine = (*jt)->direction;
            code += check(cos_x[n], cosine[0]);
            code += check(cos_y[n], cosine[1]);
            code += check(cos_z[n], cosine[2]);
            ++n;
        }
    }

    return code;
}

st_return_t call_stapi_v2_elementmap(st_context_v2_t  pcxt,
                                     st_runner_type_t runner_type)
{
    st_context* cxt = reinterpret_cast<st_context*>(pcxt);
    cxt->p_data->set_number_of_rays(1000);

    st_return_t code = st_sim_setup(pcxt, runner_type); 
    code += st_sim_run_v2(pcxt);
    code += st_sim_report(pcxt, 0);
    
    uint_fast64_t num;
    
    code += st_num_intersections(pcxt, &num);

    uint_fast64_t *els = new uint_fast64_t[num];

    code += st_elementmap(pcxt, els);

    code += check_not(els, nullptr);

    uint_fast64_t n = 0;
    for (auto it = cxt->p_results->get_ray_record_iterator(); !cxt->p_results->is_at_end(it); ++it)
    {
        ray_record_ptr rec = *it;
        for (auto jt = rec->get_interaction_record_iterator(); !rec->is_at_end(jt); ++jt)
        {
            code += check(els[n], (*jt)->element);
            ++n;
        }
    }

    return code;
}

st_return_t call_stapi_v2_stagemap(st_context_v2_t  pcxt,
                                   st_runner_type_t runner_type)
{
    st_context* cxt = reinterpret_cast<st_context*>(pcxt);
    cxt->p_data->set_number_of_rays(1000);

    st_return_t code = st_sim_setup(pcxt, runner_type); 
    code += st_sim_run_v2(pcxt);
    code += st_sim_report(pcxt, 0);
    
    uint_fast64_t num;
    
    code += st_num_intersections(pcxt, &num);

    uint_fast64_t *stages = new uint_fast64_t[num];

    code += st_stagemap(pcxt, stages);

    code += check_not(stages, nullptr);

    for (uint_fast64_t n = 0; n < num; ++n)
        code += check(stages[n], 0);
    
    return code;
}

st_return_t call_stapi_v2_raynumbers(st_context_v2_t  pcxt,
                                     st_runner_type_t runner_type)
{
    st_context* cxt = reinterpret_cast<st_context*>(pcxt);
    cxt->p_data->set_number_of_rays(1000);

    st_return_t code = st_sim_setup(pcxt, runner_type); 
    code += st_sim_run_v2(pcxt);
    code += st_sim_report(pcxt, 0);
    
    uint_fast64_t num;
    
    code += st_num_intersections(pcxt, &num);

    uint_fast64_t *raynumbers = new uint_fast64_t[num];

    code += st_raynumbers(pcxt, raynumbers);

    code += check_not(raynumbers, nullptr);

    uint_fast64_t n = 0;
    for (auto it = cxt->p_results->get_ray_record_iterator(); !cxt->p_results->is_at_end(it); ++it)
    {
        ray_record_ptr rec = *it;
        for (auto jt = rec->get_interaction_record_iterator(); !rec->is_at_end(jt); ++jt)
        {
            code += check(raynumbers[n], rec->id);
            ++n;
        }
    }

    return code;
}

st_return_t call_stapi_v2_sun_stats(st_context_v2_t  pcxt,
                                    st_runner_type_t runner_type)
{
    st_context* cxt = reinterpret_cast<st_context*>(pcxt);
    cxt->p_data->set_number_of_rays(1000);

    st_return_t code = st_sim_setup(pcxt, runner_type); 
    code += st_sim_run_v2(pcxt);
    code += st_sim_report(pcxt, 0);

    double        height;
    double        width;
    double        area;
    uint_fast64_t nrunrays;

    code += st_sun_stats(pcxt,
                         &width,
                         &height,
                         &area,
                         &nrunrays);

    SimulationResult* results =
        reinterpret_cast<SimulationResult*>(cxt->p_results);
    double _width, _height, _area;
    results->get_sun_dimensions(_width, _height);
    _area = results->get_sun_A_box(); 
    if (runner_type == st_runner_type_t::OPTIX)
    {
        code += check(area, _area);
    }
    else
    {
        code += check(width, _width);
        code += check(height, _height);
    }
    uint_fast64_t num = results->get_sun_ray_count();

    code += check(nrunrays, num);

    return code;
}

st_return_t call_stapi_v2_get_results_data(st_context_v2_t  pcxt,
                                           st_runner_type_t runner_type)
{
    st_context* cxt = reinterpret_cast<st_context*>(pcxt);
    cxt->p_data->set_number_of_rays(1000);

    st_return_t code = st_sim_setup(pcxt, runner_type); 
    code += st_sim_run_v2(pcxt);
    code += st_sim_report(pcxt, 0);
    
    uint_fast64_t num;
    
    code += st_num_intersections(pcxt, &num);

    args_results_data data;
    data.loc_x = new double[num];
    data.loc_y = new double[num];
    data.loc_z = new double[num];
    data.cos_x = new double[num];
    data.cos_y = new double[num];
    data.cos_z = new double[num];
    data.element_map = new uint_fast64_t[num];
    data.stage_map   = new uint_fast64_t[num];
    data.ray_numbers = new uint_fast64_t[num];

    code += st_get_results_data(pcxt, &data);

    code += check_not(data.loc_x, nullptr);
    code += check_not(data.loc_y, nullptr);
    code += check_not(data.loc_z, nullptr);
    code += check_not(data.cos_x, nullptr);
    code += check_not(data.cos_y, nullptr);
    code += check_not(data.cos_z, nullptr);
    code += check_not(data.element_map, nullptr);
    code += check_not(data.stage_map, nullptr);
    code += check_not(data.ray_numbers, nullptr);

    uint_fast64_t n = 0;
    glm::dvec3 loc;
    glm::dvec3 cosine;    
    for (auto it = cxt->p_results->get_ray_record_iterator(); !cxt->p_results->is_at_end(it); ++it)
    {
        ray_record_ptr rec = *it;
        for (auto jt = rec->get_interaction_record_iterator(); !rec->is_at_end(jt); ++jt)
        {
            loc = (*jt)->location;
            cosine = (*jt)->direction;
            
            code += check(data.loc_x[n], loc[0]);
            code += check(data.loc_y[n], loc[1]);
            code += check(data.loc_z[n], loc[2]);
            code += check(data.cos_x[n], cosine[0]);
            code += check(data.cos_y[n], cosine[1]);
            code += check(data.cos_z[n], cosine[2]);
            code += check(data.element_map[n], (*jt)->element);
            code += check(data.stage_map[n], 0);
            code += check(data.ray_numbers[n], rec->id);
            ++n;
        }
    }

    return code;
}
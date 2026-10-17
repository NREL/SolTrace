#include <iostream>
#include <vector>
#include <system_error>

#include "../coretrace/simulation_runner/simulation_runner.hpp"
#include "../coretrace/simulation_runner/native_runner/native_runner.hpp"
#include "../coretrace/simulation_data/simulation_data_export.hpp"
#include "../coretrace/simulation_results/simulation_result_export.hpp"

#ifdef STAPI_V2_EMBREE_SUPPORT
#include "../coretrace/simulation_runner/embree_runner/embree_runner.hpp"
#endif

#ifdef STAPI_V2_OPTIX_SUPPORT
#include "../coretrace/simulation_runner/optix_runner/optix_runner.hpp"
#endif

#include "stapi_v2.h"

/* macros for fetching pointers */
#define CONTEXT(p)                                                       \
    st_context *cxt = reinterpret_cast<st_context*>(p);                  \
    if (!cxt || cxt == nullptr) return st_return_code::CONTEXT_NOT_FOUND;

#define DATA(cxt)                                                       \
    SimulationData *data = cxt->p_data;                                 \
    if (!data || data == nullptr) return st_return_code::DATA_NOT_FOUND;

#define RUNNER(cxt)                                                           \
    SimulationRunner *runner = cxt->p_runner;                                 \
    if (!runner || runner == nullptr) return st_return_code::RUNNER_NOT_FOUND;

#define RESULT(cxt)                                                           \
    SimulationResult *result = cxt->p_results;                                \
    if (!result || result == nullptr) return st_return_code::RESULT_NOT_FOUND;

////////////////////////////////
// SolTrace Context Functions //
////////////////////////////////

/* functions for SolTrace context management */
STAPI_V2 st_return_t st_create_context(st_context_v2_t* pcxt, p_callback cb)
{
    st_context* cxt = new st_context();
    cxt->p_data     = new SimulationData();
    cxt->p_cb       = cb;
    cxt->location   = fs::path();

    *pcxt = reinterpret_cast<st_context_v2_t>(cxt);
    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_reset_context(st_context_v2_t pcxt)
{
	CONTEXT(pcxt);

    delete cxt->p_data;
    delete cxt->p_runner;
    delete cxt->p_results;

    cxt->p_data       = new SimulationData();
    cxt->runner_type  = st_runner_type_t::RUNNER_COUNT;
    cxt->p_runner     = nullptr;
    cxt->report_level = RunnerStatistics::STATISTICS_COUNT;
    cxt->p_results    = nullptr;
    cxt->location     = fs::path();

	return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_free_context(st_context_v2_t pcxt)
{
	CONTEXT(pcxt);

    delete cxt->p_data;
    cxt->p_data    = nullptr;
    delete cxt->p_runner;
    cxt->p_runner  = nullptr;
    delete cxt->p_results;
    cxt->p_results = nullptr;

    delete cxt;
	return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_locate_context(st_context_v2_t pcxt, const char *path)
{
	CONTEXT(pcxt);

    fs::path p(path);
    std::error_code ec;

    if (fs::exists(p, ec))
    {
        cxt->location = p;
        return st_return_code::SUCCESS;
    }
    else if (ec && cxt->p_cb)
    {
        cxt->p_cb("locate context", ec.message().c_str());
    }
    return st_return_code::INVALID_ARGUMENTS;
}

////////////////////////////////
// Simlulation Data Functions //
////////////////////////////////

// functions for SolTrace data management
// functions for simulation data management thru json strings
STAPI_V2 st_return_t st_read_input_json(st_context_v2_t pcxt, const char *json)
{
	CONTEXT(pcxt);
    DATA(cxt);
    
    ST_WRAP_CB_TRY_CATCH(data->import_json_string(json), cxt->p_cb);
    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_read_input_json_file(st_context_v2_t pcxt, const char *filename)
{
	CONTEXT(pcxt);
    DATA(cxt);
    
    ST_WRAP_CB_TRY_CATCH(data->import_json_file(std::string(filename)), cxt->p_cb);
    return st_return_code::SUCCESS;
}

// functions for simulation data management directly
STAPI_V2 st_return_t st_set_simulation_parameters(st_context_v2_t pcxt, args_simulation_parameters *params)
{
    CONTEXT(pcxt);
    DATA(cxt);

    SimulationParameters &sim_params = data->get_simulation_parameters();
    
    data->set_number_of_rays(params->number_of_rays);
    data->set_max_rays_traced(params->max_number_of_rays);

    data->set_tolerance(params->tolerance);
    
    data->set_latitude(params->latitude);
    data->set_longitude(params->longitude);

    data->set_include_sun_shape_errors(params->include_sun_shape_errors);
    data->set_include_optical_errors(params->include_optical_errors);
    data->set_as_power_tower(params->as_power_tower);

    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_sim_params(st_context_v2_t pcxt,
								   uint_fast64_t   raycount,
								   uint_fast64_t   maxcount,
								   bool 		   include_dynamic_group)
{
	CONTEXT(pcxt);
    DATA(cxt);

    SimulationParameters &sim_params = data->get_simulation_parameters();
    
    data->set_number_of_rays(raycount);
    data->set_max_rays_traced(maxcount);
    data->set_as_power_tower(include_dynamic_group);
    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_sim_rays(st_context_v2_t pcxt,
								 uint_fast64_t   raycount,
								 uint_fast64_t   maxcount)
{
	CONTEXT(pcxt);
    DATA(cxt);

    SimulationParameters &sim_params = data->get_simulation_parameters();
    
    data->set_number_of_rays(raycount);
    data->set_max_rays_traced(maxcount);
    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_sim_power_tower(st_context_v2_t pcxt,
								        bool 		    is_power_tower)
{
	CONTEXT(pcxt);
    DATA(cxt);

    SimulationParameters &sim_params = data->get_simulation_parameters();

    data->set_as_power_tower(is_power_tower);
    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_sim_errors(st_context_v2_t pcxt,
								   bool			   include_sun_shape,
								   bool			   include_optics)
{
	CONTEXT(pcxt);
    DATA(cxt);
    
    SimulationParameters &sim_params = data->get_simulation_parameters();
    data->set_include_sun_shape_errors((bool)include_sun_shape);
    data->set_include_optical_errors((bool)include_optics);
    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_sim_location(st_context_v2_t pcxt,
								     double		     latitude,
								     double		     longitude)
{
    CONTEXT(pcxt);
    DATA(cxt);

    SimulationParameters &sim_params = data->get_simulation_parameters();

    data->set_latitude(latitude);
    data->set_longitude(longitude);

    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_sim_tolerance(st_context_v2_t pcxt, double tolerance)
{
    CONTEXT(pcxt);
    DATA(cxt);

    SimulationParameters &sim_params = data->get_simulation_parameters();

    data->set_tolerance(tolerance);

    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_get_simulation_parameters(st_context_v2_t pcxt, args_simulation_parameters *params)
{
    CONTEXT(pcxt);
    DATA(cxt);

    SimulationParameters &sim_params = data->get_simulation_parameters();
    
    params->number_of_rays = data->get_number_of_rays();
    params->max_number_of_rays = data->get_max_number_rays_traced();

    params->tolerance = data->get_tolerance();
    
    params->latitude = data->get_latitude();
    params->longitude = data->get_longitude();

    params->include_sun_shape_errors = data->get_include_sun_shape_errors();
    params->include_optical_errors = data->get_include_optical_errors();
    params->as_power_tower = data->get_as_power_tower();

    return st_return_code::SUCCESS;
}

// functions to add/remove optical properties
STAPI_V2 st_return_t st_num_optics(st_context_v2_t pcxt, uint_fast64_t *num_optics)
{
    CONTEXT(pcxt);
    DATA(cxt);

    *num_optics = data->get_number_of_optocal_property_sets();
    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_add_optical_properties_set(st_context_v2_t 				pcxt, 
												   args_optical_properties_set 	*opt_set,
												   args_optical_properties_face *front, 
												   args_optical_properties_face *back, 
												   uint_fast64_t 				*optic_id)
{
    auto bad_type_char = [](char type) 
    {
        return type != 'g' && type != 'p' && type != 'f' && type != 'd';
    };

    if (bad_type_char(front->error_distribution_type) 
        || bad_type_char(back->error_distribution_type)) 
        return st_return_code::INVALID_ARGUMENTS;

    InteractionType inter_type = int_to_interaction(opt_set->type);
    OpticalPropertySet opt(inter_type,
                           opt_set->refraction_index_front,
                           opt_set->refraction_index_back,
                           std::string(opt_set->name));
    
    opt.set_properties(OpticalSide::Front,
                       char_to_distribution(front->error_distribution_type),
                       front->transmissivity,
                       front->reflectivity,
                       front->slope_error,
                       front->specularity_error);

    opt.set_properties(OpticalSide::Back,
                       char_to_distribution(back->error_distribution_type),
                       back->transmissivity,
                       back->reflectivity,
                       back->slope_error,
                       back->specularity_error);
    
    CONTEXT(pcxt);
    DATA(cxt);
    OpticalPropertySetReference res = data->add_optical_property_set(opt);

    if (res.id < 0) return st_return_code::DATA_INSERTION_FAILURE;

    *optic_id = res.id;
    return st_return_code::SUCCESS;
}

const OpticalPropertySetReference get_optics_ref_by_id(SimulationData *data, int_fast64_t opt_id)
{
    optical_set_ptr existing_set;
    for (auto it = data->get_optics_iterator(); !data->is_optics_at_end(it); ++it)
    {
        if (it->first == opt_id)
            existing_set = it->second;
    }

    if (!existing_set) {
        OpticalPropertySetReference ref;
        ref.id = SolTrace::Data::OPTICS_ID_TYPES::OPTICS_ID_UNASSIGNED;
        // ref.optical_property_set = nullptr;
        return ref;
    }

    return data->find_or_add_optical_property_set(*existing_set);
}

STAPI_V2 st_return_t st_get_optical_properties_set(st_context_v2_t 				pcxt, 
												   uint_fast64_t 				optic_id,
												   args_optical_properties_set 	*opt_set,
												   args_optical_properties_face *front, 
												   args_optical_properties_face *back)
{
    CONTEXT(pcxt);
    DATA(cxt);

    // check that inputs are good
    const OpticalPropertySetReference found_set = get_optics_ref_by_id(data, optic_id);
    if (found_set.id == SolTrace::Data::OPTICS_ID_TYPES::OPTICS_ID_UNASSIGNED)
        return st_return_code::DATA_VALUE_NOT_FOUND;

    auto set = found_set.optical_property_set.lock();
    if (!set) return st_return_code::DATA_VALUE_NOT_FOUND;

    opt_set->name = set->get_name().c_str();
    set->get_refraction_indices(opt_set->refraction_index_front, 
                                opt_set->refraction_index_back);
    opt_set->type = interaction_to_int(set->get_interaction_type());

    front->error_distribution_type = distribution_to_char(
        set->get_error_distribution(OpticalSide::Front));
    front->reflectivity = set->get_reflectivity(OpticalSide::Front);
    front->transmissivity = set->get_transmissivity(OpticalSide::Front);
    front->slope_error = set->get_slope_error(OpticalSide::Front);
    front->specularity_error = set->get_specularity_error(OpticalSide::Front);
    
    back->error_distribution_type = distribution_to_char(
        set->get_error_distribution(OpticalSide::Back));
    back->reflectivity = set->get_reflectivity(OpticalSide::Back);
    back->transmissivity = set->get_transmissivity(OpticalSide::Back);
    back->slope_error = set->get_slope_error(OpticalSide::Back);
    back->specularity_error = set->get_specularity_error(OpticalSide::Back);

    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_delete_optic(st_context_v2_t pcxt, st_uint_t idx)
{
    CONTEXT(pcxt);
    DATA(cxt);

    int removed = data->remove_optical_property_set(idx);

    // only warn if requested to remove optical property set that doesn't exist
    return removed > 0 ? st_return_code::SUCCESS : st_return_code::WARNING_NOT_FOUND;
}

STAPI_V2 st_return_t st_clear_optics(st_context_v2_t pcxt)
{
    CONTEXT(pcxt);
    DATA(cxt);

    data->clear_optical_property_sets();
    return st_return_code::SUCCESS;
}

// functions to add/remove elements
STAPI_V2 st_return_t st_num_elements(st_context_v2_t pcxt, uint_fast64_t *num_elements)
{
    CONTEXT(pcxt);
    DATA(cxt);

    *num_elements = data->get_number_of_elements();
    return st_return_code::SUCCESS;
}

int32_t num_aperture_params(ApertureType type)
{
    int32_t rt = -1;
    
    switch (type)
    {
        case ApertureType::CIRCLE:
        case ApertureType::EQUILATERAL_TRIANGLE:
        case ApertureType::HEXAGON:
        {
            rt = 1;
            break;
        }
        case ApertureType::RECTANGLE:
        {
            rt = 2;
            break;
        }
        case ApertureType::ANNULUS:
        case ApertureType::SINGLE_AXIS_CURVATURE_SECTION:
        {
            rt = 3;
            break;
        }
        case ApertureType::IRREGULAR_TRIANGLE:
        {
            rt = 6;
            break;
        }
        case ApertureType::IRREGULAR_QUADRILATERAL:
        {
            rt = 8;
            break;
        }
        default:
            rt = -1;
            break;
        }
    return rt;
}

int32_t num_surface_params(SurfaceType type)
{
    int32_t rt = -1;

    switch (type)
    {
        case SurfaceType::FLAT:
        {
            rt = 0;
            break;
        }
        case SurfaceType::CONE:
        case SurfaceType::CYLINDER:
        case SurfaceType::SPHERE:
        {
            rt = 1;
            break;
        }
        case SurfaceType::PARABOLA:
        {
            rt = 2;
            break;
        }
        case SurfaceType::HYPER:
        case SurfaceType::GENERAL_SPENCER_MURTY:
        case SurfaceType::TORUS:
        default:
            rt = -1; // Not implemented yet
            break;
    }

    return rt;
}

st_return_t create_geometry(st_context    *cxt,
                            args_geometry *geom,
                            aperture_ptr  &ap,
                            surface_ptr   &surf)
{
    // checking geometry
    ApertureType a_type = char_to_aperture(geom->ap);
    SurfaceType s_type = char_to_surface(geom->surf);
    int32_t num_a_params = num_aperture_params(a_type);
    int32_t num_s_params = num_surface_params(s_type);
    if (num_a_params < 0 || num_s_params < 0) return st_return_code::INVALID_ARGUMENTS;

    std::vector<double> _a_params(std::begin(geom->a_params), std::end(geom->a_params));
    std::vector<double> _s_params(std::begin(geom->s_params), std::end(geom->s_params));

    try
    {
        ap = Aperture::make_aperture_from_type(a_type, _a_params);
        surf = SolTrace::Data::make_surface_from_type(s_type, _s_params);
    }
    catch (const std::invalid_argument& e)
    {
        if (cxt->p_cb) cxt->p_cb("geometry creation", e.what());
        return st_return_code::INVALID_ARGUMENTS;
    }

    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_add_element(st_context_v2_t pcxt,
									args_element    *args,
									int_fast64_t    opt_id,
									args_geometry   *geom,
									uint_fast64_t   *element_id)
{
    CONTEXT(pcxt);
    DATA(cxt);

    // check that inputs are good
    // checking optics
    const OpticalPropertySetReference existing_set = get_optics_ref_by_id(data, opt_id);
    if (existing_set.id == SolTrace::Data::OPTICS_ID_TYPES::OPTICS_ID_UNASSIGNED)
    return st_return_code::DATA_VALUE_NOT_FOUND;
    
    aperture_ptr ap;
    surface_ptr surf;
    st_return_t code = create_geometry(cxt, geom, ap, surf);

    if (code > st_return_code::SUCCESS) return code;

    single_element_ptr el = make_element<SingleElement>();
    
    el->set_origin(args->x, args->y, args->z);
    el->set_aim_vector(args->ax, args->ay, args->az);
    el->set_zrot(args->zrot);
    if (args->enabled_flag) el->enable();
    else                    el->disable();
    if (args->virtual_flag) el->mark_virtual();
    else                    el->unmark_virtual();
    
    el->set_optical_property_set(existing_set);
    el->set_aperture(ap);
    el->set_surface(surf);

    // TODO in simulation_data.cpp saying add_element will be throwable in the future
    ST_WRAP_CB_TRY_CATCH(data->add_element(el), cxt->p_cb);
    
    *element_id = el->get_id();
    
    // to maintain grouping, changing groups goes thru SimulationData, so set group after adding the element
    if (args->group < -1) code = st_return_code::WARNING_GROUP_IGNORED;
    else                  data->set_element_group(*element_id, args->group);
    
    return code;
}

STAPI_V2 st_return_t st_get_element(st_context_v2_t pcxt,
									uint_fast64_t   idx,
									args_element    *args,
									int_fast64_t    *opt_id,
									args_geometry   *geom)
{
    CONTEXT(pcxt);
    DATA(cxt);
    
    element_ptr el = data->get_element(idx);
    if (!el) return st_return_code::DATA_VALUE_NOT_FOUND;

    // get element
    args->enabled_flag = el->is_enabled();
    args->virtual_flag = el->is_virtual();
    glm::dvec3 pos = el->get_origin_global();
    args->x = pos[0];
    args->y = pos[1];
    args->z = pos[2];
    glm::dvec3 aim = el->get_aim_vector_global();
    args->ax = aim[0];
    args->ay = aim[1];
    args->az = aim[2];
    args->zrot = el->get_zrot();
    args->group = el->get_group();
    *opt_id = el->get_optical_property_set_id();
    
    // get geometry
    aperture_ptr ap_ptr = el->get_aperture();
    ApertureType ap_type = ap_ptr->get_type();
    geom->ap = aperture_to_char(ap_type);
    surface_ptr surf_ptr = el->get_surface();
    SurfaceType surf_type = surf_ptr->get_type();
    geom->surf = surface_to_char(surf_type);
    // reset params
    for (size_t i = 0; i < 8; ++i)
    {
        geom->a_params[i] = 0;
        geom->s_params[i] = 0;
    }
    
    switch (ap_type)
    {
        case ApertureType::ANNULUS:
        {
            auto temp = std::dynamic_pointer_cast<Annulus>(ap_ptr);
            geom->a_params[0] = temp->inner_radius;
            geom->a_params[1] = temp->outer_radius;
            geom->a_params[2] = temp->arc_angle;
            break;
        }
        case ApertureType::CIRCLE:
        {
            auto temp = std::dynamic_pointer_cast<Circle>(ap_ptr);
            geom->a_params[0] = temp->diameter;
            break;
        }
        case ApertureType::HEXAGON:
        {
            auto temp = std::dynamic_pointer_cast<Hexagon>(ap_ptr);
            geom->a_params[0] = temp->circumscribe_diameter;
            break;
        }
        case ApertureType::RECTANGLE:
        {
            auto temp = std::dynamic_pointer_cast<Rectangle>(ap_ptr);
            geom->a_params[0] = temp->x_length();
            geom->a_params[1] = temp->y_length();
            break;
        }
        case ApertureType::EQUILATERAL_TRIANGLE:
        {
            auto temp = std::dynamic_pointer_cast<EquilateralTriangle>(ap_ptr);
            geom->a_params[0] = temp->circumscribe_diameter;
            break;
        }
        
        // case ApertureType::SINGLE_AXIS_CURVATURE_SECTION:

        case ApertureType::IRREGULAR_TRIANGLE:
        {
            auto temp = std::dynamic_pointer_cast<IrregularTriangle>(ap_ptr);
            geom->a_params[0] = temp->x1;
            geom->a_params[1] = temp->y1;
            geom->a_params[2] = temp->x2;
            geom->a_params[3] = temp->y2;
            geom->a_params[4] = temp->x3;
            geom->a_params[5] = temp->y3;
            break;
        }
        case ApertureType::IRREGULAR_QUADRILATERAL:
        {
            auto temp = std::dynamic_pointer_cast<IrregularQuadrilateral>(ap_ptr);
            geom->a_params[0] = temp->x1;
            geom->a_params[1] = temp->y1;
            geom->a_params[2] = temp->x2;
            geom->a_params[3] = temp->y2;
            geom->a_params[4] = temp->x3;
            geom->a_params[5] = temp->y3;
            geom->a_params[6] = temp->x4;
            geom->a_params[7] = temp->y4;
            break;
        }
        default:
            break;
    }

    switch (surf_type)
    {
        case SurfaceType::CONE:
        {
            auto temp = std::dynamic_pointer_cast<Cone>(surf_ptr);
            geom->s_params[0] = temp->half_angle;
            break;
        }
        case SurfaceType::CYLINDER:
        {
            auto temp = std::dynamic_pointer_cast<Cylinder>(surf_ptr);
            geom->s_params[0] = temp->radius;
            break;
        }
        case SurfaceType::FLAT:
        {
            break;
        }
        case SurfaceType::PARABOLA:
        {
            auto temp = std::dynamic_pointer_cast<Parabola>(surf_ptr);
            geom->s_params[0] = temp->focal_length_x;
            geom->s_params[1] = temp->focal_length_y;
            break;
        }
        case SurfaceType::SPHERE:
        {
            auto temp = std::dynamic_pointer_cast<Sphere>(surf_ptr);
            geom->s_params[0] = temp->vertex_curv;
            break;
        }

        // case SurfaceType::HYPER:
        // case SurfaceType::GENERAL_SPENCER_MURTY:
        // case SurfaceType::TORUS:

        default:
            break;
    }

    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_delete_element(st_context_v2_t pcxt, st_uint_t idx)
{
    CONTEXT(pcxt);
    DATA(cxt);

    uint_fast64_t removed = data->remove_element(idx);

    if (removed == 0) return st_return_code::WARNING_NOT_FOUND;
    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_clear_elements(st_context_v2_t pcxt)
{
    CONTEXT(pcxt);
    DATA(cxt);

    data->clear_elements();
    return st_return_code::SUCCESS;
}

// functions to modify elements
STAPI_V2 st_return_t st_element_enabled(st_context_v2_t pcxt,
										st_uint_t 		idx,
										bool 			enabled_flag)
{
    CONTEXT(pcxt);
    DATA(cxt);
    
    element_ptr el = data->get_element(idx);
    if (!el) return st_return_code::WARNING_NOT_FOUND;

    if (enabled_flag) el->enable();
    else              el->disable();

    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_element_virtual(st_context_v2_t pcxt,
										st_uint_t 		idx,
										bool 			virtual_flag)
{
    CONTEXT(pcxt);
    DATA(cxt);
    
    element_ptr el = data->get_element(idx);
    if (!el) return st_return_code::WARNING_NOT_FOUND;

    if (virtual_flag) el->mark_virtual();
    else              el->unmark_virtual();

    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_element_xyz(st_context_v2_t pcxt, 
									st_uint_t 		idx,
									double 	  		x,
									double 	  		y,
									double 	  		z)
{
    CONTEXT(pcxt);
    DATA(cxt);

    element_ptr el = data->get_element(idx);
    if (!el) return st_return_code::WARNING_NOT_FOUND;

    el->set_origin(x, y, z);
    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_element_aim(st_context_v2_t pcxt, 
									st_uint_t 		idx,
									double 	  		ax,
									double 	  		ay,
									double 	  		az)
{
    CONTEXT(pcxt);
    DATA(cxt);
    
    element_ptr el = data->get_element(idx);
    if (!el) return st_return_code::WARNING_NOT_FOUND;

    el->set_aim_vector(ax, ay, az);
    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_element_zrot(st_context_v2_t pcxt,
									 st_uint_t 		 idx,
									 double 		 zrot)
{
    CONTEXT(pcxt);
    DATA(cxt);

    element_ptr el = data->get_element(idx);
    if (!el) return st_return_code::WARNING_NOT_FOUND;

    el->set_zrot(zrot);
    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_element_aperture(st_context_v2_t pcxt,
                                         st_uint_t 	     idx,
                                         char      	     ap,
                                         double    	     params[8])
{
    CONTEXT(pcxt);
    DATA(cxt);

    element_ptr el = data->get_element(idx);
    if (!el) return st_return_code::WARNING_NOT_FOUND;
    
    ApertureType ap_type = char_to_aperture(ap);
    int32_t num_params = num_aperture_params(ap_type);
    if (num_params < 0) return st_return_code::INVALID_ARGUMENTS;

    std::vector<double> _params(params, params + num_params);
    aperture_ptr a;
    try
    {
        a = Aperture::make_aperture_from_type(ap_type, _params);
    }
    catch (const std::invalid_argument& e)
    {
        if (cxt->p_cb) cxt->p_cb("geometry creation", e.what());
        return st_return_code::INVALID_ARGUMENTS;
    }

    el->set_aperture(a);
    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_element_surface(st_context_v2_t pcxt,
                                        st_uint_t 	    idx,
                                        char      	    surf,
                                        double    	    params[8])
{
    CONTEXT(pcxt);
    DATA(cxt);
    
    element_ptr el = data->get_element(idx);
    if (!el) return st_return_code::WARNING_NOT_FOUND;
    
    SurfaceType surf_type = char_to_surface(surf);
    int32_t num_params = num_surface_params(surf_type);
    if (num_params < 0) return st_return_code::INVALID_ARGUMENTS;

    std::vector<double> _params(params, params + num_params);
    surface_ptr s;
    try
    {
        s = SolTrace::Data::make_surface_from_type(surf_type, _params);
    }
    catch (const std::invalid_argument& e)
    {
        if (cxt->p_cb) cxt->p_cb("geometry creation", e.what());
        return st_return_code::INVALID_ARGUMENTS;
    }

    el->set_surface(s);
    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_element_geometry(st_context_v2_t pcxt,
                                         st_uint_t 	     idx,
                                         args_geometry	 *geom)
{
    CONTEXT(pcxt);
    DATA(cxt);
    
    element_ptr el = data->get_element(idx);
    if (!el) return st_return_code::WARNING_NOT_FOUND;

    aperture_ptr ap;
    surface_ptr surf;
    st_return_t code = create_geometry(cxt, geom, ap, surf);

    if (code > st_return_code::SUCCESS) return code;
    
    el->set_aperture(ap);
    el->set_surface(surf);
    
    return code;
}

STAPI_V2 st_return_t st_element_optic(st_context_v2_t pcxt,
									  st_uint_t 	  idx,
									  int_fast64_t 	  opt_id)
{
    CONTEXT(pcxt);
    DATA(cxt);
    
    const OpticalPropertySetReference existing_set = get_optics_ref_by_id(data, opt_id);
    if (existing_set.id == SolTrace::Data::OPTICS_ID_TYPES::OPTICS_ID_UNASSIGNED)
        return st_return_code::DATA_VALUE_NOT_FOUND;

    element_ptr el = data->get_element(idx);
    if (!el) return st_return_code::WARNING_NOT_FOUND;

    el->set_optical_property_set(existing_set);
    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_element_group(st_context_v2_t pcxt,
									  st_uint_t 	  idx,
									  int32_t   	  group)
{
    CONTEXT(pcxt);
    DATA(cxt);
    
    if (group < -1) return st_return_code::WARNING_GROUP_IGNORED;

    element_ptr el = data->get_element(idx);
    if (!el) return st_return_code::WARNING_NOT_FOUND;

    data->set_element_group(idx, group);

    return st_return_code::SUCCESS;
}

// sun functions
std::shared_ptr<Sun> get_or_create_sun(SimulationData *data)
{
    if (!data->get_number_of_ray_sources())
        data->add_ray_source(make_ray_source<Sun>());
    
    ray_source_ptr sun_ptr = data->get_ray_source(0);
    return std::dynamic_pointer_cast<Sun>(sun_ptr);
}

STAPI_V2 st_return_t st_add_sun(st_context_v2_t pcxt,
								args_sun 		*args,
								double 		 	*angle,
								double 		 	*intensity)
{
    CONTEXT(pcxt);
    DATA(cxt);

    st_return_t code = st_return_code::SUCCESS;
    auto sun = get_or_create_sun(data);

    sun->set_position(args->x, args->y, args->z);

    if (args->npoints)
    {
        std::vector<double> v_angle(angle, 
                                    angle + args->npoints);
        std::vector<double> v_intensity(intensity, 
                                        intensity + args->npoints);
        ST_WRAP_CB_TRY_CATCH(sun->set_shape(SunShape::USER_DEFINED,
                                            0, 0, 0,
                                            v_angle, 
                                            v_intensity),
                            cxt->p_cb);
    }
    else
    {
        SunShape sun_shape = SunShape::GAUSSIAN;
        double s_hw_csr = args->sigma_halfwidth_csr;
        switch (args->shape)
        {
            /* default to gaussian shape, 
            using default to add warning */
            case 'g':
            case 'G':
                // sun_shape = SunShape::GAUSSIAN;
                break;
            case 'p':
            case 'P':
                sun_shape = SunShape::PILLBOX;
                break;
            /* currently no set_limbdarkend_distribution
            function is implemented, so ignore and 
            emit warning
            case 'l':
            case 'L':
            {
                sun_shape = SunShape::LIMBDARKENED;
                break;
            }                                           */
            case 'b':
            case 'B':
                sun_shape = SunShape::BUIE_CSR;
                break;
            /* user defined sun shape goes thru same 
            Sun::set_shape, emit warning because not 
            enough information from this signiture, 
            maybe add defaults?                      
            case 'u':
            case 'U':                             */
        default:
            // warning code default to gaussian, default sigma = 4.65
            sun_shape = SunShape::GAUSSIAN;
            s_hw_csr = 4.65;
            code = st_return_code::WARNING_SUN_SHAPE_IGNORED;
            break;
        }

        auto sun = get_or_create_sun(data);
        ST_WRAP_CB_TRY_CATCH(sun->set_shape(sun_shape, 
                                            s_hw_csr,
                                            s_hw_csr,
                                            s_hw_csr),
                             cxt->p_cb);
    }

    return code;
}

STAPI_V2 st_return_t st_get_sun(st_context_v2_t pcxt,
								args_sun 		*args,
								double 		 	**angle,
								double 		 	**intensity,
								st_uint_t 		npoints)
{
    CONTEXT(pcxt);
    DATA(cxt);

    if (!data->get_number_of_ray_sources())
        return st_return_code::DATA_VALUE_NOT_FOUND;
 
    auto sun = data->get_ray_source(0);
    glm::dvec3 pos = sun->get_position();
    
    args->npoints = 0;
    args->x = pos[0];
    args->y = pos[1];
    args->z = pos[2];

    SunShape shape = sun->get_shape();
    args->shape = sunshape_to_char(shape);

    switch (shape)
    {
        case SunShape::GAUSSIAN:
        {
            args->sigma_halfwidth_csr = sun->get_sigma();
            break;
        }
        case SunShape::PILLBOX:
        {
            args->sigma_halfwidth_csr = sun->get_half_width();
            break;
        }
        case SunShape::BUIE_CSR:
        {
            args->sigma_halfwidth_csr = sun->get_circumsolar_ratio();
            break;
        }
        case SunShape::USER_DEFINED:
        {
            std::vector<double> _angle;
            std::vector<double> _intensity;
            sun->get_user_data(_angle, _intensity);

            if (npoints == 0 || npoints != _angle.size())
                return st_return_code::INVALID_ARGUMENTS;
            
            args->npoints = _angle.size();

            std::copy(_angle.begin(), _angle.end(), *angle);
            std::copy(_intensity.begin(), _intensity.end(), *intensity);
            break;
        }
    }

    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_sun_shape(st_context_v2_t pcxt,
							      char            shape, 
							      double          sigma_halfwidth_csr)
{
    CONTEXT(pcxt);
    DATA(cxt);

    st_return_t code = st_return_code::SUCCESS; 
    SunShape sun_shape = SunShape::GAUSSIAN;

    switch (shape)
    {
        /* default to gaussian shape, 
           using default to add warning */
        case 'g':
        case 'G':
            // sun_shape = SunShape::GAUSSIAN;
            break;
        case 'p':
        case 'P':
            sun_shape = SunShape::PILLBOX;
            break;
        /* currently no set_limbdarkend_distribution
           function is implemented, so ignore and 
           emit warning
        case 'l':
        case 'L':
        {
            sun_shape = SunShape::LIMBDARKENED;
            break;
        }                                           */
        case 'b':
        case 'B':
            sun_shape = SunShape::BUIE_CSR;
            break;
        /* user defined sun shape goes thru same 
           Sun::set_shape, emit warning because not 
           enough information from this signiture, 
           maybe add defaults?                      
           case 'u':
           case 'U':                             */
    default:
        // warning code default to gaussian, default sigma = 4.65
        sun_shape = SunShape::GAUSSIAN;
        sigma_halfwidth_csr = 4.65;
        code = st_return_code::WARNING_SUN_SHAPE_IGNORED;
        break;
    }

    auto sun = get_or_create_sun(data);
    ST_WRAP_CB_TRY_CATCH(sun->set_shape(sun_shape,
                                        sigma_halfwidth_csr,
                                        sigma_halfwidth_csr,
                                        sigma_halfwidth_csr),
                         cxt->p_cb);
    return code;
}

STAPI_V2 st_return_t st_sun_xyz(st_context_v2_t pcxt,
								double          x,
								double          y,
								double          z)
{
    CONTEXT(pcxt);
    DATA(cxt);
    
    auto sun = get_or_create_sun(data);
    sun->set_position(x, y, z);
    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_sun_position(st_context_v2_t pcxt,
									 double          lat,
									 double          day,
									 double          hour,
									 double          *x,
									 double          *y,
									 double          *z)
{
    /* TODO: make st_sun_position_v2 that take a 
	computes the sun vector xyz given arguments
	lat : [deg] latitude 
	day : [] day of the year 
	hour : [hour] solar time. 12.00 corresponds to sun at maximum elevation and does not necessarily match local time

	xyz coordinate system:
		x: +west
		y: +zenith
		z: +north
	*/

	double Declination, HourAngle, Elevation, Azimuth;

	Declination = 180 / M_PI * asin(0.39795 * cos(0.98563 * M_PI / 180 * (day - 173)));
	HourAngle = 15 * (hour - 12);
	Elevation = 180 / M_PI * asin(sin(Declination * M_PI / 180) * sin(lat * M_PI / 180) + cos(Declination * M_PI / 180) * cos(HourAngle * M_PI / 180) * cos(lat * M_PI / 180));
	Azimuth = 180 / M_PI * acos((sin(M_PI / 180 * Declination) * cos(M_PI / 180 * lat) - cos(M_PI / 180 * Declination) * sin(M_PI / 180 * lat) * cos(M_PI / 180 * HourAngle)) / cos(M_PI / 180 * Elevation) + 0.0000000001);
	if (sin(HourAngle * M_PI / 180) > 0.0)
		Azimuth = 360 - Azimuth;
	*x = -sin(Azimuth * M_PI / 180) * cos(Elevation * M_PI / 180);
	*y = sin(Elevation * M_PI / 180);
	*z = cos(Azimuth * M_PI / 180) * cos(Elevation * M_PI / 180);

	return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_sun_userdata(st_context_v2_t pcxt,
									 st_uint_t 		 npoints,
									 double 		 *angle,
									 double 		 *intensity)
{
    CONTEXT(pcxt);
    DATA(cxt);

    std::vector<double> v_angle(angle, angle + npoints);
    std::vector<double> v_intensity(intensity, intensity + npoints);
    auto sun = get_or_create_sun(data);
    ST_WRAP_CB_TRY_CATCH(sun->set_shape(SunShape::USER_DEFINED,
                                        0, 0, 0,
                                        v_angle, 
                                        v_intensity),
                        cxt->p_cb);
    return st_return_code::SUCCESS;
}

st_return_code int_to_calc_method(int *arg, spcm *calc_method)
{
    if (*arg >= spcm::CALCULATOR_COUNT)
        return st_return_code::INVALID_ARGUMENTS;

    *calc_method = static_cast<spcm>(*arg);
    return st_return_code::SUCCESS;
}

void set_up_calculator(spcm				       calc,
					   args_sun_location       *loc,
					   args_sun_datetime       *dt,
                       SolarPositionCalculator *solar_position)
{
    solar_position->set_location(loc->latitude, loc->longitude, loc->timeZone, loc->altitude);
    solar_position->set_datetime(dt->year, dt->month, dt->day, dt->hour, dt->minute, dt->second);
    solar_position->set_method(calc);
}

STAPI_V2 st_return_t st_get_sun_az_zen(st_context_v2_t   pcxt,
									   int				 calc,
									   args_sun_location *loc,
									   args_sun_datetime *dt,
									   double 			 *azimuth,
									   double 			 *zenith)
{
    spcm method;
    st_return_code code = int_to_calc_method(&calc, &method);
    if (code != st_return_code::SUCCESS) return code;

    CONTEXT(pcxt);
    SolarPositionCalculator solar_position;
    
    ST_WRAP_CB_TRY_CATCH(set_up_calculator(method, loc, dt, &solar_position), cxt->p_cb);
    ST_WRAP_CB_TRY_CATCH(solar_position.get_azimuth_zenith(azimuth, zenith), cxt->p_cb);

    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_get_sun_az_el(st_context_v2_t   pcxt,
									  int				calc,
									  args_sun_location *loc,
									  args_sun_datetime *dt,
									  double 			*azimuth,
									  double 			*elevation)
{
    spcm method;
    st_return_code code = int_to_calc_method(&calc, &method);
    if (code != st_return_code::SUCCESS) return code;

    CONTEXT(pcxt);
    SolarPositionCalculator solar_position;
    
    ST_WRAP_CB_TRY_CATCH(set_up_calculator(method, loc, dt, &solar_position), cxt->p_cb);
    ST_WRAP_CB_TRY_CATCH(solar_position.get_azimuth_elevation(azimuth, elevation), cxt->p_cb);

    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_get_sun_vector(st_context_v2_t   pcxt,
									   int				 calc,
									   args_sun_location *loc,
									   args_sun_datetime *dt,
									   double 			 *sun_x,
									   double 			 *sun_y,
									   double 			 *sun_z)
{
    spcm method;
    st_return_code code = int_to_calc_method(&calc, &method);
    if (code != st_return_code::SUCCESS) return code;

    CONTEXT(pcxt);
    SolarPositionCalculator solar_position;
    
    ST_WRAP_CB_TRY_CATCH(set_up_calculator(method, loc, dt, &solar_position), cxt->p_cb);
    ST_WRAP_CB_TRY_CATCH(solar_position.get_sun_vector(sun_x, sun_y, sun_z), cxt->p_cb);

    return st_return_code::SUCCESS;
}

// functions for writing input files for SolTrace
STAPI_V2 st_return_t st_export_json_file(st_context_v2_t pcxt, const char *filename)
{
    CONTEXT(pcxt);
    DATA(cxt);

    ST_WRAP_CB_TRY_CATCH(data->export_json_file(std::string(filename)), cxt->p_cb);
    return st_return_code::SUCCESS;
}

//////////////////////////////////
// Simlulation Runner Functions //
//////////////////////////////////

// functions for SolTrace runner management
STAPI_V2 st_return_t st_get_installed_runners(uint8_t *installed)
{
    *installed = 1 << st_runner_type_t::NATIVE;
#ifdef STAPI_V2_EMBREE_SUPPORT
    *installed += 1 << st_runner_type_t::EMBREE;
#endif
#ifdef STAPI_V2_OPTIX_SUPPORT
    *installed += 1 << st_runner_type_t::OPTIX;
#endif
    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_is_runner_installed(st_runner_type_t type,
                                            bool             *installed)
{
    uint8_t runners;
    st_get_installed_runners(&runners);
    *installed = runners & (1 << type);
    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_sim_setup(st_context_v2_t  pcxt, 
								  st_runner_type_t runner_type, 
								  uint_fast64_t    num_threads,
								  st_uint_t 	   *seeds,
								  size_t		   num_seeds)
{
    CONTEXT(pcxt);
    DATA(cxt);

    delete cxt->p_runner;
    delete cxt->p_results;
    cxt->p_runner  = nullptr;
    cxt->p_results = nullptr;

    bool use_embree = runner_type == st_runner_type_t::EMBREE;
    bool use_optix = runner_type == st_runner_type_t::OPTIX;

    RunnerStatus sts;
    st_return_t  rt;
    rt = st_return_code::SUCCESS;
    SimulationRunner *runner;

#ifdef STAPI_V2_EMBREE_SUPPORT
    if (use_embree)
    {
        runner = new SolTrace::EmbreeRunner::EmbreeRunner();
    }
    else
#endif
#ifdef STAPI_V2_OPTIX_SUPPORT
        if (use_optix)
    {
        runner = new OptixRunner();
    }
    else
#endif
    {
#ifndef STAPI_V2_EMBREE_SUPPORT
        if (use_embree) 
        {
            rt = st_return_code::WARNING_FELLBACK_FROM_EMBREE;
        }
#endif
#ifndef STAPI_V2_OPTIX_SUPPORT
        if (use_optix) 
        {
            rt = st_return_code::WARNING_FELLBACK_FROM_OPTIX;
        }
#endif
        runner_type = st_runner_type_t::NATIVE;
        runner = new NativeRunner();
    }

    sts = runner->initialize();
    if (sts != RunnerStatus::SUCCESS)
    {
        delete runner;
        return st_return_code::RUNNER_INILIALIZE_FAILURE;
    }

    /* optix doesn't use threads the way native/embree does, 
       and check if user requested optix but didn't build it,
       in either case, set the number of threads for the runner */
    if (!use_optix || rt == st_return_code::WARNING_FELLBACK_FROM_OPTIX)
    {
        NativeRunner *temp_native = dynamic_cast<NativeRunner*>(runner);
        if (seeds != nullptr)
        {
            // ensures runtime_error won't be raised
            if (num_threads != num_seeds)
            {
                delete runner;
                return st_return_code::RUNNER_NUMBER_THREADS_SEEDS_MISMATCH_FAILURE;
            }

            const std::vector<st_uint_t> temp_seeds(seeds, seeds + num_seeds);
            temp_native->set_number_of_threads(num_threads, temp_seeds);
        }
        else
        {
            temp_native->set_number_of_threads(num_threads);
        }
        // by default disable stages
        temp_native->disable_stages();
    }
    // otherwise using optix runner
    else {
        // if context is located, add for ptx discovery, include ifdef guard b/c need to access OptixRunner
#ifdef STAPI_V2_OPTIX_SUPPORT
        if (!cxt->location.empty())
        {
            OptixRunner *temp_optix = dynamic_cast<OptixRunner*>(runner);
            temp_optix->set_additional_ptx_directory(cxt->location);
        }
#endif
        // if using optix runner and requested threads, emit warning that it was ignored 
        if (num_threads != DEFAULT_NUM_THREADS) rt = st_return_code::WARNING_ARGUMENT_IGNORED_BY_RUNNER;
    }

    // auto t_setup_start = std::chrono::steady_clock::now();
    
    try { sts = runner->setup_simulation(data); }
    catch (std::runtime_error &e)
    {
        if (cxt->p_cb) cxt->p_cb("set up runner", e.what());
        delete runner;
        return st_return_code::RUNNER_SETUP_ERROR;
    }
    
    // auto t_setup_end = std::chrono::steady_clock::now();
    if (sts != RunnerStatus::SUCCESS)
    {
        delete runner;
        return st_return_code::RUNNER_SETUP_FAILURE;
    }
    
    cxt->runner_type = runner_type;
    cxt->p_runner = runner;

    return rt;
}

STAPI_V2 st_return_t st_sim_run_v2(st_context_v2_t pcxt)
{
    CONTEXT(pcxt);
    RUNNER(cxt);

    if (!runner->is_ready_to_run()) return st_return_code::RUNNER_NOT_READY_TO_RUN;

    RunnerStatus sts = runner->run_simulation();

    if (sts == RunnerStatus::CANCEL) return st_return_code::CANCEL;
    else if (sts != RunnerStatus::SUCCESS) return st_return_code::FAILURE;
    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_sim_report(st_context_v2_t pcxt, int level)
{
    if (level < RunnerStatistics::RAY_RECORDS ||
        level >= RunnerStatistics::STATISTICS_COUNT)
        return st_return_code::INVALID_ARGUMENTS;
        
    CONTEXT(pcxt);
    RUNNER(cxt);

    if (!runner->is_ready_to_report()) return st_return_code::RUNNER_NOT_READY_TO_REPORT;

    cxt->report_level = (RunnerStatistics)level;
    delete cxt->p_results;
    cxt->p_results = new SimulationResult();
    RunnerStatus sts = runner->report_simulation(cxt->p_results, level);

    if (sts != RunnerStatus::SUCCESS) return st_return_code::FAILURE;
    return st_return_code::SUCCESS;
}

///////////////////////////////////
// Simlulation Results Functions //
///////////////////////////////////

STAPI_V2 st_return_t st_write_results_csv(st_context_v2_t pcxt, 
										  const char 	  *filename, 
										  int 			  precision)
{
    CONTEXT(pcxt);
    RESULT(cxt);

    if (cxt->report_level == RunnerStatistics::GROUPED_COUNTS)
        return st_return_code::RESULT_NOT_REPORTED;

    ST_WRAP_CB_TRY_CATCH(result->write_csv_file(filename, precision), cxt->p_cb);
    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_write_group_results_json(st_context_v2_t pcxt, 
										         const char 	 *filename, 
										         int 			 precision,
                                                 int             indent)
{
    CONTEXT(pcxt);
    RESULT(cxt);

    if (cxt->report_level == RunnerStatistics::RAY_RECORDS)
        return st_return_code::RESULT_NOT_REPORTED;

    ST_WRAP_CB_TRY_CATCH(result->write_group_json_file(filename, precision, indent), cxt->p_cb);
    return st_return_code::SUCCESS;
}

// functions to get results directly
uint_fast64_t number_of_interactions(SimulationResult *result)
{
    uint_fast64_t num = 0;

    for (auto it = result->get_ray_record_iterator(); !result->is_at_end(it); ++it)
        num += (*it)->get_number_of_interactions();
    
    return num;
}

STAPI_V2 st_return_t st_num_intersections(st_context_v2_t pcxt, uint_fast64_t *num_intersections)
{
    CONTEXT(pcxt);
    RESULT(cxt);
    
    if (cxt->report_level == RunnerStatistics::GROUPED_COUNTS)
        return st_return_code::RESULT_NOT_REPORTED;

    *num_intersections = number_of_interactions(result);
    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_locations(st_context_v2_t pcxt,
								  double 		  *loc_x,
								  double 		  *loc_y,
								  double 		  *loc_z)
{
    /* *loc_x, *loc_y, *loc_z must be allocated by 
       caller using st_num_intersections */
    CONTEXT(pcxt);
    RESULT(cxt);
    
    if (cxt->report_level == RunnerStatistics::GROUPED_COUNTS)
        return st_return_code::RESULT_NOT_REPORTED;

    uint_fast64_t n = 0;
    glm::dvec3 loc;
    for (auto it = result->get_ray_record_iterator(); !result->is_at_end(it); ++it)
    {
        ray_record_ptr rec = *it;
        for (auto jt = rec->get_interaction_record_iterator(); !rec->is_at_end(jt); ++jt)
        {
            loc = (*jt)->location;
            loc_x[n] = loc[0];
            loc_y[n] = loc[1];
            loc_z[n] = loc[2];
            ++n;
        }
    }

    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_cosines(st_context_v2_t pcxt,
							 	double 		    *cos_x,
							 	double 		    *cos_y,
							 	double 		    *cos_z)
{
    /* *cos_x, *cos_y, *cos_z must be allocated by 
       caller using st_num_intersections */
    CONTEXT(pcxt);
    RESULT(cxt);
    
    if (cxt->report_level == RunnerStatistics::GROUPED_COUNTS)
        return st_return_code::RESULT_NOT_REPORTED;

    uint_fast64_t n = 0;
    glm::dvec3 cosine;
    for (auto it = result->get_ray_record_iterator(); !result->is_at_end(it); ++it)
    {
        ray_record_ptr rec = *it;
        for (auto jt = rec->get_interaction_record_iterator(); !rec->is_at_end(jt); ++jt)
        {
            cosine = (*jt)->direction;
            cos_x[n] = cosine[0];
            cos_y[n] = cosine[1];
            cos_z[n] = cosine[2];
            ++n;
        }
    }

    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_elementmap(st_context_v2_t pcxt, uint_fast64_t *element_map)
{
    /* *element_map must be allocated by caller using st_num_intersections */
    CONTEXT(pcxt);
    RESULT(cxt);
    
    if (cxt->report_level == RunnerStatistics::GROUPED_COUNTS)
        return st_return_code::RESULT_NOT_REPORTED;

    uint_fast64_t n = 0;
    for (auto it = result->get_ray_record_iterator(); !result->is_at_end(it); ++it)
    {
        ray_record_ptr rec = *it;
        for (auto jt = rec->get_interaction_record_iterator(); !rec->is_at_end(jt); ++jt)
        {
            element_map[n] = (*jt)->element;
            ++n;
        }
    }

    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_stagemap(st_context_v2_t pcxt, uint_fast64_t *stage_map)
{
    /* *stage_map must be allocated by caller using st_num_intersections */
    CONTEXT(pcxt);
    RESULT(cxt);
    
    if (cxt->report_level == RunnerStatistics::GROUPED_COUNTS)
        return st_return_code::RESULT_NOT_REPORTED;

    uint_fast64_t n = number_of_interactions(result);
    // deprecating stages return array of 0s
    for (uint_fast64_t i = 0; i < n; ++i)
        stage_map[i] = 0;

    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_raynumbers(st_context_v2_t pcxt, uint_fast64_t *ray_numbers)
{
    /* *ray_numbers must be allocated by caller using st_num_intersections */
    CONTEXT(pcxt);
    RESULT(cxt);
    
    if (cxt->report_level == RunnerStatistics::GROUPED_COUNTS)
        return st_return_code::RESULT_NOT_REPORTED;

    uint_fast64_t n = 0;
    for (auto it = result->get_ray_record_iterator(); !result->is_at_end(it); ++it)
    {
        ray_record_ptr rec = *it;
        for (auto jt = rec->get_interaction_record_iterator(); !rec->is_at_end(jt); ++jt)
        {
            ray_numbers[n] = rec->id;
            ++n;
        }
    }

    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_sun_stats(st_context_v2_t pcxt,
                                  double 		  *width,
								  double 		  *height,
								  double 		  *area,
								  uint_fast64_t	  *nsunrays)
{
    // TODO: compute power per ray
    CONTEXT(pcxt);
    RESULT(cxt);

    *nsunrays = result->get_sun_ray_count();

    if (cxt->runner_type == st_runner_type_t::OPTIX)
    {
        // optix only sets sun area
        *area = result->get_sun_A_box();
    }
    else
    {
        result->get_sun_dimensions(*width, *height);
        *area = (*width) * (*height);
    }

    return st_return_code::SUCCESS;
}

STAPI_V2 st_return_t st_get_results_data(st_context_v2_t pcxt, args_results_data *data)
{
    /* all arrays in data must be allocated by caller using st_num_intersections */
    CONTEXT(pcxt);
    RESULT(cxt);
    
    if (cxt->report_level == RunnerStatistics::GROUPED_COUNTS)
        return st_return_code::RESULT_NOT_REPORTED;

    uint_fast64_t n = 0;
    glm::dvec3 loc;
    glm::dvec3 cosine;
    for (auto it = result->get_ray_record_iterator(); !result->is_at_end(it); ++it)
    {
        ray_record_ptr rec = *it;
        for (auto jt = rec->get_interaction_record_iterator(); !rec->is_at_end(jt); ++jt)
        {
            loc = (*jt)->location;
            cosine = (*jt)->direction;
            
            data->loc_x[n]       = loc[0];
            data->loc_y[n]       = loc[1];
            data->loc_z[n]       = loc[2];
            data->cos_x[n]       = cosine[0];
            data->cos_y[n]       = cosine[1];
            data->cos_z[n]       = cosine[2];
            data->element_map[n] = (*jt)->element;
            data->stage_map[n]   = 0;
            data->ray_numbers[n] = rec->id;
            ++n;
        }
    }

    return st_return_code::SUCCESS;
}

// TODO: make a way to get group results (number too)
//       make checks positive

/////////////////////////
// Batch api call work //
/////////////////////////

// could upgrade by moving non macro stuff to internal _st_* functions and fetch 
// full context at the beginning of the for loop

bool is_error(st_return_t code)
{
    return code > st_return_code::SUCCESS && code < st_return_code::WARNING_FELLBACK_FROM_EMBREE;
}

STAPI_V2 st_return_t st_batch(st_context_v2_t  pcxt,
							//   st_api_call_args **arguments, 
							  void             **arguments, 
							  st_uint_t        num_calls,
                              st_uint_t        *fail_iteration,
                              bool             verbose)
{
    CONTEXT(pcxt);

    st_return_t code = st_return_code::SUCCESS;
    
    // wrapping in try catch for runtime errors
    try {
        // for each pair of function and arguments, determine which function to call
        for (int i = 0; i < num_calls; ++i) {
            // cast the generic void* argument pointer back to its real type
            // so we can inspect the tag and reach the right payload
            st_api_call_args *call_args = (st_api_call_args*)arguments[i];
            
            if (verbose && cxt->p_cb) 
            {
                std::string msg = "batch call " + std::to_string(i);
                // TODO: magic_enum
                cxt->p_cb(msg.data(), std::to_string(call_args->type).c_str());
            }
            
            // based on the type of call, cast the generic function pointer
            // to the specific function pointer and signature.
            switch (call_args->type) {
        		// Simlulation Data Functions
                // functions for simulation data management thru json strings
                case st_api_call::CALL_ST_READ_INPUT_JSON:
                {
                    code = st_read_input_json(pcxt, call_args->payload.read_input_json_args.json);
                    break;
                }
                case st_api_call::CALL_ST_READ_INPUT_JSON_FILE:
                {
                    code = st_read_input_json_file(pcxt, call_args->payload.read_input_json_file_args.filename);
                    break;
                }
                // functions for simulation data management directly
                case st_api_call::CALL_ST_SET_SIMULATION_PARAMETERS:
                {
                    code = st_set_simulation_parameters(pcxt,
                                                        call_args->payload.set_simulation_parameters_args.params);
                    break;
                }
	            case st_api_call::CALL_ST_SIM_PARAMS:
                {
                    code = st_sim_params(pcxt,
                                         call_args->payload.sim_params_args.raycount,
                                         call_args->payload.sim_params_args.maxcount,
                                         call_args->payload.sim_params_args.include_dynamic_group);
                    break;
                }
	            case st_api_call::CALL_ST_SIM_RAYS:
                {
                    code = st_sim_rays(pcxt,
                                       call_args->payload.sim_rays_args.raycount,
                                       call_args->payload.sim_rays_args.maxcount);
                    break;
                }
	            case st_api_call::CALL_ST_SIM_POWER_TOWER:
                {
                    code = st_sim_power_tower(pcxt,
                                              call_args->payload.sim_power_tower_args.is_power_tower);
                    break;
                }
	            case st_api_call::CALL_ST_SIM_ERRORS:
                {
                    code = st_sim_errors(pcxt,
                                         call_args->payload.sim_errors_args.include_sun_shape,
                                         call_args->payload.sim_errors_args.include_optics);
                    break;
                }
                case st_api_call::CALL_ST_SIM_LOCATION:
                {
                    code = st_sim_location(pcxt,
                                         call_args->payload.sim_location_args.latitude,
                                         call_args->payload.sim_location_args.longitude);
                    break;
                }
                case st_api_call::CALL_ST_SIM_TOLERANCE:
                {
                    code = st_sim_tolerance(pcxt,
                                            call_args->payload.sim_tolerance_args.tolerance);
                    break;
                }
                // functions to add/remove optical properties
                case st_api_call::CALL_ST_NUM_OPTICS:
                {
                    code = st_num_optics(pcxt, call_args->payload.num_optics_args.num_optics);
                    break;
                }
                case st_api_call::CALL_ST_ADD_OPTICAL_PROPERIES_SET:
                {
                    code = st_add_optical_properties_set(pcxt,
                                                         call_args->payload.add_optical_properties_set_args.opt_set,
                                                         call_args->payload.add_optical_properties_set_args.front,
                                                         call_args->payload.add_optical_properties_set_args.back,
                                                         call_args->payload.add_optical_properties_set_args.optic_id);
                    break;
                }
                case st_api_call::CALL_ST_DELETE_OPTIC:
                {
                    code = st_delete_optic(pcxt, call_args->payload.delete_optic_args.idx);
                    break;
                }
                case st_api_call::CALL_ST_CLEAR_OPTICS:
                {
                    code = st_clear_optics(pcxt);
                    break;
                }
                // functions to add/remove elements
                case st_api_call::CALL_ST_NUM_ELEMENTS:
                {
                    code = st_num_elements(pcxt, call_args->payload.num_elements_args.num_elements);
                    break;
                }
	            case st_api_call::CALL_ST_ADD_ELEMENT:
                {
                    code = st_add_element(pcxt,
                                          call_args->payload.add_element_args.args,
                                          call_args->payload.add_element_args.opt_id,
                                          call_args->payload.add_element_args.geom,
                                          call_args->payload.add_element_args.element_id);
                    break;
                }
	            case st_api_call::CALL_ST_DELETE_ELEMENT:
                {
                    code = st_delete_element(pcxt, call_args->payload.delete_element_args.idx);
                    break;
                }
	            case st_api_call::CALL_ST_CLEAR_ELEMENTS:
                {
                    code = st_clear_elements(pcxt);
                    break;
                }
	            // functions to modify elements
	            case st_api_call::CALL_ST_ELEMENT_ENABLED:
                {
                    code = st_element_enabled(pcxt,
                                              call_args->payload.element_enabled_args.idx,
                                              call_args->payload.element_enabled_args.enabled_flag);
                    break;
                }
	            case st_api_call::CALL_ST_ELEMENT_VIRTUAL:
                {
                    code = st_element_virtual(pcxt,
                                              call_args->payload.element_virtual_args.idx,
                                              call_args->payload.element_virtual_args.virtual_flag);
                    break;
                }
	            case st_api_call::CALL_ST_ELEMENT_XYZ:
                {
                    code = st_element_xyz(pcxt,
                                          call_args->payload.element_xyz_args.idx,
                                          call_args->payload.element_xyz_args.x,
                                          call_args->payload.element_xyz_args.y,
                                          call_args->payload.element_xyz_args.z);
                    break;
                }
	            case st_api_call::CALL_ST_ELEMENT_AIM:
                {
                    code = st_element_aim(pcxt,
                                          call_args->payload.element_aim_args.idx,
                                          call_args->payload.element_aim_args.ax,
                                          call_args->payload.element_aim_args.ay,
                                          call_args->payload.element_aim_args.az);
                    break;
                }
	            case st_api_call::CALL_ST_ELEMENT_ZROT:
                {
                    code = st_element_zrot(pcxt,
                                           call_args->payload.element_zrot_args.idx,
                                           call_args->payload.element_zrot_args.zrot);
                    break;
                }
	            case st_api_call::CALL_ST_ELEMENT_APERTURE:
                {
                    code = st_element_aperture(pcxt,
                                               call_args->payload.element_aperture_args.idx,
                                               call_args->payload.element_aperture_args.ap,
                                               call_args->payload.element_aperture_args.params);
                    break;
                }
	            case st_api_call::CALL_ST_ELEMENT_SURFACE:
                {
                    code = st_element_surface(pcxt,
                                              call_args->payload.element_surface_args.idx,
                                              call_args->payload.element_surface_args.surf,
                                              call_args->payload.element_surface_args.params);
                    break;
                }
	            case st_api_call::CALL_ST_ELEMENT_GEOMETRY:
                {
                    code = st_element_geometry(pcxt,
                                               call_args->payload.element_geometry_args.idx,
                                               call_args->payload.element_geometry_args.geom);
                    break;
                }
	            case st_api_call::CALL_ST_ELEMENT_OPTIC:
                {
                    code = st_element_optic(pcxt,
                                            call_args->payload.element_optic_args.idx,
                                            call_args->payload.element_optic_args.opt_id);
                    break;
                }
	            case st_api_call::CALL_ST_ELEMENT_GROUP:
                {
                    code = st_element_group(pcxt,
                                            call_args->payload.element_group_args.idx,
                                            call_args->payload.element_group_args.group);
                    break;
                }
                // sun functions
                case st_api_call::CALL_ST_ADD_SUN:
                {
                    code = st_add_sun(pcxt,
                                      call_args->payload.add_sun_args.args,
                                      call_args->payload.add_sun_args.angle,
                                      call_args->payload.add_sun_args.intensity);
                    break;
                }
                case st_api_call::CALL_ST_SUN_SHAPE:
                {
                    code = st_sun_shape(pcxt,
                                        call_args->payload.sun_shape_args.shape,
                                        call_args->payload.sun_shape_args.sigma_halfwidth_csr);
                    break;
                }
                case st_api_call::CALL_ST_SUN_XYZ:
                {
                    code = st_sun_xyz(pcxt,
                                      call_args->payload.sun_xyz_args.x,
                                      call_args->payload.sun_xyz_args.y,
                                      call_args->payload.sun_xyz_args.z);
                    break;
                }
                case st_api_call::CALL_ST_SUN_POSITION:
                {
                    code = st_sun_position(pcxt,
                                           call_args->payload.sun_position_args.lat,
                                           call_args->payload.sun_position_args.day,
                                           call_args->payload.sun_position_args.hour,
                                           call_args->payload.sun_position_args.x,
                                           call_args->payload.sun_position_args.y,
                                           call_args->payload.sun_position_args.z);
                    break;
                }
                case st_api_call::CALL_ST_SUN_USERDATA:
                {
                    code = st_sun_userdata(pcxt,
                                           call_args->payload.sun_userdata_args.npoints,
                                           call_args->payload.sun_userdata_args.angle,
                                           call_args->payload.sun_userdata_args.intensity);
                    break;
                }
                case st_api_call::CALL_ST_GET_SUN_AZ_ZEN:
                {
                    code = st_get_sun_az_zen(pcxt,
                                             call_args->payload.get_sun_az_zen_args.calc,
                                             call_args->payload.get_sun_az_zen_args.loc,
                                             call_args->payload.get_sun_az_zen_args.dt,
                                             call_args->payload.get_sun_az_zen_args.azimuth,
                                             call_args->payload.get_sun_az_zen_args.zenith);
                    break;
                }
                case st_api_call::CALL_ST_GET_SUN_AZ_EL:
                { 
                    code = st_get_sun_az_el(pcxt,
                                            call_args->payload.get_sun_az_el_args.calc,
                                            call_args->payload.get_sun_az_el_args.loc,
                                            call_args->payload.get_sun_az_el_args.dt,
                                            call_args->payload.get_sun_az_el_args.azimuth,
                                            call_args->payload.get_sun_az_el_args.elevation);
                    break;
                }
                case st_api_call::CALL_ST_GET_SUN_VECTOR:
                {
                    code = st_get_sun_vector(pcxt,
                                             call_args->payload.get_sun_vector_args.calc,
                                             call_args->payload.get_sun_vector_args.loc,
                                             call_args->payload.get_sun_vector_args.dt,
                                             call_args->payload.get_sun_vector_args.sun_x,
                                             call_args->payload.get_sun_vector_args.sun_y,
                                             call_args->payload.get_sun_vector_args.sun_z);
                    break;
                }
                // functions for writing input files for SolTrace
                case st_api_call::CALL_ST_EXPORT_JSON_FILE:
                {
                    code = st_export_json_file(pcxt, call_args->payload.export_json_file_args.filename);
                    break;
                }
                // Simlulation Runner Functions
                case st_api_call::CALL_ST_SIM_SETUP:
                {
                    code = st_sim_setup(pcxt, 
                                        call_args->payload.sim_setup_args.runner_type,
                                        call_args->payload.sim_setup_args.num_threads,
                                        call_args->payload.sim_setup_args.seeds,
                                        call_args->payload.sim_setup_args.num_seeds);
                    break;
                }
                case st_api_call::CALL_ST_SIM_RUN_V2:
                {
                    code = st_sim_run_v2(pcxt);
                    break;
                }
                case st_api_call::CALL_ST_SIM_REPORT:
                {
                    code = st_sim_report(pcxt, call_args->payload.sim_report_args.level);
                    break;
                }
                // Simlulation Results Functions
                case st_api_call::CALL_ST_WRITE_RESULTS_CSV:
                {
                    code = st_write_results_csv(pcxt,
                                                call_args->payload.write_results_csv_args.filename,
                                                call_args->payload.write_results_csv_args.precision);
                    break;
                }
                // functions to get results directly
                case st_api_call::CALL_ST_NUM_INTERSECTIONS:
                {
                    code = st_num_intersections(pcxt,
                                                call_args->payload.num_intersections_args.num_intersections);
                    break;
                }
                case st_api_call::CALL_ST_LOCATIONS:
                {
                    code = st_locations(pcxt,
                                        call_args->payload.locations_args.loc_x,
                                        call_args->payload.locations_args.loc_y,
                                        call_args->payload.locations_args.loc_z);
                    break;
                }
                case st_api_call::CALL_ST_COSINES:
                {
                    code = st_cosines(pcxt,
                                      call_args->payload.cosines_args.cos_x,
                                      call_args->payload.cosines_args.cos_y,
                                      call_args->payload.cosines_args.cos_z);
                    break;
                }
                case st_api_call::CALL_ST_ELEMENTMAP:
                {
                    code = st_elementmap(pcxt,
                                         call_args->payload.elementmap_args.element_map);
                    break;
                }
                case st_api_call::CALL_ST_STAGEMAP:
                {
                    code = st_stagemap(pcxt,
                                       call_args->payload.stagemap_args.stage_map);
                    break;
                }
                case st_api_call::CALL_ST_RAYNUMBERS:
                {
                    code = st_raynumbers(pcxt,
                                         call_args->payload.raynumbers_args.ray_numbers);
                    break;
                }
                case st_api_call::CALL_ST_SUN_STATS:
                {
                    code = st_sun_stats(pcxt,
                                        call_args->payload.sun_stats_args.width,
                                        call_args->payload.sun_stats_args.height,
                                        call_args->payload.sun_stats_args.area,
                                        call_args->payload.sun_stats_args.nsunrays);
                    break;
                }
                case st_api_call::CALL_ST_GET_RESULTS_DATA:
                {
                    code = st_get_results_data(pcxt,
                                               call_args->payload.get_results_data_args.data);
                    break;
                }
                default:
                    code = st_return_code::UKNOWN_BATCH_API_CALL_FAILURE;
                    break;
            }

            if (is_error(code))
            {
                *fail_iteration = i;
                break;
            }
        }
    } catch (const std::exception& e) {
        if (cxt->p_cb) cxt->p_cb("st_batch", e.what());
    }
    return code;
}

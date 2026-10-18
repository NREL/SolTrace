/*
STCORE_API: version 2
functions for interacting with new SimulationData/Runner/Results structure through json

(x: done, -: skipped, i: in progress, T: have TODO, blank: not done)

function name			    | stapi_v2.h/cpp | gtested | stapi_v2.py | unittested | h/cpp batch | py batch | unittested | legacy.py | unittested 
-------------------------------------------------------------------------------------------------------------------------------------------------
 st_create_context 					[x]			 [ ]		[x]			   [ ]			[-]			[ ]			[-]			[-]			  [ ]
 st_reset_context 					[x]			 [ ]		[x]			   [ ]			[-]			[ ]			[-]			[-]			  [ ]
 st_free_context 					[x]			 [ ]		[x]			   [ ]			[-]			[ ]			[-]			[-]			  [ ]
 st_read_input_json 				[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_set_simulation_parameters 		[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_sim_params 						[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[x]			  [ ]
 st_sim_errors 						[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[x]			  [ ]
 st_sim_location 					[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_sim_tolerance 					[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_get_simulation_parameters 		[x]			 [x]		[x]			   [x]			[-]			[ ]			[-]			[ ]			  [ ]
 st_num_optics 						[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_add_optical_properties_set 		[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[x]			  [x]
 st_get_optical_properties_set 		[x]			 [x]		[x]			   [x]			[-]			[ ]			[-]			[ ]			  [ ]
 st_delete_optic 					[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_clear_optics 					[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_num_elements 					[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_add_element 					[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[x]			  [x]
 st_get_element 					[x]			 [x]		[x]			   [x]			[-]			[ ]			[-]			[ ]			  [ ]
 st_delete_element 					[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_clear_elements 					[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_element_enabled 				[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_element_virtual 				[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_element_xyz 					[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_element_aim 					[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_element_zrot 					[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_element_aperture 				[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_element_surface 				[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_element_optic 					[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_element_group 					[x]			 [x]		[ ]			   [ ]			[ ]			[ ]			[ ]			[ ]			  [ ]
 st_add_sun 						[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[x]			  [x]
 st_get_sun 						[T]			 [x]		[x]			   [x]			[-]			[ ]			[-]			[ ]			  [ ]
 st_sun_shape 						[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_sun_xyz 						[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_sun_position 					[T]			 [ ]		[x]			   [ ]			[x]			[ ]			[x]			[ ]			  [ ]
 st_sun_userdata 					[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_get_sun_az_zen 					[x]			 [x]		[x]			   [x]			[ ]			[ ]			[ ]			[ ]			  [ ]
 st_get_sun_az_el 					[x]			 [x]		[x]			   [x]			[ ]			[ ]			[ ]			[ ]			  [ ]
 st_get_sun_vector 					[x]			 [x]		[x]			   [x]			[ ]			[ ]			[ ]			[ ]			  [ ]
 st_get_installed_runners 			[x]			 [x]		[x]			   [ ]			[-]			[ ]			[-]			[ ]			  [ ]
 st_is_runner_installed 			[x]			 [x]		[x]			   [ ]			[-]			[ ]			[-]			[ ]			  [ ]
 st_sim_setup 						[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[x]			  [ ]
 st_sim_run_v2 						[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[x]			  [ ]
 st_sim_report 						[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_write_results_csv 				[x]			 [x]		[x]			   [ ]			[x]			[ ]			[x]			[ ]			  [ ]
 st_write_group_results_json		[x]			 [x]		[ ]			   [ ]			[ ]			[ ]			[ ]			[ ]			  [ ]
 st_num_intersections				[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[ ]			  [ ]
 st_locations 	     				[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[-]			  [ ]
 st_cosines 	     				[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[-]			  [ ]
 st_elementmap 	     				[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[-]			  [ ]
 st_stagemap 	     				[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[-]			  [ ]
 st_raynumbers 	     				[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[-]			  [ ]
 st_sun_stats 	     				[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[x]			  [ ]
 st_get_results_data   				[x]			 [x]		[x]			   [x]			[x]			[ ]			[x]			[x]			  [ ]
 st_calc_zrot_azel   				[ ]			 [ ]		[ ]			   [ ]			[ ]			[ ]			[ ]			[ ]			  [ ]
 st_calc_euler_angles 				[ ]			 [ ]		[ ]			   [ ]			[ ]			[ ]			[ ]			[ ]			  [ ]
 st_transform_to_local 				[ ]			 [ ]		[ ]			   [ ]			[ ]			[ ]			[ ]			[ ]			  [ ]
 st_transform_to_reference 			[ ]			 [ ]		[ ]			   [ ]			[ ]			[ ]			[ ]			[ ]			  [ ]
 st_matrix_vector_mult 				[ ]			 [ ]		[ ]			   [ ]			[ ]			[ ]			[ ]			[ ]			  [ ]
 st_calc_transform_matrices 		[ ]			 [ ]		[ ]			   [ ]			[ ]			[ ]			[ ]			[ ]			  [ ]
 st_matrix_transpose 				[ ]			 [ ]		[ ]			   [ ]			[ ]			[ ]			[ ]			[ ]			  [ ]
 st_batch 							[i,T]		 [ ]		[x]			   [ ]			[ ]			[ ]			[ ]			[ ]			  [ ]
-------------------------------------------------------------------------------------------------------------------------------------------------
other items
-------------------------------------------------------------------------------------------------------------------------------------------------
 STAPIv2.__check_return_code  		[-]			 [-]		[x]			   [x]			[ ]			[ ]			[ ]			[ ]			  [ ]
 st_return_code error messages		[-]			 [-]		[x]			   [x]			[ ]			[ ]			[ ]			[ ]			  [ ]
 st_return_code warning messages	[-]			 [-]		[x]			   [x]			[ ]			[ ]			[ ]			[ ]			  [ ]
 STAPIv2.generate_api_call  		[-]			 [-]		[x]			   [i]			[ ]			[ ]			[ ]			[ ]			  [ ]
 legacy.write_soltrace_input_file	[-]			 [-]		[-]			   [-]			[-]			[-]			[-]			[x]			  [ ]
 tower demo					  		[-]			 [-]		[-]			   [-]			[-]			[-]			[-]			[-]			  [i]
 heliostat template					[-]			 [-]		[-]			   [-]			[-]			[-]			[-]			[-]			  [ ]
 linear fresnel template			[-]			 [-]		[-]			   [-]			[-]			[-]			[-]			[-]			  [ ]
 parabolic trough template			[-]			 [-]		[-]			   [-]			[-]			[-]			[-]			[-]			  [ ]
 parabolic dish template			[-]			 [-]		[-]			   [-]			[-]			[-]			[-]			[-]			  [ ]
 point.__repr__						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__bool__						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__float__					[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__int__						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__iter__						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__get_set_item__ 			[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__len__						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__negation__					[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__abs__						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__add__						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__radd__						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__iadd__						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__sub__						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__rsub__						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__isub__						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__mul__						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__rmul__						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__imul__						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__floordiv__					[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__ifloordiv__				[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__truediv__					[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__itruediv__					[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__matmul__					[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__imatmul__					[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__eq__						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__neq__						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__lt__						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__gt__						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__lte__						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.__gte__						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.reduce						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.radius						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.unitize						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 point.as_list						[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
 Point.from_list					[-]			 [-]		[x]			   [x]			[-]			[-]			[-]			[-]			  [-]
-------------------------------------------------------------------------------------------------------------------------------------------------
 */

#ifndef __soltraceapi_v2_h
#define __soltraceapi_v2_h

#ifdef _STAPI_V2_DLL_
	#ifdef STAPI_V2_EXPORTS
	#define STAPI_V2 __declspec(dllexport)
	#else
	#define STAPI_V2 __declspec(dllimport)
	#endif
#else
	#define STAPI_V2
#endif

#ifndef STAPI_V2
#define STAPI_V2
#endif

#include <filesystem>

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

#ifdef __cplusplus
extern "C" {
#endif

#ifndef M_PI
	#define M_PI 3.141592653589793238462643
#endif


#define ST_WRAP_CB_TRY_CATCH(call, cb) 	  \
	try { (call); } 					  \
	catch (const std::exception& e) {     \
		if (cb) cb(#call, e.what()); 	  \
		return st_return_code::EXCEPTION; \
	} 									  \

namespace fs = std::filesystem;

using SolTrace::Runner::SimulationRunner;
using SolTrace::Runner::RunnerStatus;
using SolTrace::Runner::RunnerStatistics;
using SolTrace::NativeRunner::NativeRunner;
using SolTrace::Data::Aperture;
using SolTrace::Data::single_element_ptr;
using spcm = SolTrace::Data::SolarPositionCalculationMethod;

typedef uint32_t st_uint_t;

/* changing return code convention from v1
   to be in line with industry convention.
   i.e. 0 for success, non-zero for failure. */
typedef st_uint_t st_return_t;

typedef enum st_return_code : st_return_t {
	SUCCESS = 0,
	FAILURE,
	CANCEL,
	CONTEXT_NOT_FOUND,
	DATA_NOT_FOUND,
	RUNNER_NOT_FOUND,
	RESULT_NOT_FOUND,
	INVALID_ARGUMENTS,
	DATA_INSERTION_FAILURE,
	DATA_VALUE_NOT_FOUND,
	RUNNER_INILIALIZE_FAILURE,
	RUNNER_NUMBER_THREADS_SEEDS_MISMATCH_FAILURE,
	RUNNER_SETUP_ERROR,
	RUNNER_SETUP_FAILURE,
	RUNNER_NOT_READY_TO_RUN,
	RUNNER_NOT_READY_TO_REPORT,
	RESULT_NOT_REPORTED,
	EXCEPTION,
	UKNOWN_BATCH_API_CALL_FAILURE,

	WARNING_FELLBACK_FROM_EMBREE,
	WARNING_FELLBACK_FROM_OPTIX,
	WARNING_ARGUMENT_IGNORED_BY_RUNNER,
	WARNING_SUN_SHAPE_IGNORED,
	WARNING_GROUP_IGNORED,
	WARNING_NOT_FOUND,

	RETURN_COUNT /* sentinel (not a valid return type) */
} st_return_code;

#ifndef DEFAULT_NUM_THREADS
#define DEFAULT_NUM_THREADS 8
#endif

typedef enum st_runner_type_t : st_uint_t {
	NATIVE = 0,   /* 0 */
	OPTIX,        /* 1 */
	EMBREE,		  /* 2 */
	RUNNER_COUNT         /* sentinel (not a valid runner) */
} st_runner_type_t;

typedef int (*p_callback)(const char* loc, const char* msg);

typedef struct st_context {
	SimulationData*   p_data	   = nullptr;
	st_runner_type_t  runner_type  = st_runner_type_t::RUNNER_COUNT;
	SimulationRunner* p_runner	   = nullptr;
	RunnerStatistics  report_level = RunnerStatistics::STATISTICS_COUNT;
	SimulationResult* p_results	   = nullptr;
	p_callback		  p_cb	   	   = nullptr;
	fs::path		  location;
} st_context;

typedef void* st_context_v2_t;

////////////////////////////////
// SolTrace Context Functions //
////////////////////////////////

// functions for SolTrace context management
STAPI_V2 st_return_t st_create_context(st_context_v2_t* pcxt, p_callback cb = nullptr);
STAPI_V2 st_return_t st_reset_context(st_context_v2_t pcxt);
STAPI_V2 st_return_t st_free_context(st_context_v2_t pcxt);
STAPI_V2 st_return_t st_locate_context(st_context_v2_t pcxt, const char *path);

////////////////////////////////
// Simlulation Data Functions //
////////////////////////////////

// functions for SolTrace data management
// functions for simulation data management thru json strings
STAPI_V2 st_return_t st_read_input_json(st_context_v2_t pcxt, const char *json);
STAPI_V2 st_return_t st_read_input_json_file(st_context_v2_t pcxt, const char *filename);

// functions for simulation data management directly
typedef struct args_simulation_parameters {
	uint_fast64_t number_of_rays;
	uint_fast64_t max_number_of_rays;
	double 		  tolerance;
	double 		  latitude;
	double 		  longitude;
	bool 		  include_sun_shape_errors;
	bool 		  include_optical_errors;
	bool 		  as_power_tower;
} args_simulation_parameters;
STAPI_V2 st_return_t st_set_simulation_parameters(st_context_v2_t pcxt, args_simulation_parameters *params);
STAPI_V2 st_return_t st_sim_params(st_context_v2_t pcxt,
								   uint_fast64_t   raycount,
								   uint_fast64_t   maxcount,
								   bool			   include_dynamic_group);
STAPI_V2 st_return_t st_sim_rays(st_context_v2_t pcxt,
								 uint_fast64_t   raycount,
								 uint_fast64_t   maxcount);
STAPI_V2 st_return_t st_sim_power_tower(st_context_v2_t pcxt,
								        bool 		    is_power_tower);
STAPI_V2 st_return_t st_sim_errors(st_context_v2_t pcxt,
								   bool			   include_sun_shape,
								   bool			   include_optics);
STAPI_V2 st_return_t st_sim_location(st_context_v2_t pcxt,
								     double		     latitude,
								     double		     longitude);
STAPI_V2 st_return_t st_sim_tolerance(st_context_v2_t pcxt, double tolerance);
STAPI_V2 st_return_t st_get_simulation_parameters(st_context_v2_t pcxt, args_simulation_parameters *params);

// functions to add/remove/set optical properties
STAPI_V2 st_return_t st_num_optics(st_context_v2_t pcxt, uint_fast64_t *num_optics);
typedef struct args_optical_properties_face {
	double transmissivity;
	double reflectivity;
	double slope_error;
	double specularity_error;
	char   error_distribution_type;
} args_optical_properties_face;
typedef struct args_optical_properties_set {
	const char *name;
	double 	   refraction_index_front;
	double 	   refraction_index_back;
	st_uint_t  type; // 1 = refraction, otherwise reflection
} args_optical_properties_set;
STAPI_V2 st_return_t st_add_optical_properties_set(st_context_v2_t 				pcxt, 
												   args_optical_properties_set 	*opt_set,
												   args_optical_properties_face *front, 
												   args_optical_properties_face *back, 
												   uint_fast64_t 				*optic_id);
STAPI_V2 st_return_t st_get_optical_properties_set(st_context_v2_t 				pcxt, 
												   uint_fast64_t 				optic_id,
												   args_optical_properties_set 	*opt_set,
												   args_optical_properties_face *front, 
												   args_optical_properties_face *back);
STAPI_V2 st_return_t st_delete_optic(st_context_v2_t pcxt, st_uint_t idx);
STAPI_V2 st_return_t st_clear_optics(st_context_v2_t pcxt);

// functions to add/remove elements
STAPI_V2 st_return_t st_num_elements(st_context_v2_t pcxt, uint_fast64_t *num_elements);
typedef struct args_element {
	double  x;
	double  y;
	double  z;
	double  ax;
	double  ay;
	double  az;
	double  zrot;
	bool    enabled_flag;
	bool    virtual_flag;
	int32_t group = -1;
} args_element;
typedef struct args_geometry {
	double a_params[8];
	double s_params[8];
	char   ap;
	char   surf;
} args_geometry;
STAPI_V2 st_return_t st_add_element(st_context_v2_t pcxt,
									args_element    *args,
									int_fast64_t    opt_id,
									args_geometry	*geom,
									uint_fast64_t   *element_id);
STAPI_V2 st_return_t st_get_element(st_context_v2_t pcxt,
									uint_fast64_t   element_id,
									args_element    *args,
									int_fast64_t    *opt_id,
									args_geometry	*geom);
STAPI_V2 st_return_t st_delete_element(st_context_v2_t pcxt, st_uint_t idx);
STAPI_V2 st_return_t st_clear_elements(st_context_v2_t pcxt);
// functions to modify elements
STAPI_V2 st_return_t st_element_enabled(st_context_v2_t pcxt,
										st_uint_t 		idx,
										bool 			enabled_flag);
STAPI_V2 st_return_t st_element_virtual(st_context_v2_t pcxt,
										st_uint_t 		idx,
										bool 			virtual_flag);
STAPI_V2 st_return_t st_element_xyz(st_context_v2_t pcxt, 
									st_uint_t 		idx,
									double 	  		x,
									double 	  		y,
									double 	  		z);
STAPI_V2 st_return_t st_element_aim(st_context_v2_t pcxt, 
									st_uint_t 		idx,
									double 	  		ax,
									double 	  		ay,
									double 	  		az);
STAPI_V2 st_return_t st_element_zrot(st_context_v2_t pcxt,
									 st_uint_t 		 idx,
									 double 		 zrot);
STAPI_V2 st_return_t st_element_aperture(st_context_v2_t pcxt,
                                         st_uint_t 	     idx,
                                         char      	     ap,
                                         double    	     params[8]);
STAPI_V2 st_return_t st_element_surface(st_context_v2_t pcxt,
                                        st_uint_t 	    idx,
                                        char      	    surf,
                                        double    	    params[8]);
STAPI_V2 st_return_t st_element_geometry(st_context_v2_t pcxt,
                                         st_uint_t 	     idx,
                                         args_geometry	 *geom);
STAPI_V2 st_return_t st_element_optic(st_context_v2_t pcxt,
									  st_uint_t 	  idx,
									  int_fast64_t 	  opt_id);
STAPI_V2 st_return_t st_element_group(st_context_v2_t pcxt,
									  st_uint_t 	  idx,
									  int32_t   	  group);

// sun functions
typedef struct args_sun {
	st_uint_t npoints;
	double x;
	double y;
	double z;
	double sigma_halfwidth_csr;
	char shape;
} args_sun;
STAPI_V2 st_return_t st_add_sun(st_context_v2_t pcxt,
								args_sun 		*args,
								double 		 	*angle,
								double 		 	*intensity);
STAPI_V2 st_return_t st_get_sun(st_context_v2_t pcxt,
								args_sun 		*args,
								double 		 	**angle,
								double 		 	**intensity,
								st_uint_t 		npoints = 0);
STAPI_V2 st_return_t st_sun_shape(st_context_v2_t pcxt,
								  char   		  shape, 
								  double 		  sigma_halfwidth_csr);
STAPI_V2 st_return_t st_sun_xyz(st_context_v2_t pcxt,
								double 			x,
								double 			y,
								double 			z);
STAPI_V2 st_return_t st_sun_position(st_context_v2_t pcxt,
									 double  		 lat,
									 double  		 day,
									 double  		 hour,
									 double			 *x,
									 double			 *y,
									 double			 *z);
STAPI_V2 st_return_t st_sun_userdata(st_context_v2_t pcxt,
									 st_uint_t 		 npoints,
									 double 		 *angle,
									 double 		 *intensity);
typedef struct args_sun_location {
	double latitude;
	double longitude;
	double timeZone;
	double altitude = 0.0;
} args_sun_location;
typedef struct args_sun_datetime {
	int year;
	int month;
	int day;
	int hour   = 12;
	int minute = 0;
	int second = 0;
} args_sun_datetime;
STAPI_V2 st_return_t st_get_sun_az_zen(st_context_v2_t   pcxt,
									   int				 calc,
									   args_sun_location *loc,
									   args_sun_datetime *dt,
									   double 			 *azimuth,
									   double 			 *zenith);
STAPI_V2 st_return_t st_get_sun_az_el(st_context_v2_t   pcxt,
									  int				calc,
									  args_sun_location *loc,
									  args_sun_datetime *dt,
									  double 			*azimuth,
									  double 			*elevation);
STAPI_V2 st_return_t st_get_sun_vector(st_context_v2_t   pcxt,
									   int				 calc,
									   args_sun_location *loc,
									   args_sun_datetime *dt,
									   double 			 *sun_x,
									   double 			 *sun_y,
									   double 			 *sun_z);

// functions for writing input files for SolTrace
STAPI_V2 st_return_t st_export_json_file(st_context_v2_t pcxt, const char *filename);

//////////////////////////////////
// Simlulation Runner Functions //
//////////////////////////////////

// functions for SolTrace runner management
// NOTE: following two functions are static ie don't require pcxt
STAPI_V2 st_return_t st_get_installed_runners(uint8_t *installed);
STAPI_V2 st_return_t st_is_runner_installed(st_runner_type_t type,
                                            bool             *installed);
STAPI_V2 st_return_t st_sim_setup(st_context_v2_t  pcxt, 
								  st_runner_type_t runner_type, 
								  uint_fast64_t    num_threads = DEFAULT_NUM_THREADS,
								  st_uint_t 	   *seeds = nullptr,
								  size_t		   num_seeds = 0);
STAPI_V2 st_return_t st_sim_run_v2(st_context_v2_t pcxt);
STAPI_V2 st_return_t st_sim_report(st_context_v2_t pcxt, int level);

///////////////////////////////////
// Simlulation Results Functions //
///////////////////////////////////

// functions for dumping results to file
// TODO: add st_write_results_json once element_groups pr is merged
STAPI_V2 st_return_t st_write_results_csv(st_context_v2_t pcxt, 
										  const char 	  *filename, 
										  int 			  precision = 12);
STAPI_V2 st_return_t st_write_group_results_json(st_context_v2_t pcxt, 
										         const char 	 *filename, 
										         int 			 precision = 12,
                                                 int             indent = 4);

// functions to get results directly
STAPI_V2 st_return_t st_num_intersections(st_context_v2_t pcxt, uint_fast64_t *num_intersections);
STAPI_V2 st_return_t st_locations(st_context_v2_t pcxt,
									 double 		 *loc_x,
									 double 		 *loc_y,
									 double 		 *loc_z);
STAPI_V2 st_return_t st_cosines(st_context_v2_t pcxt,
							 	   double 		   *cos_x,
							 	   double 		   *cos_y,
							 	   double 		   *cos_z);
STAPI_V2 st_return_t st_elementmap(st_context_v2_t pcxt, uint_fast64_t *element_map);
STAPI_V2 st_return_t st_stagemap(st_context_v2_t pcxt, uint_fast64_t *stage_map);
STAPI_V2 st_return_t st_raynumbers(st_context_v2_t pcxt, uint_fast64_t *ray_numbers);
STAPI_V2 st_return_t st_sun_stats(st_context_v2_t pcxt,
								  double 		  *width,
								  double 		  *height,
								  double 		  *area,
								  uint_fast64_t	  *nsunrays);
typedef struct args_results_data {
	double 		  *loc_x;
	double 		  *loc_y;
	double 		  *loc_z;
	double 		  *cos_x;
	double 		  *cos_y;
	double 		  *cos_z;
	uint_fast64_t *element_map;
	uint_fast64_t *stage_map;
	uint_fast64_t *ray_numbers;
} args_results_data;
STAPI_V2 st_return_t st_get_results_data(st_context_v2_t pcxt, args_results_data *data);

/*
create simualtion information -> st_context_v2_t st_create_context_v2();
free simualtion information   -> int st_free_context_v2(st_context_v2_t pcxt):

- simulation parameter setup
	set threads 			  -> int st_num_threads(st_context_v2_t pcxt, int n);
	set verbose 			  -> int st_verbose(st_context_v2_t pcxt, bool v);

- simlulation data set up
	read json 			      -> int st_read_input_json(st_context_v2_t pcxt, const char *json);

  - ray sources set up
	  add sun 			  	  -> int st_sun_v2(st_context_v2_t pcxt, const char *json);

  - optical property sets
	  add optical set  		  -> int st_add_optical_set(st_context_v2_t pcxt, const char *json);
	  add optical face 		  -> int st_add_optical_face(st_context_v2_t pcxt, const char *json);
  
  - no longer continuing support for stages
  
  - elements
	  add element  			  -> int st_add_element_v2(st_context_v2_t pcxt, const char *json);
	  add elements 			  -> int st_add_elements_v2(st_context_v2_t pcxt, const char *json);

- simulation runner
	setup  					  -> int st_sim_setup(st_context_v2_t pcxt, st_runner_type_t runner_type);
	run    					  -> int st_sim_run_v2(st_context_v2_t pcxt, st_uint_t seed, int (*callback)(...), void *data);
	report 					  -> int st_sim_report(st_context_v2_t pcxt, int level);

- simulation results
	write csv  				  -> int st_write_results_csv(st_context_v2_t pcxt, const char *filename);
	write json 				  -> int st_write_results_json(st_context_v2_t pcxt, const char *filename);*/

///////////////////////
// Batch caller work //
///////////////////////

// enum defining all calls available to batch together
// NOTE: union payload of st_api_call_args must follow this order
typedef enum st_api_call : st_uint_t {
	// Simlulation Data Functions
	// functions for simulation data management thru json strings
    CALL_ST_READ_INPUT_JSON = 0,
    CALL_ST_READ_INPUT_JSON_FILE,
	// functions for simulation data management directly
	CALL_ST_SET_SIMULATION_PARAMETERS,
	CALL_ST_SIM_PARAMS,
	CALL_ST_SIM_RAYS,
	CALL_ST_SIM_POWER_TOWER,
	CALL_ST_SIM_ERRORS,
	CALL_ST_SIM_LOCATION,
	CALL_ST_SIM_TOLERANCE,
	// functions to add/remove optical properties
	CALL_ST_NUM_OPTICS,
	CALL_ST_ADD_OPTICAL_PROPERIES_SET,
	CALL_ST_DELETE_OPTIC,
	CALL_ST_CLEAR_OPTICS,
	// functions to add/remove elements
    CALL_ST_NUM_ELEMENTS,
	CALL_ST_ADD_ELEMENT,
	CALL_ST_DELETE_ELEMENT,
	CALL_ST_CLEAR_ELEMENTS,
	// functions to modify elements
	CALL_ST_ELEMENT_ENABLED,
	CALL_ST_ELEMENT_VIRTUAL,
	CALL_ST_ELEMENT_XYZ,
	CALL_ST_ELEMENT_AIM,
	CALL_ST_ELEMENT_ZROT,
	CALL_ST_ELEMENT_APERTURE,
	CALL_ST_ELEMENT_SURFACE,
	CALL_ST_ELEMENT_GEOMETRY,
	CALL_ST_ELEMENT_OPTIC,
	CALL_ST_ELEMENT_GROUP,
	// sun functions
	CALL_ST_ADD_SUN,
	CALL_ST_SUN_SHAPE,
	CALL_ST_SUN_XYZ,
	CALL_ST_SUN_POSITION,
	CALL_ST_SUN_USERDATA,
	CALL_ST_GET_SUN_AZ_ZEN,
	CALL_ST_GET_SUN_AZ_EL,
	CALL_ST_GET_SUN_VECTOR,
	// functions for writing input files for SolTrace
	CALL_ST_EXPORT_JSON_FILE,
	// Simlulation Runner Functions
	CALL_ST_SIM_SETUP,
	CALL_ST_SIM_RUN_V2,
	CALL_ST_SIM_REPORT,
	// Simlulation Results Functions
	CALL_ST_WRITE_RESULTS_CSV,
	// functions to get results directly
	CALL_ST_NUM_INTERSECTIONS,
	CALL_ST_LOCATIONS,
	CALL_ST_COSINES,
	CALL_ST_ELEMENTMAP,
	CALL_ST_STAGEMAP,
	CALL_ST_RAYNUMBERS,
	CALL_ST_SUN_STATS,
	CALL_ST_GET_RESULTS_DATA,

	API_CALL_COUNT			// sentinel
} st_api_call;

typedef struct empty_args {} empty_args;

// argument layouts — one per function signature

// Simlulation Data Functions

// functions for simulation data management thru json strings
typedef struct args_st_read_input_json {
	const char *json;
} args_st_read_input_json;

typedef struct args_st_read_input_json_file {
	const char *filename;
} args_st_read_input_json_file;

// functions for simulation data management directly
typedef struct args_st_set_simulation_parameters {
	args_simulation_parameters *params;
} args_st_set_simulation_parameters; 

typedef struct args_st_sim_params {
	uint_fast64_t raycount;
	uint_fast64_t maxcount;
	bool 		  include_dynamic_group;
} args_st_sim_params;

typedef struct args_st_sim_rays {
	uint_fast64_t raycount;
	uint_fast64_t maxcount;
} args_st_sim_rays;

typedef struct args_st_sim_power_tower {
	bool is_power_tower;
} args_st_sim_power_tower;

typedef struct args_st_sim_errors {
	bool include_sun_shape;
	bool include_optics;
} args_st_sim_errors;

typedef struct args_st_sim_location {
	double latitude;
	double longitude;
} args_st_sim_location;

typedef struct args_st_sim_tolerance {
	double tolerance;
} args_st_sim_tolerance;

typedef struct args_st_get_simulation_parameters {
	args_simulation_parameters *params;
} args_st_get_simulation_parameters; 

// functions to add/remove/set optical properties
typedef struct args_st_num_optics {
	uint_fast64_t *num_optics;
} args_st_num_optics;

typedef struct args_st_add_optical_properties_set {
	args_optical_properties_set *opt_set;
	args_optical_properties_face *front;
	args_optical_properties_face *back;
	uint_fast64_t *optic_id;
} args_st_add_optical_properties_set;

typedef struct args_st_get_optical_properties_set {
	uint_fast64_t 				 optic_id;
	args_optical_properties_set  *opt_set;
	args_optical_properties_face *front;
	args_optical_properties_face *back;
}  args_st_get_optical_properties_set;

typedef struct args_st_delete_optic {
	st_uint_t idx;
} args_st_delete_optic;

typedef empty_args args_st_clear_optics;

// functions to add/remove elements
typedef struct args_st_num_elements {
	uint_fast64_t *num_elements;
} args_st_num_elements;

typedef struct args_st_add_element {
	args_element  *args;
	int_fast64_t  opt_id;
	args_geometry *geom;
	uint_fast64_t *element_id;
} args_st_add_element;

typedef struct args_st_get_element {
	uint_fast64_t element_id;
	args_element  *args;
	int_fast64_t  *opt_id;
	args_geometry *geom;
} args_st_get_element;

typedef struct args_st_delete_element {
	st_uint_t idx;
} args_st_delete_element;

typedef empty_args args_st_clear_elements;

// functions to modify elements
typedef struct args_st_element_enabled {
	st_uint_t idx;
	bool 	  enabled_flag;
} args_st_element_enabled;

typedef struct args_st_element_virtual {
	st_uint_t idx;
	bool 	  virtual_flag;
} args_st_element_virtual;

typedef struct args_st_element_xyz {
	st_uint_t idx;
	double    x;
	double    y;
	double    z;
} args_st_element_xyz;

typedef struct args_st_element_aim {
	st_uint_t idx;
	double    ax;
	double    ay;
	double    az;
} args_st_element_aim;

typedef struct args_st_element_zrot {
	st_uint_t idx;
	double 	  zrot;
} args_st_element_zrot;

typedef struct args_st_element_aperture {
	st_uint_t idx;
	char 	  ap;
	double    *params;
} args_st_element_aperture;

typedef struct args_st_element_surface {
	st_uint_t idx;
	char 	  surf;
	double    *params;
} args_st_element_surface;

typedef struct args_st_element_geometry {
	st_uint_t 	  idx;
	args_geometry *geom;
} args_st_element_geometry;

typedef struct args_st_element_optic {
	st_uint_t 	 idx;
	int_fast64_t opt_id;
} args_st_element_optic;

typedef struct args_st_element_group {
	st_uint_t idx;
	int32_t   group;
} args_st_element_group;

// sun functions
typedef struct args_st_add_sun {
	args_sun *args;
	double* angle;
	double* intensity;
} args_st_add_sun;

typedef struct args_st_get_sun {
	args_sun  *args;
	double 	  **angle;
	double 	  **intensity;
	st_uint_t npoints = 0;
} args_st_get_sun;

typedef struct args_st_sun_shape {
	char   shape;
	double sigma_halfwidth_csr;
} args_st_sun_shape;

typedef struct args_st_sun_xyz {
	double x;
	double y;
	double z;
} args_st_sun_xyz;

typedef struct args_st_sun_position {
	double  lat;
	double  day;
	double  hour;
	double* x;
	double* y;
	double* z;
} args_st_sun_position;

typedef struct args_st_sun_userdata {
	st_uint_t npoints;
	double*	  angle;
	double*	  intensity;
} args_st_sun_userdata;

typedef struct args_st_get_sun_az_zen {
	int			  	  calc;
	args_sun_location *loc;
	args_sun_datetime *dt;
	double			  *azimuth;
	double			  *zenith;
} args_st_get_sun_az_zen;

typedef struct args_st_get_sun_az_el {
	int			  	  calc;
	args_sun_location *loc;
	args_sun_datetime *dt;
	double			  *azimuth;
	double			  *elevation;
} args_st_get_sun_az_el;

typedef struct args_st_get_sun_vector {
	int			  	  calc;
	args_sun_location *loc;
	args_sun_datetime *dt;
	double			  *sun_x;
	double			  *sun_y;
	double			  *sun_z;
} args_st_get_sun_vector;

// functions for writing input files for SolTrace
typedef struct args_st_export_json_file {
	const char *filename;
} args_st_export_json_file;

// Simlulation Runner Functions
typedef struct args_st_get_installed_runners {
	uint8_t *installed;
} args_st_get_installed_runners;

typedef struct args_st_is_runner_installed {
	st_runner_type_t type;
	bool             *installed;
} args_st_is_runner_installed;

typedef struct args_st_sim_setup {
	st_runner_type_t runner_type;
	uint_fast64_t    num_threads = DEFAULT_NUM_THREADS;
	st_uint_t 	 	 *seeds = nullptr;
	size_t		     num_seeds = 0;
} args_st_sim_setup;

typedef empty_args args_st_sim_run_v2;

typedef struct args_st_sim_report {
	int level;
} args_st_sim_report;

// Simlulation Results Functions
typedef struct args_st_write_results_csv {
	const char *filename;
	int 	   precision = 12;
} args_st_write_results_csv;

// functions to get results directly
typedef struct args_st_num_intersections {
	uint_fast64_t *num_intersections;
} args_st_num_intersections;

typedef struct args_st_locations {
	double *loc_x;
	double *loc_y;
	double *loc_z;
} args_st_locations;

typedef struct args_st_cosines {
	double *cos_x;
	double *cos_y;
	double *cos_z;
} args_st_cosines;

typedef struct args_st_elementmap {
	uint_fast64_t *element_map;
} args_st_elementmap;

typedef struct args_st_stagemap {
	uint_fast64_t *stage_map;
} args_st_stagemap;

typedef struct args_st_raynumbers {
	uint_fast64_t *ray_numbers;
} args_st_raynumbers;

typedef struct args_st_sun_stats {
	double 		  *width;
	double 		  *height;
	double 		  *area;
	uint_fast64_t *nsunrays;
} args_st_sun_stats;

typedef struct args_st_get_results_data {
	args_results_data *data;
} args_st_get_results_data;

/* Every argument layout passed through the generic arrays is one of
 * these: a tag telling us which payload is active, plus the payload
 * itself. This is what a void* in the "arguments" array actually
 * points to. */
typedef struct st_api_call_args {
	// NOTE: this must follow the order of st_api_call enum
    union {
		// Simlulation Data Functions
		// functions for simulation data management thru json strings
        args_st_read_input_json 	 read_input_json_args;
        args_st_read_input_json_file read_input_json_file_args;
		// functions for simulation data management directly
		args_st_set_simulation_parameters set_simulation_parameters_args;
		args_st_sim_params 	  			  sim_params_args;
		args_st_sim_rays 	  			  sim_rays_args;
		args_st_sim_power_tower  		  sim_power_tower_args;
		args_st_sim_errors 	  			  sim_errors_args;
		args_st_sim_location  			  sim_location_args;
		args_st_sim_tolerance 			  sim_tolerance_args;
		// functions to add/remove/set optical properties
		args_st_num_optics				   num_optics_args;
		args_st_add_optical_properties_set add_optical_properties_set_args;
		args_st_delete_optic			   delete_optic_args;
		args_st_clear_optics			   clear_optics_args;
		// functions to add/remove elements
		args_st_num_elements   num_elements_args;
		args_st_add_element    add_element_args;
		args_st_delete_element delete_element_args;
		args_st_clear_elements clear_elements_args;
		// functions to modify elements
		args_st_element_enabled  element_enabled_args;
		args_st_element_virtual  element_virtual_args;
		args_st_element_xyz      element_xyz_args;
		args_st_element_aim      element_aim_args;
		args_st_element_zrot     element_zrot_args;
		args_st_element_aperture element_aperture_args;
		args_st_element_surface  element_surface_args;
		args_st_element_geometry element_geometry_args;
		args_st_element_optic    element_optic_args;
		args_st_element_group    element_group_args;
		// sun functions
		args_st_add_sun 	   add_sun_args;
		args_st_sun_shape 	   sun_shape_args;
		args_st_sun_xyz 	   sun_xyz_args;
		args_st_sun_position   sun_position_args;
		args_st_sun_userdata   sun_userdata_args;
		args_st_get_sun_az_zen get_sun_az_zen_args;
		args_st_get_sun_az_el  get_sun_az_el_args;
		args_st_get_sun_vector get_sun_vector_args;
		// functions for writing input files for SolTrace
		args_st_export_json_file export_json_file_args;
		// Simlulation Runner Functions
		args_st_sim_setup  sim_setup_args;
		args_st_sim_run_v2 sim_run_v2_args;
		args_st_sim_report sim_report_args;
		// Simlulation Results Functions
		args_st_write_results_csv write_results_csv_args;
		// functions to get results directly
		args_st_num_intersections num_intersections_args;
		args_st_locations		  locations_args;
		args_st_cosines		      cosines_args;
		args_st_elementmap		  elementmap_args;
		args_st_stagemap		  stagemap_args;
		args_st_raynumbers		  raynumbers_args;
		args_st_sun_stats		  sun_stats_args;
		args_st_get_results_data  get_results_data_args;
    } payload;
    st_api_call type;
} st_api_call_args;

/* Generic function pointer type used to store heterogeneous function
 * pointers in the "functions" array. Each one gets cast back to its
 * real signature right before it's invoked. */
typedef st_return_t (*st_api_func_ptr)(void);

/* ------------------------------------------------------------------ */
/* The function containing the loop for batching stapi_v2 calls.      */
/*                                                                    */
/* functions : array of function pointers (stored generically)        */
/* arguments : array of void* pointers, each actually pointing at a   */
/*             st_api_call_args struct                                */
/* count     : number of calls to make                                */
/*                                                                    */
/* For each call, the void* function pointer and void* argument       */
/* pointer are cast back to their real, concrete types inside the 	  */
/* loop, based on the st_api_call tag carried in st_api_call_args.    */
/* Error codes break the loop, and warning codes are ignored.         */
/* ------------------------------------------------------------------ */
STAPI_V2 st_return_t st_batch(st_context_v2_t  pcxt,
							  void 			   **arguments, // TODO: don't need to cast as void because all are st_api_call_args
							  st_uint_t  	   count,
                              st_uint_t		   *fail_iteration,
                              bool 			   verbose);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif

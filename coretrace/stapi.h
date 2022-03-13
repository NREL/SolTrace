
/*******************************************************************************************************
*  Copyright 2018 Alliance for Sustainable Energy, LLC
*
*  NOTICE: This software was developed at least in part by Alliance for Sustainable Energy, LLC
*  ("Alliance") under Contract No. DE-AC36-08GO28308 with the U.S. Department of Energy and the U.S.
*  The Government retains for itself and others acting on its behalf a nonexclusive, paid-up,
*  irrevocable worldwide license in the software to reproduce, prepare derivative works, distribute
*  copies to the public, perform publicly and display publicly, and to permit others to do so.
*
*  Redistribution and use in source and binary forms, with or without modification, are permitted
*  provided that the following conditions are met:
*
*  1. Redistributions of source code must retain the above copyright notice, the above government
*  rights notice, this list of conditions and the following disclaimer.
*
*  2. Redistributions in binary form must reproduce the above copyright notice, the above government
*  rights notice, this list of conditions and the following disclaimer in the documentation and/or
*  other materials provided with the distribution.
*
*  3. The entire corresponding source code of any redistribution, with or without modification, by a
*  research entity, including but not limited to any contracting manager/operator of a United States
*  National Laboratory, any institution of higher learning, and any non-profit organization, must be
*  made publicly available under this license for as long as the redistribution is made available by
*  the research entity.
*
*  4. Redistribution of this software, without modification, must refer to the software by the same
*  designation. Redistribution of a modified version of this software (i) may not refer to the modified
*  version by the same designation, or by any confusingly similar designation, and (ii) must refer to
*  the underlying software originally provided by Alliance as "SolTrace". Except to comply with the 
*  foregoing, the term "SolTrace", or any confusingly similar designation may not be used to refer to 
*  any modified version of this software or any modified version of the underlying software originally 
*  provided by Alliance without the prior written consent of Alliance.
*
*  5. The name of the copyright holder, contributors, the United States Government, the United States
*  Department of Energy, or any of their employees may not be used to endorse or promote products
*  derived from this software without specific prior written permission.
*
*  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
*  IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
*  FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER,
*  CONTRIBUTORS, UNITED STATES GOVERNMENT OR UNITED STATES DEPARTMENT OF ENERGY, NOR ANY OF THEIR
*  EMPLOYEES, BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
*  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
*  DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER
*  IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF
*  THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*******************************************************************************************************/


#ifndef __soltraceapi_h
#define __soltraceapi_h

#ifdef _STCOREDLL_
	#ifdef STCORE_API_EXPORTS
	#define STCORE_API __declspec(dllexport)
	#else
	#define STCORE_API __declspec(dllimport)
	#endif
#else
	#define STCORE_API
#endif

#ifndef STCORE_API
#define STCORE_API
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef unsigned long     st_uint_t;     // unsigned integer type, at least 32 bits, could be 64 bit in the future
typedef void*             st_context_t;  // opaque reference type, 32 or 64 bit, depending on system/compiler

/* functions to create system contexts */
STCORE_API st_context_t st_create_context();
STCORE_API int st_free_context(st_context_t pcxt);

/* functions for debugging systems */
STCORE_API int st_num_messages(st_context_t pcxt);
STCORE_API const char *st_message(st_context_t pcxt, st_uint_t idx);
STCORE_API void st_clear_messages(st_context_t pcxt);

/*
STCORE_API int st_dump( st_context_t pcxt, const char *file );
STCORE_API int st_load_file( st_context_t pcxt, const char *file );
STCORE_API int st_write_output(st_context_t pcxt, const char *file);
STCORE_API void st_reset(st_context_t pcxt);
*/

/* functions to add/remove optical properties */
STCORE_API int st_num_optics(st_context_t pcxt);
STCORE_API int st_add_optic(st_context_t pcxt, const char *name);
STCORE_API int st_delete_optic(st_context_t pcxt, st_uint_t idx);
STCORE_API int st_clear_optics(st_context_t pcxt);
STCORE_API int st_optic(st_context_t pcxt, st_uint_t idx, int fb, /* 1=front,2=back */
				char dist, int optnum, int apgr, int order,
				double rreal, double rimag,
				double ref, double tra,
				double gratingab12[3],
				double rmsslope, double rmsspec,
				int userefltable, int refl_npoints,
				double *refl_angles, double *refls,
				int usetranstable, int trans_npoints,
				double* trans_angles, double* transs);

/* functions to add/remove stages */
STCORE_API int st_num_stages(st_context_t pcxt);
STCORE_API int st_add_stage(st_context_t pcxt);
STCORE_API int st_add_stages(st_context_t pcxt, st_uint_t num);
STCORE_API int st_delete_stage(st_context_t pcxt, st_uint_t idx);
STCORE_API int st_clear_stages(st_context_t pcxt);
STCORE_API int st_stage_xyz(st_context_t pcxt, st_uint_t idx, double x, double y, double z);
STCORE_API int st_stage_aim(st_context_t pcxt, st_uint_t idx, double ax, double ay, double az);
STCORE_API int st_stage_zrot(st_context_t pcxt, st_uint_t idx, double zrot);
STCORE_API int st_stage_flags(st_context_t pcxt, st_uint_t idx, int virt, int multihit, int tracethrough);

/* functions to add/remove elements */
STCORE_API int st_num_elements(st_context_t pcxt, st_uint_t stage);
STCORE_API int st_add_element(st_context_t pcxt, st_uint_t stage);
STCORE_API int st_add_elements(st_context_t pcxt, st_uint_t stage, st_uint_t num);
STCORE_API int st_delete_element(st_context_t pcxt, st_uint_t stage, st_uint_t idx);
STCORE_API int st_clear_elements(st_context_t pcxt, st_uint_t stage);

/* functions to set element properties */
STCORE_API int st_element_enabled(st_context_t pcxt, st_uint_t stage, st_uint_t idx, int enabled);
STCORE_API int st_element_xyz(st_context_t pcxt, st_uint_t stage, st_uint_t idx, double x, double y, double z);
STCORE_API int st_element_aim(st_context_t pcxt, st_uint_t stage, st_uint_t idx, double ax, double ay, double az);
STCORE_API int st_element_zrot(st_context_t pcxt, st_uint_t stage, st_uint_t idx, double zrot);
STCORE_API int st_element_aperture(st_context_t pcxt, st_uint_t stage, st_uint_t idx, char ap);
STCORE_API int st_element_aperture_params(st_context_t pcxt, st_uint_t stage, st_uint_t idx, double params[8]);
STCORE_API int st_element_surface(st_context_t pcxt, st_uint_t stage, st_uint_t idx, char surf);
STCORE_API int st_element_surface_params(st_context_t pcxt, st_uint_t stage, st_uint_t idx, double params[8]);
STCORE_API int st_element_surface_file(st_context_t pcxt, st_uint_t stage, st_uint_t idx, const char *file);
STCORE_API int st_element_interaction(st_context_t pcxt, st_uint_t stage, st_uint_t idx, int type); /* 1=refract, 2=reflect */
STCORE_API int st_element_optic(st_context_t pcxt, st_uint_t stage, st_uint_t idx, const char *name);

/* functions to configure sun geometry and shape */
STCORE_API int st_sun(st_context_t pcxt, int point_source, char shape, double sigma_halfwidth);
STCORE_API int st_sun_xyz(st_context_t pcxt,  double x, double y, double z );
STCORE_API int st_sun_userdata(st_context_t pcxt,  st_uint_t npoints, double angle[], double intensity[]);

/* function to retrieve intersection data */
STCORE_API int st_num_intersections(st_context_t pcxt);
STCORE_API int st_locations(st_context_t pcxt, double *loc_x, double *loc_y, double *loc_z);
STCORE_API int st_cosines(st_context_t pcxt, double *cos_x, double *cos_y, double *cos_z);
STCORE_API int st_elementmap(st_context_t pcxt, int *element_map);
STCORE_API int st_stagemap(st_context_t pcxt, int *stage_map);
STCORE_API int st_raynumbers(st_context_t pcxt, int *ray_numbers);
STCORE_API int st_sun_stats(st_context_t pcxt, double *xmin, double *xmax, double *ymin, double *ymax, int *nsunrays );
	
/* functions to control simulation */
STCORE_API int st_sim_params(st_context_t pcxt, int raycount, int maxcount);
STCORE_API int st_sim_errors(st_context_t pcxt, int include_sun_shape, int include_optics);
STCORE_API int st_sim_run( st_context_t pcxt, unsigned int seed, bool AsPowerTower,
						  int (*callback)(st_uint_t ntracedtotal, st_uint_t ntraced, st_uint_t ntotrace, st_uint_t curstage, st_uint_t nstages, void *data), void *data);

/*
STCORE_API int st_sim_run_data( st_context_t pcxt, unsigned int seed, std::vector<std::vector< double > > *data_s1, std::vector<std::vector< double > > *data_s2, bool save_stage_data,
						  int (*callback)(st_uint_t ntracedtotal, st_uint_t ntraced, st_uint_t ntotrace, st_uint_t curstage, st_uint_t nstages, void *data), void *data);
*/


/* utility transform/math functions */
STCORE_API void st_calc_euler_angles( double origin[3], double aimpoint[3], double zrot, double euler[3] );
STCORE_API void st_transform_to_local( double posref[3], double cosref[3], double origin[3], double rreftoloc[3][3], double posloc[3], double cosloc[3]);
STCORE_API void st_transform_to_reference( double posloc[3], double cosloc[3], double origin[3], double rloctoref[3][3], double posref[3], double cosref[3]);
STCORE_API void st_matrix_vector_mult( double m[3][3], double v[3], double mxv[3] );
STCORE_API void st_calc_transform_matrices( double euler[3], double rreftoloc[3][3], double rloctoref[3][3] );
STCORE_API void st_matrix_transpose( double input[3][3], double output[3][3] );


#ifdef __cplusplus
} /* extern "C" */
#endif

#endif


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


#include "types.h"
#include "procs.h"
#include "stapi.h"
#include "mtrand.h"

#define SYSTEM(p,r) TSystem *sys = reinterpret_cast<TSystem*>(p); if(!sys) return r;
#define SYSTEM_NR(p) TSystem *sys = reinterpret_cast<TSystem*>(p); if(!sys) return;


STCORE_API st_context_t st_create_context()
{
	return reinterpret_cast<st_context_t>(new TSystem);
}
STCORE_API int st_free_context(st_context_t pcxt)
{
	SYSTEM(pcxt,-1);
	delete sys;
	return 1;
}

/* functions to get messages out of the core */
STCORE_API int st_num_messages(st_context_t pcxt)
{
	SYSTEM(pcxt,-1);
	return sys->messages.size();
}

STCORE_API void st_clear_messages(st_context_t pcxt)
{
	SYSTEM_NR(pcxt);
	sys->messages.clear();
}

STCORE_API const char *st_message(st_context_t pcxt, st_uint_t idx)
{
	SYSTEM(pcxt,NULL);
	if (idx >= 0 && idx < sys->messages.size())
		return sys->messages[idx].c_str();
	else
		return NULL;
}

/*
STCORE_API int st_dump( st_context_t pcxt, const char *file )
{
	SYSTEM(pcxt,-1);
	return DumpSystem( file, sys )?1:-1;
}

STCORE_API int st_write_output(st_context_t pcxt, const char *file)
{
	SYSTEM(pcxt,-1);
	return WriteOutput(file, sys, true)?1:-1;
}

STCORE_API int st_load_file( st_context_t pcxt, const char *file )
{
	SYSTEM(pcxt,-1);
	sys->ClearAll();
	if (ReadInputFile(file, sys))
		return 1;
	else
		return -1;
}

STCORE_API void st_reset(st_context_t pcxt)
{
	SYSTEM_NR(pcxt);
	sys->ClearAll();
}
*/

/* functions to add/remove optics */

STCORE_API int st_num_optics(st_context_t pcxt)
{
	SYSTEM(pcxt,-1);
	return sys->OpticsList.size();
}

STCORE_API int st_add_optic(st_context_t pcxt, const char *name)
{
	SYSTEM(pcxt,-1);
	sys->OpticsList.push_back( new TOpticalPropertySet );
	sys->OpticsList[ sys->OpticsList.size()-1 ]->Name = std::string(name);
	return sys->OpticsList.size()-1;
}

STCORE_API int st_delete_optic(st_context_t pcxt, st_uint_t idx)
{
	SYSTEM(pcxt,-1);
	if (idx >= 0 && idx < sys->OpticsList.size())
	{
		delete sys->OpticsList[idx];
		sys->OpticsList.erase( sys->OpticsList.begin() + idx );
		return 1;
	}
	else return -1;
}

STCORE_API int st_clear_optics(st_context_t pcxt)
{
	SYSTEM(pcxt,-1);
	for (st_uint_t i=0;i<sys->OpticsList.size();i++)
		delete sys->OpticsList[i];
	sys->OpticsList.clear();
	return 1;
}

STCORE_API int st_optic(st_context_t pcxt, st_uint_t idx, int fb, /* 1=front,2=back */
				char dist, int optnum, int apgr, int order,
				double rreal, double rimag,
				double ref, double tra,
				double gratingab12[3],
				double rmsslope, double rmsspec,
				int userefltable, int refl_npoints,
				double* refl_angles, double* refls,
				int usetranstable, int trans_npoints,
				double* trans_angles, double* transs)
{
	SYSTEM(pcxt,-1);

	TOpticalPropertySet *set = NULL;
	if (idx >= 0 && idx < sys->OpticsList.size())
		set = sys->OpticsList[idx];

	if (!set) return -1;

	TOpticalProperties *topt = (fb==2) ? &set->Back : &set->Front;

	topt->DistributionType = dist;
	topt->OpticSurfNumber = optnum;
	topt->ApertureStopOrGratingType = apgr;
	topt->DiffractionOrder = order;
	topt->RefractiveIndex[0] = rreal;
	topt->RefractiveIndex[1] = rimag;
	topt->Reflectivity = ref;
	topt->Transmissivity = tra;

	for (st_uint_t i=0;i<4;i++)
		topt->AB12[i] = gratingab12[i];

	topt->RMSSlopeError = rmsslope;
	topt->RMSSpecError = rmsspec;

	if (userefltable && refl_npoints > 0 && refl_angles != 0 && refls != 0)
	{
		topt->UseReflectivityTable = true;
		topt->ReflectivityTable.resize( refl_npoints );
		for (int i=0;i<refl_npoints;i++)
		{
			topt->ReflectivityTable[i].angle = refl_angles[i];
			topt->ReflectivityTable[i].refl = refls[i];
		}
	}
	if (usetranstable && trans_npoints > 0 && trans_angles != 0 && transs != 0)
	{
		topt->UseTransmissivityTable = true;
		topt->TransmissivityTable.resize(trans_npoints);
		for (int i = 0; i < trans_npoints; i++)
		{
			topt->TransmissivityTable[i].angle = trans_angles[i];
			topt->TransmissivityTable[i].trans = transs[i];
		}
	}

	return 1;
}



#define STAGE(i) ((i>=0&&i<sys->StageList.size())?sys->StageList[i]:NULL)
#define GETSTAGE(i) TStage *s = STAGE(i); if (!s) return -1;

/* functions to add/remove stages */
STCORE_API int st_num_stages(st_context_t pcxt)
{
	SYSTEM(pcxt,-1);
	return sys->StageList.size();
}

STCORE_API int st_add_stage(st_context_t pcxt)
{
	SYSTEM(pcxt,-1);
	sys->StageList.push_back( new TStage );
	return sys->StageList.size()-1;
}

STCORE_API int st_add_stages(st_context_t pcxt, st_uint_t num)
{
	SYSTEM(pcxt,-1);
	if (num < 0) return -1;

	for (st_uint_t i=0;i<num;i++)
		sys->StageList.push_back( new TStage );

	return num;
}

STCORE_API int st_delete_stage(st_context_t pcxt, st_uint_t idx)
{
	SYSTEM(pcxt,-1);
	if (idx >= 0 && idx < sys->StageList.size())
	{
		delete sys->StageList[idx];
		sys->StageList.erase( sys->StageList.begin() + idx );
		return 1;
	}
	else
		return -1;
}

STCORE_API int st_clear_stages(st_context_t pcxt)
{
	SYSTEM(pcxt,-1);
	for (st_uint_t i=0;i<sys->StageList.size();i++)
		delete sys->StageList[i];
	sys->StageList.clear();
	return 1;
}

STCORE_API int st_stage_xyz(st_context_t pcxt, st_uint_t idx, double x, double y, double z)
{
	SYSTEM(pcxt,-1);
	GETSTAGE(idx);
	s->Origin[0] = x;
	s->Origin[1] = y;
	s->Origin[2] = z;
	return 1;
}

STCORE_API int st_stage_aim(st_context_t pcxt, st_uint_t idx, double ax, double ay, double az)
{
	SYSTEM(pcxt,-1);
	GETSTAGE(idx);
	s->AimPoint[0] = ax;
	s->AimPoint[1] = ay;
	s->AimPoint[2] = az;
	return 1;
}

STCORE_API int st_stage_zrot(st_context_t pcxt, st_uint_t idx, double zrot)
{
	SYSTEM(pcxt,-1);
	GETSTAGE(idx);
	s->ZRot = zrot;
	return 1;
}

STCORE_API int st_stage_flags(st_context_t pcxt, st_uint_t idx, int virt, int multihit, int tracethrough)
{
	SYSTEM(pcxt,-1);
	GETSTAGE(idx);
	s->Virtual = virt?true:false;
	s->MultiHitsPerRay = multihit?true:false;
	s->TraceThrough = tracethrough?true:false;
	return 1;
}

//==============================================================

/* functions to add/remove elements */
STCORE_API int st_num_elements(st_context_t pcxt, st_uint_t stage)
{
	SYSTEM(pcxt,-1);
	GETSTAGE(stage);
	return s->ElementList.size();
}

STCORE_API int st_add_element(st_context_t pcxt, st_uint_t stage)
{
	SYSTEM(pcxt,-1);
	GETSTAGE(stage);
	s->ElementList.push_back( new TElement );
	return s->ElementList.size()-1;
}

STCORE_API int st_add_elements(st_context_t pcxt, st_uint_t stage, st_uint_t num)
{
	SYSTEM(pcxt,-1);
	if (num < 1)
		return -1;

	GETSTAGE(stage);

	for (st_uint_t i=0;i<num;i++)
		s->ElementList.push_back( new TElement );

	return num;
}

STCORE_API int st_delete_element(st_context_t pcxt, st_uint_t stage, st_uint_t idx)
{
	SYSTEM(pcxt,-1);
	GETSTAGE(stage);
	if (idx >= 0 && idx < s->ElementList.size())
	{
		delete s->ElementList[idx];
		s->ElementList.erase( s->ElementList.begin() + idx );
		return 1;
	}
	else
		return -1;
}

STCORE_API int st_clear_elements(st_context_t pcxt, st_uint_t stage)
{
	SYSTEM(pcxt,-1);
	GETSTAGE(stage)
	for (st_uint_t i=0;i<s->ElementList.size();i++)
		delete s->ElementList[i];
	s->ElementList.clear();
	return 1;
}

#define ELEMENT(i) ((i>=0&&i<s->ElementList.size())?s->ElementList[i]:NULL)
#define GETELEMENT(i) TElement *e = ELEMENT(i); if (!e) return -1;


/* functions to set element properties */
STCORE_API int st_element_enabled(st_context_t pcxt, st_uint_t stage, st_uint_t idx, int enabled)
{
	SYSTEM(pcxt,-1);
	GETSTAGE(stage);
	GETELEMENT(idx);
	e->Enabled = (enabled?true:false);
	return 1;
}

STCORE_API int st_element_xyz(st_context_t pcxt, st_uint_t stage, st_uint_t idx, double x, double y, double z)
{
	SYSTEM(pcxt,-1);
	GETSTAGE(stage);
	GETELEMENT(idx);
	e->Origin[0] = x;
	e->Origin[1] = y;
	e->Origin[2] = z;
	return 1;
}

STCORE_API int st_element_aim(st_context_t pcxt, st_uint_t stage, st_uint_t idx, double ax, double ay, double az)
{
	SYSTEM(pcxt,-1);
	GETSTAGE(stage);
	GETELEMENT(idx);
	e->AimPoint[0] = ax;
	e->AimPoint[1] = ay;
	e->AimPoint[2] = az;
	return 1;
}

STCORE_API int st_element_zrot(st_context_t pcxt, st_uint_t stage, st_uint_t idx, double zrot)
{
	SYSTEM(pcxt,-1);
	GETSTAGE(stage);
	GETELEMENT(idx);
	e->ZRot = zrot;
	return 1;
}

STCORE_API int st_element_aperture(st_context_t pcxt, st_uint_t stage, st_uint_t idx, char ap)
{
	SYSTEM(pcxt,-1);
	GETSTAGE(stage);
	GETELEMENT(idx);
	e->ShapeIndex = ap;
	return 1;
}

STCORE_API int st_element_aperture_params(st_context_t pcxt, st_uint_t stage, st_uint_t idx, double params[8])
{
	SYSTEM(pcxt,-1);
	GETSTAGE(stage);
	GETELEMENT(idx);
	e->ParameterA = params[0];
	e->ParameterB = params[1];
	e->ParameterC = params[2];
	e->ParameterD = params[3];
	e->ParameterE = params[4];
	e->ParameterF = params[5];
	e->ParameterG = params[6];
	e->ParameterH = params[7];
	return 1;
}

STCORE_API int st_element_surface(st_context_t pcxt, st_uint_t stage, st_uint_t idx, char surf)
{
	SYSTEM(pcxt,-1);
	GETSTAGE(stage);
	GETELEMENT(idx);
	e->SurfaceIndex = surf;
	return 1;
}

STCORE_API int st_element_surface_params(st_context_t pcxt, st_uint_t stage, st_uint_t idx, double params[8])
{
	SYSTEM(pcxt,-1);
	GETSTAGE(stage);
	GETELEMENT(idx);
	return TranslateSurfaceParams( sys, e, params ) ? 1 : -1;
}

STCORE_API int st_element_surface_file(st_context_t pcxt, st_uint_t stage, st_uint_t idx, const char *file)
{
	SYSTEM(pcxt,-1);
	GETSTAGE(stage);
	GETELEMENT(idx);
	return ReadSurfaceFile( file, e, sys ) ? 1 : -1;
}

STCORE_API int st_element_interaction(st_context_t pcxt, st_uint_t stage, st_uint_t idx, int type)
{
	SYSTEM(pcxt,-1);
	GETSTAGE(stage);
	GETELEMENT(idx);
	e->InteractionType = type;
	return 1;
}

STCORE_API int st_element_optic(st_context_t pcxt, st_uint_t stage, st_uint_t idx, const char *name)
{
	SYSTEM(pcxt,-1);
	GETSTAGE(stage);
	GETELEMENT(idx);

	e->OpticName = std::string(name);
	return 1;
}


/* functions to configure sun geometry and shape */
STCORE_API int st_sun(st_context_t pcxt, int point_source, char shape, double sigma_halfwidth )
{
	SYSTEM(pcxt,-1);
	sys->Sun.PointSource = point_source?true:false;
	sys->Sun.ShapeIndex = shape;
	sys->Sun.Sigma = sigma_halfwidth;
	return 1;
}

STCORE_API int st_sun_xyz( st_context_t pcxt, double x, double y, double z )
{
	SYSTEM(pcxt,-1);
	sys->Sun.Origin[0] = x;
	sys->Sun.Origin[1] = y;
	sys->Sun.Origin[2] = z;
	return 1;
}

STCORE_API int st_sun_userdata( st_context_t pcxt, st_uint_t npoints, double angle[], double intensity[])
{
	SYSTEM(pcxt,-1);
	if (npoints > 2)
	{
		sys->Sun.MaxAngle = 0;
		sys->Sun.MaxIntensity = 0;

		sys->Sun.SunShapeAngle.resize(2*npoints-1);
		sys->Sun.SunShapeIntensity.resize(2*npoints-1);

		for (st_uint_t i=0;i<npoints;i++)
		{
			sys->Sun.SunShapeAngle[npoints+i-1] = angle[i];
			sys->Sun.SunShapeIntensity[npoints+i-1] = intensity[i];

			if (angle[i] > sys->Sun.MaxAngle) sys->Sun.MaxAngle = angle[i];
			if (intensity[i] > sys->Sun.MaxIntensity ) sys->Sun.MaxIntensity = intensity[i];
		}

		// fill negative angle side of array
		for (st_uint_t i=0;i<npoints-1;i++)
		{
			sys->Sun.SunShapeAngle[i] = -angle[npoints-i-1];
			sys->Sun.SunShapeIntensity[i] = intensity[npoints-i-1];
		}
		
		return npoints;
	}
	else
		return -1;
}


/* functions to retrieve intersection data */
STCORE_API int st_num_intersections(st_context_t pcxt)
{
	SYSTEM(pcxt,-1);
	return sys->AllRayData.Count();
}

STCORE_API int st_locations(st_context_t pcxt, double *loc_x, double *loc_y, double *loc_z)
{
	SYSTEM(pcxt,-1);

	double pos[3];
	for (st_uint_t i=0;i<sys->AllRayData.Count();i++)
	{
		if( sys->AllRayData.Query(i, pos, 0, 0, 0, 0) )
		{
			loc_x[i] = pos[0];
			loc_y[i] = pos[1];
			loc_z[i] = pos[2];
		}
		else
			return -2;
	}
	return 1;
}

STCORE_API int st_cosines(st_context_t pcxt, double *cos_x, double *cos_y, double *cos_z)
{
	SYSTEM(pcxt,-1);
	double cos[3];
	for (st_uint_t i=0;i<sys->AllRayData.Count();i++)
	{
		if (sys->AllRayData.Query(i, 0, cos, 0, 0, 0))
		{
			cos_x[i] = cos[0];
			cos_y[i] = cos[1];
			cos_z[i] = cos[2];
		}
		else
			return -2;
	}
	return 1;
}

STCORE_API int st_elementmap(st_context_t pcxt, int *element_map)
{
	SYSTEM(pcxt,-1);
	int elm;
	for (st_uint_t i=0;i<sys->AllRayData.Count();i++)
	{
		if (sys->AllRayData.Query(i, 0, 0, &elm, 0, 0))
			element_map[i] = elm;
		else
			return -2;
	}
	return 1;
}

STCORE_API int st_stagemap(st_context_t pcxt, int *stage_map)
{
	SYSTEM(pcxt,-1);
	int stage;
	for (st_uint_t i=0;i<sys->AllRayData.Count();i++)
	{
		if (sys->AllRayData.Query(i, 0, 0, 0, &stage, 0))
			stage_map[i] = stage;
		else
			return -2;
	}
	return 1;
}

STCORE_API int st_raynumbers(st_context_t pcxt, int *ray_numbers)
{
	SYSTEM(pcxt,-1);
	unsigned int ray;
	for (st_uint_t i=0;i<sys->AllRayData.Count();i++)
	{
		if (sys->AllRayData.Query(i, 0, 0, 0, 0, &ray))
			ray_numbers[i] = (int)ray;
		else
			return -2;
	}
	return 1;
}


STCORE_API int st_sun_stats( st_context_t pcxt, double *xmin, double *xmax, double *ymin, double *ymax, int *nsunrays )
{
	SYSTEM(pcxt,-1);
	if (xmin) *xmin = sys->Sun.MinXSun;
	if (xmax) *xmax = sys->Sun.MaxXSun;
	if (ymin) *ymin = sys->Sun.MinYSun;
	if (ymax) *ymax = sys->Sun.MaxYSun;
	if (nsunrays) *nsunrays = sys->SunRayCount;
	return 1;
}


/* functions to control simulation */
STCORE_API int st_sim_params(st_context_t pcxt, int raycount, int maxcount)
{
	SYSTEM(pcxt,-1);
	sys->sim_raycount = raycount;
	sys->sim_raymax = maxcount;
	return 1;
}

STCORE_API int st_sim_errors(st_context_t pcxt, int include_sun_shape, int include_optics)
{
	SYSTEM(pcxt,-1);
	sys->sim_errors_sunshape = include_sun_shape?true:false;
	sys->sim_errors_optical = include_optics?true:false;
	return 1;
}

STCORE_API int st_sim_run_data( st_context_t pcxt, unsigned int seed, 
                            bool AsPowerTower,
                            std::vector<std::vector< double > > *data_s1, 
                            std::vector<std::vector< double > > *data_s2, 
                            bool save_st_data,
						    int (*callback)(st_uint_t ntracedtotal, st_uint_t ntraced, st_uint_t ntotrace, st_uint_t curstage, st_uint_t nstages, void *data), void *cbdata
                            )
{
    /* 
    Run a simulation and save the stage 1-2 ray data for future use. Return this data in the data_s1, data_s2 structures.
    data_s1 will be formatted each vector entry:

    PosRayGlobal[0],PosRayGlobal[1],PosRayGlobal[2],CosRayGlobal[0],CosRayGlobal[1],CosRayGlobal[2],ElementNum,RayNum

    data_s2 will be formatted each vector entry:
    PosRayGlobal[0],PosRayGlobal[1],PosRayGlobal[2],CosRayGlobal[0],CosRayGlobal[1],CosRayGlobal[2],RayNum
    
    */
    SYSTEM(pcxt,-1);

	sys->AllRayData.Clear();

	if ( !InitGeometries(sys) )
		return -1;

    int rayct = sys->sim_raycount;
    if(data_s2 != 0)
        if(data_s2->size() > 0)
            rayct = data_s2->size();

	if ( !Trace(sys, seed,
		rayct, sys->sim_raymax,
		sys->sim_errors_sunshape, sys->sim_errors_optical, AsPowerTower,
		callback, cbdata, data_s1, data_s2, save_st_data) )
		return -1;


	try
	{
		for (st_uint_t i=0;i<sys->StageList.size();i++)
			sys->AllRayData.Merge( sys->StageList[i]->RayData );

		return sys->AllRayData.Count();
	}
	catch( const std::exception &e )
	{
		sys->errlog("caught exception: %s", e.what());
		return -1;
	}
}

STCORE_API int st_sim_run( st_context_t pcxt, unsigned int seed, bool AsPowerTower,
						  int (*callback)(st_uint_t ntracedtotal, st_uint_t ntraced, st_uint_t ntotrace, st_uint_t curstage, st_uint_t nstages, void *data), void *cbdata)
{
    return st_sim_run_data( pcxt, seed, AsPowerTower, 0, 0, false, callback, cbdata);
}


STCORE_API void st_calc_euler_angles( double origin[3], double aimpoint[3], double zrot, double euler[3] )
{
	double dx, dy, dz;
	euler[0] = euler[1] = euler[2] = 0;
	dx = aimpoint[0]-origin[0];
	dy = aimpoint[1]-origin[1];
	dz = aimpoint[2]-origin[2];
	double dtot = sqrt(dx*dx + dy*dy + dz*dz);
	if (dtot == 0.0)
		return;	
	dx = dx/dtot;
	dy = dy/dtot;
	dz = dz/dtot;
	euler[0] = atan2(dx,dz);
	euler[1] = asin(dy);
	euler[2] = zrot*(ACOSM1O180);
}

STCORE_API void st_transform_to_local( double posref[3], double cosref[3], double origin[3], double rreftoloc[3][3], double posloc[3], double cosloc[3])
{
	TransformToLocal(posref, cosref, origin, rreftoloc, posloc, cosloc);
}

STCORE_API void st_transform_to_reference( double posloc[3], double cosloc[3], double origin[3], double rloctoref[3][3], double posref[3], double cosref[3])
{
	TransformToReference(posloc, cosloc, origin, rloctoref, posref, cosref);
}

STCORE_API void st_matrix_vector_mult( double m[3][3], double v[3], double mxv[3] )
{
	MatrixVectorMult(m,v,mxv);
}

STCORE_API void st_calc_transform_matrices( double euler[3], double rreftoloc[3][3], double rloctoref[3][3] )
{
	CalculateTransformMatrices(euler, rreftoloc, rloctoref);
}

STCORE_API void st_matrix_transpose( double input[3][3], double output[3][3] )
{
	MatrixTranspose(input, 3, output);
}


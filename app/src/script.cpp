
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


#include <lk/env.h>

#include "soltrace.h"
#include "project.h"
#include "trace.h"
#include "sunshape.h"
#include "geometry.h"
#include "optics.h"
//#include "fluxmapform.h"

#include "script.h"

static void get_3( double v[3], lk::vardata_t &arg )
{
	for (int i=0;i<3;i++)
		v[i] = arg.index(i)->as_number();
}

static void set_3( lk::vardata_t &arg, double v[3] )
{
	arg.empty_vector();
	arg.vec_append( v[0] );
	arg.vec_append( v[1] );
	arg.vec_append( v[2] );
}

static void get_3x3( double m[3][3], lk::vardata_t &arg )
{
	for (int i=0;i<3;i++)
	{
		lk::vardata_t *v = arg.index(i);
		for (int j=0;j<3;j++)
			m[i][j] = v!=0 ? v->index(j)->as_number() : 0.0;
	}
}

static void set_3x3( lk::vardata_t &arg, double m[3][3] )
{
	arg.empty_vector();
	arg.resize(3);
	for (int i=0;i<3;i++)
	{
		arg.index(i)->empty_vector();
		arg.index(i)->resize(3);
		for (int j=0;j<3;j++)
			arg.index(i)->index(j)->assign( m[i][j] );
	}
}



static void _dot( lk::invoke_t &cxt )
{
	LK_DOC("dot", "Calculates the dot product between two 3 dimensional vectors.", "(array[3]:a, array[3]:b):real");
	double a[3],b[3], dot=0;
	get_3(a, cxt.arg(0));
	get_3(b, cxt.arg(1));
	for (int i=0;i<3;i++)
		dot += a[i]*b[i];
	cxt.result().assign(dot);
}

static void _euler( lk::invoke_t &cxt )
{
	LK_DOC("euler", "Calculates the euler vector given an origin, aim point, and z rotation.", "(array[3]:origin, array[3]:aimpoint, real:zrot):array[3]");

	double origin[3], aimpoint[3], zrot, eu[3];
	get_3( origin, cxt.arg(0) ); //	args[0] //origin
	get_3( aimpoint, cxt.arg(1) ); //	args[1] //aimpoint
	zrot = cxt.arg(2).as_number(); //	args[2] //zrot
	::st_calc_euler_angles( origin, aimpoint, zrot, eu );
	set_3(cxt.result(), eu);
}

static void _reftoloc( lk::invoke_t &cxt )
{
	LK_DOC("reftoloc", "Calculates the 3x3 transform matrix from euler angles for going from reference to local coordinates.", "(array[3]:euler):array[3][3]");

	double euler[3], rreftoloc[3][3], rloctoref[3][3];
	get_3( euler, cxt.arg(0) );
	::st_calc_transform_matrices( euler, rreftoloc, rloctoref );
	set_3x3( cxt.result(), rreftoloc );
}

static void _loctoref( lk::invoke_t &cxt )
{
	LK_DOC("loctoref", "Calculates the 3x3 transform matrix from euler angles for going from local to reference coordinates.", "(array[3]:euler):array[3][3]");

	double euler[3], rreftoloc[3][3], rloctoref[3][3];
	get_3( euler, cxt.arg(0) );
	::st_calc_transform_matrices( euler, rreftoloc, rloctoref );
	set_3x3( cxt.result(), rloctoref );
}

static void _toloc( lk::invoke_t &cxt )
{
	LK_DOC("toloc", "Transforms a reference position/cosine point description to local coordinates given the transform matrix.", "(array[3]:posref, array[3]:cosref, array[3]:origin, array[3][3]:reftoloc, array[3]:posloc, array[3]:cosloc):void");

	double posref[3], cosref[3], origin[3], rreftoloc[3][3], posloc[3], cosloc[3];

	get_3(posref, cxt.arg(0));
	get_3(cosref, cxt.arg(1));
	get_3(origin, cxt.arg(2));
	get_3x3(rreftoloc, cxt.arg(3));

	::st_transform_to_local( posref, cosref, origin, rreftoloc, posloc, cosloc);

	set_3(cxt.arg(4), posloc);
	set_3(cxt.arg(5), cosloc);
}

static void _toref( lk::invoke_t &cxt )
{
	LK_DOC("toref", "Transforms a local position/cosine point description to reference coordinates given the transform matrix.", "(array[3]:posloc, array[3]:cosloc, array[3]:origin, array[3][3]:loctoref, array[3]:posref, array[3]:cosref):void");

	double posloc[3], cosloc[3], origin[3], rloctoref[3][3], posref[3], cosref[3];

	get_3(posloc, cxt.arg(0));
	get_3(cosloc, cxt.arg(1));
	get_3(origin, cxt.arg(2));
	get_3x3(rloctoref, cxt.arg(3));

	::st_transform_to_reference(posloc,cosloc,origin,rloctoref,posref,cosref);

	set_3( cxt.arg(4), posref );
	set_3( cxt.arg(5), cosref );
}

static void _transpose( lk::invoke_t &cxt )
{
	LK_DOC("transpose", "Transposes a 3x3 matrix.", "(array[3][3]):array[3][3]");

	double m1[3][3], m2[3][3];
	get_3x3( m1, cxt.arg(0) );
	::st_matrix_transpose( m1, m2 );
	set_3x3( cxt.result(), m2 );
}

static void _matvecmult( lk::invoke_t &cxt )
{
	LK_DOC("matvecmult", "Multiples a 3x3 matrix and a vector, returning the result.", "(array[3][3]:m, array[3]:v):array[3]");

	double m[3][3], v[3], mxv[3];
	get_3x3(m, cxt.arg(0));
	get_3(v, cxt.arg(1));
	::st_matrix_vector_mult(m,v,mxv);
	set_3( cxt.result(), mxv );
}


static void _workdir( lk::invoke_t &cxt )
{
	LK_DOC2("workdir", "Gets or sets the current SolTrace working directory.",
			"Returns the current working directory for SolTrace.", "(void):string",
			"Sets the current working directory for SolTrace.", "(string:path):void");

	if (cxt.arg_count() == 1)
		MainWindow::Instance().GetTrace()->SetWorkDir( cxt.arg(0).as_string() );
	else
		cxt.result().assign( MainWindow::Instance().GetTrace()->GetWorkDir() );
}

static void _file_name( lk::invoke_t &cxt )
{
	LK_DOC("file_name", "Returns the current SolTrace project file name.", "(void):string");
	cxt.result().assign( MainWindow::Instance().GetFileName() );
}

static void _save_project( lk::invoke_t &cxt )
{
	LK_DOC("save_project", "Saves the current SolTrace project as an .stinput file.", "({string:optional file name}):boolean");
	wxString file = MainWindow::Instance().GetFileName();
	if (cxt.arg_count() > 0) file = cxt.arg(0).as_string();
	if (!file.IsEmpty()) cxt.result().assign( MainWindow::Instance().SaveProject(file, true) ? 1.0 : 0.0 );
	else cxt.result().assign( 0.0 );
}

static void _open_project( lk::invoke_t &cxt )
{
	LK_DOC("open_project", "Opens an existing SolTrace .stinput project file, without asking to save or close any currently opened project.", "(string:filename):boolean");
	wxString file = cxt.arg(0).as_string();
	MainWindow::Instance().CloseProject( true );
	cxt.result().assign( MainWindow::Instance().LoadProject( file, true ) ? 1.0 : 0.0 );
}

static void _clear_project( lk::invoke_t &cxt )
{
	LK_DOC("clear_project", "Clears the current SolTrace project, without asking to save or close first.", "(void):void");
	MainWindow::Instance().CloseProject( true );
}

static void _traceopt( lk::invoke_t &cxt )
{
	LK_DOC2("traceopt", "Two modes of operation: Gets or sets ray trace parameters. For example: traceopt( {\"seed\"=152, \"cpus\"=3} ) sets the seed value to 152 and the number of CPUs to use to 3.",
			"Sets various ray trace parameters. The argument is a table with keys {rays,maxrays,cpus,seed,include_sunshape,optical_errors,point_focus}, whose values are the corresponding integers.", "(table:parameters):void",
			"Returns a table with the following five fields filled in with their integer values: {rays,maxrays,cpus,seed,include_sunshape,optical_errors,point_focus}", "(void):table" );

	TraceForm *tf = MainWindow::Instance().GetTrace();
	size_t nrays, nmax;
	int ncpu, seed;
	bool ss, oe, pf;
	tf->GetOptions( &nrays, &nmax, &ncpu, &seed, &ss, &oe, &pf );

	if (cxt.arg_count() == 0)
	{
		lk::vardata_t &r = cxt.result();
		r.empty_hash();
		r.hash_item("rays", nrays );
		r.hash_item("maxrays", nmax );
		r.hash_item("cpus", ncpu );
		r.hash_item("seed", seed );
		r.hash_item("include_sunshape", ss ? 1.0 : 0.0 );
		r.hash_item("optical_errors", oe ? 1.0 : 0.0 );
		r.hash_item("point_focus", pf ? 1.0 : 0.0 );
	}
	else if (cxt.arg_count() == 1)
	{
		lk::vardata_t *vval = 0;
		if ( (vval = cxt.arg(0).lookup("rays")) )
			nrays = vval->deref().as_unsigned();

		if ( (vval = cxt.arg(0).lookup("maxrays")) )
			nmax = vval->deref().as_unsigned();

		if ( (vval = cxt.arg(0).lookup("cpus")) )
			ncpu = vval->deref().as_integer();

		if ( (vval = cxt.arg(0).lookup("seed")) )
			seed = vval->deref().as_integer();

		if ( (vval = cxt.arg(0).lookup("include_sunshape")) )
			ss = vval->deref().as_integer() ? true : false;

		if ( (vval = cxt.arg(0).lookup("optical_errors")) )
			oe = vval->deref().as_integer() ? true : false;
		
		if ( (vval = cxt.arg(0).lookup("point_focus")) )
			pf = vval->deref().as_integer() ? true : false;

		tf->SetOptions( nrays, nmax, ncpu, seed, ss, oe, pf );
	}
	else
	{
		cxt.error("invalid number of arguments to traceopt.  must be 0 or 1");
	}
}

static void _trace( lk::invoke_t &cxt )
{
	LK_DOC("trace", "Starts a new ray trace operation.  Returns time elapsed in milliseconds, or negative if an error occured or was canceled.  Optionally, fills the 1st argument with error messages if passed in.", "( [string:errors] ):integer");
	wxArrayString errors;
	int ms = MainWindow::Instance().GetTrace()->StartTrace( true, true, &errors );
	if ( ms < 0 && cxt.arg_count() == 1 )
		cxt.arg(0).assign( wxJoin( errors, '\n' ) );

	cxt.result().assign( ms );
}

static void _nintersect( lk::invoke_t &cxt )
{
	LK_DOC2("nintersect", "Two modes of operation: returns the total number of intersections calculated, or the number of intersections with a particular element (stagenum, elementnum)",
			"Returns the number of ray intersections in the results.", "(void):integer",
			"Returns the number of ray intersections with a particular element.", "(integer:stagenum, integer:elementnum):integer");
	if (cxt.arg_count() == 0)
	{
		cxt.result().assign( MainWindow::Instance().GetProject().Results.Length );
	}
	else if (cxt.arg_count() == 2)
	{
		Element *e = MainWindow::Instance().GetProject().GetElement(
				cxt.arg(0).as_integer(), cxt.arg(1).as_integer() );
		cxt.result().assign( e ? e->RayHits : 0.0 );
	}
	else
		cxt.error("invalid number of arguments to nintersect. must be 0 or 2");
}

static void _raydata( lk::invoke_t &cxt )
{
	LK_DOC( "raydata", "Returns the ray data as a 9 item array in stage coordinates [X,Y,Z,CosX,CosY,CosZ,Element,Stage,RayNum] for the specified intersection number. Use nintersect to get the number of intersections.", "(integer:index):array");
	Project &prj = MainWindow::Instance().GetProject();
	size_t idx = cxt.arg(0).as_unsigned();
	if (idx <= prj.Results.Length)
	{
		lk::vardata_t &r = cxt.result();

		r.empty_vector();
		r.vec_append( prj.Results.Xi[idx] );
		r.vec_append( prj.Results.Yi[idx] );
		r.vec_append( prj.Results.Zi[idx] );
		r.vec_append( prj.Results.Xc[idx] );
		r.vec_append( prj.Results.Yc[idx] );
		r.vec_append( prj.Results.Zc[idx] );
		r.vec_append( prj.Results.ElementMap[idx] );
		r.vec_append( prj.Results.StageMap[idx] );
		r.vec_append( prj.Results.RayNumbers[idx] );
	}
	else
		cxt.result().nullify();
}

static void _sundata( lk::invoke_t &cxt )
{
	LK_DOC( "sundata", "Returns the number of generated sun rays and the extents as a table with fields {nrays,xmin,xmax,ymin,ymax}.", "(void):table");
	Project &prj = MainWindow::Instance().GetProject();
	lk::vardata_t &r = cxt.result();
	r.empty_hash();
	r.hash_item( "nrays", prj.Results.SunRayCount );
	r.hash_item( "xmin",  prj.Results.SunXMin );
	r.hash_item( "xmax",  prj.Results.SunXMax );
	r.hash_item( "ymin",  prj.Results.SunYMin );
	r.hash_item( "ymax",  prj.Results.SunYMax );
}

static void _sunopt( lk::invoke_t &cxt )
{
	LK_DOC2("sunopt", "Two modes of operation.  Gets or sets the sun shape parameters using a table with fields {ptsrc:boolean, shape:character, sigma:real, halfwidth:real, x:real, y:real, z:real, lat:real, day:real, hour:real, userdata:array[array[2]]}.",
		   "Sets various sun shape parameters. The argument is a table with keys {ptsrc, shape, sigma, halfwidth, x, y, z, useldh, lat, day, hour, userdata}.", "(table:parameters):void",
		   "Gets various sun shape parameters as a field indexed table.", "(void):table");

	SunShape &sun = MainWindow::Instance().GetProject().Sun;

	if (cxt.arg_count() == 1)
	{
		lk::vardata_t &h = cxt.arg(0);
		lk::vardata_t *v = 0;

		if ( (v=h.lookup("ptsrc")) )
			sun.PointSource = v->deref().as_boolean();

		if ( (v=h.lookup("shape")) )
		{
			sun.Shape = SunShape::PILLBOX;
			wxString s = v->deref().as_string();
			if (s.length() > 0)
			{
				if (tolower(s[0]) == 'g') sun.Shape = SunShape::GAUSSIAN;
				if (tolower(s[0]) == 'u' || tolower(s[0]) == 'd') sun.Shape = SunShape::USER_DEFINED;
			}
		}

		if ( (v=h.lookup("sigma")) )
			sun.Sigma = v->deref().as_number();

		if ( (v=h.lookup("halfwidth")) )
			sun.HalfWidth = v->deref().as_number();

		if ( (v=h.lookup("x")) )
			sun.X = v->deref().as_number();

		if ( (v=h.lookup("y")) )
			sun.Y = v->deref().as_number();

		if ( (v=h.lookup("z")) )
			sun.Z = v->deref().as_number();

		if ( (v=h.lookup("useldh")) )
			sun.UseLDHSpec = v->deref().as_boolean();

		if ( (v=h.lookup("lat")) )
			sun.Latitude = v->deref().as_number();

		if ( (v=h.lookup("day")) )
			sun.Day = v->deref().as_number();

		if ( (v=h.lookup("hour")) )
			sun.Hour = v->deref().as_number();

		if ( (v=h.lookup("userdata")) )
		{
			std::vector<lk::vardata_t> *arr = v->deref().vec();
			sun.UserShapeData.resize( arr->size() );
			for (size_t i=0;i<arr->size();i++)
			{
				std::vector<lk::vardata_t> *p = arr->at(i).deref().vec();
				if (p->size() == 2)
				{
					sun.UserShapeData[i].x = p->at(0).as_number();
					sun.UserShapeData[i].y = p->at(1).as_number();
				}
			}
		}

		MainWindow::Instance().GetSunShape()->UpdateFromData();
		MainWindow::Instance().SetModified();
	}
	else if (cxt.arg_count() == 0)
	{
		lk::vardata_t &r = cxt.result();
		r.empty_hash();
		r.hash_item( "ptsrc", sun.PointSource ? 1.0 : 0.0 );

		wxString shape = "g";
		if (sun.Shape == SunShape::PILLBOX) shape = "p";
		if (sun.Shape == SunShape::USER_DEFINED) shape = "u";
		r.hash_item( "shape", shape );

		r.hash_item( "sigma", sun.Sigma );
		r.hash_item( "halfwidth", sun.HalfWidth );
		r.hash_item( "x", sun.X );
		r.hash_item( "y", sun.Y );
		r.hash_item( "z", sun.Z );
		r.hash_item( "useldh", sun.UseLDHSpec ? 1.0 : 0.0 );
		r.hash_item( "lat", sun.Latitude );
		r.hash_item( "day", sun.Day );
		r.hash_item( "hour", sun.Hour );

		lk::vardata_t ud;
		ud.empty_vector();
		ud.resize( sun.UserShapeData.size() );
		for (int i=0;i<sun.UserShapeData.size();i++)
		{
			ud.index(i)->empty_vector();
			ud.index(i)->resize(2);
			ud.index(i)->index(0)->assign( sun.UserShapeData[i].x );
			ud.index(i)->index(1)->assign( sun.UserShapeData[i].y );
		}

		r.hash_item( "userdata", ud );
	}
	else
		cxt.error("invalid number of arguments for sunopt. must be 1 or 0.");
}

static void _addoptic( lk::invoke_t &cxt )
{
	LK_DOC("addoptic", "Adds a new optical property set with the given name.", "(string:name):void");
	wxString name = cxt.arg(0).as_string();
	if (name.empty()) name = "untitled optic";
	MainWindow::Instance().GetOptics()->AddOptic( name );
}

static void _clearoptics( lk::invoke_t &cxt )
{
	LK_DOC("clearoptics", "Deletes all of the optical property sets.", "(void):void");
	MainWindow::Instance().GetOptics()->ClearOptics();
}

static void _listoptics( lk::invoke_t &cxt )
{
	LK_DOC("listoptics", "Returns a list of all the optical property sets.", "(void):array");
	cxt.result().empty_vector();

	Project &prj = MainWindow::Instance().GetProject();
	for (size_t i=0;i<prj.OpticsList.size();i++)
		cxt.result().vec_append( prj.OpticsList[i]->Name );
}

static void _opticopt( lk::invoke_t &cxt )
{
	LK_DOC2("opticopt", "Two modes of operation: gets or sets optical property information, using a table with fields {dist=string, apstop=integer, surfnum=integer, difford=integer, refl=real, trans=real, errslope=real, errspec=real, refractr=real, refracti=real, grating=array[4], refltable=array[nx2]}",
			"Set various optical properties for the given optic name and side (1=front,2=back) with a table whose fields are {dist, apstop, surfnum, difford, refl, trans, errslope, errspec, refractr, refracti, grating, refltable}", "(string:name, integer:front or back, table:properties):void",
			"Get optical property information as table with named fields for the given optic name and side (1=front,2=back)", "(string:name, integer:front or back):table");

	wxString name = cxt.arg(0).as_string();
	int side = cxt.arg(1).as_integer();

	Project &prj = MainWindow::Instance().GetProject();
	Optical *opt = 0;

	for (size_t idx=0;idx<prj.OpticsList.size();idx++)
		if (prj.OpticsList[idx]->Name.CmpNoCase(name)==0)
			opt = prj.OpticsList[idx];

	if (!opt)
	{
		cxt.result().nullify();
		return;
	}

	SurfaceOptic &o = (side==2) ? opt->Back : opt->Front;

	//{ dist=string, apstop=integer, surfnum=integer, difford=integer,
	//  refl=real, trans=real, errslope=real, errspec=real, refractr=real, refracti=real,
	//  grating=array[4]}

	if (cxt.arg_count() == 3)
	{
		lk::vardata_t &h = cxt.arg(2);
		lk::vardata_t *v = 0;

		if ( (v=h.lookup("dist")) )
		{
			char d = 'g';
			wxString param = v->deref().as_string();
			if ( param.length() > 0 && tolower( param[0] ) == 'p' ) d = 'p';
			if ( param.length() > 0 && tolower( param[0] ) == 'f' ) d = 'f';
			o.ErrorDistribution = d;
		}

		if ( (v=h.lookup("apstop")) )
			o.ApertureStopOrGratingType = v->deref().as_integer();

		if ( (v=h.lookup("surfnum")) )
			o.OpticalSurfaceNumber = v->deref().as_integer();

		if ( (v=h.lookup("difford")) )
			o.DiffractionOrder = v->deref().as_integer();

		if ( (v=h.lookup("refl")) )
			o.Reflectivity = v->deref().as_number();

		if ( (v=h.lookup("trans")) )
			o.Transmissivity = v->deref().as_number();

		if ( (v=h.lookup("errslope")) )
			o.RMSSlope = v->deref().as_number();

		if ( (v=h.lookup("errspec")) )
			o.RMSSpecularity = v->deref().as_number();

		if ( (v=h.lookup("refractr")) )
			o.RefractionIndexReal = v->deref().as_number();

		if ( (v=h.lookup("refracti")) )
			o.RefractionIndexImag = v->deref().as_number();

		if ( (v=h.lookup("grating")) )
		{
			for (int i=0;i<4;i++)
				o.GratingCoeffs[i] = v->deref().index(i)->as_number();
		}

		if ( (v=h.lookup("refltable")) )
		{
			if (v->deref().type() == lk::vardata_t::NULLVAL)
				o.UseReflectivityTable = false;
			else
			{
				o.UseReflectivityTable = true;
				size_t npoints = v->deref().length();
				if (npoints > 0)
					o.ReflectivityTable.resize( npoints );

				for (size_t i=0;i<npoints;i++)
				{
					o.ReflectivityTable[i].x = 0.0;
					o.ReflectivityTable[i].y = 0.0;

					lk::vardata_t *item = v->deref().index(i);
					if ( item != 0 )
					{
						if (lk::vardata_t *angle = item->deref().index(0))
							o.ReflectivityTable[i].x = angle->as_number();

						if (lk::vardata_t *refl = item->deref().index(1))
							o.ReflectivityTable[i].y = refl->as_number();
					}
				}
			}
		}

		MainWindow::Instance().GetOptics()->UpdateOptForms();
		MainWindow::Instance().SetModified();
	}
	else if (cxt.arg_count() == 2)
	{
		lk::vardata_t &r = cxt.result();
		r.empty_hash();

		wxString buf;
		buf += o.ErrorDistribution;
		r.hash_item( "dist", buf);
		r.hash_item( "apstop", o.ApertureStopOrGratingType );
		r.hash_item( "surfnum", o.OpticalSurfaceNumber );
		r.hash_item( "difford", o.DiffractionOrder );
		r.hash_item( "refl", o.Reflectivity );
		r.hash_item( "trans", o.Transmissivity );
		r.hash_item( "errslope", o.RMSSlope );
		r.hash_item( "errspec", o.RMSSpecularity );
		r.hash_item( "refractr", o.RefractionIndexReal );
		r.hash_item( "refracti", o.RefractionIndexImag );
		lk::vardata_t grat;
		grat.empty_vector();
		grat.resize(4);
		for (int j=0;j<4;j++)
			grat.index(j)->assign( o.GratingCoeffs[j] );
		r.hash_item( "grating", grat );

		size_t npoints = o.ReflectivityTable.size();
		if (o.UseReflectivityTable && npoints > 0)
		{
			lk::vardata_t tab;
			tab.empty_vector();
			tab.resize( npoints );
			for (size_t j=0;j<npoints;j++)
			{
				tab.index(j)->empty_vector();
				tab.index(j)->resize(2);
				tab.index(j)->index(0)->assign( o.ReflectivityTable[j].x );
				tab.index(j)->index(1)->assign( o.ReflectivityTable[j].y );
			}

			r.hash_item( "refltable", tab );
		}
		else
			r.hash_item( "refltable", lk::vardata_t() ); // set to NULL value
	}
	else
		cxt.error("invalid number of arguments to opticopt. must be 3 or 2.");
}

// tracks the currently active stage for the elementXX functions
static wxString active_stage_name;

static void _addstage( lk::invoke_t &cxt )
{
	LK_DOC("addstage", "Adds a stage to the system with the given name.  Sets the new stage as the currently active stage.", "(string:name):void");
	wxString name = cxt.arg(0).as_string();
	if (name.empty()) name = "untitled stage";

	MainWindow::Instance().GetGeometry()->NewStage( name );
	active_stage_name = name;
}

static void _clearstages(lk::invoke_t &cxt)
{
	LK_DOC("clearstages", "Clears all stages from the project.", "(void):void");
	MainWindow::Instance().GetGeometry()->ClearStages();
	active_stage_name.Clear();
}

static void _liststages( lk::invoke_t &cxt )
{
	LK_DOC("liststages", "Returns a list of all the stage names.", "(void):array");
	cxt.result().empty_vector();

	Project &prj = MainWindow::Instance().GetProject();
	for (size_t i=0;i<prj.StageList.size();i++)
		cxt.result().vec_append( prj.StageList[i]->Name );
}

static void _activestage( lk::invoke_t &cxt )
{

	LK_DOC2("activestage", "Two modes of operation: gets or sets the currently active stage that the various 'elementXX functions operate on.",
			"Set the currently active stage.  Returns the active stage name also.", "(string:stage name):string",
			"Get the currently active stage name.", "(void):string");

	if (MainWindow::Instance().GetProject().StageList.size() == 0)
		active_stage_name.Clear();

	if ( cxt.arg_count() > 0)
	{
		wxString name = cxt.arg(0).as_string();
		Project &prj = MainWindow::Instance().GetProject();
		for (size_t i=0;i<prj.StageList.size();i++)
			if (prj.StageList[i]->Name.CmpNoCase( name ) == 0)
				active_stage_name = prj.StageList[i]->Name;
	}

	cxt.result().assign( active_stage_name );
}


static void _stageopt( lk::invoke_t &cxt )
{
	LK_DOC2("stageopt", "Two modes of operation. Gets or sets stage properties using a table with fields { virtual:boolean, multihit:boolean, tracethrough:boolean, x:real, y:real, z:real, ax:real, ay:real, az:real, zrot:real }",
			"Set various stage properties for the given stage name with a table whose fields can include { virtual, multihit, tracethrough, x, y, z, ax, ay, az, zrot }", "(string:stage name, table:properties):void",
			"Get stage properties for the given stage name as a table.", "(string:stage name):table");

	wxString name = cxt.arg(0).as_string();

	Project &prj = MainWindow::Instance().GetProject();
	Stage *stage = 0;

	for (size_t idx=0;idx<prj.StageList.size();idx++)
		if (prj.StageList[idx]->Name.CmpNoCase(name)==0)
			stage = prj.StageList[idx];

	if (!stage)
	{
		cxt.result().nullify();
		return;
	}

	if (cxt.arg_count() == 2)
	{
		lk::vardata_t &h = cxt.arg(1);
		lk::vardata_t *v = 0;

		if ( (v=h.lookup("virtual")) )
			stage->Virtual = v->deref().as_boolean();

		if ( (v=h.lookup("multihit")) )
			stage->MultiHit = v->deref().as_boolean();

		if ( (v=h.lookup("tracethrough")) )
			stage->TraceThrough = v->deref().as_boolean();

		if ( (v=h.lookup("x")) )
			stage->X = v->deref().as_number();

		if ( (v=h.lookup("y")) )
			stage->Y = v->deref().as_number();

		if ( (v=h.lookup("z")) )
			stage->Z = v->deref().as_number();

		if ( (v=h.lookup("ax")) )
			stage->AX = v->deref().as_number();

		if ( (v=h.lookup("ay")) )
			stage->AY = v->deref().as_number();

		if ( (v=h.lookup("az")) )
			stage->AZ = v->deref().as_number();

		if ( (v=h.lookup("zrot")) )
			stage->ZRot = v->deref().as_number();


		if ( StageForm *sf = MainWindow::Instance().GetGeometry()->GetStageForm( stage ) )
			sf->UpdateFromData();

		MainWindow::Instance().SetModified();
	}
	else if (cxt.arg_count() == 1)
	{
		lk::vardata_t &r = cxt.result();
		r.empty_hash();
		r.hash_item( "virtual", stage->Virtual ? 1.0 : 0.0 );
		r.hash_item( "multihit", stage->MultiHit ? 1.0 : 0.0 );
		r.hash_item( "tracethrough", stage->TraceThrough ? 1.0 : 0.0 );
		r.hash_item( "x", stage->X );
		r.hash_item( "y", stage->Y );
		r.hash_item( "z", stage->Z );
		r.hash_item( "ax", stage->AX );
		r.hash_item( "ay", stage->AY );
		r.hash_item( "az", stage->AZ );
		r.hash_item( "zrot", stage->ZRot );
	}
	else
		cxt.error("invalid number of arguments to stageopt. must be 2 or 1.");

}

static bool deref_stage( Stage **pstage, StageForm **pform )
{
	Project &prj = MainWindow::Instance().GetProject();
	for (size_t i=0;i<prj.StageList.size();i++)
	{
		if (prj.StageList[i]->Name.CmpNoCase( active_stage_name ) == 0)
		{
			*pstage = prj.StageList[i];
			*pform = MainWindow::Instance().GetGeometry()->GetStageForm( prj.StageList[i] );
			return *pstage && *pform;
		}
	}

	return false;
}

static void _addelement( lk::invoke_t &cxt )
{
	LK_DOC("addelement", "Adds one (or optionally more) elements to the active stage.", "([integer:optional number of elements]):void");

	Stage *s = 0;
	StageForm *sf = 0;
	if (!deref_stage( &s, &sf )) return;

	int n = 1;
	if (cxt.arg_count() > 0)
		n = cxt.arg(0).as_integer();
	if (n < 1) n = 1;
	if (n > 1000) n = 1000;

	sf->Append( n );
	MainWindow::Instance().SetModified();
}

static void _clearelements( lk::invoke_t &cxt )
{
	LK_DOC("clearelements", "Clears all the elements from the currently active stage.", "(void):void");

	Stage *s = 0;
	StageForm *sf = 0;
	if (!deref_stage( &s, &sf )) return;

	sf->Clear();
	MainWindow::Instance().SetModified();
}

static void _nelements( lk::invoke_t &cxt )
{
	LK_DOC("nelements", "Returns the number of elements in the current stage.", "(void):integer");

	Stage *s = 0;
	StageForm *sf = 0;
	if (!deref_stage( &s, &sf )) return;

	cxt.result().assign( (double) s->ElementList.size() );
}

static lk::vardata_t apersurfvar( char index, double params[8], const wxString &file = "" )
{
	wxString buf;
	buf += tolower(index);
	lk::vardata_t r;
	r.empty_vector();
	if (!file.empty())
	{
		r.vec_append( buf );
		r.vec_append(file);
	}
	else
	{
		r.vec_append( buf );
		for (int i=0;i<8;i++)
			r.vec_append( params[i] );
	}
	return r;
}

static void apersurfparse( std::vector<lk::vardata_t> *vec, char *index, double params[8], wxString *file)
{
	if (!vec || vec->size() < 2) return;

	wxString buf = vec->at(0).as_string();
	if (buf.length() < 0) return;

	*index = buf[0];

	if (file != 0 && vec->size() == 2)
	{
		*file = vec->at(1).as_string();
		for (int i=0;i<8;i++)
			params[i] = 0.0;
	}
	else
	{
		for (size_t i=0;i<8;i++)
		{
			if (i+1 < vec->size())
				params[i] = vec->at(i+1).as_number();
			else
				params[i] = 0.0;
		}
	}
}

static void _elementopt( lk::invoke_t &cxt )
{
	LK_DOC2("elementopt", "Two modes of operation: gets or sets various element properties for the given element in the currently active stage.  The parameters are transferred via a table with named fields { en:boolean, x:real, y:real, z:real, ax:real, ay:real, az:real, zrot:real, aper:array, surf:array, interact:string, optic:string, comment:string }.",
			"Set various element properties for the specified element in the active stage. Table fields are { en, x, y, z, ax, ay, az, zrot, aper, surf, interact, optic, comment }", "(integer:element index, table:properties):void",
			"Gets the properties of the specified element in the active stage.", "(integer:element index):table");

	Stage *s = 0;
	StageForm *sf = 0;
	if (!deref_stage( &s, &sf ))
	{
		cxt.result().nullify();
		return;
	}

	Element *e = 0;
	size_t idx = cxt.arg(0).as_unsigned();
	if (idx < s->ElementList.size())
		e = s->ElementList[idx];

	if (!e)
	{
		cxt.result().nullify();
		return;
	}

	if (cxt.arg_count() == 2)
	{
		lk::vardata_t &h = cxt.arg(1);
		lk::vardata_t *v = 0;

		if ( (v=h.lookup("en")) )
			e->Enabled = v->deref().as_boolean();

		if ( (v=h.lookup("x")) )
			e->X = v->deref().as_number();

		if ( (v=h.lookup("y")) )
			e->Y = v->deref().as_number();

		if ( (v=h.lookup("z")) )
			e->Z = v->deref().as_number();

		if ( (v=h.lookup("ax")) )
			e->AX = v->deref().as_number();

		if ( (v=h.lookup("ay")) )
			e->AY = v->deref().as_number();

		if ( (v=h.lookup("az")) )
			e->AZ = v->deref().as_number();

		if ( (v=h.lookup("zrot")) )
			e->ZRot = v->deref().as_number();

		if ( (v=h.lookup("aper")) )
			apersurfparse( v->deref().vec(), &e->ApertureIndex, e->ApertureParams, 0 );

		if ( (v=h.lookup("surf")) )
			apersurfparse( v->deref().vec(), &e->SurfaceIndex, e->SurfaceParams, &e->SurfaceFile );

		if ( (v=h.lookup("interact")) )
		{
			wxString buf = v->deref().as_string();
			if ( buf == "refraction" )
				e->InteractionType = Element::REFRACTION;
			else
				e->InteractionType = Element::REFLECTION;
		}

		if ( (v=h.lookup("optic")) )
			e->OpticName = v->deref().as_string();

		if ( (v=h.lookup("comment")) )
			e->Comment = v->deref().as_string();

		sf->UpdateFromData();
		MainWindow::Instance().SetModified();
	}
	else if (cxt.arg_count() == 1)
	{
		lk::vardata_t &r = cxt.result();
		r.empty_hash();
		r.hash_item( "en", e->Enabled ? 1.0 : 0.0 );
		r.hash_item( "x", e->X );
		r.hash_item( "y", e->Y );
		r.hash_item( "z", e->Z );
		r.hash_item( "ax", e->AX );
		r.hash_item( "ay", e->AY );
		r.hash_item( "az", e->AZ );
		r.hash_item( "zrot", e->ZRot );
		r.hash_item( "aper", apersurfvar(e->ApertureIndex, e->ApertureParams) );
		r.hash_item( "surf", apersurfvar(e->SurfaceIndex, e->SurfaceParams, e->SurfaceFile) );
		r.hash_item( "interact", e->InteractionType == Element::REFLECTION ? "reflection" : "refraction" );
		r.hash_item( "optic", e->OpticName );
		r.hash_item( "comment", e->Comment );
	}
	else
		cxt.error("invalid number of arguments. must be 2 or 1.");
}

static void _elementstats( lk::invoke_t &cxt )
{
	LK_DOC("elementstats", "Returns a table of statistical information about a particular element.", "(integer:stage index, integer:element index, integer:x bins, integer:y bins, double:dni, boolean:final rays only [, double:minx, double:maxx, double:miny, double:maxy]):table");

	int stageIdx = cxt.arg(0).as_integer();
	int elementIdx = cxt.arg(1).as_integer();
	int nbinx = cxt.arg(2).as_integer();
	int nbiny = cxt.arg(3).as_integer();
	double dni = cxt.arg(4).as_number();
	bool finalonly = cxt.arg(5).as_boolean();

	double minx, maxx, miny, maxy;
	bool autoscale = true;
	if (cxt.arg_count() == 10)
	{
		autoscale = false;
		minx = cxt.arg(6).as_number();
		maxx = cxt.arg(7).as_number();
		miny = cxt.arg(8).as_number();
		maxy = cxt.arg(9).as_number();
	}

	ElementStatistics es( MainWindow::Instance().GetProject() );
	if (es.Compute( stageIdx, elementIdx,
					nbinx, nbiny,
					autoscale, finalonly,
					dni,
					minx, miny, maxx, maxy))
	{
		cxt.result().empty_hash();;
		cxt.result().hash_item( "bin_size_x", es.binszx );
		cxt.result().hash_item( "bin_size_y", es.binszy );
		cxt.result().hash_item( "power_per_ray", es.PowerPerRay );
		cxt.result().hash_item( "peak_flux", es.PeakFlux );
		cxt.result().hash_item( "peak_flux_uncertainty", es.PeakFluxUncertainty );
		cxt.result().hash_item( "ave_flux", es.AveFlux);
		cxt.result().hash_item( "ave_flux_uncertainty", es.AveFluxUncertainty );
		cxt.result().hash_item( "min_flux", es.MinFlux );
		cxt.result().hash_item( "sigma_flux", es.SigmaFlux );
		cxt.result().hash_item( "uniformity", es.Uniformity );
		cxt.result().hash_item( "radius", es.Radius );
		cxt.result().hash_item( "num_rays", es.NumberOfRays );

		lk::vardata_t cent;
		cent.empty_vector();
		cent.vec_append( es.Centroid[0] );
		cent.vec_append( es.Centroid[1] );
		cent.vec_append( es.Centroid[2] );
		cxt.result().hash_item( "centroid", cent );

		lk::vardata_t xvals;
		xvals.empty_vector();
		for (size_t i=0;i<es.xValues.size();i++)
			xvals.vec_append( es.xValues[i] );
		cxt.result().hash_item( "xvalues", xvals );

		lk::vardata_t yvals;
		yvals.empty_vector();
		for (size_t i=0;i<es.yValues.size();i++)
			yvals.vec_append( es.yValues[i] );
		cxt.result().hash_item( "yvalues", yvals );

		lk::vardata_t flux;
		flux.empty_vector();
		for (size_t r=0;r<es.fluxGrid.nrows();r++)
		{
			lk::vardata_t row;
			row.empty_vector();
			for (size_t c=0;c<es.fluxGrid.ncols();c++)
			{
				row.vec_append( es.fluxGrid.at(r,c) );
			}
			flux.vec()->push_back( row );
		}
		cxt.result().hash_item( "flux", flux );
	}

}

static void _rayhits( lk::invoke_t &cxt )
{
	LK_DOC("rayhits", "Returns the number of ray hits on an element.", "(integer:stage index, integer:element index, [boolean: final only]):integer");
	int stageIdx = cxt.arg(0).as_integer();
	int elementIdx = cxt.arg(1).as_integer();
	bool final_only = false;
	if (cxt.arg_count() > 2)
		final_only = cxt.arg(2).as_boolean();

	if (Element *e = MainWindow::Instance().GetProject().GetElement(stageIdx, elementIdx))
		cxt.result().assign( final_only ? e->FinalRayHits : e->RayHits );
}

static void _writerayfile( lk::invoke_t &cxt )
{
	LK_DOC("writerayfile", "Writes a binary ray data file with the current trace results", "(string:file):boolean");
	cxt.result().assign( MainWindow::Instance().GetProject().Results.WriteDataFile( cxt.arg(0).as_string() ) ? 1.0 : 0.0 );
}

static lk::fcall_t *soltrace_functions()
{
	static lk::fcall_t st[] = {
		_dot,
		_euler,
		_reftoloc,
		_loctoref,
		_toloc,
		_toref,
		_transpose,
		_matvecmult,

		_workdir,
		_file_name,
		_save_project,
		_open_project,
		_clear_project,
		_writerayfile,
		_trace,
		_traceopt,
		_nintersect,
		_raydata,
		_sundata,

		_sunopt,

		_addoptic,
		_clearoptics,
		_listoptics,
		_opticopt,

		_addstage,
		_clearstages,
		_liststages,
		_stageopt,
		_activestage,

		_addelement,
		_clearelements,
		_nelements,
		_elementopt,
		_elementstats,
		_rayhits,
		0 };

	return (lk::fcall_t*)st;
}


SolTraceScriptWindowFactory::SolTraceScriptWindowFactory()
{
	// nothing to do
}

SolTraceScriptWindowFactory::~SolTraceScriptWindowFactory()
{
	// nothing to do
}

wxLKScriptWindow *SolTraceScriptWindowFactory::Create()
{
	wxLKScriptWindow *sw = new SolTraceScriptWindow( &MainWindow::Instance(), wxID_ANY );
#ifdef __WXMSW__
	sw->SetIcon( wxICON( appicon ) );
#endif	
	return sw;
}

enum { ID_VARIABLES = wxID_HIGHEST+494 };

BEGIN_EVENT_TABLE( SolTraceScriptWindow, wxLKScriptWindow )
END_EVENT_TABLE()

SolTraceScriptWindow::SolTraceScriptWindow( wxWindow *parent, int id )
	: wxLKScriptWindow( parent, id )
{
	GetEditor()->RegisterLibrary( soltrace_functions(), "SolTrace Functions");
}

void SolTraceScriptWindow::OnHelp( )
{
	MainWindow::ShowHelpTopic( "macros" );
}

void SolTraceScriptWindow::OnScriptStarted()
{
	// let the SAM window be the parent for plots
	// rather than the current toplevel window so that they
	// hang around after a script window is closed
	wxLKSetToplevelParentForPlots( &MainWindow::Instance() );

	// make sure there's no current plot active
	wxLKSetPlotTarget( NULL );
}

void SolTraceScriptWindow::OnScriptStopped()
{
	MainWindow::Instance().GetTrace()->CancelTrace();
}

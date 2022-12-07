
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


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>

#include "types.h"
#include "procs.h"
#include "interpolate.h"

#define sqr(x) ((x)*(x))

#define NLINEBUF 1024
#define READLN  my_fgets(line, NLINEBUF, fp, &line_count, -1)

#ifdef WIN32
	#define strcasecmp _stricmp
#endif

char *my_fgets(char *buf, int len, FILE *fp, int *line_count, int printchars)
{
	char *ret = fgets(buf, len, fp);
	int nch = strlen(buf);
	if (nch > 0 && buf[nch-1] == '\n')
		buf[nch-1] = 0;
	if (printchars > 0)
	{
		printf("PROCESSING (%d): '", (*line_count)++);
		int i=0;
		while (i<nch && i<printchars) printf("%c", buf[i++]);
		if (i<nch) printf(" (...)");
		printf("'\n");
	}
	return ret;
}

static bool CheckInputs( TSystem *sys )
{
	for (st_uint_t i=0;i<sys->StageList.size();i++)
	{
		TStage *stage = sys->StageList[i];
		for (st_uint_t j=0;j<stage->ElementList.size();j++)
		{
			TElement *elm = stage->ElementList[j];

			if ( (elm->SurfaceIndex == 'p' || elm->SurfaceIndex == 'P')
				&& (elm->ShapeIndex == 'l' || elm->ShapeIndex == 'L')
				&& elm->VertexCurvY != 0.0 )
			{
				sys->errlog("Parabolic surface combined with single axis curvature section aperture must have Cy = 0.  Element %d Stage %d.", j+1, i+1 );
				return false;
			}

			if (elm->SurfaceIndex == 's' || elm->SurfaceIndex == 'S')
			{
				if (elm->VertexCurvX == 0.0)
				{
					sys->errlog("Spherical surface must have a curvature greater than 0.0.  Consider using a flat surface instead. Element %d, Stage %d.", j+1, i+1);
					return false;
				}

				if (elm->VertexCurvY != 0.0 && elm->VertexCurvY != elm->VertexCurvX)
				{
					sys->errlog("Spherical surface must have Cy=0 or Cy=Cx.  Currently, Cy=%lg for Element %d, Stage %d", elm->VertexCurvY, j+1, i+1);
					return false;
				}

			}

			if (elm->SurfaceIndex == 'o' || elm->SurfaceIndex == 'O')
			{
				if (elm->Kappa == 1.0 || elm->Kappa == 0.0)
				{
					sys->errlog("Kappa = 0.0 or 1.0 defines either a parabola or a sphere.  Use the specific surface definition for either of those cases. Element %d, Stage %d.", j+1, i+1);
					return false;
				}
			}

			if (elm->ShapeIndex == 'l' || elm->ShapeIndex == 'L')
			{

				if ( elm->SurfaceIndex  == 'p' || elm->SurfaceIndex == 'P' ) {
					if (elm->VertexCurvY != 0.0)
					{
						sys->errlog("Parabolic surface Cy must equal 0 when using a single axis curvature section aperture, Element %d, Stage %d", j+1, i+1);
						return false;
					}
				} else if ( elm->SurfaceIndex == 's' || elm->SurfaceIndex == 'S' ) {
					if ( elm->VertexCurvY != 0.0 )
					{
						sys->errlog("Spherical surface Cy must equal 0 when using a single axis curvature section aperture, Element %d, Stage %d", j+1, i+1);
						return false;
					}
				}
			}
		}
	}

	return true;
}


bool InitGeometries(TSystem *sys)
{
	for (st_uint_t i=0;i<sys->StageList.size();i++)
	{
		TStage *stage = sys->StageList[i];

		//calculate Euler angles [1], [2] relative to global system from origin and aimpoint.
		double dx, dy, dz, dtot;
		dx = stage->AimPoint[0]-stage->Origin[0];
		dy = stage->AimPoint[1]-stage->Origin[1];
		dz = stage->AimPoint[2]-stage->Origin[2];
		dtot = sqrt(dx*dx + dy*dy + dz*dz);
		if (dtot == 0.0)
		{
			sys->errlog("Stage %d geometry as defined not possible.\n", i);
			return false;
		}
		
		dx = dx/dtot;
		dy = dy/dtot;
		dz = dz/dtot;
		stage->Euler[0] = atan2(dx,dz);
		stage->Euler[1] = asin(dy);
		stage->Euler[2] = stage->ZRot*ACOSM1O180;
		//calculate transformation matrices for sys and place into memory
		CalculateTransformMatrices(stage->Euler, stage->RRefToLoc, stage->RLocToRef);

		for (st_uint_t j=0;j<stage->ElementList.size();j++)
		{
			TElement *elm = stage->ElementList[j];
			
			dx = elm->AimPoint[0]-elm->Origin[0];
			dy = elm->AimPoint[1]-elm->Origin[1];
			dz = elm->AimPoint[2]-elm->Origin[2];
			dtot = sqrt(dx*dx + dy*dy + dz*dz);
			if (dtot == 0.0)
			{
				sys->errlog("Stage %d element %d geometry as defined not possible.\n", i, j);
				return false;
			}
			
			dx = dx/dtot;
			dy = dy/dtot;
			dz = dz/dtot;
			elm->Euler[0] = atan2(dx,dz);
			elm->Euler[1] = asin(dy);
			elm->Euler[2] = elm->ZRot*(ACOSM1O180);
			CalculateTransformMatrices( elm->Euler, elm->RRefToLoc, elm->RLocToRef );

			double v11,v12,v21,v22,v31,v32,v41,v42,v1m,v2m,v3m,v4m,Theta1,Theta2;
			
			switch( elm->ShapeIndex )
			{
			case 'c': case 'C': // circular aperture
				elm->ApertureArea = M_PI*sqr(elm->ParameterA/2);
				break;
			case 'h': case 'H': // hexagonal aperture
				elm->ApertureArea = 5.0*sqr(elm->ParameterA/2.0)*cos(30.0*(ACOSM1O180))*sin(30.0*(ACOSM1O180));
				break;
			case 'r': case 'R': // rectangular aperture
				elm->ApertureArea = elm->ParameterA*elm->ParameterB;
				break;
			case 't': case 'T': // equilateral triangle
				elm->ApertureArea = 0.25*sqr(elm->ParameterA*cos(30.0*(ACOSM1O180)))*sqrt(3.0);
				break;
			case 'a': case 'A': // annulus
				elm->ApertureArea =  elm->ParameterC*(ACOSM1O180)*(elm->ParameterB - elm->ParameterA);
				break;
			case 'l': case 'L': // off axis section
				elm->ApertureArea = (elm->ParameterB - elm->ParameterA)*elm->ParameterC;
				break;
			case 'i': case 'I': // irregular triangle			
				v11 =  elm->ParameterA-elm->ParameterC;
				v12 =  elm->ParameterB-elm->ParameterD;
				v21 =  elm->ParameterE-elm->ParameterC;
				v22 =  elm->ParameterF-elm->ParameterD;
				v1m = sqrt(v11*v11 + v12*v12);
				v2m = sqrt(v21*v21 + v22*v22);
				Theta1 = acos((v11*v21+v12*v22)/(v1m*v2m));
				elm->ApertureArea =  0.5*v1m*v2m*sin(Theta1); //from CRC math tables 25th edition p. 139
				break;
			case 'q': case 'Q': // irregular quadrilateral
				v11 =  elm->ParameterA-elm->ParameterC;
				v12 =  elm->ParameterB-elm->ParameterD;
				v21 =  elm->ParameterE-elm->ParameterC;
				v22 =  elm->ParameterF-elm->ParameterD;
				v31 =  elm->ParameterE-elm->ParameterG;
				v32 =  elm->ParameterF-elm->ParameterH;
				v41 =  elm->ParameterA-elm->ParameterG;
				v42 =  elm->ParameterB-elm->ParameterH;
				v1m = sqrt(v11*v11 + v12*v12);
				v2m = sqrt(v21*v21 + v22*v22);
				v3m = sqrt(v31*v31 + v32*v32);
				v4m = sqrt(v41*v41 + v42*v42);
				Theta1 = acos((v11*v21+v12*v22)/(v1m*v2m));
				Theta2 = acos((v31*v41+v32*v42)/(v3m*v4m));
				elm->ApertureArea =  0.5*(v1m*v2m*sin(Theta1) + v3m*v4m*sin(Theta2)); //from CRC math tables 25th edition p. 139
				break;
			default:
				sys->errlog("invalid shape '%c' for aperture, element %d\n", elm->ShapeIndex, j);
				return false;
			}
			
			if (!stage->Virtual)
			{
				// locate element optical properties
				TOpticalPropertySet *set = NULL;
				for (st_uint_t k=0;k<sys->OpticsList.size();k++)
					if (strcasecmp( sys->OpticsList[k]->Name.c_str(), elm->OpticName.c_str() ) == 0)
						set = sys->OpticsList[k];

				if (!set)
				{
					sys->errlog("Stage %d element %d invalid optical property set name: '%s'\n", i, j, elm->OpticName.c_str());
					return false;
				}

				elm->Optics = set;

				set->Front.RefractiveIndex[2] = set->Back.RefractiveIndex[0];
				set->Front.RefractiveIndex[3] = set->Back.RefractiveIndex[1];
				set->Back.RefractiveIndex[2] = set->Front.RefractiveIndex[0];
				set->Back.RefractiveIndex[3] = set->Front.RefractiveIndex[1];
			}
			else
				elm->Optics = NULL;

			// calculate distance from aperture plane to element origin
			AperturePlane( elm );
		}
	}

	return CheckInputs( sys );
}

bool TranslateSurfaceParams( TSystem *sys, TElement *elm, double params[8])
{
	switch( elm->SurfaceIndex )
	{
	case 's': case 'S': // spherical
		elm->VertexCurvX = params[0];

		if (elm->ShapeIndex == 'l' || elm->ShapeIndex == 'L')
			elm->VertexCurvY = 0.0;
		else
			elm->VertexCurvY = elm->VertexCurvX;

		elm->Kappa = 1.0;
		elm->SurfaceType = ( elm->VertexCurvY==0.0 ) ? 7 : 1;
		break;
	case 'p': case 'P': // parabolic
		elm->VertexCurvX = params[0];
		elm->VertexCurvY = params[1];
		elm->Kappa = 0.0;
		elm->SurfaceType = ( elm->VertexCurvY==0.0 ) ? 7 : 1;					
		break;
	case 'o': case 'O': // other than sphere or parabola
		elm->VertexCurvX = params[0];
		elm->VertexCurvY = params[0];
		elm->Kappa = params[1];
		elm->SurfaceType = ( elm->VertexCurvY==0.0 ) ? 7 : 1;					
		break;
	case 'g': case 'G': // general spencer & murty
		elm->VertexCurvX = params[0];
		elm->VertexCurvY = params[1];
		elm->Kappa = params[2];
		elm->Alpha[0] = params[3];
		elm->Alpha[1] = params[4];
		elm->Alpha[2] = params[5];
		elm->Alpha[3] = params[6];
		elm->Alpha[4] = params[7];
		elm->SurfaceType = ( elm->VertexCurvY==0.0 ) ? 7 : 1;	
		break;
	case 'f': case 'F': // flat
		elm->Alpha[0] = 0;
		elm->Alpha[1] = 0;
		elm->Alpha[2] = 1;
		elm->Alpha[3] = 0;
		elm->Alpha[4] = 0;
		elm->SurfaceType = 3;
		break;
	case 'c': case 'C': // conical
		elm->ConeHalfAngle = params[0];
		elm->SurfaceType = 1;
		break;
	case 't': case 'T': // cylindrical
		if ( (elm->ShapeIndex != 'l' && elm->ShapeIndex != 'L')
		    || elm->ParameterA != 0.0 || elm->ParameterB != 0.0 )
		{
			sys->errlog("Bad aperture or parameters of element not consistent with cylindrical surface\n");
			return false;
		}
		elm->CurvOfRev = params[0];
		elm->VertexCurvX = elm->CurvOfRev;
		elm->SurfaceType = 2;
		break;
	case 'd': case 'D': // torus ('donut')
		if ( (elm->ShapeIndex != 'a' && elm->ShapeIndex != 'A')
		    || elm->ParameterA != 0.0 || elm->ParameterB != 0.0 )
		{
			sys->errlog("Bad aperture or parameters of element not consistent with toroidal surface\n");
			return false;
		}
		elm->AnnularRadius = params[0];
		elm->CrossSectionRadius = params[1];
		elm->SurfaceType = 10;
		break;
	default:
		sys->errlog("invalid surface index: %c\n", elm->SurfaceIndex);
		return false;
	}

	return true;
}

bool ReadSurfaceFile( const char *file, TElement *elm , TSystem *sys)
{
	char ext[16], *p;
	int line_count = 1;
	char line[NLINEBUF];

	p = (char*)strrchr(file, '.');
	if (!p)
	{
		sys->errlog("Could not determine surface file type: '%s'\n", file);
		return false;
	}
	strncpy( ext, p+1, 15 );
	p = ext;

	int len = strlen(ext);
	for (int i=0;i<len;i++)
		ext[i] = tolower(ext[i]);

	FILE *fp = NULL;
	fp = fopen(file, "r");
	if (!fp)
	{
		sys->errlog("Could not open surface file for read: %s\n", file);
		return false;
	}

	if (strcmp(ext, "mon")==0)
	{
		int Order = 0;
		READLN; Order = atoi( line );

		elm->BCoefficients.resize( Order+1, Order+1 );

		for (int k=0;k<=Order;k++)
		{
			for (int m=0;m<=k;m++)
			{
				READLN; 
				elm->BCoefficients.at( k, m ) = atof( line );
			}
		}

		elm->FitOrder = Order;
		elm->SurfaceIndex = 'm';
		elm->SurfaceType = 6;
	}
	else if (strcmp(ext, "sht")==0)
	{
		READLN; // skip first line (file name)
		READLN; sscanf(line, "%lg %lg %lg", 
			&elm->VSHOTRadius, &elm->VSHOTFocLen, &elm->VSHOTTarDis);

		int Order=0, NumPoints=0, idum;
		READLN; sscanf(line, "%d %d %d", &idum, &Order, &NumPoints);
		READLN; sscanf(line, "%lg %lg", &elm->VSHOTRMSSlope, &elm->VSHOTRMSScale);

		elm->FitOrder = Order;

		elm->BCoefficients.resize(Order+1, Order+1);

		for (int k=0;k<=Order;k++)
		{
			for (int m=0;m<=k;m++)
			{
				READLN;
				elm->BCoefficients.at( k, m) = atof(line);
			}
		}

		elm->VSHOTData.resize( NumPoints, 5 );
		for (int i=0;i<NumPoints;i++)
		{
			double a,b,c,d,e;
			READLN;
			sscanf(line, "%lg %lg %lg %lg %lg", &a, &b, &c, &d, &e );
			elm->VSHOTData.at(i,0) = a;
			elm->VSHOTData.at(i,1) = b;
			elm->VSHOTData.at(i,2) = c;
			elm->VSHOTData.at(i,3) = d;
			elm->VSHOTData.at(i,4) = e;
		}

		elm->SurfaceIndex = 'v';
		elm->SurfaceType = 5;

	}
	else if (strcmp(ext, "ply")==0)
	{
		int Order = 0;
		READLN; Order = atoi(line);

		elm->PolyCoeffs.resize( Order + 1 );

		for (int k=0;k<=Order;k++)
		{
			READLN;
			elm->PolyCoeffs[k] = atof( line );
		}

		elm->FitOrder = Order;
		elm->SurfaceIndex = 'r';
		elm->SurfaceType = 8;
	}
	else if (strcmp(ext, "csi")==0)
	{
		int NPoints = 0;
		READLN; NPoints = atoi(line);

		elm->CubicSplineXData.resize(NPoints);
		elm->CubicSplineYData.resize(NPoints);
		elm->CubicSplineY2Data.resize(NPoints);

		double x, y;
		for (int k=0;k<NPoints;k++)
		{
			READLN; sscanf(line, "%lg %lg", &x, &y);
			elm->CubicSplineXData[k] = x;
			elm->CubicSplineYData[k] = y;
		}
		READLN; sscanf(line, "%lg %lg", &x, &y);
		elm->CubicSplineDYDXbc1 = x;
		elm->CubicSplineDYDXbcN = y;

		spline( elm->CubicSplineXData, elm->CubicSplineYData,
			NPoints, x, y, elm->CubicSplineY2Data );

		if ( ((elm->ShapeIndex != 'a') && (elm->ShapeIndex != 'A') && (elm->ShapeIndex != 'l') && (elm->ShapeIndex != 'L'))
			|| (elm->ParameterA < elm->CubicSplineXData[0])
			|| (elm->ParameterB > elm->CubicSplineXData[elm->CubicSplineXData.size()-1])
			)
		{
			sys->errlog("Error: Element uses cubic spline interpolation:"
				"\tMake sure aperture is type 'a' or 'l' and that "
				"1st and 2nd parameters are not < and not > the "
				"1st and last X values in file respectively.\n");
			fclose(fp);
			return false;
		}
		elm->SurfaceIndex = 'i';
		elm->SurfaceType = 9;
	}
	else if (strcmp(ext, "fed") == 0)
	{
		int NumPoints = 0;
		READLN; // skip FE file name
		READLN; NumPoints = atoi(line);

		//elm->FEData.x.resize(NumPoints, VectDoub(2, 0.));
		//elm->FEData.y.resize(NumPoints, 0.);
		elm->FEData.nodes.resize(NumPoints, VectDoub(3));

		double xmax, xmin, ymax, ymin;
		xmax = ymax = -std::numeric_limits<double>::max();
		xmin = ymin = std::numeric_limits<double>::max();

		
		for (int i=0;i<NumPoints;i++)
		{
			double x,y,z;
			READLN; sscanf(line, "%lg %lg %lg", &x, &y, &z);
			elm->FEData.nodes.at(i).at(0) = x;
			elm->FEData.nodes.at(i).at(1) = y;
			elm->FEData.nodes.at(i).at(2) = z;

			//track largest/smallest
			xmax = x > xmax ? x : xmax;
			ymax = y > ymax ? y : ymax;
			xmin = x < xmin ? x : xmin;
			ymin = y < ymin ? y : ymin;
		}

		KDLayoutData node_ld;
		node_ld.xlim[0] = xmin;
		node_ld.xlim[1] = xmax;
		node_ld.ylim[0] = ymin;
		node_ld.ylim[1] = ymax;
		double rapprox = 1.5*std::sqrt(((xmax - xmin) * (ymax - ymin)) / (double)NumPoints);
		node_ld.min_unit_dx = node_ld.min_unit_dy = rapprox;

		elm->FEData.create_mesh(node_ld);

		//Load node objects into the mesh
		for (int i = 0; i < NumPoints; i++)
		{
			VectDoub* v = &elm->FEData.nodes.at(i);
			elm->FEData.add_object((void*)v, v->at(0), v->at(1));
		}
		elm->FEData.add_neighborhood_data();


		elm->SurfaceIndex = 'e';
		elm->SurfaceType = 4;
	}
	else
	{
		fclose(fp);
		sys->errlog("Surface file type extension unknown: '%s'\n", ext);
		return false;
	}
	
	elm->SurfaceFile = file;

	fclose(fp);
	return true;
}

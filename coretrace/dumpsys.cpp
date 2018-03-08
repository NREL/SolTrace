
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

#include "types.h"
#include "procs.h"


void DumpVec3( FILE *fp, double vec[3], char *caption )
{
	fprintf(fp, "%s %lg %lg %lg\n", caption, vec[0], vec[1], vec[2]);
}

void DumpMat3x3( FILE *fp, double mat[3][3], char *caption )
{
	fprintf(fp, "%s\n", caption);
	for (int i=0;i<3;i++)
		fprintf(fp, "\t%lg %lg %lg\n", mat[i][0], mat[i][1], mat[i][2]);	
}


void DumpSun(FILE *fp, TSun *sun)
{
	fprintf(fp, "--- BEGIN SUN ---\n");
	fprintf(fp, "SHAPE %c SIGMA %lg POINT %d\n", sun->ShapeIndex, sun->Sigma, sun->PointSource?1:0);
	if (sun->SunShapeAngle.size() != sun->SunShapeIntensity.size())
	{
		printf("sun error, unequal shape data array lengths\n");
		return;
	}
	
	fprintf(fp, "DATA %d\n", sun->SunShapeAngle.size());
	for (int i=0;i<sun->SunShapeAngle.size();i++)
		fprintf(fp, "\t%lg %lg\n", sun->SunShapeAngle[i], sun->SunShapeIntensity[i]);
	//fprintf(fp, "MAXANGLE %lg MAXINTENSITY %lg\n", sun->MaxAngle, sun->MaxIntensity);
	
	DumpVec3(fp, sun->Origin, "ORIGIN");
	//DumpVec3(fp, sun->Euler, "EULER");
	//DumpMat3x3(fp, sun->RRefToLoc, "RREFTOLOC");
	//DumpMat3x3(fp, sun->RLocToRef, "RLOCTOREF");
	
	fprintf(fp, "MAXRAD %lg XCM %lg YCM %lg\n", sun->MaxRad, sun->Xcm, sun->Ycm);
	fprintf(fp, "MINX %lg MAXX %lg MINY %lg MAXY %lg\n", sun->MinXSun, sun->MaxXSun, sun->MinYSun, sun->MaxYSun);
	
	fprintf(fp, "--- END SUN ---\n\n");
		
}

void DumpOptics(FILE *fp, TOpticalProperties *opt, char *name)
{
	fprintf(fp, "\tOPTICAL PROPERTIES: %s\n", name);
	fprintf(fp, "\t  OPTICSURF %d APERSTOPGRATING %d DIFFRACTORD %d\n",
			opt->OpticSurfNumber, opt->ApertureStopOrGratingType, opt->DiffractionOrder);
	fprintf(fp, "\t  REFLECTIVITY %lg TRANSMISSIVITY %lg SLOPEERR %lg SPECERR %lg DISTTYPE %c\n",
			opt->Reflectivity, opt->Transmissivity, opt->RMSSlopeError, opt->RMSSpecError, opt->DistributionType);
	fprintf(fp, "\t  REFRACTIVEINDEX %lg %lg %lg %lg\n", opt->RefractiveIndex[0], opt->RefractiveIndex[1],
			 opt->RefractiveIndex[2],  opt->RefractiveIndex[3] );
	fprintf(fp, "\t  AB12 %lg %lg %lg %lg\n", opt->AB12[0], opt->AB12[1],
			opt->AB12[2], opt->AB12[3] );
}

void DumpElement(FILE *fp, TElement *elm, int idx_s, int idx_e)
{
	fprintf(fp, "--- BEGIN ELEMENT %d STAGE %d ---\n", idx_e, idx_s);
	fprintf(fp, "\tENABLED %d XYZ %lg %lg %lg AIM %lg %lg %lg ZR %lg\n", elm->Enabled?1:0,
		elm->Origin[0],elm->Origin[1],elm->Origin[2],
		elm->AimPoint[0],elm->AimPoint[1],elm->AimPoint[2],
		elm->ZRot
		);

	//DumpVec3(fp, elm->Euler, "EULER");
	//DumpMat3x3(fp, elm->RRefToLoc, "RREFTOLOC");
	//DumpMat3x3(fp, elm->RLocToRef, "RLOCTOREF");
	
	fprintf(fp, "\tAPERTURE %c PARAMS %lg %lg %lg %lg %lg %lg %lg %lg\n", elm->ShapeIndex,
		elm->ParameterA, elm->ParameterB, elm->ParameterC, elm->ParameterD, 
		elm->ParameterE, elm->ParameterF, elm->ParameterG, elm->ParameterH );
	//fprintf(fp, "APERTURE AREA %lg ZAPER %lg\n", elm->ApertureArea, elm->ZAperture);
	
	fprintf(fp, "\tSURFACE INDEX %c TYPE %d FILE '%s'\n", elm->SurfaceIndex, elm->SurfaceType, elm->SurfaceFile.c_str());
	fprintf(fp, "\tVCURVX %lg VCURVY %lg KAPPA %lg\n", elm->VertexCurvX, elm->VertexCurvY, elm->Kappa);
	fprintf(fp, "\tALPHA %lg %lg %lg %lg %lg\n", elm->Alpha[0], elm->Alpha[1], elm->Alpha[2], elm->Alpha[3], elm->Alpha[4]);
	fprintf(fp, "\tANNULARRAD %lg CROSSSECRAD %lg CONEHALFANG %lg CURVOFREV %lg\n", elm->AnnularRadius, elm->CrossSectionRadius,
		elm->ConeHalfAngle, elm->CurvOfRev);
	
	fprintf(fp, "\tFITORDER %d\n", elm->FitOrder);
	
	fprintf(fp, "\tOPTICNAME '%s' INTERACTIONTYPE %d\n", elm->OpticName.c_str(), elm->InteractionType);
	
	fprintf(fp, "--- END ELEMENT %d STAGE %d ---\n\n", idx_e, idx_s);
	
	
}

void DumpStage(FILE *fp, TStage *stage, int nstage)
{
	
	
	fprintf(fp, "BEGIN STAGE %d ELEMENTS %d XYZ %lg %lg %lg AIM %lg %lg %lg ZR %lg MULTIHITS %d VIRTUAL %d\n", nstage,
		stage->ElementList.size(),
		stage->Origin[0],stage->Origin[1],stage->Origin[2],
		stage->AimPoint[0],stage->AimPoint[1],stage->AimPoint[2],
		stage->ZRot,
		stage->MultiHitsPerRay?1:0, stage->Virtual?1:0);
	//DumpMat3x3(fp, sys->RRefToLoc, "RREFTOLOC");
	//DumpMat3x3(fp, sys->RLocToRef, "RLOCTOREF");
	
	for (int i=0;i<stage->ElementList.size();i++)
		DumpElement( fp, stage->ElementList[i], 0, i );
		
	fprintf(fp, "END STAGE %d.\n", nstage);
	
}

bool DumpSystem(const char *file, TSystem *sys)
{
	FILE *fp = fopen(file, "w");
	if (!fp)
	{
		printf("could not dump system to: %s\n", file);
		return false;
	}
	
	fprintf(fp, "CONCENTRATOR SYSTEM DUMP\n\n");
	DumpSun( fp, &sys->Sun );

	fprintf(fp, "\nOPTICAL SETS\n");
	for (int i=0;i<sys->OpticsList.size();i++)
	{
		fprintf(fp,"OPTIC SET %d '%s'\n", i, sys->OpticsList[i]->Name.c_str());
		DumpOptics(fp, &sys->OpticsList[i]->Front, "FRONT");
		DumpOptics(fp, &sys->OpticsList[i]->Back, "BACK");
	}

	fprintf(fp, "\nSTAGE LIST %d\n", sys->StageList.size());

	for (int i=0;i<sys->StageList.size();i++)
		DumpStage(fp, sys->StageList[i], i);
		
	fprintf(fp, "\nEND.\n");
	
	fclose(fp);
	return true;
}


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


#include <math.h>
#include <stdio.h>

#include "types.h"
#include "procs.h"

#define sqr(x) (x*x)

bool SunToPrimaryStage(
				TSystem *System,
				TStage *Stage,
				TSun *Sun,
				double PosSunStage[3])
{

/*{Purpose: To compute the sun position within primary sage and the maximum radius of a cicle seen from sun which encircles
all elements  within the primary stage.  Used for genenating rays from sun.   modified on 09/26/05 to establish rectangular
region of interest - more efficient.
     Input - Sys     = Primary stage
             Sun       = Sun description block
     Output - PosSunStage  = position of sun in primary stage coordinate system
              Sun.MaxRad =    maximum radius of a cicle seen from sun which encircles
                          all elements within the primary stage  relative to center of mass of all elements in that stage
              Sun.Xcm
              Sun.Ycm    =  center of mass of all elements in primary stage as seen in sun coordinate system
}*/


	double dx=0, dy=0, dz=0, dtot=0;
	double CosSunGlob[3] = {0.0, 0.0, 0.0};
	double PosSunGlob[3] = {0.0, 0.0, 0.0};
	double CosSunStage[3] = {0.0, 0.0, 0.0};
	
	st_uint_t i = 0;
	double x=0,y=0,radius=0;
	double Origin[3] = {0.0, 0.0, 0.0};
	double CosDum[3] = {0.0, 0.0, 0.0};
	double PosLoc[3] = {0.0, 0.0, 0.0};
	double CosLoc[3] = {0.0, 0.0, 0.0};
	double RRefToLoc[3][3] = {  {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0} };
	double radius1 = 0.0,radius2 = 0.0,radius3 = 0.0,radius4 = 0.0,radiustemp = 0.0;
	double Xsum = 0.0,Ysum = 0.0, xminsun = 0.0,yminsun = 0.0,xmaxsun = 0.0,ymaxsun = 0.0;
	double XLegofRadius = 0.0;
	
	//PosSunGlob[0] = Sun.Origin[0];//Position of sun coord. system origin in global system
	//PosSunGlob[1] = Sun.Origin[1]; //changed 5/1/00 to place sun at primary stage origin; direction vector
	//PosSunGlob[2] = Sun.Origin[2]; //calculated below from difference between entered sun position and global
	PosSunGlob[0] = Stage->Origin[0]; //origin
	PosSunGlob[1] = Stage->Origin[1];
	PosSunGlob[2] = Stage->Origin[2];

	//First calculate direction cosines of sun z-axis in global coord. system
	dx = 0.0 - Sun->Origin[0];  //changed 5/1/00 to tie the sun direction to global coordinate system origin
	dy = 0.0 - Sun->Origin[1];  //for any stage; not different for each stage.  this represents reality because
	dz = 0.0 - Sun->Origin[2];  //sun is essentially inifinitely far away.
	dtot = sqrt(dx*dx + dy*dy + dz*dz);

	if (dtot == 0.0)
	{
		// flag error somehow?
		System->errlog("error calculating sun position in primary stage, dtot = 0.0\n");
		return false;
	}

	dx = dx/dtot; //unit vector in global coord.system of sun coord. system z axis.
	dy = dy/dtot;
	dz = dz/dtot;

	CosSunGlob[0] = dx;  //direction cosines of sun Z-axis in global system.
	CosSunGlob[1] = dy;
	CosSunGlob[2] = dz;

	//Transform sun direction vector to Stage system; CosSunStage is dir cosines of sun ray in Stage coord. system
	//PosSunStage is position of sun coord. system origin in Stage system
	TransformToLocal(PosSunGlob, CosSunGlob, Stage->Origin, Stage->RRefToLoc, PosSunStage, CosSunStage);

	Sun->Euler[0] = atan2(CosSunStage[0],CosSunStage[2]);   //Euler angles relating sun to Stage system
	Sun->Euler[1] = asin(CosSunStage[1]);
	Sun->Euler[2] = 0.0;

/*     {Now we have the Euler angles from Stage to the sun coordinate system.  We have to now transform the
      element locations in the stage system to the sun coordinate system and find the smallest circle in the
      xy plane of the sun system that completely encompasses the projected images of the elements onto that plane}*/
	
	Origin[0] = 0.0;  //Origin of transformed system and stage system the same
	Origin[1] = 0.0;
	Origin[2] = 0.0;
	
	CosDum[0] = 0.0; //direction cosines not important; only interested in point locations
	CosDum[1] = 0.0;
	CosDum[2] = 1.0;
	
	Sun->MaxRad = 0.0;
	Sun->Xcm = 0.0;
	Sun->Ycm = 0.0;
	Sun->MaxXSun = -1.0e20;
	Sun->MinXSun =  1.0e20;
	Sun->MaxYSun = -1.0e20;
	Sun->MinYSun = 1.0e20;


	CalculateTransformMatrices(Sun->Euler, RRefToLoc, Sun->RLocToRef);

     //{Now calculate center of mass of projected distribution. Added 09/26/05}
	Xsum = 0.0;
	Ysum = 0.0;
	for (i=0;i<Stage->ElementList.size();i++)
	{
		if ( !Stage->ElementList[i]->Enabled ) continue;
		TransformToLocal(Stage->ElementList[i]->Origin, CosDum, Origin, RRefToLoc, PosLoc, CosLoc);
		//Now have PosLoc which is the projected position of element[i] in xy plane of sun coord. system
		Xsum = Xsum + PosLoc[0];
		Ysum = Ysum + PosLoc[1];
	}
	Sun->Xcm = Xsum/Stage->ElementList.size(); //center of mass of distribution of element locations as projected in sun coord.
	Sun->Ycm = Ysum/Stage->ElementList.size(); //system.   Added 09/26/05

	size_t nelements = 0;

	for (i=0;i<Stage->ElementList.size();i++)
	{
		if ( !Stage->ElementList[i]->Enabled ) continue;

		TransformToLocal(Stage->ElementList[i]->Origin, CosDum, Origin, RRefToLoc, PosLoc, CosLoc);
		//Now have PosLoc which is the projected position of element[i] in xy plane of sun coord. system
		x = PosLoc[0] - Sun->Xcm;   //changes origin to center of mass of all elements  09/26/05
		y = PosLoc[1] - Sun->Ycm;
		radius = sqrt(x*x + y*y);
		
		xminsun = PosLoc[0];
		xmaxsun = PosLoc[0];
		yminsun = PosLoc[1];
		ymaxsun = PosLoc[1];

        //save the projected position of the element on the sun coordinate plane
        Stage->ElementList[i]->PosSunCoords[0] = PosLoc[0];
        Stage->ElementList[i]->PosSunCoords[1] = PosLoc[1];
        Stage->ElementList[i]->PosSunCoords[2] = PosLoc[2];

       //Add radius of element circle of interest - different radius for each shape: circular, hexagonal, rectangular, triangular, annular, off-axis rectangle
		if (Stage->ElementList[i]->ShapeIndex == 'c' || Stage->ElementList[i]->ShapeIndex == 'C')
		{
			radius = radius + Stage->ElementList[i]->ParameterA/2.0;
			xminsun = xminsun - Stage->ElementList[i]->ParameterA/2.0;
			yminsun = yminsun - Stage->ElementList[i]->ParameterA/2.0;
			xmaxsun = xmaxsun + Stage->ElementList[i]->ParameterA/2.0;
			ymaxsun = ymaxsun + Stage->ElementList[i]->ParameterA/2.0;
		}
		if (Stage->ElementList[i]->ShapeIndex == 'h' || Stage->ElementList[i]->ShapeIndex == 'H')
		{
			radius = radius + Stage->ElementList[i]->ParameterA/2.0;
			xminsun = xminsun - Stage->ElementList[i]->ParameterA/2.0;
			yminsun = yminsun - Stage->ElementList[i]->ParameterA/2.0;
			xmaxsun = xmaxsun + Stage->ElementList[i]->ParameterA/2.0;
			ymaxsun = ymaxsun + Stage->ElementList[i]->ParameterA/2.0;
		}
		if (Stage->ElementList[i]->ShapeIndex == 'r' || Stage->ElementList[i]->ShapeIndex == 'R')
		{
			radius = radius + sqrt((sqr(Stage->ElementList[i]->ParameterA) + sqr(Stage->ElementList[i]->ParameterB))/4.0);
			xminsun = xminsun - sqrt((sqr(Stage->ElementList[i]->ParameterA) + sqr(Stage->ElementList[i]->ParameterB))/4.0);
			yminsun = yminsun - sqrt((sqr(Stage->ElementList[i]->ParameterA) + sqr(Stage->ElementList[i]->ParameterB))/4.0);
			xmaxsun = xmaxsun + sqrt((sqr(Stage->ElementList[i]->ParameterA) + sqr(Stage->ElementList[i]->ParameterB))/4.0);
			ymaxsun = ymaxsun + sqrt((sqr(Stage->ElementList[i]->ParameterA) + sqr(Stage->ElementList[i]->ParameterB))/4.0);
		}
		if (Stage->ElementList[i]->ShapeIndex == 't' || Stage->ElementList[i]->ShapeIndex == 'T')
		{
			radius = radius + Stage->ElementList[i]->ParameterA/(2.0*cos(30.0*(ACOSM1O180)));
			xminsun = xminsun - Stage->ElementList[i]->ParameterA/(2.0*cos(30.0*(ACOSM1O180)));
			yminsun = yminsun - Stage->ElementList[i]->ParameterA/(2.0*cos(30.0*(ACOSM1O180)));
			xmaxsun = xmaxsun + Stage->ElementList[i]->ParameterA/(2.0*cos(30.0*(ACOSM1O180)));
			ymaxsun = ymaxsun + Stage->ElementList[i]->ParameterA/(2.0*cos(30.0*(ACOSM1O180)));
		}
		if (Stage->ElementList[i]->ShapeIndex == 'a' || Stage->ElementList[i]->ShapeIndex == 'A')
		{
			radius = radius + Stage->ElementList[i]->ParameterB;
			xminsun = xminsun - Stage->ElementList[i]->ParameterB;
			yminsun = yminsun - Stage->ElementList[i]->ParameterB;
			xmaxsun = xmaxsun + Stage->ElementList[i]->ParameterB;
			ymaxsun = ymaxsun + Stage->ElementList[i]->ParameterB;
		}
		
		if (Stage->ElementList[i]->ShapeIndex == 'l' || Stage->ElementList[i]->ShapeIndex == 'L')
		{
			if (fabs(Stage->ElementList[i]->ParameterB) >= fabs(Stage->ElementList[i]->ParameterA))          //change made on 02-12-09  replaced above with following block
				XLegofRadius = Stage->ElementList[i]->ParameterB;
			else
				XLegofRadius = Stage->ElementList[i]->ParameterA;
				
			radius =  radius  + sqrt(sqr(XLegofRadius) + 0.25*Stage->ElementList[i]->ParameterC*Stage->ElementList[i]->ParameterC);
			xminsun = xminsun - sqrt(sqr(XLegofRadius) + 0.25*Stage->ElementList[i]->ParameterC*Stage->ElementList[i]->ParameterC);
			yminsun = yminsun - sqrt(sqr(XLegofRadius) + 0.25*Stage->ElementList[i]->ParameterC*Stage->ElementList[i]->ParameterC);
			xmaxsun = xmaxsun + sqrt(sqr(XLegofRadius) + 0.25*Stage->ElementList[i]->ParameterC*Stage->ElementList[i]->ParameterC);
			ymaxsun = ymaxsun + sqrt(sqr(XLegofRadius) + 0.25*Stage->ElementList[i]->ParameterC*Stage->ElementList[i]->ParameterC);
		}
		
       //****************************************************************
		if (Stage->ElementList[i]->ShapeIndex == 'i' || Stage->ElementList[i]->ShapeIndex == 'I')
		{
			radius1 = sqrt(sqr(Stage->ElementList[i]->ParameterA) + sqr(Stage->ElementList[i]->ParameterB));
			radius2 = sqrt(sqr(Stage->ElementList[i]->ParameterC) + sqr(Stage->ElementList[i]->ParameterD));
			radius3 = sqrt(sqr(Stage->ElementList[i]->ParameterE) + sqr(Stage->ElementList[i]->ParameterF));
			radiustemp = radius1;
			if (radius2 > radiustemp) radiustemp = radius2;
			if (radius3 > radiustemp) radiustemp = radius3;
			radius = radius + radiustemp;
			xminsun = xminsun - radiustemp;
			yminsun = yminsun - radiustemp;
			xmaxsun = xmaxsun + radiustemp;
			ymaxsun = ymaxsun + radiustemp;
		}
		if (Stage->ElementList[i]->ShapeIndex == 'q' || Stage->ElementList[i]->ShapeIndex == 'Q')
		{
			radius1 = sqrt(sqr(Stage->ElementList[i]->ParameterA) + sqr(Stage->ElementList[i]->ParameterB));
			radius2 = sqrt(sqr(Stage->ElementList[i]->ParameterC) + sqr(Stage->ElementList[i]->ParameterD));
			radius3 = sqrt(sqr(Stage->ElementList[i]->ParameterE) + sqr(Stage->ElementList[i]->ParameterF));
			radius4 = sqrt(sqr(Stage->ElementList[i]->ParameterG) + sqr(Stage->ElementList[i]->ParameterH));
			radiustemp = radius1;
			if ( radius2 > radiustemp ) radiustemp = radius2;
			if ( radius3 > radiustemp ) radiustemp = radius3;
			if ( radius4 > radiustemp ) radiustemp = radius4;
			radius = radius + radiustemp;
			xminsun = xminsun - radiustemp;
			yminsun = yminsun - radiustemp;
			xmaxsun = xmaxsun + radiustemp;
			ymaxsun = ymaxsun + radiustemp;
		}

		if ( radius > Sun->MaxRad ) Sun->MaxRad = radius;     //establishes a circular region
		
		if ( xminsun < Sun->MinXSun ) Sun->MinXSun = xminsun;  //restablishes a rectangular region instead of a circular region
		if ( xmaxsun > Sun->MaxXSun ) Sun->MaxXSun = xmaxsun;  // Added 09/26/05
		if ( yminsun < Sun->MinYSun ) Sun->MinYSun = yminsun;
		if ( ymaxsun > Sun->MaxYSun ) Sun->MaxYSun = ymaxsun;


		nelements++;
	}

	if (nelements == 0)
		System->errlog("error calculating sun position in primary stage because no elements were enabled");

	return (nelements > 0);
}
//End of Procedure--------------------------------------------------------------

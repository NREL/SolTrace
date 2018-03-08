
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

#include "types.h"
#include "procs.h"

#define sqr(x) (x*x)

void SpencerandMurtySurfaceClosedForm( 
			TElement *Element,
			double PosLoc[3],
			double CosLoc[3],
			double PosXYZ[3],
			double DFXYZ[3],
			double *PathLength,
			int *ErrorFlag )
{
//handles spheres, parabolas, hyperbolas, ellipsoids
	double t1=0,t2=0,A=0,B=0,C=0,slopemag=0;
	double LineFocusFactor=0;

	*ErrorFlag = 0;
	
	if (Element->VertexCurvY == 0.0)
		LineFocusFactor = 0.0; //line focus section with curvature in one direction
	else
		LineFocusFactor = 1.0; // not line focus but surface of revolution
	
	A = 0.5*Element->VertexCurvX*(CosLoc[0]*CosLoc[0]+LineFocusFactor*CosLoc[1]*CosLoc[1]+CosLoc[2]*CosLoc[2]*Element->Kappa);
	B = Element->VertexCurvX*(PosLoc[0]*CosLoc[0]+LineFocusFactor*PosLoc[1]*CosLoc[1]+PosLoc[2]*CosLoc[2]*Element->Kappa)-CosLoc[2];
	C = 0.5*Element->VertexCurvX*(PosLoc[0]*PosLoc[0]+LineFocusFactor*PosLoc[1]*PosLoc[1]+PosLoc[2]*PosLoc[2]*Element->Kappa)-PosLoc[2];

	if (sqr(B) > 4.0*A*C)
	{
		t1 = (-B + sqrt(sqr(B)-4.0*A*C))/(2.0*A);
		t2 = (-B - sqrt(sqr(B)-4.0*A*C))/(2.0*A);
		if (t2 > 0)    //initial ray location outside surface
		{
			PosXYZ[0] = PosLoc[0] + t2*CosLoc[0];
			PosXYZ[1] = PosLoc[1] + t2*CosLoc[1];
			PosXYZ[2] = PosLoc[2] + t2*CosLoc[2];
			*PathLength = t2;
			if (PosXYZ[2] > Element->ZAperture)  //makes sure to get shortest ray path on valid side of closed surface;not sure what's going on here 06-12-07
			{
				PosXYZ[0] = PosLoc[0] + t1*CosLoc[0];
				PosXYZ[1] = PosLoc[1] + t1*CosLoc[1];
				PosXYZ[2] = PosLoc[2] + t1*CosLoc[2];
				*PathLength = t1;
			}
			goto Label_100;
		}
		if (t2 == 0)   //initial ray location at surface
		{
			PosXYZ[0] = PosLoc[0] + t1*CosLoc[0];
			PosXYZ[1] = PosLoc[1] + t1*CosLoc[1];
			PosXYZ[2] = PosLoc[2] + t1*CosLoc[2];
			*PathLength = t1;
			goto Label_100;
		}
		if (t2 < 0 && t1 > 0)     //initial ray location inside surface
		{
			PosXYZ[0] = PosLoc[0] + t1*CosLoc[0];
			PosXYZ[1] = PosLoc[1] + t1*CosLoc[1];
			PosXYZ[2] = PosLoc[2] + t1*CosLoc[2];
			*PathLength = t1;
			goto Label_100;
		}
		if (t1 <= 0)
		{
			*PathLength = t1; //ray heading away from surface
			*ErrorFlag = 1;
			return;
		}
	}
	else
	{
		*PathLength = 0.0; //ray tangent or missed
		*ErrorFlag = 1;
		return;
	}

Label_100:
	slopemag = sqrt(sqr(-Element->VertexCurvX*PosXYZ[0])+sqr(-Element->VertexCurvY*PosXYZ[1])+sqr(1.0-Element->Kappa*Element->VertexCurvX*PosXYZ[2]));

	DFXYZ[0] = -Element->VertexCurvX*PosXYZ[0]/slopemag;
	DFXYZ[1] = -Element->VertexCurvY*PosXYZ[1]/slopemag;
	DFXYZ[2] = (1.0-Element->Kappa*Element->VertexCurvX*PosXYZ[2])/slopemag;
}
//end of procedure--------------------------------------------------------------

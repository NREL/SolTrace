
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


void NewZStartforCubicSplineSurf(
			double CRadius,
			double PosLoc[3],
			double CosLoc[3],
			char AperShapeIndex,
			double *NewZStart,
			double *PLength,
			int *EFlag )
{
//This procedure calculates a new ZStart starting point for rays which strike cubic spline surfaces so obliquely that they intersect the
//default ZA planes (i.e. the origin or the aperture) either inside or outside the cubic spline data set.  The new plane is found by
//finding the intersection point of the ray with a virtual cylinder (or plane) with radius (or distance from origin) = inner or outer limit of cubic spline data set.
// Z value of this point becomes the new ZStart plane for the Newton-Raphson interation process.  If EFlag = 0 then success, if EFlag =
//1, ray heading away from virtual cylinder (or plane). If EFlag = 2, then ray totally misses virtual cylinder (or plane).  These last two should never happen.

	double t1 = 0.0,t2 = 0.0,A = 0.0,B = 0.0,C = 0.0;
	*EFlag = 0;
	
	switch (AperShapeIndex)
	{
	case 'a':
	case 'A'://rotationally symmetric annulus
	
			A = CosLoc[0]*CosLoc[0] + CosLoc[1]*CosLoc[1];
			B = 2.0*(PosLoc[0]*CosLoc[0] + PosLoc[1]*CosLoc[1]);
			C = PosLoc[0]*PosLoc[0] + PosLoc[1]*PosLoc[1] - CRadius*CRadius ;
						
			if (B*B > 4.0*A*C)
			{
				t1 = (-B + sqrt(B*B-4.0*A*C))/(2.0*A);
				t2 = (-B - sqrt(B*B-4.0*A*C))/(2.0*A);
				if (t2 > 0)    //initial ray location outside surface
				{
					*NewZStart = PosLoc[2] + t2*CosLoc[2];
					*PLength = t2;
					return;
				}
				
				if (t2 == 0)   //initial ray location at surface
				{
					*NewZStart = PosLoc[2] + t1*CosLoc[2];
					*PLength = t1;
					return;
				}
				
				if (t2 < 0 && t1 > 0) //initial ray location inside surface
				{
					*NewZStart = PosLoc[2] + t1*CosLoc[2];
					*PLength = t1;
					return;
				}
				
				if (t1 <= 0)
				{
					*PLength = t1; //ray heading away from surface
					*EFlag = 1;
					return;
				}
			}
			else
			{
			  *PLength = 0.0; //ray tangent or missed
			  *EFlag = 2;
				return;
			}
		break;

	case 'l':
	case 'L':  //single axis curvature section
	
			t1 = (CRadius-PosLoc[0])/CosLoc[0];
			if (t1 >= 0.0)
			{
				*NewZStart = PosLoc[2] + t1*CosLoc[2];
				*PLength = t1;
				return;
			}
			else
			{
				*PLength = 0.0;
				*EFlag = 1;
				return;
			}
		break;
	}
}
//end of procedure--------------------------------------------------------------

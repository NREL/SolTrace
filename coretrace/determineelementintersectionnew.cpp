
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
#include <math.h>
#include "types.h"
#include "procs.h"

void DetermineElementIntersectionNew(
			TElement *Element,
			double PosRayIn[3],
			double CosRayIn[3],
			double PosRayOut[3],
			double CosRayOut[3],
			double DFXYZ[3],
			double *PathLength,
			int *ErrorFlag,
			int *Intercept,
			int *BacksideFlag )
{
	double r = 0.0, Ro = 0.0, Ri = 0.0, XL = 0.0, x = 0.0, y = 0.0;
	double /*SLOP30 = 0.0, SLOP60 = 0.0,*/ Y1 = 0.0, Y2 = 0.0, Y3 = 0.0, Y4 = 0.0;
	double P1x = 0.0, P1y = 0.0, P2x = 0.0, P2y = 0.0, P3x = 0.0, P3y = 0.0, P4x = 0.0, P4y = 0.0;
	//double Tn;
	int in_quad = 0;
   //ZAperPlane: real;

	*ErrorFlag = 0;
	double SLOP30 = 0.57735026918962573; //tan(30.0*(acos(-1.0)/180.0));
	double SLOP60 = 1.7320508075688767; //tan(60.0*(acos(-1.0)/180.0));

	//AperturePlane(Element);           <------- calculated now in ODConcentrator
	//ZAperPlane = Element->ZAperture;

	//find intersection with surface first
	Intersect(PosRayIn, CosRayIn, Element, PosRayOut, CosRayOut, DFXYZ, PathLength, ErrorFlag);
	if (*ErrorFlag > 0 || *PathLength < 0)
	{
		*Intercept = 0;
		PosRayOut[0] = 0.0;
		PosRayOut[1] = 0.0;
		PosRayOut[2] = 0.0;
		CosRayOut[0] = 0.0;
		CosRayOut[1] = 0.0;
		CosRayOut[2] = 0.0;
		DFXYZ[0] = 0.0;
		DFXYZ[1] = 0.0;
		DFXYZ[2] = 0.0;
		*BacksideFlag = 0;
		*PathLength = 0.0;
		goto Label_100;
	}

	x = PosRayOut[0];
	y = PosRayOut[1];
	r = sqrt(x*x + y*y);



	switch (Element->ShapeIndex)
	{
	case 'c':
	case 'C': //circular aperture
			Ro = Element->ParameterA/2.0;

			if (r > Ro) //ray falls outsideside circular aperture
			{
			   *Intercept = 0;
			   PosRayOut[0] = 0.0;
			   PosRayOut[1] = 0.0;
			   PosRayOut[2] = 0.0;
			   CosRayOut[0] = 0.0;
			   CosRayOut[1] = 0.0;
			   CosRayOut[2] = 0.0;
			   DFXYZ[0] = 0.0;
			   DFXYZ[1] = 0.0;
			   DFXYZ[2] = 0.0;
			   *PathLength = 0.0;
			   *ErrorFlag = 0;
			   *BacksideFlag = 0;
			   goto Label_100;
			}
			else
			{
				if (DOT(CosRayIn, DFXYZ) < 0)
					*BacksideFlag = 0;
				else
					*BacksideFlag = 1;
				*Intercept = 1;
				goto Label_100;
			}
		break;
	
	case 'h':
	case 'H': //hexagonal aperture
			Ro = Element->ParameterA/2.0;

			if (r > Ro) //ray falls outside circular circumference aperture
			{
			   *Intercept = 0;
			   PosRayOut[0] = 0.0;
			   PosRayOut[1] = 0.0;
			   PosRayOut[2] = 0.0;
			   CosRayOut[0] = 0.0;
			   CosRayOut[1] = 0.0;
			   CosRayOut[2] = 0.0;
			   DFXYZ[0] = 0.0;
			   DFXYZ[1] = 0.0;
			   DFXYZ[2] = 0.0;
			   *PathLength = 0.0;
			   *ErrorFlag = 0;
			   *BacksideFlag = 0;
			   goto Label_100;
			}
			Ri = Ro*cos(30.0*(ACOSM1O180));

			if ( r <= Ri ) //ray falls inside inscribed circle
			{
				if ( DOT(CosRayIn, DFXYZ) < 0 )
					*BacksideFlag = 0;
				else
					*BacksideFlag = 1;
				*Intercept = 1;
				goto Label_100;
			}

			XL = sqrt(Ro*Ro - Ri*Ri); //otherwise break hexagon into 3 sections
			if ( (x <= Ro) && (x > XL) )  //1st section
			{
				Y1 = SLOP60*(x-Ro);
				Y2 = -Y1;
				if ( (y >= Y1) && (y <= Y2) )
				{
					if ( DOT(CosRayIn, DFXYZ) < 0 )
						*BacksideFlag = 0;
					else
						*BacksideFlag = 1;
						
					*Intercept = 1;
					goto Label_100;
				}
				*Intercept = 0;
				PosRayOut[0] = 0.0;
				PosRayOut[1] = 0.0;
				PosRayOut[2] = 0.0;
				CosRayOut[0] = 0.0;
				CosRayOut[1] = 0.0;
				CosRayOut[2] = 0.0;
				DFXYZ[0] = 0.0;
				DFXYZ[1] = 0.0;
				DFXYZ[2] = 0.0;
				*PathLength = 0.0;
				*BacksideFlag = 0;
				*ErrorFlag = 0;
				goto Label_100;
			}
			
			if ( (x <= XL) && (x >= -XL) )    //2nd section
			{
				if ( (y >= -Ri) && (y <= Ri) )
				{
					if ( DOT(CosRayIn, DFXYZ) < 0 )
						*BacksideFlag = 0;
					else
						*BacksideFlag = 1;
					*Intercept = 1;
					goto Label_100;
				}
				*Intercept = 0;
				PosRayOut[0] = 0.0;
				PosRayOut[1] = 0.0;
				PosRayOut[2] = 0.0;
				CosRayOut[0] = 0.0;
				CosRayOut[1] = 0.0;
				CosRayOut[2] = 0.0;
				DFXYZ[0] = 0.0;
				DFXYZ[1] = 0.0;
				DFXYZ[2] = 0.0;
				*PathLength = 0.0;
				*BacksideFlag = 0;
				*ErrorFlag = 0;
				goto Label_100;
			}
			
			if ( (x < -XL) && (x >= -Ro) )    //3rd section
			{
				Y3 = SLOP60*(x+Ro);
				Y4 = -Y3;
				if ( (y >= Y4) && (y <= Y3) )
				{
					if ( DOT(CosRayIn, DFXYZ) < 0 )
						*BacksideFlag = 0;
					else
						*BacksideFlag = 1;
					*Intercept = 1;
					goto Label_100;
				}
				*Intercept = 0;
				PosRayOut[0] = 0.0;
				PosRayOut[1] = 0.0;
				PosRayOut[2] = 0.0;
				CosRayOut[0] = 0.0;
				CosRayOut[1] = 0.0;
				CosRayOut[2] = 0.0;
				DFXYZ[0] = 0.0;
				DFXYZ[1] = 0.0;
				DFXYZ[2] = 0.0;
				*PathLength = 0.0;
				*BacksideFlag = 0;
				*ErrorFlag = 0;
				goto Label_100;
			}
		break;

	case 't':
	case 'T': //Triangular aperture
			Ro = Element->ParameterA/2.0;

			if ( r > Ro ) //ray falls outside circular circumference aperture
			{
				*Intercept = 0;
				PosRayOut[0] = 0.0;
				PosRayOut[1] = 0.0;
				PosRayOut[2] = 0.0;
				CosRayOut[0] = 0.0;
				CosRayOut[1] = 0.0;
				CosRayOut[2] = 0.0;
				DFXYZ[0] = 0.0;
				DFXYZ[1] = 0.0;
				DFXYZ[2] = 0.0;
				*PathLength = 0.0;
				*BacksideFlag = 0;
				*ErrorFlag = 0;
				goto Label_100;
			}
			Ri = Ro*sin(30.0*(ACOSM1O180));

			if ( r <= Ri )  //ray falls inside inscribed circle
			{
				if ( DOT(CosRayIn, DFXYZ) < 0 )
					*BacksideFlag = 0;
				else
					*BacksideFlag = 1;
				*Intercept = 1;
				goto Label_100;
			}

			if ( (x <= Ro) && (x > 0.0) )  //1st section
			{
				Y1 = -SLOP60*(x-Ri/cos(30.0*(ACOSM1O180)));
				Y2 = -Ri;
				if ( (y <= Y1) && (y >= Y2) )
				{
					if ( DOT(CosRayIn, DFXYZ) < 0 )
						*BacksideFlag = 0;
					else
						*BacksideFlag = 1;
					*Intercept = 1;
					goto Label_100;
				}
				*Intercept = 0;
				PosRayOut[0] = 0.0;
				PosRayOut[1] = 0.0;
				PosRayOut[2] = 0.0;
				CosRayOut[0] = 0.0;
				CosRayOut[1] = 0.0;
				CosRayOut[2] = 0.0;
				DFXYZ[0] = 0.0;
				DFXYZ[1] = 0.0;
				DFXYZ[2] = 0.0;
				*PathLength = 0.0;
				*BacksideFlag = 0;
				*ErrorFlag = 0;
				goto Label_100;
			}
			if ( (x >= -Ro) && (x <= 0.0) )  //2nd section
			{
				Y3 = SLOP60*(x+Ri/cos(30.0*(ACOSM1O180)));
				Y4 = -Ri;
				if ( (y >= Y4) && (y <= Y3) )
				{
					if ( DOT(CosRayIn, DFXYZ) < 0 )
						*BacksideFlag = 0;
					else
						*BacksideFlag = 1;
					*Intercept = 1;
					goto Label_100;
				}
				*Intercept = 0;
				PosRayOut[0] = 0.0;
				PosRayOut[1] = 0.0;
				PosRayOut[2] = 0.0;
				CosRayOut[0] = 0.0;
				CosRayOut[1] = 0.0;
				CosRayOut[2] = 0.0;
				DFXYZ[0] = 0.0;
				DFXYZ[1] = 0.0;
				DFXYZ[2] = 0.0;
				*PathLength = 0.0;
				*BacksideFlag = 0;
				*ErrorFlag = 0;
				goto Label_100;
			}
		break;

	case 'r':
	case 'R': //Rectangular aperture
               
			if ( (x > Element->ParameterA/2.0) || (x < -Element->ParameterA/2.0) )
			{
				*Intercept = 0;
				PosRayOut[0] = 0.0;
				PosRayOut[1] = 0.0;
				PosRayOut[2] = 0.0;
				CosRayOut[0] = 0.0;
				CosRayOut[1] = 0.0;
				CosRayOut[2] = 0.0;
				DFXYZ[0] = 0.0;
				DFXYZ[1] = 0.0;
				DFXYZ[2] = 0.0;
				*PathLength = 0.0;
				*BacksideFlag = 0;
				*ErrorFlag = 0;
				goto Label_100;
			}
			
			if ( (y > Element->ParameterB/2.0) || (y < -Element->ParameterB/2.0) )
			{
				*Intercept = 0;
				PosRayOut[0] = 0.0;
				PosRayOut[1] = 0.0;
				PosRayOut[2] = 0.0;
				CosRayOut[0] = 0.0;
				CosRayOut[1] = 0.0;
				CosRayOut[2] = 0.0;
				DFXYZ[0] = 0.0;
				DFXYZ[1] = 0.0;
				DFXYZ[2] = 0.0;
				*PathLength = 0.0;
				*BacksideFlag = 0;
				*ErrorFlag = 0;
				goto Label_100;
			}

			if ( DOT(CosRayIn, DFXYZ) < 0 ) //successfully falls on rectangle
				*BacksideFlag = 0;
			else
				*BacksideFlag = 1;
				
			*Intercept = 1;
			goto Label_100;

		break;

	case 'a':
	case 'A': //Annulus or torus contour
	
			if ( (Element->ParameterA == 0.0) && (Element->ParameterB == 0.0) ) goto Label_5; //torus

			if ( (r < Element->ParameterA) || (r > Element->ParameterB) )
			{
				*Intercept = 0;
				PosRayOut[0] = 0.0;
				PosRayOut[1] = 0.0;
				PosRayOut[2] = 0.0;
				CosRayOut[0] = 0.0;
				CosRayOut[1] = 0.0;
				CosRayOut[2] = 0.0;
				DFXYZ[0] = 0.0;
				DFXYZ[1] = 0.0;
				DFXYZ[2] = 0.0;
				*PathLength = 0.0;
				*BacksideFlag = 0;
				*ErrorFlag = 0;
				goto Label_100;
			}
Label_5:
			if ( x >= 0.0 )
			{
				if ( (asin(y/r) > Element->ParameterC*(ACOSM1O180)/2.0) || (asin(y/r) < -Element->ParameterC*(ACOSM1O180)/2.0) )
				{
					*Intercept = 0;
					PosRayOut[0] = 0.0;
					PosRayOut[1] = 0.0;
					PosRayOut[2] = 0.0;
					CosRayOut[0] = 0.0;
					CosRayOut[1] = 0.0;
					CosRayOut[2] = 0.0;
					DFXYZ[0] = 0.0;
					DFXYZ[1] = 0.0;
					DFXYZ[2] = 0.0;
					*PathLength = 0.0;
					*BacksideFlag = 0;
					*ErrorFlag = 0;
					goto Label_100;
				}

				if ( DOT(CosRayIn, DFXYZ) < 0 ) //successfully falls on annular section
					*BacksideFlag = 0;
				else
					*BacksideFlag = 1;
				*Intercept = 1;
				goto Label_100;
			}
			
			if ( x < 0.0 )
			{
				if ( (y >= 0) && ((acos(y/r)+M_PI/2.0) > Element->ParameterC*(ACOSM1O180)/2.0) )
				{
					*Intercept = 0;
					PosRayOut[0] = 0.0;
					PosRayOut[1] = 0.0;
					PosRayOut[2] = 0.0;
					CosRayOut[0] = 0.0;
					CosRayOut[1] = 0.0;
					CosRayOut[2] = 0.0;
					DFXYZ[0] = 0.0;
					DFXYZ[1] = 0.0;
					DFXYZ[2] = 0.0;
					*PathLength = 0.0;
					*BacksideFlag = 0;
					*ErrorFlag = 0;
					goto Label_100;
				}
				else if ( (y < 0) && ((-acos(-y/r)-M_PI/2.0) < -Element->ParameterC*(ACOSM1O180)/2.0) )
				{
					*Intercept = 0;
					PosRayOut[0] = 0.0;
					PosRayOut[1] = 0.0;
					PosRayOut[2] = 0.0;
					CosRayOut[0] = 0.0;
					CosRayOut[1] = 0.0;
					CosRayOut[2] = 0.0;
					DFXYZ[0] = 0.0;
					DFXYZ[1] = 0.0;
					DFXYZ[2] = 0.0;
					*PathLength = 0.0;
					*BacksideFlag = 0;
					*ErrorFlag = 0;
					goto Label_100;
				}

				if ( DOT(CosRayIn, DFXYZ) < 0 ) //successfully falls on annular section
					*BacksideFlag = 0;
				else
					*BacksideFlag = 1;
				*Intercept = 1;
				goto Label_100;
			}
		break;

	case 'l':
	case 'L': //off axis aperture section of line focus trough  or cylinder
			if ( (Element->ParameterA == 0.0) && (Element->ParameterB == 0.0) ) goto Label_10; //for cylinder, only need to check for limits on y

			if ( (x < Element->ParameterA) || (x > Element->ParameterB) )
			{
				*Intercept = 0;
				PosRayOut[0] = 0.0;
				PosRayOut[1] = 0.0;
				PosRayOut[2] = 0.0;
				CosRayOut[0] = 0.0;
				CosRayOut[1] = 0.0;
				CosRayOut[2] = 0.0;
				DFXYZ[0] = 0.0;
				DFXYZ[1] = 0.0;
				DFXYZ[2] = 0.0;
				*PathLength = 0.0;
				*BacksideFlag = 0;
				*ErrorFlag = 0;
				goto Label_100;
			}

Label_10:
			if ( (y < -Element->ParameterC/2.0) || (y > Element->ParameterC/2.0) )
			{
				*Intercept = 0;
				PosRayOut[0] = 0.0;
				PosRayOut[1] = 0.0;
				PosRayOut[2] = 0.0;
				CosRayOut[0] = 0.0;
				CosRayOut[1] = 0.0;
				CosRayOut[2] = 0.0;
				DFXYZ[0] = 0.0;
				DFXYZ[1] = 0.0;
				DFXYZ[2] = 0.0;
				*PathLength = 0.0;
				*BacksideFlag = 0;
				*ErrorFlag = 0;
				goto Label_100;
			}

			if ( DOT(CosRayIn, DFXYZ) < 0 ) //successfully falls on line focus or cylindrical section
				*BacksideFlag = 0;
			else
				*BacksideFlag = 1;
				
			*Intercept = 1;
			goto Label_100;
		break;

	case 'i':
	case 'I': //irregular triangle
			P1x = Element->ParameterA;
			P1y = Element->ParameterB;
			P2x = Element->ParameterC;
			P2y = Element->ParameterD;
			P3x = Element->ParameterE;
			P3y = Element->ParameterF;

			if (!intri( P1x, P1y, P2x, P2y, P3x, P3y, x, y ))
			{
				*Intercept = 0;
				PosRayOut[0] = 0.0;
				PosRayOut[1] = 0.0;
				PosRayOut[2] = 0.0;
				CosRayOut[0] = 0.0;
				CosRayOut[1] = 0.0;
				CosRayOut[2] = 0.0;
				DFXYZ[0] = 0.0;
				DFXYZ[1] = 0.0;
				DFXYZ[2] = 0.0;
				*PathLength = 0.0;
				*ErrorFlag = 0;
				*BacksideFlag = 0;
				goto Label_100;
			}

			if ( DOT(CosRayIn, DFXYZ) < 0 )
				*BacksideFlag = 0;
			else
				*BacksideFlag = 1;
			*Intercept = 1;
			goto Label_100;
		break;

	case 'q':
	case 'Q': //irregular quadrilateral
			P1x = Element->ParameterA;
			P1y = Element->ParameterB;
			P2x = Element->ParameterC;
			P2y = Element->ParameterD;
			P3x = Element->ParameterE;
			P3y = Element->ParameterF;
			P4x = Element->ParameterG;
			P4y = Element->ParameterH;

			in_quad = inquad(P1x, P1y, P2x, P2y, P3x, P3y, P4x, P4y, x, y);

			if (!in_quad)
			{
				*Intercept = 0;
				PosRayOut[0] = 0.0;
				PosRayOut[1] = 0.0;
				PosRayOut[2] = 0.0;
				CosRayOut[0] = 0.0;
				CosRayOut[1] = 0.0;
				CosRayOut[2] = 0.0;
				DFXYZ[0] = 0.0;
				DFXYZ[1] = 0.0;
				DFXYZ[2] = 0.0;
				*PathLength = 0.0;
				*ErrorFlag = 0;
				*BacksideFlag = 0;
				goto Label_100;
			}

			if ( DOT(CosRayIn, DFXYZ) < 0 )
				*BacksideFlag = 0;
			else
				*BacksideFlag = 1;
			*Intercept = 1;
			goto Label_100;
		break;
	} //end select case
	
Label_100:
	if ( *BacksideFlag )   //if hit on backside of element then slope of surface is reversed
	{
		DFXYZ[0] = -DFXYZ[0];
		DFXYZ[1] = -DFXYZ[1];
		DFXYZ[2] = -DFXYZ[2];
	}
}
//End of Procedure--------------------------------------------------------------

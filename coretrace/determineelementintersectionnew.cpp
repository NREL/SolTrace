
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
		*Intercept = false;
		PosRayOut[0] = 0.0;
		PosRayOut[1] = 0.0;
		PosRayOut[2] = 0.0;
		CosRayOut[0] = 0.0;
		CosRayOut[1] = 0.0;
		CosRayOut[2] = 0.0;
		DFXYZ[0] = 0.0;
		DFXYZ[1] = 0.0;
		DFXYZ[2] = 0.0;
		*BacksideFlag = false;
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
			   *Intercept = false;
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
			   *BacksideFlag = false;
			   goto Label_100;
			}
			else
			{
				if (DOT(CosRayIn, DFXYZ) < 0)
					*BacksideFlag = false;
				else
					*BacksideFlag = true;
				*Intercept = true;
				goto Label_100;
			}
		break;
	
	case 'h':
	case 'H': //hexagonal aperture
			Ro = Element->ParameterA/2.0;

			if (r > Ro) //ray falls outside circular circumference aperture
			{
			   Intercept = false;
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
			   *BacksideFlag = false;
			   goto Label_100;
			}
			Ri = Ro*cos(30.0*(ACOSM1O180));

			if ( r <= Ri ) //ray falls inside inscribed circle
			{
				if ( DOT(CosRayIn, DFXYZ) < 0 )
					*BacksideFlag = false;
				else
					*BacksideFlag = true;
				*Intercept = true;
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
						*BacksideFlag = false;
					else
						*BacksideFlag = true;
						
					*Intercept = true;
					goto Label_100;
				}
				*Intercept = false;
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
				*BacksideFlag = false;
				*ErrorFlag = 0;
				goto Label_100;
			}
			
			if ( (x <= XL) && (x >= -XL) )    //2nd section
			{
				if ( (y >= -Ri) && (y <= Ri) )
				{
					if ( DOT(CosRayIn, DFXYZ) < 0 )
						*BacksideFlag = false;
					else
						*BacksideFlag = true;
					*Intercept = true;
					goto Label_100;
				}
				*Intercept = false;
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
				*BacksideFlag = false;
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
						*BacksideFlag = false;
					else
						*BacksideFlag = true;
					*Intercept = true;
					goto Label_100;
				}
				*Intercept = false;
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
				*BacksideFlag = false;
				*ErrorFlag = 0;
				goto Label_100;
			}
		break;

	case 't':
	case 'T': //Triangular aperture
			Ro = Element->ParameterA/2.0;

			if ( r > Ro ) //ray falls outside circular circumference aperture
			{
				*Intercept = false;
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
				*BacksideFlag = false;
				*ErrorFlag = 0;
				goto Label_100;
			}
			Ri = Ro*sin(30.0*(ACOSM1O180));

			if ( r <= Ri )  //ray falls inside inscribed circle
			{
				if ( DOT(CosRayIn, DFXYZ) < 0 )
					*BacksideFlag = false;
				else
					*BacksideFlag = true;
				*Intercept = true;
				goto Label_100;
			}

			if ( (x <= Ro) && (x > 0.0) )  //1st section
			{
				Y1 = -SLOP60*(x-Ri/cos(30.0*(ACOSM1O180)));
				Y2 = -Ri;
				if ( (y <= Y1) && (y >= Y2) )
				{
					if ( DOT(CosRayIn, DFXYZ) < 0 )
						*BacksideFlag = false;
					else
						*BacksideFlag = true;
					*Intercept = true;
					goto Label_100;
				}
				*Intercept = false;
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
				*BacksideFlag = false;
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
						*BacksideFlag = false;
					else
						*BacksideFlag = true;
					*Intercept = true;
					goto Label_100;
				}
				*Intercept = false;
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
				*BacksideFlag = false;
				*ErrorFlag = 0;
				goto Label_100;
			}
		break;

	case 'r':
	case 'R': //Rectangular aperture
               
			if ( (x > Element->ParameterA/2.0) || (x < -Element->ParameterA/2.0) )
			{
				*Intercept = false;
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
				*BacksideFlag = false;
				*ErrorFlag = 0;
				goto Label_100;
			}
			
			if ( (y > Element->ParameterB/2.0) || (y < -Element->ParameterB/2.0) )
			{
				*Intercept = false;
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
				*BacksideFlag = false;
				*ErrorFlag = 0;
				goto Label_100;
			}

			if ( DOT(CosRayIn, DFXYZ) < 0 ) //successfully falls on rectangle
				*BacksideFlag = false;
			else
				*BacksideFlag = true;
				
			*Intercept = true;
			goto Label_100;

		break;

	case 'a':
	case 'A': //Annulus or torus contour
	
			if ( (Element->ParameterA == 0.0) && (Element->ParameterB == 0.0) ) goto Label_5; //torus

			if ( (r < Element->ParameterA) || (r > Element->ParameterB) )
			{
				*Intercept = false;
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
				*BacksideFlag = false;
				*ErrorFlag = 0;
				goto Label_100;
			}
Label_5:
			if ( x >= 0.0 )
			{
				if ( (asin(y/r) > Element->ParameterC*(ACOSM1O180)/2.0) || (asin(y/r) < -Element->ParameterC*(ACOSM1O180)/2.0) )
				{
					*Intercept = false;
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
					*BacksideFlag = false;
					*ErrorFlag = 0;
					goto Label_100;
				}

				if ( DOT(CosRayIn, DFXYZ) < 0 ) //successfully falls on annular section
					*BacksideFlag = false;
				else
					*BacksideFlag = true;
				*Intercept = true;
				goto Label_100;
			}
			
			if ( x < 0.0 )
			{
				if ( (y >= 0) && ((acos(y/r)+M_PI/2.0) > Element->ParameterC*(ACOSM1O180)/2.0) )
				{
					*Intercept = false;
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
					*BacksideFlag = false;
					*ErrorFlag = 0;
					goto Label_100;
				}
				else if ( (y < 0) && ((-acos(-y/r)-M_PI/2.0) < -Element->ParameterC*(ACOSM1O180)/2.0) )
				{
					*Intercept = false;
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
					*BacksideFlag = false;
					*ErrorFlag = 0;
					goto Label_100;
				}

				if ( DOT(CosRayIn, DFXYZ) < 0 ) //successfully falls on annular section
					*BacksideFlag = false;
				else
					*BacksideFlag = true;
				*Intercept = true;
				goto Label_100;
			}
		break;

	case 'l':
	case 'L': //off axis aperture section of line focus trough  or cylinder
			if ( (Element->ParameterA == 0.0) && (Element->ParameterB == 0.0) ) goto Label_10; //for cylinder, only need to check for limits on y

			if ( (x < Element->ParameterA) || (x > Element->ParameterB) )
			{
				*Intercept = false;
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
				*BacksideFlag = false;
				*ErrorFlag = 0;
				goto Label_100;
			}

Label_10:
			if ( (y < -Element->ParameterC/2.0) || (y > Element->ParameterC/2.0) )
			{
				*Intercept = false;
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
				*BacksideFlag = false;
				*ErrorFlag = 0;
				goto Label_100;
			}

			if ( DOT(CosRayIn, DFXYZ) < 0 ) //successfully falls on line focus or cylindrical section
				*BacksideFlag = false;
			else
				*BacksideFlag = true;
				
			*Intercept = true;
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
				*Intercept = false;
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
				*BacksideFlag = false;
				goto Label_100;
			}

			if ( DOT(CosRayIn, DFXYZ) < 0 )
				*BacksideFlag = false;
			else
				*BacksideFlag = true;
			*Intercept = true;
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
				*Intercept = false;
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
				*BacksideFlag = false;
				goto Label_100;
			}

			if ( DOT(CosRayIn, DFXYZ) < 0 )
				*BacksideFlag = false;
			else
				*BacksideFlag = true;
			*Intercept = true;
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

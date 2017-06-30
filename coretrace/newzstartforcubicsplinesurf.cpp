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

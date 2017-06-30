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

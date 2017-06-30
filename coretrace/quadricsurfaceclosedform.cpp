#include <math.h>

#include "types.h"
#include "procs.h"

#define sqr(x) (x*x)

void QuadricSurfaceClosedForm(
			TElement *Element,
			double PosLoc[3],
			double CosLoc[3],
			double PosXYZ[3],
			double DFXYZ[3],
			double *PathLength,
			int *ErrorFlag)
{
	double Xdelta = 0.0,Ydelta = 0.0,Zdelta = 0.0;
	double Xc=0,Yc=0,Zc=0,Kx=0,Ky=0,Kz=0;
	double r = 0.0,r2 = 0.0,a2=0,b2=0,c2=0;
	double t1 = 0.0,t2 = 0.0,A=0,B=0,C=0,slopemag = 0.0;

	*ErrorFlag = 0;

	switch( Element->SurfaceIndex )
	{
	case 's':
	case 'S': // sphere
			a2 = 1;
			b2 = 1;
			c2 = 1;
			Kx = 1;
			Ky = 1;
			Kz = 1;
			r = 1.0/Element->VertexCurvX;
			r2 = r*r;
			Xc = 0.0;
			Yc = 0.0;
			Zc = r;
			
			Xdelta = PosLoc[0] - Xc;
			Ydelta = PosLoc[1] - Yc;
			Zdelta = PosLoc[2] - Zc;
			
			A = CosLoc[0]*CosLoc[0]*Kx/a2 + CosLoc[1]*CosLoc[1]*Ky/b2 + CosLoc[2]*CosLoc[2]*Kz/c2;
			B = 2.0*(Kx*Xdelta*CosLoc[0]/a2 + Ky*Ydelta*CosLoc[1]/b2 + Kz*Zdelta*CosLoc[2]/c2);
			C = Kx*Xdelta*Xdelta/a2 + Ky*Ydelta*Ydelta/a2 + Kz*Zdelta*Zdelta/c2 - r2;
		break;

	case 'p':
	case 'P':   //parabola
			a2 = 4.0*(1.0/Element->VertexCurvX)/2.0;
			b2 = a2;
			c2 = 1.0;
			Xc = 0.0;
			Yc = 0.0;
			Zc = 0.0;
			
			Xdelta = PosLoc[0] - Xc;
			Ydelta = PosLoc[1] - Yc;
			Zdelta = PosLoc[2] - Zc;
			
			A = sqr(CosLoc[0])/a2 + sqr(CosLoc[1])/b2;
			B = 2.0*CosLoc[0]*Xdelta/a2 + 2.0*CosLoc[1]*Ydelta/b2 - CosLoc[2]/c2;
			C = sqr(Xdelta)/a2 + sqr(Ydelta)/b2 - Zdelta/c2;
		break;
		
	case 'o':
	case 'O':   //other
		break;

	case 't':
	case 'T':   //cylinder
			a2 = 1;
			b2 = 1;
			c2 = 1;
			Kx = 1;
			Ky = 0;
			Kz = 1;
			r = 1.0/Element->CurvOfRev;
			r2 = r*r;
			Xc = 0.0;
			Yc = 0.0;
			Zc = r;
			
			Xdelta = PosLoc[0] - Xc;
			Ydelta = PosLoc[1] - Yc;
			Zdelta = PosLoc[2] - Zc;
			
			A = CosLoc[0]*CosLoc[0]*Kx/a2 + CosLoc[1]*CosLoc[1]*Ky/b2 + CosLoc[2]*CosLoc[2]*Kz/c2;
			B = 2.0*(Kx*Xdelta*CosLoc[0]/a2 + Ky*Ydelta*CosLoc[1]/b2 + Kz*Zdelta*CosLoc[2]/c2);
			C = Kx*Xdelta*Xdelta/a2 + Ky*Ydelta*Ydelta/a2 + Kz*Zdelta*Zdelta/c2 - r2;
		break;

	case 'c':
	case 'C':   //cone
		break;

	case 'f':
	case 'F':   //flat
		break;
	}

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

			//*************************************************************************************************************
			//makes sure to get shortest ray path on valid side of surface; 10-05-10    for open surface of parabola
			//if cylinder, then PosXYZ[3] will always be less than or equal to Element.Zaperture so never passes this test.
			// Test for  cylinder follows below.
			if (PosXYZ[2] > Element->ZAperture)
			{
				PosXYZ[0] = PosLoc[0] + t1*CosLoc[0];
				PosXYZ[1] = PosLoc[1] + t1*CosLoc[1];
				PosXYZ[2] = PosLoc[2] + t1*CosLoc[2];
				*PathLength = t1;
			}

			// Remember at this point, intersection is being found on an INFINITELY long cylinder.
			//if 1st intersection on INFINITELY long cylinder is from the outside, t2, check to make sure
			//intersection is within the finite
			// length of the actual cylinder geometry, if not then 2nd intersection on the inside, t1,
			//is valid one to use.  This means ray could
			// enter from the open end  of the cylinder and hit on the inside.  The final test for this is in the
			//calling routine:  DetermineElementIntersectionNew
			// Wendelin 10-05-10
			if ((Element->SurfaceIndex == 't') || (Element->SurfaceIndex == 'T'))
			{
				if ((PosXYZ[1] < -Element->ParameterC/2.0) || (PosXYZ[1] > Element->ParameterC/2.0))
				{
					PosXYZ[0] = PosLoc[0] + t1*CosLoc[0];
					PosXYZ[1] = PosLoc[1] + t1*CosLoc[1];
					PosXYZ[2] = PosLoc[2] + t1*CosLoc[2];
					*PathLength = t1;
				}
			}
           //***********************************************************************************************************

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
	slopemag = sqrt(sqr(2.0*Kx*(PosXYZ[0] - Xc)/a2)+sqr(2.0*Ky*(PosXYZ[1] - Yc)/b2)+sqr(2.0*Kz*(PosXYZ[2] - Zc)/c2));
	DFXYZ[0] = -(2.0*Kx*(PosXYZ[0] - Xc)/a2)/slopemag;
	DFXYZ[1] = -(2.0*Ky*(PosXYZ[1] - Yc)/b2)/slopemag;
	DFXYZ[2] = -(2.0*Kz*(PosXYZ[2] - Zc)/c2)/slopemag;
}
//end of procedure--------------------------------------------------------------

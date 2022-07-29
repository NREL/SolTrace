
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


#include <stdlib.h>
#include <math.h>

#include "types.h"
#include "procs.h"

void Surface(
			double PosXYZ[3],
			TElement *Element,
			double *FXYZ,
			double DFXYZ[3],
			int *ErrorFlag )
{
/*{Purpose: To compute the surface equation and it's derivatives for various
geometric surfaces.
    Input - PosXYZ[3] = X, Y, Z coordinate position
            Element.SurfaceType = Surface type flag
                          = 1 for rotationally symmetric surfaces
                          = 2 for torics and cylinders
                          = 3 for plane surfaces
                          = 4 for surface interpolated from finite element data points
                          = 5 for surface interpolated from VSHOT data points
                          = 6 for surface described by Zernike monomials
                          = 7 for single axis parabolic curvature surfaces
                          = 8 for rotationally symmetric polynomial description
                          = 9 for       "          "     cubic spline interpolation
            Element.Alpha = Sensitivity coefficients which specify deviation from conic
                    of revolution
            Element.VertexCurvX = Vertex Curvature of surface
            Element.Kappa = Surface specifier
                 < 0         ==> Hyperboloid
                 = 0         ==> Paraboloid
                 > 0 and < 1 ==> Hemelipsoid of revolution about major axis
                 = 1         ==> Hemisphere
                 > 1         ==> Hemelipsoid of revolution about minor axis
            Element.ConeHalfAngle = Half-angle of cone for cones or revolution or axicons
            Element.CurvOfRev = Curvature of revolution

    Output - FXYZ = Surface equation
             DFXYZ[3] = Derivatives of surface equation
             ErrorFlag = Error Flag
                         = 0  ==> no errors
                         > 0  ==> interpolation error
}*/

	int i=0,k=0;
	double X=0.0,Y=0.0,Z=0.0;
	double Rho=0.0, Rho2=0.0, Rho2i=0.0;
	double DFDX=0, DFDY=0, DFDZ=0;
	double Sum1=0.0, Sum2=0.0, Term=0.0;
	double Y2=0.0, Y2J=0.0;
	double FY=0.0, ZZ=0.0;
	double zm=0.0, zr=0.0;
	double dzrdx=0.0, dzrdy=0.0, delzx=0.0, delzy=0.0;
	double density=0.0, delta=0.0;
	double zx=0.0, zy=0.0;
	double dzdRho=0.0, dRhodx=0.0, dRhody=0.0;

//label 100, 130,160, 850, 990;

	//Initialize variables
	X = PosXYZ[0];
	Y = PosXYZ[1];
	Z = PosXYZ[2];
	*ErrorFlag = 0;
     
	//===SurfaceType = 1, 7  Rotationally Symmetric surfaces and single axis curvature sections===========================
	if (Element->SurfaceType == 1 || Element->SurfaceType == 7)
	{
		if (Element->SurfaceType == 1)
			Rho2 = X*X + Y*Y;    //rotationally symmetric
		else
			Rho2 = X*X;         //single axis curvature depends only on x

		if (Element->ConeHalfAngle != 0.0) 
			goto Label_160;
		
		//wendelin 5-18-11 changes to allow different vertex curvature in the x and y directions for the parabola; this block of code
		//is a subset of the more general form below therefore it has been commented out.  It also assumes VertexCurvY = either VertexCurvX or zero
		//and doesn't allow different nonzero values for the parabolic case.   Not using the alpha parameters for the general case for now.
/*
		for (i=0;i<5;i++)
			if (Element->Alpha[i] != 0.0)
				goto Label_130;

		*FXYZ = Z - 0.5*Element->VertexCurvX*(Rho2 + Element->Kappa*Z*Z);
		DFDX = -Element->VertexCurvX*X;
		DFDY = -Element->VertexCurvY*Y; //VertexCurvY = VertexCurvX if rotationally symmetric or 0 if single axis curved
		DFDZ = 1.0 - Element->Kappa*Element->VertexCurvX*Z;
		goto Label_990;
*/
		Sum1 = 0.0;
		Sum2 = 0.0;
		Rho2i = 1.0;

		// wendelin 5-18-11
		/*
		for (i=0;i<5;i++)
		{
			Sum1 = i*Element->Alpha[i]*Rho2i + Sum1;
			Rho2i = Rho2i*Rho2;
			Sum2 = Element->Alpha[i]*Rho2i + Sum2;
		}*/

		//wendelin 5-18-11 changes to allow different vertex curvature in the x and y directions for the parabola only
		// Term = sqrt(1.0 - Element->Kappa*Element->VertexCurvX*Element->VertexCurvX*Rho2);
		Term = sqrt(1.0 - Element->Kappa*(Element->VertexCurvX*Element->VertexCurvX*X*X+Element->VertexCurvY*Element->VertexCurvY*Y*Y));   //new
		//*FXYZ = Z - Element->VertexCurvX*Rho2/(1.0 + Term) - Sum2;
		*FXYZ = Z - (Element->VertexCurvX*X*X+Element->VertexCurvY*Y*Y)/(1.0 + Term) - Sum2;   //new

		DFDX = -X*(Element->VertexCurvX/Term + 2.0*Sum1);
		DFDY = -Y*(Element->VertexCurvY/Term + 2.0*Sum1); //VertexCurvY = VertexCurvX if rotationally symmetric or 0 if single axis curved
		DFDZ = 1.0;
		goto Label_990;

Label_160:    
		*FXYZ = Z - sqrt(Rho2)/tan(Element->ConeHalfAngle*(ACOSM1O180));
		DFDX = -X/(sqrt(Rho2)*tan(Element->ConeHalfAngle*(ACOSM1O180)));
		DFDY = -Y/(sqrt(Rho2)*tan(Element->ConeHalfAngle*(ACOSM1O180)));
		DFDZ = 1.0;
		goto Label_990;
	}

//===SurfaceType = 2, Toroidal or Cylindrical surfaces========================== //not currently used
	if (Element->SurfaceType == 2)
	{
		Sum1 = 0.0;
		Sum2 = 0.0;
		Y2 = Y*Y;
		Y2J = 1.0;
		
		for (i=0;i<5;i++)
		{
			Sum1 = i*Element->Alpha[i]*Y2J*Y + Sum1;
			Y2J = Y2J*Y2;
			Sum2 = Element->Alpha[i]*Y2J + Sum2;
		}
		
		Term = sqrt(1.0 - Element->Kappa*Element->VertexCurvX*Element->VertexCurvX*Y2);
		FY = Element->VertexCurvX*Y2/(1.0 + Term) + Sum2;
		*FXYZ = Z - FY - 0.5*Element->CurvOfRev*(X*X + Z*Z - FY*FY);
		DFDX = -Element->CurvOfRev*X;
		DFDY = (Element->CurvOfRev*FY - 1.0)*(Element->VertexCurvX*Y/Term + 2.0*Sum1);
		DFDZ = 1.0 - Element->CurvOfRev*Z;
		goto Label_990;
	}

//===SurfaceType = 3, Plane Surfaces============================================
/*     {The equation of a plane is: kx + ly + mz = p,  where k,l,m are the direction
     cosines of the normal to the plane, and p is the distance from the origin
     to the plane.  In this case, these parameters are contained in the Alpha array.}*/
	if (Element->SurfaceType == 3)
	{
		DFDX = Element->Alpha[0];
		DFDY = Element->Alpha[1];
		DFDZ = Element->Alpha[2];
		*FXYZ = DFDX*X + DFDY*Y + DFDZ*Z - Element->Alpha[3];
		goto Label_990;
	}

//===SurfaceType = 4, Surface specified by finite element data==================
	if (Element->SurfaceType == 4)
	{
		Rho2 = X*X + Y*Y;
		if (Rho2 == 0.0)
		{
			//FXYZ := Z - ZA[1];  ZA not defined yet
			*FXYZ = Z;
			DFDX = 0.0;
			DFDY = 0.0;
			DFDZ = 1;
			goto Label_990;
		}
		
		//Interpolate to find the z
		density = Element->FEData.nodes.size()/Element->ApertureArea;
		delta = 0.001/sqrt(density);
		FEInterpKD(X, Y, &Element->FEData, delta, &zr, &dzrdx, &dzrdy);
		
		PosXYZ[2] = zr;
		*FXYZ = Z - zr;
		DFDX = dzrdx;
		DFDY = dzrdy;
		DFDZ = 1.0;
		goto Label_990;
	}

//===SurfaceType = 5, VSHOT data================================================
	if (Element->SurfaceType == 5)
	{
		Rho2 = X*X + Y*Y;
		if (Rho2 == 0.0)
		{
			*FXYZ = Z;
			DFDX = 0.0;
			DFDY = 0.0;
			DFDZ = 1.0;
			goto Label_990;
		}
		// evaluate z, dz/dx and dz/dy from the monomial fit at x,y
		EvalMono(X, Y, Element->BCoefficients, Element->FitOrder, 0.0, 0.0, &zm); //the 0.0's are values for DeltaX and DeltaY; **[need to look at this further]**

		//Interpolate to find the slope residuals
		density = Element->VSHOTData.nrows()/Element->ApertureArea;
		
		/*
		if (Element->ShapeIndex == 'l' || Element->ShapeIndex == 'L')       //interpolation scheme for single axis curvature surfaces
			VSHOTInterpolateNew(X, Y, density, Element->VSHOTData, Element->VSHOTData.nrows(), &delzx, &delzy);
		else
			VSHOTInterpolate(X, Y, density, Element->VSHOTData, Element->VSHOTData.nrows(), &delzx, &delzy);
		*/

		::VSHOTInterpolateModShepard(X, Y, density, Element->VSHOTData, Element->VSHOTData.nrows(), &delzx, &delzy, ErrorFlag);

		if ( *ErrorFlag != 0 ) return;


		//Evaluate "real" z (i.e. the best estimate for z comes from the monomial fit)
		zr = zm;
		
		//Now evaluate the slopes  -  what we want here is the measured slope which is the best value to use
		//dzrdx := dzmdx + delzx;     //fit slope + (fit slope - meas. slope) =  wrong value
		//dzrdy := dzmdy + delzy;
		
		//dzrdx := dzmdx - delzx;       //fit slope - (fit slope - meas. slope) = meas. slope  (this is what we want)
		//dzrdy := dzmdy - delzy;
		dzrdx = delzx;                //if VSHOTInterpolate returns interpolated measured slopes and not slope RESIDUALS
		dzrdy = delzy;                // These values are angles of the slope in radians. Need to convert to dimensionless dz/dy and dz/dx so take tangent of angle
		
		dzrdx = tan(dzrdx);
		dzrdy = tan(dzrdy);
		
		PosXYZ[2] = zr;
		*FXYZ = Z - zr;
		DFDX = dzrdx;
		DFDY = dzrdy;
		//change sign of derivatives to agree with SurfaceType = 1
		DFDX = -DFDX;
		DFDY = -DFDY;
		DFDZ = 1.0;
		goto Label_990;
	}

//===SurfaceType = 6, Zernike monomials=========================================
	if (Element->SurfaceType == 6)
	{
		k = Element->FitOrder + 1;
		ZZ = 0.0;
		DFDX = 0.0;
		DFDY = 0.0;
		
		// evaluate z from the monomial expression at x,y
		EvalMono(X, Y, Element->BCoefficients, Element->FitOrder, 0.0, 0.0, &ZZ); //the 0.0's are values for DeltaX and DeltaY; **[need to look at this further]**
		MonoSlope(Element->BCoefficients, Element->FitOrder, X, Y, &DFDX, &DFDY);
		
		PosXYZ[2] = ZZ;
		*FXYZ = Z - ZZ;
		//{change sign of derivatives to agree with SurfaceType = 1}
		DFDX = -DFDX;
		DFDY = -DFDY;
		DFDZ = 1.0;
		goto Label_990;
	}

//===SurfaceType = 8, rotationally symmetric polynomial surface=============================
	if (Element->SurfaceType == 8)
	{
		ZZ = 0.0;
		DFDX = 0.0;
		DFDY = 0.0;
		
		double yval = Y;
		if ( Element->ShapeIndex == 'l' || Element->ShapeIndex == 'L' )
			yval = 0.0;

		// evaluate z & slopes from the polynomial expression at r = sqrt(x^2+y^2)
		EvalPoly(X, yval, Element->PolyCoeffs, Element->FitOrder, &ZZ);
		PolySlope(Element->PolyCoeffs, Element->FitOrder, X, yval, &DFDX, &DFDY);
		
		PosXYZ[2] = ZZ;
		*FXYZ = Z - ZZ;
		//{change sign of derivatives to agree with SurfaceType = 1}
		DFDX = -DFDX;
		DFDY = -DFDY;
		DFDZ = 1.0;
		goto Label_990;
	}
     
//===SurfaceType = 9, rotationally symmetric cubic spline interpolation surface==============
	if (Element->SurfaceType == 9)
	{
		ZZ = 0.0;
		DFDX = 0.0;
		DFDY = 0.0;
		
		Rho = sqrt(X*X+Y*Y);
		dRhodx = X/Rho;
		dRhody = Y/Rho;
		
		if (Element->ShapeIndex == 'l'|| Element->ShapeIndex == 'L')  //x dimension only for single axis curvature
		{
			Rho = X;
			dRhodx = 1.0;
			dRhody = 0.0;
		}
		
		//evaluate z & slopes using cubic spline interpolation
		if (!splint(Element->CubicSplineXData,
				Element->CubicSplineYData,
				Element->CubicSplineY2Data,
				Element->CubicSplineXData.size(),
				Rho,&ZZ,&dzdRho))
		{
			*ErrorFlag = 3;
			return;
		}
						  
		DFDX = dzdRho*dRhodx;
		DFDY = dzdRho*dRhody;
		
		PosXYZ[2] = ZZ;
		*FXYZ = Z - ZZ;
		//{change sign of derivatives to agree with SurfaceType = 1}
		DFDX = -DFDX;
		DFDY = -DFDY;
		DFDZ = 1.0;
		goto Label_990;
	}

 //the following surfacetype is now handled above in the general case

//===SurfaceType = 7, single axis curvature parabolic or spherical surface=============================
     /*{if SurfaceType = 7 then
     begin
       if (SurfaceIndex = 'p') or (SurfaceIndex = 'P') then
       begin
          FXYZ := Z - X*X*VertexCurvX/2.0;
          DFDX := -X*VertexCurvX;
          DFDY := 0.0;
          DFDZ := 1.0;
       end;
       if (SurfaceIndex = 's') or (SurfaceIndex = 'S') then
       begin
        FXYZ := Z - 0.5*VertexCurvX*(X*X + Z*Z);
        DFDX := -VertexCurvX*X;
        DFDY := 0.0;
        DFDZ := 1.0 - VertexCurvX*Z;
       end;
     end;}*/

//{****Offload derivatives of the surface equation}
Label_990:
	DFXYZ[0] = DFDX;
	DFXYZ[1] = DFDY;
	DFXYZ[2] = DFDZ;

}
//end of procedure--------------------------------------------------------------

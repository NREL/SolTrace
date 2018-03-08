
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

void SurfaceZatXYPair(
			double PosXYZ[3],
			TElement *Element,
			double *FXYZ,
			int *ErrorFlag )
{
/*{Purpose: To compute the Z value of the surface equation at an X,Y pair.
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

    Output - FXYZ = Z value of Surface equation
             ErrorFlag = Error Flag
                         = 0  ==> no errors
                         > 0  ==> interpolation error
}*/
	int i=0;
	double X=0.0,Y=0.0,Z=0.0;
	double Rho2=0.0, Rho=0.0;
	double Sum1=0.0, ZZ=0.0,  zm=0.0;

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

        Rho = sqrt(Rho2);

        if (Element->ConeHalfAngle != 0.0) goto Label_160;

		//wendelin 5-18-11

		//if (Element->Kappa*Element->VertexCurvX*Element->VertexCurvX*Rho2 > 1.0)  //xy pair cannot be found on closed surface   06-10-07
		if ( Element->Kappa*(Element->VertexCurvX*Element->VertexCurvX*X*X+Element->VertexCurvY*Element->VertexCurvY*Y*Y) > 1.0 )  //xy pair cannot be found on closed surface   06-10-07
        {
			*FXYZ = 0.0;
			return;
		}

		//wendelin 5-18-11
		// *FXYZ = Element->VertexCurvX*Rho2/(1.0+sqrt(1.0-Element->Kappa*Element->VertexCurvX*Element->VertexCurvX*Rho2));
		*FXYZ = (Element->VertexCurvX*X*X+Element->VertexCurvY*Y*Y)
				/ (1.0+sqrt(1.0-Element->Kappa*(Element->VertexCurvX*Element->VertexCurvX*X*X+Element->VertexCurvY*Y*Y)));

/*        for (i=0;i<5;i++)
             if (Element->Alpha[i] != 0.0) goto Label_130;
			 */
             
        return;

		Sum1 = 0.0;
		for (i=0;i<5;i++)
             Sum1 = Element->Alpha[i]*pow(Rho,2*(i+1)) + Sum1;

        *FXYZ += Sum1;
        return;

Label_160:
		*FXYZ = sqrt(Rho2)/tan(Element->ConeHalfAngle*(ACOSM1O180));
        return;
	}

//===SurfaceType = 3, Plane Surfaces============================================
    /* {The equation of a plane is: kx + ly + mz = p,  where k,l,m are the direction
     cosines of the normal to the plane, and p is the distance from the origin
     to the plane.  In this case, these parameters are contained in the Alpha array.}
     {if SurfaceType = 3 then
     begin
        DFDX = Alpha[1];
        DFDY = Alpha[2];
        DFDZ = Alpha[3];
        FXYZ = DFDX*X + DFDY*Y + DFDZ*Z - Alpha[4];
        return;
     end;}*/

//===SurfaceType = 4, Surface specified by finite element data==================
    /*{if SurfaceType = 4 then
     begin
        Rho2 = X*X + Y*Y;
        if Rho2 = 0.0 then
        begin
             //FXYZ = Z - ZA[1];  ZA not defined yet
             FXYZ = Z;
             DFDX = 0.0;
             DFDY = 0.0;
             DFDZ = 1;
             return;
        end;
          //Interpolate to find the z
          Density = FENumPoints/ApertureArea;
          delta = 0.1/sqrt(density);
          FEInterpNew(X, Y, Density, FEData, FENumPoints, zr);

          //Now evaluate the slopes
          FEInterpNew(X+delta, Y, Density, FEData, FENumPoints, zx);
          FEInterpNew(X, Y+delta, Density, FEData, FENumPoints, zy);
          dzrdx = (zx-zr)/delta;
          dzrdy = (zy-zr)/delta;

          PosXYZ[3] = zr;
          FXYZ = z - zr;
          DFDX = dzrdx;
          DFDY = dzrdy;
          //change sign of derivatives to agree with SurfaceType = 1
          DFDX = -DFDX;
          DFDY = -DFDY;
          DFDZ = 1.0;
          return;
     end;}*/

//===SurfaceType = 5, VSHOT data================================================
	if (Element->SurfaceType == 5)
	{
		Rho2 = X*X + Y*Y;
		if (Rho2 == 0.0)
		{
			*FXYZ = 0.0;
			return;
		}
		// evaluate z, dz/dx and dz/dy from the monomial fit at x,y
		EvalMono(X, Y, Element->BCoefficients, Element->FitOrder, 0.0, 0.0, &zm); //the 0.0's are values for DeltaX and DeltaY; **[need to look at this further]**
		*FXYZ = zm;
		return;
	}

//===SurfaceType = 6, Zernike monomials=========================================
	if (Element->SurfaceType == 6)
	{
          // evaluate z from the monomial expression at x,y
		EvalMono(X, Y, Element->BCoefficients, Element->FitOrder, 0.0, 0.0, &ZZ); //the 0.0's are values for DeltaX and DeltaY; **[need to look at this further]**
		*FXYZ = ZZ;
		return;
	}

//===SurfaceType = 8, rotationally symmetric polynomial surface=============================
	if (Element->SurfaceType == 8)
	{
		// evaluate z & slopes from the polynomial expression at r = sqrt(x^2+y^2)

		double yval = Y;
		if ( Element->ShapeIndex == 'l' || Element->ShapeIndex == 'L' )
			yval = 0.0;

		EvalPoly(X, yval, Element->PolyCoeffs, Element->FitOrder, &ZZ);
		*FXYZ = ZZ;
		return;
	}
//===SurfaceType = 9, rotationally symmetric cubic spline interpolation surface==============
     /*if (Element->SurfaceType == 9)
	 {
		ZZ = 0.0;
		DFDX = 0.0;
		DFDY = 0.0;

		Rho = sqrt(X*X+Y*Y);
		dRhodx = X/Rho;
		dRhody = Y/Rho;
		//evaluate z & slopes using cubic spline interpolation
		splint(Element->CubicSplineXData,
			Element->CubicSplineYData,
			Element->CubicSplineY2Data,
			Element->CubicSplineXData.length(),
			Rho,
			&ZZ,&dzdRho);

		DFDX = dzdRho*dRhodx;
		DFDY = dzdRho*dRhody;

		PosXYZ[2] = ZZ;
		*FXYZ = Z - ZZ;
		//change sign of derivatives to agree with SurfaceType = 1
		DFDX = -DFDX;
		DFDY = -DFDY;
		return;
	 }*/

 //the following surfacetype is now handled above in the general case

//===SurfaceType = 7, single axis curvature parabolic or spherical surface=============================
     /*{if SurfaceType = 7 then
     begin
       if (SurfaceIndex = 'p') or (SurfaceIndex = 'P') then
       begin
          FXYZ = Z - X*X*VertexCurvX/2.0;
          DFDX = -X*VertexCurvX;
          DFDY = 0.0;
          DFDZ = 1.0;
       end;
       if (SurfaceIndex = 's') or (SurfaceIndex = 'S') then
       begin
        FXYZ = Z - 0.5*VertexCurvX*(X*X + Z*Z);
        DFDX = -VertexCurvX*X;
        DFDY = 0.0;
        DFDZ = 1.0 - VertexCurvX*Z;
       end;
     end;}*/
}

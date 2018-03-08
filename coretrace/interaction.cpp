
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

inline double sqr(double x) { return (x)*(x); }

void Interaction(
			MTRand &myrng,
			double PosXYZ[3],
			double CosKLM[3],
			double DFXYZ[3],
			int InteractionType,
			TOpticalProperties *Opticl,
			double Wavelength,
			double PosOut[3],
			double CosOut[3],
			int *ErrorFlag )
{
/* {Purpose: To compute the direction cosines of the ray due to optical interaction
           at the intersection point of the ray with the surface
     Input - PosXYZ[2] = X,Y,Z coordinates of intersection point.
             DFXYZ     = direction numbers for the surface normal at the
                         intersection point (partial derivatives with respect
                         to X,Y,Z of surface equation)
             InteractionType = Optical interaction type indicator
                       = 1, refraction
                       = 2, reflection
                       = 3, aperture stop
                       = 4, diffraction, transmission grating
                       = 5, diffraction, reflection grating
             CosKLM[2] = direction cosines of incident ray
             Opticl    = record of optical properties
                   .RefractiveIndex[4] = Refractive index of incident and outgoing medium
                                   [0] = real part of incident medium refractive index
                                   [1] = imaginary part of ""
                                   [2] = real part of outgoing medium refractive index
                                   [4] = imaginary part of ""
                   .ApertureStopOrGratingType
                                       for InteractionType = 3, aperture stop
                                           = 1, slit
                                           = 2, elliptical
                                       for InteractionType = 4,5 grating
                                           = 1, planes parallel to Y-Z plane
                                           = 2, concentric cylinders centered about Z-axis
                   .DiffractionOrder = integral order of diffraction for InteractionTypes=4,5, grating
                   .AB12[4] = coefficients of polynomial specifying grating spacing for InteractionTypes=4,5
                        [0] = lower X limit, ApertureStopOrGratingType = 1
                              semi-X axis, ApertureStopOrGratingType = 2
                        [1] = lower Y limit, ApertureStopOrGratingType = 1
                              semi-Y axis, ApertureStopOrGratingType = 2
                        [2] = upper X limit, ApertureStopOrGratingType = 1
                              unused, ApertureStopOrGratingType = 2
                        [4] = upper Y limit, ApertureStopOrGratingType = 1
                              unused, ApertureStopOrGratingType = 2
             Wavelength = wavelength of ray

     Output - PosOut[2] = position of ray after optical interaction
              CosOut[2] = direction cosines of ray after optical interaction
              ErrorFlag = Error flag indicating successful interaction
}*/

   int i = 0;
   double CosUVW[3] = {0.0, 0.0, 0.0};
   int NIter = 0, IType = 0, NMord = 0;
   double Epsilon = 0.0, Refr1 = 0.0, Refr2 = 0.0, RMU = 0.0, RM2 = 0.0, D2 = 0.0, B = 0.0, A = 0.0, A2 = 0.0;
   double Gamn = 0.0, Gamn1 = 0.0;
   double X = 0.0,Y = 0.0,A1 = 0.0,B1 = 0.0,Ellips = 0.0,B2 = 0.0;
   double RK = 0.0,RL = 0.0,RM = 0.0,Denom,U=0,V=0,W=0;
   double Varr=0,GFactr=0,Rho2 = 0.0,Rho = 0.0,Term = 0.0,G = 0.0,D = 0.0,XX = 0.0,Ordiff = 0.0,RLamda = 0.0;
   double Rave = 0.0, Rs = 0.0, Rp = 0.0, UnitDFXYZ[3] = {0.0,0.0,0.0}, IncidentAngle = 0.0;

   NIter = 10;
	Epsilon = 0.000005;
	
	*ErrorFlag = 0;
	for (i=0;i<3;i++)
		PosOut[i] = PosXYZ[i];		 

	switch (InteractionType)
	{

/*{  InteractionType = 1, Refraction
===============================================================================}*/
	case 1:
			Refr1 = Opticl->RefractiveIndex[0];
			Refr2 = Opticl->RefractiveIndex[2];
			RMU = Refr1/Refr2;
			D2 = DOT(DFXYZ,DFXYZ);
			B = (RMU*RMU - 1.0)/D2;
			A = RMU*DOT(CosKLM, DFXYZ) /D2;
			A2 = A*A;
			if (B > A2)     //Total internal reflection
			{
				A = DOT(CosKLM, DFXYZ)/DOT(DFXYZ, DFXYZ);
				for (i=0; i<3; i++)
					CosOut[i] = CosKLM[i] - 2.0*A*DFXYZ[i];
				return;
			}

			//fresnel equations
			UnitDFXYZ[0] = -DFXYZ[0]/sqrt(DOT(DFXYZ,DFXYZ));   //unit surface normals
			UnitDFXYZ[1] = -DFXYZ[1]/sqrt(DOT(DFXYZ,DFXYZ));
			UnitDFXYZ[2] = -DFXYZ[2]/sqrt(DOT(DFXYZ,DFXYZ));
			IncidentAngle = acos(DOT(CosKLM,UnitDFXYZ));
			Rs = sqr(((Refr1*cos(IncidentAngle)-Refr2*sqrt(1-sqr(Refr1*sin(IncidentAngle)/Refr2))))/
					  ((Refr1*cos(IncidentAngle)+Refr2*sqrt(1-sqr(Refr1*sin(IncidentAngle)/Refr2)))));
			Rp = sqr(((Refr1*sqrt(1-sqr(Refr1*sin(IncidentAngle)/Refr2)))-Refr2*cos(IncidentAngle))/
					  ((Refr1*sqrt(1-sqr(Refr1*sin(IncidentAngle)/Refr2)))+Refr2*cos(IncidentAngle)));
			Rave = (Rp + Rs)/2.0;    //average of s and p polarized light; equal parts of both = non-polarized
			if (Rave < myrng())   //transmitted through surface
			{
				Gamn = -B/(2.0*A);

				//Begin Newton-Raphson loop to converge on correct root.
				for (i=1;i<NIter;i++)
				{
					Gamn1 = (Gamn*Gamn - B)/(2.0*(Gamn + A));
					if (fabs(Gamn - Gamn1) < Epsilon) goto Label_Converge;
					Gamn = Gamn1;
				}
				//Failed to converge
				*ErrorFlag = 12;
				return;
			//Have converged on Gamma, Compute direction cosines of refracted ray.
Label_Converge:
				for (i=0;i<3;i++)
					CosOut[i] = RMU*CosKLM[i] + Gamn1*DFXYZ[i];
			}
			else      //reflected from surface
			{
				A = DOT(CosKLM, DFXYZ)/DOT(DFXYZ, DFXYZ);
				for (i=0;i<3;i++)
					CosOut[i] = CosKLM[i] - 2.0*A*DFXYZ[i];
			}
			return;
		break;


/*{  InteractionType = 2, Reflection
===============================================================================}*/
	case 2:
			A = DOT(CosKLM, DFXYZ)/DOT(DFXYZ, DFXYZ);
			//Compute direction cosines for reflected ray
			for (i=0;i<3;i++)
				CosOut[i] = CosKLM[i] - 2.0*A*DFXYZ[i];

			return;
		break;


/*{  InteractionType = 3, Aperture Stop
===============================================================================}*/
	case 3:
			X = PosXYZ[0];
			Y = PosXYZ[1];
			IType = Opticl->ApertureStopOrGratingType;
			A1 = Opticl->AB12[0];
			B1 = Opticl->AB12[1];
			
			if (IType == 1)    //Slit Aperture
			{
				A2 = Opticl->AB12[2];
				B2 = Opticl->AB12[3];
				if (X < A1 || X > A2)
				{
					*ErrorFlag = 31;
					goto RayMissesAperture;
				}
				if (Y >= B1 && Y <= B2) return;
				
				*ErrorFlag = 31;
				goto RayMissesAperture;
			}
			
			if (IType == 2)      //Elliptical Aperture
			{
				Ellips = X*X/(A1*A1) + Y*Y/(B1*B1);
				if (Ellips <= 1.0) return;
				*ErrorFlag = 32;
			}

RayMissesAperture:   //Ray misses aperture
			for (i=0;i<3;i++)
				CosOut[i] = 0.0;
				
			return;

		break;


/*{  InteractionType = 4,5; Diffraction
===============================================================================}*/
	case 4:
	case 5:
			IType = Opticl->ApertureStopOrGratingType;
			NMord = Opticl->DiffractionOrder;
			Refr1 = Opticl->RefractiveIndex[0];
			Refr2 = Opticl->RefractiveIndex[2];
			RMU = Refr1/Refr2;
			D2 = DOT(DFXYZ, DFXYZ);
			RK = DFXYZ[0];
			RL = DFXYZ[1];
			RM = DFXYZ[2];
			X = PosXYZ[0];
			Y = PosXYZ[1];

			if (IType == 1)     //Parallel plane grating
			{
				Denom = RL*RL + RM*RM;
				U = 1.0/sqrt(1.0 + RK*RK/Denom);
				V = -RK*RL*U/Denom;
				W = -RK*RM*U/Denom;
				Varr = X;
				GFactr = 1.0/U;
				goto CompDiffInt;
			}

			if (IType == 2)   //Concentric Cylinder Grating
			{
				Rho2 = X*X + Y*Y;
				Rho = sqrt(Rho2);
				RM2 = RM*RM;
				Term = RL*X - RK*Y;
				G = sqrt(D2*(RM2*Rho2 + Term*Term));
				U = (RM2*X + RL*Term)/G;
				V = (RM2*Y - RK*Term)/G;
				W = -RM*(RK*X + RL*Y)/G;
				Varr = Rho;
				GFactr = Rho/(X*U + Y*V);
			}
CompDiffInt:         //Compute interaction due to diffraction
			CosUVW[0] = U;
			CosUVW[1] = V;
			CosUVW[2] = W;
			
			D = 0.0;
			XX = 1.0;

			for (i=0;i<4;i++)
			{
				D = D + Opticl->AB12[i]*XX;
				XX = XX*Varr;
			}

			D = D*GFactr;
			Ordiff = NMord;
			RLamda = Ordiff*Wavelength/(Refr2*D);
			B = (RMU*RMU - 1.0 + RLamda*RLamda - 2.0*RMU*RLamda*DOT(CosKLM, CosUVW))/D2;
			A = RMU*DOT(CosKLM, DFXYZ) /D2;
			A2 = A*A;
			if (B > A2)     //Total internal reflection
			{
				for (i=0;i<3;i++)
					CosOut[i] = 0.0;
				*ErrorFlag = 11;
				return;
			}
			
			Gamn = -B/(2.0*A);
			if (InteractionType == 5)
				Gamn = -Gamn - 2.0*A;

//Begin Newton-Raphson loop to converge on correct root.
			i=0;
			while(i++<NIter)
			{
				Gamn1 = (Gamn*Gamn - B)/(2.0*(Gamn + A));
				if (fabs(Gamn - Gamn1) < Epsilon) goto CompDCos;
				Gamn = Gamn1;
			}
//Failed to converge
			*ErrorFlag = 12;
			return;
//Have converged on Gamn1. Compute direction cosines of diffracted ray.
CompDCos:
			for (i=0;i<3;i++)
				CosOut[i] = RMU*CosKLM[i] - RLamda*CosUVW[i] + Gamn1*DFXYZ[i];

		break;
		
	default:
		break;
	}
}
//End of Procedure--------------------------------------------------------------

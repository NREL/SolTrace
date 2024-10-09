
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

#define RANGEN myrng

#define sqr(x) (x*x)

void Errors (
			MTRand &myrng,
			double CosIn[3],
			int Source,
			TSun *Sun,
			TElement *Element,
			TOpticalProperties *OptProperties, 
			double CosOut[3],
			double DFXYZ[3] ) 
{
/*{Purpose:  To add error terms to the perturbed ray at the surface in question

           Input - Seed    = Seed for RNG
                   CosIn   = Direction cosine vector of ray to which errors will be applied.
                              If Source below is 1 (i.e. sunshape) then this ray vector is before interaction with element surface
                              If Source below is 2 (i.e. surface error) then this ray vector is after interaction with element surface
                                (i.e. reflected ray or transmitted ray)

                   Source  = Source indicator flag
                           = 1 for Sunshape error (Can be gaussian, pillbox or profile data distribution)
                           = 2 for surface errors (Can be gaussian or pillbox distribution)
                   Sun     = Sun data record
                   Element = Element data record
                   DFXYZ   = surface normal vector at interaction point

           Output - CosOut  = Output direction cosine vector of ray after error terms have been included
                   }*/
	double Origin[3] = {0.0, 0.0, 0.0};
	double Euler[3] = {0.0, 0.0, 0.0};
	double PosIn[3] = {0.0, 0.0, 0.0};
	double PosOut[3] = {0.0, 0.0, 0.0};
	char dist = 'g';
	double delop=0, delop3=0, thetax=0, thetay=0, ttheta=0, theta2=0, phi=0, theta=0, stest=0;
	st_uint_t i;
	double RRefToLoc[3][3] = { {0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0} };
	double RLocToRef[3][3] = { {0.0,0.0,0.0}, {0.0,0.0,0.0}, {0.0,0.0,0.0} };


	if (CosIn[2] == 0.0)
	{
		if (CosIn[0] == 0.0)
		{
			Euler[0] = 0.0;
			Euler[1] = M_PI/2.0;
			goto Label_9;
		}
		else
		{
			Euler[0] = M_PI/2.0;
			goto Label_8;
		}
	}
     
	Euler[0] = atan2(CosIn[0],CosIn[2]);
	
Label_8:
	Euler[1] = atan2(CosIn[1],sqrt(CosIn[0]*CosIn[0]+CosIn[2]*CosIn[2]));
	
Label_9:
	Euler[2] = 0.0;

	CalculateTransformMatrices(Euler, RRefToLoc, RLocToRef);

     //g,p,d
	if (Source == 1)
	{
		dist = Sun->ShapeIndex;  //sun
		delop = Sun->Sigma/1000.0;
	}
	
	if (Source == 2)
	{
		dist = OptProperties->DistributionType; //errors
		//delop = sqrt(4.0*sqr(OptProperties->RMSSlopeError)+sqr(OptProperties->RMSSpecError))/1000.0;
		delop = OptProperties->RMSSpecError/1000.0;
	}

	unsigned int maxcall = 0;

Label_50:
	switch(dist)
	{
	case 'g':
	case 'G': //gaussian distribution
            thetax = myrng.randNorm(0., delop);
            thetay = myrng.randNorm(0., delop);

            theta2 = thetax * thetax + thetay * thetay;

		break;
	
	case 'p':
	case 'P':  //pillbox distribution
Label_200:
			thetax = 2.0*delop*RANGEN() - delop;
			thetay = 2.0*delop*RANGEN() - delop;
			theta2 = thetax*thetax + thetay*thetay;
			if (theta2 > (delop*delop)) goto Label_200;
		break;

	case 'd':
	case 'D': //sunshape data  (for sunshape only)
Label_300:
			thetax = 2.0*Sun->MaxAngle*RANGEN() - Sun->MaxAngle;
			thetay = 2.0*Sun->MaxAngle*RANGEN() - Sun->MaxAngle;
			theta2 = thetax*thetax + thetay*thetay;
			theta = sqrt(theta2);  //wendelin 1-9-12  do the test once on theta NOT individually on thetax and thetay as before

			i = 0;
			while ( i < Sun->SunShapeAngle.size()-1
				&& Sun->SunShapeAngle[i] < theta)
				i++;
			
			if (i == 0)
				stest = Sun->SunShapeIntensity[0];
			else  //change from average interpolation between data points to linear interpolation  12-20-11 wendelin
                stest = Sun->SunShapeIntensity[i-1]
                        + (Sun->SunShapeIntensity[i] - Sun->SunShapeIntensity[i-1])*(theta - Sun->SunShapeAngle[i-1])/
						  (Sun->SunShapeAngle[i] - Sun->SunShapeAngle[i-1]);
				//stest = (Sun->SunShapeIntensity[i] + Sun->SunShapeIntensity[i-1])/2.0;

			if (RANGEN() > (stest/Sun->MaxIntensity)) goto Label_300;

			if (theta2 > (Sun->MaxAngle*Sun->MaxAngle)) goto Label_300;
			theta2 = theta2/1000000.0;
		break;

	case 'f':  //gray diffuse distribution
	case 'F':
		theta2 = pow(asin(sqrt(RANGEN())),2);
		break;

	}

     /*{Transform to local coordinate system of ray to set up rotation matrices for coord and inverse
       transforms}*/
	TransformToLocal(PosIn, CosIn, Origin, RRefToLoc, PosOut, CosOut);

	// {Generate errors in terms of direction cosines in local ray coordinate system}
	theta = sqrt(theta2);

	//phi = atan2(thetay, thetax); //This function appears to  present irregularities that bias results incorrectly for small values of thetay or thetax
	phi = RANGEN()*2.0*3.1415926535897932385; // Therefore have chosen to randomize phi rather than calculate from randomized theta components
                                                      //  obtained from the distribution. The two approaches are equivalent save for this issue with
                                                      //  arctan2.      wendelin 01-12-11 


	CosOut[0] = sin (theta)*cos(phi);
	CosOut[1] = sin (theta)*sin(phi);
	CosOut[2] = cos(theta);

	for (i=0;i<3;i++)
	{
		PosIn[i] = PosOut[i];
		CosIn[i] = CosOut[i];
	}

     //{Transform perturbed ray back to element system}
     TransformToReference(PosIn, CosIn, Origin, RLocToRef, PosOut, CosOut);

     /*{If reflection error applicaton and new ray direction (after errors) physically goes through opaque surface,
     then go back and get new perturbation 06-12-07}*/
     if ( (Source == 2) && (Element->InteractionType == 2) 
			&& (DOT(CosOut, DFXYZ) < 0)
			&& maxcall++ < 50000)
			goto Label_50;
}
//End of Procedure--------------------------------------------------------------

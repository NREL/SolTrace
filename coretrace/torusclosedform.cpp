
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


void TorusClosedForm(
			TElement *Element,
			double PosLoc[3],
			double CosLoc[3],
			double PosXYZ[3],
			double DFXYZ[3],
			double *PathLength,
			int *ErrorFlag)
{
	double Xo=0.0,Yo=0.0,Zo=0.0,Epsilon=0.0,Eta=0.0,Rho=0.0,Rs=0.0,Ra=0.0,X=0.0,Y=0.0,Z=0.0,Fx=0.0,Fy=0.0,Fz=0.0;
	double amatrix[5][5];
	double rvector[5];
	double imagroot1=0.0,imagroot2=0.0;
	int nn = 0;

	Rs = Element->CrossSectionRadius;
	Ra = Element->AnnularRadius;
	Xo = PosLoc[0];
	Yo = PosLoc[1];
	Zo = PosLoc[2];
	Epsilon = CosLoc[0];
	Eta = CosLoc[1];
	Rho = CosLoc[2];
	nn = 4;
	*ErrorFlag = 0;

	for (int i=0;i<5;i++)
	{
		rvector[i] = 0.0;
		for (int j=0;j<5;j++)
			amatrix[i][j] = 0.0;
	}

	amatrix[nn][4] = pow(Epsilon,4)+2.0*Epsilon*Epsilon*(Eta*Eta+Rho*Rho)+
						pow(Eta,4)+2.0*Eta*Eta*Rho*Rho+pow(Rho,4);

	amatrix[nn][3] = 4.0*(Epsilon*Epsilon+Eta*Eta+Rho*Rho)*(Epsilon*Xo+Eta*Yo+
						Rho*Zo-Rho*Rs);
	
	amatrix[nn][2] = Xo*Xo*(6.0*Epsilon*Epsilon+2.0*Eta*Eta+2.0*Rho*Rho)+
						8.0*Epsilon*Xo*(Eta*Yo+Rho*Zo-Rho*Rs)+
						2.0*Yo*Yo*(Epsilon*Epsilon+3.0*Eta*Eta+Rho*Rho)+
						8.0*Eta*Rho*Yo*(Zo-Rs)+
						(Epsilon*Epsilon+Eta*Eta+3.0*Rho*Rho)*(2.0*Zo*Zo-4.0*Rs*Zo)-
						2.0*Ra*Ra*(Epsilon*Epsilon+Eta*Eta-Rho*Rho)+4.0*Rho*Rho*Rs*Rs;

	amatrix[nn][1] = 4.0*(Xo*Xo*(Epsilon*Xo+Eta*Yo+Rho*Zo-Rho*Rs)+Yo*Yo*(Epsilon*Xo+
						Eta*Yo+Rho*Zo-Rho*Rs)+Zo*Zo*(Epsilon*Xo+Eta*Yo+Rho*Zo-3.0*Rho*Rs)-
						2.0*Epsilon*Rs*Xo*Zo-Epsilon*Ra*Ra*Xo-2.0*Eta*Rs*Yo*Zo-
						Eta*Ra*Ra*Yo+Rho*Ra*Ra*(Zo-Rs)+2.0*Rho*Rs*Rs*Zo);

	amatrix[nn][0] = pow(Xo,4)+2.0*Xo*Xo*(Yo*Yo+Zo*Zo-2.0*Rs*Zo-Ra*Ra)+
						pow(Yo,4)+2.0*Yo*Yo*(Zo*Zo-2.0*Rs*Zo-Ra*Ra)+
						pow(Zo,4)-4.0*Rs*Zo*Zo*Zo+2.0*Ra*Ra*Zo*Zo+4.0*Rs*Rs*Zo*Zo
						-4.0*Ra*Ra*Rs*Zo+pow(Ra,4);

    Root_432(nn,amatrix,rvector,&imagroot1,&imagroot2);

 // {pass results}
    if (imagroot1==0.0 && imagroot2==0.0)
    {
      //four real roots rvector[1-4] which are the pathlengths
		piksrt(4,rvector); //sort them
		if (rvector[4] <= 0.0)  //ray heading away from surface
		{
			*PathLength = rvector[4];
			*ErrorFlag = 1;
			return;
		}
		*PathLength = rvector[4];
		if (rvector[3] <= 0.0 )
			goto Label_10;
		else
			*PathLength = rvector[3];
			
		if (rvector[2] <= 0.0 )
			goto Label_10;
		else
			*PathLength = rvector[2];
			
		if (rvector[1] <= 0.0 )
			goto Label_10;
		else
			*PathLength = rvector[1];
			
		goto Label_10;
	}


	if (imagroot1==0.0 && imagroot2!=0.0)
	{
		//two real roots rvector[1-2] which are the pathlengths
		if (rvector[1] <= rvector[2])
			*PathLength = rvector[1];
		else
			*PathLength = rvector[2];
		goto Label_10;
	}

    if (imagroot1!=0.0 && imagroot2==0.0)
    {
      //two real roots rvector[3-4] which are the pathlengths
      if (rvector[3] <= rvector[4])
		*PathLength = rvector[3];
      else
		*PathLength = rvector[4];
		
      goto Label_10;
    }

    *PathLength = 0.0;        
    
Label_10:
	if (*PathLength==0.0)   //ray missed torus completely
	{
		*ErrorFlag = 1;
		return;
	}
	
	X = Xo+*PathLength*Epsilon;
	Y = Yo+*PathLength*Eta;
	Z = Zo+*PathLength*Rho;
	Fx = -2.0*X*(Ra-sqrt(X*X+Y*Y))/sqrt(X*X+Y*Y);
	Fy = -2.0*Y*(Ra-sqrt(X*X+Y*Y))/sqrt(X*X+Y*Y);
	Fz = 2.0*(Z-Rs);
	PosXYZ[0] = X;
	PosXYZ[1] = Y;
	PosXYZ[2] = Z;
	DFXYZ[0] = -Fx;
	DFXYZ[1] = -Fy;
	DFXYZ[2] = -Fz;
}
//End of Procedure--------------------------------------------------------------

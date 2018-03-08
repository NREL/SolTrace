
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
#include <stdlib.h>

#include "types.h"
#include "procs.h"


void MatrixVectorMult(double M[3][3], double V[3], double MxV[3])
{
/*{Purpose: To perform multiplication of a matrix (3,3) and a vector (3) to result
          in new vector (3)

          Input -
                Matrix = Matrix
                Vector = Vector
          Output -
                 MxV = Matrix*Vector = new vector}*/

	MxV[0] = M[0][0]*V[0] + M[0][1]*V[1] + M[0][2]*V[2];
	MxV[1] = M[1][0]*V[0] + M[1][1]*V[1] + M[1][2]*V[2];
	MxV[2] = M[2][0]*V[0] + M[2][1]*V[1] + M[2][2]*V[2];
/*
	int i, j;
	for (i=0;i<3;i++)
	{
		MxV[i] = 0.0;
		for (j=0;j<3;j++)
			MxV[i] = Matrix[i][j]*Vector[j] + MxV[i];
	}*/
}
//End of procedure-------------------------------------------------------------


double DOT(double A[3], double B[3])
{
//{Purpose: To compute the dot product of 2 N-dimensional vectors, A and B
  //        Input -
  //              A(N) = First input vector
  //              B(N) = Second input vector
  //              N = dimension of vectors
  //        Output -
  //             Result of DOT = dot product of A and B}

	return (A[0]*B[0] + A[1]*B[1] + A[2]*B[2]);
}


void MatrixTranspose(double InputMatrix[3][3], int NumRowsCols, double OutputMatrix[3][3])
{
/*{Purpose: To obtain the transpose of a square matrix InputMatrix(NumRowsCols,NumRowscols)
         Input -
               InputMatrix = Input Matrix
               NumRowsCols = Number of rows/columns in InputMatrix
         Output -
                OutputMatrix = Output Matrix OutputMatrix(i,j) = InputMatrix(j,i)}*/
	for (int i=0;i<NumRowsCols;i++)
		for (int j=0;j<NumRowsCols;j++)
			OutputMatrix[i][j] = InputMatrix[j][i];
}
//end of procedure--------------------------------------------------------------


//end of procedure--------------------------------------------------------------

void TransformToLocal(double PosRef[3], double CosRef[3], double Origin[3], 
		double RRefToLoc[3][3], 
		double PosLoc[3], double CosLoc[3])
{
/*{Purpose:  To perform coordinate transformation from reference system to local
           system.
           Input -
                 PosRef = X,Y,Z coordinates of ray point in reference system
                 CosRef = Direction cosines of ray in reference system
                 Origin = X,Y,Z coordinates of origin of local system as measured
                          in reference system
                 RRefToLoc = Rotation matrices required for coordinate transform
                             from reference to local
           Output -
                 PosLoc = X,Y,Z coordinates of ray point in local system
                 CosLoc = Direction cosines of ray in local system }*/

	double PosDum[3];

/*{Multiply the position vector and the direction cosine vector by the transformation
 matrix to get the new vectors in the local system.  The position vector is first
 referenced to the origin of the local frame.}*/
	for (int i=0;i<3;i++)
		PosDum[i] = PosRef[i] - Origin[i];

	MatrixVectorMult(RRefToLoc, PosDum, PosLoc);
	MatrixVectorMult(RRefToLoc, CosRef, CosLoc);
}
//end of procedure--------------------------------------------------------------


void TransformToReference(double PosLoc[3], double CosLoc[3], double Origin[3], 
		double RLocToRef[3][3], 
		double PosRef[3], double CosRef[3])
{
/*{Purpose:  To perform coordinate transformation from local system to reference
           system.
           Input -
                 PosLoc = X,Y,Z coordinates of ray point in local system
                 CosLoc = Direction cosines of ray in Loc system
                 Origin = X,Y,Z coordinates of origin of local system as measured
                          in reference system
                 RLocToRef = Rotation matrices required for coordinate transform
                             from local to reference (inverse of reference to
                             local transformation)
           Output -
                 PosRef = X,Y,Z coordinates of ray point in reference system
                 CosRef = Direction cosines of ray in reference system}*/

	double PosDum[3];
	
/*{Use previously calculated RLocToRef matrix (in TransformToLocal) to obtain the
 inverse transformation back to Reference system.}*/
     MatrixVectorMult(RLocToRef, PosLoc, PosDum);
     MatrixVectorMult(RLocToRef, CosLoc, CosRef);

	for (int i=0;i<3;i++)
		PosRef[i] = PosDum[i] + Origin[i];
}
//end of procedure--------------------------------------------------------------
void CalculateTransformMatrices(double Euler[3], double RRefToLoc[3][3], double RLocToRef[3][3])
{
/*{input:  Euler = Euler angles
 output: RRefToLoc = Transformation matrix from Reference to Local system
         RLocToRef = ""             ""     ""   Local to Reference system (transpose of above)} */
	double Alpha = 0.0;
	double Beta = 0.0;
	double Gamma = 0.0;
	double CosAlpha = 0.0;
	double CosBeta = 0.0;
	double CosGamma = 0.0;
	double SinAlpha = 0.0;
	double SinBeta = 0.0;
	double SinGamma = 0.0;
	
/*{Offload Alpha, Beta and Gamma: the three Euler rotations from Ref frame to Local
 frame. Also calculate their cosines.}*/
	Alpha = Euler[0];
	Beta  = Euler[1];
	Gamma = Euler[2];
	CosAlpha = cos(Alpha);
	CosBeta  = cos(Beta);
	CosGamma = cos(Gamma);
	SinAlpha = sin(Alpha);
	SinBeta  = sin(Beta);
	SinGamma = sin(Gamma);

/*{Fill in elements of the transformation matrix as per Spencer and Murty paper
 page 673 equation (2)}*/
	RRefToLoc[0][0] = CosAlpha*CosGamma + SinAlpha*SinBeta*SinGamma;
	RRefToLoc[0][1] = -CosBeta*SinGamma;
	RRefToLoc[0][2] = -SinAlpha*CosGamma + CosAlpha*SinBeta*SinGamma;
	RRefToLoc[1][0] = CosAlpha*SinGamma - SinAlpha*SinBeta*CosGamma;
	RRefToLoc[1][1] = CosBeta*CosGamma;
	RRefToLoc[1][2] = -SinAlpha*SinGamma - CosAlpha*SinBeta*CosGamma;
	RRefToLoc[2][0] = SinAlpha*CosBeta;
	RRefToLoc[2][1] = SinBeta;
	RRefToLoc[2][2] = CosAlpha*CosBeta;

/*{Transpose the matrix to get the inverse transformation matrix for going back
 to reference system.  This is used by the TransformToReference procedure.}*/
	MatrixTranspose(RRefToLoc, 3, RLocToRef);
}

//end of procedure--------------------------------------------------------------
void EvalPoly(double ax, double ay, std::vector<double> &Coeffs, int POrder, double *az) //the 0.0's are values for DeltaX and DeltaY; **[need to look at this further]**
{
    *az = 0.0;
	double r = sqrt(ax*ax + ay*ay);
	for (int i=0;i<=POrder;i++)
      *az = *az + Coeffs[i]*pow(r,i);
}

//end of procedure--------------------------------------------------------------
void PolySlope( std::vector<double> &Coeffs, int POrder, double ax, double ay, double *dzdx, double *dzdy)
{
	double dzdr = 0.0;
	double r = sqrt(ax*ax + ay*ay);
	double drdx = ax/r;
	double drdy = ay/r;
	for (int i=1;i<=POrder;i++)
		dzdr = dzdr + i*Coeffs[i]*pow(r, i-1);
    
    *dzdx = dzdr*drdx;
    *dzdy = dzdr*drdy;
}
//end of procedure--------------------------------------------------------------


bool splint( std::vector<double> &xa,
			std::vector<double> &ya,
			std::vector<double> &y2a,
			int n,
			double x,
			double *y,
			double *dydx )
{
	int klo = 0,khi = 0,k = 0;
	double h = 0.0,b = 0.0,a = 0.0;

	klo = 0;
	khi = n-1;
	while( khi-klo > 1)
	{
		k = (khi+klo) / 2;
		if( xa[k] > x )
			khi = k;
		else
			klo = k;
	}
	
	h = xa[khi]-xa[klo];
	if (h != 0.0)
	{
	
		a = (xa[khi]-x)/h;
		b = (x-xa[klo])/h;
		*y = a*ya[klo]+b*ya[khi]+
			((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;
		*dydx = (ya[khi]-ya[klo])/(xa[khi]-xa[klo])-
			(3.0*a*a-1.0)*(xa[khi]-xa[klo])*y2a[klo]/6.0 + (3.0*b*b-1.0)*(xa[khi]-xa[klo])*y2a[khi]/6.0;
	}
	else
	{
		return false;
	}

	return true;
}

//end of procedure--------------------------------------------------------------

void spline( std::vector<double> &x, 
			std::vector<double> &y,
			int n,
			double yp1, double ypn,
			std::vector<double> &y2 )
{
	int i=0,k=0;
	double p = 0.0,qn = 0.0,sig = 0.0,un = 0.0;
	std::vector<double> u( x.size() );

	if (yp1 > 0.99e30)
	{
		y2[0] = 0.0;
		u[0] = 0.0;
	}
	else
	{
		y2[0] = -0.5;
		u[0] = (3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
	}
	for (i=2;i<n-1;i++)
	{
		sig = (x[i]-x[i-1])/(x[i+1]-x[i-1]);
		p = sig*y2[i-1]+2.0;
		y2[i] = (sig-1.0)/p;
		u[i] = (y[i+1]-y[i])/(x[i+1]-x[i])-(y[i]-y[i-1])/(x[i]-x[i-1]);
		u[i] = (6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
	}
	
	if (ypn > 0.99e30 )
	{
		qn = 0.0;
		un = 0.0;
	}
	else
	{
		qn = 0.5;
		un = (3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
	}
	y2[n-1] = (un-qn*u[n-2])/(qn*y2[n-2]+1.0);
	for (k=n-2;k>=0;k--)
	  y2[k] = y2[k]*y2[k+1]+u[k];

}
//end of procedure--------------------------------------------------------------

void piksrt(int n, double arr[5] )
{
	int i;
	for (int j=1;j<n;j++)
	{
		double a = arr[j];
		for (i = j-1; i >= 0; i--)
		{
			if (arr[i] <= a ) goto Label_10;
			arr[i+1] = arr[i];
		}
		i = 0;
Label_10:
		arr[i+1] = a;
	}
}
//end of procedure--------------------------------------------------------------

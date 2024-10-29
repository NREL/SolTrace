
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
#include <algorithm>

#include "types.h"
#include "procs.h"
#include "interpolate.h"

void EvalMono(double ax, double ay, HPM2D &B, int order, double DeltaX, double DeltaY, double *az)
{
/*{This procedure evaluates monomial expressions.
 Input -   ax = location in x direction
           ay = location in y direction
           b = array of coefficients
           order = order of monomial
           DeltaX = offset of dataset in x direction
           DeltaY = offset of dataset in y direction

 Output -   az = result of evaluating series at point ax, ay */
 
 /* {It should be pointed out that the form of the Zernike equation as presented in Malacara is
        Z(x,y) = Sum Sum Bij x^j y^i-j.  In the VSHOT system, this was rotated such that
        Z(x,y) = Sum Sum Bij x^i-j y^j.  Soltrace was made to be consistent with this form. However
        as of 02/14/03, the VSHOT was made to be consistent with the correct form in Malacara (i.e. the first
        equation above. Subsequently, SolTRACE was changed to reflect this. Thus the following  makes sense.
        The important terms B(1,0) B(1,1) B(2,0) and B(2,2) are
        such that B(1,0) = tilt about the x-axis
                  B(1,1) = tilt about the y-axis
                  B(2,0) = best focal distance in y-direction
                  B(2,2) = best focal distance in x-direction} */

	double XP = 0.0, YP = 0.0, z = 0.0, XX = 0.0, YY = 0.0;
	int i=0,j=0;

	//initialize sum
	XP = ax - DeltaX;
	YP = ay - DeltaY;
	z = 0.0;

	//sum over i from 0 to order of monomials
	for (i = 0;i<=order;i++)
	{
		//sum over j from 0 to i
		for (j = 0; j <= i;j++)
		{
			if (j == 0)
				XX = 1.0;
			else
				XX = pow(XP, j);

			if (i-j == 0)
				YY = 1.0;
			else
				YY = pow(YP, (i-j));

			z = z + B.at(i,j)*XX*YY;
		}
	}

	*az = z;
}
//End of Procedure--------------------------------------------------------------

void MonoSlope(HPM2D &B, int order, double sxp, double syp, double *dzdx, double *dzdy)
{
/* {It should be pointed out that the form of the Zernike equation as presented in Malacara is
        Z(x,y) = Sum Sum Bij x^j y^i-j.  In the VSHOT system, this was rotated such that
        Z(x,y) = Sum Sum Bij x^i-j y^j.  However
        as of 02/14/03, the VSHOT was made to be consistent with the correct form in Malacara (i.e. the first
        equation above. Subsequently, SolTRACE was changed to reflect this.
        Thus the following derivatives make sense.}*/
	
    if (sxp == 0.0) sxp = sxp + 1.0e-10;
    if (syp == 0.0) syp = syp + 1.0e-10;
    
    *dzdx = 0.0;
    *dzdy = 0.0;
    
	for (int i = 1; i <= order; i++)
    {
		for (int j = 0; j<=i; j++)
		{
            //dzdx = dzdx + (i-j)*B[i,j]*power(sxp, (i-j-1))*power(syp,j);
            //dzdy = dzdy + j*B[i,j]*power(sxp, (i-j))*power(syp, (j-1));
			*dzdx = *dzdx +   j   * B.at(i,j) * pow(sxp, (j-1)) * pow(syp, (i-j));
			*dzdy = *dzdy + (i-j) * B.at(i,j) * pow(sxp, j)     * pow(syp, (i-j-1));
		}
	}
}

//End of Procedure--------------------------------------------------------------

void FEInterpGM(double Xray, double Yray, GaussMarkov* gm, double* zr)
{
	/*
	
	*/
	VectDoub xy = { Xray, Yray };

	*zr = gm->interp(xy);
}

struct __r_el
{
	double r;
	VectDoub* v;
	__r_el() {};
	__r_el(double _r, VectDoub* _v) { r = _r; v = _v; };
};

bool __rcomp(__r_el &a, __r_el &b)
{
	return a.r < b.r;
}

void FEInterpKD(double Xray, double Yray, FEDataObj* kd, double step, double* zr, double* dzrdx, double* dzrdy)
{
	/*
	Interpolate using a local kriging model of the surface. The model is constructed using adjacent nodes to the 
	ray position {Xray,Yray} determined using a K-d tree. The derivatives of the surface in x and y are also
	computed numerically using the kriging model.
	*/

	// Retrieve ND nearest neighbors and linearly interpolate 
	std::vector<void*> nearby_nodes;
	bool ok = kd->get_all_data_at_loc(nearby_nodes, Xray, Yray);
	if (!ok || nearby_nodes.size() < 3)
		throw std::runtime_error("Interpolation failed to find nearby nodes to ray hit in FEA surface.");
	
	int nn = (int)nearby_nodes.size();
	
	std::vector<__r_el> ndist;

	//sort to find nearest neighbors
	for (int i = 0; i < nn; i++)
	{
		VectDoub* v = static_cast<VectDoub*>(nearby_nodes.at(i));

		double r = std::sqrt(std::pow(Xray - v->at(0), 2) + std::pow(Yray - v->at(1), 2));

		ndist.push_back(__r_el(r, v));
	}

	// Find the nearest points
	std::sort(ndist.begin(), ndist.end(), __rcomp);

	// build a local model of the surface using kriging. Use nearest 7 points max (number selected based on limited testing)
	GaussMarkov gm;
	for (int i = 0; i < (nn > 7 ? 7 : nn); i++)
	{
		gm.x.push_back({ ndist.at(i).v->at(0), ndist.at(i).v->at(1) });
		gm.y.push_back(ndist.at(i).v->at(2));
	}

	gm.setup(1.999,0);

	VectDoub xy = { Xray,Yray };
	*zr = gm.interp(xy);

	// evaluate the slopes using numerical derivative
	xy.at(0) += step;
	double za = gm.interp(xy);
	xy.at(0) -= step;
	xy.at(1) += step;
	double zb = gm.interp(xy);
	*dzrdx = (*zr - za) / step;
	*dzrdy = (*zr - zb) / step;

	//--> the following was also tried. Barycentric interpolation doesn't work very well for low-quality meshes.
	//keep this code for now as a reference.
	//VectDoub* p1 = ndist.at(0).v;
	//VectDoub* p2 = ndist.at(1).v;
	//VectDoub* p3 = ndist.at(2).v;
	//
	//double det = (p2->at(1) - p3->at(1)) * (p1->at(0) - p3->at(0)) + (p3->at(0) - p2->at(0)) * (p1->at(1) - p3->at(1));
	//double u = ((p2->at(1) - p3->at(1)) * (Xray - p3->at(0)) + (p3->at(0) - p2->at(0)) * (Yray - p3->at(1))) / det;
	//double v = ((p3->at(1) - p1->at(1)) * (Xray - p3->at(0)) + (p1->at(0) - p3->at(0)) * (Yray - p3->at(1))) / det;
	//double w = (1 - u - v);
	//*zr = u * p1->at(2) + v * p2->at(2) + w * p3->at(2);
	//// slopes
	//VectDoub r_12 = { p2->at(0) - p1->at(0), p2->at(1) - p1->at(1), p2->at(2) - p1->at(2) };
	//VectDoub r_13 = { p3->at(0) - p1->at(0), p3->at(1) - p1->at(1), p3->at(2) - p1->at(2) };
	//double a = r_12.at(1) * r_13.at(2) - r_12.at(2) * r_13.at(1);
	//double b = r_12.at(2) * r_13.at(0) - r_12.at(0) * r_13.at(2);
	////double c = r_12.at(0) * r_13.at(1) - r_12.at(1) * r_13.at(0);
	//double tiny = 1e-19;
	//*dzrdx = a / (det == 0 ? tiny : det);
	//*dzrdy = b / (det == 0 ? tiny : det);
	//<---

	return;
}


void FEInterpNew(double Xray, double Yray, double Density,
			HPM2D &FEData, int NumFEPoints,
			double *zr)
{

/*{Interpolation scheme for random finite element data. Given a location defined by
xray and yray, interpolates to find the best guess for the corresponding residual z.
    Input - Xray = x coordinate of incoming ray
            Yray = y coordinate of incoming ray
            Density = density of surface data points (np/pi*radish^2)
            FEData = array of FE data
                        (in the form of FEData[1] = x location of data point
                                        FEData[2] = y location of data point
                                        FEData[3] = z location of data point )
            NumFEPoints = number of FE datum points

    Output - zr = best  z
}*/
	double SUMRR2 = 0.0, sumwtr = 0.0, XX = 0.0, YY = 0.0, R2 = 0.0;
	
    //initialize
    SUMRR2 = 0.0;
    sumwtr = 0.0;
	for (int i=0;i<NumFEPoints;i++)
    {
        XX = FEData.at(i,0) - Xray;
        YY = FEData.at(i,1) - Yray;

        //test to see if ray coincides with data point
        if (XX == 0.0 && YY == 0.0)
        {
            *zr = FEData.at(i,2);
            return;
		}

        R2 = XX*XX + YY*YY;
        SUMRR2 = SUMRR2 + 1.0/R2;
        sumwtr = sumwtr + FEData.at(i,2)/R2;
	}

    //zr = sumwtr/(Density + SUMRR2);
    *zr = sumwtr/(SUMRR2);
}
//End of Procedure--------------------------------------------------------------

void VSHOTInterpolateModShepard( double Xray, double Yray, double Density,
			HPM2D &VSHOTData, int NumVSHOTPoints,
			double *zx, double *zy,
			int *ErrorFlag)
{
/*{Interpolation scheme for VSHOT data. This method is a modified version of Shepard's method.
Shepard, Donald (1968). "A two-dimensional interpolation function for irregularly-spaced data". Proceedings of the 1968 ACM National Conference. pp. 517–524
Given a location defined by xray and yray, interpolates to find the best guess for the corresponding measured slope.
This  version assumes a dual curvature or axisymmetric optical element.  Weighting is a function of the actual radial distance from the ray
intersection point to VSHOT data points within a radius defined by the density of points and a chosen limited number of nearest neighbors.
The  slopes in both the x and y directions should be affected by all points equally.

	Input - Xray = x coordinate of incoming ray
			Yray = y coordinate of incoming ray
			Density = density of surface data points (np/pi*radish^2)
			VSHOTData = array of VSHOT data
			NumVSHOTPoints = number of VSHOT datum points

	Output - zx = best slope in x-direction
			 zy = best slope in y-direction}

var
	SUMRR2: real;
	sumwtx: real;
	sumwty: real;
	i: integer;
	XX: real;
	YY: real;
	R2: real;
	ROIRadius2: real;*/

	double SUMRR2 = 0.0, sumwtx = 0.0, sumwty = 0.0;
	int i=0;
	double XX = 0.0, YY = 0.0, R2 = 0.0, ROIRadius2 = 0.0;
	double *vpt = 0;

	*ErrorFlag = 0;

	//initialize
	ROIRadius2 = 30.0/Density/3.14159;
	SUMRR2 = 0.0;
	sumwtx = 0.0;
	sumwty = 0.0;
	for (i = 0; i < NumVSHOTPoints; i++)
	{
		vpt = VSHOTData.data()+VSHOTData.ncols()*i; // pointer arithmetic for top performance

		XX = vpt[0] - Xray;
		YY = vpt[1] - Yray;

		//test to see if ray coincides with data point
		if (XX == 0.0 && YY == 0.0)
		{
			*zx = vpt[2];     //measured slopes, not slope residuals
			*zy = vpt[3];
			return;
		}

		R2 = XX*XX + YY*YY;
		if (R2 < ROIRadius2)
		{
			SUMRR2 = SUMRR2 + 1.0/R2;
			sumwtx = sumwtx + vpt[2]/R2;  //measured slopes, not slope residuals
			sumwty = sumwty + vpt[3]/R2;
		}
	}

	if ( SUMRR2 == 0.0 )
	{
		*zx = 0.0;
		*zy = 0.0;
		//*ErrorFlag = 1;
		return;
	}

	*zx = sumwtx/(SUMRR2);
	*zy = sumwty/(SUMRR2);
}
//End of Procedure--------------------------------------------------------------


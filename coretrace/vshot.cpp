#include <math.h>
#include <stdlib.h>

#include "types.h"
#include "procs.h"


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
void FEInterpolate(double Xray, double Yray, double Delta, double Density, 
			double *FEData[5], int NumFEPoints, 
			double *z, double *zx, double *zy)
{
/* {Interpolation scheme for random finite element data. Given a location defined by
xray and yray, interpolates to find the best guess for the corresponding residual z.
    Input - Xray = x coordinate of incoming ray
            Yray = y coordinate of incoming ray
            Delta = change in x-y coordinates used to calculate the corresponding slopes.
            Density = density of surface data points (np/pi*radish^2)
            FEData = array of FE data
                        (in the form of FEData[1] = x location of data point
                                        FEData[2] = y  location of data point
                                        FEData[3] = Residual z (displacement in z)) ???
            NumFEPoints = number of FE datum points

    Output - z = best residual z
             zx = best residual z at point Delta away in x-direction
             zy = best residual z at point delta away in y-direction}*/

	double SUMRR2 = 0.0, SUMWHT = 0.0, sumrx2 = 0.0, sumry2 = 0.0, sumwtx = 0.0, sumwty = 0.0, Delta2 = 0.0;
	int i=0;
	double XX = 0.0, YY = 0.0, rx2 = 0.0, ry2 = 0.0, R2 = 0.0;
	
    //initialize
    Delta2 = Delta*Delta;
    SUMRR2 = 0.0;
    SUMWHT = 0.0;
    sumrx2 = 0.0;
    sumry2 = 0.0;
    sumwtx = 0.0;
    sumwty = 0.0;

	for (i = 0;i< NumFEPoints;i++)
	{
		XX = FEData[i][0] - Xray;
		YY = FEData[i][1] - Yray;

		//test to see if ray coincides with datum point
		if (XX == 0.0 && YY == 0.0)
		{
			*z = FEData[i][2];
			*zx = *z;
			*zy = *z;
			return;
		}

		R2 = XX*XX + YY*YY;
		rx2 = R2 - 2.0*Delta*XX + Delta2;
		ry2 = R2 - 2.0*Delta*YY + Delta2;
		SUMRR2 = SUMRR2 + 1.0/R2;
		sumrx2 = sumrx2 + 1.0/rx2;
		sumry2 = sumry2 + 1.0/ry2;
		SUMWHT = SUMWHT + FEData[i][2]/R2;
		sumwtx = sumwtx + FEData[i][2]/rx2;
		sumwty = sumwty + FEData[i][2]/ry2;
	}

    *z = SUMWHT/(Density + SUMRR2);
    *zx = sumwtx/(Density + sumrx2);
    *zy = sumwty/(Density + sumry2);
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


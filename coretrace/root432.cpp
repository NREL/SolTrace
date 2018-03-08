
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

class CR432
{
public:

	double a[5][5];
	double r[5];
	int i,j,n;
	double aa,b,c,d,ii,im,k,l,m,q,rr,s,sw;


	CR432() {
		for (i=0;i<5;i++)
		{
			r[i] = 0.0;
			for (j=0;j<5;j++)
				a[i][j] = 0.0;
		}

		i=j=n=0;
		aa=b=c=d=ii=im=k=l=m=q=rr=s=sw=0.0;
	}

	double f(double x)
	{
		return  x * x * x + a[3][2] * x * x + a[3][1] * x + a[3][0];
	}

	//*****************************************************************************

	/*{******************************************
	* This Subroutine calculates the roots of *
	*    X^2 + A(2,1)*X + A(2,0) = 0          *
	*    W = Determinant of equation          *
	******************************************}*/
	void Root_2()
	{
		double q1=0.0, q2=0.0, w=0.0;
		w = a[2][1] * a[2][1] - 4 * a[2][0];
		q1 = -a[2][1] / 2.0;
		if (w < 0)
		{
			q2 = sqrt(-w) / 2.0;
			im = q2;
			r[1] = q1; r[2] = q1;
		}
		else if (w==0)
		{
			r[1] = q1; r[2] = q1; im = 0.0;
		}
		else if (w > 0)
		{
			q2 = sqrt(w) / 2.0; im = 0.0;
			r[1] = q1 + q2; r[2] = q1 - q2;
		}
	}
	//*****************************************************************************

	/*{*******************************************
	* This subroutine calculates the roots of  *
	* X^3 + A(3,2)*X^2 + A(3,1)*X + A(3,0) = 0 *
	*******************************************}*/
	void Root_3()
	{
		double am=0, er=0, te=0, tt=0, xa=0, xb=0, xc=0, y1=0, y2=0;
		
	  //{one root equals zero}
		if (a[3][0] == 0)
		{
			xc = 0; goto Label_500;
		}
		//{looking for maximum coefficient in absolute magnitude}
		am = fabs(a[3][0]);
		for (i=1;i<=3;i++)
		{
			tt = fabs(a[3][i]);
			if (am < tt) am = tt;
		}
		//{Define interval where a real root exists
		// according to sign of A(3,0)  }
		if (a[3][0] > 0)
		{
			aa = -am - 1;
			b = 0;
			goto Label_100;
		}
		aa = 0; b = am + 1;

	Label_100:
		//{Searching for xc = real root in interval (a,b)
		// (by Bisection method)  }

		//{Define segment (xa,xb) including the root}
		xa = aa; xb = b; y1 = f(aa); y2 = f(b);

		//{Desired precision}
		er = 0.000001;
		if (y1 == 0)
		{
			xc = aa; goto Label_500;
		}
		
		if (y2 == 0 )
		{
			xc = b; goto Label_500;
		}
	//  {xc : middle of segment}
	Label_200:
		xc = (xa + xb) / 2.0; te = f(xc);
		if (te == 0) goto Label_500;
		if (xb - xa < er) goto Label_500;
		if (f(xa) * te > 0)
		{
			xa = xc; goto Label_200;
		}
		xb = xc; goto Label_200;
		
	//	{r[3] is a real root}
	Label_500:
		r[3] = xc;
		if (sw == -1) return;
		//{Calculates the roots of remaining quadratic equation
		//Define the equation coefficients }
		a[2][1] = a[3][2] + xc;
		a[2][0] = a[3][1] + a[3][2] * xc + xc * xc;
		Root_2();
	}
	//*****************************************************************************

	//{Root search main subroutine}
	void Root_4()
	{
		//{Normalization of coefficients }
		for (i=0;i<n;i++)
			a[n][i] = a[n][i] / a[n][n];
			
		//{Branching according to degree n}
		if (n==2)
		{
			Root_2();
			return;
		}
		else if (n==3)
		{
			Root_3();
			return;
		}
		else if (n==4)
		{
			aa = a[4][3];
			b = a[4][2];
			c = a[4][1];
			d = a[4][0];
			q = b - (3.0 * aa * aa / 8.0);
			rr = c - (aa * b / 2) + (aa * aa * aa / 8.0);
			s = d - (aa * c / 4) + (aa * aa * b / 16.0) - (3 * aa * aa * aa * aa / 256.0);
			//{Define coefficients of cubic equation}
			a[3][2] = q / 2.0;
			a[3][1] = (q * q - 4 * s) / 16.0;
			a[3][0] = -(rr * rr / 64.0);
			//{Calculate a real root of this equation}
			if ((rr != 0) || (a[3][1] >= 0)) goto Label_100;
			//{Particular case when this equation is of 2nd order}
			a[2][1] = a[3][2]; a[2][0] = a[3][1];
			Root_2();
			//{One takes the positive root}
			r[3] = r[1];
			goto Label_200;
			//{Calling Root_3 with sw=-1 to calculate one root only}
	Label_100:
			sw = -1;
			Root_3();
			//{real root of above equation}
	Label_200:
			k = sqrt(r[3]);
			//{Calculate L and M if k=0}
			if (k == 0)
			{
				rr = sqrt(q * q - 4 * s); goto Label_300;
			}
			q = q + (4 * r[3]);
			rr = rr / (2 * k);
	Label_300:
			l = (q - rr) / 2.0; m = (q + rr) / 2.0;
			//{Solving two equations of degree 2}
			a[2][1] = 2 * k; a[2][0] = l;
			//{1st equation}
			Root_2();
			//{Transferring solutions in r[3], r[4], ii}
			r[3] = r[1] - (a[4][3] / 4); r[4] = r[2] - (a[4][3] / 4.0); ii = im;
			a[2][1] = -2 * k; a[2][0] = m;
			//{2nd equation}
			Root_2();
			r[2] = r[2] - (a[4][3] / 4.0); r[1] = r[1] - (a[4][3] / 4.0);
		}
	}
	//*****************************************************************************

	void Root_432(int order, double Coeffs[5][5], double RealRoots[5], double *ImRoot1, double *ImRoot2)
	{
	//{initialize roots vector}
		RealRoots[0] = 0.0;
		RealRoots[1] = 0.0;
		RealRoots[2] = 0.0;
		RealRoots[3] = 0.0;
		RealRoots[4] = 0.0;

		n = order;
		if( n < 2  ) n = 2;
		if( n > 4  ) n = 4;
		
		for (i=0;i<5;i++)
			for (j=0;j<5;j++)
				a[i][j] = Coeffs[i][j];
			
	  //{calling root search main subroutine}
		Root_4();
	  
	  //{passing results}
		for (i=0;i<5;i++)
			RealRoots[i] = r[i];
		
		*ImRoot1 = im;
		*ImRoot2 = ii;
	}
	//*****************************************************************************
};

void Root_432(int order, double Coeffs[5][5], double RealRoots[5], double *ImRoot1, double *ImRoot2)
{
	CR432 calc;
	calc.Root_432( order, Coeffs, RealRoots, ImRoot1, ImRoot2 );
}

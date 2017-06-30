#include <math.h>
#include "types.h"
#include "procs.h"

#define RANGEN myrng
#define sqr(x) (x*x)

void SurfaceNormalErrors( MTRand &myrng, double CosIn[3],
						 TOpticalProperties *OptProperties,
						 double CosOut[3] ) throw(nanexcept)
{

/*{Purpose:  To add error terms to the surface normal vector at the surface in question

           Input - Seed    = Seed for RNG
                   CosIn   = Direction cosine vector of surface normal to which errors will be applied.
                   Element = Element data record
                   DFXYZ   = surface normal vector at interaction point

           Output - CosOut  = Output direction cosine vector of surface normal after error terms have been included
                   }*/

	int i=0;
	double Origin[3] = { 0.0, 0.0, 0.0 },
		Euler[3] = { 0.0, 0.0, 0.0 };
	double PosIn[3] = { 0.0, 0.0, 0.0 },
		PosOut[3] = { 0.0, 0.0, 0.0 };
	char dist = ' ';
	double delop = 0.0, delop3 = 0.0, thetax = 0.0,
		thetay = 0.0, ttheta = 0.0, theta2 = 0.0,
		phi = 0.0, theta = 0.0;
	double RRefToLoc[3][3] = { {0.0, 0.0, 0.0},
							   {0.0, 0.0, 0.0},
							   {0.0, 0.0, 0.0} };
	double RLocToRef[3][3] = { {0.0, 0.0, 0.0},
							   {0.0, 0.0, 0.0},
							   {0.0, 0.0, 0.0} };

	if ( CosIn[2] == 0.0 )
	{
		if ( CosIn[0] == 0.0 )
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

	CalculateTransformMatrices( Euler, RRefToLoc, RLocToRef );

	dist = OptProperties->DistributionType;
	delop = OptProperties->RMSSlopeError/1000.0;


	int nninner = 0;
	switch( dist )
	{
	case 'g':
	case 'G':
		//gaussian distribution
		delop3 = 3.0*delop;

		do
		{
			nninner = 0;

			do
			{
				thetax = 2.0*delop3*RANGEN() - delop3;
				if (delop == 0.0)       //handles case were specific optical errors are set to zero  06-11-07
					ttheta = 1.0;
				else
					ttheta = 1.0/exp(thetax*thetax/(2.0*delop*delop));

			} while( RANGEN() > ttheta );

			do
			{
				thetay = 2.0*delop3*RANGEN() - delop3;
				if (delop == 0.0)      //handles case were specific optical errors are set to zero  06-11-07
					ttheta = 1.0;
				else
					ttheta = 1.0/exp(thetay*thetay/(2.0*delop*delop));

			} while ( RANGEN() > ttheta );

			theta2 = thetax*thetax + thetay*thetay;

		} while ( theta2 > (delop3*delop3) );

		break;

	case 'p':
	case 'P':
		//pillbox distribution
		do
		{
			thetax = 2.0*delop*RANGEN() - delop;
			thetay = 2.0*delop*RANGEN() - delop;
			theta2 = thetax*thetax + thetay*thetay;
		}
		while ( theta2 > (delop*delop) );

		break;
	}

    /* {Transform to local coordinate system of ray to set up rotation matrices for coord and inverse
       transforms} */

	TransformToLocal(PosIn, CosIn, Origin, RRefToLoc, PosOut, CosOut);

	/* {Generate errors in terms of direction cosines in local ray coordinate system} */
	theta = sqrt(theta2);

	//phi = atan2(thetay, thetax); //This function appears to  present irregularities that bias results incorrectly for small values of thetay or thetax
	phi = RANGEN()*2.0*3.1415926535897932385; // Therefore have chosen to randomize phi rather than calculate from randomized theta components
                                                      //  obtained from the distribution. The two approaches are equivalent save for this issue with
                                                      //  arctan2.      wendelin 01-12-11 

	CosOut[0] = sin(theta)*cos(phi);
	CosOut[1] = sin(theta)*sin(phi);
	CosOut[2] = cos(theta);

	for (i=0;i<3;i++)
	{
		PosIn[i]=PosOut[i];
		CosIn[i]=CosOut[i];
	}

	/*{Transform perturbed ray back to element system}*/
	TransformToReference(PosIn, CosIn, Origin, RLocToRef, PosOut, CosOut);
}
//End of Procedure--------------------------------------------------------------

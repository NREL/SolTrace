#include <math.h>

#include "types.h"
#include "procs.h"

#define RANGEN myrng

void GenerateRay(
			MTRand &myrng,
			double PosSunStage[3],
			double Origin[3],
			double RLocToRef[3][3],
			TSun *Sun,
			double PosRayGlobal[3],
			double CosRayGlobal[3],
            double PosRaySun[3]
            )
{
/*{This procedure generates a randomly located ray in the x-y plane of the sun coordinate system in
 the z direction of the sun coord. system, checks to see that the ray is within the region of interest
 defined by the spatial extent of the elements of Stage as seen from the sun
 and ultimately transforms that ray to the global coord. system.   The z-axis of the sun coord. system points
 towards the Stage coord. system origin.

 Input
       - Seed = Seed for random number generator
       - Sun = Sun data record of type TSun
       - Origin = Primary Stage origin
       - RLocToRef = transformation matrix from local to reference frame
 Output
       - PosRayGlobal = Position of ray in Global coordinate system
       - CosRayGlobal = Direction cosines of ray in Global coordinate system} */

	double XRaySun = 0.0, YRaySun = 0.0, ZRaySun = 0.0;
	double CosRaySun[3] = { 0.0, 0.0, 0.0 };
	double PosRayStage[3] = { 0.0, 0.0, 0.0 };
	double CosRayStage[3] = { 0.0, 0.0, 0.0 };
	int NegPosSign = 0;
    PosRaySun[0] = 0.;
    PosRaySun[1] = 0.;
    PosRaySun[2] = 0.;

      //ZRaySun := 0.0;  //Origin of rays in xy plane of sun coord system.
	ZRaySun = -10000.0;  //changed 5/1/00.  rays originate from well bebind the sun coordinate system xy
                            // plane which has been translated to primary stage origin.         This value has been reduced signficantly because of numerical issues in tracing rays from sun
                            // to the closer form solution for a cylinder.  It used to 1e6 and has been reduced to 1e4, which should still be sufficient.   10-26-09 Wendelin


  //{Generate random rays inside of region of interest or from point source}
	
        if (Sun->PointSource) //fixed this on 3-18-13
	{
		PosRayGlobal[0] = Sun->Origin[0];
		PosRayGlobal[1] = Sun->Origin[1];
		PosRayGlobal[2] = Sun->Origin[2];
		
		if (RANGEN() <= 0.5)
			NegPosSign = -1;
		else
			NegPosSign = 1;
			
		CosRayGlobal[0] = NegPosSign*RANGEN();     //random direction for x part of ray vector
		
		if (RANGEN() <= 0.5)
			NegPosSign = -1;
		else
			NegPosSign = 1;
			
		CosRayGlobal[1] = NegPosSign*RANGEN();    //random direction for y part of ray vector
		
                if (RANGEN() <= 0.5)
                        NegPosSign = -1;
                else
                        NegPosSign = 1;

                CosRayGlobal[2] = NegPosSign*RANGEN();   //random direction for z part of ray vector

                double CosRayGMag = sqrt(CosRayGlobal[0]*CosRayGlobal[0]+CosRayGlobal[1]*CosRayGlobal[1]+CosRayGlobal[2]*CosRayGlobal[2]);

                CosRayGlobal[0] = CosRayGlobal[0]/CosRayGMag;  // obtain unit vector by dividing by magnitude
                CosRayGlobal[1] = CosRayGlobal[1]/CosRayGMag;
                CosRayGlobal[2] = CosRayGlobal[2]/CosRayGMag;
        }
	else
	{
        //following changed on 09/26/05 to more efficiently generate rays relative to element center of mass in primary stage
        /*{XRaySun := 2.0*MaxRad*ran3(Seed) - MaxRad;  //ran3 produces results independent of platform.
        YRaySun := 2.0*MaxRad*ran3(Seed) - MaxRad;
        if (XRaySun*XRaySun + YRaySun*YRaySun) > MaxRad*MaxRad then goto GENRAY;
        XRaySun := Xcm + XRaySun;  //adjust location of generated rays about element center of mass
        YRaySun := Ycm + YRaySun;}*/

		XRaySun = Sun->MinXSun + (Sun->MaxXSun - Sun->MinXSun)*RANGEN();     //uses a rectangular region of interest about the primary
		YRaySun = Sun->MinYSun + (Sun->MaxYSun - Sun->MinYSun)*RANGEN();     //stage. Added 09/26/05
		
		
		//{Offload ray location and direction cosines into sun array}
		PosRaySun[0] = XRaySun;
		PosRaySun[1] = YRaySun;
		PosRaySun[2] = ZRaySun;
		CosRaySun[0] = 0.0;
		CosRaySun[1] = 0.0;
		CosRaySun[2] = 1.0;

		//{Transform ray locations and dir cosines into Stage system}
		TransformToReference(PosRaySun, CosRaySun, PosSunStage, Sun->RLocToRef, PosRayStage, CosRayStage);
		
		//{Transform ray locations and dir cosines into global system}
		TransformToReference(PosRayStage, CosRayStage, Origin, RLocToRef, PosRayGlobal, CosRayGlobal);
	}
}
//End of Procedure--------------------------------------------------------------


bool LoadExistingStage0Ray(
    int index,
    std::vector<std::vector< double > > *raydat,
	double PosRayGlobal[3],
    double CosRayGlobal[3],
    st_uint_t &ElementNum,
    st_uint_t &RayNum
    )
{
    /* 
    Load an existing saved ray entering stage 0
    */

    std::vector<double> *theray = &raydat->at(index);
    
    if(theray->size() != 8)
        return false;

    for(int i=0; i<3; i++)
    {
        PosRayGlobal[i] = theray->at(i);
        CosRayGlobal[i] = theray->at(i+3);
    }
    ElementNum = (st_uint_t)theray->at(6);
    RayNum = (st_uint_t)theray->at(7);

    return true;
}

bool LoadExistingStage1Ray(
    int index,
    std::vector<std::vector< double > > *raydat,
	double PosRayGlobal[3],
    double CosRayGlobal[3],
    int &raynum
    )
{
    /* 
    Load an existing saved ray leaving stage 1
    */

    std::vector<double> *theray = &raydat->at(index);
    
    if(theray->size() != 7)
        return false;

    for(int i=0; i<3; i++)
    {
        PosRayGlobal[i] = theray->at(i);
        CosRayGlobal[i] = theray->at(i+3);
    }
    raynum = (int)theray->at(6);

    return true;
}


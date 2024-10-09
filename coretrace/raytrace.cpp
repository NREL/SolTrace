
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



#include <vector>
#include <string>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <ctime>
//#define WITH_DEBUG_TIMER
#ifdef WITH_DEBUG_TIMER
    #include <chrono>    //comment out for production
#endif

#include "types.h"
#include "procs.h"
#include "treemesh.h"


void time(const char *message, ofstream *fout)
{
#ifdef WITH_DEBUG_TIMER
    (*fout) << message << chrono::duration_cast< chrono::milliseconds >( chrono::system_clock::now().time_since_epoch() ).count() << "\n"; 
#endif
}

inline void CopyVec3( double dest[3], const std::vector<double> &src )
{
	dest[0] = src[0];
	dest[1] = src[1];
	dest[2] = src[2];
}

inline void CopyVec3( std::vector<double> &dest, double src[3] )
{
	dest[0] = src[0];
	dest[1] = src[1];
	dest[2] = src[2];
}

inline void CopyVec3( double dest[3], double src[3] )
{
	dest[0] = src[0];
	dest[1] = src[1];
	dest[2] = src[2];
}

#define ZeroVec(x) x[0]=x[1]=x[2]=0.0

class GlobalRay
{
public:
	GlobalRay() {
		Num = 0;
		for (int i=0;i<3;i++) Pos[i]=Cos[i]=0.0;
	}

	double Pos[3];
	double Cos[3];
	st_uint_t Num;
};

//structure to store element address and projected polar coordinate size
struct eprojdat
{
    TElement* el_addr;
    double d_proj;
    double az;
    double zen;

    eprojdat(){};
    eprojdat(TElement* e, double d, double a, double z)
    {
        el_addr = e;
        d_proj = d;
        az = a;
        zen = z;
    };
};

//Comparison function for sorting vector of eprojdat
static bool eprojdat_compare(const eprojdat &A, const eprojdat &B)
{
    return A.d_proj > B.d_proj;
};

bool Trace(TSystem *System, unsigned int seed,
		   st_uint_t NumberOfRays, 
		   st_uint_t MaxNumberOfRays,
		   bool IncludeSunShape, 
		   bool IncludeErrors,
           bool AsPowerTower,
		   int (*callback)(st_uint_t ntracedtotal, st_uint_t ntraced, st_uint_t ntotrace, st_uint_t curstage, st_uint_t nstages, void *data),
		   void *cbdata,
           std::vector< std::vector< double > > *st0data,
           std::vector< std::vector< double > > *st1in,
           bool save_st_data)
{
    
    bool PT_override = false;        //override speed improvements (use as compiled option for benchmarking old version)
    
    //don't try to use the element filtering method if: 
    if( System->StageList.size() > 0 
		&& (System->StageList[0]->ElementList.size() < 10    //the first stage contains only a few elements
			|| System->StageList.size() == 1)                //there's only one stage
      )
        PT_override = true;         

    bool load_st_data = st0data != 0 && st1in != 0;
    if(load_st_data)
    {
        load_st_data = st0data->size() > 0 && st1in->size() > 0;
    }
	bool StageHit = false;
	st_uint_t LastElementNumber = 0, LastRayNumber = 0;
	st_uint_t MultipleHitCount = 0;
	double LastPathLength = 0.0, PathLength = 0.0;

	double PosSunStage[3] = { 0.0, 0.0, 0.0 };
	double PosRayOutElement[3] = { 0.0, 0.0, 0.0 };
	double CosRayOutElement[3] = { 0.0, 0.0, 0.0 };
	double CosIn[3] = { 0.0, 0.0, 0.0 };
	double CosOut[3] = { 0.0, 0.0, 0.0 };
	double PosRayGlob[3] = { 0.0, 0.0, 0.0 };
	double CosRayGlob[3] = { 0.0, 0.0, 0.0 };
	double PosRayStage[3] = { 0.0, 0.0, 0.0 };
	double CosRayStage[3] = { 0.0, 0.0, 0.0 };
	double PosRayElement[3] = { 0.0, 0.0, 0.0 };
	double CosRayElement[3] = { 0.0, 0.0, 0.0 };
	double PosRaySurfElement[3] = { 0.0, 0.0, 0.0 };
	double CosRaySurfElement[3] = { 0.0, 0.0, 0.0 };
	double LastPosRaySurfElement[3] = { 0.0, 0.0, 0.0 };
	double LastCosRaySurfElement[3] = { 0.0, 0.0, 0.0 };
	double LastPosRaySurfStage[3] = { 0.0, 0.0, 0.0 };
	double LastCosRaySurfStage[3] = { 0.0, 0.0, 0.0 };
	double PosRaySurfStage[3] = { 0.0, 0.0, 0.0 };
	double CosRaySurfStage[3] = { 0.0, 0.0, 0.0 };
	double DFXYZ[3] = { 0.0, 0.0, 0.0 };
	double LastDFXYZ[3] = { 0.0, 0.0, 0.0 };
	double IncidentAngle = 0.0;
	double UnitLastDFXYZ[3] = { 0.0, 0.0, 0.0 };
	int ErrorFlag = 0, InterceptFlag = 0, HitBackSide = 0, LastHitBackSide = 0;

	std::vector<GlobalRay> IncomingRays;
	st_uint_t StageDataArrayIndex=0;
	bool PreviousStageHasRays = false;
	st_uint_t PreviousStageDataArrayIndex = 0;
	st_uint_t LastRayNumberInPreviousStage = NumberOfRays;

	ZeroVec( LastPosRaySurfStage );
	ZeroVec( LastCosRaySurfStage );

	//bool aspowertower_ok = false;
	bool in_multi_hit_loop = false;


	try
	{
		TOpticalProperties *optics=NULL;
		PreviousStageHasRays = false;

		int k = 0;
		TElement *optelm = 0;
		TRayData::ray_t *p_ray = 0;
		TStage *Stage;

		System->SunRayCount=0;
		st_uint_t RayNumber = 1;
		MTRand myrng(seed);
		st_uint_t RaysTracedTotal = 0;

        if( load_st_data && st0data->size() < 1)
        {
            System->errlog("empty stage 0 data array provided to Trace()");
            return false;
        }
        if( load_st_data && st1in->size() < 1 )
        {
            System->errlog("empty stage 1 data array provided to Trace()");
            return false;
        }

		if (NumberOfRays < 1)
		{
			System->errlog("invalid number of rays: %d", NumberOfRays);
			return false;
		}

		if (System->StageList.size() < 1)
		{
			System->errlog("no stages defined.");
			return false;
		}

		try
		{
			IncomingRays.resize( NumberOfRays );
		} catch (std::exception &e) {
			System->errlog("Incoming rays resize exception: %d, '%s'", NumberOfRays, e.what());
			return false;
		}


		if (!SunToPrimaryStage(System, System->StageList[0], &System->Sun, PosSunStage))
			return false;

#ifdef WITH_DEBUG_TIMER
        ofstream fout("C:\\Users\\mwagner\\Documents\\NREL\\Dev\\SolTraceWX\\log.txt");
        fout.clear();
#else 
        ofstream fout;
#endif

        time("Initialize:\t", &fout);

        /* 
        Calculate hash tree for sun incoming plane.
        
        Calculate hash tree for reflection to receiver plane (polar coordinates).
        */
        st_hash_tree sun_hash;
        st_hash_tree rec_hash;
        double reccm_helio[3];  //receiver centroid in heliostat field coordinates
        if(! PT_override )
        {
            //Calculate the center of mass of the receiver stage (StageList[1]) in heliostat stage coordinates.
            double reccm[] = {0., 0., 0.};
            int nelrec=0;
            if(AsPowerTower)
            {
                for(st_uint_t j=0; j<System->StageList[1]->ElementList.size(); j++)
                {
                    TElement* el = System->StageList[1]->ElementList.at(j);

                    if(! el->Enabled)
                        continue;
                
                    nelrec++;

                    for(int jj=0; jj<3; jj++)
                        reccm[jj] += el->Origin[jj]; 
                }
                for(int jj=0; jj<3; jj++)
                    reccm[jj] /= (double)nelrec;    //average
            

                //Transform to reference 
                double dum1[] = {0., 0., 1.};
                double dum2[3];
                double reccm_global[3];
                TransformToReference(reccm, dum1, System->StageList[1]->Origin, System->StageList[1]->RLocToRef, reccm_global, dum2);

                //Transform to local (heliostat). reccm_helio is the x,y,z position of the receiver centroid in heliostat stage coordinates.
                TransformToLocal(reccm_global, dum1, System->StageList[0]->Origin, System->StageList[0]->RRefToLoc, reccm_helio, dum2);
            }
            //Create an array that stores the element address and the projected size in polar coordinates
            vector<eprojdat> el_proj_dat;
            el_proj_dat.reserve( System->StageList[0]->ElementList.size() );

            //calculate the smallest zone size. This should be on the order of the largest element in the stage. 
            //load stage 0 elements into the mesh
            double d_elm_max = -9.e9;

            time("Calculating element sizes:\t", &fout);

            for( st_uint_t i=0; i<System->StageList[0]->ElementList.size(); i++)
            {
                TElement* el = System->StageList[0]->ElementList.at(i);

                el->element_number = i+1;   //use index for element number

                double d_elm;

                switch (el->ShapeIndex)
                {
                //circular aperture
                case 'c':
                case 'C':
                //hexagonal aperture
                case 'h':
                case 'H':
                //triangular aperture
                case 't':
                case 'T':
                    d_elm =  el->ParameterA;
                    break;
                //rectangular aperture
                case 'r':
                case 'R':
                    d_elm =  sqrt(el->ParameterA*el->ParameterA + el->ParameterB*el->ParameterB);
                    break;
                //annular aperture
                case 'a':
                case 'A':
                    d_elm =  el->ParameterB;
                    break;
                case 'l':
                case 'L': 
                    //off axis aperture section of line focus trough  or cylinder
                    d_elm =  sqrt(el->ParameterB*el->ParameterB*4. + el->ParameterC*el->ParameterC);
                    break;
                //Irregular triangle
                case 'i':
                case 'I':
                //irregular quadrilateral
                case 'q':
                case 'Q':
                {
                    double xmax = fmax( el->ParameterA, fmax( el->ParameterC, el->ParameterE ) );
                    double xmin = fmin( el->ParameterA, fmin( el->ParameterC, el->ParameterE ) );
                    double ymax = fmax( el->ParameterB, fmax( el->ParameterD, el->ParameterF ) );
                    double ymin = fmin( el->ParameterB, fmin( el->ParameterD, el->ParameterF ) );

                    if( el->ShapeIndex == 'q' || el->ShapeIndex == 'Q' )
                    {
                        xmax = fmax(xmax, el->ParameterG);
                        xmin = fmin(xmin, el->ParameterG);
                        ymax = fmax(ymax, el->ParameterH);
                        ymin = fmin(ymin, el->ParameterH);
                    }

                    double dx = xmax - xmin;
                    double dy = ymax - ymin; 

                    d_elm =  sqrt(dx*dx + dy*dy);
                
                    break;
                }
                default:
                    break;
                }

                d_elm_max = fmax(d_elm_max, d_elm);
                
                if(AsPowerTower)
                {
                    //Calculate the distance from the receiver to the element and the max projected size
                    double dX[3];
                    for(int jj=0; jj<3; jj++)
                        dX[jj] = el->Origin[jj] - reccm_helio[jj];  //vector from receiver to heliostat (not unitized)
                    double r_elm = 0.;
                    for(int jj=0; jj<3; jj++)
                        r_elm += dX[jj]*dX[jj];     
                    r_elm = sqrt(r_elm);            //vector length
                    double d_elm_proj = d_elm / r_elm;  //Projected size of the element from the view of the receiver (radians)
                
                    //calculate az,zen coordinate
                    double az,zen;
                    az = atan2(dX[0]/r_elm, dX[1]/r_elm);       //Az coordinate of the heliostat from the receiver's perspective
                    zen = asin(dX[2]/r_elm);                    //Zen coordinate """"

                    el_proj_dat.push_back( eprojdat(el, d_elm_proj, az, zen) );
                }
            }

            if(AsPowerTower)
            {
                time("Sorting polar mesh entries:\t", &fout);

                //Sort the polar projections by size, largest to smallest
                std::sort(el_proj_dat.begin(), el_proj_dat.end(), eprojdat_compare);
            }

            //set up the layout data object that provides configuration details for the hash tree
            KDLayoutData sun_ld;
            sun_ld.xlim[0] = System->Sun.MinXSun;
            sun_ld.xlim[1] = System->Sun.MaxXSun;
            sun_ld.ylim[0] = System->Sun.MinYSun;
            sun_ld.ylim[1] = System->Sun.MaxYSun;
            sun_ld.min_unit_dx = d_elm_max;
            sun_ld.min_unit_dy = d_elm_max;

            sun_hash.create_mesh( sun_ld );
            time("Adding solar mesh elements:\t", &fout);
            
           //load stage 0 elements into the mesh
            for( st_uint_t i=0; i<System->StageList[0]->ElementList.size(); i++)
            {
                TElement* el = System->StageList[0]->ElementList.at(i);
                sun_hash.add_object( (void*)el, el->PosSunCoords[0], el->PosSunCoords[1] );
            }

            //calculate and associate neighbors with each zone
            time("Adding solar mesh neighbors:\t", &fout);
            sun_hash.add_neighborhood_data();

            if(AsPowerTower)
            {
                //Set things up for the polar coordinate tree
                KDLayoutData rec_ld;
                rec_ld.xlim[0] = -M_PI;
                rec_ld.xlim[1] = M_PI;
                rec_ld.ylim[0] = -M_PI/2.;
                rec_ld.ylim[1] = M_PI/2.;
                //use smallest element to set the minimum size
                rec_ld.min_unit_dx = rec_ld.min_unit_dy = el_proj_dat.back().d_proj; //radians at equator
            
                rec_hash.create_mesh( rec_ld );
                time("Adding polar mesh elements:\t", &fout);

                //load stage 0 elements into the receiver mesh in the order of largest projection to smallest
                for( int i=0; i<el_proj_dat.size(); i++)
                {
                    eprojdat* D = &el_proj_dat.at(i);

                    //Calculate the angular span of the element
                    double angspan[2];
                    double adjmult = 1.5;
                    angspan[0] = D->d_proj/cos(fabs(D->zen))*adjmult;   //azimuthal span
                    angspan[0] = fmin(angspan[0], 2.*M_PI);     //limit to circumference 
                    angspan[1] = D->d_proj/M_PI*adjmult;    //zenithal span
                    rec_hash.add_object( (void*)D->el_addr,  D->az, D->zen, angspan);     
                }
                time("Adding polar mesh neighbors:\t",&fout);
                //associate neighbors with each zone
                rec_hash.add_neighborhood_data(); 
            }
        }

//#define WRITE_NODE_FILE
#ifdef WRITE_NODE_FILE
        {
        //Write out to a file for debugging
        ofstream fout2("C:\\Users\\mwagner\\Documents\\NREL\\Dev\\SolTraceWX\\meshxy.txt");
        fout2.clear();

        fout2 << "node,xmin,xmax,ymin,ymax\n";

        vector<st_opt_element>* all_nodes = sun_hash.get_all_nodes();
        for(int i=0; i<all_nodes->size(); i++)
        {
            double* xr = all_nodes->at(i).get_xr();
            double* yr = all_nodes->at(i).get_yr();

            fout2 << all_nodes->at(i).get_address() << "," << xr[0] << "," << xr[1] << "," << yr[0] << "," << yr[1] << "\n";
        }

        fout2.close();

        ofstream fout3("C:\\Users\\mwagner\\Documents\\NREL\\Dev\\SolTraceWX\\meshpolar.txt");
        fout3.clear();

        fout3 << "node,xmin,xmax,ymin,ymax\n";

        all_nodes = rec_hash.get_all_nodes();
        for(int i=0; i<all_nodes->size(); i++)
        {
            double* xr = all_nodes->at(i).get_xr();
            double* yr = all_nodes->at(i).get_yr();

            fout3 << all_nodes->at(i).get_address() << "," << xr[0] << "," << xr[1] << "," << yr[0] << "," << yr[1] << "\n";
        }

        fout3.close();
        }
#endif


        //declare items used within the loop
        vector<void*> sunint_elements;
        vector<void*> reflint_elements;
        bool has_elements;

        time("Starting stage calculations:\t", &fout);
#ifdef WITH_DEBUG_TIMER
        fout.close();
#endif

        //use the callbacks based on elapsed time rather than fixed rays processed. 

        clock_t startTime = clock();     //start timer
        int rays_per_callback_estimate = 50;    //starting rough estimate for how often to check the clock

		for (st_uint_t i=0;i<System->StageList.size();i++)
		{

			if (i > 0 && PreviousStageHasRays == false)
			{
				// no rays to pass through from previous stage
				// so nothing to trace in this stage
				goto Label_EndStageLoop;
			}
            
			Stage = System->StageList[i];

			LastElementNumber = 0;
			LastRayNumber = 0;
			LastHitBackSide = 0;

			StageDataArrayIndex = 0;
			PreviousStageDataArrayIndex = 0;


            //if loading stage 0 data, construct appropriate arrays here
            if(i==0 && load_st_data)
            {
                double rpos[3],rcos[3];
                //Stage 0 data
                for(int j=0; j<st0data->size(); j++)   
                {
                    
                    LoadExistingStage0Ray(j, st0data, 
                        rpos, rcos,
                        LastElementNumber, LastRayNumber);


                    p_ray = Stage->RayData.Append( 
                        rpos, rcos,
						LastElementNumber, 1,
						LastRayNumber );

                }

                //stage 1 data
                for(int j=0; j<st1in->size(); j++)
                {
                    int rnum;
                    LoadExistingStage1Ray(j, st1in, rpos, rcos, rnum);
                    CopyVec3(IncomingRays[j].Pos, rpos);
                    CopyVec3(IncomingRays[j].Cos, rcos);
                    IncomingRays[j].Num = rnum;
                }

                PreviousStageHasRays = true;
                PreviousStageDataArrayIndex = st1in->size()-1;
                System->SunRayCount = LastRayNumber;
                goto Label_EndStageLoop;
            }
            

            
Label_StartRayLoop:
			MultipleHitCount = 0;
            sunint_elements.clear();

            has_elements = true;
			if ( i == 0 )
			{


                // we are in the first stage, so 
				// generate a new sun ray in global coords
                double PosRaySun[3];
				GenerateRay(myrng, PosSunStage, Stage->Origin,
							Stage->RLocToRef, &System->Sun,
							PosRayGlob, CosRayGlob, PosRaySun);
				    System->SunRayCount++;


				if (System->SunRayCount > MaxNumberOfRays)
				{
					System->errlog("generated sun rays reached maximum count: %d", MaxNumberOfRays);
					return false;
				}

                /* 
                Find the list of elements that could potentially interact with this ray. If empty, continue
                */
                if(! PT_override) //AsPowerTower)
                    has_elements = sun_hash.get_all_data_at_loc( sunint_elements, PosRaySun[0], PosRaySun[1] );

			}
			else
			{
				// we are in a subsequent stage, so trace using an incoming ray
				// saved from the previous stages
                RayNumber = IncomingRays[StageDataArrayIndex].Num;
				CopyVec3( PosRayGlob, IncomingRays[StageDataArrayIndex].Pos );
				CopyVec3( CosRayGlob, IncomingRays[StageDataArrayIndex].Cos );
				StageDataArrayIndex++;
				
			}

			// transform the global incoming ray to local stage coordinates
			TransformToLocal(PosRayGlob, CosRayGlob, 
				Stage->Origin, Stage->RRefToLoc, 
				PosRayStage, CosRayStage);


			// CheckForCancelAndUpdateProgressBar
			if (callback != 0
				&& RaysTracedTotal++ % rays_per_callback_estimate == 0)
			{
                if( RaysTracedTotal > 1 )
                {
                    //update how often to call this
                    double msec_per_ray = 1000.*( clock() - startTime ) / CLOCKS_PER_SEC / (double)(RaysTracedTotal > 0 ? RaysTracedTotal : 1);
                    //set the new callback estimate to be about 50 ms
                    rays_per_callback_estimate = (int)( 200. / msec_per_ray );
                    //limit to something reasonable
                    rays_per_callback_estimate = rays_per_callback_estimate < 5 ? 5 : rays_per_callback_estimate;
                }

                //do the callback
				if ( ! (*callback)( RaysTracedTotal, RayNumber,
									LastRayNumberInPreviousStage, i+1,
									System->StageList.size(), cbdata ))
					return true;
			}
            
            in_multi_hit_loop = false;
            
Label_MultiHitLoop:
			LastPathLength = 1e99;
			StageHit = false;

            st_uint_t nintelements;
            if( i==0 && !PT_override)
            {
                if( in_multi_hit_loop )
                {
                    if( AsPowerTower )
                    {
                        //>=Second time through - checking for first stage multiple element interactions
                    
                        //get ray position in receiver polar coordinates
                        double raypvec[3];
                        for(int jj=0; jj<3; jj++)
                            raypvec[jj] = PosRayStage[jj] - reccm_helio[jj];
                        double raypvecmag = sqrt(raypvec[0]*raypvec[0] + raypvec[1]*raypvec[1] + raypvec[2]*raypvec[2]);
                        double raypol[2];
                        raypol[0] = atan2(raypvec[0], raypvec[1]);
                        raypol[1] = asin(raypvec[2]/raypvecmag);
                        //get elements in the vicinity of the ray's polar coordinates
                        reflint_elements.clear();
                        rec_hash.get_all_data_at_loc( reflint_elements, raypol[0], raypol[1]);
                        nintelements = reflint_elements.size();
                        has_elements = nintelements > 0;

                    }
                    else
                    {
                        nintelements = Stage->ElementList.size();
                    }
                }
                else
                {
                    //First time through - checking for sun ray intersections
                    if( has_elements )
                        nintelements = sunint_elements.size();
                    else
                        nintelements = 0;
                }
            }
            else
                nintelements = Stage->ElementList.size();

            for( st_uint_t j=0; j<nintelements; j++)
			{
                TElement *Element; // = Stage->ElementList[j];
                if( i == 0 && !PT_override )
                {
                    if( in_multi_hit_loop )
                    {
                        if( AsPowerTower )
                            Element = (TElement*)reflint_elements.at(j);
                        else
                            Element = (TElement*)Stage->ElementList[j];
                    }
                    else
                        Element = (TElement*)sunint_elements.at(j);
                }
                else
				    Element = Stage->ElementList[j];

				if (!Element->Enabled)
					continue;

				//  {Transform ray to element[j] coord system of Stage[i]}
				TransformToLocal( PosRayStage, CosRayStage,
								  Element->Origin, Element->RRefToLoc,
								  PosRayElement, CosRayElement);

				ErrorFlag = 0;
				HitBackSide = 0;
				InterceptFlag = 0;

				// increment position by tiny amount to get off the element if tracing to the same element
				PosRayElement[0] = PosRayElement[0] + 1.0e-5*CosRayElement[0];
				PosRayElement[1] = PosRayElement[1] + 1.0e-5*CosRayElement[1];
				PosRayElement[2] = PosRayElement[2] + 1.0e-5*CosRayElement[2];

				// {Determine if ray intersects element[j]; if so, Find intersection point with surface of element[j] }
				DetermineElementIntersectionNew(Element, PosRayElement, CosRayElement,
					PosRaySurfElement, CosRaySurfElement, DFXYZ, 
					&PathLength, &ErrorFlag, &InterceptFlag, &HitBackSide);



				if (InterceptFlag)
				{
				  //{If hit multiple elements, this loop determines which one hit first.
				  //Also makes sure that correct part of closed surface is hit. Also, handles wavy, but close to flat zernikes and polynomials correctly.}
				  //if (PathLength < LastPathLength) and (PosRaySurfElement[2] <= Element->ZAperture) then
					if (PathLength < LastPathLength)
					{
						if (PosRaySurfElement[2] <= Element->ZAperture 
							|| Element->SurfaceIndex == 'm'
							|| Element->SurfaceIndex == 'M'
							|| Element->SurfaceIndex == 'r'
							|| Element->SurfaceIndex == 'R') 
						{
							StageHit = true;
							LastPathLength = PathLength;
							CopyVec3( LastPosRaySurfElement, PosRaySurfElement );
							CopyVec3( LastCosRaySurfElement, CosRaySurfElement );
							CopyVec3( LastDFXYZ, DFXYZ );
							LastElementNumber = ( i == 0 && !PT_override )? Element->element_number : j+1;    //mjw change from j index to element id
							LastRayNumber = RayNumber;
							TransformToReference(PosRaySurfElement, CosRaySurfElement, 
								Element->Origin, Element->RLocToRef, 
								PosRaySurfStage, CosRaySurfStage);

							CopyVec3( LastPosRaySurfStage, PosRaySurfStage );
							CopyVec3( LastCosRaySurfStage, CosRaySurfStage );
							LastHitBackSide = HitBackSide;
						}
					}
				}			
			}

			//  {Logic for ray which misses stage element - Note that all rays eventually satisfy this
			//  condition because rays are continually traced until they no longer hit the stage}
Label_StageHitLogic:

			if ( !StageHit )
			{
				if ( i == 0 ) // first stage only
				{
					if (MultipleHitCount == 0)
						goto Label_StartRayLoop; // ray misses 1st stage completely so get a new sun ray
					else
					{

						// at least one hit on stage, so move on to next ray
						CopyVec3( IncomingRays[PreviousStageDataArrayIndex].Pos, PosRayGlob );
						CopyVec3( IncomingRays[PreviousStageDataArrayIndex].Cos, CosRayGlob );
						IncomingRays[PreviousStageDataArrayIndex].Num = RayNumber;

						if (RayNumber == NumberOfRays)
							goto Label_EndStageLoop;

						PreviousStageDataArrayIndex++;
						PreviousStageHasRays = true;

						RayNumber++;
						goto Label_StartRayLoop;
					} 
				}
				else
				{
					// stages beyond first stage
					if (Stage->TraceThrough || MultipleHitCount > 0)
					{

						CopyVec3( IncomingRays[PreviousStageDataArrayIndex].Pos, PosRayGlob );
						CopyVec3( IncomingRays[PreviousStageDataArrayIndex].Cos, CosRayGlob );
						IncomingRays[PreviousStageDataArrayIndex].Num = RayNumber;

						if (RayNumber == LastRayNumberInPreviousStage)
							goto Label_EndStageLoop;

						PreviousStageDataArrayIndex++;
						PreviousStageHasRays = true;

						if (MultipleHitCount == 0)
							goto Label_FlagMiss;

						goto Label_StartRayLoop;
					}
Label_FlagMiss:
					LastElementNumber = 0;
					LastRayNumber = RayNumber;
					CopyVec3(LastPosRaySurfStage, PosRayStage);
					CopyVec3(LastCosRaySurfStage, CosRayStage);
				}
			} // end of not stagehit logic

			p_ray = Stage->RayData.Append( LastPosRaySurfStage,
								  LastCosRaySurfStage,
								  LastElementNumber,
								  i+1,
								  LastRayNumber );

			if (!p_ray)
			{
				System->errlog("Failed to save ray data at index %d", Stage->RayData.Count()-1);
				return false;
			}

			if (LastElementNumber == 0) // {If missed all elements}
			{
				if (RayNumber == LastRayNumberInPreviousStage)
				{
					if ( !Stage->TraceThrough )
					{
						PreviousStageHasRays = false;
						if (PreviousStageDataArrayIndex > 0)
						{
							PreviousStageHasRays = true;
							PreviousStageDataArrayIndex--; // last ray was previous one
						}
					}
					goto Label_EndStageLoop;
				}
				else
				{
					if (i == 0) RayNumber++; // generate new sun ray
					goto Label_StartRayLoop;
				}
			}

			MultipleHitCount++;

			if ( Stage->Virtual )
			{
				CopyVec3(PosRayOutElement, LastPosRaySurfElement);
				CopyVec3(CosRayOutElement, LastCosRaySurfElement);
				goto Label_TransformBackToGlobal;
			}

			// {Otherwise trace ray through interaction}
			// {Determine if backside or frontside properties should be used}
		
			// trace through the interaction
			optelm = Stage->ElementList[ p_ray->element - 1 ];
			optics = 0;
			
			if (LastHitBackSide)
				optics = &optelm->Optics->Back;
			else
				optics = &optelm->Optics->Front;


			double TestValue;
			switch(optelm->InteractionType )
			{
			case 1: // refraction
				if (optics->UseTransmissivityTable)
				{
					int npoints = optics->TransmissivityTable.size();
					int m = 0;
					UnitLastDFXYZ[0] = -LastDFXYZ[0] / sqrt(DOT(LastDFXYZ, LastDFXYZ));
					UnitLastDFXYZ[1] = -LastDFXYZ[1] / sqrt(DOT(LastDFXYZ, LastDFXYZ));
					UnitLastDFXYZ[2] = -LastDFXYZ[2] / sqrt(DOT(LastDFXYZ, LastDFXYZ));
					IncidentAngle = acos(DOT(LastCosRaySurfElement, UnitLastDFXYZ))*1000.;  //[mrad]
					if (IncidentAngle >= optics->TransmissivityTable[npoints - 1].angle)
					{
						TestValue = optics->TransmissivityTable[npoints - 1].trans;
					}
					else
					{
						while (optics->TransmissivityTable[m].angle < IncidentAngle)
							m++;

						if (m == 0)
							TestValue = optics->TransmissivityTable[m].trans;
						else
							TestValue = (optics->TransmissivityTable[m].trans + optics->TransmissivityTable[m - 1].trans) / 2.0;
					}
				}
				else
					TestValue = optics->Transmissivity; 
				break;
			case 2: // reflection

				if ( optics->UseReflectivityTable )
				{
					int npoints = optics->ReflectivityTable.size();
					int m = 0;
					UnitLastDFXYZ[0] = -LastDFXYZ[0]/sqrt(DOT(LastDFXYZ,LastDFXYZ));
					UnitLastDFXYZ[1] = -LastDFXYZ[1]/sqrt(DOT(LastDFXYZ,LastDFXYZ));
					UnitLastDFXYZ[2] = -LastDFXYZ[2]/sqrt(DOT(LastDFXYZ,LastDFXYZ));
					IncidentAngle = acos(DOT(LastCosRaySurfElement,UnitLastDFXYZ)) * 1000.;  //[mrad]
					if (IncidentAngle >= optics->ReflectivityTable[ npoints-1 ].angle )
					{
						TestValue = optics->ReflectivityTable[ npoints-1 ].refl;
					}
					else
					{
						while ( optics->ReflectivityTable[m].angle < IncidentAngle )
							m++;

						if (m == 0)
							TestValue = optics->ReflectivityTable[m].refl;
						else
							TestValue = (optics->ReflectivityTable[m].refl + optics->ReflectivityTable[m-1].refl)/2.0;
					}
				}
				else
					TestValue = optics->Reflectivity;
				break;
			default:
				System->errlog("Bad optical interaction type = %d (stage %d)",i,optelm->InteractionType);
				return false;
			}


		//  {Apply MonteCarlo probability of absorption. Limited for now, but can make more complex later on if desired}
			if (TestValue <= myrng())
			{
				// ray was fully absorbed, so indicate by negating the element number
				p_ray->element = 0 - p_ray->element;

				if (RayNumber == LastRayNumberInPreviousStage)
				{
					PreviousStageHasRays = false;
					if (PreviousStageDataArrayIndex > 0)
					{
						PreviousStageDataArrayIndex--;
						PreviousStageHasRays = true;
					}
					goto Label_EndStageLoop;
				}
				else
				{
					if (i == 0)
					{
						if (RayNumber == NumberOfRays)
							goto Label_EndStageLoop;
						else
							RayNumber++;
					}

					goto Label_StartRayLoop;
				}
			}

Label_TransformBackToGlobal:
			k = abs( p_ray->element ) - 1;

			if ( !Stage->Virtual )
			{
				if (IncludeSunShape && i == 0 && MultipleHitCount == 1)//change to account for first hit only in primary stage 8-11-31
				{
					// Apply sunshape to UNPERTURBED ray at intersection point
					//only apply sunshape error once for primary stage
					CopyVec3(CosIn, LastCosRaySurfElement);
					Errors(myrng, CosIn, 1, &System->Sun,
						   Stage->ElementList[k], optics, CosOut, LastDFXYZ);  //sun shape
					CopyVec3(LastCosRaySurfElement, CosOut);
				}

				//{Determine interaction at surface and direction of perturbed ray}
				ErrorFlag = 0;

				// {Apply surface normal errors to surface normal before interaction ray at intersection point - Wendelin 11-23-09}
				if( IncludeErrors )
				{
					CopyVec3( CosIn, CosRayOutElement );
					SurfaceNormalErrors(myrng, LastDFXYZ, optics, CosOut);  //surface normal errors
					CopyVec3( LastDFXYZ, CosOut );
				}

				Interaction( myrng, LastPosRaySurfElement, LastCosRaySurfElement, LastDFXYZ,
					Stage->ElementList[k]->InteractionType, optics, 630.0, 
					PosRayOutElement, CosRayOutElement, &ErrorFlag);

				// {Apply specularity optical error to PERTURBED (i.e. after interaction) ray at intersection point}
				if( IncludeErrors )
				{
					if (optics->DistributionType == 'F' || optics->DistributionType == 'f')
						CopyVec3(CosIn, LastDFXYZ);  // Apply diffuse errors relative to surface normal
					else
						CopyVec3(CosIn, CosRayOutElement); // Apply all other errors relative to the specularly-reflected direction

					Errors(myrng, CosIn, 2, &System->Sun,
						   Stage->ElementList[k], optics, CosOut, LastDFXYZ);  //optical errors
					CopyVec3(CosRayOutElement, CosOut);
				}
			}

			// { Transform ray back to stage coord system and trace through stage again}
			TransformToReference(PosRayOutElement, CosRayOutElement, 
					Stage->ElementList[k]->Origin, Stage->ElementList[k]->RLocToRef, 
					PosRayStage, CosRayStage);
			TransformToReference(PosRayStage, CosRayStage, 
					Stage->Origin, Stage->RLocToRef, 
					PosRayGlob, CosRayGlob);

			if (!Stage->MultiHitsPerRay)
			{
				StageHit = false;
				goto Label_StageHitLogic;
			}
			else
            {
                in_multi_hit_loop = true;
				goto Label_MultiHitLoop;
            }

Label_EndStageLoop:

			if(i==0 && save_st_data)
            {
                //if flagged save the stage 0 incoming rays data
                TRayData *raydat = &Stage->RayData;
                st_uint_t nray0 = raydat->Count();
        
                for(st_uint_t ii=0; ii<nray0; ii++)
                {
                    TRayData::ray_t *rr = raydat->Index(ii,false);

                    std::vector<double> ray(8);
                    for(int j=0; j<3; j++)
                        ray[j] = rr->pos[j];
                    for(int j=0; j<3; j++)
                        ray[j+3] = rr->cos[j];
                    ray[6] = rr->element;
                    ray[7] = rr->raynum;
                    st0data->push_back(ray);
                }
            }

            if(i==1 && save_st_data)
            {
                //if flagged, save the stage 1 incoming rays data to the data structure passed into the algorithm
                for(int ir=0; ir<StageDataArrayIndex; ir++)
                {
                    st1in->push_back(std::vector<double>(7));
                    for(int jr=0; jr<3; jr++)
                    {
                        st1in->back().at(jr) = IncomingRays[ir].Pos[jr];
                        st1in->back().at(jr+3) = IncomingRays[ir].Cos[jr];
                    }
                    st1in->back().at(6) = IncomingRays[ir].Num;
                }
            }
            

            if (!PreviousStageHasRays)
			{
				LastRayNumberInPreviousStage = 0;
				continue; // no rays to carry forward
			}

			if (PreviousStageDataArrayIndex < IncomingRays.size())
			{
				LastRayNumberInPreviousStage = IncomingRays[PreviousStageDataArrayIndex].Num;
				if (LastRayNumberInPreviousStage == 0)
				{
					size_t pp = IncomingRays[PreviousStageDataArrayIndex-1].Num;
					System->errlog("LastRayNumberInPreviousStage=0, stage %d, PrevIdx=%d, CurIdx=%d, pp=%d", i+1,
										PreviousStageDataArrayIndex, StageDataArrayIndex, pp);
					return false;
				}
			}
			else
			{
				System->errlog("Invalid PreviousStageDataArrayIndex: %u, @ stage %d",
							   PreviousStageDataArrayIndex, i+1);
				return false;
			}
		}

		return true;
	}
	catch( const std::exception &e )
	{
		System->errlog("trace error: %s", e.what());
		return false;
	}
}

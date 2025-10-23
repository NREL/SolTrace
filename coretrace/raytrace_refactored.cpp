
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
#include <cmath>
#include <algorithm>
#include <ctime>

#include "types.h"
#include "procs.h"
#include "treemesh.h"

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

class GlobalRay_refactored
{
public:
	GlobalRay_refactored() {
		Num = 0;
		for (int i=0;i<3;i++) Pos[i]=Cos[i]=0.0;
	}

	double Pos[3];
	double Cos[3];
	st_uint_t Num;
};

void FindElementHit(
	// stage info
	const int i, const TStage* Stage, const bool PT_override, const bool AsPowerTower,
	
	// element info
	const int nintelements, const vector<void*>& sunint_elements, const vector<void*>& reflint_elements,
	
	// ray info
	const int RayNumber, const bool in_multi_hit_loop,
	double(&PosRayStage)[3], double(&CosRayStage)[3],

	// outputs
	double(&LastPosRaySurfElement)[3], double(&LastCosRaySurfElement)[3], double(&LastDFXYZ)[3],
	st_uint_t& LastElementNumber, st_uint_t& LastRayNumber, 
	double(&LastPosRaySurfStage)[3], double(&LastCosRaySurfStage)[3],

	int& ErrorFlag, int& LastHitBackSide,
	bool& StageHit)
{
	// Initialize Variables
	double LastPathLength = 1e99;
	int HitBackSide = 0;
	int InterceptFlag = 0;
	double DFXYZ[3] = { 0.0, 0.0, 0.0 };
	double PosRayElement[3] = { 0.0, 0.0, 0.0 };
	double CosRayElement[3] = { 0.0, 0.0, 0.0 };
	double PosRaySurfStage[3] = { 0.0, 0.0, 0.0 };
	double CosRaySurfStage[3] = { 0.0, 0.0, 0.0 };
	double PosRaySurfElement[3] = { 0.0, 0.0, 0.0 };
	double CosRaySurfElement[3] = { 0.0, 0.0, 0.0 };
	StageHit = false;


	for (st_uint_t j = 0; j < nintelements; j++)
	{
		TElement* Element; // = Stage->ElementList[j];
		if (i == 0 && !PT_override)
		{
			if (in_multi_hit_loop)
			{
				if (AsPowerTower)
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
		TransformToLocal(PosRayStage, CosRayStage,
			Element->Origin, Element->RRefToLoc,
			PosRayElement, CosRayElement);

		ErrorFlag = 0;
		HitBackSide = 0;
		InterceptFlag = 0;
		double PathLength = 0;

		// increment position by tiny amount to get off the element if tracing to the same element
		PosRayElement[0] = PosRayElement[0] + 1.0e-5 * CosRayElement[0];
		PosRayElement[1] = PosRayElement[1] + 1.0e-5 * CosRayElement[1];
		PosRayElement[2] = PosRayElement[2] + 1.0e-5 * CosRayElement[2];

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
					CopyVec3(LastPosRaySurfElement, PosRaySurfElement);
					CopyVec3(LastCosRaySurfElement, CosRaySurfElement);
					CopyVec3(LastDFXYZ, DFXYZ);
					LastElementNumber = (i == 0 && !PT_override) ? Element->element_number : j + 1;    //mjw change from j index to element id
					LastRayNumber = RayNumber;
					TransformToReference(PosRaySurfElement, CosRaySurfElement,
						Element->Origin, Element->RLocToRef,
						PosRaySurfStage, CosRaySurfStage);

					CopyVec3(LastPosRaySurfStage, PosRaySurfStage);
					CopyVec3(LastCosRaySurfStage, CosRaySurfStage);
					LastHitBackSide = HitBackSide;
				}
			}
		}
	}
}

void ProcessInteraction(
	// system info
	TSystem* System, MTRand &myrng, const bool IncludeSunShape, TOpticalProperties* optics,
	const bool IncludeErrors,

	// stage info
	const int i, const TStage* Stage, const int k,
	
	// ray info
	const st_uint_t MultipleHitCount,
	double(&LastDFXYZ)[3],
	
	// Outputs
	double(&LastCosRaySurfElement)[3], int& ErrorFlag,
	double(&CosRayOutElement)[3], double(&LastPosRaySurfElement)[3],
	double(&PosRayOutElement)[3], int& myrng_counter)
{
	// Initialize
	double CosIn[3] = { 0.0, 0.0, 0.0 };
	double CosOut[3] = { 0.0, 0.0, 0.0 };

	if (!Stage->Virtual)
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
		if (IncludeErrors)
		{
			CopyVec3(CosIn, CosRayOutElement);
			SurfaceNormalErrors(myrng, LastDFXYZ, optics, CosOut);  //surface normal errors
			myrng_counter++;
			CopyVec3(LastDFXYZ, CosOut);
		}

		Interaction(myrng, LastPosRaySurfElement, LastCosRaySurfElement, LastDFXYZ,
			Stage->ElementList[k]->InteractionType, optics, 630.0,
			PosRayOutElement, CosRayOutElement, &ErrorFlag);
		myrng_counter++;

		// {Apply specularity optical error to PERTURBED (i.e. after interaction) ray at intersection point}
		if (IncludeErrors)
		{
			if (optics->DistributionType == 'F' || optics->DistributionType == 'f')
				CopyVec3(CosIn, LastDFXYZ);  // Apply diffuse errors relative to surface normal
			else
				CopyVec3(CosIn, CosRayOutElement); // Apply all other errors relative to the specularly-reflected direction

			Errors(myrng, CosIn, 2, &System->Sun,
				Stage->ElementList[k], optics, CosOut, LastDFXYZ);  //optical errors
			myrng_counter++;
			CopyVec3(CosRayOutElement, CosOut);
		}
	}
}

// PT Optimization Methods and Structs

struct eprojdat
{
	TElement* el_addr;
	double d_proj;
	double az;
	double zen;

	eprojdat() {};
	eprojdat(TElement* e, double d, double a, double z)
	{
		el_addr = e;
		d_proj = d;
		az = a;
		zen = z;
	};
};

//Comparison function for sorting vector of eprojdat
static bool eprojdat_compare_refactored(const eprojdat& A, const eprojdat& B)
{
	return A.d_proj > B.d_proj;
};

void SetupPTOptimizations(
	// system info
	TSystem* System, const bool AsPowerTower,

	// outputs
	st_hash_tree &sun_hash, st_hash_tree &rec_hash, double(&reccm_helio)[3]
	)
{
	//Calculate the center of mass of the receiver stage (StageList[1]) in heliostat stage coordinates.
	double reccm[] = { 0., 0., 0. };
	int nelrec = 0;
	if (AsPowerTower)
	{
		for (st_uint_t j = 0; j < System->StageList[1]->ElementList.size(); j++)
		{
			TElement* el = System->StageList[1]->ElementList.at(j);

			if (!el->Enabled)
				continue;

			nelrec++;

			for (int jj = 0; jj < 3; jj++)
				reccm[jj] += el->Origin[jj];
		}
		for (int jj = 0; jj < 3; jj++)
			reccm[jj] /= (double)nelrec;    //average


		//Transform to reference 
		double dum1[] = { 0., 0., 1. };
		double dum2[3];
		double reccm_global[3];
		TransformToReference(reccm, dum1, System->StageList[1]->Origin, System->StageList[1]->RLocToRef, reccm_global, dum2);

		//Transform to local (heliostat). reccm_helio is the x,y,z position of the receiver centroid in heliostat stage coordinates.
		TransformToLocal(reccm_global, dum1, System->StageList[0]->Origin, System->StageList[0]->RRefToLoc, reccm_helio, dum2);
	}
	//Create an array that stores the element address and the projected size in polar coordinates
	vector<eprojdat> el_proj_dat;
	el_proj_dat.reserve(System->StageList[0]->ElementList.size());

	//calculate the smallest zone size. This should be on the order of the largest element in the stage. 
	//load stage 0 elements into the mesh
	double d_elm_max = -9.e9;

	for (st_uint_t i = 0; i < System->StageList[0]->ElementList.size(); i++)
	{
		TElement* el = System->StageList[0]->ElementList.at(i);

		el->element_number = i + 1;   //use index for element number

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
				d_elm = el->ParameterA;
				break;
				//rectangular aperture
			case 'r':
			case 'R':
				d_elm = sqrt(el->ParameterA * el->ParameterA + el->ParameterB * el->ParameterB);
				break;
				//annular aperture
			case 'a':
			case 'A':
				d_elm = el->ParameterB;
				break;
			case 'l':
			case 'L':
				//off axis aperture section of line focus trough  or cylinder
				d_elm = sqrt(el->ParameterB * el->ParameterB * 4. + el->ParameterC * el->ParameterC);
				break;
				//Irregular triangle
			case 'i':
			case 'I':
				//irregular quadrilateral
			case 'q':
			case 'Q':
			{
				double xmax = fmax(el->ParameterA, fmax(el->ParameterC, el->ParameterE));
				double xmin = fmin(el->ParameterA, fmin(el->ParameterC, el->ParameterE));
				double ymax = fmax(el->ParameterB, fmax(el->ParameterD, el->ParameterF));
				double ymin = fmin(el->ParameterB, fmin(el->ParameterD, el->ParameterF));

				if (el->ShapeIndex == 'q' || el->ShapeIndex == 'Q')
				{
					xmax = fmax(xmax, el->ParameterG);
					xmin = fmin(xmin, el->ParameterG);
					ymax = fmax(ymax, el->ParameterH);
					ymin = fmin(ymin, el->ParameterH);
				}

				double dx = xmax - xmin;
				double dy = ymax - ymin;

				d_elm = sqrt(dx * dx + dy * dy);

				break;
			}
			default:
				break;
		}

		d_elm_max = fmax(d_elm_max, d_elm);

		if (AsPowerTower)
		{
			//Calculate the distance from the receiver to the element and the max projected size
			double dX[3];
			for (int jj = 0; jj < 3; jj++)
				dX[jj] = el->Origin[jj] - reccm_helio[jj];  //vector from receiver to heliostat (not unitized)
			double r_elm = 0.;
			for (int jj = 0; jj < 3; jj++)
				r_elm += dX[jj] * dX[jj];
			r_elm = sqrt(r_elm);            //vector length
			double d_elm_proj = d_elm / r_elm;  //Projected size of the element from the view of the receiver (radians)

			//calculate az,zen coordinate
			double az, zen;
			az = atan2(dX[0] / r_elm, dX[1] / r_elm);       //Az coordinate of the heliostat from the receiver's perspective
			zen = asin(dX[2] / r_elm);                    //Zen coordinate """"

			el_proj_dat.push_back(eprojdat(el, d_elm_proj, az, zen));
		}
	}

	if (AsPowerTower)
	{
		//Sort the polar projections by size, largest to smallest
		std::sort(el_proj_dat.begin(), el_proj_dat.end(), eprojdat_compare_refactored);
	}

	//set up the layout data object that provides configuration details for the hash tree
	KDLayoutData sun_ld;
	sun_ld.xlim[0] = System->Sun.MinXSun;
	sun_ld.xlim[1] = System->Sun.MaxXSun;
	sun_ld.ylim[0] = System->Sun.MinYSun;
	sun_ld.ylim[1] = System->Sun.MaxYSun;
	sun_ld.min_unit_dx = d_elm_max;
	sun_ld.min_unit_dy = d_elm_max;

	sun_hash.create_mesh(sun_ld);

	//load stage 0 elements into the mesh
	for (st_uint_t i = 0; i < System->StageList[0]->ElementList.size(); i++)
	{
		TElement* el = System->StageList[0]->ElementList.at(i);
		sun_hash.add_object((void*)el, el->PosSunCoords[0], el->PosSunCoords[1]);
	}

	//calculate and associate neighbors with each zone
	sun_hash.add_neighborhood_data();

	if (AsPowerTower)
	{
		//Set things up for the polar coordinate tree
		KDLayoutData rec_ld;
		rec_ld.xlim[0] = -M_PI;
		rec_ld.xlim[1] = M_PI;
		rec_ld.ylim[0] = -M_PI / 2.;
		rec_ld.ylim[1] = M_PI / 2.;
		//use smallest element to set the minimum size
		rec_ld.min_unit_dx = rec_ld.min_unit_dy = el_proj_dat.back().d_proj; //radians at equator

		rec_hash.create_mesh(rec_ld);

		//load stage 0 elements into the receiver mesh in the order of largest projection to smallest
		for (int i = 0; i < el_proj_dat.size(); i++)
		{
			eprojdat* D = &el_proj_dat.at(i);

			//Calculate the angular span of the element
			double angspan[2];
			double adjmult = 1.5;
			angspan[0] = D->d_proj / cos(fabs(D->zen)) * adjmult;   //azimuthal span
			angspan[0] = fmin(angspan[0], 2. * M_PI);     //limit to circumference 
			angspan[1] = D->d_proj / M_PI * adjmult;    //zenithal span
			rec_hash.add_object((void*)D->el_addr, D->az, D->zen, angspan);
		}

		//associate neighbors with each zone
		rec_hash.add_neighborhood_data();
	}
}

st_uint_t GetPTElements(
	// system info
	const bool AsPowerTower,
	
	// Stage info
	const TStage* Stage, const int i,

	// Ray info
	const bool in_multi_hit_loop, const double(&PosRayStage)[3],
	const double(&reccm_helio)[3], st_hash_tree& rec_hash,

	const vector<void*> &sunint_elements,

	// Outputs
	vector<void*> &reflint_elements,
	bool& has_elements
	)
{
	st_uint_t nintelements = 0;

	if (i == 0)
	{
		if (in_multi_hit_loop)
		{
			if (AsPowerTower)
			{
				//>=Second time through - checking for first stage multiple element interactions

				//get ray position in receiver polar coordinates
				double raypvec[3];
				for (int jj = 0; jj < 3; jj++)
					raypvec[jj] = PosRayStage[jj] - reccm_helio[jj];
				double raypvecmag = sqrt(raypvec[0] * raypvec[0] + raypvec[1] * raypvec[1] + raypvec[2] * raypvec[2]);
				double raypol[2];
				raypol[0] = atan2(raypvec[0], raypvec[1]);
				raypol[1] = asin(raypvec[2] / raypvecmag);
				//get elements in the vicinity of the ray's polar coordinates
				reflint_elements.clear();
				rec_hash.get_all_data_at_loc(reflint_elements, raypol[0], raypol[1]);
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
			if (has_elements)
				nintelements = sunint_elements.size();
			else
				nintelements = 0;
		}
	}
	else
		nintelements = Stage->ElementList.size();

	return nintelements;
}

// Trace method

bool Trace_refactored(TSystem* System, unsigned int seed,
	st_uint_t NumberOfRays,
	st_uint_t MaxNumberOfRays,
	bool IncludeSunShape,
	bool IncludeErrors,
	bool AsPowerTower,
	int (*callback)(st_uint_t ntracedtotal, st_uint_t ntraced, st_uint_t ntotrace, st_uint_t curstage, st_uint_t nstages, void* data),
	void* cbdata)
{
	// Determine if PT optimizations should be applied
	bool PT_override = false;
	if (System->StageList.size() > 0
		&& (System->StageList[0]->ElementList.size() < 10    //the first stage contains only a few elements
			|| System->StageList.size() == 1)                //there's only one stage
		)
	{
		PT_override = true;
	}

	// Initialize variables
	MTRand myrng(seed);
	int myrng_counter = 0;

	// Initialize Internal State Variables
	st_uint_t RayNumber = 1;						// Ray Number of current ray
	bool PreviousStageHasRays = false;
	st_uint_t LastRayNumberInPreviousStage = NumberOfRays;

	// Check Inputs
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

	// Define IncomingRays
	std::vector<GlobalRay_refactored> IncomingRays;	// Vector of rays from previous stage, going into next stage
	try
	{
		IncomingRays.resize(NumberOfRays);
	}
	catch (std::exception& e) {
		System->errlog("Incoming rays resize exception: %d, '%s'", NumberOfRays, e.what());
		return false;
	}

	// Initialize Sun
	double PosSunStage[3] = { 0.0, 0.0, 0.0 };
	if (!SunToPrimaryStage(System, System->StageList[0], &System->Sun, PosSunStage))
		return false;

	// Calculate hash tree for reflection to receiver plane(polar coordinates).
	st_hash_tree sun_hash;
	st_hash_tree rec_hash;
	double reccm_helio[3];  //receiver centroid in heliostat field coordinates
	if (!PT_override)
	{
		SetupPTOptimizations(System, AsPowerTower, sun_hash, rec_hash, reccm_helio);
	}

	// Start the clock
	clock_t startTime = clock();
	int rays_per_callback_estimate = 50;
	st_uint_t RaysTracedTotal = 0;

	// Loop through stages
	for (st_uint_t i = 0; i < System->StageList.size(); i++)
	{
		// Check if previous stage has rays
		bool StageHasRays = true;
		if (i > 0 && PreviousStageHasRays == false)
		{
			StageHasRays = false;
		}

		// Get Current Stage
		TStage* Stage = System->StageList[i];

		// Initialize stage variables
		st_uint_t StageDataArrayIndex = 0;
		st_uint_t PreviousStageDataArrayIndex = 0;

		// Loop through rays
		while (StageHasRays)
		{
			// Initialize Global Coordinates
			double PosRayGlob[3] = { 0.0, 0.0, 0.0 };
			double CosRayGlob[3] = { 0.0, 0.0, 0.0 };

			// Initialize Stage Coordinates
			double PosRayStage[3] = { 0.0, 0.0, 0.0 };
			double CosRayStage[3] = { 0.0, 0.0, 0.0 };

			// Initialize PT Optimization variables
			bool has_elements = true;
			vector<void*> sunint_elements;

			// Get Ray
			if (i == 0)
			{
				// Make ray (if first stage)
				double PosRaySun[3];
				GenerateRay(myrng, PosSunStage, Stage->Origin,
					Stage->RLocToRef, &System->Sun,
					PosRayGlob, CosRayGlob, PosRaySun);
				myrng_counter++;
				System->SunRayCount++;

				// If using PT optimizations, check if stage has elements that could interact with ray
				if (!PT_override)
				{
					has_elements = sun_hash.get_all_data_at_loc(sunint_elements, PosRaySun[0], PosRaySun[1]);
				}
			}
			else
			{
				// Get ray from previous stage
				RayNumber = IncomingRays[StageDataArrayIndex].Num;
				CopyVec3(PosRayGlob, IncomingRays[StageDataArrayIndex].Pos);
				CopyVec3(CosRayGlob, IncomingRays[StageDataArrayIndex].Cos);
				StageDataArrayIndex++;
			}

			// transform the global incoming ray to local stage coordinates
			TransformToLocal(PosRayGlob, CosRayGlob,
				Stage->Origin, Stage->RRefToLoc,
				PosRayStage, CosRayStage);

			// Update callback
			if (callback != 0
				&& RaysTracedTotal++ % rays_per_callback_estimate == 0)
			{
				if (RaysTracedTotal > 1)
				{
					//update how often to call this
					double msec_per_ray = 1000. * (clock() - startTime) / CLOCKS_PER_SEC / (double)(RaysTracedTotal > 0 ? RaysTracedTotal : 1);
					//set the new callback estimate to be about 50 ms
					rays_per_callback_estimate = (int)(200. / msec_per_ray);
					//limit to something reasonable
					rays_per_callback_estimate = rays_per_callback_estimate < 5 ? 5 : rays_per_callback_estimate;
				}

				//do the callback
				if (!(*callback)(RaysTracedTotal, RayNumber,
					LastRayNumberInPreviousStage, i + 1,
					System->StageList.size(), cbdata))
					return true;
			}


			// Initialize internal variables for ray intersection tracing
			bool RayInStage = true;
			bool in_multi_hit_loop = false;
			double LastPosRaySurfElement[3] = { 0.0, 0.0, 0.0 };
			double LastCosRaySurfElement[3] = { 0.0, 0.0, 0.0 };
			double LastPosRaySurfStage[3] = { 0.0, 0.0, 0.0 };
			double LastCosRaySurfStage[3] = { 0.0, 0.0, 0.0 };
			double LastDFXYZ[3] = { 0.0, 0.0, 0.0 };
			st_uint_t LastElementNumber = 0;
			st_uint_t LastRayNumber = 0;
			int ErrorFlag;
			int LastHitBackSide;
			bool StageHit;
			int MultipleHitCount = 0;
			double PosRayOutElement[3] = { 0.0, 0.0, 0.0 };
			double CosRayOutElement[3] = { 0.0, 0.0, 0.0 };

			// Start Loop to trace ray until it leaves stage
			bool RayIsAbsorbed = false;
			while (RayInStage)
			{
				// Set number of elements to search through
				st_uint_t nintelements = 0;
				vector<void*> reflint_elements;
				if (!PT_override)	// if using opt AND first stage
				{
					nintelements = GetPTElements(AsPowerTower, Stage, i, in_multi_hit_loop, PosRayStage,
						reccm_helio, rec_hash, sunint_elements,
						reflint_elements, has_elements);
				}
				else
				{
					nintelements = Stage->ElementList.size();
				}

				// Find the element the ray hits
				FindElementHit(i, Stage, PT_override, AsPowerTower,
					nintelements, sunint_elements, reflint_elements,
					RayNumber, in_multi_hit_loop,
					PosRayStage, CosRayStage,

					LastPosRaySurfElement, LastCosRaySurfElement, LastDFXYZ,
					LastElementNumber, LastRayNumber,
					LastPosRaySurfStage, LastCosRaySurfStage,
					ErrorFlag, LastHitBackSide, StageHit);

				// Breakout if ray left stage
				if (!StageHit)
				{
					RayInStage = false;
					break;
				}

				// Add ray to Stage RayData
				TRayData::ray_t* p_ray = Stage->RayData.Append(LastPosRaySurfStage,
					LastCosRaySurfStage,
					LastElementNumber,
					i + 1,
					LastRayNumber);

				// Check p_ray saved correctly
				if (!p_ray)
				{
					System->errlog("Failed to save ray data at index %d", Stage->RayData.Count() - 1);
					return false;
				}

				// Skipping LastElementNumber == 0 check

				// Increment MultipleHitCount
				MultipleHitCount++;

				// Get optics and check for absorption
				TOpticalProperties* optics = 0;
				if (Stage->Virtual)
				{
					// If stage is virtual, there is no interaction
					CopyVec3(PosRayOutElement, LastPosRaySurfElement);
					CopyVec3(CosRayOutElement, LastCosRaySurfElement);
				}
				else
				{
					// trace through the interaction
					TElement* optelm = Stage->ElementList[p_ray->element - 1];

					if (LastHitBackSide)
						optics = &optelm->Optics->Back;
					else
						optics = &optelm->Optics->Front;

					double TestValue;
					double UnitLastDFXYZ[3] = { 0.0, 0.0, 0.0 };
					double IncidentAngle = 0;
					switch (optelm->InteractionType)
					{
						case 1: // refraction
							if (optics->UseTransmissivityTable)
							{
								int npoints = optics->TransmissivityTable.size();
								int m = 0;

								UnitLastDFXYZ[0] = -LastDFXYZ[0] / sqrt(DOT(LastDFXYZ, LastDFXYZ));
								UnitLastDFXYZ[1] = -LastDFXYZ[1] / sqrt(DOT(LastDFXYZ, LastDFXYZ));
								UnitLastDFXYZ[2] = -LastDFXYZ[2] / sqrt(DOT(LastDFXYZ, LastDFXYZ));
								IncidentAngle = acos(DOT(LastCosRaySurfElement, UnitLastDFXYZ)) * 1000.;  //[mrad]
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

							if (optics->UseReflectivityTable)
							{
								int npoints = optics->ReflectivityTable.size();
								int m = 0;
								UnitLastDFXYZ[0] = -LastDFXYZ[0] / sqrt(DOT(LastDFXYZ, LastDFXYZ));
								UnitLastDFXYZ[1] = -LastDFXYZ[1] / sqrt(DOT(LastDFXYZ, LastDFXYZ));
								UnitLastDFXYZ[2] = -LastDFXYZ[2] / sqrt(DOT(LastDFXYZ, LastDFXYZ));
								IncidentAngle = acos(DOT(LastCosRaySurfElement, UnitLastDFXYZ)) * 1000.;  //[mrad]
								if (IncidentAngle >= optics->ReflectivityTable[npoints - 1].angle)
								{
									TestValue = optics->ReflectivityTable[npoints - 1].refl;
								}
								else
								{
									while (optics->ReflectivityTable[m].angle < IncidentAngle)
										m++;

									if (m == 0)
										TestValue = optics->ReflectivityTable[m].refl;
									else
										TestValue = (optics->ReflectivityTable[m].refl + optics->ReflectivityTable[m - 1].refl) / 2.0;
								}
							}
							else
								TestValue = optics->Reflectivity;
							break;
						default:
							System->errlog("Bad optical interaction type = %d (stage %d)", i, optelm->InteractionType);
							return false;
					}

					//  {Apply MonteCarlo probability of absorption. Limited for now, but can make more complex later on if desired}
					if (TestValue <= myrng())
					{
						myrng_counter++;
						// ray was fully absorbed, so indicate by negating the element number
						p_ray->element = 0 - p_ray->element;
						RayIsAbsorbed = true;
						break;
					}

				}

				// Process Interaction
				int k = abs(p_ray->element) - 1;
				ProcessInteraction(System, myrng, IncludeSunShape, optics, IncludeErrors,
					i, Stage, k,
					MultipleHitCount, LastDFXYZ,
					LastCosRaySurfElement, ErrorFlag,
					CosRayOutElement, LastPosRaySurfElement,
					PosRayOutElement, myrng_counter);

				// Transform ray back to stage coordinate system
				TransformToReference(PosRayOutElement, CosRayOutElement,
					Stage->ElementList[k]->Origin, Stage->ElementList[k]->RLocToRef,
					PosRayStage, CosRayStage);
				TransformToReference(PosRayStage, CosRayStage,
					Stage->Origin, Stage->RLocToRef,
					PosRayGlob, CosRayGlob);

				// Break out if multiple hits are not allowed
				if (!Stage->MultiHitsPerRay)
				{
					StageHit = false;
					break;
				}
				else
				{
					in_multi_hit_loop = true;
				}

			}

			// Handle if Ray was absorbed
			if (RayIsAbsorbed)
			{
				// ray was fully absorbed
				if (RayNumber == LastRayNumberInPreviousStage)
				{
					PreviousStageHasRays = false;
					if (PreviousStageDataArrayIndex > 0)
					{
						PreviousStageDataArrayIndex--;
						PreviousStageHasRays = true;
					}
					//goto Label_EndStageLoop;
					break;
				}
				else
				{
					if (i == 0)
					{
						if (RayNumber == NumberOfRays)
							//goto Label_EndStageLoop;
							break;
						else
							RayNumber++;
					}

					// Next ray in loop
					continue;
				}
			}

			// !StageHit logic goes here
			if (StageHit == true)
			{
				// This shouldn't happen...
			}

			// Ray has left the stage
			bool FlagMiss = false;
			if (i == 0)
			{
				if (MultipleHitCount == 0)
				{
					// Ray in first stage missed stage entirely
					// Generate new ray
					continue;
				}
				else
				{
					// Ray hit an element, so save it for next stage
					CopyVec3(IncomingRays[PreviousStageDataArrayIndex].Pos, PosRayGlob);
					CopyVec3(IncomingRays[PreviousStageDataArrayIndex].Cos, CosRayGlob);
					IncomingRays[PreviousStageDataArrayIndex].Num = RayNumber;

					// Is Ray the last in the stage?
					if (RayNumber == NumberOfRays)
					{
						StageHasRays = false;
						break;
					}

					PreviousStageDataArrayIndex++;
					PreviousStageHasRays = true;

					// Move on to next ray
					RayNumber++;
					continue;
				}
			}
			else
			{
				// After the first stage
				// Ray hit element OR is traced through stage
				if (Stage->TraceThrough || MultipleHitCount > 0)
				{
					// Ray is saved for the next stage
					CopyVec3(IncomingRays[PreviousStageDataArrayIndex].Pos, PosRayGlob);
					CopyVec3(IncomingRays[PreviousStageDataArrayIndex].Cos, CosRayGlob);
					IncomingRays[PreviousStageDataArrayIndex].Num = RayNumber;

					// Check if ray is last in stage
					if (RayNumber == LastRayNumberInPreviousStage)
					{
						StageHasRays = false;
						break;
					}

					PreviousStageDataArrayIndex++;
					PreviousStageHasRays = true;

					if (MultipleHitCount == 0)
					{
						FlagMiss = true;
					}

					// Go to next ray
					continue;

				}
				// Ray missed stage entirely and is not traced
				else
				{
					FlagMiss = true;
				}

				// Handle FlagMiss condition (
				if (FlagMiss == true)
				{
					LastElementNumber = 0;
					LastRayNumber = RayNumber;
					CopyVec3(LastPosRaySurfStage, PosRayStage);
					CopyVec3(LastCosRaySurfStage, CosRayStage);

					// Copying this here to handle FlagMiss condition
					TRayData::ray_t* p_ray = Stage->RayData.Append(LastPosRaySurfStage,
						LastCosRaySurfStage,
						LastElementNumber,
						i + 1,
						LastRayNumber);

					if (RayNumber == LastRayNumberInPreviousStage)
					{
						if (!Stage->TraceThrough)
						{
							PreviousStageHasRays = false;
							if (PreviousStageDataArrayIndex > 0)
							{
								PreviousStageHasRays = true;
								PreviousStageDataArrayIndex--; // last ray was previous one
							}
						}
						
						// Exit stage
						StageHasRays = false;
						break;
					}
					else
					{
						if (i == 0) RayNumber++; // generate new sun ray
						
						// Start new ray
						continue;
					}
				}
			}
		}

		// EndStage section...

		// skipping save_st_data logic

		if (!PreviousStageHasRays)
		{
			LastRayNumberInPreviousStage = 0;
			continue;	// No rays to carry forward
		}

		if (PreviousStageDataArrayIndex < IncomingRays.size())
		{
			LastRayNumberInPreviousStage = IncomingRays[PreviousStageDataArrayIndex].Num;
			if (LastRayNumberInPreviousStage == 0)
			{
				size_t pp = IncomingRays[PreviousStageDataArrayIndex - 1].Num;
				System->errlog("LastRayNumberInPreviousStage=0, stage %d, PrevIdx=%d, CurIdx=%d, pp=%d", i + 1,
					PreviousStageDataArrayIndex, StageDataArrayIndex, pp);
				return false;
			}
		}
		else
		{
			System->errlog("Invalid PreviousStageDataArrayIndex: %u, @ stage %d",
				PreviousStageDataArrayIndex, i + 1);
			return false;
		}

	}

	return true;
}





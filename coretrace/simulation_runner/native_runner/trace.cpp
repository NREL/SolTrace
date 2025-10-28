
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

#include "trace.hpp"

// Standard library headers
#include <cmath>
#include <limits>
#include <vector>

// SimulationData headers
#include "constants.hpp"
#include "matvec.hpp"
#include "simulation_data_export.hpp"

// NativeRunner headers
#include "find_element_hit.hpp"
#include "generate_ray.hpp"
#include "native_runner_types.hpp"
#include "process_interaction.hpp"
#include "pt_optimizations.hpp"
#include "sun_to_primary_stage.hpp"
#include "treemesh.hpp"

namespace SolTrace::NativeRunner
{

	// #define   Order 3
	// #define   NumIterations 20
	// #define   Epsilon 0.000001

	using SolTrace::Result::RayEvent;

	// Trace method
	bool trace_native(
		TSystem *System,
		unsigned int seed,
		uint_fast64_t NumberOfRays,
		uint_fast64_t MaxNumberOfRays,
		bool IncludeSunShape,
		bool IncludeErrors,
		bool AsPowerTower)
	{
		// Determine if PT optimizations should be applied
		bool PT_override = false;
		if (System->StageList.size() > 0 &&
			(System->StageList[0]->ElementList.size() < 10 || System->StageList.size() == 1))
		{
			PT_override = true;
		}

		// Initialize variables
		// std::cout << "Seed: " << seed << std::endl;
		MTRand myrng(seed);
		int myrng_counter = 0;

		// Initialize Internal State Variables
		uint_fast64_t RayNumber = 1; // Ray Number of current ray
		bool PreviousStageHasRays = false;
		uint_fast64_t LastRayNumberInPreviousStage = NumberOfRays;

		// Define IncomingRays
		std::vector<GlobalRay_refactored> IncomingRays; // Vector of rays from previous stage, going into next stage
		IncomingRays.resize(NumberOfRays);

		// Initialize Sun
		double PosSunStage[3] = {0.0, 0.0, 0.0};
		if (!SunToPrimaryStage(System,
							   System->StageList[0].get(),
							   &System->Sun,
							   PosSunStage))
			return false;

		// Calculate hash tree for reflection to receiver plane(polar coordinates).
		st_hash_tree sun_hash;
		st_hash_tree rec_hash;
		double reccm_helio[3]; // receiver centroid in heliostat field coordinates
		if (!PT_override)
		{
			SetupPTOptimizations(System, AsPowerTower, sun_hash,
								 rec_hash, reccm_helio);
		}

		// Start the clock
		clock_t startTime = clock();
		// int rays_per_callback_estimate = 50;
		uint_fast64_t RaysTracedTotal = 0;

		// Initialize stage variables
		uint_fast64_t StageDataArrayIndex = 0;
		uint_fast64_t PreviousStageDataArrayIndex = 0;
		uint_fast64_t n_rays_active = NumberOfRays;

		// Loop through stages
		for (uint_fast64_t i = 0; i < System->StageList.size(); i++)
		{
			// std::cout << "Processing stage " << i << "..." << std::endl;
			// Check if previous stage has rays
			bool StageHasRays = true;
			if (i > 0 && PreviousStageHasRays == false)
			{
				StageHasRays = false;
			}

			// Get Current Stage
			tstage_ptr Stage = System->StageList[i];

			// Initialize stage variables
			StageDataArrayIndex = 0;
			PreviousStageDataArrayIndex = 0;

			// Loop through rays
			while (StageHasRays)
			{
				// Initialize Global Coordinates
				double PosRayGlob[3] = {0.0, 0.0, 0.0};
				double CosRayGlob[3] = {0.0, 0.0, 0.0};

				// Initialize Stage Coordinates
				double PosRayStage[3] = {0.0, 0.0, 0.0};
				double CosRayStage[3] = {0.0, 0.0, 0.0};

				// Initialize PT Optimization variables
				bool has_elements = true;
				std::vector<void *> sunint_elements;

				// Get Ray
				if (i == 0)
				{
					// TODO: This function seems to ignore the MaxNumberOfRays
					// argument. Should fix that.

					// Make ray (if first stage)
					double PosRaySun[3];
					GenerateRay(myrng, PosSunStage, Stage->Origin,
								Stage->RLocToRef, &System->Sun,
								PosRayGlob, CosRayGlob, PosRaySun);
					myrng_counter++;
					System->SunRayCount++;

					// If using PT optimizations, check if stage has elements
					// that could interact with ray
					if (!PT_override)
					{
						has_elements =
							sun_hash.get_all_data_at_loc(sunint_elements,
														 PosRaySun[0],
														 PosRaySun[1]);
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

				// Initialize internal variables for ray intersection tracing
				bool RayInStage = true;
				bool in_multi_hit_loop = false;
				double LastPosRaySurfElement[3] = {0.0, 0.0, 0.0};
				double LastCosRaySurfElement[3] = {0.0, 0.0, 0.0};
				double LastPosRaySurfStage[3] = {0.0, 0.0, 0.0};
				double LastCosRaySurfStage[3] = {0.0, 0.0, 0.0};
				double LastDFXYZ[3] = {0.0, 0.0, 0.0};
				uint_fast64_t LastElementNumber = 0;
				uint_fast64_t LastRayNumber = 0;
				int ErrorFlag;
				int LastHitBackSide;
				bool StageHit;
				int MultipleHitCount = 0;
				double PosRayOutElement[3] = {0.0, 0.0, 0.0};
				double CosRayOutElement[3] = {0.0, 0.0, 0.0};

				// Start Loop to trace ray until it leaves stage
				bool RayIsAbsorbed = false;
				while (RayInStage)
				{
					// Set number of elements to search through
					uint_fast64_t nintelements = 0;
					std::vector<void *> reflint_elements;
					if (!PT_override) // if using opt AND first stage
					{
						nintelements = GetPTElements(AsPowerTower, Stage, i,
													 in_multi_hit_loop, PosRayStage,
													 reccm_helio, rec_hash,
													 sunint_elements,
													 reflint_elements, has_elements);
					}
					else
					{
						nintelements = Stage->ElementList.size();
					}

					// Find the element the ray hits
					FindElementHit(i, Stage, PT_override, AsPowerTower,
								   nintelements, sunint_elements,
								   reflint_elements,
								   RayNumber, in_multi_hit_loop,
								   PosRayStage, CosRayStage,
								   LastPosRaySurfElement,
								   LastCosRaySurfElement,
								   LastDFXYZ,
								   LastElementNumber, LastRayNumber,
								   LastPosRaySurfStage, LastCosRaySurfStage,
								   ErrorFlag, LastHitBackSide, StageHit);

					// Breakout if ray left stage
					if (!StageHit)
					{
						RayInStage = false;
						break;
					}

					// Increment MultipleHitCount
					MultipleHitCount++;

					if (i == 0 && MultipleHitCount == 1)
					{
						System->RayData.Append(PosRayGlob,
											   CosRayGlob,
											   ELEMENT_NULL,
											   i + 1,
											   LastRayNumber,
											   RayEvent::CREATE);
					}

					// Get optics and check for absorption
					const OpticalProperties *optics = 0;
					RayEvent rev = RayEvent::VIRTUAL;
					if (Stage->Virtual)
					{
						// If stage is virtual, there is no interaction
						CopyVec3(PosRayOutElement, LastPosRaySurfElement);
						CopyVec3(CosRayOutElement, LastCosRaySurfElement);
					}
					else
					{
						// trace through the interaction
						telement_ptr optelm = Stage->ElementList[LastElementNumber - 1];

						if (LastHitBackSide)
							optics = &optelm->Optics.Back;
						else
							optics = &optelm->Optics.Front;

						double TestValue;
						double UnitLastDFXYZ[3] = {0.0, 0.0, 0.0};
						double IncidentAngle = 0;
						// switch (optelm->InteractionType)
						switch (optics->my_type)
						{
						case InteractionType::REFRACTION: // refraction
							// TODO: Implement transmissivity table?
							// if (optics->UseTransmissivityTable)
							// {
							// 	int npoints = optics->TransmissivityTable.size();
							// 	int m = 0;

							// 	UnitLastDFXYZ[0] = -LastDFXYZ[0] / sqrt(DOT(LastDFXYZ, LastDFXYZ));
							// 	UnitLastDFXYZ[1] = -LastDFXYZ[1] / sqrt(DOT(LastDFXYZ, LastDFXYZ));
							// 	UnitLastDFXYZ[2] = -LastDFXYZ[2] / sqrt(DOT(LastDFXYZ, LastDFXYZ));
							// 	IncidentAngle = acos(DOT(LastCosRaySurfElement, UnitLastDFXYZ)) * 1000.; //[mrad]
							// 	if (IncidentAngle >= optics->TransmissivityTable[npoints - 1].angle)
							// 	{
							// 		TestValue = optics->TransmissivityTable[npoints - 1].trans;
							// 	}
							// 	else
							// 	{
							// 		while (optics->TransmissivityTable[m].angle < IncidentAngle)
							// 			m++;

							// 		if (m == 0)
							// 			TestValue = optics->TransmissivityTable[m].trans;
							// 		else
							// 			TestValue = (optics->TransmissivityTable[m].trans + optics->TransmissivityTable[m - 1].trans) / 2.0;
							// 	}
							// }
							// else
							// 	TestValue = optics->Transmissivity;
							TestValue = optics->transmitivity;
							rev = RayEvent::TRANSMIT;
							break;
						case InteractionType::REFLECTION: // reflection
							// TODO: Implement reflectivity table?
							// if (optics->UseReflectivityTable)
							// {
							// 	int npoints = optics->ReflectivityTable.size();
							// 	int m = 0;
							// 	UnitLastDFXYZ[0] = -LastDFXYZ[0] / sqrt(DOT(LastDFXYZ, LastDFXYZ));
							// 	UnitLastDFXYZ[1] = -LastDFXYZ[1] / sqrt(DOT(LastDFXYZ, LastDFXYZ));
							// 	UnitLastDFXYZ[2] = -LastDFXYZ[2] / sqrt(DOT(LastDFXYZ, LastDFXYZ));
							// 	IncidentAngle = acos(DOT(LastCosRaySurfElement, UnitLastDFXYZ)) * 1000.; //[mrad]
							// 	if (IncidentAngle >= optics->ReflectivityTable[npoints - 1].angle)
							// 	{
							// 		TestValue = optics->ReflectivityTable[npoints - 1].refl;
							// 	}
							// 	else
							// 	{
							// 		while (optics->ReflectivityTable[m].angle < IncidentAngle)
							// 			m++;

							// 		if (m == 0)
							// 			TestValue = optics->ReflectivityTable[m].refl;
							// 		else
							// 			TestValue = (optics->ReflectivityTable[m].refl + optics->ReflectivityTable[m - 1].refl) / 2.0;
							// 	}
							// }
							// else
							// 	TestValue = optics->Reflectivity;
							TestValue = optics->reflectivity;
							rev = RayEvent::REFLECT;
							break;
						default:
							System->errlog(
								"Bad optical interaction type = %d (stage %d)",
								i, optics->my_type);
							return false;
						}

						// Apply MonteCarlo probability of absorption. Limited
						// for now, but can make more complex later on if desired
						// if (TestValue <= myrng())
						double flip = myrng();
						// if (RayNumber == 36)
						// {
						// 	std::cout << "RAY NUMBER 36!" << std::endl;
						// }
						// if (i > 0)
						// {
						// 	std::cout << "Ray Number: " << RayNumber
						// 			  << "\nStage: " << i + 1
						// 			  << "\nMultiHit: " << Stage->MultiHitsPerRay
						// 			  << "\nmyrng_counter: " << myrng_counter
						// 			  << "\nTestValue: " << TestValue
						// 			  << "\nflip: " << flip
						// 			  << "\nRayIsAbsorbed: " << (TestValue < flip)
						// 			  << "\nPosRayGlob: " << Vector3d(PosRayGlob)
						// 			  << "\nCosRayGlob: " << Vector3d(CosRayGlob)
						// 			  << "\nDFXYZ: " << Vector3d(LastDFXYZ)
						// 			  << std::endl;
						// }
						if (TestValue <= flip)
						{
							myrng_counter++;
							// ray was fully absorbed
							RayIsAbsorbed = true;
							break;
						}
					}

					// Process Interaction
					int_fast64_t k = LastElementNumber - 1;
					ProcessInteraction(System, myrng, IncludeSunShape,
									   optics,
									   IncludeErrors,
									   i, Stage, // k,
									   MultipleHitCount, LastDFXYZ,
									   LastCosRaySurfElement, ErrorFlag,
									   CosRayOutElement, LastPosRaySurfElement,
									   PosRayOutElement, myrng_counter);

					// Transform ray back to stage coordinate system
					TransformToReference(PosRayOutElement, CosRayOutElement,
										 Stage->ElementList[k]->Origin,
										 Stage->ElementList[k]->RLocToRef,
										 PosRayStage, CosRayStage);
					TransformToReference(PosRayStage, CosRayStage,
										 Stage->Origin, Stage->RLocToRef,
										 PosRayGlob, CosRayGlob);

					System->RayData.Append(PosRayGlob,
										   CosRayGlob,
										   LastElementNumber,
										   i + 1,
										   LastRayNumber,
										   rev);

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
					TransformToReference(LastPosRaySurfStage,
										 LastCosRaySurfStage,
										 Stage->Origin,
										 Stage->RLocToRef,
										 PosRayGlob,
										 CosRayGlob);

					System->RayData.Append(PosRayGlob,
										   CosRayGlob,
										   LastElementNumber,
										   i + 1,
										   LastRayNumber,
										   RayEvent::ABSORB);

					n_rays_active--;

					// ray was fully absorbed
					if (RayNumber == LastRayNumberInPreviousStage)
					{
						PreviousStageHasRays = false;
						if (PreviousStageDataArrayIndex > 0)
						{
							PreviousStageDataArrayIndex--;
							PreviousStageHasRays = true;
						}
						break;
					}
					else
					{
						if (i == 0)
						{
							if (RayNumber == NumberOfRays)
								break;
							else
								RayNumber++;
						}

						// Next ray in loop
						continue;
					}
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
						CopyVec3(IncomingRays[PreviousStageDataArrayIndex].Pos,
								 PosRayGlob);
						CopyVec3(IncomingRays[PreviousStageDataArrayIndex].Cos,
								 CosRayGlob);
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
						CopyVec3(IncomingRays[PreviousStageDataArrayIndex].Pos,
								 PosRayGlob);
						CopyVec3(IncomingRays[PreviousStageDataArrayIndex].Cos,
								 CosRayGlob);
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
						LastRayNumber = RayNumber;

						System->RayData.Append(PosRayGlob,
											   CosRayGlob,
											   ELEMENT_NULL,
											   i + 1,
											   LastRayNumber,
											   RayEvent::EXIT);

						n_rays_active--;

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
							if (i == 0)
								RayNumber++; // generate new sun ray

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
				continue; // No rays to carry forward
			}

			if (PreviousStageDataArrayIndex < IncomingRays.size())
			{
				LastRayNumberInPreviousStage = IncomingRays[PreviousStageDataArrayIndex].Num;
				if (LastRayNumberInPreviousStage == 0)
				{
					return false;
				}
			}
			else
			{
				return false;
			}
		}

		// Close out any remaining rays as misses
		unsigned idx = System->StageList.size() - 1;
		tstage_ptr Stage = System->StageList[idx];
		for (uint_fast64_t k = 0; k < n_rays_active; ++k)
		{
			GlobalRay_refactored ray = IncomingRays[k];
			System->RayData.Append(ray.Pos,
								   ray.Cos,
								   ELEMENT_NULL,
								   idx + 1,
								   ray.Num,
								   RayEvent::EXIT);
		}

		return true;
	}

} // namespace SolTrace::NativeRunner

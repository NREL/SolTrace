
#include "simulation_data_export.hpp"

#include "determine_element_intersection_new.hpp"
#include "find_element_hit.hpp"
#include "native_runner_types.hpp"

namespace SolTrace::NativeRunner
{

	void FindElementHit(
		// stage info
		const int i,
		const tstage_ptr Stage,
		const bool PT_override,
		const bool AsPowerTower,
		// element info
		const int nintelements,
		const std::vector<void *> &sunint_elements,
		const std::vector<void *> &reflint_elements,
		// ray info
		const int RayNumber,
		const bool in_multi_hit_loop,
		double (&PosRayStage)[3],
		double (&CosRayStage)[3],
		// outputs
		double (&LastPosRaySurfElement)[3],
		double (&LastCosRaySurfElement)[3],
		double (&LastDFXYZ)[3],
		uint_fast64_t &LastElementNumber,
		uint_fast64_t &LastRayNumber,
		double (&LastPosRaySurfStage)[3],
		double (&LastCosRaySurfStage)[3],
		int &ErrorFlag,
		int &LastHitBackSide,
		bool &StageHit)
	{
		// Initialize Variables
		double LastPathLength = 1e99;
		int HitBackSide = 0;
		int InterceptFlag = 0;
		double DFXYZ[3] = {0.0, 0.0, 0.0};
		double PosRayElement[3] = {0.0, 0.0, 0.0};
		double CosRayElement[3] = {0.0, 0.0, 0.0};
		double PosRaySurfStage[3] = {0.0, 0.0, 0.0};
		double CosRaySurfStage[3] = {0.0, 0.0, 0.0};
		double PosRaySurfElement[3] = {0.0, 0.0, 0.0};
		double CosRaySurfElement[3] = {0.0, 0.0, 0.0};
		StageHit = false;

		for (uint_fast64_t j = 0; j < nintelements; j++)
		{
			TElement *Element; // = Stage->ElementList[j];
			if (i == 0 && !PT_override)
			{
				if (in_multi_hit_loop)
				{
					if (AsPowerTower)
					{
						Element = (TElement *)reflint_elements.at(j);
					}
					else
					{
						Element = Stage->ElementList[j].get();
					}
				}
				else
				{
					Element = (TElement *)sunint_elements.at(j);
				}
			}
			else
			{
				Element = Stage->ElementList[j].get();
			}

			// if (!Element->Enabled)
			// 	continue;

			//  {Transform ray to element[j] coord system of Stage[i]}
			TransformToLocal(PosRayStage, CosRayStage,
							 Element->Origin, Element->RRefToLoc,
							 PosRayElement, CosRayElement);

			ErrorFlag = 0;
			HitBackSide = 0;
			InterceptFlag = 0;
			double PathLength = 0;

			// increment position by tiny amount to get off the element
			// if tracing to the same element
			PosRayElement[0] = PosRayElement[0] + 1.0e-5 * CosRayElement[0];
			PosRayElement[1] = PosRayElement[1] + 1.0e-5 * CosRayElement[1];
			PosRayElement[2] = PosRayElement[2] + 1.0e-5 * CosRayElement[2];

			// {Determine if ray intersects element[j]; if so, Find intersection
			// point with surface of element[j] }
			DetermineElementIntersectionNew(Element,
											PosRayElement,
											CosRayElement,
											PosRaySurfElement,
											CosRaySurfElement,
											DFXYZ,
											&PathLength,
											&ErrorFlag,
											&InterceptFlag,
											&HitBackSide);

			if (InterceptFlag)
			{
				// {If hit multiple elements, this loop determines which one hit
				// first. Also makes sure that correct part of closed surface is
				// hit. Also, handles wavy, but close to flat zernikes and
				// polynomials correctly.}
				if (PathLength < LastPathLength)
				{
					// if (PosRaySurfElement[2] <= Element->ZAperture ||
					// 	Element->SurfaceIndex == 'm' ||
					// 	Element->SurfaceIndex == 'M' ||
					// 	Element->SurfaceIndex == 'r' ||
					// 	Element->SurfaceIndex == 'R')
					// TODO: Is this the correct thing to do?
					if (PosRaySurfElement[2] <= Element->ZAperture)
					{
						StageHit = true;
						LastPathLength = PathLength;
						CopyVec3(LastPosRaySurfElement, PosRaySurfElement);
						CopyVec3(LastCosRaySurfElement, CosRaySurfElement);
						CopyVec3(LastDFXYZ, DFXYZ);
						// LastElementNumber = ((i == 0 && !PT_override) ? Element->element_number : j + 1); // mjw change from j index to element id
						LastElementNumber = Element->element_number;
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

} // namespace SolTrace::NativeRunner

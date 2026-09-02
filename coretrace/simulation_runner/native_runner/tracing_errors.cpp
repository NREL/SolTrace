
#include "tracing_errors.hpp"

#include "simulation_data_export.hpp"

namespace SolTrace::NativeRunner {

	// NOTES: SurfaceNormalErrors() and Errors() differ in that
	// SurfaceNormalErrors() uses the slope error and applies it to the surface
	// normal, while Errors() uses the specularity error, applies it to the ray
	// direction, has a diffuse option, rejects perturbations that pass through
	// an opaque surface, and also handles sun shape errors.
	//
	// Remaining cleanup:
	//     - Split sun shape and surface error handling into separate functions.
	//     - Fix the mrad->rad conversion, which wrongly scales the diffuse case.
	//     - Reduce the random number calls. I.e., sample theta directly rather than thetax
	//       and thetay -> this will break tests because the number of RNG calls will change.

namespace {

// Rotation taking the frame whose +z axis is `axis` back to the reference frame.
glm::dmat3 BuildRayFrame(const glm::dvec3& axis)
{
	glm::dvec3 Euler(0.0, 0.0, 0.0);

	if (axis.z == 0.0)
	{
		if (axis.x == 0.0)
		{
			Euler.x = 0.0;
			Euler.y = PI / 2.0;
		}
		else
		{
			Euler.x = PI / 2.0;
			Euler.y = atan2(axis.y, sqrt(axis.x * axis.x + axis.z * axis.z));
		}
	}
	else
	{
		Euler.x = atan2(axis.x, axis.z);
		Euler.y = atan2(axis.y, sqrt(axis.x * axis.x + axis.z * axis.z));
	}

	Euler.z = 0.0;

	glm::dmat3 RRefToLoc(0.0);
	glm::dmat3 RLocToRef(0.0);
	Data::CalculateTransformMatrices(Euler, RRefToLoc, RLocToRef);

	return RLocToRef;
}

// Draws a uniform azimuth and returns the direction at polar angle `theta` from `axis`.
glm::dvec3 PerturbAboutAxis(MTRand& myrng, const glm::dvec3& axis, double theta)
{
	const glm::dmat3 RLocToRef = BuildRayFrame(axis);

	// phi = atan2(thetay, thetax); //This function appears to  present irregularities that bias results incorrectly for small values of thetay or thetax
	const double phi = myrng() * 2.0 * PI; // Therefore have chosen to randomize phi rather than calculate from randomized theta components
	                                       //  obtained from the distribution. The two approaches are equivalent save for this issue with arctan2.      wendelin 01-12-11

	const glm::dvec3 Origin(0.0, 0.0, 0.0);
	const glm::dvec3 CosLoc(sin(theta) * cos(phi),
	                        sin(theta) * sin(phi),
	                        cos(theta));

	glm::dvec3 PosRef(0.0, 0.0, 0.0);
	glm::dvec3 CosRef(0.0, 0.0, 0.0);

	/*{Transform perturbed ray back to element system}*/
	Data::TransformToReference(Origin, CosLoc, Origin, RLocToRef, PosRef, CosRef);

	return CosRef;
}

} // namespace

void SurfaceNormalErrors(MTRand &myrng,
                         const glm::dvec3 &CosIn,
                         const SolTrace::Data::OpticalPropertySet* OptProperties,
						 const bool LastHitBackSide,
                         glm::dvec3 &CosOut)
{

	/*{Purpose:  To add error terms to the surface normal vector at the surface in question

			   Input - Seed    = Seed for RNG
					   CosIn   = Direction cosine vector of surface normal to which errors will be applied.
					   Element = Element data record
					   DFXYZ   = surface normal vector at interaction point

			   Output - CosOut  = Output direction cosine vector of surface normal after error terms have been included
					   }*/

	const OpticalSide side = LastHitBackSide == false ? OpticalSide::Front : OpticalSide::Back;

	const double delop = OptProperties->get_slope_error(side) / 1000.0;

	double thetax = 0.0, thetay = 0.0, theta2 = 0.0;

	switch (OptProperties->get_error_distribution(side))
	{
	case DistributionType::GAUSSIAN:		// case 'g':
		// gaussian distribution
		thetax = myrng.randNorm(0., delop);
		thetay = myrng.randNorm(0., delop);
		theta2 = thetax * thetax + thetay * thetay;
		break;
	case DistributionType::PILLBOX:			// case 'p':
		// pillbox distribution
		do
		{
			thetax = 2.0 * delop * myrng() - delop;
			thetay = 2.0 * delop * myrng() - delop;
			theta2 = thetax * thetax + thetay * thetay;
		} while (theta2 > (delop * delop));
		break;
	default:
		// TODO: Need an error here.
		break;
	}

	CosOut = PerturbAboutAxis(myrng, CosIn, sqrt(theta2));
}

void Errors(
    MTRand& myrng,
    const glm::dvec3& CosIn,
    int Source,
    TSun* Sun,
    const SolTrace::Data::OpticalPropertySet* OptProperties,
	const bool LastHitBackSide,
    glm::dvec3& CosOut,
    const glm::dvec3& DFXYZ)
{
	/*{Purpose:  To add error terms to the perturbed ray at the surface in question

			   Input - Seed    = Seed for RNG
					   CosIn   = Direction cosine vector of ray to which errors will be applied.
								  If Source below is 1 (i.e. sunshape) then this ray vector is before interaction with element surface
								  If Source below is 2 (i.e. surface error) then this ray vector is after interaction with element surface
									(i.e. reflected ray or transmitted ray)

					   Source  = Source indicator flag
							   = 1 for Sunshape error (Can be gaussian, pillbox or profile data distribution)
							   = 2 for surface errors (Can be gaussian or pillbox distribution)
					   Sun     = Sun data record
					   Element = Element data record
					   DFXYZ   = surface normal vector at interaction point

			   Output - CosOut  = Output direction cosine vector of ray after error terms have been included
					   }*/

	const OpticalSide side = LastHitBackSide == false ? OpticalSide::Front : OpticalSide::Back;

    double delop = 0.0, thetax = 0.0, thetay = 0.0, theta2 = 0.0, theta = 0.0, stest = 0.0;

    unsigned int maxcall = 0;
	// g,p,d
	if (Source == 1)  // sun error
	{
		delop = Sun->Sigma;

		switch (Sun->ShapeIndex)
		{
		case SunShape::GAUSSIAN:			// case 'g':
			thetax = myrng.randNorm(0., delop);
			thetay = myrng.randNorm(0., delop);

			theta2 = thetax * thetax + thetay * thetay;
			break;

		case SunShape::PILLBOX:				// case 'p':
			do
			{
				thetax = 2.0 * delop * myrng() - delop;
				thetay = 2.0 * delop * myrng() - delop;
				theta2 = thetax * thetax + thetay * thetay;
			} while (theta2 > (delop * delop));
			//theta = delop * sqrt(myrng()); // Wang et al. 2010 Solar Energy 195 461-474
			//theta2 = theta * theta;
			break;
		case SunShape::LIMBDARKENED:
			do {
				thetax = 2.0 * Sun->MaxAngle * myrng() - Sun->MaxAngle;
				thetay = 2.0 * Sun->MaxAngle * myrng() - Sun->MaxAngle;
				theta2 = thetax * thetax + thetay * thetay;
				theta = sqrt(theta2);

				stest = 1.0 - 0.5138 * std::pow((theta / Sun->MaxAngle), 4);
			} while ((myrng() > (stest / Sun->MaxIntensity)) || (theta2 > (Sun->MaxAngle * Sun->MaxAngle)));
			break;
		case SunShape::BUIE_CSR:
			// This sun model has long tails so this might take more iterations
			// TODO: add an option to set the max angle (thereby reducing the tail)
			do 
			{
				thetax = 2.0 * Sun->MaxAngle * myrng() - Sun->MaxAngle;
				thetay = 2.0 * Sun->MaxAngle * myrng() - Sun->MaxAngle;
				theta2 = thetax * thetax + thetay * thetay;
				theta = sqrt(theta2);

				if (std::abs(theta) <= 4.65) // within solar disc
					stest = cos(0.326 * theta) / cos(0.308 * theta);
				else // within circumsolar region
					stest = std::exp(Sun->buie_kappa) * std::pow(std::abs(theta), Sun->buie_gamma);

			} while ((myrng() > (stest / Sun->MaxIntensity)) || (theta2 > (Sun->MaxAngle * Sun->MaxAngle)));
			break;
		case SunShape::USER_DEFINED:
			do
			{
				thetax = 2.0 * Sun->MaxAngle * myrng() - Sun->MaxAngle;
				thetay = 2.0 * Sun->MaxAngle * myrng() - Sun->MaxAngle;
				theta2 = thetax * thetax + thetay * thetay;
				theta = sqrt(theta2); // wendelin 1-9-12  do the test once on theta NOT individually on thetax and thetay as before

				size_t idx = 0;
				while (idx < Sun->SunShapeAngle.size() - 1 && Sun->SunShapeAngle[idx] < theta)
					idx++;

				if (idx == 0)
                    stest = Sun->SunShapeIntensity[0];
				else // linear interpolation (switched from average) 12-20-11 wendelin
					stest = Sun->SunShapeIntensity[idx - 1] + (Sun->SunShapeIntensity[idx] - Sun->SunShapeIntensity[idx - 1]) * (theta - Sun->SunShapeAngle[idx - 1]) /
					(Sun->SunShapeAngle[idx] - Sun->SunShapeAngle[idx - 1]);

			} while ((myrng() > (stest / Sun->MaxIntensity)) || (theta2 > (Sun->MaxAngle * Sun->MaxAngle)));
			break;
		default:
			// TODO: This shouldn't throw here...
			throw std::invalid_argument("Unsupported sun shape in Errors function.");
		}
	}

	if (Source == 2)	// surface error
	{
		delop = OptProperties->get_specularity_error(side);

	Label_50:
		switch (OptProperties->get_error_distribution(side))
		{
		case DistributionType::GAUSSIAN:			// case 'g':
			thetax = myrng.randNorm(0., delop);
			thetay = myrng.randNorm(0., delop);

			theta2 = thetax * thetax + thetay * thetay;
			break;

		case DistributionType::PILLBOX:				// case 'p':
			do
			{
				thetax = 2.0 * delop * myrng() - delop;
				thetay = 2.0 * delop * myrng() - delop;
				theta2 = thetax * thetax + thetay * thetay;
			} while (theta2 > (delop * delop));
			break;

		case DistributionType::DIFFUSE:
			theta2 = pow(asin(sqrt(myrng())), 2);
			break;

		default:
			// TODO: Add error message here.
			break;
		}
	}

	theta = sqrt(theta2) / 1.e3; // convert from mrad to rad

	CosOut = PerturbAboutAxis(myrng, CosIn, theta);

    // TODO: Remove goto, should we always do dot product check? // We could move this out of the function and into the caller.

    /*{If reflection error application and new ray direction (after errors) physically goes through opaque surface,
    then go back and get new perturbation 06-12-07}*/		
	if ((Source == 2) &&
		(OptProperties->get_interaction_type() == InteractionType::REFLECTION) &&
        (glm::dot(CosOut, DFXYZ) < 0) &&
		maxcall++ < 50000)
	{
		goto Label_50;
	}
}
// End of Procedure--------------------------------------------------------------

} // namespace SolTrace::NativeRunner


#include "tracing_errors.hpp"

#include "simulation_data_export.hpp"

namespace SolTrace::NativeRunner {

    // TODO: These two functions are basically the same. Refactor to avoid code duplication.
	// NOTES:
	// SurfaceNormalErrors()						Errors()
    //     - uses slope error							- uses specularity error
    //     - applies to surface normal (input)		    - applies to ray direction (input)
    //     - no diffuse option  						- has diffuse option
	//													- has an dot product check (goto) for surface reflection (requires DFXYZ)
    //													- handles sun shape errors
	//      
    // Plan: Break up surface and sun shape error handling into separate functions
	//     - remove the special case of diffuse surfaces before Errors()
    //     - Reduce the random number calls. I.e., sample theta directly rather than thetax and thetay -> this will break tests because the number of RNG calls will change.

void SurfaceNormalErrors(MTRand &myrng,
						 double CosIn[3],
						 const OpticalProperties *OptProperties,
						 double CosOut[3]) noexcept(false) // throw(nanexcept)
{

	/*{Purpose:  To add error terms to the surface normal vector at the surface in question

			   Input - Seed    = Seed for RNG
					   CosIn   = Direction cosine vector of surface normal to which errors will be applied.
					   Element = Element data record
					   DFXYZ   = surface normal vector at interaction point

			   Output - CosOut  = Output direction cosine vector of surface normal after error terms have been included
					   }*/

	int i = 0;
	double Origin[3] = {0.0, 0.0, 0.0},
		   Euler[3] = {0.0, 0.0, 0.0};
	double PosIn[3] = {0.0, 0.0, 0.0},
		   PosOut[3] = {0.0, 0.0, 0.0};
	DistributionType dist;
	double delop = 0.0, thetax = 0.0,
		   thetay = 0.0, theta2 = 0.0,
		   phi = 0.0, theta = 0.0;
	double RRefToLoc[3][3] = {{0.0, 0.0, 0.0},
							  {0.0, 0.0, 0.0},
							  {0.0, 0.0, 0.0}};
	double RLocToRef[3][3] = {{0.0, 0.0, 0.0},
							  {0.0, 0.0, 0.0},
							  {0.0, 0.0, 0.0}};

	if (CosIn[2] == 0.0)
	{
		if (CosIn[0] == 0.0)
		{
			Euler[0] = 0.0;
			Euler[1] = PI / 2.0;
		}
		else
		{
			Euler[0] = PI / 2.0;
			Euler[1] = atan2(CosIn[1], sqrt(CosIn[0] * CosIn[0] + CosIn[2] * CosIn[2]));
		}
	}
	else
	{
		Euler[0] = atan2(CosIn[0], CosIn[2]);
		Euler[1] = atan2(CosIn[1], sqrt(CosIn[0] * CosIn[0] + CosIn[2] * CosIn[2]));
	}

	Euler[2] = 0.0;

	CalculateTransformMatrices(Euler, RRefToLoc, RLocToRef);

	// TODO: Add distribution type to optical properties
	// dist = OptProperties->DistributionType;
	dist = OptProperties->error_distribution_type;
	// delop = OptProperties->RMSSlopeError / 1000.0;
	delop = OptProperties->slope_error / 1000.0;

	switch (dist)
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

	/* {Transform to local coordinate system of ray to set up rotation matrices for coord and inverse
	   transforms} */

	TransformToLocal(PosIn, CosIn, Origin, RRefToLoc, PosOut, CosOut);

	/* {Generate errors in terms of direction cosines in local ray coordinate system} */
	theta = sqrt(theta2);
	// phi = atan2(thetay, thetax); //This function appears to  present irregularities that bias results incorrectly for small values of thetay or thetax
	phi = myrng() * 2.0 * PI; // Therefore have chosen to randomize phi rather than calculate from randomized theta components
	                          //  obtained from the distribution. The two approaches are equivalent save for this issue with arctan2.      wendelin 01-12-11

	CosOut[0] = sin(theta) * cos(phi);
	CosOut[1] = sin(theta) * sin(phi);
	CosOut[2] = cos(theta);

	for (i = 0; i < 3; i++)
	{
		PosIn[i] = PosOut[i];
		CosIn[i] = CosOut[i];
	}

	/*{Transform perturbed ray back to element system}*/
	TransformToReference(PosIn, CosIn, Origin, RLocToRef, PosOut, CosOut);
}

void Errors(
	MTRand &myrng,
	double CosIn[3],
	int Source,
	TSun *Sun,
	// telement_ptr Element,
	const OpticalProperties *OptProperties,
	// TElement *Element,
	// TOpticalProperties *OptProperties,
	double CosOut[3],
	double DFXYZ[3])
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

	double Origin[3] = {0.0, 0.0, 0.0};
	double Euler[3] = {0.0, 0.0, 0.0};
	double PosIn[3] = {0.0, 0.0, 0.0};
	double PosOut[3] = {0.0, 0.0, 0.0};
	double delop = 0, thetax = 0, thetay = 0, theta2 = 0, phi = 0, theta = 0, stest = 0;
	uint_fast64_t i;
	double RRefToLoc[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
	double RLocToRef[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};

	if (CosIn[2] == 0.0)
	{
		if (CosIn[0] == 0.0)
		{
			Euler[0] = 0.0;
			Euler[1] = PI / 2.0;
		}
		else
		{
			Euler[0] = PI / 2.0;
			Euler[1] = atan2(CosIn[1], sqrt(CosIn[0] * CosIn[0] + CosIn[2] * CosIn[2]));
		}
	}
	else
	{
		Euler[0] = atan2(CosIn[0], CosIn[2]);
		Euler[1] = atan2(CosIn[1], sqrt(CosIn[0] * CosIn[0] + CosIn[2] * CosIn[2]));
	}

	Euler[2] = 0.0;

	CalculateTransformMatrices(Euler, RRefToLoc, RLocToRef);

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

				i = 0;
				while (i < Sun->SunShapeAngle.size() - 1 && Sun->SunShapeAngle[i] < theta)
					i++;

				if (i == 0)
					stest = Sun->SunShapeIntensity[0];
				else // linear interpolation (switched from average) 12-20-11 wendelin
					stest = Sun->SunShapeIntensity[i - 1] + (Sun->SunShapeIntensity[i] - Sun->SunShapeIntensity[i - 1]) * (theta - Sun->SunShapeAngle[i - 1]) /
					(Sun->SunShapeAngle[i] - Sun->SunShapeAngle[i - 1]);

			} while ((myrng() > (stest / Sun->MaxIntensity)) || (theta2 > (Sun->MaxAngle * Sun->MaxAngle)));
			break;
		default:
			throw std::invalid_argument("Unsupported sun shape in Errors function.");
		}
	}

	if (Source == 2)	// surface error
	{
		// dist = OptProperties->DistributionType; // errors
		// // delop = sqrt(4.0*sqr(OptProperties->RMSSlopeError)+sqr(OptProperties->RMSSpecError))/1000.0;
		delop = OptProperties->specularity_error;

	Label_50:
		switch (OptProperties->error_distribution_type)
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

	// {Transform to local coordinate system of ray to set up rotation matrices for coordinate and inverse transforms}
	TransformToLocal(PosIn, CosIn, Origin, RRefToLoc, PosOut, CosOut);

	// {Generate errors in terms of direction cosines in local ray coordinate system}
	theta = sqrt(theta2);
	theta = theta / 1.e3; // convert from mrad to rad


	// phi = atan2(thetay, thetax); //This function appears to  present irregularities that bias results incorrectly for small values of thetay or thetax
	phi = myrng() * 2.0 * PI; // Therefore have chosen to randomize phi rather than calculate from randomized theta components
							  //  obtained from the distribution. The two approaches are equivalent save for this issue with arctan2.      wendelin 01-12-11

	CosOut[0] = sin(theta) * cos(phi);
	CosOut[1] = sin(theta) * sin(phi);
	CosOut[2] = cos(theta);

	for (i = 0; i < 3; i++)
	{
		PosIn[i] = PosOut[i];
		CosIn[i] = CosOut[i];
	}

	//{Transform perturbed ray back to element system}
	TransformToReference(PosIn, CosIn, Origin, RLocToRef, PosOut, CosOut);

	// TODO: Remove goto, should we always do dot product check? // We could move this out of the function and into the caller.

	/*{If reflection error application and new ray direction (after errors) physically goes through opaque surface,
    then go back and get new perturbation 06-12-07}*/		
	if ((Source == 2) &&
		(OptProperties->my_type == InteractionType::REFLECTION) &&
		(DOT(CosOut, DFXYZ) < 0) &&
		maxcall++ < 50000)
	{
		goto Label_50;
	}
}
// End of Procedure--------------------------------------------------------------

} // namespace SolTrace::NativeRunner

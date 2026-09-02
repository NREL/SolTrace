
#include "tracing_errors.hpp"

#include "simulation_data_export.hpp"

namespace SolTrace::NativeRunner {

	// NOTES: ApplySlopeError() applies the slope error to the surface normal.
	// ApplySpecularityError() applies the specularity error to the ray
	// direction, has a diffuse option, and rejects perturbations that pass
	// through an opaque surface. ApplySunShape() applies the sun shape to the
	// incoming ray.
	//
	// Remaining cleanup:
	//     - Reduce the random number calls. I.e., sample theta directly rather than thetax
	//       and thetay.
	//     - Validate the sun shape and distribution type during setup so the
	//       sampling paths do not need to throw.

namespace {

constexpr unsigned int kMaxRejectionAttempts = 50000;

constexpr double kMradPerRad = 1000.0;

// Angular radius of the solar disc.
constexpr double kSolarDiscHalfAngleMrad = 4.65;

// Limb darkening profile coefficient.
constexpr double kLimbDarkeningCoeff = 0.5138;

// Buie (2003) solar disc profile, cos(a*theta)/cos(b*theta).
constexpr double kBuieDiscCoeffA = 0.326;
constexpr double kBuieDiscCoeffB = 0.308;

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

// Polar angle of a normally distributed perturbation.
double SampleGaussianAngle(MTRand& myrng, double sigma)
{
	const double thetax = myrng.randNorm(0., sigma);
	const double thetay = myrng.randNorm(0., sigma);
	return sqrt(thetax * thetax + thetay * thetay);
}

// Polar angle of a perturbation drawn uniformly over a disc of radius `half_width`.
double SampleDiscAngle(MTRand& myrng, double half_width)
{
	double theta2 = 0.0;
	do
	{
		const double thetax = 2.0 * half_width * myrng() - half_width;
		const double thetay = 2.0 * half_width * myrng() - half_width;
		theta2 = thetax * thetax + thetay * thetay;
	} while (theta2 > (half_width * half_width));

	return sqrt(theta2);
}

// Polar angle for a radial intensity profile, by rejection sampling a
// (theta_x, theta_y) box out to `max_angle`.
template <typename IntensityFn>
double SampleProfileAngle(MTRand& myrng, double max_angle, double max_intensity,
                          IntensityFn intensity)
{
	double theta2 = 0.0;
	double theta = 0.0;
	double stest = 0.0;

	do
	{
		const double thetax = 2.0 * max_angle * myrng() - max_angle;
		const double thetay = 2.0 * max_angle * myrng() - max_angle;
		theta2 = thetax * thetax + thetay * thetay;
		theta = sqrt(theta2); // wendelin 1-9-12  do the test once on theta NOT individually on thetax and thetay as before

		stest = intensity(theta);

	} while ((myrng() > (stest / max_intensity)) || (theta2 > (max_angle * max_angle)));

	return theta;
}

// Sun shape perturbation angle, in mrad. The sun profiles are tabulated in mrad,
// so sampling happens in those units and the caller converts.
double SampleSunAngleMrad(MTRand& myrng, const TSun& Sun)
{
	switch (Sun.ShapeIndex)
	{
	case SunShape::GAUSSIAN:			// case 'g':
		return SampleGaussianAngle(myrng, Sun.Sigma);

	case SunShape::PILLBOX:				// case 'p':
		//theta = delop * sqrt(myrng()); // Wang et al. 2010 Solar Energy 195 461-474
		return SampleDiscAngle(myrng, Sun.Sigma);

	case SunShape::LIMBDARKENED:
		return SampleProfileAngle(myrng, Sun.MaxAngle, Sun.MaxIntensity,
			[&](double theta) {
				return 1.0 - kLimbDarkeningCoeff * std::pow((theta / Sun.MaxAngle), 4);
			});

	case SunShape::BUIE_CSR:
		// This sun model has long tails so this might take more iterations
		// TODO: add an option to set the max angle (thereby reducing the tail)
		return SampleProfileAngle(myrng, Sun.MaxAngle, Sun.MaxIntensity,
			[&](double theta) {
				if (std::abs(theta) <= kSolarDiscHalfAngleMrad) // within solar disc
					return cos(kBuieDiscCoeffA * theta) / cos(kBuieDiscCoeffB * theta);
				// within circumsolar region
				return std::exp(Sun.buie_kappa) * std::pow(std::abs(theta), Sun.buie_gamma);
			});

	case SunShape::USER_DEFINED:
		return SampleProfileAngle(myrng, Sun.MaxAngle, Sun.MaxIntensity,
			[&](double theta) {
				size_t i = 0;
				while (i < Sun.SunShapeAngle.size() - 1 && Sun.SunShapeAngle[i] < theta)
					i++;

				if (i == 0)
					return Sun.SunShapeIntensity[0];

				// linear interpolation (switched from average) 12-20-11 wendelin
				return Sun.SunShapeIntensity[i - 1]
					+ (Sun.SunShapeIntensity[i] - Sun.SunShapeIntensity[i - 1])
					* (theta - Sun.SunShapeAngle[i - 1])
					/ (Sun.SunShapeAngle[i] - Sun.SunShapeAngle[i - 1]);
			});

	default:
		// TODO: This shouldn't throw here...
		throw std::invalid_argument("Unsupported sun shape.");
	}
}

// Surface error perturbation angle, in radians.
double SampleSurfaceErrorAngle(MTRand& myrng,
                               const SolTrace::Data::OpticalPropertySet& OptProperties,
                               const OpticalSide side)
{
	// delop = sqrt(4.0*sqr(OptProperties->RMSSlopeError)+sqr(OptProperties->RMSSpecError))/1000.0;
	const double delop = OptProperties.get_specularity_error(side) / kMradPerRad;

	switch (OptProperties.get_error_distribution(side))
	{
	case DistributionType::GAUSSIAN:			// case 'g':
		return SampleGaussianAngle(myrng, delop);

	case DistributionType::PILLBOX:				// case 'p':
		return SampleDiscAngle(myrng, delop);

	case DistributionType::DIFFUSE:
		// Gray diffuse (Lambertian) surface: cosine-weighted over the
		// hemisphere, and already in radians.
		return asin(sqrt(myrng()));

	default:
		// TODO: Add error message here.
		return 0.0;
	}
}

} // namespace

glm::dvec3 ApplySlopeError(MTRand& myrng,
                           const glm::dvec3& CosIn,
                           const SolTrace::Data::OpticalPropertySet& OptProperties,
                           const bool LastHitBackSide)
{
	/*{Purpose:  To add error terms to the surface normal vector at the surface in question

			   Input - myrng   = RNG
					   CosIn   = Direction cosine vector of surface normal to which errors
								 will be applied.
					   OptProperties = record of optical properties

			   Returns the surface normal after the slope error has been applied.
					   }*/

	const OpticalSide side = LastHitBackSide == false ? OpticalSide::Front : OpticalSide::Back;

	const double delop = OptProperties.get_slope_error(side) / kMradPerRad;

	double theta = 0.0;

	switch (OptProperties.get_error_distribution(side))
	{
	case DistributionType::GAUSSIAN:		// case 'g':
		theta = SampleGaussianAngle(myrng, delop);
		break;
	case DistributionType::PILLBOX:			// case 'p':
		theta = SampleDiscAngle(myrng, delop);
		break;
	default:
		break;
	}

	return PerturbAboutAxis(myrng, CosIn, theta);
}

glm::dvec3 ApplySunShape(MTRand& myrng, const glm::dvec3& CosIn, const TSun& Sun)
{
	/*{Purpose:  To apply the sun shape to the unperturbed ray at the surface in question

			   Input - myrng   = RNG
					   CosIn   = Direction cosine vector of the ray before interaction
								 with the element surface
					   Sun     = Sun data record

			   Returns the ray direction after the sun shape has been applied.
					   }*/

	const double theta = SampleSunAngleMrad(myrng, Sun) / kMradPerRad;

	return PerturbAboutAxis(myrng, CosIn, theta);
}

glm::dvec3 ApplySpecularityError(MTRand& myrng,
                                 const glm::dvec3& CosIn,
                                 const SolTrace::Data::OpticalPropertySet& OptProperties,
                                 const bool LastHitBackSide,
                                 const glm::dvec3& DFXYZ)
{
	/*{Purpose:  To add error terms to the perturbed ray at the surface in question

			   Input - myrng   = RNG
					   CosIn   = Direction cosine vector of the ray after interaction with
								 the element surface (i.e. reflected or transmitted ray)
					   OptProperties = record of optical properties
					   DFXYZ   = surface normal vector at interaction point

			   Returns the ray direction after error terms have been included.
					   }*/

	const OpticalSide side = LastHitBackSide == false ? OpticalSide::Front : OpticalSide::Back;

	const bool reflecting =
		OptProperties.get_interaction_type() == InteractionType::REFLECTION;

	glm::dvec3 CosOut(0.0, 0.0, 0.0);

	for (unsigned int attempt = 0; attempt < kMaxRejectionAttempts; ++attempt)
	{
		const double theta = SampleSurfaceErrorAngle(myrng, OptProperties, side);

		CosOut = PerturbAboutAxis(myrng, CosIn, theta);

		/*{If reflection error application and new ray direction (after errors) physically goes through opaque surface,
		then go back and get new perturbation 06-12-07}*/
		if (!reflecting || glm::dot(CosOut, DFXYZ) >= 0.0)
			break;
	}

	return CosOut;
}
// End of Procedure--------------------------------------------------------------

} // namespace SolTrace::NativeRunner

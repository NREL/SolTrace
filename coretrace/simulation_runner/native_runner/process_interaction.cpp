#include "process_interaction.hpp"

// SimulationData headers
#include "simulation_data_export.hpp"

// NativeRunner headers
#include "tracing_errors.hpp"

namespace SolTrace::NativeRunner
{

void ProcessInteraction(
    TSystem*                                  System,
    MTRand&                                   myrng,
    const bool                                IncludeSunShape,
    const SolTrace::Data::OpticalPropertySet* optics,
    const bool                                LastHitBackSide,
    const bool                                IncludeErrors,
    // stage info
    const int i,
    const tstage_ptr Stage,
    const uint_fast64_t MultipleHitCount,
    glm::dvec3&         LastDFXYZ,
    // Outputs
    glm::dvec3& LastCosRaySurfElement,
    int&        ErrorFlag,
    glm::dvec3& CosRayOutElement,
    glm::dvec3& LastPosRaySurfElement,
    glm::dvec3& PosRayOutElement)
{
    // Initialize
    glm::dvec3 CosIn(0.0, 0.0, 0.0);
    glm::dvec3 CosOut(0.0, 0.0, 0.0);

    if (!Stage->Virtual)
    {
        // change to account for first hit only in primary stage 8-11-31
        if (IncludeSunShape && i == 0 && MultipleHitCount == 1)
        {
            // Apply sunshape to UNPERTURBED ray at intersection point
            // only apply sunshape error once for primary stage
            CosIn = LastCosRaySurfElement;
            SampleSunShape(myrng, CosIn, &System->Sun, CosOut);
            LastCosRaySurfElement = CosOut;
        }

        //{Determine interaction at surface and direction of perturbed ray}
        ErrorFlag = 0;

        // {Apply surface normal errors to surface normal before interaction
        // ray at intersection point - Wendelin 11-23-09}
        if (IncludeErrors)
        {
            // surface normal errors
            SurfaceNormalErrors(
                myrng, LastDFXYZ, optics, LastHitBackSide, CosOut);
            LastDFXYZ = CosOut;
        }

        Interaction(myrng,
                    LastPosRaySurfElement,
                    LastCosRaySurfElement,
                    LastDFXYZ, // Stage->ElementList[k]->InteractionType,
                    optics,
                    LastHitBackSide,
                    630.0,
                    PosRayOutElement,
                    CosRayOutElement,
                    &ErrorFlag);

        // {Apply specularity optical error to PERTURBED (i.e. after
        // interaction) ray at intersection point}
        if (IncludeErrors)
        {

            const OpticalSide side = LastHitBackSide == false
                                         ? OpticalSide::Front
                                         : OpticalSide::Back;

            CosIn = optics->get_error_distribution(side) ==
                            DistributionType::DIFFUSE
                        ? LastDFXYZ
                        : CosRayOutElement;

            ApplySurfaceError(
                myrng, CosIn, optics, LastHitBackSide, LastDFXYZ, CosOut);
            CosRayOutElement = CosOut;
        }
    }
}

inline double sqr(double x)
{ return (x) * (x); }

void Interaction(MTRand&                                   myrng,
                 const glm::dvec3&                         PosXYZ,
                 const glm::dvec3&                         CosKLM,
                 const glm::dvec3&                         DFXYZ,
                 const SolTrace::Data::OpticalPropertySet* Opticl,
                 const bool                                LastHitBackSide,
                 double                                    Wavelength,
                 glm::dvec3&                               PosOut,
                 glm::dvec3&                               CosOut,
                 int*                                      ErrorFlag)
{
    /* {Purpose: To compute the direction cosines of the ray due to optical
                 interaction at the intersection point of the ray with the
                 surface

        Input - PosXYZ[2] = X,Y,Z coordinates of intersection point.
                DFXYZ     = direction numbers for the surface normal at the
                            intersection point (partial derivatives with
                            respect to X,Y,Z of surface equation)
                InteractionType = Optical interaction type indicator
                          = 1, refraction
                          = 2, reflection
                          = 3, aperture stop
                          = 4, diffraction, transmission grating
                          = 5, diffraction, reflection grating
                CosKLM[2] = direction cosines of incident ray
                Opticl    = record of optical properties
                  .RefractiveIndex[4] = Refractive index of incident and
                                        outgoing medium
                                  [0] = real part of incident medium
                                        refractive index
                                  [1] = imaginary part of ""
                                  [2] = real part of outgoing medium
                                        refractive index
                                  [4] = imaginary part of ""
                  .ApertureStopOrGratingType
                        for InteractionType = 3, aperture stop
                            = 1, slit
                            = 2, elliptical
                        for InteractionType = 4,5 grating
                            = 1, planes parallel to Y-Z plane
                            = 2, concentric cylinders centered about Z-axis
                  .DiffractionOrder = integral order of diffraction for
                                      InteractionTypes = 4,5, grating
                  .AB12[4] = coefficients of polynomial specifying grating
                             spacing for InteractionTypes = 4,5
                       [0] = lower X limit, ApertureStopOrGratingType = 1
                             semi-X axis,   ApertureStopOrGratingType = 2
                       [1] = lower Y limit, ApertureStopOrGratingType = 1
                             semi-Y axis,   ApertureStopOrGratingType = 2
                       [2] = upper X limit, ApertureStopOrGratingType = 1
                             unused,        ApertureStopOrGratingType = 2
                       [4] = upper Y limit, ApertureStopOrGratingType = 1
                             unused,        ApertureStopOrGratingType = 2
                Wavelength = wavelength of ray

        Output - PosOut[2] = position of ray after optical interaction
                 CosOut[2] = direction cosines of ray after optical
                             interaction
                 ErrorFlag = Error flag indicating successful interaction
    } */

    int        i = 0;
    glm::dvec3 CosUVW(0.0, 0.0, 0.0);

    int    NIter = 0, IType = 0, NMord = 0;
    double Epsilon = 0.0, Refr1 = 0.0, Refr2 = 0.0, RMU = 0.0, RM2 = 0.0;
    double D2 = 0.0, B = 0.0, A = 0.0, A2 = 0.0;
    double Gamn = 0.0, Gamn1 = 0.0;
    double X = 0.0, Y = 0.0, A1 = 0.0, B1 = 0.0, Ellips = 0.0, B2 = 0.0;
    double RK = 0.0, RL = 0.0, RM = 0.0, Denom, U = 0, V = 0, W = 0;
    double Varr = 0, GFactr = 0, Rho2 = 0.0, Rho = 0.0, Term = 0.0, G = 0.0,
           D = 0.0, XX = 0.0, Ordiff = 0.0, RLamda = 0.0;
    double Rave = 0.0, Rs = 0.0, Rp = 0.0;

    glm::dvec3 UnitDFXYZ(0.0, 0.0, 0.0);
    double     IncidentAngle = 0.0;

    NIter   = 10;
    Epsilon = 0.000005;

    *ErrorFlag = 0;

    PosOut = PosXYZ;

    switch (Opticl->get_interaction_type())
    {

        /*{  InteractionType = 1, Refraction
        ===============================================================================}*/
    case InteractionType::REFRACTION:
    {
        // TODO: Check that this grabs the correct/savem values
        // as the commented out bit immediately below.
        // Refr1 = Opticl->RefractiveIndex[0];
        // Refr2 = Opticl->RefractiveIndex[2];

        double refrac_front, refrac_back;
        Opticl->get_refraction_indices(refrac_front, refrac_back);

        if (!LastHitBackSide)
        {
            Refr1 = refrac_front;
            Refr2 = refrac_back;
        }
        else
        {
            Refr1 = refrac_back;
            Refr2 = refrac_front;
        }

        RMU = Refr1 / Refr2;
        D2  = glm::dot(DFXYZ, DFXYZ);
        B   = (RMU * RMU - 1.0) / D2;
        A   = RMU * glm::dot(CosKLM, DFXYZ) / D2;
        A2  = A * A;
        if (B > A2) // Total internal reflection
        {
            A      = glm::dot(CosKLM, DFXYZ) / glm::dot(DFXYZ, DFXYZ);
            CosOut = CosKLM - 2.0 * A * DFXYZ;
            return;
        }

        // fresnel equations
        UnitDFXYZ = -glm::normalize(DFXYZ);

        IncidentAngle = acos(glm::dot(CosKLM, UnitDFXYZ));
        Rs = sqr(((Refr1 * cos(IncidentAngle) -
                   Refr2 * sqrt(1 - sqr(Refr1 * sin(IncidentAngle) / Refr2)))) /
                 ((Refr1 * cos(IncidentAngle) +
                   Refr2 * sqrt(1 - sqr(Refr1 * sin(IncidentAngle) / Refr2)))));
        Rp = sqr(((Refr1 * sqrt(1 - sqr(Refr1 * sin(IncidentAngle) / Refr2))) -
                  Refr2 * cos(IncidentAngle)) /
                 ((Refr1 * sqrt(1 - sqr(Refr1 * sin(IncidentAngle) / Refr2))) +
                  Refr2 * cos(IncidentAngle)));
        Rave = (Rp + Rs) / 2.0; // average of s and p polarized light; equal
                                // parts of both = non-polarized
        if (Rave < myrng())     // transmitted through surface
        {
            Gamn = -B / (2.0 * A);

            // Begin Newton-Raphson loop to converge on correct root.
            bool converged = false;
            for (i = 1; i < NIter; i++)
            {
                Gamn1 = (Gamn * Gamn - B) / (2.0 * (Gamn + A));
                if (fabs(Gamn - Gamn1) < Epsilon)
                {
                    converged = true;
                    break;
                }

                Gamn = Gamn1;
            }
            // Failed to converge
            if (converged == false)
            {
                *ErrorFlag = 12;
                return;
            }

            // Have converged on Gamma, Compute direction cosines of refracted
            // ray. Label_Converge:
            CosOut = RMU * CosKLM + Gamn1 * DFXYZ;
            // for (i = 0; i < 3; i++)
            //     CosOut[i] = RMU * CosKLM[i] + Gamn1 * DFXYZ[i];
        }
        else // reflected from surface
        {
            A      = glm::dot(CosKLM, DFXYZ) / glm::dot(DFXYZ, DFXYZ);
            CosOut = CosKLM - 2.0 * A * DFXYZ;
            // for (i = 0; i < 3; i++)
            //     CosOut[i] = CosKLM[i] - 2.0 * A * DFXYZ[i];
        }
        return;
        break;
    }

        /*{  InteractionType = 2, Reflection
        ===============================================================================}*/
    case InteractionType::REFLECTION:
    {
        A = glm::dot(CosKLM, DFXYZ) / glm::dot(DFXYZ, DFXYZ);
        // Compute direction cosines for reflected ray
        CosOut = CosKLM - 2.0 * A * DFXYZ;
        // for (i = 0; i < 3; i++)
        //     CosOut[i] = CosKLM[i] - 2.0 * A * DFXYZ[i];

        return;
        break;
    }

    // 	/*{  InteractionType = 3, Aperture Stop
    // 	===============================================================================}*/
    // case 3:
    // {
    // 	X = PosXYZ[0];
    // 	Y = PosXYZ[1];
    // 	IType = Opticl->ApertureStopOrGratingType;
    // 	A1 = Opticl->AB12[0];
    // 	B1 = Opticl->AB12[1];

    // 	bool ray_missed_aperture = false;
    // 	if (IType == 1) // Slit Aperture
    // 	{
    // 		A2 = Opticl->AB12[2];
    // 		B2 = Opticl->AB12[3];
    // 		if (X < A1 || X > A2)
    // 		{
    // 			*ErrorFlag = 31;
    // 			ray_missed_aperture = true;
    // 		}
    // 		else
    // 		{
    // 			if (Y >= B1 && Y <= B2)
    // 				return;

    // 			*ErrorFlag = 31;
    // 			ray_missed_aperture = true;
    // 		}
    // 	}

    // 	else if (IType == 2) // Elliptical Aperture
    // 	{
    // 		Ellips = X * X / (A1 * A1) + Y * Y / (B1 * B1);
    // 		if (Ellips <= 1.0)
    // 			return;
    // 		*ErrorFlag = 32;
    // 	}

    // 	// RayMissesAperture:
    // 	// Ray misses aperture
    // 	if (ray_missed_aperture == true)
    // 	{
    // 		for (i = 0; i < 3; i++)
    // 			CosOut[i] = 0.0;
    // 	}
    // 	return;

    // 	break;
    // }

    // 	/*{  InteractionType = 4,5; Diffraction
    // 	===============================================================================}*/
    // case 4:
    // case 5:
    // {
    // 	IType = Opticl->ApertureStopOrGratingType;
    // 	NMord = Opticl->DiffractionOrder;
    // 	Refr1 = Opticl->RefractiveIndex[0];
    // 	Refr2 = Opticl->RefractiveIndex[2];
    // 	RMU = Refr1 / Refr2;
    // 	D2 = DOT(DFXYZ, DFXYZ);
    // 	RK = DFXYZ[0];
    // 	RL = DFXYZ[1];
    // 	RM = DFXYZ[2];
    // 	X = PosXYZ[0];
    // 	Y = PosXYZ[1];

    // 	if (IType == 1) // Parallel plane grating
    // 	{
    // 		Denom = RL * RL + RM * RM;
    // 		U = 1.0 / sqrt(1.0 + RK * RK / Denom);
    // 		V = -RK * RL * U / Denom;
    // 		W = -RK * RM * U / Denom;
    // 		Varr = X;
    // 		GFactr = 1.0 / U;
    // 	}

    // 	else if (IType == 2) // Concentric Cylinder Grating
    // 	{
    // 		Rho2 = X * X + Y * Y;
    // 		Rho = sqrt(Rho2);
    // 		RM2 = RM * RM;
    // 		Term = RL * X - RK * Y;
    // 		G = sqrt(D2 * (RM2 * Rho2 + Term * Term));
    // 		U = (RM2 * X + RL * Term) / G;
    // 		V = (RM2 * Y - RK * Term) / G;
    // 		W = -RM * (RK * X + RL * Y) / G;
    // 		Varr = Rho;
    // 		GFactr = Rho / (X * U + Y * V);
    // 	}
    // 	// CompDiffInt:         //Compute interaction due to diffraction
    // 	CosUVW[0] = U;
    // 	CosUVW[1] = V;
    // 	CosUVW[2] = W;

    // 	D = 0.0;
    // 	XX = 1.0;

    // 	for (i = 0; i < 4; i++)
    // 	{
    // 		D = D + Opticl->AB12[i] * XX;
    // 		XX = XX * Varr;
    // 	}

    // 	D = D * GFactr;
    // 	Ordiff = NMord;
    // 	RLamda = Ordiff * Wavelength / (Refr2 * D);
    // 	B = (RMU * RMU - 1.0 + RLamda * RLamda - 2.0 * RMU * RLamda *
    // DOT(CosKLM, CosUVW)) / D2; 	A = RMU * DOT(CosKLM, DFXYZ) / D2;
    // A2 = A * A; 	if (B > A2) // Total internal reflection
    // 	{
    // 		for (i = 0; i < 3; i++)
    // 			CosOut[i] = 0.0;
    // 		*ErrorFlag = 11;
    // 		return;
    // 	}

    // 	Gamn = -B / (2.0 * A);
    // 	if (InteractionType == 5)
    // 		Gamn = -Gamn - 2.0 * A;

    // 	// Begin Newton-Raphson loop to converge on correct root.
    // 	i = 0;
    // 	bool converged = false;
    // 	while (i++ < NIter)
    // 	{
    // 		Gamn1 = (Gamn * Gamn - B) / (2.0 * (Gamn + A));
    // 		if (fabs(Gamn - Gamn1) < Epsilon)
    // 		{
    // 			converged = true;
    // 			break;
    // 		}
    // 		Gamn = Gamn1;
    // 	}
    // 	// Failed to converge
    // 	if (converged == false)
    // 	{
    // 		*ErrorFlag = 12;
    // 		return;
    // 	}
    // 	// Have converged on Gamn1. Compute direction cosines of diffracted ray.
    // 	// CompDCos:
    // 	for (i = 0; i < 3; i++)
    // 		CosOut[i] = RMU * CosKLM[i] - RLamda * CosUVW[i] + Gamn1 *
    // DFXYZ[i];

    // 	break;
    // }
    default: break;
    }
    return;
}

} // namespace SolTrace::NativeRunner

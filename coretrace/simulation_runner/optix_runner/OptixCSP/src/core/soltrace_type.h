#pragma once

#include <ostream>

namespace OptixCSP
{

	enum class ApertureType
	{
		RECTANGLE,
		CIRCLE,
		TRIANGLE,
		QUADRILATERAL,
		HEXAGON,
		ANNULUS
	};

	// types for both scene building and pipeline assembly
	enum class SurfaceType
	{
		FLAT,
		PARABOLIC,
		MESH,
		CYLINDER,
		SPHERICAL
	};

	// mapping of the surface type combined with the aperture type
	// for lookup in the sbt mapping
	struct SurfaceApertureMap
	{
		SurfaceType surfaceType;
		ApertureType apertureType;

		SurfaceApertureMap() : surfaceType(SurfaceType::FLAT),
							   apertureType(ApertureType::RECTANGLE)
		{
		}

		SurfaceApertureMap(SurfaceType surf, ApertureType ap) : surfaceType(surf),
																apertureType(ap)
		{
		}

		// TODO: might not need this, since i always compare
		// surface and aperture types separately ....
		bool operator==(const SurfaceApertureMap &map) const
		{
			return (surfaceType == map.surfaceType) && (apertureType == map.apertureType);
		}

		bool operator<(const SurfaceApertureMap &b) const
		{
			return surfaceType < b.surfaceType ||
				   (surfaceType == b.surfaceType && apertureType < b.apertureType);
		}

		friend std::ostream &operator<<(std::ostream &os, const SurfaceApertureMap &sam);
	};

	inline std::ostream &operator<<(std::ostream &os, const SurfaceApertureMap &sam)
	{
		os << "SurfApMap -- Surface: " << static_cast<int>(sam.surfaceType)
		   << " Aperture: " << static_cast<int>(sam.apertureType);
		return os;
	}
}

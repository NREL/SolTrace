#ifndef SOLTRACE_EMBREE_HELPER_H
#define SOLTRACE_EMBREE_HELPER_H

#include <embree4/rtcore.h>

#include <native_runner_types.hpp>

#include <cstdint>
#include <limits>
#include <string>

namespace SolTrace::EmbreeRunner
{

	struct RayIntersectPayload
	{
		RTCRayQueryContext context; // Embree built-in context (MUST come first)

		int LastHitBackSide = -1;
		double LastDFXYZ[3] = {0.0, 0.0, 0.0};

		// double LastPosRaySurfGlob[3] = { 0.0, 0.0, 0.0 };
		// double LastCosRaySurfGlob[3] = { 0.0, 0.0, 0.0 };

		double LastPosRaySurfStage[3] = {0.0, 0.0, 0.0};
		double LastCosRaySurfStage[3] = {0.0, 0.0, 0.0};

		double LastPosRaySurfElement[3] = {0.0, 0.0, 0.0};
		double LastCosRaySurfElement[3] = {0.0, 0.0, 0.0};

		int_fast64_t element_number = -1;
		int ErrorFlag = 0;

		double PosRayGlobIn[3] = {0.0, 0.0, 0.0};
		double CosRayGlobIn[3] = {0.0, 0.0, 0.0};

		double LastPathLength = std::numeric_limits<double>::infinity();
	};

	inline unsigned int embree_mask(int_fast64_t stage_index)
	{
		return 1u << stage_index;
	}

	void error_function(void *userPtr,
						RTCError error,
						const char *str);

	RTCScene make_scene(RTCDevice &device,
						SolTrace::NativeRunner::TSystem &system);

	bool validate_intersect(double (&LastPosRaySurfElement1)[3],
							double (&LastCosRaySurfElement1)[3],
							double (&LastDFXYZ1)[3],
							uint_fast64_t &LastElementNumber1,
							uint_fast64_t &LastRayNumber1,
							double (&LastPosRaySurfStage1)[3],
							double (&LastCosRaySurfStage1)[3],
							int &ErrorFlag1,
							int &LastHitBackSide1,
							bool &StageHit1,
							double (&LastPosRaySurfElement2)[3],
							double (&LastCosRaySurfElement2)[3],
							double (&LastDFXYZ2)[3],
							uint_fast64_t &LastElementNumber2,
							uint_fast64_t &LastRayNumber2,
							double (&LastPosRaySurfStage2)[3],
							double (&LastCosRaySurfStage2)[3],
							int &ErrorFlag2,
							int &LastHitBackSide2,
							bool &StageHit2);

} // namespace SolTrace::EmbreeRunner

#endif

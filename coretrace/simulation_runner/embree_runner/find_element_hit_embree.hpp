#ifndef SOLTRACE_FIND_ELEMENT_HIT_H
#define SOLTRACE_FIND_ELEMENT_HIT_H

#include <cstdint>

#include <embree4/rtcore.h>

namespace SolTrace::EmbreeRunner
{
    void FindElementHit_embree(
        // Embree args
        const RTCScene &scene,
        // Ray info
        const int i,
        const uint_fast64_t RayNumber,
        const double (&PosRayGlob)[3],
        const double (&CosRayGlob)[3],
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
        bool &StageHit);
} // namespace SolTrace::EmbreeRunner

#endif

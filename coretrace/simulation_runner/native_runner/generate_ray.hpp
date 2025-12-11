#ifndef SOLTRACE_GENERATE_RAY_H
#define SOLTRACE_GENERATE_RAY_H

#include "mtrand.hpp"
#include "native_runner_types.hpp"

namespace SolTrace::NativeRunner {

void GenerateRay(MTRand &myrng,
                 const double PosSunStage[3],
                 double Origin[3],
                 double RLocToRef[3][3],
                 TSun *Sun,
                 double PosRayGlobal[3],
                 double CosRayGlobal[3],
                 double PosRaySun[3]);

} // namespace SolTrace::NativeRunner

#endif

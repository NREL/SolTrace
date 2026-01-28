#ifndef SOLTRACE_SUN_TO_PRIMARY_STAGE_H
#define SOLTRACE_SUN_TO_PRIMARY_STAGE_H

#include "native_runner_types.hpp"
#include "trace_logger.hpp"

namespace SolTrace::NativeRunner
{

    bool SunToPrimaryStage(
        trace_logger_ptr logger,
        TSystem *System,
        TStage *Stage,
        TSun *Sun,
        double PosSunStage[3]);

} // namespace SolTrace::NativeRunner

#endif

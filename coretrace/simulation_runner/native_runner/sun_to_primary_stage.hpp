#ifndef SOLTRACE_SUN_TO_PRIMARY_STAGE_H
#define SOLTRACE_SUN_TO_PRIMARY_STAGE_H

#include "native_runner_types.hpp"
#include "thread_manager.hpp"

namespace SolTrace::NativeRunner
{

    bool SunToPrimaryStage(
        thread_manager_ptr manager,
        TSystem *System,
        TStage *Stage,
        TSun *Sun,
        double PosSunStage[3]);

} // namespace SolTrace::NativeRunner

#endif

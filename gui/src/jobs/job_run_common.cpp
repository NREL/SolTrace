#include "job_run_common.h"
#include "support/math_utility.h"

#include "analysis/rays/ray_volume_raster.h"

#define SECTION(VALUE, TEXT)                                                   \
    promise.setProgressValueAndText(VALUE, TEXT);                              \
    promise.suspendIfRequested();                                              \
    if (promise.isCanceled()) {                                                \
        promise.emplaceResult("Cancelled at " TEXT);                           \
        return;                                                                \
    }

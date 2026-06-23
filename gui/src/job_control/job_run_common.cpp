#include "job_run_common.h"
#include "utilities/math_utility.h"

#include "analysis/ray_volume_raster.h"

#define SECTION(VALUE, TEXT)                                                   \
    promise.setProgressValueAndText(VALUE, TEXT);                              \
    promise.suspendIfRequested();                                              \
    if (promise.isCanceled()) {                                                \
        promise.emplaceResult("Cancelled at " TEXT);                           \
        return;                                                                \
    }

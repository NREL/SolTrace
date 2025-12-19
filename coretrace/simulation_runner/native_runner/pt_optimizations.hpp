
#ifndef SOLTRACE_PT_OPTIMIZATIONS_H
#define SOLTRACE_PT_OPTIMIZATIONS_H

#include "native_runner_types.hpp"
#include "treemesh.hpp"

namespace SolTrace::NativeRunner {

struct eprojdat
{
    TElement *el_addr;
    double d_proj;
    double az;
    double zen;

    eprojdat() {};
    eprojdat(TElement *e, double d, double a, double z)
    {
        el_addr = e;
        d_proj = d;
        az = a;
        zen = z;
    };
};

// Comparison function for sorting vector of eprojdat
static inline bool eprojdat_compare_refactored(const eprojdat &A, const eprojdat &B)
{
    return A.d_proj > B.d_proj;
};

void SetupPTOptimizations(TSystem *System,
                          const bool AsPowerTower,
                          st_hash_tree &sun_hash,
                          st_hash_tree &rec_hash,
                          double (&reccm_helio)[3]);

uint_fast64_t GetPTElements(const bool AsPowerTower,
                            const tstage_ptr Stage,
                            const int i,
                            const bool in_multi_hit_loop, 
                            const double (&PosRayStage)[3],
                            const double (&reccm_helio)[3], 
                            st_hash_tree *rec_hash,
                            const std::vector<void *> &sunint_elements,
                            std::vector<void *> &reflint_elements,
                            bool &has_elements);

} // namespace SolTrace::NativeRunner

#endif

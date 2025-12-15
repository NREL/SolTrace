#ifndef SOLTRACE_PROCESS_INTERACTION_H
#define SOLTRACE_PROCESS_INTERACTION_H

// SimulationData headers
#include "optical_properties.hpp"

// NativeRunner headers
#include "mtrand.hpp"
#include "native_runner_types.hpp"

namespace SolTrace::NativeRunner {

void ProcessInteraction(
    // system info
    TSystem *System,
    MTRand &myrng,
    const bool IncludeSunShape,
    const SolTrace::Data::OpticalProperties *optics,
    const bool IncludeErrors,
    // stage info
    const int i,
    // const TStage *Stage,
    const tstage_ptr Stage,
    // const telement_ptr Elem,
    // const int k,
    // ray info
    const uint_fast64_t MultipleHitCount,
    double (&LastDFXYZ)[3],
    // Outputs
    double (&LastCosRaySurfElement)[3],
    int &ErrorFlag,
    double (&CosRayOutElement)[3],
    double (&LastPosRaySurfElement)[3],
    double (&PosRayOutElement)[3]);

void Interaction(MTRand &myrng,
                 const double PosXYZ[3],
                 const double CosKLM[3],
                 const double DFXYZ[3],
                 const SolTrace::Data::OpticalProperties *Opticl,
                 double Wavelength,
                 double PosOut[3],
                 double CosOut[3],
                 int *ErrorFlag);

} // namespace SolTrace::NativeRunner

#endif

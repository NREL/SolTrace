#ifndef SOLTRACE_PROCESS_INTERACTION_H
#define SOLTRACE_PROCESS_INTERACTION_H

// SimulationData headers
#include "optical_properties.hpp"

// NativeRunner headers
#include "mtrand.hpp"
#include "native_runner_types.hpp"

namespace SolTrace::NativeRunner {

void ProcessInteraction(
    MTRand& myrng,
    const SolTrace::Data::OpticalPropertySet* optics,
    const bool LastHitBackSide,
    const bool IncludeErrors,
    // stage info
    const int i,
    // const TStage *Stage,
    const tstage_ptr Stage,
    // const telement_ptr Elem,
    // const int k,
    // ray info
    const uint_fast64_t MultipleHitCount,
    glm::dvec3& LastDFXYZ,
    // Outputs
    glm::dvec3& LastCosRaySurfElement,
    int& ErrorFlag,
    glm::dvec3& CosRayOutElement,
    glm::dvec3& LastPosRaySurfElement,
    glm::dvec3& PosRayOutElement);

void Interaction(MTRand& myrng,
                 const glm::dvec3& PosXYZ,
                 const glm::dvec3& CosKLM,
                 const glm::dvec3& DFXYZ,
                 const SolTrace::Data::OpticalPropertySet* Opticl,
                 const bool LastHitBackSide,
                 double Wavelength,
                 glm::dvec3& PosOut,
                 glm::dvec3& CosOut,
                 int* ErrorFlag);


} // namespace SolTrace::NativeRunner

#endif

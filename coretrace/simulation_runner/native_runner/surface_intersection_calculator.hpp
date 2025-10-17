#ifndef SOLTRACE_SURFACE_INTERSECTION_CALCULATOR_H
#define SOLTRACE_SURFACE_INTERSECTION_CALCULATOR_H

#include <memory>

#include "aperture.hpp"

namespace SolTrace::NativeRunner {

class SurfaceIntersectionCalculator
{
public:
    SurfaceIntersectionCalculator() {}
    ~SurfaceIntersectionCalculator() {}

    // TODO: Should probably move these to using Vector3d rather
    // than arrays of doubles...
    // TODO: Not sure what purpose CosKLM serves here. Investigate...
    virtual int intersect(const double PosLoc[3],
                          const double CosLoc[3],
                          double PosXYZ[3],
                          double CosKLM[3],
                          double DFXYZ[3],
                          double ZAperture,
                          double *PathLength) = 0;

    virtual double compute_z_aperture(SolTrace::Data::aperture_ptr ap) = 0;

private:
};

using calculator_ptr = typename std::shared_ptr<SurfaceIntersectionCalculator>;

} // namespace SolTrace::NativeRunner

#endif

#ifndef SOLTRACE_FLAT_CALCULATOR_H
#define SOLTRACE_FLAT_CALCULATOR_H

#include "aperture.hpp"
#include "surface.hpp"
#include "surface_intersection_calculator.hpp"

namespace SolTrace::NativeRunner
{

    class FlatCalculator : public SurfaceIntersectionCalculator
    {
    public:
        FlatCalculator(SolTrace::Data::surface_ptr surf,
                       SolTrace::Data::aperture_ptr ap);
        virtual ~FlatCalculator();
        virtual int intersect(const double PosLoc[3],
                              const double CosLoc[3],
                              double PosXYZ[3],
                              double CosKLM[3],
                              double DFXYZ[3],
                              double *PathLength);

        virtual double compute_z_aperture(SolTrace::Data::aperture_ptr ap);

    private:
        SolTrace::Data::aperture_ptr aper;
    };

} // namespace SolTrace::NativeRunner

#endif

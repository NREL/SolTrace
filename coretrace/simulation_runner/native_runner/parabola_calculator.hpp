#ifndef SOLTRACE_PARABOLA_CALCULATOR_H
#define SOLTRACE_PARABOLA_CALCULATOR_H

#include "surface_intersection_calculator.hpp"

#include "aperture.hpp"
#include "surface.hpp"

namespace SolTrace::NativeRunner
{

    class ParabolaCalculator : public SurfaceIntersectionCalculator
    {
    public:
        ParabolaCalculator(SolTrace::Data::surface_ptr surf,
                           SolTrace::Data::aperture_ptr ap);
        virtual ~ParabolaCalculator();
        virtual int intersect(const double PosLoc[3],
                              const double CosLoc[3],
                              double PosXYZ[3],
                              double CosKLM[3],
                              double DFXYZ[3],
                              double *PathLength);

        void surface_normal(const double PosXYZ[3], double DFXYZ[3]);

        virtual double compute_z_aperture(SolTrace::Data::aperture_ptr ap);

    private:
        double cx;
        double cy;

        SolTrace::Data::aperture_ptr aper;
    };

} // namespace SolTrace::NativeRunner

#endif

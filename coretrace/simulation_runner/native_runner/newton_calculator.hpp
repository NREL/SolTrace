#ifndef SOLTRACE_NEWTON_CALCULATOR_H
#define SOLTRACE_NEWTON_CALCULATOR_H

#include <cstdint>

#include "aperture.hpp"
#include "surface_intersection_calculator.hpp"

// #include "surface.hpp"

namespace SolTrace::NativeRunner
{

    class NewtonCalculator : public SurfaceIntersectionCalculator
    {
    public:
        NewtonCalculator(SolTrace::Data::aperture_ptr ap,
                         double tol = 1e-6,
                         uint_fast64_t max_iters = 20);
        virtual ~NewtonCalculator() {}

        virtual int intersect(const double PosLoc[3],
                              const double CosLoc[3],
                              double PosXYZ[3],
                              double CosKLM[3],
                              double DFXYZ[3],
                              double *PathLength);

        // For x = PosXYZ[0], y = PosXYZ[1], make a guess at
        // value of z and place in PosXYZ[2].
        virtual void set_zstart(double PosXYZ[3]) = 0;
        virtual void surface_and_jacobian(const double PosXYZ[3],
                                          double *F,
                                          double DFXYZ[3]) = 0;

        inline double get_tolerance() const
        {
            return this->tolerance;
        }
        inline uint_fast64_t get_max_iters() const
        {
            return this->max_iters;
        }

    private:
        SolTrace::Data::aperture_ptr aper;
        double tolerance;
        uint_fast64_t max_iters;
    };

} // namespace SolTrace::NativeRunner

#endif

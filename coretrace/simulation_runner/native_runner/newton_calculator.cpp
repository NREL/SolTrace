
#include "newton_calculator.hpp"

#include <cmath>
#include <cstdint>

#include "simulation_data_export.hpp"
// #include "vector3d.hpp"

namespace SolTrace::NativeRunner
{

    NewtonCalculator::NewtonCalculator(aperture_ptr ap,
                                       double tol,
                                       uint_fast64_t max_iters)
        : SurfaceIntersectionCalculator(),
          aper(ap),
          tolerance(tol),
          max_iters(max_iters)
    {
    }

    int NewtonCalculator::intersect(const double PosLoc[3],
                                    const double CosLoc[3],
                                    double PosXYZ[3],
                                    double CosKLM[3],
                                    double DFXYZ[3],
                                    double *PathLength)
    {
        int sts = 1;
        uint_fast64_t count = 0;

        // double x0 = PosLoc[0], y0 = PosLoc[1], z0 = PosLoc[2];
        // double mx = CosLoc[0], my = CosLoc[1], mz = CosLoc[2];
        Vector3d u0(PosLoc);
        Vector3d dt(CosLoc);

        double t0 = 0.0;
        double delta = 0.0;
        // double res;
        double fvalue;
        double fprime;
        Vector3d v0(PosLoc);
        Vector3d dv;

        *PathLength = 0.0;
        ZeroVec3(PosXYZ);
        ZeroVec3(CosKLM);
        ZeroVec3(DFXYZ);

        // this->set_zstart(v0.data);
        this->surface_and_jacobian(v0.data, &fvalue, dv.data);
        fprime = dot_product(dv, dt);

        // std::cout << "Tolerance: " << this->tolerance << std::endl;

        while (fabs(fvalue) > this->tolerance && count <= this->max_iters)
        {
            // std::cout << "**** ITER: " << count << " ****"
            //           << "\nt = " << t0
            //           << "\ndelta = " << delta
            //           << "\nf(t) = " << fvalue
            //           << "\nf'(t) = " << fprime
            //           << "\nu0 = " << v0
            //           << std::endl;
            delta = fvalue / fprime;
            t0 -= delta;
            // Updates v0 to (x0, y0, z0) + t0 * (mx, my, mz)
            vector_add(1.0, u0, t0, dt, v0);
            this->surface_and_jacobian(v0.data, &fvalue, dv.data);
            fprime = dot_product(dv, dt);
            ++count;
        }

        if (fabs(fvalue) <= this->tolerance &&
            this->aper->is_in(v0[0], v0[1]))
        {
            sts = 0;
            *PathLength = t0;
            CopyVec3(CosKLM, CosLoc);
            CopyVec3(PosXYZ, v0.data);
            CopyVec3(DFXYZ, dv.data);
        }

        return sts;
    }

} // namespace SolTrace::NativeRunner

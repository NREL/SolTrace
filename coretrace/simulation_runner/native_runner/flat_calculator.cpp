
#include "flat_calculator.hpp"

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <sstream>

#include "matvec.hpp"
#include "simulation_data_export.hpp"
#include "surface.hpp"

namespace SolTrace::NativeRunner
{

    FlatCalculator::FlatCalculator(surface_ptr surf, aperture_ptr ap)
        : SurfaceIntersectionCalculator(),
          aper(ap)
    {
        if (surf == nullptr)
        {
            throw std::invalid_argument("FlatCalculator: Surface pointer cannot be null");
        }

        auto flat = std::dynamic_pointer_cast<Flat>(surf);
        if (flat == nullptr)
        {
            throw std::invalid_argument("FlatCalculator: Surface must be of type Flat");
        }

        if (ap == nullptr)
        {
            throw std::invalid_argument("FlatCalculator: Aperture pointer cannot be null");
        }
    }

    FlatCalculator::~FlatCalculator() {}

    int FlatCalculator::intersect(const double PosLoc[3],
                                  const double CosLoc[3],
                                  double PosXYZ[3],
                                  double CosKLM[3],
                                  double DFXYZ[3],
                                  double *PathLength)
    {
        int sts = 0;
        // TODO: Make sure that intersection should be computed
        // in local coordinate frame.
        double x0 = PosLoc[0];
        double y0 = PosLoc[1];
        double z0 = PosLoc[2];
        double mx = CosLoc[0];
        double my = CosLoc[1];
        double mz = CosLoc[2];
        double t;

        ZeroVec3(PosXYZ);
        ZeroVec3(DFXYZ);
        ZeroVec3(CosKLM);

        if (fabs(mz) < 1e-12)
        {
            // Ray is parallel to the flat plane.
            // It does not intersect.
            sts = 1;
            *PathLength = 0.0;
        }
        else
        {
            t = -z0 / mz;
            x0 += t * mx;
            y0 += t * my;

            if (t > 0.0 && this->aper->is_in(x0, y0))
            {
                CopyVec3(CosKLM, CosLoc);

                PosXYZ[0] = x0;
                PosXYZ[1] = y0;
                PosXYZ[2] = 0.0;

                DFXYZ[0] = DFXYZ[1] = 0.0;
                DFXYZ[2] = 1.0;

                assert(fabs(z0 + t * mz) < 1e-12);

                *PathLength = t;
            }
            else
            {
                sts = 1;
                *PathLength = 0.0;
            }
        }

        return sts;
    }

    double FlatCalculator::compute_z_aperture(aperture_ptr ap)
    {
        return 0.0;
    }

} // namespace SolTrace::NativeRunner

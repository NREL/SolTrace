
#include "parabola_calculator.hpp"

#include <cassert>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <sstream>

#include "matvec.hpp"
#include "simulation_data_export.hpp"
#include "surface.hpp"
// #include "vector3d.hpp"

namespace SolTrace::NativeRunner
{

    ParabolaCalculator::ParabolaCalculator(surface_ptr surf, aperture_ptr ap)
        : aper(ap)
    {
        if (surf == nullptr)
        {
            throw std::invalid_argument("ParabolaCalculator: Surface pointer cannot be null");
        }

        auto para = std::dynamic_pointer_cast<Parabola>(surf);
        if (para == nullptr)
        {
            throw std::invalid_argument("ParabolaCalculator: Surface must be of type Parabola");
        }

        // TODO: Check that this is the correct thing to do
        // this->cx = para->vertex_x_curv;
        // this->cy = para->vertex_y_curv;

        double fx = para->focal_length_x;
        double fy = para->focal_length_y;

        // Validate focal lengths
        if (std::isnan(fx) || std::isnan(fy))
        {
            throw std::invalid_argument("ParabolaCalculator: Focal lengths cannot be NaN");
        }

        if (std::isinf(fx) && std::isinf(fy))
        {
            throw std::invalid_argument("ParabolaCalculator: Both focal lengths cannot be infinite");
        }

        if (ap == nullptr)
        {
            throw std::invalid_argument("ParabolaCalculator: Aperture pointer cannot be null");
        }

        this->cx = fabs(fx) < 1e-12 ? 0.0 : 0.5 / fx;
        this->cy = fabs(fy) < 1e-12 ? 0.0 : 0.5 / fy;

        // std::cout << "cx: " << cx
        //           << "\ncy: " << cy
        //           << std::endl;

        return;
    }

    ParabolaCalculator::~ParabolaCalculator()
    {
    }

    int ParabolaCalculator::intersect(const double PosLoc[3],
                                      const double CosLoc[3],
                                      double PosXYZ[3],
                                      double CosKLM[3],
                                      double DFXYZ[3],
                                      double *PathLength)
    {
        int sts = 0;

        // std::cout << "Computing parabola intersection" << std::endl;

        double x0 = PosLoc[0], y0 = PosLoc[1], z0 = PosLoc[2];
        double mx = CosLoc[0], my = CosLoc[1], mz = CosLoc[2];
        double cx = this->cx, cy = this->cy;
        double t1, t2;
        double a, b, c;
        Vector3d p1, p2;

        c = 0.5 * (cx * x0 * x0 + cy * y0 * y0) - z0;
        b = (x0 * cx * mx + y0 * cy * my - mz);
        a = 0.5 * (cx * mx * mx + cy * my * my);

        ZeroVec3(PosXYZ);
        ZeroVec3(CosKLM);
        ZeroVec3(DFXYZ);

        // std::cout << "a: " << a
        //           << " b: " << b
        //           << " c: " << c
        //           << std::endl;

        if (fabs(a) < 1e-12)
        {
            // This should only happen if mx == my == 0.0
            t1 = -c / b;
            AddVec3(1.0, PosLoc, t1, CosLoc, p1.data);

            if (t1 > 0.0 && this->aper->is_in(p1[0], p1[1]))
            {
                *PathLength = t1;
                // SetVec3(PosXYZ, x0 + t1 * mx, y0 + t1 * my, z0 + t1 * mz);
                // AddVec3(1.0, PosLoc, t1, CosLoc, PosXYZ);
                CopyVec3(PosXYZ, p1.data);
            }
            else
            {
                // Intersection is behind the ray -- same as
                // no solution.
                sts = 1;
                *PathLength = 0.0;
            }
        }
        else
        {
            double scratch = b * b - 4.0 * a * c;
            if (scratch < 0.0)
            {
                // Only imaginary solutions so the ray does
                // NOT intersect the surface
                sts = 1;
                *PathLength = 0.0;
            }
            else
            {
                // TODO: In quadricsurfaceclosedform.cpp checks for
                // PosXYZ[2] > Element->ZAperture. Do we need to do
                // this check here?

                scratch = sqrt(scratch);
                // Negative root
                t1 = -0.5 * (b + scratch) / a;
                // Positive root
                t2 = -0.5 * (b - scratch) / a;

                AddVec3(1.0, PosLoc, t1, CosLoc, p1.data);
                AddVec3(1.0, PosLoc, t2, CosLoc, p2.data);

                // std::cout << "P1: " << p1
                //           << "\nP2: " << p2
                //           << std::endl;

                if (t1 > 0.0 && this->aper->is_in(p1[0], p1[1]))
                {
                    *PathLength = t1;
                    // SetVec3(PosXYZ, x0 + t1 * mx, y0 + t1 * my, z0 + t1 * mz);
                    // AddVec3(1.0, PosLoc, t1, CosLoc, PosXYZ);
                    CopyVec3(PosXYZ, p1.data);
                }
                else if (t2 > 0.0 && this->aper->is_in(p2[0], p2[1]))
                {
                    *PathLength = t2;
                    // SetVec3(PosXYZ, x0 + t2 * mx, y0 + t2 * my, z0 + t2 * mz);
                    // AddVec3(1.0, PosLoc, t2, CosLoc, PosXYZ);
                    CopyVec3(PosXYZ, p2.data);
                }
                else
                {
                    sts = 1;
                    *PathLength = 0.0;
                }
            }
        }

        if (sts == 0)
        {
            this->surface_normal(PosXYZ, DFXYZ);
            CopyVec3(CosKLM, CosLoc);
        }

        // Vector3d grad(DFXYZ);
        // Vector3d pos(PosXYZ);
        // std::cout << "Surface normal: " << grad
        //           << "\nPosition: " << pos
        //           << std::endl;

        // std::cout << "Exiting with status: " << sts << std::endl;

        return sts;
    }

    void ParabolaCalculator::surface_normal(const double PosXYZ[3],
                                            double DFXYZ[3])
    {
        // TODO: Need to default to returning the surface normal of
        // whatever is the "front". Is that the inside of the parabola
        // (which is used currently) or the outside?
        double cx = this->cx;
        double cy = this->cy;
        DFXYZ[0] = -cx * PosXYZ[0];
        DFXYZ[1] = -cy * PosXYZ[1];
        DFXYZ[2] = 1.0;
        return;
    }

    double ParabolaCalculator::compute_z_aperture(aperture_ptr ap)
    {
        // z(x,y) = 1/2 * (cx * x^2 + cy * y^2)
        double r = ap->radius_circumscribed_circle();
        double c = std::max(this->cx, this->cy);
        double zmax = 0.5 * c * r * r;
        return zmax;
    }

} // namespace SolTrace::NativeRunner

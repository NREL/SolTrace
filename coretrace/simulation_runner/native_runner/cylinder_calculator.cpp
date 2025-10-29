
#include "cylinder_calculator.hpp"

#include <cassert>
#include <memory>
#include <stdexcept>
#include <sstream>

#include "aperture.hpp"
#include "matvec.hpp"
#include "simulation_data_export.hpp"
#include "surface.hpp"

namespace SolTrace::NativeRunner
{

    CylinderCalculator::CylinderCalculator(surface_ptr surf, aperture_ptr ap)
        : SurfaceIntersectionCalculator(),
          aper(ap)
    {
        if (surf == nullptr)
        {
            throw std::invalid_argument("CylinderCalculator: Surface pointer cannot be null");
        }

        if (ap == nullptr)
        {
            throw std::invalid_argument("CylinderCalculator: Aperture pointer cannot be null");
        }

        auto cylinder = std::dynamic_pointer_cast<Cylinder>(surf);
        if (cylinder == nullptr)
        {
            throw std::invalid_argument("CylinderCalculator: Surface must be of type Cylinder");
        }

        auto rect = std::dynamic_pointer_cast<Rectangle>(ap);
        if (rect == nullptr)
        {
            throw std::invalid_argument("CylinderCalculator: Aperture must be of type Rectangle");
        }

        // Validate cylinder radius
        if (cylinder->radius <= 0.0)
        {
            std::stringstream ss;
            ss << "CylinderCalculator: Cylinder radius must be positive, got: " << cylinder->radius;
            throw std::invalid_argument(ss.str());
        }

        if (std::isnan(cylinder->radius) || std::isinf(cylinder->radius))
        {
            throw std::invalid_argument("CylinderCalculator: Cylinder radius cannot be NaN or infinite");
        }

        // Validate rectangle dimensions
        if (rect->x_length <= 0.0 || rect->y_length <= 0.0)
        {
            std::stringstream ss;
            ss << "CylinderCalculator: Rectangle dimensions must be positive, got: ("
               << rect->x_length << ", " << rect->y_length << ")";
            throw std::invalid_argument(ss.str());
        }

        // Validate that rectangle x_length matches cylinder diameter
        double expected_x_length = 2.0 * cylinder->radius;
        if (fabs(expected_x_length - rect->x_length) > 1e-8)
        {
            std::stringstream ss;
            ss << "CylinderCalculator: Rectangle x_length (" << rect->x_length
               << ") must equal cylinder diameter (" << expected_x_length << ")";
            throw std::invalid_argument(ss.str());
        }

        this->radius = cylinder->radius;
        // this->length_y = rect->y_length;
        // this->length_x = rect->x_length;
        // this->length_x = 2.0 * this->radius;

        // std::cout << "radius: " << radius
        // << "  length y: " << length_y
        // << "\naperture: (" << rect->x_coord << ", " << rect->y_coord << ") -- ("
        // << rect->x_length << ", " << rect->y_length << ")"
        // << std::endl;

        return;
    }

    CylinderCalculator::~CylinderCalculator()
    {
    }

    int CylinderCalculator::intersect(const double PosLoc[3],
                                      const double CosLoc[3],
                                      double PosXYZ[3],
                                      double CosKLM[3],
                                      double DFXYZ[3],
                                      double *PathLength)
    {
        // std::cout << "Computing cylinder intersection" << std::endl;

        int sts = 0;

        double x0 = PosLoc[0], y0 = PosLoc[1], z0 = PosLoc[2];
        double mx = CosLoc[0], my = CosLoc[1], mz = CosLoc[2];
        double r = this->radius;
        // double yu = 0.5 * this->length_y;
        // double yl = -yu;
        double t1, t2;
        double a, b, c;
        // double y1, y2;

        Vector3d p1, p2;

        c = x0 * x0 + z0 * z0 - 2.0 * z0 * r;
        b = 2.0 * (x0 * mx + (z0 - r) * mz);
        a = mx * mx + mz * mz;

        ZeroVec3(PosXYZ);
        ZeroVec3(CosKLM);
        ZeroVec3(DFXYZ);

        if (fabs(a) < 1e-12)
        {
            // Ray is parallel to the cylinder -- no solution
            sts = 1;
            *PathLength = 0.0;
        }
        else
        {
            double scratch = b * b - 4.0 * a * c;
            if (scratch < 0.0)
            {
                // Only imaginary solutions to quadratic equation -- no solution
                sts = 1;
                *PathLength = 0.0;
            }
            else
            {
                scratch = sqrt(scratch);
                // Negative root
                t1 = -0.5 * (b + scratch) / a;
                // Positive root
                t2 = -0.5 * (b - scratch) / a;

                // y1 = y0 + t1 * my;
                // y2 = y0 + t2 * my;
                AddVec3(1.0, PosLoc, t1, CosLoc, p1.data);
                AddVec3(1.0, PosLoc, t2, CosLoc, p2.data);

                if (t1 > 0.0 && this->aper->is_in(p1[0], p1[1]))
                {
                    // First intersection point is on the surface of
                    // the finite length cylinder
                    *PathLength = t1;
                    // SetVec3(PosXYZ, x0 + t1 * mx, y1, z0 + t1 * mz);
                    // AddVec3(1.0, PosLoc, t1, CosLoc, PosXYZ);
                    // this->surface_normal(PosXYZ, DFXYZ);
                    CopyVec3(PosXYZ, p1.data);
                }
                else if (t2 > 0.0 && this->aper->is_in(p2[0], p2[1]))
                {
                    // First intersection point is outside of surface of finite
                    // length cylinder. Second intersection point is on the
                    // surface of finite cylinder.
                    *PathLength = t2;
                    // SetVec3(PosXYZ, x0 + t2 * mx, y2, z0 + t2 * mz);
                    // AddVec3(1.0, PosLoc, t2, CosLoc, PosXYZ);
                    // this->surface_normal(PosXYZ, DFXYZ);
                    CopyVec3(PosXYZ, p2.data);
                }
                else
                {
                    // Both intersection points are not on the surface of the
                    // finite length cylinder.
                    sts = 1;
                    *PathLength = 0.0;
                }
            }
        }

        if (sts == 0)
        {
            CopyVec3(CosKLM, CosLoc);
            this->surface_normal(PosXYZ, DFXYZ);
        }

        // std::cout << "Exiting with status: " << sts << std::endl;

        return sts;
    }

    void CylinderCalculator::surface_normal(const double PosXYZ[3],
                                            double DFXYZ[3])
    {
        // TODO: Need to default to returning the surface normal of
        // whatever is the "front". Is that the inside of the cylinder
        // (which is used currently) or the outside?
        DFXYZ[0] = 2.0 * PosXYZ[0];
        DFXYZ[1] = 0.0;
        DFXYZ[2] = 2.0 * (PosXYZ[2] - this->radius);
        return;
    }

    double CylinderCalculator::compute_z_aperture(aperture_ptr ap)
    {
        double zmax = 2.0 * this->radius;
        return zmax;
    }

} // namespace SolTrace::NativeRunner

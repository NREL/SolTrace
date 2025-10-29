
#include "sphere_calculator.hpp"

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <sstream>

#include "aperture.hpp"
#include "matvec.hpp"
#include "simulation_data_export.hpp"
#include "surface.hpp"

// #include "vector3d.hpp"

namespace SolTrace::NativeRunner
{

    SphereCalculator::SphereCalculator(surface_ptr surf, aperture_ptr ap)
        : aper(ap)
    {
        if (surf == nullptr)
        {
            throw std::invalid_argument("SphereCalculator: Surface pointer cannot be null");
        }

        auto sph = std::dynamic_pointer_cast<Sphere>(surf);
        if (sph == nullptr)
        {
            throw std::invalid_argument("SphereCalculator: Surface must be of type Sphere");
        }

        // Validate vertex curvature
        if (sph->vertex_curv <= 0.0)
        {
            std::stringstream ss;
            ss << "SphereCalculator: Vertex curvature must be positive, got: " << sph->vertex_curv;
            throw std::invalid_argument(ss.str());
        }

        if (std::isnan(sph->vertex_curv) || std::isinf(sph->vertex_curv))
        {
            throw std::invalid_argument("SphereCalculator: Vertex curvature cannot be NaN or infinite");
        }

        if (ap == nullptr)
        {
            throw std::invalid_argument("SphereCalculator: Aperture pointer cannot be null");
        }

        this->curvature = sph->vertex_curv;
        this->radius = 1.0 / this->curvature;

        // std::cout << "Curvature: " << this->curvature
        //           << "\nRadius: " << this->radius
        //           << std::endl;

        return;
    }

    int SphereCalculator::intersect(const double PosLoc[3],
                                    const double CosLoc[3],
                                    double PosXYZ[3],
                                    double CosKLM[3],
                                    double DFXYZ[3],
                                    double *PathLength)
    {
        // std::cout << "SphereCalculator" << std::endl;

        int sts = 0;

        double x0 = PosLoc[0], y0 = PosLoc[1], z0 = PosLoc[2];
        double alpha = CosLoc[0], beta = CosLoc[1], gamma = CosLoc[2];
        double r = this->radius;
        double t1, t2;
        // double x1, y1, z1, x2, y2, z2;
        double a, b, c;
        double scratch;

        Vector3d p1, p2;

        c = x0 * x0 + y0 * y0 + z0 * z0 - 2.0 * r * z0;
        b = 2.0 * (alpha * x0 + beta * y0 + gamma * (z0 - r));
        a = alpha * alpha + beta * beta + gamma * gamma;

        ZeroVec3(PosXYZ);
        ZeroVec3(CosKLM);
        ZeroVec3(DFXYZ);

        scratch = b * b - 4.0 * a * c;

        // std::cout << "Scratch: " << scratch << std::endl;

        if (scratch < 0.0)
        {
            // No intersection
            sts = 1;
            *PathLength = 0.0;
        }
        else
        {
            scratch = sqrt(scratch);
            t1 = -0.5 * (b + scratch) / a;
            t2 = 0.5 * (scratch - b) / a;

            // z1 = z0 + gamma * t1;
            // z2 = z0 + gamma * t2;
            AddVec3(1.0, PosLoc, t1, CosLoc, p1.data);
            AddVec3(1.0, PosLoc, t2, CosLoc, p2.data);

            // std::cout << "P1: " << p1
            //           << "\nP2: " << p2
            //           << std::endl;

            if (t1 > 0.0 && p1[2] <= r && this->aper->is_in(p1[0], p1[1]))
            {
                // SetVec3(PosXYZ, x0 + t1 * alpha, y0 + t1 * beta, z1);
                // AddVec3(1.0, PosLoc, t1, CosLoc, PosXYZ);
                CopyVec3(PosXYZ, p1.data);
                *PathLength = t1;
            }
            else if (t2 > 0.0 && p2[2] <= r && this->aper->is_in(p2[0], p2[1]))
            {
                // SetVec3(PosXYZ, x0 + t2 * alpha, y0 + t2 * beta, z2);
                // AddVec3(1.0, PosLoc, t2, CosLoc, PosXYZ);
                CopyVec3(PosXYZ, p2.data);
                *PathLength = t2;
            }
            else
            {
                sts = 1;
                *PathLength = 0.0;
            }
        }

        if (sts == 0)
        {
            this->surface_normal(PosXYZ, DFXYZ);
            CopyVec3(CosKLM, CosLoc);
        }

        // if (sts == 0)
        // {
        //     Vector3d p0(PosLoc), v0(CosLoc);
        //     Vector3d p1(PosXYZ), n0(DFXYZ);

        //     std::cout << "Ray Position: " << p0
        //               << "\nRay Direction: " << v0
        //               << "\nIntersection: " << p1
        //               << "\nDistance: " << *PathLength
        //               << "\nNormal: " << n0
        //               << std::endl;
        // }

        return sts;
    }

    void SphereCalculator::surface_normal(const double PosXYZ[3],
                                          double DFXYZ[3])
    {
        double x0 = PosXYZ[0];
        double y0 = PosXYZ[1];
        double z0 = PosXYZ[2];
        double r = this->radius;

        DFXYZ[0] = -2.0 * x0;
        DFXYZ[1] = -2.0 * y0;
        DFXYZ[2] = -2.0 * (z0 - r);

        return;
    }

    double SphereCalculator::compute_z_aperture(aperture_ptr ap)
    {
        double zmax = this->radius;
        double R = this->radius;
        double r = ap->radius_circumscribed_circle();

        // std::cout << "Radius: " << R
        //           << "\nAperture Radius: " << r
        //           << "\nSQRT: " << sqrt(R * R - r * r) << std::endl;

        if (r < this->radius)
        {
            zmax = R - sqrt(R * R - r * r);
        }

        return zmax;
    }

} // namespace SolTrace::NativeRunner

#include <gtest/gtest.h>

#include <cmath>

#include <newton_calculator.hpp>
#include <vector3d.hpp>

#include <common.hpp>

using SolTrace::NativeRunner::NewtonCalculator;

class ParabolaNewton : public NewtonCalculator
{
public:
    ParabolaNewton(double cx,
                   double cy,
                   aperture_ptr ap,
                   double tol,
                   uint_fast64_t max_iters)
        : NewtonCalculator(ap, tol, max_iters),
          cx(cx),
          cy(cy)
    {
    }
    virtual ~ParabolaNewton() {}
    virtual void set_zstart(double PosXYZ[3])
    {
        PosXYZ[2] = 0.0;
        return;
    }
    virtual void surface_and_jacobian(const double PosXYZ[3],
                                      double *F,
                                      double DFXYZ[3])
    {
        double x0 = PosXYZ[0];
        double y0 = PosXYZ[1];
        double z0 = PosXYZ[2];
        *F = z0 - 0.5 * (cx * x0 * x0 + cy * y0 * y0);
        DFXYZ[0] = -cx * x0;
        DFXYZ[1] = -cy * y0;
        DFXYZ[2] = 1.0;
        return;
    }
    double cx;
    double cy;

    // Unused -- here so that this is not a pure virtual class
    virtual double compute_z_aperture(aperture_ptr ap)
    {
        return 0.0;
    }
};

TEST(NewtonCalculator, Case1)
{
    // Case: newton's method converges
    // Ray location
    Vector3d x0(-2.0, 1.0, 0.0);
    // Ray direction
    Vector3d m(2.0, 0.5, 3.0);
    // Intersection point
    const double T = 2.0 / 3.0;
    const double TOL = 1e-12;
    // const double NEWTON_TOL = 1e-12;
    const uint_fast64_t MAX_ITERS = 20;
    // Parabola constants
    const double cx = 1.0;
    const double cy = 2.0;

    // Solution values
    double t;
    Vector3d xt;
    Vector3d mt;
    Vector3d gradf;

    auto rect = create_rectangle_aperture(20.0, 20.0);
    ParabolaNewton pnewt(cx, cy, rect, TOL, MAX_ITERS);
    int sts = pnewt.intersect(x0.data, m.data,
                              xt.data, mt.data, 
                              gradf.data, &t);
    EXPECT_EQ(sts, 0);
    EXPECT_NEAR(t, T, TOL);
    EXPECT_NEAR(xt[0], m[0] * T + x0[0], TOL);
    EXPECT_NEAR(xt[1], m[1] * T + x0[1], TOL);
    EXPECT_NEAR(xt[2], m[2] * T + x0[2], TOL);
    EXPECT_TRUE(is_identical(mt, m));
    EXPECT_NEAR(gradf[0], -cx * (m[0] * T + x0[0]), TOL);
    EXPECT_NEAR(gradf[1], -cy * (m[1] * T + x0[1]), TOL);
    EXPECT_NEAR(gradf[2], 1.0, TOL);
    EXPECT_NEAR(xt[2] - 0.5 * (cx * xt[0] * xt[0] + cy * xt[1] * xt[1]), 0.0, TOL);
}

TEST(NewtonCalculator, Case2)
{
    // Case: newton's method diverges
    Vector3d zero;
    zero.zero();
    // Ray location
    Vector3d x0(-2.0, 1.0, 0.0);
    // Ray direction
    Vector3d m(1.0, 1.0, 1.0);

    // Parabola constants and newton constants
    const double TOL = 1e-12;
    const uint_fast64_t MAX_ITERS = 20;
    const double cx = 1.0;
    const double cy = 2.0;

    // Solution values
    double t;
    Vector3d xt;
    Vector3d mt;
    Vector3d gradf;

    auto rect = create_rectangle_aperture(20.0, 20.0);
    ParabolaNewton pnewt(cx, cy, rect, TOL, MAX_ITERS);
    int sts = pnewt.intersect(x0.data, m.data,
                              xt.data, mt.data, gradf.data, &t);
    EXPECT_EQ(sts, 1);
    EXPECT_EQ(t, 0.0);
    EXPECT_TRUE(is_identical(xt, zero));
    EXPECT_TRUE(is_identical(mt, zero));
    EXPECT_TRUE(is_identical(gradf, zero));
}

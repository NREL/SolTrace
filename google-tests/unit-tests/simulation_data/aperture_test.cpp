
#include <gtest/gtest.h>

#include <aperture.hpp>
#include <constants.hpp>

#include "common.hpp"

// Bounding box test helper functions
std::vector<double> generate_grid(double a, double b, uint_fast64_t npoints)
{
    double h = (b - a) / (npoints - 1);
    std::vector<double> grid;
    grid.resize(npoints);

    for (uint_fast64_t k = 0; k < npoints; ++k)
    {
        grid[k] = k * h + a;
    }

    return grid;
}

bool is_in_box(double x, double y, const double xbox[2], const double ybox[2])
{
    return xbox[0] <= x && x <= xbox[1] &&
           ybox[0] <= y && y <= ybox[1];
}

void bounding_box_test(aperture_ptr ap,
                       double ax, double bx, uint_fast64_t nx,
                       double ay, double by, uint_fast64_t ny,
                       uint_fast64_t &nhit_outside,
                       uint_fast64_t &nhit_inside)
{
    std::vector<double> xgrid = generate_grid(ax, bx, nx);
    std::vector<double> ygrid = generate_grid(ay, by, ny);

    double xbox[2] = {0.0, 0.0};
    double ybox[2] = {0.0, 0.0};
    ap->bounding_box(xbox[0], xbox[1], ybox[0], ybox[1]);

    // std::cout << "Grid bounds: [ax, bx] X [ay, by] = "
    //           << "[" << ax << ", " << bx << "] X ["
    //           << ay << ", " << by << "]"
    //           << std::endl;
    // std::cout << "Bounding box: [ax, bx] X [ay, by] = "
    //           << "[" << xbox[0] << ", " << xbox[1] << "] X ["
    //           << ybox[0] << ", " << ybox[1] << "]"
    //           << std::endl;

    nhit_outside = 0;
    nhit_inside = 0;

    bool hit = false;
    bool in_box = false;
    for (auto xi : xgrid)
    {
        for (auto yj : ygrid)
        {
            // std::cout << "Checking grid point (x,y) = (" << xi << ", " << yj << ")" << std::endl;
            hit = ap->is_in(xi, yj);
            if (hit)
            {
                // std::cout << "Grid point within aperture" << std::endl;
                in_box = is_in_box(xi, yj, xbox, ybox);
                if (in_box)
                {
                    ++nhit_inside;
                }
                else
                {
                    ++nhit_outside;
                }
            }
        }
    }

    return;
}

TEST(Aperture, ApertureBase)
{
    struct TestAperture : public Aperture
    {
        double my_value;
        TestAperture(double mv, ApertureType at) : Aperture(at), my_value(mv)
        {
        }
        virtual ~TestAperture() {}
        virtual double aperture_area() const override { return my_value; }
        virtual double diameter_circumscribed_circle() const override { return 2.0; }
        virtual void bounding_box(double &xmin,
                                  double &xmax,
                                  double &ymin,
                                  double &ymax) const override
        {
            return;
        }
        virtual bool is_in(double x, double y) const override { return false; }
        virtual aperture_ptr make_copy() const override
        {
            return make_aperture<TestAperture>(*this);
        }
        virtual void write_json(nlohmann::ordered_json &jnode) const override
        {
            jnode["my_value"] = my_value;
        }
        virtual std::tuple<std::vector<double>, std::vector<int>> triangulation() const override
        {
            return std::make_tuple(std::vector<double>{}, std::vector<int>{});
        }
    };

    TestAperture ta1(1.2, ApertureType::CIRCLE);
    TestAperture ta2(5.3, ApertureType::RECTANGLE);

    ta1 = ta2;
    EXPECT_EQ(ta1.my_type, ta2.my_type);
    EXPECT_EQ(ta1.my_value, ta2.my_value);

    EXPECT_EQ(ta1.radius_circumscribed_circle(), 1.0);

    auto ta3 = ta1.make_copy();
    EXPECT_EQ(ta3->aperture_area(), ta1.aperture_area());
    EXPECT_EQ(ta3->get_type(), ta1.get_type());
}

TEST(Aperture, Annulus)
{
    constexpr double RO = 5.0;
    const double RI = 1.0;
    const double ARC1 = 90.0;
    const double ARC2 = 360.0;

    // Inside both
    const double X1 = 2.0;
    const double Y1 = 1.0;
    // Outside both in center
    const double X2 = 0.5;
    const double Y2 = 0.5;
    // Outside both
    const double X3 = 1.0;
    const double Y3 = -5.0;
    // Inside ann2 but outside ann1
    const double X4 = -2.0;
    const double Y4 = -3.0;
    // Inside ann2 but outside ann1
    const double X5 = -2.0;
    const double Y5 = 3.0;
    // Inside both
    const double X6 = 2.0;
    const double Y6 = -1.0;

    auto ann1 = make_aperture<Annulus>(RI, RO, ARC1);
    auto ann2 = make_aperture<Annulus>(RI, RO, ARC2);

    EXPECT_EQ(ann1->diameter_circumscribed_circle(), 2 * RO);
    EXPECT_EQ(ann2->diameter_circumscribed_circle(), 2 * RO);

    EXPECT_EQ(ann1->radius_circumscribed_circle(), RO);
    EXPECT_EQ(ann2->radius_circumscribed_circle(), RO);

    EXPECT_EQ(ann2->aperture_area(), PI * (RO * RO - RI * RI));
    EXPECT_EQ(ann1->aperture_area(), 0.25 * ann2->aperture_area());

    EXPECT_TRUE(ann1->is_in(X1, Y1));
    EXPECT_FALSE(ann1->is_in(X2, Y2));
    EXPECT_FALSE(ann1->is_in(X3, Y3));
    EXPECT_FALSE(ann1->is_in(X4, Y4));
    EXPECT_FALSE(ann1->is_in(X5, Y5));
    EXPECT_TRUE(ann1->is_in(X6, Y6));

    EXPECT_TRUE(ann2->is_in(X1, Y1));
    EXPECT_FALSE(ann2->is_in(X2, Y2));
    EXPECT_FALSE(ann2->is_in(X3, Y3));
    EXPECT_TRUE(ann2->is_in(X4, Y4));
    EXPECT_TRUE(ann2->is_in(X5, Y5));
    EXPECT_TRUE(ann2->is_in(X6, Y6));

    aperture_ptr a1 = ann1->make_copy();
    EXPECT_EQ(a1->diameter_circumscribed_circle(),
              ann1->diameter_circumscribed_circle());
    EXPECT_EQ(a1->radius_circumscribed_circle(),
              ann1->radius_circumscribed_circle());
    EXPECT_EQ(a1->aperture_area(), ann1->aperture_area());
    EXPECT_TRUE(a1->is_in(X1, Y1));
    EXPECT_FALSE(a1->is_in(X2, Y2));
    EXPECT_FALSE(a1->is_in(X3, Y3));
    EXPECT_FALSE(a1->is_in(X4, Y4));
    EXPECT_FALSE(a1->is_in(X5, Y5));

    const uint_fast64_t NPOINTS = 1001;
    constexpr double AX = -RO - 1.0;
    constexpr double BX = RO + 1.0;
    constexpr double AY = -RO - 1.0;
    constexpr double BY = RO + 1.0;
    constexpr double GRID_AREA = (BX - AX) * (BY - AY);
    uint_fast64_t nhit_out = 0;
    uint_fast64_t nhit_in = 0;
    bounding_box_test(ann1,
                      AX, BX, NPOINTS,
                      AY, BY, NPOINTS,
                      nhit_out, nhit_in);

    double frac = ann1->aperture_area() / GRID_AREA;

    EXPECT_EQ(nhit_out, 0);
    EXPECT_NEAR((double)(nhit_in) / (NPOINTS * NPOINTS), frac, 1e-3);

    std::cout << "Number of hits outside of bounding box: " << nhit_out << std::endl;
    std::cout << "Number of hits inside bounding box: " << nhit_in << std::endl;
    std::cout << "Expected fraction of hits: " << frac << std::endl;
}

TEST(Aperture, Circle)
{
    constexpr double D = 2.0;
    const double X1 = 0.5;
    const double Y1 = -0.5;
    const double X2 = 1.0;
    const double Y2 = 1.5;
    auto cir = make_aperture<Circle>(D);

    EXPECT_EQ(cir->diameter_circumscribed_circle(), D);
    EXPECT_EQ(cir->radius_circumscribed_circle(), 0.5 * D);
    EXPECT_EQ(cir->aperture_area(), PI * 0.25 * D * D);

    EXPECT_TRUE(cir->is_in(X1, Y1));
    EXPECT_FALSE(cir->is_in(X2, Y2));

    aperture_ptr ap = cir->make_copy();
    EXPECT_EQ(ap->diameter_circumscribed_circle(),
              cir->diameter_circumscribed_circle());
    EXPECT_EQ(ap->radius_circumscribed_circle(),
              cir->radius_circumscribed_circle());
    EXPECT_EQ(ap->aperture_area(), cir->aperture_area());
    EXPECT_TRUE(ap->is_in(X1, Y1));
    EXPECT_FALSE(ap->is_in(X2, Y2));

    const uint_fast64_t NPOINTS = 1001;
    constexpr double AX = -0.5 * D - 1.0;
    constexpr double BX = 0.5 * D + 1.0;
    constexpr double AY = -0.5 * D - 1.0;
    constexpr double BY = 0.5 * D + 1.0;
    constexpr double GRID_AREA = (BX - AX) * (BY - AY);
    uint_fast64_t nhit_out = 0;
    uint_fast64_t nhit_in = 0;
    bounding_box_test(cir,
                      AX, BX, NPOINTS,
                      AY, BY, NPOINTS,
                      nhit_out, nhit_in);

    double frac = cir->aperture_area() / GRID_AREA;

    EXPECT_EQ(nhit_out, 0);
    EXPECT_NEAR((double)(nhit_in) / (NPOINTS * NPOINTS), frac, 1e-3);

    std::cout << "Number of hits outside of bounding box: " << nhit_out << std::endl;
    std::cout << "Number of hits inside bounding box: " << nhit_in << std::endl;
    std::cout << "Expected fraction of hits: " << frac << std::endl;
}

TEST(Aperture, EqualateralTriangle)
{
    const double TOL = 1e-12;
    constexpr double D = 2.0;
    constexpr double R = 0.5 * D;
    const double S = sqrt(3.0) * R; // Side length of triangle
    const double AREA = sqrt(27.0) * R * R / 4.0;

    // Inside inscribed circle
    const double X1 = 0.1;
    const double Y1 = 0.1;
    // Inside but outside inscribed circle on left
    const double X2 = -0.375 * S;
    const double Y2 = -0.375 * R;
    // Outside circumscribed circle
    const double X3 = 1.0;
    const double Y3 = -1.0;
    // Outside but inside circumscribed circle
    const double X4 = -R / sqrt(3.0) - 0.1;
    const double Y4 = 0.1;
    // Inside but outside inscribed circle on right
    const double X5 = 0.375 * S;
    const double Y5 = -0.375 * R;

    auto et = make_aperture<EqualateralTriangle>(D);

    EXPECT_EQ(et->diameter_circumscribed_circle(), D);
    EXPECT_EQ(et->radius_circumscribed_circle(), 0.5 * D);
    EXPECT_NEAR(et->aperture_area(), AREA, TOL);

    EXPECT_TRUE(et->is_in(X1, Y1));
    EXPECT_TRUE(et->is_in(X2, Y2));
    EXPECT_FALSE(et->is_in(X3, Y3));
    EXPECT_FALSE(et->is_in(X4, Y4));
    EXPECT_TRUE(et->is_in(X5, Y5));

    aperture_ptr ap = et->make_copy();
    EXPECT_EQ(ap->diameter_circumscribed_circle(),
              et->diameter_circumscribed_circle());
    EXPECT_EQ(ap->radius_circumscribed_circle(),
              et->radius_circumscribed_circle());
    EXPECT_EQ(ap->aperture_area(), et->aperture_area());
    EXPECT_TRUE(ap->is_in(X1, Y1));
    EXPECT_FALSE(ap->is_in(X3, Y3));

    const uint_fast64_t NPOINTS = 1001;
    constexpr double AX = -R - 1.0;
    constexpr double BX = R + 1.0;
    constexpr double AY = -R - 1.0;
    constexpr double BY = R + 1.0;
    constexpr double GRID_AREA = (BX - AX) * (BY - AY);
    uint_fast64_t nhit_out = 0;
    uint_fast64_t nhit_in = 0;
    bounding_box_test(et,
                      AX, BX, NPOINTS,
                      AY, BY, NPOINTS,
                      nhit_out, nhit_in);

    double frac = et->aperture_area() / GRID_AREA;

    EXPECT_EQ(nhit_out, 0);
    EXPECT_NEAR((double)(nhit_in) / (NPOINTS * NPOINTS), frac, 1e-3);

    std::cout << "Number of hits outside of bounding box: " << nhit_out << std::endl;
    std::cout << "Number of hits inside bounding box: " << nhit_in << std::endl;
    std::cout << "Expected fraction of hits: " << frac << std::endl;
}

TEST(Aperture, Hexagon)
{
    const double TOL = 1e-12;
    const double D = 2.0;
    const double R = 0.5 * D;
    const double S = sqrt(3.0) * R; // Side length of hexagon
    const double AREA = 0.5 * sqrt(27.0) * R * R;

    const double X1 = 1.0;
    const double Y1 = 1.0;
    const double X2 = -0.5;
    const double Y2 = 0.5;
    const double X3 = 0.9;
    const double Y3 = 0.1;
    const double X4 = 0.9;
    const double Y4 = 0.25;
    const double X5 = -0.45;
    const double Y5 = -0.8;
    const double X6 = 0.1;
    const double Y6 = 0.95;

    auto hex = make_aperture<Hexagon>(D);

    EXPECT_EQ(hex->diameter_circumscribed_circle(), D);
    EXPECT_EQ(hex->radius_circumscribed_circle(), R);
    EXPECT_NEAR(hex->aperture_area(), AREA, TOL);

    // Outside Circumscribed, Inside Inscribed
    EXPECT_FALSE(hex->is_in(X1, Y1));
    EXPECT_TRUE(hex->is_in(X2, Y2));
    // Left side inside (outside inscribed circle), outside above, below
    EXPECT_TRUE(hex->is_in(-X3, Y3));
    EXPECT_FALSE(hex->is_in(-X4, Y4));
    EXPECT_FALSE(hex->is_in(-X4, -Y4));
    // Right side inside (outside inscribed circle), outside above, below
    EXPECT_TRUE(hex->is_in(X3, Y3));
    EXPECT_FALSE(hex->is_in(X4, Y4));
    EXPECT_FALSE(hex->is_in(X4, -Y4));
    // Center inside (outside inscribed circle), outside above, below
    EXPECT_TRUE(hex->is_in(X5, Y5));
    EXPECT_FALSE(hex->is_in(X6, Y6));
    EXPECT_FALSE(hex->is_in(-X6, -Y6));

    aperture_ptr ap = hex->make_copy();
    EXPECT_EQ(ap->diameter_circumscribed_circle(),
              hex->diameter_circumscribed_circle());
    EXPECT_EQ(ap->radius_circumscribed_circle(),
              hex->radius_circumscribed_circle());
    EXPECT_EQ(ap->aperture_area(), hex->aperture_area());
    EXPECT_TRUE(ap->is_in(X2, Y2));
    EXPECT_FALSE(ap->is_in(X4, Y4));

    const uint_fast64_t NPOINTS = 1001;
    const double AX = -R - 1.0;
    const double BX = R + 1.0;
    const double AY = -R - 1.0;
    const double BY = R + 1.0;
    const double GRID_AREA = (BX - AX) * (BY - AY);
    uint_fast64_t nhit_out = 0;
    uint_fast64_t nhit_in = 0;
    bounding_box_test(hex,
                      AX, BX, NPOINTS,
                      AY, BY, NPOINTS,
                      nhit_out, nhit_in);

    double frac = hex->aperture_area() / GRID_AREA;

    EXPECT_EQ(nhit_out, 0);
    EXPECT_NEAR((double)(nhit_in) / (NPOINTS * NPOINTS), frac, 1e-3);

    std::cout << "Number of hits outside of bounding box: " << nhit_out << std::endl;
    std::cout << "Number of hits inside bounding box: " << nhit_in << std::endl;
    std::cout << "Expected fraction of hits: " << frac << std::endl;
}

TEST(Aperture, Rectangle)
{
    /**** Common Constants ****/
    const double TOL = 1e-12;
    const double D = 2.0;
    const double LY = 1.0;
    const double LX = sqrt(D * D - LY * LY);
    const double AREA = LY * LX;

    /**** Rectangle that is at the origin ****/
    // Inside
    const double X1 = -0.5 * LX;
    const double Y1 = 0.5 * LY;
    // Outside left
    const double X2 = -2.0 * LX;
    const double Y2 = Y1;
    // Outside right
    const double X3 = 2.0 * LX;
    const double Y3 = -Y1;
    // Outside top
    const double X4 = X1;
    const double Y4 = 1.5 * LY;
    // Outside bottom
    const double X5 = -X1;
    const double Y5 = -1.5 * LY;

    auto rect = make_aperture<Rectangle>(LX, LY);

    EXPECT_NEAR(rect->diameter_circumscribed_circle(), D, TOL);
    EXPECT_NEAR(rect->radius_circumscribed_circle(), 0.5 * D, TOL);
    EXPECT_NEAR(rect->aperture_area(), AREA, TOL);

    EXPECT_TRUE(rect->is_in(X1, Y1));
    EXPECT_FALSE(rect->is_in(X2, Y2));
    EXPECT_FALSE(rect->is_in(X3, Y3));
    EXPECT_FALSE(rect->is_in(X4, Y4));
    EXPECT_FALSE(rect->is_in(X5, Y5));

    aperture_ptr ap = rect->make_copy();
    EXPECT_EQ(ap->diameter_circumscribed_circle(),
              rect->diameter_circumscribed_circle());
    EXPECT_EQ(ap->radius_circumscribed_circle(),
              rect->radius_circumscribed_circle());
    EXPECT_EQ(ap->aperture_area(), rect->aperture_area());
    EXPECT_TRUE(ap->is_in(X1, Y1));
    EXPECT_FALSE(ap->is_in(X3, Y3));

    /**** Rectangle that is shifted from the origin ****/
    const double XL = -1.0;
    const double YL = -2.0;

    // Inside
    const double X6 = 0.5 * LX + XL;
    const double Y6 = 0.5 * LY + YL;
    // Outside left
    const double X7 = -2.0 * LX + XL;
    const double Y7 = Y1 + YL;
    // Outside right
    const double X8 = 2.0 * LX + XL;
    const double Y8 = -Y1 + YL;
    // Outside top
    const double X9 = X1 + XL;
    const double Y9 = 1.5 * LY + YL;
    // Outside bottom
    const double X10 = -X1 + XL;
    const double Y10 = -1.5 * LY + YL;

    auto rect_shift = make_aperture<Rectangle>(LX, LY, XL, YL);

    EXPECT_NEAR(rect_shift->diameter_circumscribed_circle(), D, TOL);
    EXPECT_NEAR(rect_shift->radius_circumscribed_circle(), 0.5 * D, TOL);
    EXPECT_NEAR(rect_shift->aperture_area(), AREA, TOL);

    EXPECT_TRUE(rect_shift->is_in(X6, Y6));
    EXPECT_FALSE(rect_shift->is_in(X7, Y7));
    EXPECT_FALSE(rect_shift->is_in(X8, Y8));
    EXPECT_FALSE(rect_shift->is_in(X9, Y9));
    EXPECT_FALSE(rect_shift->is_in(X10, Y10));

    aperture_ptr ap_shift = rect_shift->make_copy();
    EXPECT_EQ(ap_shift->diameter_circumscribed_circle(),
              rect_shift->diameter_circumscribed_circle());
    EXPECT_EQ(ap_shift->radius_circumscribed_circle(),
              rect_shift->radius_circumscribed_circle());
    EXPECT_EQ(ap_shift->aperture_area(), rect_shift->aperture_area());
    EXPECT_TRUE(ap_shift->is_in(X6, Y6));
    EXPECT_FALSE(ap_shift->is_in(X8, Y8));

    const uint_fast64_t NPOINTS = 1001;
    const double AX = -LX - 1.0;
    const double BX = LX + 1.0;
    const double AY = -LY - 1.0;
    const double BY = LY + 1.0;
    const double GRID_AREA = (BX - AX) * (BY - AY);
    uint_fast64_t nhit_out = 0;
    uint_fast64_t nhit_in = 0;
    bounding_box_test(rect,
                      AX, BX, NPOINTS,
                      AY, BY, NPOINTS,
                      nhit_out, nhit_in);

    double frac = rect->aperture_area() / GRID_AREA;

    EXPECT_EQ(nhit_out, 0);
    EXPECT_NEAR((double)(nhit_in) / (NPOINTS * NPOINTS), frac, 1e-3);

    std::cout << "Number of hits outside of bounding box: " << nhit_out << std::endl;
    std::cout << "Number of hits inside bounding box: " << nhit_in << std::endl;
    std::cout << "Expected fraction of hits: " << frac << std::endl;
}

TEST(Aperture, IrregularTriangle)
{
    const double TOL = 1e-12;
    const double x1 = 0.0, x2 = 1.0, x3 = 2.0 * x2;
    const double y1 = 0.0, y2 = 2.0, y3 = y1;
    auto tri = make_aperture<IrregularTriangle>(x1, y1, x2, y2, x3, y3);

    EXPECT_NEAR(tri->aperture_area(), 0.5 * y2 * (x3 - x1), TOL);

    EXPECT_TRUE(tri->is_in(1.0, 1.0));
    EXPECT_FALSE(tri->is_in(1.5, 2.0));

    auto ap = tri->make_copy();
    EXPECT_NEAR(ap->aperture_area(), 0.5 * y2 * (x3 - x1), TOL);
    EXPECT_TRUE(ap->is_in(1.0, 1.0));
    EXPECT_FALSE(ap->is_in(1.5, 2.0));

    const uint_fast64_t NPOINTS = 1001;
    const double AX = x1 - 1.0;
    const double BX = x3 + 1.0;
    const double AY = y1 - 1.0;
    const double BY = y2 + 1.0;
    const double GRID_AREA = (BX - AX) * (BY - AY);
    uint_fast64_t nhit_out = 0;
    uint_fast64_t nhit_in = 0;
    bounding_box_test(tri,
                      AX, BX, NPOINTS,
                      AY, BY, NPOINTS,
                      nhit_out, nhit_in);

    double frac = tri->aperture_area() / GRID_AREA;

    EXPECT_EQ(nhit_out, 0);
    EXPECT_NEAR((double)(nhit_in) / (NPOINTS * NPOINTS), frac, 1e-3);

    std::cout << "Number of hits outside of bounding box: " << nhit_out << std::endl;
    std::cout << "Number of hits inside bounding box: " << nhit_in << std::endl;
    std::cout << "Expected fraction of hits: " << frac << std::endl;
}

TEST(Aperture, IrregularQuadrilateral)
{
    const double TOL = 1e-12;
    // Parallelogram
    const double x1 = 0.0, x2 = 3.0, x3 = (x2 - x1) + 1.0, x4 = x3 - x2 + x1;
    const double y1 = 0.0, y2 = y1, y3 = 2.0, y4 = y3;
    auto quad = make_aperture<IrregularQuadrilateral>(
        x1, y1, x2, y2, x3, y3, x4, y4);

    // EXPECT_NEAR(quad->aperture_area(), (y2 - y1) * (x3 - x1), TOL);
    EXPECT_NEAR(quad->aperture_area(), (x3 - x4) * (y3 - y1), TOL);

    EXPECT_TRUE(quad->is_in(3.0, 1.0));
    EXPECT_TRUE(quad->is_in(1.0, 1.5));
    EXPECT_FALSE(quad->is_in(4.0, 1.0));

    auto ap = quad->make_copy();
    EXPECT_NEAR(ap->aperture_area(), (x3 - x4) * (y3 - y1), TOL);
    EXPECT_TRUE(ap->is_in(3.0, 1.0));
    EXPECT_TRUE(ap->is_in(1.0, 1.5));
    EXPECT_FALSE(ap->is_in(4.0, 1.0));

    const uint_fast64_t NPOINTS = 1001;
    const double AX = x1 - 1.0;
    const double BX = x2 + 1.0;
    const double AY = y1 - 1.0;
    const double BY = y3 + 1.0;
    const double GRID_AREA = (BX - AX) * (BY - AY);
    uint_fast64_t nhit_out = 0;
    uint_fast64_t nhit_in = 0;
    bounding_box_test(quad,
                      AX, BX, NPOINTS,
                      AY, BY, NPOINTS,
                      nhit_out, nhit_in);

    double frac = quad->aperture_area() / GRID_AREA;

    EXPECT_EQ(nhit_out, 0);
    EXPECT_NEAR((double)(nhit_in) / (NPOINTS * NPOINTS), frac, 1e-3);

    std::cout << "Number of hits outside of bounding box: " << nhit_out << std::endl;
    std::cout << "Number of hits inside bounding box: " << nhit_in << std::endl;
    std::cout << "Expected fraction of hits: " << frac << std::endl;
}

TEST(Aperture, MakeApertureFromType)
{
    const double TOL = 1e-12;

    // Test Circle creation
    std::vector<double> circle_args = {2.0}; // diameter
    auto circle_ap = Aperture::make_aperture_from_type(ApertureType::CIRCLE, circle_args);
    ASSERT_TRUE(circle_ap != nullptr);
    EXPECT_EQ(circle_ap->get_type(), ApertureType::CIRCLE);
    EXPECT_EQ(circle_ap->diameter_circumscribed_circle(), 2.0);
    EXPECT_NEAR(circle_ap->aperture_area(), PI, TOL);

    // Test Rectangle creation
    std::vector<double> rect_args = {3.0, 2.0}; // width, height
    auto rect_ap = Aperture::make_aperture_from_type(ApertureType::RECTANGLE, rect_args);
    ASSERT_TRUE(rect_ap != nullptr);
    EXPECT_EQ(rect_ap->get_type(), ApertureType::RECTANGLE);
    EXPECT_NEAR(rect_ap->aperture_area(), 6.0, TOL);
    EXPECT_TRUE(rect_ap->is_in(0.0, 0.0)); // Center should be inside

    // Test Annulus creation
    std::vector<double> annulus_args = {1.0, 3.0, 180.0}; // inner_radius, outer_radius, arc_angle
    auto annulus_ap = Aperture::make_aperture_from_type(ApertureType::ANNULUS, annulus_args);
    ASSERT_TRUE(annulus_ap != nullptr);
    EXPECT_EQ(annulus_ap->get_type(), ApertureType::ANNULUS);
    EXPECT_EQ(annulus_ap->diameter_circumscribed_circle(), 6.0);
    EXPECT_NEAR(annulus_ap->aperture_area(), 0.5 * PI * (9.0 - 1.0), TOL);

    // Test Hexagon creation
    std::vector<double> hex_args = {4.0}; // circumscribe_diameter
    auto hex_ap = Aperture::make_aperture_from_type(ApertureType::HEXAGON, hex_args);
    ASSERT_TRUE(hex_ap != nullptr);
    EXPECT_EQ(hex_ap->get_type(), ApertureType::HEXAGON);
    EXPECT_EQ(hex_ap->diameter_circumscribed_circle(), 4.0);

    // Test Equilateral Triangle creation
    std::vector<double> tri_args = {2.0}; // circumscribe_diameter
    auto tri_ap = Aperture::make_aperture_from_type(ApertureType::EQUILATERAL_TRIANGLE, tri_args);
    ASSERT_TRUE(tri_ap != nullptr);
    EXPECT_EQ(tri_ap->get_type(), ApertureType::EQUILATERAL_TRIANGLE);
    EXPECT_EQ(tri_ap->diameter_circumscribed_circle(), 2.0);

    // Test Irregular Triangle creation
    std::vector<double> irregular_tri_args = {0.0, 0.0, 1.0, 2.0, 2.0, 0.0}; // x1,y1, x2,y2, x3,y3
    auto irregular_tri_ap = Aperture::make_aperture_from_type(ApertureType::IRREGULAR_TRIANGLE, irregular_tri_args);
    ASSERT_TRUE(irregular_tri_ap != nullptr);
    EXPECT_EQ(irregular_tri_ap->get_type(), ApertureType::IRREGULAR_TRIANGLE);
    EXPECT_NEAR(irregular_tri_ap->aperture_area(), 2.0, TOL);

    // Test Irregular Quadrilateral creation
    std::vector<double> quad_args = {0.0, 0.0, 3.0, 0.0, 4.0, 2.0, 1.0, 2.0}; // x1,y1, x2,y2, x3,y3, x4,y4
    auto quad_ap = Aperture::make_aperture_from_type(ApertureType::IRREGULAR_QUADRILATERAL, quad_args);
    ASSERT_TRUE(quad_ap != nullptr);
    EXPECT_EQ(quad_ap->get_type(), ApertureType::IRREGULAR_QUADRILATERAL);

    // Test insufficient arguments - should return null pointer
    std::vector<double> insufficient_args = {1.0}; // Not enough args for annulus
    auto null_ap = Aperture::make_aperture_from_type(ApertureType::ANNULUS, insufficient_args);
    EXPECT_TRUE(null_ap == nullptr);

    // Test empty arguments
    std::vector<double> empty_args;
    auto null_ap2 = Aperture::make_aperture_from_type(ApertureType::CIRCLE, empty_args);
    EXPECT_TRUE(null_ap2 == nullptr);
}

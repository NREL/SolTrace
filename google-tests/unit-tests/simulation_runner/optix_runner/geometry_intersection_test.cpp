#include <gtest/gtest.h>

#include <constants.hpp>
#include <optix_runner.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result_export.hpp>
#include <simulation_runner.hpp>

#include <functional>
#include <iostream>
#include <limits>

using SolTrace::Data::D2R;
using SolTrace::Runner::RunnerStatus;

const double Z_ELEM = 50.0;
const double Z_BACKSTOP = Z_ELEM - 0.5 * Z_ELEM;
const double TOL = 1e-6;
const uint_fast64_t NRAYS = 10000;

// Adapted from CspElement::set_bounding_box_local.
// Computes the global axis-aligned bounding box for a local bounding box
// transformed by a z-axis rotation (in degrees) and origin offset.
static void compute_global_bounding_box(
    const glm::dvec3 &lower_local,
    const glm::dvec3 &upper_local,
    double zrot_deg, // z-axis rotation in degrees (local-to-global)
    const glm::dvec3 &origin,
    glm::dvec3 &lower_global,
    glm::dvec3 &upper_global)
{
    const double cr = cos(zrot_deg * D2R);
    const double sr = sin(zrot_deg * D2R);
    // glm is column-major: columns are (cr, sr, 0), (-sr, cr, 0), (0, 0, 1)
    const glm::dmat3 rotation(cr, sr, 0.0,
                              -sr, cr, 0.0,
                              0.0, 0.0, 1.0);
    const glm::dvec3 c0 = lower_local;
    const glm::dvec3 c7 = upper_local;
    const glm::dvec3 corners[8] = {
        rotation * glm::dvec3(c0[0], c0[1], c0[2]) + origin,
        rotation * glm::dvec3(c0[0], c0[1], c7[2]) + origin,
        rotation * glm::dvec3(c0[0], c7[1], c0[2]) + origin,
        rotation * glm::dvec3(c7[0], c0[1], c0[2]) + origin,
        rotation * glm::dvec3(c0[0], c7[1], c7[2]) + origin,
        rotation * glm::dvec3(c7[0], c0[1], c7[2]) + origin,
        rotation * glm::dvec3(c7[0], c7[1], c0[2]) + origin,
        rotation * glm::dvec3(c7[0], c7[1], c7[2]) + origin,
    };

    lower_global = glm::dvec3(std::numeric_limits<double>::max());
    upper_global = glm::dvec3(std::numeric_limits<double>::lowest());
    for (const auto &c : corners)
    {
        lower_global = glm::min(lower_global, c);
        upper_global = glm::max(upper_global, c);
    }
}

element_id set_default_sd(SimulationData &sd,
                          surface_ptr surf,
                          aperture_ptr ap,
                          double rotation)
{
    sd.clear();

    // Sun
    auto sun = make_ray_source<Sun>();
    sun->set_position(0, 0, 100);
    sd.add_ray_source(sun);

    // Make target element
    element_ptr el = make_element<SingleElement>();
    el->set_origin(0, 0, Z_ELEM);
    el->set_aim_vector(0, 0, 100); // Face up towards sun
    el->set_surface(surf);
    el->set_aperture(ap);

    OpticalPropertySet optics(InteractionType::REFLECTION, 0, 0);
    optics.set_ideal_absorption(OpticalSide::Both);
    auto opt_ref = sd.add_optical_property_set(optics);
    el->set_optical_property_set(opt_ref);

    el->set_name("el");
    el->set_zrot(rotation);

    // Add element to stage
    element_id id = sd.add_element(el);

    // Back stop element that is bigger than the created element so that the
    // testing element casts a shadow on this big thing.
    element_ptr stop = make_element<SingleElement>();
    double xlb, xub, ylb, yub, zlb, zub;
    ap->bounding_box(xlb, xub, ylb, yub);
    surf->bounding_box(xlb, xub, ylb, yub, zlb, zub);

    // Transform local AABB into global AABB accounting for zrot
    glm::dvec3 lower_global, upper_global;
    compute_global_bounding_box(glm::dvec3(xlb, ylb, zlb),
                                glm::dvec3(xub, yub, zub),
                                rotation,
                                glm::dvec3(0.0, 0.0, Z_ELEM),
                                lower_global, upper_global);

    const double sx = std::max(fabs(lower_global[0]), fabs(upper_global[0])) + 2.0;
    const double sy = std::max(fabs(lower_global[1]), fabs(upper_global[1])) + 2.0;
    stop->set_origin(0, 0, Z_BACKSTOP);
    stop->set_aim_vector(0, 0, 100);
    stop->set_surface(make_surface<Flat>());
    stop->set_aperture(make_aperture<Rectangle>(2.0 * sx, 2.0 * sy));
    stop->set_optical_property_set(opt_ref);
    sd.add_element(stop);

    // Set parameters
    SimulationParameters &params = sd.get_simulation_parameters();
    params.number_of_rays = NRAYS;
    params.max_number_of_rays = params.number_of_rays * 100;
    params.include_optical_errors = false;
    params.include_sun_shape_errors = false;
    params.seed = 123;

    return id;
}

template <typename ApertureT>
static void run_geometry_intersection_test(
    const surface_ptr &surf,
    const std::shared_ptr<ApertureT> &aper,
    double rotation_deg,
    const std::function<double(double, double)> &surface_z,
    bool use_local_coordinates = true)
{
    uint_fast64_t fpos = 0, fneg = 0, hits = 0, misses = 0;

    SimulationData sd;
    element_id test_elid = set_default_sd(sd, surf, aper, rotation_deg);
    SimulationResult result;

    OptixRunner runner;
    RunnerStatus sts = runner.initialize();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.setup_simulation(&sd);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.run_simulation();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.report_simulation(&result, 0);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    ASSERT_EQ(result.get_number_of_records(),
              sd.get_simulation_parameters().number_of_rays);

    const double cos_rot = cos(rotation_deg * D2R);
    const double sin_rot = sin(rotation_deg * D2R);

    for (int i = 0; i < (int)result.get_number_of_records(); ++i)
    {
        auto rr = result[i];
        ASSERT_GE(rr->get_number_of_interactions(), 2);
        glm::dvec3 p0, p1;
        rr->get_position(0, p0);
        rr->get_position(1, p1);
        auto id = rr->get_element(1);
        EXPECT_NEAR(p0[0], p1[0], TOL) << "ray " << i;
        EXPECT_NEAR(p0[1], p1[1], TOL) << "ray " << i;

        const double lx = use_local_coordinates
                              ? p1[0] * cos_rot - p1[1] * sin_rot
                              : p1[0];
        const double ly = use_local_coordinates
                              ? p1[0] * sin_rot + p1[1] * cos_rot
                              : p1[1];

        if (id == test_elid)
        {
            EXPECT_NEAR(p1[2], Z_ELEM + surface_z(lx, ly), TOL * Z_ELEM)
                << "ray " << i;
            EXPECT_TRUE(aper->is_in(lx, ly));
            ++hits;
            if (!aper->is_in(lx, ly))
                ++fpos;
        }
        else
        {
            EXPECT_NEAR(p1[2], Z_BACKSTOP, TOL * Z_ELEM);
            EXPECT_FALSE(aper->is_in(lx, ly));
            ++misses;
            if (aper->is_in(lx, ly))
                ++fneg;
        }
    }

    EXPECT_GT(hits, 0u);
    EXPECT_GT(misses, 0u);
    std::cout << "hits: " << hits << ", misses: " << misses
              << ", false positives: " << fpos
              << ", false negatives: " << fneg << std::endl;
}

template <typename ApertureT>
static void run_flat_geometry_intersection_test(
    const std::shared_ptr<ApertureT> &aper,
    double rotation_deg,
    bool use_local_coordinates = true)
{
    run_geometry_intersection_test(
        make_surface<Flat>(),
        aper,
        rotation_deg,
        [](double, double) { return 0.0; },
        use_local_coordinates);
}

template <typename ApertureT>
static void run_parabolic_geometry_intersection_test(
    const std::shared_ptr<Parabola> &surf,
    const std::shared_ptr<ApertureT> &aper,
    double rotation_deg)
{
    run_geometry_intersection_test(
        surf,
        aper,
        rotation_deg,
        [surf](double lx, double ly)
        {
            return lx * lx / (4.0 * surf->focal_length_x) +
                   ly * ly / (4.0 * surf->focal_length_y);
        });
}

template <typename ApertureT>
static void run_cylindrical_geometry_intersection_test(
    const std::shared_ptr<Cylinder> &surf,
    const std::shared_ptr<ApertureT> &aper,
    double rotation_deg)
{
    run_geometry_intersection_test(
        surf,
        aper,
        rotation_deg,
        [surf](double lx, double)
        {
            return sqrt(surf->radius * surf->radius - lx * lx);
        });
}

template <typename ApertureT>
static void run_spherical_geometry_intersection_test(
    const std::shared_ptr<Sphere> &surf,
    const std::shared_ptr<ApertureT> &aper,
    double rotation_deg)
{
    run_geometry_intersection_test(
        surf,
        aper,
        rotation_deg,
        [surf](double lx, double ly) { return surf->z(lx, ly); });
}

TEST(OptixRunner, FlatRectangle)
{
    const double XL = 10.0;
    const double YL = 5.0;
    const double ROT_DEG = -10.0;
    auto aper = make_aperture<Rectangle>(XL, YL);

    run_flat_geometry_intersection_test(aper, ROT_DEG);
}

TEST(OptixRunner, FlatEquilateralTriangle)
{
    const double d = 4.0;
    const double ROT_DEG = 90.0;
    auto aper = make_aperture<EquilateralTriangle>(d);

    run_flat_geometry_intersection_test(aper, ROT_DEG);
}

TEST(OptixRunner, FlatTriangle)
{
    // Same triangle as FlatTriangle but vertices supplied in CCW order:
    // (0,0)->(2,0)->(1,2). The optix_runner's cross-product sign check sees
    // a positive area and does NOT swap, so the CCW path is taken directly.
    // Results should be identical to FlatTriangle.
    const double x1 = 0.0, x2 = 2.0, x3 = 1.0;
    const double y1 = 0.0, y2 = 0.0, y3 = 2.0;
    const double ROT_DEG = 110.0;
    auto aper = make_aperture<IrregularTriangle>(x1, y1, x2, y2, x3, y3);

    run_flat_geometry_intersection_test(aper, ROT_DEG);
}

TEST(OptixRunner, FlatTriangle_CW)
{
    // Vertices are in CW order (signed area = -2). The optix_runner's
    // cross-product sign check detects this and swaps p1<->p2 to produce
    // CCW winding before building the geometry.
    const double x1 = 0.0, x2 = 1.0, x3 = 2.0 * x2;
    const double y1 = 0.0, y2 = 2.0, y3 = y1;
    const double ROT_DEG = 110.0;
    auto aper = make_aperture<IrregularTriangle>(x1, y1, x2, y2, x3, y3);

    run_flat_geometry_intersection_test(aper, ROT_DEG);
}

TEST(OptixRunner, FlatQuadrilateral)
{
    // Parallelogram
    const double x1 = 0.0, x2 = 3.0, x3 = (x2 - x1) + 1.0, x4 = x3 - x2 + x1;
    const double y1 = 0.0, y2 = y1, y3 = 2.0, y4 = y3;
    const double ROT_DEG = -45.0;
    auto aper = make_aperture<IrregularQuadrilateral>(
        x1, y1, x2, y2, x3, y3, x4, y4);

    run_flat_geometry_intersection_test(aper, ROT_DEG);
}

TEST(OptixRunner, FlatQuadrilateral_CW)
{
    // Same parallelogram as FlatQuadrilateral but vertices supplied in
    // clockwise winding order: (0,0)->(1,2)->(4,2)->(3,0).
    // The x1-x3 diagonal is already interior for this ordering, so
    // ensure_valid_diagonal() leaves the vertices unchanged (CW).
    // The optix_runner's shoelace sign check detects the negative area and
    // swaps p1<->p3 to restore CCW winding before building the geometry.
    // Results should be identical to FlatQuadrilateral.
    const double x1 = 0.0, x2 = 1.0, x3 = 4.0, x4 = 3.0;
    const double y1 = 0.0, y2 = 2.0, y3 = 2.0, y4 = 0.0;
    const double ROT_DEG = -45.0;
    auto aper = make_aperture<IrregularQuadrilateral>(
        x1, y1, x2, y2, x3, y3, x4, y4);

    run_flat_geometry_intersection_test(aper, ROT_DEG);
}

TEST(OptixRunner, ParabolaRectangle)
{
    constexpr double CX = 0.5;
    constexpr double CY = 1.0;
    constexpr double FX = 0.5 / CX;
    constexpr double FY = 0.5 / CY;
    const double XL = 10.0, YL = 5.0;
    const double ROT_DEG = -135.0;
    auto surf = make_surface<Parabola>(FX, FY);
    auto aper = make_aperture<Rectangle>(XL, YL);

    run_parabolic_geometry_intersection_test(surf, aper, ROT_DEG);
}

TEST(OptixRunner, Cylinder)
{
    const double R = 5.0;
    const double YL = 3.0; // Total cylinder length
    const double ROT_DEG = 25.0;
    auto surf = make_surface<Cylinder>(R);
    auto aper = make_aperture<Rectangle>(2 * R, YL);

    run_cylindrical_geometry_intersection_test(surf, aper, ROT_DEG);
}

TEST(OptixRunner, FlatCircle)
{
    const double R = 5.0;
    const double ROT_DEG = 10.0; // Should make no difference
    auto aper = make_aperture<Circle>(2 * R);

    run_flat_geometry_intersection_test(aper, ROT_DEG, false);
}

TEST(OptixRunner, FlatHexagon)
{
    const double S = 5.0;
    const double ROT_DEG = 30.0;
    auto aper = make_aperture<Hexagon>(2 * S);

    run_flat_geometry_intersection_test(aper, ROT_DEG);
}

TEST(OptixRunner, FlatAnnulus_FullArc)
{
    const double R0 = 5.0;
    const double R1 = 180.0;
    const double ARC = 360.0;
    const double ROT_DEG = -15.0; // Should make no difference
    auto aper = make_aperture<Annulus>(R0, R1, ARC);

    run_flat_geometry_intersection_test(aper, ROT_DEG, false);
}

TEST(OptixRunner, FlatAnnulus_PartialArc)
{
    const double R0 = 5.0;
    const double R1 = 180.0;
    const double ARC = 90.0;
    const double ROT_DEG = -15.0;
    auto aper = make_aperture<Annulus>(R0, R1, ARC);

    run_flat_geometry_intersection_test(aper, ROT_DEG);
}

TEST(OptixRunner, ParabolicHexagon)
{
    const double S = 5.0;
    const double FOCAL_X = 10.0;
    const double FOCAL_Y = 30.0;
    const double ROT_DEG = 30.0;
    auto surf = make_surface<Parabola>(FOCAL_X, FOCAL_Y);
    auto aper = make_aperture<Hexagon>(2 * S);

    run_parabolic_geometry_intersection_test(surf, aper, ROT_DEG);
}

TEST(OptixRunner, ParabolicTriangle)
{
    const double d = 4.0;
    const double FOCAL_X = 10.0;
    const double FOCAL_Y = 30.0;
    const double ROT_DEG = 90.0;
    auto surf = make_surface<Parabola>(FOCAL_X, FOCAL_Y);
    auto aper = make_aperture<EquilateralTriangle>(d);

    run_parabolic_geometry_intersection_test(surf, aper, ROT_DEG);
}

TEST(OptixRunner, ParabolicAnnulus)
{
    const double R0 = 2.0;
    const double R1 = 5.0;
    const double ARC = 180.0;
    const double FOCAL_X = 10.0;
    const double FOCAL_Y = 30.0;
    const double ROT_DEG = -15.0;
    auto surf = make_surface<Parabola>(FOCAL_X, FOCAL_Y);
    auto aper = make_aperture<Annulus>(R0, R1, ARC);

    run_parabolic_geometry_intersection_test(surf, aper, ROT_DEG);
}

TEST(OptixRunner, ParabolicQuadrilateral)
{
    // Parallelogram
    const double x1 = 0.0, x2 = 3.0, x3 = (x2 - x1) + 1.0, x4 = x3 - x2 + x1;
    const double y1 = 0.0, y2 = y1, y3 = 2.0, y4 = y3;
    const double FOCAL_X = 10.0;
    const double FOCAL_Y = 30.0;
    const double ROT_DEG = -45.0;
    auto surf = make_surface<Parabola>(FOCAL_X, FOCAL_Y);
    auto aper = make_aperture<IrregularQuadrilateral>(x1, y1, x2, y2, x3, y3, x4, y4);

    run_parabolic_geometry_intersection_test(surf, aper, ROT_DEG);
}

TEST(OptixRunner, ParabolicCircle)
{
    const double R = 5.0;
    const double FOCAL_X = 10.0;
    const double FOCAL_Y = 30.0;
    const double ROT_DEG = -45.0;
    auto surf = make_surface<Parabola>(FOCAL_X, FOCAL_Y);
    auto aper = make_aperture<Circle>(2 * R);

    run_parabolic_geometry_intersection_test(surf, aper, ROT_DEG);
}

TEST(OptixRunner, SphericalRectangle)
{
    const double XL = 10.0;
    const double YL = 5.0;
    const double SPHERE_R = 20.0;
    const double ROT_DEG = -10.0;
    auto surf = make_surface<Sphere>(1.0 / SPHERE_R);
    auto aper = make_aperture<Rectangle>(XL, YL);

    run_spherical_geometry_intersection_test(surf, aper, ROT_DEG);
}

TEST(OptixRunner, SphericalTriangle)
{
    const double x1 = 0.0, x2 = 2.0, x3 = 1.0;
    const double y1 = 0.0, y2 = 0.0, y3 = 2.0;
    const double SPHERE_R = 20.0;
    const double ROT_DEG = 110.0;
    auto surf = make_surface<Sphere>(1.0 / SPHERE_R);
    auto aper = make_aperture<IrregularTriangle>(x1, y1, x2, y2, x3, y3);

    run_spherical_geometry_intersection_test(surf, aper, ROT_DEG);
}

TEST(OptixRunner, SphericalQuadrilateral)
{
    const double x1 = 0.0, x2 = 3.0, x3 = (x2 - x1) + 1.0, x4 = x3 - x2 + x1;
    const double y1 = 0.0, y2 = y1, y3 = 2.0, y4 = y3;
    const double SPHERE_R = 20.0;
    const double ROT_DEG = -45.0;
    auto surf = make_surface<Sphere>(1.0 / SPHERE_R);
    auto aper = make_aperture<IrregularQuadrilateral>(
        x1, y1, x2, y2, x3, y3, x4, y4);

    run_spherical_geometry_intersection_test(surf, aper, ROT_DEG);
}

TEST(OptixRunner, SphericalCircle)
{
    const double R = 5.0;
    const double SPHERE_R = 20.0;
    const double ROT_DEG = 10.0;
    auto surf = make_surface<Sphere>(1.0 / SPHERE_R);
    auto aper = make_aperture<Circle>(2 * R);

    run_spherical_geometry_intersection_test(surf, aper, ROT_DEG);
}

TEST(OptixRunner, SphericalHexagon)
{
    const double S = 5.0;
    const double SPHERE_R = 20.0;
    const double ROT_DEG = 30.0;
    auto surf = make_surface<Sphere>(1.0 / SPHERE_R);
    auto aper = make_aperture<Hexagon>(2 * S);

    run_spherical_geometry_intersection_test(surf, aper, ROT_DEG);
}

TEST(OptixRunner, SphericalAnnulus)
{
    const double R0 = 2.0;
    const double R1 = 5.0;
    const double ARC = 180.0;
    const double SPHERE_R = 20.0;
    const double ROT_DEG = -15.0;
    auto surf = make_surface<Sphere>(1.0 / SPHERE_R);
    auto aper = make_aperture<Annulus>(R0, R1, ARC);

    run_spherical_geometry_intersection_test(surf, aper, ROT_DEG);
}

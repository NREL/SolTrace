#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include <determine_interaction_type.hpp>
#include <mtrand.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result_export.hpp>
#include <trace_logger.hpp>

using SolTrace::NativeRunner::determine_interaction_type;
using SolTrace::NativeRunner::make_trace_logger;
using SolTrace::NativeRunner::MTRand;

namespace
{
    // Fixed local surface normal; LastCosRaySurfElement is built below to give
    // an exact, known incidence angle against this normal.
    glm::dvec3 surface_normal()
    {
        return glm::dvec3(0.0, 0.0, 1.0);
    }

    // Ray direction (in the same local frame as surface_normal()) that hits
    // the front side at exactly incidence_angle_rad.
    glm::dvec3 ray_direction_for_angle(double incidence_angle_rad)
    {
        return glm::dvec3(std::sin(incidence_angle_rad), 0.0, -std::cos(incidence_angle_rad));
    }
}

TEST(DetermineInteractionType, ZeroReflectivityTableValueAlwaysAbsorbs)
{
    OpticalPropertySet optics(InteractionType::REFLECTION, "zero_reflectivity_table");
    optics.set_errors(OpticalSide::Both, DistributionType::NONE, 0.0, 0.0);

    constexpr double angle_mrad = 500.0;
    std::vector<AngularTablePoint> table = {{0.0, 1.0}, {angle_mrad, 0.0}, {1000.0, 1.0}};
    optics.set_reflectivity_table(OpticalSide::Front, table);

    MTRand myrng(1);
    auto logger = make_trace_logger();
    RayEvent rev;

    bool good = determine_interaction_type(logger, 0, 0, myrng, &optics,
                                           surface_normal(),
                                           ray_direction_for_angle(angle_mrad * 1.0e-3),
                                           /*LastHitBackSide=*/false, rev);

    EXPECT_TRUE(good);
    EXPECT_EQ(rev, RayEvent::ABSORB);
}

TEST(DetermineInteractionType, OneReflectivityTableValueNeverAbsorbs)
{
    OpticalPropertySet optics(InteractionType::REFLECTION, "one_reflectivity_table");
    optics.set_errors(OpticalSide::Both, DistributionType::NONE, 0.0, 0.0);

    constexpr double angle_mrad = 500.0;
    std::vector<AngularTablePoint> table = {{0.0, 0.0}, {angle_mrad, 1.0}, {1000.0, 0.0}};
    optics.set_reflectivity_table(OpticalSide::Front, table);

    MTRand myrng(1);
    auto logger = make_trace_logger();
    RayEvent rev;

    bool good = determine_interaction_type(logger, 0, 0, myrng, &optics,
                                           surface_normal(),
                                           ray_direction_for_angle(angle_mrad * 1.0e-3),
                                           /*LastHitBackSide=*/false, rev);

    EXPECT_TRUE(good);
    EXPECT_EQ(rev, RayEvent::REFLECT);
}

TEST(DetermineInteractionType, ZeroTransmissivityTableValueAlwaysAbsorbs)
{
    OpticalPropertySet optics(InteractionType::REFRACTION, "zero_transmissivity_table");
    optics.set_errors(OpticalSide::Both, DistributionType::NONE, 0.0, 0.0);

    constexpr double angle_mrad = 200.0;
    std::vector<AngularTablePoint> table = {{0.0, 1.0}, {angle_mrad, 0.0}, {1000.0, 1.0}};
    optics.set_transmissivity_table(OpticalSide::Front, table);

    MTRand myrng(1);
    auto logger = make_trace_logger();
    RayEvent rev;

    bool good = determine_interaction_type(logger, 0, 0, myrng, &optics,
                                           surface_normal(),
                                           ray_direction_for_angle(angle_mrad * 1.0e-3),
                                           /*LastHitBackSide=*/false, rev);

    EXPECT_TRUE(good);
    EXPECT_EQ(rev, RayEvent::ABSORB);
}

TEST(DetermineInteractionType, OneTransmissivityTableValueYieldsTransmit)
{
    OpticalPropertySet optics(InteractionType::REFRACTION, "one_transmissivity_table");
    optics.set_errors(OpticalSide::Both, DistributionType::NONE, 0.0, 0.0);

    constexpr double angle_mrad = 200.0;
    std::vector<AngularTablePoint> table = {{0.0, 0.0}, {angle_mrad, 1.0}, {1000.0, 0.0}};
    optics.set_transmissivity_table(OpticalSide::Front, table);

    MTRand myrng(1);
    auto logger = make_trace_logger();
    RayEvent rev;

    bool good = determine_interaction_type(logger, 0, 0, myrng, &optics,
                                           surface_normal(),
                                           ray_direction_for_angle(angle_mrad * 1.0e-3),
                                           /*LastHitBackSide=*/false, rev);

    EXPECT_TRUE(good);
    EXPECT_EQ(rev, RayEvent::TRANSMIT);
}

TEST(DetermineInteractionType, FrontAndBackSidesUseIndependentTables)
{
    OpticalPropertySet optics(InteractionType::REFLECTION, "front_back_tables");
    optics.set_errors(OpticalSide::Both, DistributionType::NONE, 0.0, 0.0);

    constexpr double angle_mrad = 400.0;
    std::vector<AngularTablePoint> front_table = {{0.0, 1.0}, {angle_mrad, 1.0}, {1000.0, 1.0}};
    std::vector<AngularTablePoint> back_table = {{0.0, 0.0}, {angle_mrad, 0.0}, {1000.0, 0.0}};
    optics.set_reflectivity_table(OpticalSide::Front, front_table);
    optics.set_reflectivity_table(OpticalSide::Back, back_table);

    auto logger = make_trace_logger();
    const glm::dvec3 dir = ray_direction_for_angle(angle_mrad * 1.0e-3);

    MTRand front_rng(1);
    RayEvent front_rev;
    bool front_good = determine_interaction_type(logger, 0, 0, front_rng, &optics,
                                                 surface_normal(), dir,
                                                 /*LastHitBackSide=*/false, front_rev);

    MTRand back_rng(1);
    RayEvent back_rev;
    bool back_good = determine_interaction_type(logger, 0, 0, back_rng, &optics,
                                                surface_normal(), dir,
                                                /*LastHitBackSide=*/true, back_rev);

    EXPECT_TRUE(front_good);
    EXPECT_EQ(front_rev, RayEvent::REFLECT);
    EXPECT_TRUE(back_good);
    EXPECT_EQ(back_rev, RayEvent::ABSORB);
}

TEST(DetermineInteractionType, DisablingTableFallsBackToScalarValue)
{
    OpticalPropertySet optics(InteractionType::REFLECTION, "disabled_table_fallback");
    optics.set_errors(OpticalSide::Both, DistributionType::NONE, 0.0, 0.0);
    optics.set_reflectivity(OpticalSide::Front, 1.0);

    constexpr double angle_mrad = 500.0;
    std::vector<AngularTablePoint> table = {{0.0, 0.0}, {angle_mrad, 0.0}, {1000.0, 0.0}};
    optics.set_reflectivity_table(OpticalSide::Front, table);
    optics.disable_reflectivity_table(OpticalSide::Front);

    MTRand myrng(1);
    auto logger = make_trace_logger();
    RayEvent rev;

    bool good = determine_interaction_type(logger, 0, 0, myrng, &optics,
                                           surface_normal(),
                                           ray_direction_for_angle(angle_mrad * 1.0e-3),
                                           /*LastHitBackSide=*/false, rev);

    EXPECT_TRUE(good);
    EXPECT_EQ(rev, RayEvent::REFLECT);
}

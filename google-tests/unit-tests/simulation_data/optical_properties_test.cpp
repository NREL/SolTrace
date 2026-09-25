#include <gtest/gtest.h>

#include <cmath>
#include <limits>
#include <vector>

#include <nlohmann/json.hpp>

#include <simulation_data_export.hpp>


#include "common.hpp"

TEST(OpticalProperties, DangleOpticalPointer)
{
    SimulationData sd;

    // Make optical property set
    auto opt_set = OpticalPropertySet();
    opt_set.set_ideal_reflection(OpticalSide::Both);

    // Add to simulation data
    auto opt_ref = sd.add_optical_property_set(opt_set);

    // Make element
    auto test_element = make_element<SingleElement>();
    test_element->set_aperture(make_aperture<Circle>(1));
    test_element->set_surface(make_surface<Flat>());

    // Attach optical property reference to element
    test_element->set_optical_property_set(opt_ref);

    // Add element to simulation data
    ASSERT_NO_THROW(sd.add_element(test_element));

    sd.clear();

    EXPECT_EQ(sd.get_optical_property_set(*test_element), nullptr);
}

TEST(OpticalProperties, MissingOptical)
{
    SimulationData sd;

    // Make element
    auto test_element = make_element<SingleElement>();
    test_element->set_aperture(make_aperture<Circle>(1));
    test_element->set_surface(make_surface<Flat>());

    // Add element to simulation data
    EXPECT_THROW(sd.add_element(test_element), std::invalid_argument);
}

TEST(OpticalProperties, InvalidErrorSettingsThrow)
{
    OpticalPropertySet optics(InteractionType::REFLECTION, "invalid_errors");
    optics.set_ideal_reflection(OpticalSide::Both);

    EXPECT_THROW(optics.set_errors(
                     OpticalSide::Both, DistributionType::UNKNOWN, 1.0, 1.0),
                 std::invalid_argument);
    EXPECT_THROW(optics.set_errors(
                     OpticalSide::Both, DistributionType::GAUSSIAN, -1.0, 1.0),
                 std::invalid_argument);
    EXPECT_THROW(optics.set_errors(OpticalSide::Both,
                                   DistributionType::GAUSSIAN,
                                   1.0,
                                   std::numeric_limits<double>::quiet_NaN()),
                 std::invalid_argument);
}

TEST(OpticalProperties, MutateOptics)
{
    SimulationData sd;

    // Make optical property set
    auto opt_set = OpticalPropertySet();
    opt_set.set_ideal_reflection(OpticalSide::Both);

    // Add to simulation data
    auto opt_ref = sd.add_optical_property_set(opt_set);

    // Make element 1
    auto test_element_1 = make_element<SingleElement>();
    test_element_1->set_aperture(make_aperture<Circle>(1));
    test_element_1->set_surface(make_surface<Flat>());

    // Attach optical property reference to element 1
    test_element_1->set_optical_property_set(opt_ref);

    // Make element 2
    auto test_element_2 = make_element<SingleElement>();
    test_element_2->set_aperture(make_aperture<Circle>(1));
    test_element_2->set_surface(make_surface<Flat>());

    // Attach same optical properties to elemetn 2
    test_element_2->set_optical_property_set(opt_ref);

    // Modify optical properties directly from simulation data
    auto* optics_ptr = sd.get_mutable_optical_property_set(*test_element_1);
    optics_ptr->set_reflectivity(OpticalSide::Front, 0.5);

    // Confirm reflectivity is modified for both
    EXPECT_EQ(test_element_1->get_optical_property_set()->get_reflectivity(
                  OpticalSide::Front),
              0.5);
    EXPECT_EQ(test_element_2->get_optical_property_set()->get_reflectivity(
                  OpticalSide::Front),
              0.5);

    EXPECT_EQ(test_element_1->get_optical_property_set(),
              test_element_2->get_optical_property_set());
}

TEST(OpticalProperties, NullPointer)
{
    // Make optics
    OpticalPropertySet opt_set = OpticalPropertySet();
    opt_set.set_ideal_reflection(OpticalSide::Both);
    auto opt_ptr = std::make_shared<OpticalPropertySet>(opt_set);

    // Make optics reference
    OpticalPropertySetReference ref = { 0, opt_ptr };

    // Make element
    auto test_element = make_element<SingleElement>();
    test_element->set_aperture(make_aperture<Circle>(1));
    test_element->set_surface(make_surface<Flat>());
    test_element->set_optical_property_set(ref);

    // Add element to simulation data
    SimulationData sd;
    sd.add_element(test_element);

    // Delete optics pointer
    opt_ptr.reset();

    // Try to get pointer
    EXPECT_EQ(test_element->get_optical_property_set(), nullptr);
}

TEST(OpticalProperties, ElementCanAccessSimulationOwnedOptics)
{
    SimulationData sd;

    OpticalPropertySet opt_set;
    opt_set.set_ideal_reflection(OpticalSide::Both);

    auto opt_ref = sd.add_optical_property_set(opt_set);

    auto test_element = make_element<SingleElement>();
    test_element->set_aperture(make_aperture<Circle>(1));
    test_element->set_surface(make_surface<Flat>());
    test_element->set_optical_property_set(opt_ref);

    EXPECT_NE(test_element->get_optical_property_set(), nullptr);
    EXPECT_EQ(test_element->get_optical_property_set()->get_reflectivity(
                  OpticalSide::Front),
              1.0);
}

TEST(OpticalProperties, RemovingElementDoesNotRemoveOptics)
{
    SimulationData sd;

    OpticalPropertySet opt_set;
    opt_set.set_ideal_reflection(OpticalSide::Both);

    auto opt_ref = sd.add_optical_property_set(opt_set);

    auto test_element = make_element<SingleElement>();
    test_element->set_aperture(make_aperture<Circle>(1));
    test_element->set_surface(make_surface<Flat>());
    test_element->set_optical_property_set(opt_ref);

    auto id = sd.add_element(test_element);
    ASSERT_TRUE(SolTrace::Data::Element::is_success(id));

    sd.remove_element(id);

    EXPECT_NE(sd.get_optical_property_set(*test_element), nullptr);
}

TEST(OpticalProperties, FindOrAddDeduplicatesEquivalentOptics)
{
    SimulationData sd;

    OpticalPropertySet opt_set_1;
    opt_set_1.set_ideal_reflection(OpticalSide::Both);

    OpticalPropertySet opt_set_2;
    opt_set_2.set_ideal_reflection(OpticalSide::Both);

    auto ref_1 = sd.find_or_add_optical_property_set(opt_set_1);
    auto ref_2 = sd.find_or_add_optical_property_set(opt_set_2);

    EXPECT_EQ(ref_1.id, ref_2.id);
    EXPECT_EQ(ref_1.optical_property_set.lock(),
              ref_2.optical_property_set.lock());
}

TEST(OpticalProperties, FindOrAddKeepsDifferentOpticsSeparate)
{
    SimulationData sd;

    OpticalPropertySet opt_set_1;
    opt_set_1.set_ideal_reflection(OpticalSide::Both);

    OpticalPropertySet opt_set_2;
    opt_set_2.set_ideal_reflection(OpticalSide::Both);
    opt_set_2.set_reflectivity(OpticalSide::Front, 0.5);

    auto ref_1 = sd.find_or_add_optical_property_set(opt_set_1);
    auto ref_2 = sd.find_or_add_optical_property_set(opt_set_2);

    EXPECT_NE(ref_1.id, ref_2.id);
    EXPECT_NE(ref_1.optical_property_set.lock(),
              ref_2.optical_property_set.lock());
}

TEST(OpticalProperties, AddElementWithExpiredOpticsThrows)
{
    auto test_element = make_element<SingleElement>();
    test_element->set_aperture(make_aperture<Circle>(1));
    test_element->set_surface(make_surface<Flat>());

    {
        SimulationData temporary_sd;

        OpticalPropertySet opt_set;
        opt_set.set_ideal_reflection(OpticalSide::Both);

        auto opt_ref = temporary_sd.add_optical_property_set(opt_set);
        test_element->set_optical_property_set(opt_ref);
    }

    SimulationData sd;

    EXPECT_THROW(sd.add_element(test_element), std::invalid_argument);
}

TEST(OpticalProperties, ReplaceElementDoesNotInvalidateOptics)
{
    SimulationData sd;

    OpticalPropertySet opt_set;
    opt_set.set_ideal_reflection(OpticalSide::Both);

    auto opt_ref = sd.add_optical_property_set(opt_set);

    auto old_element = make_element<SingleElement>();
    old_element->set_aperture(make_aperture<Circle>(1));
    old_element->set_surface(make_surface<Flat>());
    old_element->set_optical_property_set(opt_ref);

    auto id = sd.add_element(old_element);
    ASSERT_TRUE(SolTrace::Data::Element::is_success(id));

    auto new_element = make_element<SingleElement>();
    new_element->set_aperture(make_aperture<Circle>(2));
    new_element->set_surface(make_surface<Flat>());
    new_element->set_optical_property_set(opt_ref);

    EXPECT_TRUE(sd.replace_element(id, new_element));

    EXPECT_NE(new_element->get_optical_property_set(), nullptr);
    EXPECT_NE(sd.get_optical_property_set(*new_element), nullptr);
}

TEST(OpticalProperties, AngularTableInterpolatesExactAtTablePoints)
{
    OpticalPropertySet optics(InteractionType::REFLECTION,
                              "table_exact_points");
    optics.set_errors(OpticalSide::Both, DistributionType::NONE, 0.0, 0.0);

    std::vector<AngularTablePoint> table = { { 0.0, 0.9 },
                                             { 500.0, 0.5 },
                                             { 1000.0, 0.1 } };
    optics.set_reflectivity_table(OpticalSide::Front, table);

    EXPECT_DOUBLE_EQ(
        optics.get_reflectivity(OpticalSide::Front, std::cos(0.0 * 1.0e-3)),
        0.9);
    EXPECT_DOUBLE_EQ(
        optics.get_reflectivity(OpticalSide::Front, std::cos(500.0 * 1.0e-3)),
        0.5);
    EXPECT_DOUBLE_EQ(
        optics.get_reflectivity(OpticalSide::Front, std::cos(1000.0 * 1.0e-3)),
        0.1);
}

TEST(OpticalProperties, AngularTableInterpolation)
{
    OpticalPropertySet optics(InteractionType::REFLECTION,
                              "table_cosine_interp");
    optics.set_errors(OpticalSide::Both, DistributionType::NONE, 0.0, 0.0);

    constexpr double               angle0      = 0.0;
    constexpr double               angle1      = 1500.0;
    constexpr double               query_angle = 750.0;
    std::vector<AngularTablePoint> table = { { angle0, 1.0 }, { angle1, 0.0 } };
    optics.set_reflectivity_table(OpticalSide::Front, table);

    const double c0 = std::cos(angle0 * 1.0e-3);
    const double c1 = std::cos(angle1 * 1.0e-3);
    const double cq = std::cos(query_angle * 1.0e-3);
    const double expected_cosine_interp =
        1.0 + (cq - c0) / (c1 - c0) * (0.0 - 1.0);
    const double naive_angle_interp =
        0.5; // midpoint if interpolating linearly in angle instead

    const double actual = optics.get_reflectivity(OpticalSide::Front, cq);

    EXPECT_NEAR(actual, expected_cosine_interp, 1e-12);
    EXPECT_GT(std::fabs(actual - naive_angle_interp), 0.1);
}

TEST(OpticalProperties,
     AngularTableInterpolationSatisfiesSecondDerivativeErrorBound)
{
    // g(u) stands in for whatever function is tabulated against u = cos(theta),
    // since that's the variable the table is actually interpolated linearly in.
    // The standard normal pdf is used only because it is smooth with a simple,
    // exact second derivative -- this isn't meant to model Fresnel physics,
    // just to check the cosine-space interpolation formula against the classic
    // linear-interpolation error bound |g(u) - p1(u)| <= (h^2 / 8) * max|g''|.
    auto g = [](double u)
    { return std::exp(-0.5 * u * u) / std::sqrt(2.0 * PI); };

    // g''(u) = (u^2 - 1) * g(u); its only critical points are u = 0 and
    // u = +-sqrt(3), and |g''(0)| = g(0) is the largest of the three, so it
    // bounds |g''(u)| for every u -- no per-interval sampling needed.
    const double max_abs_g2 = g(0.0);

    constexpr double          kMradToRad  = 1.0e-3;
    const std::vector<double> angles_mrad = { 0.0,   200.0,  400.0,  600.0,
                                              800.0, 1000.0, 1200.0, 1400.0 };

    std::vector<AngularTablePoint> table;
    for (double angle : angles_mrad)
    {
        table.push_back(
            AngularTablePoint { angle, g(std::cos(angle * kMradToRad)) });
    }

    OpticalPropertySet optics(InteractionType::REFLECTION, "table_error_bound");
    optics.set_errors(OpticalSide::Both, DistributionType::NONE, 0.0, 0.0);
    optics.set_reflectivity_table(OpticalSide::Front, table);

    constexpr int kQueriesPerInterval = 25;

    for (std::size_t k = 1; k < angles_mrad.size(); ++k)
    {
        const double theta0 = angles_mrad[k - 1] * kMradToRad;
        const double theta1 = angles_mrad[k] * kMradToRad;
        const double h = std::cos(theta0) -
                         std::cos(theta1); // u decreases as theta increases

        const double error_bound = (h * h / 8.0) * max_abs_g2;

        for (int i = 1; i < kQueriesPerInterval; ++i)
        {
            const double t = static_cast<double>(i) / kQueriesPerInterval;
            const double query_theta      = theta0 + t * (theta1 - theta0);
            const double query_angle_mrad = query_theta / kMradToRad;

            const double actual = optics.get_reflectivity(
                OpticalSide::Front, std::cos(query_theta));
            const double exact = g(std::cos(query_theta));

            EXPECT_LE(std::fabs(actual - exact), error_bound)
                << "interval [" << angles_mrad[k - 1] << ", " << angles_mrad[k]
                << "] mrad, query angle " << query_angle_mrad << " mrad";
        }
    }
}

TEST(OpticalProperties, AngularTableClampsOutsideRange)
{
    OpticalPropertySet optics(InteractionType::REFRACTION, "table_clamp");
    optics.set_errors(OpticalSide::Both, DistributionType::NONE, 0.0, 0.0);


    std::vector<AngularTablePoint> table = { { 100.0, 0.8 }, { 900.0, 0.2 } };
    optics.set_transmissivity_table(OpticalSide::Front, table);

    EXPECT_DOUBLE_EQ(
        optics.get_transmissivity(OpticalSide::Front, std::cos(0.0 * 1.0e-3)),
        0.8);
    EXPECT_DOUBLE_EQ(optics.get_transmissivity(OpticalSide::Front,
                                               std::cos(2000.0 * 1.0e-3)),
                     0.2);
}

TEST(OpticalProperties, SinglePointAngularTableReturnsConstantValue)
{
    OpticalPropertySet optics(InteractionType::REFLECTION,
                              "table_single_point");
    optics.set_errors(OpticalSide::Both, DistributionType::NONE, 0.0, 0.0);

    std::vector<AngularTablePoint> table = { { 300.0, 0.42 } };
    optics.set_reflectivity_table(OpticalSide::Front, table);

    EXPECT_DOUBLE_EQ(
        optics.get_reflectivity(OpticalSide::Front, std::cos(0.0 * 1.0e-3)),
        0.42);
    EXPECT_DOUBLE_EQ(
        optics.get_reflectivity(OpticalSide::Front, std::cos(300.0 * 1.0e-3)),
        0.42);
    EXPECT_DOUBLE_EQ(
        optics.get_reflectivity(OpticalSide::Front, std::cos(5000.0 * 1.0e-3)),
        0.42);
}

TEST(OpticalProperties, EnablingEmptyAngularTableThrows)
{
    OpticalPropertySet optics(InteractionType::REFLECTION, "table_empty");
    optics.set_errors(OpticalSide::Both, DistributionType::NONE, 0.0, 0.0);

    EXPECT_THROW(optics.enable_reflectivity_table(OpticalSide::Front),
                 std::invalid_argument);
}

TEST(OpticalProperties, AngularTableRejectsNonAscendingAngles)
{
    OpticalPropertySet optics(InteractionType::REFLECTION,
                              "table_non_ascending");
    optics.set_errors(OpticalSide::Both, DistributionType::NONE, 0.0, 0.0);

    std::vector<AngularTablePoint> descending_table = { { 500.0, 0.5 },
                                                        { 100.0, 0.9 } };
    EXPECT_THROW(
        optics.set_reflectivity_table(OpticalSide::Front, descending_table),
        std::invalid_argument);

    std::vector<AngularTablePoint> duplicate_angle_table = { { 100.0, 0.9 },
                                                             { 100.0, 0.5 } };
    EXPECT_THROW(optics.set_reflectivity_table(OpticalSide::Front,
                                               duplicate_angle_table),
                 std::invalid_argument);
}

TEST(OpticalProperties, AngularTableRejectsInvalidValues)
{
    OpticalPropertySet optics(InteractionType::REFLECTION,
                              "table_invalid_values");
    optics.set_errors(OpticalSide::Both, DistributionType::NONE, 0.0, 0.0);

    std::vector<AngularTablePoint> negative_angle_table = { { -1.0, 0.9 } };
    EXPECT_THROW(
        optics.set_reflectivity_table(OpticalSide::Front, negative_angle_table),
        std::invalid_argument);

    std::vector<AngularTablePoint> nan_value_table = {
        { 100.0, std::numeric_limits<double>::quiet_NaN() }
    };
    EXPECT_THROW(
        optics.set_reflectivity_table(OpticalSide::Front, nan_value_table),
        std::invalid_argument);
}

TEST(OpticalProperties, DisablingAngularTableFallsBackToScalarValue)
{
    OpticalPropertySet optics(InteractionType::REFLECTION, "table_disable");
    optics.set_errors(OpticalSide::Both, DistributionType::NONE, 0.0, 0.0);
    optics.set_reflectivity(OpticalSide::Front, 0.75);

    std::vector<AngularTablePoint> table = { { 0.0, 0.9 }, { 1000.0, 0.1 } };
    optics.set_reflectivity_table(OpticalSide::Front, table);
    EXPECT_DOUBLE_EQ(
        optics.get_reflectivity(OpticalSide::Front, std::cos(0.0 * 1.0e-3)),
        0.9);

    optics.disable_reflectivity_table(OpticalSide::Front);
    EXPECT_DOUBLE_EQ(
        optics.get_reflectivity(OpticalSide::Front, std::cos(0.0 * 1.0e-3)),
        0.75);
    EXPECT_DOUBLE_EQ(optics.get_reflectivity(OpticalSide::Front), 0.75);

    // Re-enabling re-uses the previously-set table without needing to supply it
    // again.
    optics.enable_reflectivity_table(OpticalSide::Front);
    EXPECT_DOUBLE_EQ(
        optics.get_reflectivity(OpticalSide::Front, std::cos(0.0 * 1.0e-3)),
        0.9);
}

TEST(OpticalProperties, UsesAngularTableReflectsPerSideAndPerTableFlags)
{
    OpticalPropertySet optics(InteractionType::REFRACTION, "table_flags");
    optics.set_errors(OpticalSide::Both, DistributionType::NONE, 0.0, 0.0);

    EXPECT_FALSE(optics.uses_angular_table(OpticalSide::Front));
    EXPECT_FALSE(optics.uses_angular_table(OpticalSide::Back));

    std::vector<AngularTablePoint> table = { { 0.0, 0.9 }, { 1000.0, 0.1 } };
    optics.set_transmissivity_table(OpticalSide::Front, table);

    EXPECT_TRUE(optics.uses_angular_table(OpticalSide::Front));
    EXPECT_FALSE(optics.uses_angular_table(OpticalSide::Back));

    optics.disable_transmissivity_table(OpticalSide::Front);
    EXPECT_FALSE(optics.uses_angular_table(OpticalSide::Front));
}

TEST(OpticalProperties, AngularTableRoundTripsThroughJson)
{
    OpticalPropertySet optics(InteractionType::REFRACTION,
                              "table_json_roundtrip");
    optics.set_ideal_transmission(1.0, 1.5);

    std::vector<AngularTablePoint> table = { { 0.0, 1.0 },
                                             { 500.0, 0.8 },
                                             { 1000.0, 0.2 } };
    optics.set_transmissivity_table(OpticalSide::Front, table);

    nlohmann::ordered_json jnode;
    optics.write_json(jnode);

    OpticalPropertySet roundtripped(jnode);

    EXPECT_EQ(optics, roundtripped);
    EXPECT_DOUBLE_EQ(roundtripped.get_transmissivity(OpticalSide::Front,
                                                     std::cos(500.0 * 1.0e-3)),
                     0.8);
}

TEST(OpticalProperties, MissingTableKeysInJsonDefaultToDisabled)
{
    nlohmann::ordered_json jnode;
    jnode["my_type"] =
        SolTrace::Data::InteractionTypeMap.at(InteractionType::REFLECTION);
    jnode["my_name"]                = "no_table_keys";
    jnode["refraction_index_front"] = 1.0;
    jnode["refraction_index_back"]  = 1.0;

    nlohmann::ordered_json jside;
    jside["error_distribution_type"] =
        SolTrace::Data::DistributionTypeMap.at(DistributionType::NONE);
    jside["transmissivity"]    = 0.0;
    jside["reflectivity"]      = 0.6;
    jside["slope_error"]       = 0.0;
    jside["specularity_error"] = 0.0;

    jnode["front"] = jside;
    jnode["back"]  = jside;

    OpticalPropertySet optics(jnode);

    EXPECT_FALSE(optics.uses_angular_table(OpticalSide::Front));
    EXPECT_FALSE(optics.uses_angular_table(OpticalSide::Back));
    EXPECT_DOUBLE_EQ(
        optics.get_reflectivity(OpticalSide::Front, std::cos(123.0 * 1.0e-3)),
        0.6);
}

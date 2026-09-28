#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include <constants.hpp>
#include <optical_properties.hpp>
#include <optix_runner.hpp>
#include <simulation_data.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result_export.hpp>

using SolTrace::Runner::RunnerStatus;

using SolTrace::Data::AngularTablePoint;
using SolTrace::Data::D2R;
using SolTrace::Data::DistributionType;
using SolTrace::Data::InteractionType;
using SolTrace::Data::OpticalPropertySet;
using SolTrace::Data::OpticalSide;

namespace
{

// Sun directly overhead (rays travel in -Z). Plate is tilted about the X axis
// so the ray's angle of incidence on the plate equals incidence_angle_deg.
void make_tilted_plate_sd(SimulationData& sd,
                          element_ptr&    plate,
                          double          incidence_angle_deg)
{
    sd.clear();

    auto sun = make_ray_source<Sun>();
    sun->set_position(0, 0, 100);
    sd.add_ray_source(sun);

    auto stage = make_stage(0);
    stage->set_origin(0, 0, 0);
    stage->set_aim_vector(0, 0, 1);
    stage->set_name("stage");

    plate = make_element<SingleElement>();
    plate->set_origin(0, 0, 50);

    const double theta = incidence_angle_deg * D2R;
    plate->set_aim_vector(
        0, 100.0 * std::sin(theta), 50.0 + 100.0 * std::cos(theta));
    plate->set_surface(make_surface<Flat>());
    plate->set_aperture(make_aperture<Rectangle>(5, 5));
    plate->set_name("plate");

    // Placeholder optics; tests mutate this in place via
    // get_mutable_optical_property_set() once the element is registered.
    OpticalPropertySet default_optics(
        InteractionType::REFRACTION, 1.0, 1.0, "PlateOptics");
    default_optics.set_properties(
        OpticalSide::Both, DistributionType::NONE, 1.0, 0.0, 0.0, 0.0);
    auto default_optics_ref = sd.add_optical_property_set(default_optics);
    plate->set_optical_property_set(default_optics_ref);

    stage->add_element(plate);
    sd.add_stage(stage);

    SimulationParameters& params    = sd.get_simulation_parameters();
    params.number_of_rays           = 20000;
    params.max_number_of_rays       = params.number_of_rays * 100;
    params.include_optical_errors   = false;
    params.include_sun_shape_errors = false;
    params.seed                     = 123;
}

void count_hits(const SimulationResult& result,
                int&                    absorbed_count,
                int&                    transmitted_count,
                int&                    reflected_count)
{
    absorbed_count = transmitted_count = reflected_count = 0;
    int n_records = result.get_number_of_records();

    for (int i = 0; i < n_records; i++)
    {
        ray_record_ptr rec            = result[i];
        int            n_interactions = rec->get_number_of_interactions();
        for (int j = 0; j < n_interactions; j++)
        {
            RayEvent rev = rec->get_event(j);

            if (rev == RayEvent::ABSORB) absorbed_count++;
            else if (rev == RayEvent::TRANSMIT)
                transmitted_count++;
            else if (rev == RayEvent::REFLECT)
                reflected_count++;
        }
    }
}

// Independent reference implementation of the unpolarized Fresnel reflectance
// (textbook form), used as ground truth for the OptiX shader's Monte Carlo
// split.
double
expected_fresnel_reflectance(double n1, double n2, double incidence_angle_deg)
{
    const double theta_i = incidence_angle_deg * D2R;
    const double sin_t   = (n1 / n2) * std::sin(theta_i);
    if (sin_t >= 1.0) return 1.0; // Total internal reflection

    const double cos_i = std::cos(theta_i);
    const double cos_t = std::sqrt(1.0 - sin_t * sin_t);

    const double rs =
        std::pow((n1 * cos_i - n2 * cos_t) / (n1 * cos_i + n2 * cos_t), 2);
    const double rp =
        std::pow((n1 * cos_t - n2 * cos_i) / (n1 * cos_t + n2 * cos_i), 2);
    return 0.5 * (rs + rp);
}

double deg_to_mrad(double deg) { return deg * D2R * 1000.0; }

// Independent reference implementation of the angular table lookup (linear
// interpolation in cos(angle), clamped at the ends), used as ground truth for
// the OptiX shader's lookup_angular_table().
double interpolate_angular_table(const std::vector<AngularTablePoint>& table,
                                 double incidence_angle_deg)
{
    const double ci = std::cos(incidence_angle_deg * D2R);

    std::vector<double> cos_cache(table.size());
    for (std::size_t i = 0; i < table.size(); ++i)
        cos_cache[i] = std::cos(table[i].angle * 1.0e-3);

    if (ci >= cos_cache.front()) return table.front().value;
    if (ci <= cos_cache.back()) return table.back().value;

    std::size_t k = 1;
    while (cos_cache[k] > ci)
        ++k;

    const double c0 = cos_cache[k - 1];
    const double c1 = cos_cache[k];
    const double v0 = table[k - 1].value;
    const double v1 = table[k].value;
    return v0 + (ci - c0) / (c1 - c0) * (v1 - v0);
}

} // namespace

TEST(MaterialInteraction, FresnelSplitAtFixedAngle)
{
    const double incidence_angle_deg = 30.0;
    const double n_incident          = 1.0;
    const double n_transmit          = 1.5;

    SimulationData sd;
    element_ptr    plate;
    make_tilted_plate_sd(sd, plate, incidence_angle_deg);

    OpticalPropertySet* plate_optics = sd.get_mutable_optical_property_set(*plate);
    ASSERT_NE(plate_optics, nullptr);
    plate_optics->set_refraction_indices(n_incident, n_transmit);
    plate_optics->set_transmissivity(OpticalSide::Both, 1.0);

    OptixRunner runner;
    ASSERT_EQ(runner.initialize(), RunnerStatus::SUCCESS);
    ASSERT_EQ(runner.setup_simulation(&sd), RunnerStatus::SUCCESS);
    ASSERT_EQ(runner.run_simulation(), RunnerStatus::SUCCESS);

    SimulationResult result;
    ASSERT_EQ(runner.report_simulation(&result, 0), RunnerStatus::SUCCESS);

    int absorbed, transmitted, reflected;
    count_hits(result, absorbed, transmitted, reflected);
    EXPECT_EQ(absorbed, 0); // transmissivity = 1, so nothing should be absorbed

    const int total = transmitted + reflected;
    ASSERT_GT(total, 0);
    const double measured_reflectance = static_cast<double>(reflected) / total;
    const double expected             = expected_fresnel_reflectance(
        n_incident, n_transmit, incidence_angle_deg);

    EXPECT_NEAR(measured_reflectance, expected, 0.02);
}

TEST(MaterialInteraction, ReflectanceIncreasesTowardGrazingAngle)
{
    const double              n_incident = 1.0;
    const double              n_transmit = 1.5;
    const std::vector<double> angles_deg = { 0.0, 40.0, 60.0, 80.0 };

    double previous_reflectance = -1.0;
    double first_reflectance    = -1.0;
    for (double angle : angles_deg)
    {
        SimulationData sd;
        element_ptr    plate;
        make_tilted_plate_sd(sd, plate, angle);

        OpticalPropertySet* plate_optics =
            sd.get_mutable_optical_property_set(*plate);
        ASSERT_NE(plate_optics, nullptr);
        plate_optics->set_refraction_indices(n_incident, n_transmit);
        plate_optics->set_transmissivity(OpticalSide::Both, 1.0);

        OptixRunner runner;
        ASSERT_EQ(runner.initialize(), RunnerStatus::SUCCESS);
        ASSERT_EQ(runner.setup_simulation(&sd), RunnerStatus::SUCCESS);
        ASSERT_EQ(runner.run_simulation(), RunnerStatus::SUCCESS);

        SimulationResult result;
        ASSERT_EQ(runner.report_simulation(&result, 0), RunnerStatus::SUCCESS);

        int absorbed, transmitted, reflected;
        count_hits(result, absorbed, transmitted, reflected);
        const int total = transmitted + reflected;
        ASSERT_GT(total, 0);
        const double reflectance = static_cast<double>(reflected) / total;

        EXPECT_GE(reflectance, previous_reflectance - 0.02)
            << "angle=" << angle;
        if (first_reflectance < 0.0) first_reflectance = reflectance;
        previous_reflectance = reflectance;
    }

    // Reflectance at near-grazing incidence should be well above normal
    // incidence.
    EXPECT_GT(previous_reflectance, first_reflectance + 0.2);
}

TEST(MaterialInteraction, CombinedTransmissivityAndFresnel)
{
    const double incidence_angle_deg = 20.0;
    const double n_incident          = 1.0;
    const double n_transmit          = 1.5;
    const double transmissivity      = 0.7;

    SimulationData sd;
    element_ptr    plate;
    make_tilted_plate_sd(sd, plate, incidence_angle_deg);

    OpticalPropertySet* plate_optics = sd.get_mutable_optical_property_set(*plate);
    ASSERT_NE(plate_optics, nullptr);
    plate_optics->set_refraction_indices(n_incident, n_transmit);
    plate_optics->set_transmissivity(OpticalSide::Both, transmissivity);

    OptixRunner runner;
    ASSERT_EQ(runner.initialize(), RunnerStatus::SUCCESS);
    ASSERT_EQ(runner.setup_simulation(&sd), RunnerStatus::SUCCESS);
    ASSERT_EQ(runner.run_simulation(), RunnerStatus::SUCCESS);

    SimulationResult result;
    ASSERT_EQ(runner.report_simulation(&result, 0), RunnerStatus::SUCCESS);

    int absorbed, transmitted, reflected;
    count_hits(result, absorbed, transmitted, reflected);
    const int total = absorbed + transmitted + reflected;
    ASSERT_GT(total, 0);

    // Absorption is governed purely by transmissivity, independent of angle.
    const double absorbed_fraction = static_cast<double>(absorbed) / total;
    EXPECT_NEAR(absorbed_fraction, 1.0 - transmissivity, 0.02);

    // Among rays that survive absorption, the reflect/transmit split still
    // follows Fresnel.
    const int survived = transmitted + reflected;
    ASSERT_GT(survived, 0);
    const double measured_reflectance =
        static_cast<double>(reflected) / survived;
    const double expected = expected_fresnel_reflectance(
        n_incident, n_transmit, incidence_angle_deg);

    EXPECT_NEAR(measured_reflectance, expected, 0.03);
}

TEST(AngularTableInteraction, ReflectivityTableMatchesInterpolation)
{
    // Breakpoints at 0/30/60 deg; 75 deg exercises clamping beyond the
    // table's range.
    const std::vector<AngularTablePoint> table = {
        { deg_to_mrad(0.0), 1.0 },
        { deg_to_mrad(30.0), 0.7 },
        { deg_to_mrad(60.0), 0.3 }
    };
    const std::vector<double> angles_deg = { 0.0, 15.0, 30.0, 45.0, 60.0, 75.0 };

    for (double angle : angles_deg)
    {
        SimulationData sd;
        element_ptr    plate;
        make_tilted_plate_sd(sd, plate, angle);

        OpticalPropertySet* plate_optics =
            sd.get_mutable_optical_property_set(*plate);
        ASSERT_NE(plate_optics, nullptr);
        plate_optics->set_interaction_type(InteractionType::REFLECTION);
        plate_optics->set_reflectivity_table(OpticalSide::Both, table);

        OptixRunner runner;
        ASSERT_EQ(runner.initialize(), RunnerStatus::SUCCESS);
        ASSERT_EQ(runner.setup_simulation(&sd), RunnerStatus::SUCCESS);
        ASSERT_EQ(runner.run_simulation(), RunnerStatus::SUCCESS);

        SimulationResult result;
        ASSERT_EQ(runner.report_simulation(&result, 0), RunnerStatus::SUCCESS);

        int absorbed, transmitted, reflected;
        count_hits(result, absorbed, transmitted, reflected);
        const int total = absorbed + reflected;
        ASSERT_GT(total, 0);

        const double measured_reflectance =
            static_cast<double>(reflected) / total;
        const double expected = interpolate_angular_table(table, angle);

        EXPECT_NEAR(measured_reflectance, expected, 0.03) << "angle=" << angle;
    }
}

TEST(AngularTableInteraction, TransmissivityTableGovernsAbsorption)
{
    // Matched refraction indices (mu=1) neutralize refract()'s internal
    // Fresnel split, isolating the transmissivity table's absorb/transmit
    // split.
    const double n_incident = 1.0;
    const double n_transmit = 1.0;

    const std::vector<AngularTablePoint> table = {
        { deg_to_mrad(0.0), 1.0 },
        { deg_to_mrad(20.0), 0.75 },
        { deg_to_mrad(50.0), 0.2 }
    };
    const std::vector<double> angles_deg = { 0.0, 20.0, 35.0, 65.0 };

    for (double angle : angles_deg)
    {
        SimulationData sd;
        element_ptr    plate;
        make_tilted_plate_sd(sd, plate, angle);

        OpticalPropertySet* plate_optics =
            sd.get_mutable_optical_property_set(*plate);
        ASSERT_NE(plate_optics, nullptr);
        plate_optics->set_refraction_indices(n_incident, n_transmit);
        plate_optics->set_transmissivity_table(OpticalSide::Both, table);

        OptixRunner runner;
        ASSERT_EQ(runner.initialize(), RunnerStatus::SUCCESS);
        ASSERT_EQ(runner.setup_simulation(&sd), RunnerStatus::SUCCESS);
        ASSERT_EQ(runner.run_simulation(), RunnerStatus::SUCCESS);

        SimulationResult result;
        ASSERT_EQ(runner.report_simulation(&result, 0), RunnerStatus::SUCCESS);

        int absorbed, transmitted, reflected;
        count_hits(result, absorbed, transmitted, reflected);
        const int total = absorbed + transmitted + reflected;
        ASSERT_GT(total, 0);

        const double absorbed_fraction = static_cast<double>(absorbed) / total;
        const double expected_transmissivity =
            interpolate_angular_table(table, angle);

        EXPECT_NEAR(absorbed_fraction, 1.0 - expected_transmissivity, 0.03)
            << "angle=" << angle;
        EXPECT_NEAR(static_cast<double>(reflected) / total, 0.0, 0.01)
            << "angle=" << angle;
    }
}

namespace
{

// Ideal flat mirror tilted by tilt_angle_deg, with an ideal absorbing target
// placed exactly where the law of reflection (angle in = angle out) predicts
// the reflected beam should land.
void make_mirror_and_target_sd(SimulationData& sd,
                               element_ptr&    mirror,
                               element_ptr&    target,
                               double          tilt_angle_deg)
{
    sd.clear();

    auto sun = make_ray_source<Sun>();
    sun->set_position(0, 0, 100);
    sd.add_ray_source(sun);

    auto stage = make_stage(0);
    stage->set_origin(0, 0, 0);
    stage->set_aim_vector(0, 0, 1);
    stage->set_name("stage");

    const double theta     = tilt_angle_deg * D2R;
    const double mirror_x = 0.0, mirror_y = 0.0, mirror_z = 50.0;

    mirror = make_element<SingleElement>();
    mirror->set_origin(mirror_x, mirror_y, mirror_z);
    mirror->set_aim_vector(mirror_x,
                           mirror_y + 100.0 * std::sin(theta),
                           mirror_z + 100.0 * std::cos(theta));
    mirror->set_surface(make_surface<Flat>());
    mirror->set_aperture(make_aperture<Rectangle>(5, 5));
    mirror->set_name("mirror");

    OpticalPropertySet mirror_optics(
        InteractionType::REFLECTION, 0.0, 0.0, "IdealMirror");
    mirror_optics.set_ideal_reflection(OpticalSide::Both);
    auto mirror_optics_ref = sd.add_optical_property_set(mirror_optics);
    mirror->set_optical_property_set(mirror_optics_ref);

    // Incident ray (0,0,-1) reflects to (0, sin(2*theta), cos(2*theta)).
    const double rx       = 0.0;
    const double ry       = std::sin(2.0 * theta);
    const double rz       = std::cos(2.0 * theta);
    const double distance = 50.0;

    const double target_x = mirror_x + distance * rx;
    const double target_y = mirror_y + distance * ry;
    const double target_z = mirror_z + distance * rz;

    target = make_element<SingleElement>();
    target->set_origin(target_x, target_y, target_z);
    target->set_aim_vector(
        target_x - 100.0 * rx, target_y - 100.0 * ry, target_z - 100.0 * rz);
    target->set_surface(make_surface<Flat>());
    target->set_aperture(make_aperture<Rectangle>(10, 10));
    target->set_name("target");

    OpticalPropertySet target_optics(
        InteractionType::REFLECTION, 0.0, 0.0, "IdealAbsorber");
    target_optics.set_ideal_absorption(OpticalSide::Both);
    auto target_optics_ref = sd.add_optical_property_set(target_optics);
    target->set_optical_property_set(target_optics_ref);

    stage->add_element(mirror);
    stage->add_element(target);
    sd.add_stage(stage);

    SimulationParameters& params    = sd.get_simulation_parameters();
    params.number_of_rays           = 20000;
    params.max_number_of_rays       = params.number_of_rays * 100;
    params.include_optical_errors   = false;
    params.include_sun_shape_errors = false;
    params.seed                     = 321;
}

} // namespace

TEST(MaterialInteraction, ReflectionFollowsLawOfReflection)
{
    const std::vector<double> tilt_angles_deg = { 10.0, 25.0, 40.0 };

    for (double tilt_angle_deg : tilt_angles_deg)
    {
        SimulationData sd;
        element_ptr    mirror, target;
        make_mirror_and_target_sd(sd, mirror, target, tilt_angle_deg);

        OptixRunner runner;
        ASSERT_EQ(runner.initialize(), RunnerStatus::SUCCESS);
        ASSERT_EQ(runner.setup_simulation(&sd), RunnerStatus::SUCCESS);
        ASSERT_EQ(runner.run_simulation(), RunnerStatus::SUCCESS);

        SimulationResult result;
        ASSERT_EQ(runner.report_simulation(&result, 0), RunnerStatus::SUCCESS);

        int absorbed, transmitted, reflected;
        count_hits(result, absorbed, transmitted, reflected);
        ASSERT_GT(reflected, 0);

        // A wrong sign/formula in reflect() sends the beam elsewhere and this
        // capture fraction collapses toward zero.
        const double capture_fraction =
            static_cast<double>(absorbed) / static_cast<double>(reflected);
        EXPECT_GT(capture_fraction, 0.95) << "tilt_angle_deg=" << tilt_angle_deg;
    }
}

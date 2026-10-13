#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include <native_runner.hpp>
#include <optix_runner.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result_export.hpp>

#include "common.hpp"

using SolTrace::NativeRunner::NativeRunner;
using SolTrace::Runner::RunnerStatus;

namespace
{
    // Builds a scene with a sun directly overhead (collimated rays along -z
    // once sun-shape errors are disabled) and a single flat mirror tilted by
    // incidence_angle_rad from normal incidence, with a Gaussian slope error
    // applied to its surface normal. Returns the mirror's id.
    element_id build_tilted_mirror_scene_with_slope_error(
        SimulationData &sd,
        double incidence_angle_rad,
        double slope_error_mrad,
        const std::vector<AngularTablePoint> &reflectivity_table)
    {
        auto sun = SolTrace::Data::make_ray_source<Sun>();
        sun->set_position(0.0, 0.0, 100.0);
        sun->set_shape(SolTrace::Data::SunShape::PILLBOX, -1.0, 1.0, 0.0);
        sd.add_ray_source(sun);

        OpticalPropertySet mirror_optics(InteractionType::REFLECTION, "NativeOptixParityOptics");
        mirror_optics.set_errors(OpticalSide::Front, DistributionType::GAUSSIAN, slope_error_mrad, 0.0);
        mirror_optics.set_errors(OpticalSide::Back, DistributionType::NONE, 0.0, 0.0);
        mirror_optics.set_reflectivity_table(OpticalSide::Front, reflectivity_table);
        auto optics_ref = sd.add_optical_property_set(mirror_optics);

        auto mirror = SolTrace::Data::make_element<SingleElement>();
        mirror->set_aperture(SolTrace::Data::make_aperture<Circle>(5.0));
        mirror->set_surface(SolTrace::Data::make_surface<Flat>());
        mirror->set_reference_frame_geometry(
            glm::dvec3(0.0, 0.0, 0.0),
            glm::dvec3(std::sin(incidence_angle_rad), 0.0, std::cos(incidence_angle_rad)),
            0.0);
        mirror->set_optical_property_set(optics_ref);

        return sd.add_element(mirror);
    }

    // Nominal incidence angle sits inside the table's fully-reflective
    // plateau; only slope-error-perturbed rays can reach the falloff zone
    // starting at 350 mrad. A runner that evaluates the angular table
    // BEFORE applying slope error would see the fixed 300 mrad angle on
    // every ray and never absorb anything.
    constexpr double kIncidenceAngleMrad = 300.0;
    constexpr double kIncidenceAngleRad  = kIncidenceAngleMrad * 1.0e-3;
    constexpr double kSlopeErrorMrad     = 50.0;
    constexpr uint_fast64_t kNumRays     = 20000;

    const std::vector<AngularTablePoint> kFalloffTable = {
        {0.0, 1.0}, {350.0, 1.0}, {450.0, 0.0}, {1500.0, 0.0}};

    void configure_params(SimulationParameters &params)
    {
        params.include_optical_errors = true;
        params.include_sun_shape_errors = false;
        params.number_of_rays = kNumRays;
        params.max_number_of_rays = kNumRays * 10;
        params.seed = 1;
    }
} // namespace

// Regression guard for the fix that made all runners apply slope error to
// the surface normal BEFORE evaluating an angular reflectivity/transmissivity
// table, so the table sees the same perturbed incidence angle used for the
// interaction itself. Prior to that fix, native/Embree looked up the table
// against the pre-perturbation angle while OptiX used the post-perturbation
// angle, so identical scenes could diverge whenever both errors and an
// angular table were enabled together. OptiX already applied the correct
// (post-perturbation) ordering, so this guards against a future regression
// on either side rather than a currently-known bug.
TEST(NativeOptixParity, SlopeErrorAppliedBeforeAngularTableLookupInBothRunners)
{
    SimulationData sd_native;
    configure_params(sd_native.get_simulation_parameters());
    element_id native_mirror_id = build_tilted_mirror_scene_with_slope_error(
        sd_native, kIncidenceAngleRad, kSlopeErrorMrad, kFalloffTable);

    SimulationData sd_optix;
    configure_params(sd_optix.get_simulation_parameters());
    element_id optix_mirror_id = build_tilted_mirror_scene_with_slope_error(
        sd_optix, kIncidenceAngleRad, kSlopeErrorMrad, kFalloffTable);

    NativeRunner native_runner;
    ASSERT_EQ(native_runner.initialize(), RunnerStatus::SUCCESS);
    ASSERT_EQ(native_runner.setup_simulation(&sd_native), RunnerStatus::SUCCESS);
    ASSERT_EQ(native_runner.run_simulation(), RunnerStatus::SUCCESS);
    SimulationResult native_result;
    ASSERT_EQ(native_runner.report_simulation(&native_result, 0), RunnerStatus::SUCCESS);

    OptixRunner optix_runner;
    ASSERT_EQ(optix_runner.initialize(), RunnerStatus::SUCCESS);
    ASSERT_EQ(optix_runner.setup_simulation(&sd_optix), RunnerStatus::SUCCESS);
    ASSERT_EQ(optix_runner.run_simulation(), RunnerStatus::SUCCESS);
    SimulationResult optix_result;
    ASSERT_EQ(optix_runner.report_simulation(&optix_result, 0), RunnerStatus::SUCCESS);

    const int_fast64_t native_absorbed = count_element_event(native_result, native_mirror_id, RayEvent::ABSORB);
    const int_fast64_t native_reflected = count_element_event(native_result, native_mirror_id, RayEvent::REFLECT);
    const int_fast64_t optix_absorbed = count_element_event(optix_result, optix_mirror_id, RayEvent::ABSORB);
    const int_fast64_t optix_reflected = count_element_event(optix_result, optix_mirror_id, RayEvent::REFLECT);

    ASSERT_EQ(native_absorbed + native_reflected, static_cast<int_fast64_t>(kNumRays));
    ASSERT_EQ(optix_absorbed + optix_reflected, static_cast<int_fast64_t>(kNumRays));

    // A regression back to pre-perturbation table lookups would make this
    // exactly zero, since the unperturbed 300 mrad angle never reaches the
    // table's falloff zone (350-450 mrad).
    EXPECT_GT(native_absorbed, 0);
    EXPECT_GT(optix_absorbed, 0);

    const double native_fraction = static_cast<double>(native_reflected) / static_cast<double>(kNumRays);
    const double optix_fraction = static_cast<double>(optix_reflected) / static_cast<double>(kNumRays);

    // Both runners should land on the same reflected fraction (well below
    // 1.0) since both perturb the normal before the table lookup. A 5%
    // absolute tolerance comfortably covers Monte Carlo noise between the
    // two independent RNG streams while still catching an ordering
    // regression, which would show up as one runner reporting ~1.0.
    EXPECT_NEAR(native_fraction, optix_fraction, 0.05);
}

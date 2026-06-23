#include "database/database.h"

#include "composite_element.hpp"
#include "native_runner.hpp"
#include "simulation_data.hpp"
#include "simulation_result.hpp"
#include "simulation_runner.hpp"
#include "sun.hpp"

#include <gtest/gtest.h>

#include <glm/geometric.hpp>
#include <glm/mat3x3.hpp>
#include <glm/vec3.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <map>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

namespace {

namespace SD = SolTrace::Data;

constexpr double kTransformTolerance = 1.0e-7;
constexpr double kValueTolerance     = 1.0e-12;

struct ElementSnapshot {
    std::string name;
    bool        enabled = false;

    glm::dvec3 origin { 0.0 };
    glm::dvec3 x_axis { 0.0 };
    glm::dvec3 y_axis { 0.0 };
    glm::dvec3 z_axis { 0.0 };
    glm::dvec3 aim_axis { 0.0 };

    std::string aperture_json;
    std::string surface_json;
    std::string optics_json;
};

std::string dump_json(nlohmann::ordered_json const& node) {
    return node.dump();
}

std::string dump_aperture(SD::Element const& element) {
    auto aperture = element.get_aperture();
    if (!aperture) return "null";

    nlohmann::ordered_json node;
    aperture->write_json(node);
    return dump_json(node);
}

std::string dump_surface(SD::Element const& element) {
    auto surface = element.get_surface();
    if (!surface) return "null";

    nlohmann::ordered_json node;
    surface->write_json(node);
    return dump_json(node);
}

std::string dump_optics(std::shared_ptr<const SD::OpticalPropertySet> optics) {
    if (!optics) return "null";

    nlohmann::ordered_json node;
    optics->write_json(node);
    return dump_json(node);
}

glm::dvec3 normalized_or_zero(glm::dvec3 const& vector) {
    auto length = glm::length(vector);
    if (length == 0.0) return glm::dvec3 { 0.0 };
    return vector / length;
}

double vec_error(glm::dvec3 const& actual, glm::dvec3 const& expected) {
    return glm::length(actual - expected);
}

std::string vec_to_string(glm::dvec3 const& vector) {
    std::ostringstream out;
    out << std::setprecision(12) << "(" << vector.x << ", " << vector.y << ", "
        << vector.z << ")";
    return out.str();
}

std::optional<double> cylinder_radius(SD::Element const& element) {
    auto cylinder =
        std::dynamic_pointer_cast<SD::Cylinder const>(element.get_surface());
    if (!cylinder) return std::nullopt;
    if (!(cylinder->radius > 0.0) || !std::isfinite(cylinder->radius)) {
        return std::nullopt;
    }
    return cylinder->radius;
}

ElementSnapshot snapshot_element(SD::Element const& element,
                                 bool legacy_stinput_cylinder_origin = false) {
    auto local_to_global = element.get_local_to_global();
    auto origin          = element.get_origin_global();
    auto z_axis =
        normalized_or_zero(local_to_global * glm::dvec3 { 0.0, 0.0, 1.0 });

    if (legacy_stinput_cylinder_origin) {
        if (auto radius = cylinder_radius(element); radius) {
            origin += z_axis * *radius;
        }
    }

    return ElementSnapshot {
        .name    = element.get_name(),
        .enabled = element.is_enabled(),
        .origin  = origin,
        .x_axis =
            normalized_or_zero(local_to_global * glm::dvec3 { 1.0, 0.0, 0.0 }),
        .y_axis =
            normalized_or_zero(local_to_global * glm::dvec3 { 0.0, 1.0, 0.0 }),
        .z_axis = z_axis,
        .aim_axis =
            normalized_or_zero(element.get_aim_vector_global() - origin),
        .aperture_json = dump_aperture(element),
        .surface_json  = dump_surface(element),
        .optics_json   = dump_optics(element.get_optical_property_set()),
    };
}

std::int64_t quantize(double value) {
    return static_cast<std::int64_t>(std::llround(value * 1.0e6));
}

auto sort_key(ElementSnapshot const& item) {
    return std::tuple {
        quantize(item.origin.x), quantize(item.origin.y),
        quantize(item.origin.z), item.aperture_json,
        item.surface_json,       item.name,
    };
}

std::vector<ElementSnapshot>
collect_single_element_snapshots(SD::SimulationData const& data,
                                 bool legacy_stinput_cylinder_origin = false) {
    std::vector<ElementSnapshot> snapshots;
    snapshots.reserve(data.get_number_of_elements());

    for (auto iter = data.get_const_iterator(); !data.is_at_end(iter); ++iter) {
        auto const& element = *iter->second;
        if (!element.is_single()) continue;

        snapshots.push_back(
            snapshot_element(element, legacy_stinput_cylinder_origin));
    }

    std::sort(
        snapshots.begin(), snapshots.end(), [](auto const& a, auto const& b) {
            return sort_key(a) < sort_key(b);
        });

    return snapshots;
}

std::string mismatch_diagnostic(ElementSnapshot const& actual,
                                ElementSnapshot const& expected) {
    std::ostringstream out;
    out << std::setprecision(12) << "name actual/expected: " << actual.name
        << " / " << expected.name
        << "\norigin error: " << vec_error(actual.origin, expected.origin)
        << "\n  actual origin:   " << vec_to_string(actual.origin)
        << "\n  expected origin: " << vec_to_string(expected.origin)
        << "\naim error: " << vec_error(actual.aim_axis, expected.aim_axis)
        << "\n  actual aim:      " << vec_to_string(actual.aim_axis)
        << "\n  expected aim:    " << vec_to_string(expected.aim_axis)
        << "\nx axis error: " << vec_error(actual.x_axis, expected.x_axis)
        << "\n  actual x axis:   " << vec_to_string(actual.x_axis)
        << "\n  expected x axis: " << vec_to_string(expected.x_axis)
        << "\ny axis error: " << vec_error(actual.y_axis, expected.y_axis)
        << "\n  actual y axis:   " << vec_to_string(actual.y_axis)
        << "\n  expected y axis: " << vec_to_string(expected.y_axis)
        << "\nz axis error: " << vec_error(actual.z_axis, expected.z_axis)
        << "\n  actual z axis:   " << vec_to_string(actual.z_axis)
        << "\n  expected z axis: " << vec_to_string(expected.z_axis);
    return out.str();
}

std::string max_error_summary(std::vector<ElementSnapshot> const& actual,
                              std::vector<ElementSnapshot> const& expected) {
    auto count = std::min(actual.size(), expected.size());

    struct MaxError {
        double value = 0.0;
        size_t index = 0;
    };

    auto update = [](MaxError& max_error, double value, size_t index) {
        if (value > max_error.value) {
            max_error.value = value;
            max_error.index = index;
        }
    };

    MaxError origin;
    MaxError aim;
    MaxError x_axis;
    MaxError y_axis;
    MaxError z_axis;

    for (size_t i = 0; i < count; ++i) {
        update(origin, vec_error(actual[i].origin, expected[i].origin), i);
        update(aim, vec_error(actual[i].aim_axis, expected[i].aim_axis), i);
        update(x_axis, vec_error(actual[i].x_axis, expected[i].x_axis), i);
        update(y_axis, vec_error(actual[i].y_axis, expected[i].y_axis), i);
        update(z_axis, vec_error(actual[i].z_axis, expected[i].z_axis), i);
    }

    std::ostringstream out;
    out << std::setprecision(12) << "max origin error: " << origin.value
        << " at " << origin.index << "\nmax aim error: " << aim.value << " at "
        << aim.index << "\nmax x axis error: " << x_axis.value << " at "
        << x_axis.index << "\nmax y axis error: " << y_axis.value << " at "
        << y_axis.index << "\nmax z axis error: " << z_axis.value << " at "
        << z_axis.index;
    return out.str();
}

void expect_vec_near(glm::dvec3 const& actual,
                     glm::dvec3 const& expected,
                     double            tolerance = kTransformTolerance) {
    EXPECT_NEAR(actual.x, expected.x, tolerance);
    EXPECT_NEAR(actual.y, expected.y, tolerance);
    EXPECT_NEAR(actual.z, expected.z, tolerance);
}

void expect_scalar_near_or_both_nan(double actual,
                                    double expected,
                                    double tolerance = kValueTolerance) {
    if (std::isnan(actual) && std::isnan(expected)) return;

    EXPECT_NEAR(actual, expected, tolerance);
}

void expect_sun_near(SD::RaySource const& actual,
                     SD::RaySource const& expected) {
    EXPECT_EQ(actual.get_shape(), expected.get_shape());
    EXPECT_EQ(actual.get_gen_type(), expected.get_gen_type());
    expect_vec_near(actual.get_position(), expected.get_position());

    auto& actual_mutable   = const_cast<SD::RaySource&>(actual);
    auto& expected_mutable = const_cast<SD::RaySource&>(expected);

    expect_scalar_near_or_both_nan(actual_mutable.get_sigma(),
                                   expected_mutable.get_sigma());
    expect_scalar_near_or_both_nan(actual_mutable.get_half_width(),
                                   expected_mutable.get_half_width());
    expect_scalar_near_or_both_nan(actual_mutable.get_circumsolar_ratio(),
                                   expected_mutable.get_circumsolar_ratio());

    std::vector<double> actual_angle;
    std::vector<double> actual_intensity;
    std::vector<double> expected_angle;
    std::vector<double> expected_intensity;

    actual_mutable.get_user_data(actual_angle, actual_intensity);
    expected_mutable.get_user_data(expected_angle, expected_intensity);

    ASSERT_EQ(actual_angle.size(), expected_angle.size());
    ASSERT_EQ(actual_intensity.size(), expected_intensity.size());

    for (size_t i = 0; i < actual_angle.size(); ++i) {
        EXPECT_NEAR(actual_angle[i], expected_angle[i], kValueTolerance);
    }

    for (size_t i = 0; i < actual_intensity.size(); ++i) {
        EXPECT_NEAR(
            actual_intensity[i], expected_intensity[i], kValueTolerance);
    }
}

void expect_simulation_parameters_equal(SD::SimulationData const& actual,
                                        SD::SimulationData const& expected) {
    auto const& actual_params   = actual.get_simulation_parameters();
    auto const& expected_params = expected.get_simulation_parameters();

    EXPECT_EQ(actual_params.number_of_rays, expected_params.number_of_rays);
    EXPECT_EQ(actual_params.max_number_of_rays,
              expected_params.max_number_of_rays);
    EXPECT_NEAR(
        actual_params.tolerance, expected_params.tolerance, kValueTolerance);
    EXPECT_NEAR(
        actual_params.latitude, expected_params.latitude, kValueTolerance);
    EXPECT_NEAR(
        actual_params.longitude, expected_params.longitude, kValueTolerance);
    EXPECT_EQ(actual_params.seed, expected_params.seed);
    EXPECT_EQ(actual_params.include_sun_shape_errors,
              expected_params.include_sun_shape_errors);
    EXPECT_EQ(actual_params.include_optical_errors,
              expected_params.include_optical_errors);
}

void expect_snapshots_near(ElementSnapshot const& actual,
                           ElementSnapshot const& expected,
                           size_t                 index) {
    SCOPED_TRACE(mismatch_diagnostic(actual, expected));

    {
        SCOPED_TRACE("single element snapshot " + std::to_string(index));

        EXPECT_EQ(actual.name, expected.name);
        EXPECT_EQ(actual.enabled, expected.enabled);
        EXPECT_EQ(actual.aperture_json, expected.aperture_json);
        EXPECT_EQ(actual.surface_json, expected.surface_json);
        EXPECT_EQ(actual.optics_json, expected.optics_json);
    }

    {
        SCOPED_TRACE("single element snapshot origin " + std::to_string(index));
        expect_vec_near(actual.origin, expected.origin);
    }

    {
        SCOPED_TRACE("single element snapshot axis " + std::to_string(index));

        expect_vec_near(actual.aim_axis, expected.aim_axis);
        expect_vec_near(actual.x_axis, expected.x_axis);
        expect_vec_near(actual.y_axis, expected.y_axis);
        expect_vec_near(actual.z_axis, expected.z_axis);
    }
}

std::filesystem::path power_tower_surround_path() {
    return std::filesystem::path(SOLTRACE_REPO_ROOT) /
           "gui/assets/examples/Power-tower-surround.stinput";
}

void configure_tiny_deterministic_run(SD::SimulationData& data) {
    data.set_seed(12345);
    data.set_number_of_rays(5);
    data.set_max_rays_traced(500);

    auto& params                    = data.get_simulation_parameters();
    params.include_sun_shape_errors = false;
    params.include_optical_errors   = false;
}

std::unique_ptr<SolTrace::Result::SimulationResult>
run_native_trace(SD::SimulationData& data, bool disable_stages = false) {
    SolTrace::NativeRunner::NativeRunner runner;
    runner.set_number_of_threads(1);
    if (disable_stages) { runner.disable_stages(); }

    EXPECT_EQ(runner.initialize(), SolTrace::Runner::RunnerStatus::SUCCESS);
    EXPECT_EQ(runner.setup_simulation(&data),
              SolTrace::Runner::RunnerStatus::SUCCESS);
    EXPECT_EQ(runner.run_simulation(), SolTrace::Runner::RunnerStatus::SUCCESS);

    auto result = std::make_unique<SolTrace::Result::SimulationResult>();
    EXPECT_EQ(runner.report_simulation(result.get(), 0),
              SolTrace::Runner::RunnerStatus::SUCCESS);

    return result;
}

std::vector<SolTrace::Result::ray_record_ptr>
collect_ray_records(SolTrace::Result::SimulationResult const& result) {
    std::vector<SolTrace::Result::ray_record_ptr> records;
    records.reserve(result.get_number_of_records());

    for (auto iter = result.get_ray_record_iterator(); !result.is_at_end(iter);
         ++iter) {
        records.push_back(*iter);
    }

    return records;
}

std::map<SolTrace::Result::ray_id, SolTrace::Result::ray_record_ptr>
collect_ray_records_by_id(SolTrace::Result::SimulationResult const& result) {
    std::map<SolTrace::Result::ray_id, SolTrace::Result::ray_record_ptr>
        records;

    for (auto iter = result.get_ray_record_iterator(); !result.is_at_end(iter);
         ++iter) {
        records.emplace((*iter)->id, *iter);
    }

    return records;
}

void expect_interaction_near(
    SolTrace::Result::InteractionRecord const& actual,
    SolTrace::Result::InteractionRecord const& expected,
    SolTrace::Result::ray_id                   ray_id,
    size_t                                     event_index) {
    SCOPED_TRACE("ray id " + std::to_string(ray_id) + ", event " +
                 std::to_string(event_index));

    EXPECT_EQ(actual.event, expected.event);
    expect_vec_near(actual.location, expected.location, 1.0e-6);
    expect_vec_near(actual.direction, expected.direction, 1.0e-9);
}

void expect_ray_record_near(SolTrace::Result::RayRecord const& actual,
                            SolTrace::Result::RayRecord const& expected,
                            SolTrace::Result::ray_id           ray_id) {
    SCOPED_TRACE("ray record id " + std::to_string(ray_id));

    EXPECT_EQ(actual.id, expected.id);
    ASSERT_EQ(actual.interactions.size(), expected.interactions.size());

    for (size_t i = 0; i < expected.interactions.size(); ++i) {
        expect_interaction_near(
            *actual.interactions[i], *expected.interactions[i], ray_id, i);
    }
}

std::string sun_box_summary(SolTrace::Result::SimulationResult& actual,
                            SolTrace::Result::SimulationResult& expected) {
    double actual_width    = 0.0;
    double actual_height   = 0.0;
    double expected_width  = 0.0;
    double expected_height = 0.0;

    actual.get_sun_dimensions(actual_width, actual_height);
    expected.get_sun_dimensions(expected_width, expected_height);

    std::ostringstream out;
    out << std::setprecision(12) << "actual sun box: width=" << actual_width
        << ", height=" << actual_height << ", area=" << actual.get_sun_A_box()
        << "\nexpected sun box: width=" << expected_width
        << ", height=" << expected_height
        << ", area=" << expected.get_sun_A_box();
    return out.str();
}

::testing::AssertionResult export_succeeded(
    Result<std::shared_ptr<db::DatabaseExport>, QString> const& result) {
    if (result) return ::testing::AssertionSuccess();

    return ::testing::AssertionFailure()
           << result.get_failure().toStdString();
}

void expect_sun_box_near(SolTrace::Result::SimulationResult& actual,
                         SolTrace::Result::SimulationResult& expected,
                         bool compare_sun_ray_count = true) {
    SCOPED_TRACE(sun_box_summary(actual, expected));

    if (compare_sun_ray_count) {
        EXPECT_EQ(actual.get_sun_ray_count(), expected.get_sun_ray_count());
    }

    double actual_width    = 0.0;
    double actual_height   = 0.0;
    double expected_width  = 0.0;
    double expected_height = 0.0;

    actual.get_sun_dimensions(actual_width, actual_height);
    expected.get_sun_dimensions(expected_width, expected_height);

    EXPECT_NEAR(actual_width, expected_width, 1.0e-7);
    EXPECT_NEAR(actual_height, expected_height, 1.0e-7);
    EXPECT_NEAR(actual.get_sun_A_box(), expected.get_sun_A_box(), 1.0e-4);
}

} // namespace

TEST(DatabaseRoundTrip, PowerTowerSurroundExportsEquivalentGlobalSimData) {
    SD::SimulationData original;
    ASSERT_TRUE(
        original.import_from_file(power_tower_surround_path().string()));

    db::Database database("round-trip");
    database.import(original, true /* legacy_stinput_cylinder_origins */);

    auto exported_result = database.export_to_simdata();
    ASSERT_TRUE(export_succeeded(exported_result));

    auto exported = exported_result.get_success();
    ASSERT_NE(exported, nullptr);
    ASSERT_NE(exported->data, nullptr);

    expect_simulation_parameters_equal(*exported->data, original);

    ASSERT_EQ(exported->data->get_number_of_ray_sources(),
              original.get_number_of_ray_sources());
    expect_sun_near(*exported->data->get_ray_source(),
                    *original.get_ray_source());

    auto actual_snapshots   = collect_single_element_snapshots(*exported->data);
    auto expected_snapshots = collect_single_element_snapshots(
        original, true /* legacy_stinput_cylinder_origin */);

    ASSERT_EQ(actual_snapshots.size(), expected_snapshots.size());

    SCOPED_TRACE(max_error_summary(actual_snapshots, expected_snapshots));

    for (size_t i = 0; i < expected_snapshots.size(); ++i) {
        expect_snapshots_near(actual_snapshots[i], expected_snapshots[i], i);
    }
}

TEST(DatabaseRoundTrip, PowerTowerSurroundNativeTraceMatchesOriginalSimData) {
    SD::SimulationData original;
    ASSERT_TRUE(
        original.import_from_file(power_tower_surround_path().string()));
    configure_tiny_deterministic_run(original);

    SD::SimulationData source_for_database;
    ASSERT_TRUE(source_for_database.import_from_file(
        power_tower_surround_path().string()));

    db::Database database("round-trip-trace");
    database.import(source_for_database,
                    true /* legacy_stinput_cylinder_origins */);

    auto exported_result = database.export_to_simdata();
    ASSERT_TRUE(export_succeeded(exported_result));

    auto exported = exported_result.get_success();
    ASSERT_NE(exported, nullptr);
    ASSERT_NE(exported->data, nullptr);

    configure_tiny_deterministic_run(*exported->data);

    SD::SimulationData source_for_expected;
    ASSERT_TRUE(source_for_expected.import_from_file(
        power_tower_surround_path().string()));

    db::Database expected_database("round-trip-trace-expected");
    expected_database.import(source_for_expected,
                             true /* legacy_stinput_cylinder_origins */);

    auto expected_exported_result = expected_database.export_to_simdata();
    ASSERT_TRUE(export_succeeded(expected_exported_result));

    auto expected_exported = expected_exported_result.get_success();
    ASSERT_NE(expected_exported, nullptr);
    ASSERT_NE(expected_exported->data, nullptr);
    configure_tiny_deterministic_run(*expected_exported->data);

    auto staged_result   = run_native_trace(original);
    auto expected_result = run_native_trace(*expected_exported->data);
    auto actual_result   = run_native_trace(*exported->data);

    expect_sun_box_near(
        *actual_result, *staged_result, false /* compare_sun_ray_count */);
    expect_sun_box_near(*actual_result, *expected_result);

    auto expected_records = collect_ray_records_by_id(*expected_result);
    auto actual_records   = collect_ray_records_by_id(*actual_result);

    ASSERT_FALSE(expected_records.empty());
    ASSERT_EQ(actual_records.size(), expected_records.size());

    size_t compared = 0;
    for (auto const& [ray_id, expected_record] : expected_records) {
        auto actual_iter = actual_records.find(ray_id);
        ASSERT_NE(actual_iter, actual_records.end())
            << "missing actual ray id " << ray_id;

        expect_ray_record_near(*actual_iter->second, *expected_record, ray_id);

        compared++;
        if (compared >= 5) break;
    }
}

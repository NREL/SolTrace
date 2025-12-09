#include <gtest/gtest.h>

#include <aperture.hpp>
#include <surface.hpp>
#include <constants.hpp>
#include <native_runner.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result_export.hpp>
#include <json_helpers.hpp>

#include "common.hpp"

namespace SolTrace { namespace Data { /* forward declare in case of missing header resolution */ surface_ptr make_surface_from_json(const nlohmann::ordered_json&); } }

void get_default_element_base(nlohmann::ordered_json& jnode)
{
    jnode["active"] = true;
    jnode["virtual_flag"] = false;
    jnode["my_id"] = 1;
    jnode["my_name"] = "";
    jnode["stage"] = 0;
    jnode["origin"] = Vector3d(0, 0, 0).data;
    jnode["aim"] = Vector3d(0, 0, 0).data;
    jnode["zrot"] = 0;
}

TEST(io_json, json_round_trip)
{
    namespace fs = std::filesystem;
    using json = nlohmann::ordered_json;

    // Build paths
    const fs::path project_root(PROJECT_DIR);
    const fs::path sample_path = project_root / "High Flux Solar Furnace.stinput";
    const fs::path output_path_1 = project_root / "json_test_1.json";
    const fs::path output_path_2 = project_root / "json_test_2.json";

    ASSERT_TRUE(fs::exists(sample_path)) << "Sample .stinput not found: " << sample_path.string();

    // Load original simulation from .stinput
    SimulationData sd_original;
    ASSERT_TRUE(sd_original.import_from_file(sample_path.string())) << "Failed to import stinput";

    // Export to JSON
    ASSERT_NO_THROW(sd_original.export_json_file(output_path_1.string())) << "Failed to export first JSON";

    // Import from first JSON
    SimulationData sd_round_trip;
    ASSERT_NO_THROW(sd_round_trip.import_json_file(output_path_1.string())) << "Failed to import JSON";

    // Re-export second JSON
    ASSERT_NO_THROW(sd_round_trip.export_json_file(output_path_2.string())) << "Failed to export second JSON";

    // Parse both JSON files
    std::ifstream ifs1(output_path_1);
    std::ifstream ifs2(output_path_2);
    ASSERT_TRUE(ifs1.is_open()) << "Cannot open first JSON file";
    ASSERT_TRUE(ifs2.is_open()) << "Cannot open second JSON file";

    json root1;
    json root2;
    try {
        root1 = json::parse(ifs1);
        root2 = json::parse(ifs2);
    }
    catch (const std::exception& e) {
        FAIL() << "JSON parse failure: " << e.what();
    }

    // Basic schema/version checks
    ASSERT_TRUE(root1.contains("schema_version"));
    ASSERT_TRUE(root2.contains("schema_version"));
    EXPECT_EQ(root1["schema_version"], root2["schema_version"]);

    // Simulation parameters
    ASSERT_TRUE(root1.contains("simulation_parameters"));
    ASSERT_TRUE(root2.contains("simulation_parameters"));
    EXPECT_EQ(root1["simulation_parameters"], root2["simulation_parameters"]);

    // Ray sources
    ASSERT_TRUE(root1.contains("ray_sources"));
    ASSERT_TRUE(root2.contains("ray_sources"));
    EXPECT_EQ(root1["ray_sources"], root2["ray_sources"]);

    // Elements
    ASSERT_TRUE(root1.contains("elements"));
    ASSERT_TRUE(root2.contains("elements"));
    EXPECT_EQ(root1["elements"], root2["elements"]);

    // Global sanity checks against SimulationData instances
    EXPECT_EQ(root1["number_of_elements"], root2["number_of_elements"]);
    EXPECT_EQ(sd_original.get_number_of_elements(), sd_round_trip.get_number_of_elements());
    EXPECT_EQ(sd_original.get_number_of_ray_sources(), sd_round_trip.get_number_of_ray_sources());

    // Full structural equality
    ASSERT_TRUE(root1 == root2);

    // Conditional cleanup: remove only if test passed so far.
    if (!::testing::Test::HasFailure()) {
        std::error_code ec;
        fs::remove(output_path_1, ec);
        fs::remove(output_path_2, ec);
    }
}

TEST(io_json, large_field_comparison)
{
    namespace fs = std::filesystem;
    using json = nlohmann::ordered_json;

    // Build paths
    const fs::path project_root(PROJECT_DIR);
    const fs::path sample_path = project_root / "Power-tower-surround_singlefacet_large_afternoon.stinput";
    const fs::path output_path = project_root / "number_elements_test.json";

    ASSERT_TRUE(fs::exists(sample_path)) << "Sample .stinput not found: " << sample_path.string();

    // Load original simulation from .stinput
    SimulationData sd_original;
    ASSERT_TRUE(sd_original.import_from_file(sample_path.string())) << "Failed to import stinput";

    // Export to JSON
    ASSERT_NO_THROW(sd_original.export_json_file(output_path.string())) << "Failed to export first JSON";

    // Import from first JSON
    SimulationData sd_round_trip;
    ASSERT_NO_THROW(sd_round_trip.import_json_file(output_path.string())) << "Failed to import JSON";

    // Compare number of elements
    int N_elements_original = sd_original.get_number_of_elements();
    int N_elements_round_trip = sd_round_trip.get_number_of_elements();
    ASSERT_TRUE(N_elements_original == N_elements_round_trip) << "Element number is not equal";

    // Compare number of sources
    int N_sources_original = sd_original.get_number_of_ray_sources();
    int N_sources_round_trip = sd_round_trip.get_number_of_ray_sources();
    ASSERT_TRUE(N_sources_original == N_sources_round_trip) << "Ray sources number is not equal";

    // Compare number of rays
    int N_rays_original = sd_original.get_number_of_rays();
    int N_rays_round_trip = sd_round_trip.get_number_of_rays();
    ASSERT_TRUE(N_rays_original == N_rays_round_trip) << "Ray number is not equal";

    // Conditional cleanup: remove only if test passed so far.
    if (!::testing::Test::HasFailure()) {
        std::error_code ec;
        fs::remove(output_path, ec);
    }
}

TEST(io_json, precision_comparison)
{
    namespace fs = std::filesystem;
    using json = nlohmann::ordered_json;

    // Build paths
    const fs::path project_root(PROJECT_DIR);
    const fs::path output_path = project_root / "precision_comparison_test.json";

    // Make simulation data
    SimulationData sd_original;
    SimulationParameters& simpar_original = sd_original.get_simulation_parameters();
    simpar_original.latitude = SolTrace::Data::PI / 2.;
    simpar_original.longitude = 1.5816981651658435135814384351384351385143845;
    simpar_original.tolerance = 1e-90;
    simpar_original.seed = static_cast<int>(15654681468168136541ull);

    // Add sun
    auto sun = make_ray_source<Sun>();
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    sun->set_shape(SolTrace::Data::SunShape::GAUSSIAN, 1.0, nan, nan, {}, {});
    sd_original.add_ray_source(sun);

    // Export to JSON
    ASSERT_NO_THROW(sd_original.export_json_file(output_path.string())) << "Failed to export first JSON";

    // Import from first JSON
    SimulationData sd_round_trip;
    ASSERT_NO_THROW(sd_round_trip.import_json_file(output_path.string())) << "Failed to import JSON";

    // Compare sim parameters
    SimulationParameters& simpar_round_trip = sd_round_trip.get_simulation_parameters();
    EXPECT_DOUBLE_EQ(simpar_original.latitude, simpar_round_trip.latitude);
    EXPECT_DOUBLE_EQ(simpar_original.longitude, simpar_round_trip.longitude);
    EXPECT_DOUBLE_EQ(simpar_original.tolerance, simpar_round_trip.tolerance);
    EXPECT_EQ(simpar_original.seed, simpar_round_trip.seed);

    // Compare sun nan parameters
    double half_width_original = sd_original.get_ray_source()->get_half_width();
    double half_width_round_trip = sd_round_trip.get_ray_source()->get_half_width();
    double csr_original = sd_original.get_ray_source()->get_circumsolar_ratio();
    double csr_round_trip = sd_round_trip.get_ray_source()->get_circumsolar_ratio();

    EXPECT_TRUE(std::isnan(half_width_original) && std::isnan(half_width_round_trip));
    EXPECT_TRUE(std::isnan(csr_original) && std::isnan(csr_round_trip));

    // Conditional cleanup: remove only if test passed so far.
    if (!::testing::Test::HasFailure()) {
        std::error_code ec;
        fs::remove(output_path, ec);
    }
}

TEST(io_json, empty_case)
{
    namespace fs = std::filesystem;
    using json = nlohmann::ordered_json;

    // Build paths
    const fs::path project_root(PROJECT_DIR);
    const fs::path output_path = project_root / "empty_case.json";

    // Make simulation data
    SimulationData sd_original;

    // Export
    ASSERT_NO_THROW(sd_original.export_json_file(output_path.string())) << "Failed to export first JSON";

    // Import from first JSON
    SimulationData sd_round_trip;
    ASSERT_NO_THROW(sd_round_trip.import_json_file(output_path.string())) << "Failed to import JSON";

    // Check number of elements
    EXPECT_TRUE(sd_round_trip.get_number_of_elements() == 0);
    EXPECT_TRUE(sd_round_trip.get_number_of_ray_sources() == 0);

    // Conditional cleanup: remove only if test passed so far.
    if (!::testing::Test::HasFailure()) {
        std::error_code ec;
        fs::remove(output_path, ec);
    }
}

TEST(io_json, invalid_sun_shape)
{
    namespace fs = std::filesystem;

    // Build paths
    const fs::path project_root(PROJECT_DIR);
    const fs::path output_path = project_root / "invalid_sun_shape.json";

    // Make simulation data
    SimulationData sd;
    auto sun = SolTrace::Data::make_ray_source<Sun>();
    sun->set_shape(SunShape::GAUSSIAN, 0.01, 0.0, 0.0);
    sd.add_ray_source(sun);

    // Export JSON
    sd.export_json_file(output_path.string());

    // Tamper JSON
    nlohmann::ordered_json root;
    {
        std::ifstream ifs(output_path);
        ifs >> root;
    }
    root["ray_sources"]["0"]["my_shape"] = "Ellipse"; // invalid
    {
        std::ofstream ofs(output_path, std::ios::trunc);
        ofs << root.dump(SolTrace::Data::kJsonIndentSpaces);
    }

    SimulationData sd2;
    EXPECT_THROW(sd2.import_json_file(output_path.string()), std::runtime_error);

    // Conditional cleanup: remove only if test passed so far.
    if (!::testing::Test::HasFailure()) {
        std::error_code ec;
        fs::remove(output_path, ec);
    }
}

TEST(io_json, invalid_source_type)
{
    namespace fs = std::filesystem;

    // Build paths
    const fs::path project_root(PROJECT_DIR);
    const fs::path output_path = project_root / "invalid_source_type.json";

    SimulationData sd;
    auto sun = SolTrace::Data::make_ray_source<Sun>();
    sun->set_shape(SunShape::GAUSSIAN, 0.01, 0.0, 0.0);
    sd.add_ray_source(sun);

    // Export JSON
    sd.export_json_file(output_path.string());

    // Tamper JSON
    nlohmann::ordered_json root;
    {
        std::ifstream ifs(output_path.string());
        ifs >> root;
    }
    root["ray_sources"]["0"]["source_type"] = "The Moon"; // invalid
    {
        std::ofstream ofs(output_path.string(), std::ios::trunc);
        ofs << root.dump(2);
    }

    SimulationData sd2;
    EXPECT_THROW(sd2.import_json_file(output_path.string()), std::runtime_error);

    // Conditional cleanup: remove only if test passed so far.
    if (!::testing::Test::HasFailure()) {
        std::error_code ec;
        fs::remove(output_path, ec);
    }
}

TEST(io_json, multi_ray_source)
{
    namespace fs = std::filesystem;
    using json = nlohmann::ordered_json;

    // Build paths
    const fs::path project_root(PROJECT_DIR);
    const fs::path output_path = project_root / "multi_ray_source.json";

    // Make simulation data
    SimulationData sd_original;

    // Add sun
    auto sun = make_ray_source<Sun>();
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    sun->set_shape(SolTrace::Data::SunShape::GAUSSIAN, 1.0, nan, nan, {}, {});
    sd_original.add_ray_source(sun);

    auto sun2 = make_ray_source<Sun>();
    sun2->set_shape(SolTrace::Data::SunShape::GAUSSIAN, 1.0, nan, nan, {}, {});
    sd_original.add_ray_source(sun2);

    // Export to JSON
    ASSERT_NO_THROW(sd_original.export_json_file(output_path.string())) << "Failed to export first JSON";

    // Import from first JSON
    SimulationData sd_round_trip;
    ASSERT_NO_THROW(sd_round_trip.import_json_file(output_path.string())) << "Failed to import JSON";

    // Check ray sources
    EXPECT_EQ(sd_original.get_number_of_ray_sources(), sd_round_trip.get_number_of_ray_sources());
    EXPECT_EQ(sd_round_trip.get_number_of_ray_sources(), 2);

    // Conditional cleanup: remove only if test passed so far.
    if (!::testing::Test::HasFailure()) {
        std::error_code ec;
        fs::remove(output_path, ec);
    }
}

TEST(io_json, performance_comparison)
{
    namespace fs = std::filesystem;
    using json = nlohmann::ordered_json;
    using SolTrace::Runner::RunnerStatus;
    using SolTrace::NativeRunner::NativeRunner;
    using SolTrace::NativeRunner::TRayData;
    using SolTrace::NativeRunner::tstage_ptr;
    using SolTrace::NativeRunner::TSystem;

    // Build paths
    const fs::path project_root(PROJECT_DIR);
    const fs::path sample_path = project_root / "High Flux Solar Furnace.stinput";
    const fs::path output_path_1 = project_root / "json_test_1.json";
    const fs::path output_path_2 = project_root / "json_test_2.json";

    ASSERT_TRUE(fs::exists(sample_path)) << "Sample .stinput not found: " << sample_path.string();

    // Load original simulation from .stinput
    SimulationData sd_original;
    ASSERT_TRUE(sd_original.import_from_file(sample_path.string())) << "Failed to import stinput";

    // Export to JSON
    ASSERT_NO_THROW(sd_original.export_json_file(output_path_1.string())) << "Failed to export first JSON";

    // Import from first JSON
    SimulationData sd_round_trip;
    ASSERT_NO_THROW(sd_round_trip.import_json_file(output_path_1.string())) << "Failed to import JSON";

    // Run original case
    NativeRunner runner;
    RunnerStatus sts;
    sts = runner.initialize();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS) << "runner.initialize() failed";
    sts = runner.setup_simulation(&sd_original);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS) << "runner.setup_simulation() failed";
    sts = runner.run_simulation();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS) << "runner.run_simulation() failed";
    SimulationResult result_original;
    sts = runner.report_simulation(&result_original, 0);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS) << "runner.report_simulation() failed";

    // Run round trip case
    NativeRunner runner_round_trip;
    sts = runner_round_trip.initialize();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS) << "runner_round_trip.initialize() failed";
    sts = runner_round_trip.setup_simulation(&sd_round_trip);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS) << "runner_round_trip.setup_simulation() failed";
    sts = runner_round_trip.run_simulation();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS) << "runner_round_trip.run_simulation() failed";
    SimulationResult result_round_trip;
    sts = runner_round_trip.report_simulation(&result_round_trip, 0);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS) << "runner_round_trip.report_simulation() failed";
    // Compare number of records
    ASSERT_EQ(result_original.get_number_of_records(), result_round_trip.get_number_of_records());

    // Loop through each result set and compare
    for (int i = 0; i < result_original.get_number_of_records(); ++i)
    {
        ray_record_ptr rr_o = result_original[i];
        ray_record_ptr rr_r = result_round_trip[i];

        ASSERT_EQ(rr_o->get_number_of_interactions(),
            rr_r->get_number_of_interactions());

        for (int_fast64_t k = 0; k < rr_o->get_number_of_interactions(); ++k)
        {
            // Event & element IDs
            ASSERT_EQ(rr_o->get_event(k), rr_r->get_event(k));
            ASSERT_EQ(rr_o->get_element(k), rr_r->get_element(k));

            // Positions
            Vector3d pos_o; Vector3d pos_r;
            rr_o->get_position(k, pos_o);
            rr_r->get_position(k, pos_r);
            EXPECT_DOUBLE_EQ(pos_o[0], pos_r[0]);
            EXPECT_DOUBLE_EQ(pos_o[1], pos_r[1]);
            EXPECT_DOUBLE_EQ(pos_o[2], pos_r[2]);

            // Directions
            Vector3d dir_o; Vector3d dir_r;
            rr_o->get_direction(k, dir_o);
            rr_r->get_direction(k, dir_r);
            EXPECT_DOUBLE_EQ(dir_o[0], dir_r[0]);
            EXPECT_DOUBLE_EQ(dir_o[1], dir_r[1]);
            EXPECT_DOUBLE_EQ(dir_o[2], dir_r[2]);
        }

    }

    // Conditional cleanup: remove only if test passed so far.
    if (!::testing::Test::HasFailure()) {
        std::error_code ec;
        fs::remove(output_path_1, ec);
        fs::remove(output_path_2, ec);
    }
}

TEST(io_json, missing_key)
{
    namespace fs = std::filesystem;

    // Build paths
    const fs::path project_root(PROJECT_DIR);
    const fs::path output_path = project_root / "missing_key.json";

    // Make simulation data
    SimulationData sd;

    // Export JSON
    sd.export_json_file(output_path.string());

    // Tamper JSON
    nlohmann::ordered_json root;
    {
        std::ifstream ifs(output_path);
        ifs >> root;
    }
    root["simulation_parameters"].erase("include_sun_shape_errors");
    {
        std::ofstream ofs(output_path, std::ios::trunc);
        ofs << root.dump(SolTrace::Data::kJsonIndentSpaces);
    }

    SimulationData sd2;
    EXPECT_THROW(sd2.import_json_file(output_path.string()), nlohmann::json_abi_v3_11_3::detail::out_of_range);

    // Conditional cleanup: remove only if test passed so far.
    if (!::testing::Test::HasFailure()) {
        std::error_code ec;
        fs::remove(output_path, ec);
    }
}

TEST(io_json, apertures_read)
{
    using json = nlohmann::ordered_json;

    // ANNULUS
    json jannulus;
    jannulus["aperture_type"] = SolTrace::Data::ApertureTypeMap.at(ApertureType::ANNULUS);
    jannulus["inner_radius"] = 1;
    jannulus["outer_radius"] = 2;
    jannulus["arc_angle"] = 1;
    EXPECT_NO_THROW(Aperture::make_aperture_from_json(jannulus));
    auto ann_ptr = Aperture::make_aperture_from_json(jannulus);
    auto ann_cast = dynamic_cast<Annulus*>(ann_ptr.get());
    ASSERT_TRUE(ann_cast != nullptr);
    EXPECT_DOUBLE_EQ(1, ann_cast->inner_radius);
    EXPECT_DOUBLE_EQ(2, ann_cast->outer_radius);
    EXPECT_DOUBLE_EQ(1, ann_cast->arc_angle);

    // CIRCLE
    json jcircle;
    jcircle["aperture_type"] = SolTrace::Data::ApertureTypeMap.at(ApertureType::CIRCLE);
    jcircle["diameter"] = 1;
    EXPECT_NO_THROW(Aperture::make_aperture_from_json(jcircle));
    auto circ_ptr = Aperture::make_aperture_from_json(jcircle);
    auto circ_cast = dynamic_cast<Circle*>(circ_ptr.get());
    ASSERT_TRUE(circ_cast != nullptr);
    EXPECT_DOUBLE_EQ(1, circ_cast->diameter);

    // HEXAGON
    json jhexagon;
    jhexagon["aperture_type"] = SolTrace::Data::ApertureTypeMap.at(ApertureType::HEXAGON);
    jhexagon["circumscribe_diameter"] = 3;
    EXPECT_NO_THROW(Aperture::make_aperture_from_json(jhexagon));
    auto hex_ptr = Aperture::make_aperture_from_json(jhexagon);
    auto hex_cast = dynamic_cast<Hexagon*>(hex_ptr.get());
    ASSERT_TRUE(hex_cast != nullptr);
    EXPECT_DOUBLE_EQ(3, hex_cast->circumscribe_diameter);

    // RECTANGLE
    json jrectangle;
    jrectangle["aperture_type"] = SolTrace::Data::ApertureTypeMap.at(ApertureType::RECTANGLE);
    jrectangle["x_length"] = 4;
    jrectangle["y_length"] = 5;
    jrectangle["x_coord"] = -2;
    jrectangle["y_coord"] = -2.5;
    EXPECT_NO_THROW(Aperture::make_aperture_from_json(jrectangle));
    auto rect_ptr = Aperture::make_aperture_from_json(jrectangle);
    auto rect_cast = dynamic_cast<Rectangle*>(rect_ptr.get());
    ASSERT_TRUE(rect_cast != nullptr);
    EXPECT_DOUBLE_EQ(4, rect_cast->x_length);
    EXPECT_DOUBLE_EQ(5, rect_cast->y_length);
    EXPECT_DOUBLE_EQ(-2, rect_cast->x_coord);
    EXPECT_DOUBLE_EQ(-2.5, rect_cast->y_coord);

    // EQUILATERAL_TRIANGLE
    json jtriangle_eq;
    jtriangle_eq["aperture_type"] = SolTrace::Data::ApertureTypeMap.at(ApertureType::EQUILATERAL_TRIANGLE);
    jtriangle_eq["circumscribe_diameter"] = 6;
    EXPECT_NO_THROW(Aperture::make_aperture_from_json(jtriangle_eq));
    auto eq_ptr = Aperture::make_aperture_from_json(jtriangle_eq);
    auto eq_cast = dynamic_cast<EqualateralTriangle*>(eq_ptr.get());
    ASSERT_TRUE(eq_cast != nullptr);
    EXPECT_DOUBLE_EQ(6, eq_cast->circumscribe_diameter);

    // IRREGULAR_TRIANGLE
    json jtriangle_ir;
    jtriangle_ir["aperture_type"] = SolTrace::Data::ApertureTypeMap.at(ApertureType::IRREGULAR_TRIANGLE);
    jtriangle_ir["x1"] = 0; jtriangle_ir["y1"] = 0;
    jtriangle_ir["x2"] = 1; jtriangle_ir["y2"] = 0;
    jtriangle_ir["x3"] = 0; jtriangle_ir["y3"] = 1;
    EXPECT_NO_THROW(Aperture::make_aperture_from_json(jtriangle_ir));
    auto irtri_ptr = Aperture::make_aperture_from_json(jtriangle_ir);
    auto irtri_cast = dynamic_cast<IrregularTriangle*>(irtri_ptr.get());
    ASSERT_TRUE(irtri_cast != nullptr);
    EXPECT_DOUBLE_EQ(0, irtri_cast->x1);
    EXPECT_DOUBLE_EQ(0, irtri_cast->y1);
    EXPECT_DOUBLE_EQ(1, irtri_cast->x2);
    EXPECT_DOUBLE_EQ(0, irtri_cast->y2);
    EXPECT_DOUBLE_EQ(0, irtri_cast->x3);
    EXPECT_DOUBLE_EQ(1, irtri_cast->y3);

    // IRREGULAR_QUADRILATERAL
    json jquad_ir;
    jquad_ir["aperture_type"] = SolTrace::Data::ApertureTypeMap.at(ApertureType::IRREGULAR_QUADRILATERAL);
    jquad_ir["x1"] = 0; jquad_ir["y1"] = 0;
    jquad_ir["x2"] = 2; jquad_ir["y2"] = 0;
    jquad_ir["x3"] = 2; jquad_ir["y3"] = 1;
    jquad_ir["x4"] = 0; jquad_ir["y4"] = 1;
    EXPECT_NO_THROW(Aperture::make_aperture_from_json(jquad_ir));
    auto irquad_ptr = Aperture::make_aperture_from_json(jquad_ir);
    auto irquad_cast = dynamic_cast<IrregularQuadrilateral*>(irquad_ptr.get());
    ASSERT_TRUE(irquad_cast != nullptr);
    EXPECT_DOUBLE_EQ(0, irquad_cast->x1);
    EXPECT_DOUBLE_EQ(0, irquad_cast->y1);
    EXPECT_DOUBLE_EQ(2, irquad_cast->x2);
    EXPECT_DOUBLE_EQ(0, irquad_cast->y2);
    EXPECT_DOUBLE_EQ(2, irquad_cast->x3);
    EXPECT_DOUBLE_EQ(1, irquad_cast->y3);
    EXPECT_DOUBLE_EQ(0, irquad_cast->x4);
    EXPECT_DOUBLE_EQ(1, irquad_cast->y4);

    // No aperture type
    json jmissing;
    jmissing["radius"] = 0;
    EXPECT_THROW(Aperture::make_aperture_from_json(jmissing), std::invalid_argument);

    // Unknown aperture type
    json junknown;
    junknown["aperture_type"] = SolTrace::Data::ApertureTypeMap.at(ApertureType::APERTURE_UNKNOWN);
    EXPECT_THROW(Aperture::make_aperture_from_json(junknown), std::invalid_argument);

    int x = 0;
}

TEST(io_json, apertures_write)
{
    using json = nlohmann::ordered_json;

    // ANNULUS
    json jannulus;
    double ri = 1;
    double ro = 2;
    double arc = 1;
    auto annulus = make_aperture<Annulus>(ri, ro, arc);
    ASSERT_NO_THROW(annulus->write_json(jannulus));
    EXPECT_DOUBLE_EQ(ri, jannulus["inner_radius"]);
    EXPECT_DOUBLE_EQ(ro, jannulus["outer_radius"]);
    EXPECT_DOUBLE_EQ(arc, jannulus["arc_angle"]);
    EXPECT_TRUE(jannulus["aperture_type"] == SolTrace::Data::ApertureTypeMap.at(ApertureType::ANNULUS));

    // CIRCLE
    json jcircle;
    double cdiam = 1;
    auto circle = make_aperture<Circle>(cdiam);
    ASSERT_NO_THROW(circle->write_json(jcircle));
    EXPECT_DOUBLE_EQ(cdiam, jcircle["diameter"]);
    EXPECT_TRUE(jcircle["aperture_type"] == SolTrace::Data::ApertureTypeMap.at(ApertureType::CIRCLE));

    // HEXAGON
    json jhexagon;
    double hdiam = 3;
    auto hexagon = make_aperture<Hexagon>(hdiam);
    ASSERT_NO_THROW(hexagon->write_json(jhexagon));
    EXPECT_DOUBLE_EQ(hdiam, jhexagon["circumscribe_diameter"]);
    EXPECT_TRUE(jhexagon["aperture_type"] == SolTrace::Data::ApertureTypeMap.at(ApertureType::HEXAGON));

    // RECTANGLE
    json jrectangle;
    double xlen = 4;
    double ylen = 5;
    double xl = -2;
    double yl = -2.5;
    auto rectangle = make_aperture<Rectangle>(xlen, ylen, xl, yl);
    ASSERT_NO_THROW(rectangle->write_json(jrectangle));
    EXPECT_DOUBLE_EQ(xlen, jrectangle["x_length"]);
    EXPECT_DOUBLE_EQ(ylen, jrectangle["y_length"]);
    EXPECT_DOUBLE_EQ(xl, jrectangle["x_coord"]);
    EXPECT_DOUBLE_EQ(yl, jrectangle["y_coord"]);
    EXPECT_TRUE(jrectangle["aperture_type"] == SolTrace::Data::ApertureTypeMap.at(ApertureType::RECTANGLE));

    // EQUILATERAL_TRIANGLE
    json jtriangle_eq;
    double eqdiam = 6;
    auto triangle_eq = make_aperture<EqualateralTriangle>(eqdiam);
    ASSERT_NO_THROW(triangle_eq->write_json(jtriangle_eq));
    EXPECT_DOUBLE_EQ(eqdiam, jtriangle_eq["circumscribe_diameter"]);
    EXPECT_TRUE(jtriangle_eq["aperture_type"] == SolTrace::Data::ApertureTypeMap.at(ApertureType::EQUILATERAL_TRIANGLE));

    // IRREGULAR_TRIANGLE
    json jtriangle_ir;
    double t1x = 0; double t1y = 0;
    double t2x = 1; double t2y = 0;
    double t3x = 0; double t3y = 1;
    auto triangle_ir = make_aperture<IrregularTriangle>(t1x, t1y, t2x, t2y, t3x, t3y);
    ASSERT_NO_THROW(triangle_ir->write_json(jtriangle_ir));
    EXPECT_DOUBLE_EQ(t1x, jtriangle_ir["x1"]);
    EXPECT_DOUBLE_EQ(t1y, jtriangle_ir["y1"]);
    EXPECT_DOUBLE_EQ(t2x, jtriangle_ir["x2"]);
    EXPECT_DOUBLE_EQ(t2y, jtriangle_ir["y2"]);
    EXPECT_DOUBLE_EQ(t3x, jtriangle_ir["x3"]);
    EXPECT_DOUBLE_EQ(t3y, jtriangle_ir["y3"]);
    EXPECT_TRUE(jtriangle_ir["aperture_type"] == SolTrace::Data::ApertureTypeMap.at(ApertureType::IRREGULAR_TRIANGLE));

    // IRREGULAR_QUADRILATERAL
    json jquad_ir;
    double q1x = 0; double q1y = 0;
    double q2x = 2; double q2y = 0;
    double q3x = 2; double q3y = 1;
    double q4x = 0; double q4y = 1;
    auto quad_ir = make_aperture<IrregularQuadrilateral>(q1x, q1y, q2x, q2y, q3x, q3y, q4x, q4y);
    ASSERT_NO_THROW(quad_ir->write_json(jquad_ir));
    EXPECT_DOUBLE_EQ(q1x, jquad_ir["x1"]);
    EXPECT_DOUBLE_EQ(q1y, jquad_ir["y1"]);
    EXPECT_DOUBLE_EQ(q2x, jquad_ir["x2"]);
    EXPECT_DOUBLE_EQ(q2y, jquad_ir["y2"]);
    EXPECT_DOUBLE_EQ(q3x, jquad_ir["x3"]);
    EXPECT_DOUBLE_EQ(q3y, jquad_ir["y3"]);
    EXPECT_DOUBLE_EQ(q4x, jquad_ir["x4"]);
    EXPECT_DOUBLE_EQ(q4y, jquad_ir["y4"]);
    EXPECT_TRUE(jquad_ir["aperture_type"] == SolTrace::Data::ApertureTypeMap.at(ApertureType::IRREGULAR_QUADRILATERAL));
}

TEST(io_json, surface_read)
{
    using json = nlohmann::ordered_json;

    // CONE
    json jcone;
    jcone["surface_type"] = SolTrace::Data::SurfaceTypeMap.at(SolTrace::Data::CONE);
    jcone["half_angle"] = 1;
    EXPECT_NO_THROW(SolTrace::Data::make_surface_from_json(jcone));
    auto cone_ptr = SolTrace::Data::make_surface_from_json(jcone);
    auto cone_cast = dynamic_cast<SolTrace::Data::Cone*>(cone_ptr.get());
    ASSERT_TRUE(cone_cast != nullptr);
    EXPECT_DOUBLE_EQ(1, cone_cast->half_angle);

    // CYLINDER
    json jcyl;
    jcyl["surface_type"] = SolTrace::Data::SurfaceTypeMap.at(SolTrace::Data::CYLINDER);
    jcyl["radius"] = 2;
    EXPECT_NO_THROW(SolTrace::Data::make_surface_from_json(jcyl));
    auto cyl_ptr = SolTrace::Data::make_surface_from_json(jcyl);
    auto cyl_cast = dynamic_cast<SolTrace::Data::Cylinder*>(cyl_ptr.get());
    ASSERT_TRUE(cyl_cast != nullptr);
    EXPECT_DOUBLE_EQ(2, cyl_cast->radius);

    // FLAT
    json jflat;
    jflat["surface_type"] = SolTrace::Data::SurfaceTypeMap.at(SolTrace::Data::FLAT);
    EXPECT_NO_THROW(SolTrace::Data::make_surface_from_json(jflat));
    auto flat_ptr = SolTrace::Data::make_surface_from_json(jflat);
    auto flat_cast = dynamic_cast<SolTrace::Data::Flat*>(flat_ptr.get());
    ASSERT_TRUE(flat_cast != nullptr);

    // PARABOLA
    json jpara;
    jpara["surface_type"] = SolTrace::Data::SurfaceTypeMap.at(SolTrace::Data::PARABOLA);
    jpara["focal_length_x"] = 3;
    jpara["focal_length_y"] = 4;
    EXPECT_NO_THROW(SolTrace::Data::make_surface_from_json(jpara));
    auto para_ptr = SolTrace::Data::make_surface_from_json(jpara);
    auto para_cast = dynamic_cast<SolTrace::Data::Parabola*>(para_ptr.get());
    ASSERT_TRUE(para_cast != nullptr);
    EXPECT_DOUBLE_EQ(3, para_cast->focal_length_x);
    EXPECT_DOUBLE_EQ(4, para_cast->focal_length_y);

    // SPHERE
    json jsphere;
    jsphere["surface_type"] = SolTrace::Data::SurfaceTypeMap.at(SolTrace::Data::SPHERE);
    jsphere["vertex_curv"] = 5;
    EXPECT_NO_THROW(SolTrace::Data::make_surface_from_json(jsphere));
    auto sphere_ptr = SolTrace::Data::make_surface_from_json(jsphere);
    auto sphere_cast = dynamic_cast<SolTrace::Data::Sphere*>(sphere_ptr.get());
    ASSERT_TRUE(sphere_cast != nullptr);
    EXPECT_DOUBLE_EQ(5, sphere_cast->vertex_curv);

    // No surface type
    json jmissing;
    jmissing["radius"] = 0;
    EXPECT_THROW(SolTrace::Data::make_surface_from_json(jmissing), std::invalid_argument);

    // Unknown surface type
    json junknown;
    junknown["surface_type"] = SolTrace::Data::SurfaceTypeMap.at(SolTrace::Data::SURFACE_UNKNOWN);
    EXPECT_THROW(SolTrace::Data::make_surface_from_json(junknown), std::invalid_argument);

    int x = 0;
}

TEST(io_json, stage_read_fail)
{
    using json = nlohmann::ordered_json;
    json jstage;
    get_default_element_base(jstage);
    jstage["is_stage"] = false;
    jstage["elements"] = json::object(); // Empty node

    // Try to make stage
    EXPECT_THROW(make_stage(jstage), std::invalid_argument);

}
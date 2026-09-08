#include <gtest/gtest.h>

#include <algorithm>

#include <aperture.hpp>
#include <surface.hpp>
#include <constants.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result_export.hpp>
#include <embree_runner.hpp>

#include <json_helpers.hpp>

#include "common.hpp"

using SolTrace::EmbreeRunner::EmbreeRunner;
using SolTrace::NativeRunner::TRayData;
using SolTrace::NativeRunner::TSystem;
using SolTrace::Result::RayEvent;
using SolTrace::Runner::RunnerStatus;

class grouped_results_EmbreeRunner_helper {
    public:
        static void append(EmbreeRunner &runner, uint_fast64_t raynum, uint_fast64_t element_id, RayEvent hit_type, glm::dvec3 &hit_point) 
        {
            glm::dvec3 cos = {0.0, 0.0, 0.0};
            runner.tsys.RayData.Append(0, hit_point, cos, element_id, 0, raynum, hit_type);
            if (hit_type == RayEvent::CREATE) ++runner.tsys.SunRayCount;
        }
        static void setUp(EmbreeRunner &runner) 
        {
            runner.tsys.RayData.SetUp(1, 20000);
            runner.tsys.Sun.MaxXSun = 1;
            runner.tsys.Sun.MaxYSun = 1;
        }
};

void create_sun_record(EmbreeRunner &runner, uint_fast64_t raynum)
{
    glm::dvec3 pos = {0.0, 0.0, 1000.0};
    grouped_results_EmbreeRunner_helper::append(runner, raynum, -1, RayEvent::CREATE, pos);
}

void create_hit_record(EmbreeRunner &runner, uint_fast64_t raynum, uint_fast64_t element_id, RayEvent hit_type, glm::dvec3 &hit_point)
{
    grouped_results_EmbreeRunner_helper::append(runner, raynum, element_id, hit_type, hit_point);
}


TEST(grouped_results, counts_test) {
    namespace fs = std::filesystem;

    // Build paths
    const fs::path project_root(PROJECT_DIR);
    const std::string input_str = project_root.string() + "/field_test.json";
    const std::string output_str = project_root.string() + "/field_out_embree.json";

    SimulationData sd;
    ASSERT_NO_THROW(sd.import_json_file(input_str));

    SimulationParameters &params = sd.get_simulation_parameters();
    params.number_of_rays = static_cast<uint_fast64_t>(100);
    params.max_number_of_rays = params.number_of_rays * 100;
    params.seed = 608;
    
    EmbreeRunner runner;
    runner.set_number_of_threads(1);
    
    RunnerStatus sts = runner.initialize();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS) << "runner.initialize() failed";
    
    runner.disable_stages();
    sts = runner.setup_simulation(&sd);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS) << "runner.setup_simulation() failed";
    grouped_results_EmbreeRunner_helper::setUp(runner);
    
    // Check groups
    std::vector<std::set<uint_fast64_t>> groups = runner.get_groups();
    EXPECT_EQ(groups.size(), 5);
    ASSERT_EQ(runner.get_num_groups(), 5) << "Number of groups in system does not match expected";
    ASSERT_EQ(runner.get_group(26), -1) << "Element 26 should be ungrouped";
    ASSERT_EQ(runner.get_group(27), 0) << "Element 27 should be in group 0";
    ASSERT_EQ(runner.get_group(127), 4) << "Element 127 should be in group 4";

    // conjure up some hit records to check the counting algorithm
    // heliostat elements span 2 - > 126
    // receiver element is 127
    // number of reflection events per heliostat element is element id - 1 
    // -> expect triangular number of facets + group number * number of facets ^ 2 reflection events in group
    // for each heliostat increment the number of absorptions by 1 ie ungrouped: 1 per facet, group 0: 2 per facet, etc
    // number of absorptions on receiver per heliostat is (group number + 2) * number of facets
    // -> expect 50 absorptions from group 0, 75 from group 1, 100 from group 2, 125 from group 3, and 150 from group 4
    // therefore total absorptions 50 + 75 + 100 + 125 + 25 (ungrouped) = 375
    uint_fast64_t rec_id = 127;
    SolTrace::Data::element_ptr rec = sd.get_element(rec_id);
    glm::dvec3 rec_origin = rec->get_origin_global();
    uint_fast64_t raynum = 0;
    
    for (uint_fast64_t i = 2; i < rec_id; ++i)
    {
        uint_fast64_t facet_mod = (i - 1) % 25;
        SolTrace::Data::element_ptr el = sd.get_element(i);
        glm::dvec3 origin = el->get_origin_global();

        for (uint_fast64_t j = 1; j < i; ++j)
        {
            // add create record
            create_sun_record(runner, raynum);

            // add reflection record
            create_hit_record(
                runner,
                raynum,
                i - 1, // simulation id
                RayEvent::REFLECT,
                origin
            );
            
            if (facet_mod == j % 25)
            {
                create_hit_record(
                    runner,
                    raynum,
                    rec_id - 1, // simulation id
                    RayEvent::ABSORB,
                    rec_origin
                );
            }
        }
    }

    SimulationResult result;
    sts = runner.report_simulation(&result, SolTrace::Runner::RunnerStatistics::GROUPED_COUNTS);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS) << "runner.report_simulation() failed";

    std::vector<GroupResult> grouped_results = result.get_grouped_results();

    uint_fast64_t tri_num = 25 * 26 / 2; // triangular number of facets
    ASSERT_EQ(grouped_results[0].reflect_count, tri_num + 1 * 25 * 25);
    ASSERT_EQ(grouped_results[1].reflect_count, tri_num + 2 * 25 * 25);
    ASSERT_EQ(grouped_results[2].reflect_count, tri_num + 3 * 25 * 25);
    ASSERT_EQ(grouped_results[3].reflect_count, tri_num + 4 * 25 * 25);
    
    ASSERT_EQ(grouped_results[4].reflect_count,  0);
    ASSERT_EQ(grouped_results[4].absorb_count, 375);
    ASSERT_EQ(grouped_results[4].absorb_previous_group[0],  50);
    ASSERT_EQ(grouped_results[4].absorb_previous_group[1],  75);
    ASSERT_EQ(grouped_results[4].absorb_previous_group[2], 100);
    ASSERT_EQ(grouped_results[4].absorb_previous_group[3], 125);

    // should look like unit-tests\simulation_data\grouped_elements_io\field_out_reference.json
    ASSERT_NO_THROW(result.write_group_json_file(output_str));
}
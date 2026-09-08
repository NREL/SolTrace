#include <gtest/gtest.h>

#include <algorithm>

#include <aperture.hpp>
#include <surface.hpp>
#include <constants.hpp>
#include <optix_runner.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result_export.hpp>
#include <json_helpers.hpp>

#include "common.hpp"

OptixCSP::HitRecord create_sun_record()
{
    OptixCSP::HitRecord hr;
    hr.element_id = -1;
    hr.hit_type = static_cast<uint8_t>(OptixCSP::HitType::HIT_CREATE);
    hr.hit_point = float4();
    hr.hit_point.w = 0.0f;
    hr.hit_point.x = 0.0f;
    hr.hit_point.y = 0.0f;
    hr.hit_point.z = 1000.0f;
    return hr;
}

OptixCSP::HitRecord create_hit_record(int32_t element_id, uint8_t hit_type, const glm::dvec3 &hit_point)
{
    OptixCSP::HitRecord hr;
    hr.element_id = element_id;
    hr.hit_type = hit_type;
    hr.hit_point = float4();
    hr.hit_point.w = 0.0f;
    hr.hit_point.x = static_cast<float>(hit_point.x);
    hr.hit_point.y = static_cast<float>(hit_point.y);
    hr.hit_point.z = static_cast<float>(hit_point.z);
    return hr;
}

class grouped_results_SolTraceSystem_helper 
{
    public:
        static void set_hit_records(OptixCSP::SolTraceSystem *sys, const std::vector<OptixCSP::HitRecord> &hit_records) 
        {
            sys->m_hit_records = hit_records;
        }
        static void set_sun(OptixCSP::SolTraceSystem* sys)
        {
            sys->m_n_sun_rays = 1;
        }
};

TEST(grouped_results, counts_test) 
{
    using SolTrace::Runner::RunnerStatus;
    namespace fs = std::filesystem;

    // Build paths
    const fs::path project_root(PROJECT_DIR);
    const std::string input_str = project_root.string() + "/field_test.json";
    const std::string output_str = project_root.string() + "/field_out.json";

    SimulationData sd;
    ASSERT_NO_THROW(sd.import_json_file(input_str));

    SimulationParameters &params = sd.get_simulation_parameters();
    params.number_of_rays = static_cast<uint_fast64_t>(100);
    params.max_number_of_rays = params.number_of_rays * 100;
    params.seed = 608;
    
    OptixRunner runner;
    
    RunnerStatus sts = runner.initialize();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS) << "runner.initialize() failed";
    
    sts = runner.setup_simulation(&sd);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS) << "runner.setup_simulation() failed";
    
    // Check groups
    std::vector<std::set<int32_t>> groups = runner.get_groups();
    EXPECT_EQ(groups.size(), 5);
    ASSERT_EQ(runner.get_num_groups(), 5) << "Number of groups in system does not match expected";
    ASSERT_EQ(runner.get_group(26), -1) << "Element 26 should be ungrouped";
    ASSERT_EQ(runner.get_group(27), 0) << "Element 27 should be in group 0";
    ASSERT_EQ(runner.get_group(127), 4) << "Element 127 should be in group 4";

    // conjure up some hit records to check the counting algorithm
    OptixCSP::SolTraceSystem *sys = runner.get_optix_system();
    std::vector<OptixCSP::HitRecord> hit_records;

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
    
    for (uint_fast64_t i = 2; i < rec_id; ++i)
    {
        uint_fast64_t facet_mod = (i - 1) % 25;
        for (uint_fast64_t j = 1; j < i; ++j)
        {
            // add create record
            hit_records.push_back(create_sun_record());

            // add reflection record
            SolTrace::Data::element_ptr el = sd.get_element(i);
            glm::dvec3 origin = el->get_origin_global();

            hit_records.push_back(create_hit_record(
                static_cast<int32_t>(i),
			    static_cast<uint8_t>(OptixCSP::HitType::HIT_REFLECT),
			    origin
            ));
            
            if (facet_mod == j % 25)
            {
                hit_records.push_back(create_hit_record(
                    static_cast<int32_t>(rec_id),
                    static_cast<uint8_t>(OptixCSP::HitType::HIT_ABSORB),
                    rec_origin
                ));
            }
        }
    }

    grouped_results_SolTraceSystem_helper::set_hit_records(sys, hit_records);
    grouped_results_SolTraceSystem_helper::set_sun(sys);

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

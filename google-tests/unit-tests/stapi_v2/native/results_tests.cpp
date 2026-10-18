#include <gtest/gtest.h>

#include "across_builds.hpp"


TEST(results_tests, results_write_csv)
{
    SETUP_TEST_CXT();

    LOAD_TEST_JSON();

    code = call_stapi_v2_write_results_csv(pcxt, st_runner_type_t::NATIVE, "temp_char.csv");
    EXPECT_EQ(code, st_return_code::SUCCESS);
    
    CLEANUP_TEST_CXT();
}

TEST(results_tests, results_write_group_json)
{
    SETUP_TEST_CXT();

    LOAD_TEST_JSON();

    code = call_stapi_v2_write_group_results_json(pcxt, st_runner_type_t::NATIVE, "temp_char.csv");
    EXPECT_EQ(code, st_return_code::SUCCESS);

    CLEANUP_TEST_CXT();
}

TEST(results_tests, results_locations)
{
    SETUP_TEST_CXT();

    LOAD_TEST_JSON();

    code = call_stapi_v2_locations(pcxt, st_runner_type_t::NATIVE);
    EXPECT_EQ(code, st_return_code::SUCCESS);
    
    CLEANUP_TEST_CXT();
}

TEST(results_tests, results_cosines)
{
    SETUP_TEST_CXT();

    LOAD_TEST_JSON();

    code = call_stapi_v2_cosines(pcxt, st_runner_type_t::NATIVE);
    EXPECT_EQ(code, st_return_code::SUCCESS);
    
    CLEANUP_TEST_CXT();
}

TEST(results_tests, results_elementmap)
{
    SETUP_TEST_CXT();

    LOAD_TEST_JSON();

    code = call_stapi_v2_elementmap(pcxt, st_runner_type_t::NATIVE);
    EXPECT_EQ(code, st_return_code::SUCCESS);
    
    CLEANUP_TEST_CXT();
}

TEST(results_tests, results_stagemap)
{
    SETUP_TEST_CXT();

    LOAD_TEST_JSON();

    code = call_stapi_v2_stagemap(pcxt, st_runner_type_t::NATIVE);
    EXPECT_EQ(code, st_return_code::SUCCESS);
    
    CLEANUP_TEST_CXT();
}

TEST(results_tests, results_raynumbers)
{
    SETUP_TEST_CXT();

    LOAD_TEST_JSON();

    code = call_stapi_v2_raynumbers(pcxt, st_runner_type_t::NATIVE);
    EXPECT_EQ(code, st_return_code::SUCCESS);
    
    CLEANUP_TEST_CXT();
}

TEST(results_tests, results_sun_stats)
{
    SETUP_TEST_CXT();

    LOAD_TEST_JSON();

    code = call_stapi_v2_sun_stats(pcxt, st_runner_type_t::NATIVE);
    EXPECT_EQ(code, st_return_code::SUCCESS);
    
    CLEANUP_TEST_CXT();
}

TEST(results_tests, results_get_results_data)
{
    SETUP_TEST_CXT();

    LOAD_TEST_JSON();

    code = call_stapi_v2_get_results_data(pcxt, st_runner_type_t::NATIVE);
    EXPECT_EQ(code, st_return_code::SUCCESS);
    
    CLEANUP_TEST_CXT();
}

#include <gtest/gtest.h>

#include "across_builds.hpp"

TEST(runner_tests, runners_installed)
{
    SETUP_TEST_CXT();

    uint8_t installed;
    code = st_get_installed_runners(&installed);
    EXPECT_EQ(code, st_return_code::SUCCESS);
    EXPECT_EQ(installed, 1 << st_runner_type_t::NATIVE);

    CLEANUP_TEST_CXT();
}

TEST(runner_tests, runners_is_installed)
{
    SETUP_TEST_CXT();

    bool installed;
    code = st_is_runner_installed(st_runner_type_t::NATIVE, &installed);
    EXPECT_EQ(code, st_return_code::SUCCESS);
    EXPECT_EQ(installed, true);
    
    code = st_is_runner_installed(st_runner_type_t::EMBREE, &installed);
    EXPECT_EQ(code, st_return_code::SUCCESS);
    EXPECT_EQ(installed, false);
    
    code = st_is_runner_installed(st_runner_type_t::OPTIX, &installed);
    EXPECT_EQ(code, st_return_code::SUCCESS);
    EXPECT_EQ(installed, false);

    CLEANUP_TEST_CXT();
}

TEST(runner_tests, runner_setup)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_sim_setup(pcxt);
    EXPECT_EQ(code, 
        st_return_code::WARNING_FELLBACK_FROM_EMBREE 
        + st_return_code::WARNING_FELLBACK_FROM_OPTIX
        + st_return_code::RUNNER_NUMBER_THREADS_SEEDS_MISMATCH_FAILURE);

    CLEANUP_TEST_CXT();
}

TEST(runner_tests, runner_run_native)
{
    SETUP_TEST_CXT();

    LOAD_TEST_JSON();

    code = call_stapi_v2_sim_run_v2(pcxt, st_runner_type_t::NATIVE);
    EXPECT_EQ(code, st_return_code::SUCCESS);

    CLEANUP_TEST_CXT();
}

TEST(runner_tests, runner_not_found)
{
    SETUP_TEST_CXT();

    code = st_sim_run_v2(pcxt);
    EXPECT_EQ(code, st_return_code::RUNNER_NOT_FOUND);

    CLEANUP_TEST_CXT();
}

TEST(runner_tests, runner_report)
{
    SETUP_TEST_CXT();

    LOAD_TEST_JSON();

    code = call_stapi_v2_sim_report(pcxt, st_runner_type_t::NATIVE);
    EXPECT_EQ(code, st_return_code::RUNNER_NOT_READY_TO_REPORT
                    + 2 * st_return_code::INVALID_ARGUMENTS);

    // code = call_stapi_v2_sim_report(pcxt, st_runner_type_t::EMBREE);
    // EXPECT_EQ(code, st_return_code::RUNNER_NOT_READY_TO_REPORT
    //                 + 2 * st_return_code::INVALID_ARGUMENTS);
    
    // code = call_stapi_v2_sim_report(pcxt, st_runner_type_t::OPTIX);
    // EXPECT_EQ(code, st_return_code::RUNNER_NOT_READY_TO_REPORT
    //                 + 2 * st_return_code::INVALID_ARGUMENTS);

    CLEANUP_TEST_CXT();
}

// TEST(runner_tests, runner_report_native)
// {
//     SETUP_TEST_CXT();

//     LOAD_TEST_JSON();

//     code = call_stapi_v2_sim_run_v2(pcxt, st_runner_type_t::NATIVE);
//     EXPECT_EQ(code, st_return_code::SUCCESS);
    
//     code = st_sim_report(pcxt, 0);
//     EXPECT_EQ(code, st_return_code::SUCCESS);

//     CLEANUP_TEST_CXT();
// }
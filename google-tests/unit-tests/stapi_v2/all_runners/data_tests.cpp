#include <gtest/gtest.h>

#include "across_builds.hpp"

TEST(data_tests, data_set_simulation_params)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_set_simulation_parameters(pcxt);
    EXPECT_EQ(code, st_return_code::SUCCESS);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_params)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_sim_params(pcxt);
    EXPECT_EQ(code, st_return_code::SUCCESS);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_rays)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_sim_rays(pcxt);
    EXPECT_EQ(code, st_return_code::SUCCESS);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_power_tower)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_sim_power_tower(pcxt);
    EXPECT_EQ(code, st_return_code::SUCCESS);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_errors)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_sim_errors(pcxt);
    EXPECT_EQ(code, st_return_code::SUCCESS);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_location)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_sim_location(pcxt);
    EXPECT_EQ(code, st_return_code::SUCCESS);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_tolerance)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_sim_tolerance(pcxt);
    EXPECT_EQ(code, st_return_code::SUCCESS);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_get_simulation_params)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_get_simulation_parameters(pcxt);
    EXPECT_EQ(code, st_return_code::SUCCESS);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_add_optics)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_add_optics(pcxt);
    EXPECT_EQ(code, 2 * st_return_code::INVALID_ARGUMENTS);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_get_optic)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_get_optic(pcxt);
    EXPECT_EQ(code, st_return_code::DATA_VALUE_NOT_FOUND);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_remove_optics)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_remove_optics(pcxt);
    EXPECT_EQ(code, st_return_code::SUCCESS);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_add_elements)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_add_elements(pcxt);
    EXPECT_EQ(code, st_return_code::DATA_VALUE_NOT_FOUND
                    + 4 * st_return_code::INVALID_ARGUMENTS);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_get_element)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_get_element(pcxt);
    EXPECT_EQ(code, st_return_code::DATA_VALUE_NOT_FOUND);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_remove_elements)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_remove_elements(pcxt);
    EXPECT_EQ(code, st_return_code::WARNING_NOT_FOUND);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_toggle_element)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_toggle_element(pcxt);
    EXPECT_EQ(code, 2 * st_return_code::WARNING_NOT_FOUND);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_element_xyz)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_element_xyz(pcxt);
    EXPECT_EQ(code, st_return_code::WARNING_NOT_FOUND);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_element_aim)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_element_aim(pcxt);
    EXPECT_EQ(code, st_return_code::WARNING_NOT_FOUND);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_element_zrot)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_element_zrot(pcxt);
    EXPECT_EQ(code, st_return_code::WARNING_NOT_FOUND);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_element_aperture)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_element_aperture(pcxt);
    EXPECT_EQ(code, 2 * st_return_code::INVALID_ARGUMENTS);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_element_surface)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_element_surface(pcxt);
    EXPECT_EQ(code, 2 * st_return_code::INVALID_ARGUMENTS);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_element_geometry)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_element_geometry(pcxt);
    EXPECT_EQ(code, 4 * st_return_code::INVALID_ARGUMENTS);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_element_optic)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_element_optic(pcxt);
    EXPECT_EQ(code, st_return_code::DATA_VALUE_NOT_FOUND);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_element_group)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_element_group(pcxt);
    EXPECT_EQ(code, 2 * st_return_code::WARNING_GROUP_IGNORED);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_add_sun)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_add_sun(pcxt);
    EXPECT_EQ(code, 3 * st_return_code::WARNING_SUN_SHAPE_IGNORED
                    + st_return_code::EXCEPTION);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_get_sun)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_get_sun(pcxt);
    EXPECT_EQ(code, st_return_code::DATA_VALUE_NOT_FOUND +
                    2 * st_return_code::INVALID_ARGUMENTS);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_sun_shape)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_sun_shape(pcxt);
    EXPECT_EQ(code, 3 * st_return_code::WARNING_SUN_SHAPE_IGNORED);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_sun_xyz)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_sun_xyz(pcxt);
    EXPECT_EQ(code, st_return_code::SUCCESS);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_sun_userdata)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_sun_userdata(pcxt);
    EXPECT_EQ(code, st_return_code::EXCEPTION);

    CLEANUP_TEST_CXT();
}

TEST(data_tests, data_sun_calculator)
{
    SETUP_TEST_CXT();

    code = call_stapi_v2_solar_calculator(pcxt);
    EXPECT_EQ(code, st_return_code::INVALID_ARGUMENTS);

    CLEANUP_TEST_CXT();
}
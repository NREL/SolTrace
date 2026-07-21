#include "isolated_heliostat_test_template.hpp"

#include <embree_runner.hpp>

using EmbreeRunnerType = SolTrace::EmbreeRunner::EmbreeRunner;

using IsolatedHeliostatSimulationEmbree = IsolatedHeliostatSimulation<EmbreeRunnerType>;

static const int N_threads = static_cast<int>(std::max(1u, std::min(std::thread::hardware_concurrency(), 10u)));

//DEFAULT IS FACET FOCUS TO SLANT RANGE

//task 4a: r f 6x5, canting by band, facet focusing by band

/*TEST_F(IsolatedHeliostatSimulationEmbree, multiFacet1332_BlockingShading4a)
{
    this->runner.set_number_of_threads(N_threads);

    // Centerline aimpoints
    set_flat_facets();
    glm::dvec3 origin = {0,0,171.035};
    set_rec_origin(origin);
    std::vector<int> active {1332};
    std::vector<int> blocking {1266, 1304, 1306};
    create_active_heliostats(active);
    create_blocking_heliostats(blocking);

    assign_canted_banded(false);
    assign_canted_banded(true);

    setup_simData();
    update_from_hour("12");
    simulate_check_outputs("4a", "1");
    //simulate_check_outputs("7b", "1");
}*/

//task 4b: r p 6x5, canting by band, facet focusing by band
TEST_F(IsolatedHeliostatSimulationEmbree, multiFacet8993_BlockingShading4b)
{
    this->runner.set_number_of_threads(N_threads);

    // Centerline aimpoints
    set_slope_error(0.0);
    glm::dvec3 origin = {0,0,171.035};
    set_rec_origin(origin);
    std::vector<int> active {8993};
    std::vector<int> blocking {9100, 9102, 9208};
    create_active_heliostats(active);
    create_blocking_heliostats(blocking);

    assign_canted_banded(false);
    assign_canted_banded(true);
    
    assign_focal_lengths_banded(false);
    assign_focal_lengths_banded(true);

    setup_simData();
    update_from_hour("12");
    simulate_check_outputs("4b", "1");
    //simulate_check_outputs("7b", "1");

}

//task 4c: r f 6x5, canting by band, facet focusing by band

TEST_F(IsolatedHeliostatSimulationEmbree, multiFacet5473_BlockingShading4c)
{
    this->runner.set_number_of_threads(N_threads);

    // Centerline aimpoints
    glm::dvec3 origin = {0,0,171.035};
    set_rec_origin(origin);
    set_flat_facets();
    std::vector<int> active {5473};
    std::vector<int> blocking {5573};
    create_active_heliostats(active);
    create_blocking_heliostats(blocking);

    assign_canted_banded(false);
    assign_canted_banded(true);

    setup_simData();
    update_from_hour("12");
    simulate_check_outputs("4c", "1");
    save_flux_map_to_file("embree_test_4c.csv");
    //simulate_check_outputs("7c", "1");

}

//task 4d: r p 1, facet focusing by slant range
TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet8993_4dLongerAimpoint8)
{
    this->runner.set_number_of_threads(N_threads);
    
    // Centerline aimpoint
    set_slope_error(0.0);
    std::vector<int> active {8993};
    create_active_heliostats(active);
    setup_simData();
    
    simulate_check_outputs("4d", "1", "8");
    save_flux_map_to_file("embree_test_4d_8.csv");

    // update_from_hour("8");
    // simulate_check_outputs("7b", "1");
    
    
}

TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet8993_4dLongerAimpoint12)
{
    this->runner.set_number_of_threads(N_threads);
    
    // Centerline aimpoint
    set_slope_error(0.0);
    std::vector<int> active {8993};
    create_active_heliostats(active);
    set_no_sunShape();
    setup_simData();
    
    simulate_check_outputs("4d", "1", "12");
    // update_from_hour("12");
    // simulate_check_outputs("7b", "1");
}

/*TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet8993_4dLongerAimpoint_blocking)
{
    this->runner.set_number_of_threads(N_threads);
    
    // Centerline aimpoints
    std::vector<int> active {8993};
    std::vector<int> blocking {9100, 9102, 9208};
    create_active_heliostats(active);
    create_blocking_heliostats(blocking);
    setup_simData();
    
    simulate_check_outputs("4d", "1", "8");
    simulate_check_outputs("4d", "1", "12");
}*/

/*TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet8993_4dDoctoredAimpoint)
{
    this->runner.set_number_of_threads(N_threads);
    
    // Centerline aimpoints
    set_helio_dim(10,10);
    std::vector<int> active {8993};
    create_active_heliostats(active);
    setup_simData();
    
    simulate_check_outputs("4d", "1", "8");
    simulate_check_outputs("4d", "1", "12");
}*/

//task 4e: r p 1, facet focusing by slant range
TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet8993_TargetCoordSystemE4e)
{
    this->runner.set_number_of_threads(N_threads);

    // Centerline aimpoints
    glm::dvec3 shift = {-1,0,0};
    shift_aimpoint(shift);
    set_slope_error(0.0);
    std::vector<int> active {8993};
    create_active_heliostats(active);
    
    //  TODO: function that shifts aimpoint, no sun shape, no slope error

    set_no_sunShape();
    setup_simData();

    // update_from_hour("8");
    // simulate_check_outputs("7b", "1");
    // update_from_hour("12");
    // simulate_check_outputs("7b", "1");
    
    simulate_check_outputs("4e", "1", "8");
    save_flux_map_to_file("embree_test_4e_8.csv");
    simulate_check_outputs("4e", "1", "12");
    save_flux_map_to_file("embree_test_4e_12.csv");
}

//task 4f: r p 1, facet focusing by slant range
TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet8993_TargetCoordSystemW4f)
{
    this->runner.set_number_of_threads(N_threads);
    // Centerline aimpoints
    glm::dvec3 shift = {1,0,0};
    shift_aimpoint(shift);
    set_slope_error(0.0);
    std::vector<int> active {8993};
    create_active_heliostats(active);
    
    //  TODO: function that shifts aimpoint, no sun shape, no slope error
    
    set_no_sunShape();
    setup_simData();

    // update_from_hour("8");
    // simulate_check_outputs("7b", "1");
    // update_from_hour("12");
    // simulate_check_outputs("7b", "1");
    
    simulate_check_outputs("4f", "1", "8");
    save_flux_map_to_file("embree_test_4f_8.csv");
    simulate_check_outputs("4f", "1", "12");

}

//task 4g: r p 1, facet focusing by slant range
TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet8993_TargetCoordSystemU4g)
{
    this->runner.set_number_of_threads(N_threads);
    // Centerline aimpoints
    glm::dvec3 shift = {0,0,1};
    shift_aimpoint(shift);
    set_slope_error(0.0);
    std::vector<int> active {8993};
    create_active_heliostats(active);
    
    //  TODO: function that shifts aimpoint, no sun shape, no slope error
    
    set_no_sunShape();
    setup_simData();

    // update_from_hour("8");
    // simulate_check_outputs("7b", "1");
    // update_from_hour("12");
    // simulate_check_outputs("7b", "1");

    simulate_check_outputs("4g", "1", "8");
    save_flux_map_to_file("embree_test_4g_8.csv");
    simulate_check_outputs("4g", "1", "12");

}

//task 4h: r p 1, facet focusing by slant rangeTEST_F(IsolatedHeliostatSimulationEmbree, singleFacet8993_TargetCoordSystemD)
TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet8993_TargetCoordSystemD4h)
{
    this->runner.set_number_of_threads(N_threads);
    // Centerline aimpoints
    glm::dvec3 shift = {0,0,-1};
    shift_aimpoint(shift);
    set_slope_error(0.0);
    std::vector<int> active {8993};
    create_active_heliostats(active);
    
    //  TODO: function that shifts aimpoint, no sun shape, no slope error
    
    set_no_sunShape();
    setup_simData();

    // update_from_hour("8");
    // simulate_check_outputs("7b", "1");
    // update_from_hour("12");
    // simulate_check_outputs("7b", "1");

    simulate_check_outputs("4h", "1", "8");
    save_flux_map_to_file("embree_test_4h_8.csv");
    simulate_check_outputs("4h", "1", "12");

}

//task 4i: r f 1, facet focusing by slant range
TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet8993_4i_8)
{
    this->runner.set_number_of_threads(N_threads);
    set_slope_error(0.0);
    set_flat_facets();
    std::vector<int> active {8993};
    create_active_heliostats(active);
    
    //  TODO: function that shifts aimpoint, no sun shape, no slope error, FLAT??????
    set_no_sunShape();
    setup_simData();
    
    // update_from_hour("8");
    // simulate_check_outputs("7b", "1");

    simulate_check_outputs("4i", "1", "8");
    save_flux_map_to_file("embree_test_4i_8.csv");

}

TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet8993_4i_12)
{
    this->runner.set_number_of_threads(N_threads);
    set_slope_error(0.0);
    set_flat_facets();
    std::vector<int> active {8993};
    create_active_heliostats(active);
    
    //  TODO: function that shifts aimpoint, no sun shape, no slope error, FLAT??????
    set_no_sunShape();
    setup_simData();

    // update_from_hour("12");
    // simulate_check_outputs("7b", "1");

    simulate_check_outputs("4i", "1", "12");

}
//task 5a: r p 6x5, canting by band, facet focusing by band
TEST_F(IsolatedHeliostatSimulationEmbree, multiFacet1332_BlockingShading5a)
{
    this->runner.set_number_of_threads(N_threads);

    // Centerline aimpoints
    glm::dvec3 origin = {0,0,171.035};
    set_rec_origin(origin);
    std::vector<int> active {1332};
    std::vector<int> blocking {1266, 1304, 1306, 1334, 5321};
    create_active_heliostats(active);
    create_blocking_heliostats(blocking);

    assign_canted_banded(false);
    assign_canted_banded(true);
    
    assign_focal_lengths_banded(false);
    assign_focal_lengths_banded(true);

    setup_simData();
    
    //update_simulation_geometry(74.95, 26.26);
    update_from_hour("8");
    simulate_check_outputs("5a", "1");
    save_flux_map_to_file("embree_test_5a.csv");

    //simulate_check_outputs("7b", "1");

    

}

//task 5b: NO FILE r  1, 

//task 6a: r f 6x5, canting by band, flat
TEST_F(IsolatedHeliostatSimulationEmbree, multiFacet1332_CantingAccuracy6a)
{
    this->runner.set_number_of_threads(N_threads);

    // Centerline aimpoints
    set_flat_facets();
    glm::dvec3 origin = {0,0,171.035};
    set_rec_origin(origin);
    std::vector<int> active {1332};
    create_active_heliostats(active);

    assign_canted_banded(true);

    setup_simData();
    update_from_hour("12");
    simulate_check_outputs("6a", "1");
    //simulate_check_outputs("7b", "1");

}

//task 6b: NO FILE r  1, facet focusing by slant range

//task 6c: NO FILE r  1

//task 7a: r p 6x5, canting by slant range, facet focusing by slant range
TEST_F(IsolatedHeliostatSimulationEmbree, multiFacet5473_CantingFocusingAccuracy7a) //0 slope error
{
    this->runner.set_number_of_threads(N_threads);
    // Centerline aimpoints
    glm::dvec3 shift = {0,0,0.825};//due to truncation error in legacy
    shift_aimpoint(shift);

    set_slope_error(0.0);
    glm::dvec3 origin = {0,0,171.035};
    set_rec_origin(origin);
    std::vector<int> active {5473};
    create_active_heliostats(active);

    assign_canted_slant(true);

    setup_simData();
    update_from_hour("12");
    simulate_check_outputs("7a", "1");
    //simulate_check_outputs("5a", "1");
    save_flux_map_to_file("embree_test_7a.csv");

}

//task 7b: r p 1, facet focusing by slant range
TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet5473_7b) //0 slope error
{
    this->runner.set_number_of_threads(N_threads);
    // Centerline aimpoints
    glm::dvec3 shift = {0,0,0.825};//due to truncation error in legacy
    shift_aimpoint(shift);

    set_slope_error(0.0);
    glm::dvec3 origin = {0,0,171.035};
    set_rec_origin(origin);
    std::vector<int> active {5473};
    create_active_heliostats(active);
    setup_simData();
    update_from_hour("12");
    simulate_check_outputs("7b", "1");
    //simulate_check_outputs("5a", "1");
    save_flux_map_to_file("embree_test_7b_12.csv");

}

TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet5473_7bLongerAimpoint_8) //0 slope error
{
    this->runner.set_number_of_threads(N_threads);
    // Centerline aimpoints
    set_slope_error(0.0);
    std::vector<int> active {5473};
    create_active_heliostats(active);
    setup_simData();
    simulate_check_outputs("7b", "1", "8");
    //simulate_check_outputs("5a", "1");

}

/*TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet5473_7bDebug) //0 slope error
{
    this->runner.set_number_of_threads(N_threads);
    // Centerline aimpoints
    set_slope_error(0.0);
    glm::dvec3 origin = {-0.2273,7.7217,171.035};
    set_rec_origin(origin);
    std::vector<int> active {5473};
    create_active_heliostats(active);
    setup_simData();
    update_from_hour("12");
    simulate_check_outputs("7b", "1");

}*/

/*TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet5473_7bDebugFartherAimpoint) //0 slope error
{
    this->runner.set_number_of_threads(N_threads);
    // Centerline aimpoints
    set_slope_error(0.0);
    glm::dvec3 origin = {0,0,171.035};
    set_rec_origin(origin);
    std::vector<int> active {5473};
    create_active_heliostats(active);
    setup_simData();
    update_from_hour("12");
    simulate_check_outputs("7b", "1");

}*/

//task 7c: r p 1, facet focusing by slant range
TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet5473_BlockingShading7c) //0 slope error
{
    this->runner.set_number_of_threads(N_threads);
    // Centerline aimpoints
    glm::dvec3 shift = {0,0,0.825};//due to truncation error in legacy
    shift_aimpoint(shift);

    set_slope_error(0.0);
    glm::dvec3 origin = {0,0,171.035};
    set_rec_origin(origin);
    std::vector<int> active {5473};
    std::vector<int> blocking {5573};
    create_active_heliostats(active);
    create_blocking_heliostats(blocking);

    setup_simData();
    update_from_hour("12");
    simulate_check_outputs("7c", "1");
    //simulate_check_outputs("5a", "1");
    save_flux_map_to_file("embree_test_7c.csv");

}

//task 7d: NO FILE r p 1 (guess)
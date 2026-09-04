#include "isolated_heliostat_test_template.hpp"

#include <embree_runner.hpp>

using EmbreeRunnerType = SolTrace::EmbreeRunner::EmbreeRunner;

using IsolatedHeliostatSimulationEmbree = IsolatedHeliostatSimulation<EmbreeRunnerType>;

static const int N_threads = static_cast<int>(std::max(1u, std::min(std::thread::hardware_concurrency(), 10u)));

/*
DEFAULT IS FACET FOCUS TO SLANT RANGE
Tests were based off of Phase III of the round robin paper. 
Some stinput files provided had differences from the paper, and some tests needed further adjustments to match the result fluxmaps.
Edits are noted by each test.
*/

//task 4a: ommited because it is missing fluxmap file for result comparison
//stinput differences: receiver origin at {0,0,171.035}, flat facets

/*TEST_F(IsolatedHeliostatSimulationEmbree, multiFacet1332_BlockingShading4a)
{
    this->runner.set_number_of_threads(N_threads);

    set_flat_facets();

    glm::dvec3 origin = {0,0,171.035};
    set_rec_origin(origin);
    
    std::vector<int> active {1332};
    std::vector<int> blocking {1266, 1304, 1306};
    create_active_heliostats(active);
    create_blocking_heliostats(blocking);

    assign_canted_banded();

    setup_simData();
    update_from_hour("12");
    simulate_check_outputs("4a", "1");
}*/

//task 4b: 
//stinput differences: receiver origin {0,0,171.035}
//edits: no slope error
TEST_F(IsolatedHeliostatSimulationEmbree, multiFacet8993_BlockingShading4b)
{
    this->runner.set_number_of_threads(N_threads);

    set_slope_error(0.0);

    glm::dvec3 origin = {0,0,171.035};
    set_rec_origin(origin);

    std::vector<int> active {8993};
    std::vector<int> blocking {9100, 9102, 9208};
    create_active_heliostats(active);
    create_blocking_heliostats(blocking);

    assign_canted_banded(); 
    assign_focal_lengths_banded();

    setup_simData();
    update_from_hour("12");
    simulate_check_outputs("4b", "1");

}

//task 4c: 
// omitted due to significant differences in blocking heliostats between paper and stinput file. 
// Running with input file scenario would require a large amount of rays

/*TEST_F(IsolatedHeliostatSimulationEmbree, multiFacet5473_BlockingShading4c)
{
    this->runner.set_number_of_threads(N_threads);

    glm::dvec3 origin = {0,0,171.035};
    set_rec_origin(origin);
    
    std::vector<int> active {5473};
    std::vector<int> blocking {5573};
    create_active_heliostats(active);
    create_blocking_heliostats(blocking);

    assign_canted_banded();
    assign_focal_lengths_banded();

    setup_simData();
    update_from_hour("12");
    simulate_check_outputs("4c", "1");

}*/


//task 4d_8: 
//edits: no slope error
TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet8993_AimingAccuracy4d_8)
{
    this->runner.set_number_of_threads(N_threads);
    
    set_slope_error(0.0);

    std::vector<int> active {8993};
    create_active_heliostats(active);

    setup_simData();
    simulate_check_outputs("4d", "1", "8");
    
}

//task 4d_12
//edits: no slope error, no sunshape
TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet8993_AimingAccuracy4d_12)
{
    this->runner.set_number_of_threads(N_threads);
    
    set_slope_error(0.0);
    
    std::vector<int> active {8993};
    create_active_heliostats(active);

    set_no_sunShape();
    setup_simData();
    simulate_check_outputs("4d", "1", "12");
  
}

//task 4e:
TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet8993_TargetCoordSystemE4e)
{
    this->runner.set_number_of_threads(N_threads);

    glm::dvec3 shift = {-1,0,0};
    shift_aimpoint(shift);

    set_slope_error(0.0);

    std::vector<int> active {8993};
    create_active_heliostats(active);

    set_no_sunShape();
    setup_simData();

    simulate_check_outputs("4e", "1", "8");
    simulate_check_outputs("4e", "1", "12");
}

//task 4f: 
TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet8993_TargetCoordSystemW4f)
{
    this->runner.set_number_of_threads(N_threads);
    
    glm::dvec3 shift = {1,0,0};
    shift_aimpoint(shift);

    set_slope_error(0.0);

    std::vector<int> active {8993};
    create_active_heliostats(active);
    
    set_no_sunShape();
    setup_simData();
    
    simulate_check_outputs("4f", "1", "8");
    simulate_check_outputs("4f", "1", "12");

}

//task 4g:
TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet8993_TargetCoordSystemU4g)
{
    this->runner.set_number_of_threads(N_threads);
    
    glm::dvec3 shift = {0,0,1};
    shift_aimpoint(shift);

    set_slope_error(0.0);

    std::vector<int> active {8993};
    create_active_heliostats(active);
    
    set_no_sunShape();
    setup_simData();

    simulate_check_outputs("4g", "1", "8");
    simulate_check_outputs("4g", "1", "12");

}

//task 4h:
TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet8993_TargetCoordSystemD4h)
{
    this->runner.set_number_of_threads(N_threads);
    
    glm::dvec3 shift = {0,0,-1};
    shift_aimpoint(shift);

    set_slope_error(0.0);

    std::vector<int> active {8993};
    create_active_heliostats(active);
    
    set_no_sunShape();
    setup_simData();

    simulate_check_outputs("4h", "1", "8");
    simulate_check_outputs("4h", "1", "12");

}

//task 4i:
//stinput differences: flat facets 
TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet8993_TargetCoordSystem4i)
{
    this->runner.set_number_of_threads(N_threads);

    set_slope_error(0.0);
    set_flat_facets();

    std::vector<int> active {8993};
    create_active_heliostats(active);
    
    set_no_sunShape();
    setup_simData();

    simulate_check_outputs("4i", "1", "8");
    simulate_check_outputs("4i", "1", "12");

}

//task 5a: 
//stinput differences: additional blocking heliostat 5321, receiver origin {0,0,171.035}
TEST_F(IsolatedHeliostatSimulationEmbree, multiFacet1332_BlockingShading5a)
{
    this->runner.set_number_of_threads(N_threads);

    glm::dvec3 origin = {0,0,171.035};
    set_rec_origin(origin);

    std::vector<int> active {1332};
    std::vector<int> blocking {1266, 1304, 1306, 1334, 5321};
    create_active_heliostats(active);
    create_blocking_heliostats(blocking);

    assign_canted_banded();
    assign_focal_lengths_banded();

    setup_simData();

    update_from_hour("8");
    simulate_check_outputs("5a", "1");

}

//task 5b: NO FILE

//task 6a:
//stinput differences: receiver origin {0,0,171.035}
TEST_F(IsolatedHeliostatSimulationEmbree, multiFacet1332_CantingAccuracy6a)
{
    this->runner.set_number_of_threads(N_threads);

    set_flat_facets();

    glm::dvec3 origin = {0,0,171.035};
    set_rec_origin(origin);

    std::vector<int> active {1332};
    create_active_heliostats(active);

    assign_canted_banded();

    setup_simData();
    update_from_hour("12");
    simulate_check_outputs("6a", "1");

}

//task 6b: NO FILE

//task 6c: NO FILE

//task 7a:
//stinput differences: reveiver origin {0,0,171.035}
//edits: shift in aimpoint due to truncation error in legacy when round robin study generated fluxmaps
TEST_F(IsolatedHeliostatSimulationEmbree, multiFacet5473_CantingFocusingAccuracy7a) 
{
    this->runner.set_number_of_threads(N_threads);
    
    glm::dvec3 shift = {0,0,0.825};//due to truncation error in legacy
    shift_aimpoint(shift);

    set_slope_error(0.0);

    glm::dvec3 origin = {0,0,171.035};
    set_rec_origin(origin);

    std::vector<int> active {5473};
    create_active_heliostats(active);

    assign_canted_slant();

    setup_simData();
    update_from_hour("12");
    simulate_check_outputs("7a", "1");

}

//task 7b:
//stinput differences: reveiver origin {0,0,171.035}
//edits: shift in aimpoint due to truncation error in legacy when round robin study generated fluxmaps
TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet5473_AimingAccuracy7b)
{
    this->runner.set_number_of_threads(N_threads);
    
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

}

//task 7c:
//stinput differences: reveiver origin {0,0,171.035}
//edits: shift in aimpoint due to truncation error in legacy when round robin study generated fluxmaps
TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet5473_BlockingShading7c) 
{
    this->runner.set_number_of_threads(N_threads);
    
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

}
//7d: named 7b_8 in round robin files but matches 7d case
//stinput differences: reveiver origin {0,0,171.035}
TEST_F(IsolatedHeliostatSimulationEmbree, singleFacet5473_AimingAccuracy7d)
{
    this->runner.set_number_of_threads(N_threads);
    
    set_slope_error(0.0);

    std::vector<int> active {5473};
    create_active_heliostats(active);

    setup_simData();
    simulate_check_outputs("7b", "1", "8");

}
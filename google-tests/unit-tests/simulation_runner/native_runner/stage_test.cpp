#include <gtest/gtest.h>

#include <chrono>
#include <future>
#include <thread>

#include <native_runner.hpp>
#include <native_runner_types.hpp>
#include <optical_properties.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result_export.hpp>
#include <simulation_runner.hpp>

using SolTrace::NativeRunner::NativeRunner;
using SolTrace::Runner::RunnerStatus;

const int NRays = 10000;
const int MaxRays = 10000 * 10;

static SimulationData create_two_flat_elements_simulation(const bool separateStages, 
    element_ptr& plate1, element_ptr& plate2)
{
    SimulationData sd;

    SolTrace::Data::OpticalPropertySet reflective_optics(SolTrace::Data::InteractionType::REFLECTION, "reflective_plate_optics");
    reflective_optics.set_ideal_reflection(SolTrace::Data::OpticalSide::Both);
    reflective_optics.set_errors(SolTrace::Data::OpticalSide::Both, SolTrace::Data::DistributionType::NONE, 0.0, 0.0);
    auto reflective_optics_ref = sd.add_optical_property_set(reflective_optics);

    SolTrace::Data::OpticalPropertySet absorber_optics(SolTrace::Data::InteractionType::REFLECTION, "absorbing_plate_optics");
    absorber_optics.set_ideal_absorption(SolTrace::Data::OpticalSide::Both);
    auto absorber_optics_ref = sd.add_optical_property_set(absorber_optics);

    // Sun
    auto sun = make_ray_source<Sun>();
    sun->set_position(0.0, 0.0, 100.0);
    sd.add_ray_source(sun);

    // Plate 1: Angled 45 degrees, receives rays from sun and reflects toward plate 2
    // Position at origin, tilted to reflect rays in +Y direction
    plate1 = make_element<SingleElement>();
    plate1->set_origin(0, 0, 50);
    plate1->set_aim_vector(0, 50, 100);  // Tilted 45 degrees toward +Y
    plate1->set_surface(make_surface<Flat>());
    plate1->set_aperture(make_aperture<Rectangle>(10, 10));
    plate1->set_optical_property_set(reflective_optics_ref);
    plate1->set_name("plate1");

    // Plate 2: Positioned to receive reflected rays from plate 1
    // At Y=50, tilted 45 degrees to face plate 1
    plate2 = make_element<SingleElement>();
    plate2->set_origin(0, 50, 50);
    plate2->set_aim_vector(0, 0, 100);  // Tilted 45 degrees toward -Y (facing plate 1)
    plate2->set_surface(make_surface<Flat>());
    plate2->set_aperture(make_aperture<Rectangle>(5, 5));  // Smaller to absorb rays
    plate2->set_optical_property_set(absorber_optics_ref);
    plate2->set_name("plate2");    

    if (separateStages)
    {
        // Stage 0 with el1
        auto st0 = make_stage(0);
        st0->set_origin(0.0, 0.0, 0.0);
        st0->set_aim_vector(0.0, 0.0, 1.0);
        st0->set_name("Stage_0");
        st0->add_element(plate1);

        // Stage 1 with el2
        auto st1 = make_stage(1);
        st1->set_origin(0.0, 0.0, 0.0);
        st1->set_aim_vector(0.0, 0.0, 1.0);
        st1->set_name("Stage_1");
        st1->add_element(plate2);

        sd.add_stage(st0);
        sd.add_stage(st1);
    }
    else
    {
        // Single stage containing both elements
        auto st0 = make_stage(0);
        st0->set_origin(0.0, 0.0, 0.0);
        st0->set_aim_vector(0.0, 0.0, 1.0);
        st0->set_name("Stage_0");
        st0->add_element(plate1);
        st0->add_element(plate2);
        sd.add_stage(st0);
    }

    // Parameters
    SimulationParameters& params = sd.get_simulation_parameters();
    params.number_of_rays = NRays;
    params.max_number_of_rays = MaxRays;
    params.include_optical_errors = false;
    params.include_sun_shape_errors = false;
    params.seed = 12345;

    return sd;
}

static void count_hits(const RayHistoryResult& result, element_ptr& plate1, 
    element_ptr& plate2, std::vector<int>& plate1_hits, std::vector<int>& plate2_hits)
{
    plate1_hits = { 0,0,0 };
    plate2_hits = { 0,0,0 };

    const int n_records = result.get_number_of_records();
    for (int i = 0; i < n_records; i++)
    {
        ray_record_ptr rec = result[i];
        int n_interactions = rec->get_number_of_interactions();
        for (int j = 0; j < n_interactions; j++)
        {
            // Get ray event
            RayEvent rev = rec->get_event(j);
            EXPECT_FALSE(rev == RayEvent::UNKNOWN);

            // Get element id
            int el_id = rec->get_element(j);

            // Check hit
            // Plate1
            if (el_id == plate1->get_id())
            {
                // Add to count
                plate1_hits[j]++;

                // Confirm it is a reflect
                EXPECT_EQ(rev, RayEvent::REFLECT);
            }
            else if (el_id == plate2->get_id())
            {
                // Add to count
                plate2_hits[j]++;

                // Confirm is is an absorb
                EXPECT_EQ(rev, RayEvent::ABSORB);
            }
            else
            {
                // Must be create or exit
                EXPECT_TRUE(rev == RayEvent::CREATE || rev == RayEvent::EXIT);
            }
        }
    }
}

TEST(StageTest, StageOn)
{
    // Make simulation data
    element_ptr plate1, plate2;
    SimulationData sd = create_two_flat_elements_simulation(true, plate1, plate2);

    // Run simulation
    NativeRunner runner;

    RunnerStatus sts = runner.initialize();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.setup_simulation(&sd);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.run_simulation();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    // Collect results
    RayHistoryResult result;
    sts = runner.report_simulation(&result, 0);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    // Check order of hits
    int n_records = result.get_number_of_records();
    EXPECT_EQ(n_records, NRays);
    for (int i = 0; i < n_records; i++)
    {
        ray_record_ptr rec = result[i];
        int n_interactions = rec->get_number_of_interactions();
        for (int j = 0; j < n_interactions; j++)
        {
            // Get ray event
            RayEvent rev = rec->get_event(j);
            EXPECT_FALSE(rev == RayEvent::UNKNOWN);

            // Get element id
            int el_id = rec->get_element(j);

            if (j == 0)
            {
                EXPECT_EQ(rev, RayEvent::CREATE);   // Always created first
            } 
            else if (j == 1)
            {
                EXPECT_EQ(rev, RayEvent::REFLECT);  // Must hit the ideal reflective surface next
                EXPECT_EQ(el_id, plate1->get_id());
            }
            else if (j == 2)
            {
                if (rev != RayEvent::EXIT)
                {
                    EXPECT_EQ(rev, RayEvent::ABSORB);   // Must be absorbed, if not exit
                    EXPECT_EQ(el_id, plate2->get_id());
                }
            }
                
        }
    }

}

TEST(StageTest, StageOff)
{
    // Make simulation data
    element_ptr plate1, plate2;
    SimulationData sd = create_two_flat_elements_simulation(true, plate1, plate2);

    // Make native runner
    NativeRunner runner;

    // Turn off stages
    runner.disable_stages();

    // Run simulation
    RunnerStatus sts = runner.initialize();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.setup_simulation(&sd);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.run_simulation();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    // Collect results
    RayHistoryResult result;
    sts = runner.report_simulation(&result, 0);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    // Check order of hits
    int n_records = result.get_number_of_records();
    EXPECT_EQ(n_records, NRays);
    std::vector<int> plate1_hits, plate2_hits;
    count_hits(result, plate1, plate2, plate1_hits, plate2_hits);

    // Check counts
    EXPECT_EQ(plate1_hits[0], 0);   // First record should always be create
    EXPECT_EQ(plate2_hits[0], 0);   // ^
    EXPECT_GT(plate2_hits[1], 0);   // Plate 2 must have some direct hits
    EXPECT_GT(plate1_hits[1], 0);   // Plate 1 must have some direct hits
    EXPECT_GT(plate2_hits[2], 0);   // Plate 2 must have some hits reflected from 1
    EXPECT_EQ(plate1_hits[2], 0);   // Plate 1 should have no hits reflected from 2
}

TEST(StageTest, OffComparison)
{
    // Make case predefined without stages
    element_ptr plate1_nostage, plate2_nostage;
    SimulationData sd_nostage = create_two_flat_elements_simulation(false, plate1_nostage, plate2_nostage);
    NativeRunner runner_nostage;
    RunnerStatus sts = runner_nostage.initialize();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner_nostage.setup_simulation(&sd_nostage);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner_nostage.run_simulation();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    // Get results
    RayHistoryResult result_nostage;
    sts = runner_nostage.report_simulation(&result_nostage, 0);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    std::vector<int> plate1_hits_nostage = { 0,0,0 };
    std::vector<int> plate2_hits_nostage = { 0,0,0 };
    count_hits(result_nostage, plate1_nostage, plate2_nostage, plate1_hits_nostage, plate2_hits_nostage);

    // Make case with stages, but turn them off
    element_ptr plate1_stageoff, plate2_stageoff;
    SimulationData sd_stageoff = create_two_flat_elements_simulation(true, plate1_stageoff, plate2_stageoff);
    NativeRunner runner_stageoff;
    runner_stageoff.disable_stages();   // turn off stages
    sts = runner_stageoff.initialize();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner_stageoff.setup_simulation(&sd_stageoff);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner_stageoff.run_simulation();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    // Get results
    RayHistoryResult result_stageoff;
    sts = runner_stageoff.report_simulation(&result_stageoff, 0);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    std::vector<int> plate1_hits_stageoff = { 0,0,0 };
    std::vector<int> plate2_hits_stageoff = { 0,0,0 };
    count_hits(result_stageoff, plate1_stageoff, plate2_stageoff, plate1_hits_stageoff, plate2_hits_stageoff);

    // Compare results
    for (int i = 0; i < 3; i++)
    {
        EXPECT_EQ(plate1_hits_nostage[i], plate1_hits_stageoff[i]);
        EXPECT_EQ(plate2_hits_nostage[i], plate2_hits_stageoff[i]);
    }
        
}

TEST(StageTest, DuplicateStageNumbersThrows)
{
    SimulationData sd;
    auto sun = make_ray_source<Sun>();
    sun->set_position(0.0, 0.0, 100.0);
    sd.add_ray_source(sun);

    // Two stages with the same index
    auto st0a = make_stage(0);
    auto st0b = make_stage(0);

    // Add one enabled element to each stage so stages are materialized
    auto elA = make_element<SingleElement>();
    elA->set_origin(0, 0, 50);
    elA->set_aim_vector(0, 0, 100);
    elA->set_surface(make_surface<Flat>());
    elA->set_aperture(make_aperture<Rectangle>(1, 1));
    {
        SolTrace::Data::OpticalPropertySet optics(SolTrace::Data::InteractionType::REFLECTION, "duplicate_stage_A_optics");
        optics.set_ideal_reflection(SolTrace::Data::OpticalSide::Both);
        optics.set_errors(SolTrace::Data::OpticalSide::Both, SolTrace::Data::DistributionType::NONE, 0.0, 0.0);
        elA->set_optical_property_set(sd.add_optical_property_set(optics));
    }
    st0a->add_element(elA);

    auto elB = make_element<SingleElement>();
    elB->set_origin(0, 10, 50);
    elB->set_aim_vector(0, 0, 100);
    elB->set_surface(make_surface<Flat>());
    elB->set_aperture(make_aperture<Rectangle>(1, 1));
    {
        SolTrace::Data::OpticalPropertySet optics(SolTrace::Data::InteractionType::REFLECTION, "duplicate_stage_B_optics");
        optics.set_ideal_reflection(SolTrace::Data::OpticalSide::Both);
        optics.set_errors(SolTrace::Data::OpticalSide::Both, SolTrace::Data::DistributionType::NONE, 0.0, 0.0);
        elB->set_optical_property_set(sd.add_optical_property_set(optics));
    }
    st0b->add_element(elB);

    sd.add_stage(st0a);
    sd.add_stage(st0b);

    NativeRunner runner;
    RunnerStatus sts = runner.initialize();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    EXPECT_THROW(runner.setup_simulation(&sd), std::runtime_error);
}

TEST(StageTest, ElementBeforeStageThrowsWhenStagesEnabled)
{
    SimulationData sd;
    auto sun = make_ray_source<Sun>();
    sun->set_position(0.0, 0.0, 100.0);
    sd.add_ray_source(sun);

    auto el = make_element<SingleElement>();
    el->set_origin(0,0,50);
    el->set_aim_vector(0,0,100);
    el->set_surface(make_surface<Flat>());
    el->set_aperture(make_aperture<Rectangle>(1,1));
    {
        SolTrace::Data::OpticalPropertySet optics(SolTrace::Data::InteractionType::REFLECTION, "element_before_stage_optics");
        optics.set_ideal_reflection(SolTrace::Data::OpticalSide::Both);
        optics.set_errors(SolTrace::Data::OpticalSide::Both, SolTrace::Data::DistributionType::NONE, 0.0, 0.0);
        el->set_optical_property_set(sd.add_optical_property_set(optics));
    }
    sd.add_element(el);

    auto st0 = make_stage(0);

    auto elB = make_element<SingleElement>();
    elB->set_origin(0, 10, 50);
    elB->set_aim_vector(0, 0, 100);
    elB->set_surface(make_surface<Flat>());
    elB->set_aperture(make_aperture<Rectangle>(1, 1));
    {
        SolTrace::Data::OpticalPropertySet optics(SolTrace::Data::InteractionType::REFLECTION, "element_before_stage_B_optics");
        optics.set_ideal_reflection(SolTrace::Data::OpticalSide::Both);
        optics.set_errors(SolTrace::Data::OpticalSide::Both, SolTrace::Data::DistributionType::NONE, 0.0, 0.0);
        elB->set_optical_property_set(sd.add_optical_property_set(optics));
    }
    st0->add_element(elB);

    sd.add_stage(st0);

    NativeRunner runner;
    RunnerStatus sts = runner.initialize();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    EXPECT_THROW(runner.setup_simulation(&sd), std::runtime_error);
}

TEST(StageTest, ElementBeforeStageSucceedsWhenStagesDisabled)
{
    SimulationData sd;
    auto sun = make_ray_source<Sun>();
    sun->set_position(0.0, 0.0, 100.0);
    sd.add_ray_source(sun);

    auto el = make_element<SingleElement>();
    el->set_origin(0,0,50);
    el->set_aim_vector(0,0,100);
    el->set_surface(make_surface<Flat>());
    el->set_aperture(make_aperture<Rectangle>(1,1));
    {
        SolTrace::Data::OpticalPropertySet optics(SolTrace::Data::InteractionType::REFLECTION, "element_before_stage_disabled_optics");
        optics.set_ideal_reflection(SolTrace::Data::OpticalSide::Both);
        optics.set_errors(SolTrace::Data::OpticalSide::Both, SolTrace::Data::DistributionType::NONE, 0.0, 0.0);
        el->set_optical_property_set(sd.add_optical_property_set(optics));
    }
    sd.add_element(el);

    auto st0 = make_stage(0);

    auto elB = make_element<SingleElement>();
    elB->set_origin(0, 10, 50);
    elB->set_aim_vector(0, 0, 100);
    elB->set_surface(make_surface<Flat>());
    elB->set_aperture(make_aperture<Rectangle>(1, 1));
    {
        SolTrace::Data::OpticalPropertySet optics(SolTrace::Data::InteractionType::REFLECTION, "element_before_stage_disabled_B_optics");
        optics.set_ideal_reflection(SolTrace::Data::OpticalSide::Both);
        optics.set_errors(SolTrace::Data::OpticalSide::Both, SolTrace::Data::DistributionType::NONE, 0.0, 0.0);
        elB->set_optical_property_set(sd.add_optical_property_set(optics));
    }
    st0->add_element(elB);

    sd.add_stage(st0);

    NativeRunner runner;
    runner.disable_stages();

    RunnerStatus sts = runner.initialize();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.setup_simulation(&sd);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.run_simulation();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
}

TEST(StageTest, NoStagesCreatesSingleInternalStage)
{
    SimulationData sd;
    sd.get_simulation_parameters().number_of_rays = NRays;
    sd.get_simulation_parameters().max_number_of_rays = MaxRays;
    auto sun = make_ray_source<Sun>();
    sun->set_position(0.0, 0.0, 100.0);
    sd.add_ray_source(sun);

    auto el = make_element<SingleElement>();
    el->set_origin(0,0,50);
    el->set_aim_vector(0,0,100);
    el->set_surface(make_surface<Flat>());
    el->set_aperture(make_aperture<Rectangle>(2,2));
    {
        SolTrace::Data::OpticalPropertySet optics(SolTrace::Data::InteractionType::REFLECTION, "no_stage_internal_optics");
        optics.set_ideal_reflection(SolTrace::Data::OpticalSide::Both);
        optics.set_errors(SolTrace::Data::OpticalSide::Both, SolTrace::Data::DistributionType::NONE, 0.0, 0.0);
        el->set_optical_property_set(sd.add_optical_property_set(optics));
    }
    sd.add_element(el);

    NativeRunner runner;
    RunnerStatus sts = runner.initialize();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.setup_simulation(&sd);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.run_simulation();
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);

    RayHistoryResult result;
    sts = runner.report_simulation(&result, 0);
    ASSERT_EQ(sts, RunnerStatus::SUCCESS);
    EXPECT_EQ(result.get_number_of_records(), NRays);
}

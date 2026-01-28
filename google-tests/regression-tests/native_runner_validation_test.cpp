#include <gtest/gtest.h>

#include <chrono>
#include <sstream>

#include <native_runner.hpp>
#include <simulation_data_export.hpp>
#include <simulation_result_export.hpp>

#include "common.hpp"
#include "split_csv.h"

using SolTrace::Runner::RunnerStatus;

using SolTrace::NativeRunner::NativeRunner;
using SolTrace::NativeRunner::TRayData;
using SolTrace::NativeRunner::tstage_ptr;
using SolTrace::NativeRunner::TSystem;

using SolTrace::Result::RayEvent;

bool get_runner_element_and_stage(const NativeRunner *runner,
                                  element_id id,
                                  int_fast64_t &re,
                                  int_fast64_t &rs)
{
    // Assumes that `id` is a unique identifier and that an element
    // can belong only to a single stage
    bool found = false;
    const TSystem *sys = runner->get_system();
    for (unsigned stage_idx = 0;
         stage_idx < sys->StageList.size();
         ++stage_idx)
    {
        tstage_ptr stage = sys->StageList[stage_idx];
        for (unsigned element_idx = 0;
             element_idx < stage->ElementList.size();
             ++element_idx)
        {
            SolTrace::NativeRunner::telement_ptr tel =
                stage->ElementList[element_idx];
            if (tel->sim_data_id == id)
            {
                // NativeRunner stage and element ids are 1 based
                // so we add 1 to the vector index here.
                found = true;
                re = element_idx + 1;
                rs = stage_idx + 1;
            }
        }

        if (found)
            break;
    }

    return found;
}

int_fast64_t count_element_event(const SimulationResult &res, element_id el, RayEvent rev)
{
    int_fast64_t count = 0;

    for (auto ray_idx = 0;
         ray_idx < res.get_number_of_records();
         ++ray_idx)
    {
        auto rr = res[ray_idx];
        for (auto event_idx = 0;
             event_idx < rr->get_number_of_interactions();
             ++event_idx)
        {
            if (rr->get_event(event_idx) == rev &&
                rr->get_element(event_idx) == el)
            {
                ++count;
            }
        }
    }

    return count;
}

TEST(NativeRunner, ValidationTest1)
{
    // Pulling in path variable from CMake and creating path to .stinput sample file
    std::string sample_path = std::string(PROJECT_DIR) + std::string("/High Flux Solar Furnace.stinput");

    // Path to .csv exported from Soltrace as ground truth
    std::string ground_csv_path = PROJECT_DIR + std::string("/hfsf_example_raydata.csv");

    std::ifstream csv_file(ground_csv_path);
    std::vector<std::vector<std::string>> ground_raydata = split_csv(ground_csv_path);

    // Create Simuluation Data
    SimulationData sd;

    // Constants
    const uint_fast64_t NRAYS = 10000;
    const double TOL = 1e-4;

    // Read Input File
    bool success = sd.import_from_file(sample_path);
    EXPECT_TRUE(success);
    EXPECT_TRUE(sd.get_number_of_elements() > 0);
    EXPECT_TRUE(sd.get_number_of_ray_sources() > 0);

    // Parameters
    SimulationParameters &params = sd.get_simulation_parameters();
    params.include_optical_errors = false;
    params.include_sun_shape_errors = false;
    params.max_number_of_rays = NRAYS * 100;
    params.number_of_rays = NRAYS;
    params.seed = 1;

    // Run Ray Trace
    NativeRunner runner;
    runner.disable_point_focus();
    runner.disable_power_tower();
    RunnerStatus sts;
    sts = runner.initialize();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.setup_simulation(&sd);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);

    auto stage_list = runner.get_system()->StageList;
    tstage_ptr tstage = nullptr;
    for (auto iter = stage_list.begin();
         iter != stage_list.end();
         ++iter)
    {
        tstage = *iter;
        tstage->MultiHitsPerRay = false;
        tstage->TraceThrough = false;
    }

    sts = runner.run_simulation();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);

    // const TSystem *sys = runner.get_system();
    // const TRayData *ray_data = &(runner.get_system()->RayData);
    // size_t nrdata = ray_data->Count();

    // ray_data->Print();

    SimulationResult result;
    sts = runner.report_simulation(&result, 0);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    EXPECT_EQ(result.get_number_of_records(), NRAYS);

    Vector3d point, cosines;
    Vector3d pos_stage, dir_stage;
    Vector3d temp;
    int_fast64_t element;
    int_fast64_t stage;
    uint_fast64_t rayidx;
    uint_fast64_t iidx;
    element_ptr el = nullptr;
    int_fast64_t run_element, run_stage;

    // size_t total_lines = ground_raydata[0].size();
    size_t total_lines = std::numeric_limits<size_t>::max();
    size_t i;
    for (i = 0; i < 9; ++i)
    {
        total_lines = std::min(total_lines, ground_raydata[i].size());
    }

    // Compare saved CSV values to runner values
    for (i = 1; i < total_lines; ++i)
    {
        element = stoi(ground_raydata[6][i]);
        stage = stoi(ground_raydata[7][i]);
        // Legacy SolTrace and CSV file had 1-based ray IDs. SimulationResult
        // has 0-based ray ID's so subtract 1 here.
        rayidx = stoul(ground_raydata[8][i]) - 1;

        const ray_record_ptr rr = result[rayidx];
        if (element > 0)
        {
            bool found = false;
            for (uint_fast64_t idx = 0; idx < rr->get_number_of_interactions(); ++idx)
            {
                run_element = -1;
                run_stage = -1;
                bool sts = get_runner_element_and_stage(&runner,
                                                        rr->get_element(idx),
                                                        run_element,
                                                        run_stage);
                if (sts && element == run_element && stage == run_stage)
                {
                    iidx = idx;
                    found = true;
                    break;
                }
            }

            // assert(found);
            EXPECT_TRUE(found);

            if (found)
            {
                EXPECT_NE(rr->get_event(iidx), RayEvent::CREATE);
                EXPECT_NE(rr->get_event(iidx), RayEvent::ABSORB);
                EXPECT_NE(rr->get_event(iidx), RayEvent::EXIT);

                if (rr->get_event(iidx) == RayEvent::ABSORB)
                {
                    break;
                }
            }
            else
            {
                std::cout << "CSV Line: " << i + 1
                          << "Ray Record: " << *rr
                          << std::endl;
                break;
            }

            rr->get_position(iidx, point);
            // Legacy SolTrace stored the incoming ray direction whereas
            // NativeRunner/SimulationResult stores the exit direction so
            // we take the direction for the previous ray event.
            rr->get_direction(iidx - 1, cosines);
        }
        else
        {
            iidx = rr->get_number_of_interactions() - 1;
            if (element == 0)
            {
                // Ray miss -- check only that it was noted
                EXPECT_EQ(rr->get_event(iidx), RayEvent::EXIT);
                continue;
            }
            else
            {
                EXPECT_EQ(rr->get_event(iidx), RayEvent::ABSORB);
                get_runner_element_and_stage(&runner,
                                             rr->get_element(iidx),
                                             run_element,
                                             run_stage);
                EXPECT_EQ(run_element, abs(element));
                EXPECT_EQ(run_stage, stage);

                rr->get_position(iidx, point);
                // Legacy SolTrace stored the incoming ray direction whereas
                // NativeRunner/SimulationResult stores the exit direction so
                // we take the direction for the previous ray event.
                rr->get_direction(iidx - 1, cosines);
            }
        }

        el = sd.get_element(rr->get_element(iidx));
        EXPECT_NE(el, nullptr);

        if (el == nullptr)
        {
            std::cout << "CSV Line: " << i + 1
                      << "\nRay Number: " << rayidx + 1
                      << "\nElement: " << rr->get_element(iidx)
                      << " CSV Element: " << element
                      << " Runner Element: " << run_element
                      << "\nCSV Stage: " << stage
                      << " Runner Stage: " << run_stage
                      << std::endl;
            break;
        }

        // Runner and SimulationResult store everything in global
        // coordinate whereas the CSV file is in stage coordinates
        // as per legacy SolTrace
        el->convert_global_to_reference(pos_stage, point);

        // See previous comment about coordinates
        el->convert_vector_global_to_reference(dir_stage, cosines);

        EXPECT_NEAR(pos_stage[0], stod(ground_raydata[0][i]), TOL);
        EXPECT_NEAR(pos_stage[1], stod(ground_raydata[1][i]), TOL);
        EXPECT_NEAR(pos_stage[2], stod(ground_raydata[2][i]), TOL);

        EXPECT_NEAR(dir_stage[0], stod(ground_raydata[3][i]), TOL);
        EXPECT_NEAR(dir_stage[1], stod(ground_raydata[4][i]), TOL);
        EXPECT_NEAR(dir_stage[2], stod(ground_raydata[5][i]), TOL);

        // if (fabs(pos_stage[0] - stod(ground_raydata[0][i]) > TOL))
        // {
        //     Vector3d pos_csv(stod(ground_raydata[0][i]),
        //                      stod(ground_raydata[1][i]),
        //                      stod(ground_raydata[2][i]));
        //     Vector3d dir_csv(stod(ground_raydata[3][i]),
        //                      stod(ground_raydata[4][i]),
        //                      stod(ground_raydata[5][i]));

        //     Vector3d pos_loc;
        //     Vector3d dir_loc;

        //     Vector3d csv_glob;
        //     el->convert_stage_to_local(temp, pos_csv);
        //     el->convert_local_to_global(csv_glob, temp);

        //     el->convert_global_to_local(pos_loc, point);
        //     el->convert_vector_global_to_local(dir_loc, cosines);

        //     std::cout << "CSV Line: " << i + 1
        //               << "\nRay Number: " << rayidx + 1
        //               << "\nElement: " << rr->get_element(iidx)
        //               << " CSV Element: " << element
        //               << " Runner Element: " << run_element
        //               << "\nCSV Stage: " << stage
        //               << " Runner Stage: " << run_stage
        //               << std::endl;

        //     std::cout << "CSV Pos: " << pos_csv
        //               << "\nCSV Global: " << csv_glob
        //               << "\nPos Global: " << point
        //               << "\nPos Stage: " << pos_stage
        //               << "\nPos Local: " << pos_loc
        //               << "\nCSV Dir: " << dir_csv
        //               << "\nDir Global: " << cosines
        //               << "\nDir Stage: " << dir_stage
        //               << "\nDir Local: " << dir_loc
        //               << std::endl;

        //     std::cout << "Ray Record: "
        //               << *rr
        //               << std::endl;

        //     std::cout << "Element Global Origin: " << el->get_origin_global()
        //               << "\nElement Local to Global: " << el->get_local_to_global()
        //               << "\nElement Ref to Local: " << el->get_local_to_reference()
        //               << std::endl;

        //     break;
        // }
    }
}

TEST(NativeRunner, ValidationTest2)
{
    // Pulling in path variable from CMake and creating path to .stinput sample file
    std::string sample_path = std::string(PROJECT_DIR) + std::string("/Power-tower-surround_singlefacet.stinput");

    // Path to .csv exported from Soltrace as ground truth
    std::string ground_csv_path = PROJECT_DIR + std::string("/powertower_example_raydata.csv");

    std::ifstream csv_file(ground_csv_path);
    std::vector<std::vector<std::string>> ground_raydata = split_csv(ground_csv_path);

    // Create Simuluation Data
    SimulationData sd;

    // Constants
    const uint_fast64_t NRAYS = 50000;
    const double TOL = 1e-4;

    // Read Input File
    bool success = sd.import_from_file(sample_path);
    EXPECT_TRUE(success);
    EXPECT_TRUE(sd.get_number_of_elements() > 0);
    EXPECT_TRUE(sd.get_number_of_ray_sources() > 0);

    std::cout << "Num Elements: " << sd.get_number_of_elements() << std::endl;

    // Parameters
    SimulationParameters &params = sd.get_simulation_parameters();
    params.include_optical_errors = false;
    params.include_sun_shape_errors = false;
    params.max_number_of_rays = NRAYS * 100;
    params.number_of_rays = NRAYS;
    params.seed = 1;

    // for (auto cit = sd.get_const_iterator();
    //      !sd.is_at_end(cit);
    //      ++cit)
    // {
    //     auto elem = cit->second;
    //     if (elem->is_stage())
    //     {
    //         stage_ptr st = dynamic_pointer_cast<StageElement>(elem);
    //         assert(st != nullptr);
    //         std::cout << "Stage: " << st->get_stage()
    //                   << "\nNumber of elements: "
    //                   << st->get_number_of_elements()
    //                   << std::endl;
    //         if (st->get_stage() == 1)
    //         {
    //             auto el = st->get_const_iterator()->second;

    //             auto surf = dynamic_pointer_cast<Cylinder>(el->get_surface());
    //             assert(surf != nullptr);
    //             auto ap = dynamic_pointer_cast<Rectangle>(el->get_aperture());
    //             assert(ap != nullptr);

    //             std::cout << "------------\n"
    //                       << "Element ID: " << el->get_id()
    //                       << "\nElement name: " << el->get_name()
    //                       << "\nIs Stage: " << el->is_stage()
    //                       << "\nIs Composite: " << el->is_composite()
    //                       << "\nIs Single: " << el->is_single()
    //                       << "\nOrigin (ref): " << el->get_origin_ref()
    //                       << "\nOrigin (stage): " << el->get_origin_stage()
    //                       << "\nOrigin (global): " << el->get_origin_global()
    //                       << "\nAim (ref): " << el->get_aim_vector_ref()
    //                       << "\nAim (stage): " << el->get_aim_vector_stage()
    //                       << "\nAim (global): " << el->get_aim_vector_global()
    //                       << "\nRefToLoc: " << el->get_reference_to_local()
    //                       << "\nCylinder Radius: " << surf->radius
    //                       << "\nAperture Dimension: " << ap->x_length << " x " << ap->y_length
    //                       << " at (" << ap->x_coord << ", " << ap->y_coord << ")"
    //                       << "\n** Front Optical Properties **\n"
    //                       << *(el->get_front_optical_properties())
    //                       << "\n** Back Optical Properties **\n"
    //                       << *(el->get_back_optical_properties())
    //                       << "\n------------"
    //                       << std::endl;
    //         }
    //     }
    // }

    // return;

    // Run Ray Trace
    NativeRunner runner;
    runner.disable_point_focus();
    runner.disable_power_tower();
    RunnerStatus sts;
    sts = runner.initialize();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    sts = runner.setup_simulation(&sd);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);

    // auto stage_list = runner.get_system()->StageList;
    // stage_list[1]->MultiHitsPerRay = false;
    // tstage_ptr tstage = nullptr;
    // for (auto iter = stage_list.begin();
    //      iter != stage_list.end();
    //      ++iter)
    // {
    //     tstage = *iter;
    //     tstage->MultiHitsPerRay = false;
    //     tstage->TraceThrough = false;
    // }

    auto t0 = std::chrono::high_resolution_clock::now();
    sts = runner.run_simulation();
    auto t1 = std::chrono::high_resolution_clock::now();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);

    std::chrono::duration<double, std::milli> dur = t1 - t0;
    EXPECT_TRUE(dur.count() < 75000.0);

    std::cout << "Time: " << dur.count() << " ms" << std::endl;

    // const TSystem *sys = runner.get_system();
    // const TRayData *ray_data = &(runner.get_system()->RayData);
    // size_t nrdata = ray_data->Count();

    // ray_data->Print();

    SimulationResult result;
    sts = runner.report_simulation(&result, 0);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    EXPECT_EQ(result.get_number_of_records(), NRAYS);

    element_id absorber_id = 6285;
    int_fast64_t nabsorbed = count_element_event(result, absorber_id, RayEvent::ABSORB);
    int_fast64_t nreflect = count_element_event(result, absorber_id, RayEvent::REFLECT);
    int_fast64_t nevents = nabsorbed + nreflect;

    std::cout << "Total: " << nevents
              << "\nAbsorb: " << nabsorbed << " ("
              << static_cast<double>(nabsorbed) / nevents << ")"
              << "\nReflect: " << nreflect << " ("
              << static_cast<double>(nreflect) / nevents << ")"
              << std::endl;

    // result.write_csv_file("native_runner_result_dump.csv");

    Vector3d point, cosines;
    Vector3d pos_stage, dir_stage;
    Vector3d temp;
    int_fast64_t element;
    int_fast64_t stage;
    uint_fast64_t rayidx;
    uint_fast64_t iidx;
    element_ptr el = nullptr;
    int_fast64_t run_element, run_stage;
    double RTOL = 0.0;

    // size_t total_lines = ground_raydata[0].size();
    size_t total_lines = std::numeric_limits<size_t>::max();
    size_t i;
    for (i = 0; i < 9; ++i)
    {
        total_lines = std::min(total_lines, ground_raydata[i].size());
    }

    // Compare saved CSV values to runner values
    for (i = 1; i < total_lines; ++i)
    {
        element = stoi(ground_raydata[6][i]);
        stage = stoi(ground_raydata[7][i]);
        // Legacy SolTrace and CSV file had 1-based ray IDs. SimulationResult
        // has 0-based ray ID's so subtract 1 here.
        rayidx = stoul(ground_raydata[8][i]) - 1;

        const ray_record_ptr rr = result[rayidx];
        if (element > 0)
        {
            bool found = false;
            for (uint_fast64_t idx = 0; idx < rr->get_number_of_interactions(); ++idx)
            {
                run_element = -1;
                run_stage = -1;
                bool sts = get_runner_element_and_stage(&runner,
                                                        rr->get_element(idx),
                                                        run_element,
                                                        run_stage);
                if (sts && element == run_element && stage == run_stage)
                {
                    iidx = idx;
                    found = true;
                    break;
                }
            }

            // assert(found);
            EXPECT_TRUE(found);

            if (found)
            {
                EXPECT_NE(rr->get_event(iidx), RayEvent::CREATE);
                EXPECT_NE(rr->get_event(iidx), RayEvent::ABSORB);
                EXPECT_NE(rr->get_event(iidx), RayEvent::EXIT);
            }

            if (!found ||
                rr->get_event(iidx) == RayEvent::CREATE ||
                rr->get_event(iidx) == RayEvent::ABSORB ||
                rr->get_event(iidx) == RayEvent::EXIT)
            {
                std::cout << "CSV Line: " << i + 1
                          << "\nRay Record: " << *rr
                          << std::endl;
                break;
            }

            rr->get_position(iidx, point);
            // Legacy SolTrace stored the incoming ray direction whereas
            // NativeRunner/SimulationResult stores the exit direction so
            // we take the direction for the previous ray event.
            rr->get_direction(iidx - 1, cosines);
        }
        else
        {
            iidx = rr->get_number_of_interactions() - 1;
            if (element == 0)
            {
                // Ray miss -- check only that it was noted
                EXPECT_EQ(rr->get_event(iidx), RayEvent::EXIT);
                continue;
            }
            else
            {
                EXPECT_EQ(rr->get_event(iidx), RayEvent::ABSORB);

                get_runner_element_and_stage(&runner,
                                             rr->get_element(iidx),
                                             run_element,
                                             run_stage);
                EXPECT_EQ(run_element, abs(element));
                EXPECT_EQ(run_stage, stage);

                if (run_element != abs(element) || run_stage != stage)
                {
                    std::cout << "CSV Line: " << i + 1
                              << "\nRay Number: " << rayidx + 1
                              << "\nElement: " << rr->get_element(iidx)
                              << " CSV Element: " << element
                              << " Runner Element: " << run_element
                              << "\nCSV Stage: " << stage
                              << " Runner Stage: " << run_stage
                              << "\nRay Record: " << *rr
                              << std::endl;
                    break;
                }

                rr->get_position(iidx, point);
                // Legacy SolTrace stored the incoming ray direction whereas
                // NativeRunner/SimulationResult stores the exit direction so
                // we take the direction for the previous ray event.
                rr->get_direction(iidx - 1, cosines);
            }
        }

        el = sd.get_element(rr->get_element(iidx));
        EXPECT_NE(el, nullptr);

        if (el == nullptr)
        {
            std::cout << "CSV Line: " << i + 1
                      << "\nRay Number: " << rayidx + 1
                      << "\nElement: " << rr->get_element(iidx)
                      << " CSV Element: " << element
                      << " Runner Element: " << run_element
                      << "\nCSV Stage: " << stage
                      << " Runner Stage: " << run_stage
                      << "\nRay Record: " << *rr
                      << std::endl;
            break;
        }

        // Runner and SimulationResult store everything in global
        // coordinate whereas the CSV file is in stage coordinates
        // as per legacy SolTrace
        el->convert_global_to_reference(pos_stage, point);

        // See previous comment about coordinates
        el->convert_vector_global_to_reference(dir_stage, cosines);

        RTOL = vector_norm(pos_stage) * TOL;
        EXPECT_NEAR(pos_stage[0], stod(ground_raydata[0][i]), RTOL);
        EXPECT_NEAR(pos_stage[1], stod(ground_raydata[1][i]), RTOL);
        EXPECT_NEAR(pos_stage[2], stod(ground_raydata[2][i]), RTOL);

        EXPECT_NEAR(dir_stage[0], stod(ground_raydata[3][i]), TOL);
        EXPECT_NEAR(dir_stage[1], stod(ground_raydata[4][i]), TOL);
        EXPECT_NEAR(dir_stage[2], stod(ground_raydata[5][i]), TOL);

        if (fabs(pos_stage[0] - stod(ground_raydata[0][i])) > RTOL)
        {
            std::cout << "CSV Line: " << i + 1
                      << "\nRay Number: " << rayidx + 1
                      << "\nElement: " << rr->get_element(iidx)
                      << " CSV Element: " << element
                      << " Runner Element: " << run_element
                      << "\nCSV Stage: " << stage
                      << " Runner Stage: " << run_stage
                      << "\nRay Record: " << *rr
                      << std::endl;
            break;
        }

        // -------- Test Debugging Code -------- //
        // if (fabs(pos_stage[0] - stod(ground_raydata[0][i]) > TOL))
        // {
        //     Vector3d pos_csv(stod(ground_raydata[0][i]),
        //                      stod(ground_raydata[1][i]),
        //                      stod(ground_raydata[2][i]));
        //     Vector3d dir_csv(stod(ground_raydata[3][i]),
        //                      stod(ground_raydata[4][i]),
        //                      stod(ground_raydata[5][i]));

        //     Vector3d pos_loc;
        //     Vector3d dir_loc;

        //     Vector3d csv_glob;
        //     el->convert_stage_to_local(temp, pos_csv);
        //     el->convert_local_to_global(csv_glob, temp);

        //     el->convert_global_to_local(pos_loc, point);
        //     el->convert_vector_global_to_local(dir_loc, cosines);

        //     std::cout << "CSV Line: " << i + 1
        //               << "\nRay Number: " << rayidx + 1
        //               << "\nElement: " << rr->get_element(iidx)
        //               << " CSV Element: " << element
        //               << " Runner Element: " << run_element
        //               << "\nCSV Stage: " << stage
        //               << " Runner Stage: " << run_stage
        //               << std::endl;

        //     std::cout << "CSV Pos: " << pos_csv
        //               << "\nCSV Global: " << csv_glob
        //               << "\nPos Global: " << point
        //               << "\nPos Stage: " << pos_stage
        //               << "\nPos Local: " << pos_loc
        //               << "\nCSV Dir: " << dir_csv
        //               << "\nDir Global: " << cosines
        //               << "\nDir Stage: " << dir_stage
        //               << "\nDir Local: " << dir_loc
        //               << std::endl;

        //     std::cout << "Ray Record: "
        //               << *rr
        //               << std::endl;

        //     std::cout << "Element Global Origin: " << el->get_origin_global()
        //               << "\nElement Local to Global: " << el->get_local_to_global()
        //               << "\nElement Ref to Local: " << el->get_local_to_reference()
        //               << std::endl;

        //     break;
        // }
    }
}

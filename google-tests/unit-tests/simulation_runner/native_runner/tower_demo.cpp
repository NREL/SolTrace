#include <gtest/gtest.h>

// #include <constants.hpp>
// #include <element.hpp>
#include <native_runner.hpp>
// #include <sun.hpp>
// #include <simulation_data.hpp>
#include <simulation_data_export.hpp>
// #include <single_element.hpp>
// #include <stage_element.hpp>
// #include <vector3d.hpp>

#include <cmath>
#include <iostream>
#include <sstream>

#include "common.hpp"
#include "count_absorbed_native.h"

using SolTrace::Runner::RunnerStatus;
using SolTrace::NativeRunner::NativeRunner;

TEST(TowerDemo, NativeRunnerWithStages)
{
    SimulationData sd;

    // Sun
    auto sun = make_ray_source<Sun>();
    sun->set_position(0.0, 0.0, 100.0);
    sd.add_ray_source(sun);

    // Absorber -- Flat
    auto absorber = make_element<SingleElement>();
    absorber->set_origin(0.0, 0.0, 10.0);
    absorber->set_aim_vector(0.0, 5.0, 0.0);
    absorber->set_zrot(0.0);
    absorber->compute_coordinate_rotations();
    absorber->set_surface(make_surface<Flat>()); // surface(nullptr)
    absorber->set_aperture(make_aperture<Rectangle>(2.0, 2.0));
    OpticalProperties *foptics = absorber->get_front_optical_properties();
    foptics->my_type = InteractionType::REFLECTION;
    foptics->reflectivity = 0.0;
    absorber->set_name("Absorber");

    // // Absorber -- Cylindrical -- MAY NOT WORK DUE TO UNIMPLEMENTED CODE!
    // auto absorber = make_element<SingleElement>();
    // const double r = 5.0;
    // absorber->set_origin(0.0, -r, 0.0);
    // absorber->set_zrot(0.0);
    // absorber->set_aim_vector(0.0, 100.0, -(r + 50.0));
    // // These are the default options but to show off the constructor
    // OpticalProperties op(REFLECTION, 0.0, 0.0, 0.0, 0.0);
    // // Alternative to getting the pointer and setting values
    // absorber->set_front_optical_properties(op);
    // absorber->set_surface(make_surface<Cylinder>());
    // absorber->set_aperture(make_aperture<SingleAxisCurvatureSection>());

    // Make stage 1 -- second stage -- these can be added to SimulationData
    // in any order but should be numbered in the desired order
    auto st1 = make_stage(1);
    // Origin is initialized to zero but set it explicitly
    st1->set_origin(0.0, 0.0, 0.0);
    // Set aim vector so stage and global coordinates are identical
    st1->set_aim_vector(0.0, 0.0, 1.0);
    st1->set_zrot(0.0);
    st1->compute_coordinate_rotations();
    st1->add_element(absorber);
    // Optional -- to help the user identify things
    st1->set_name("Stage 1--Absorber");

    // Make stage 0 -- this will be the first stage if the runner uses stages
    auto st0 = make_stage(0);
    st0->set_origin(0.0, 0.0, 0.0);
    st0->set_aim_vector(0.0, 0.0, 1.0);
    st0->set_zrot(0.0);
    st0->compute_coordinate_rotations();
    st0->set_name("Stage 0--Reflectors");

    Vector3d rvec, svec, avec;
    Vector3d aim, pos;

    for (int i = -1; i < 2; ++i)
    {
        auto el = make_element<SingleElement>();
        foptics = el->get_front_optical_properties();
        foptics->reflectivity = 1.0;

        pos.set_values(5 * sin(i * PI / 2.0),
                       5 * cos(i * PI / 2.0),
                       0.0);
        el->set_origin(pos);
        vector_add(1.0, absorber->get_origin_global(),
                   -1.0, pos,
                   rvec);
        make_unit_vector(rvec);
        svec = sun->get_position();
        make_unit_vector(svec);
        vector_add(0.5, rvec, 0.5, svec, avec);

        vector_add(1.0, pos, 100.0, avec, aim);
        el->set_aim_vector(aim);

        // TODO: Set zrot as in python file?
        el->set_zrot(30.0 * i);
        // // Could also be set in radians
        // el->set_zrot_radians(PI / 3.0 * i);
        el->compute_coordinate_rotations();

        el->set_surface(make_surface<Flat>());
        el->set_aperture(make_aperture<Rectangle>(1.0, 1.95));

        std::stringstream name;
        name << "Reflector " << i;
        el->set_name(name.str());

        st0->add_element(el);
    }

    // Stages must have all elements present before adding to SimulationData!!
    sd.add_stage(st0);
    sd.add_stage(st1);

    // Set parameters
    SimulationParameters &params = sd.get_simulation_parameters();
    params.number_of_rays = 100000;
    // params.number_of_rays = 10; // Above takes too long
    params.max_number_of_rays = params.number_of_rays * 100;
    params.include_optical_errors = false;
    params.include_sun_shape_errors = false;
    params.seed = 12345;

    // // We can go over all the elements added
    // for (auto iter = sd.get_iterator();
    //      !sd.is_at_end(iter);
    //      ++iter)
    // {
    //     // iter is a iterator over the storing container which is a map
    //     // so that the iterator gives the key value pair
    //     element_id id = iter->first;
    //     // `element_ptr` is a std::shared_pointer to an Element
    //     element_ptr el = iter->second;
    //     std::cout << "------------\n"
    //               << "Element ID: " << id
    //               << "\nElement name: " << el->get_name()
    //               << "\nIs Stage: " << el->is_stage()
    //               << "\nIs Composite: " << el->is_composite()
    //               << "\nIs Single: " << el->is_single()
    //               // Below are all the same in this case
    //               << "\nOrigin (local): " << el->get_origin_local()
    //               << "\nOrigin (stage): " << el->get_origin_stage()
    //               << "\nOrigin (global): " << el->get_origin_global()
    //               << "\n";
    //     // NOTE: The above iteration includes StageElements (e.g. stages),
    //     // CompositeElements (which are collections of elements), as well as
    //     // SingleElements (actual physical elements that interact with rays).
    //     // We can distinguish between them by calling the various `is_XXXX`
    //     // methods as seen above. For runners that put everything in global
    //     // coordinates, anything that is not a SingleElement (so `is_single`
    //     // returns true) can be ignored.
    // }

    NativeRunner runner;
    RunnerStatus sts = runner.initialize();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    // Setup runs but is not complete
    sts = runner.setup_simulation(&sd);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    // Run simulation runs but returns RunnerStatus::ERROR
    sts = runner.run_simulation();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);

    // TODO: Do some post processing tests here

    // const TSystem *sys = runner.get_system();
    // // auto ray_data = sys->AllRayData;
    // sys->AllRayData.Print();
}

TEST(TowerDemo, NativeRunnerWithoutStages)
{
    SimulationData sd;

    // Sun
    auto sun = make_ray_source<Sun>();
    sun->set_position(0.0, 0.0, 100.0);
    sd.add_ray_source(sun);

    // Absorber -- Flat
    auto absorber = make_element<SingleElement>();
    absorber->set_origin(0.0, 0.0, 10.0);
    absorber->set_aim_vector(0.0, 5.0, 0.0);
    absorber->set_zrot(0.0);
    absorber->compute_coordinate_rotations();
    absorber->set_surface(make_surface<Flat>()); // surface(nullptr)
    absorber->set_aperture(make_aperture<Rectangle>(2.0, 2.0));
    OpticalProperties *foptics = absorber->get_front_optical_properties();
    foptics->my_type = InteractionType::REFLECTION;
    foptics->reflectivity = 0.0;
    absorber->set_name("Absorber");
    sd.add_element(absorber);

    Vector3d rvec, svec, avec;
    Vector3d aim, pos;

    for (int i = -1; i < 2; ++i)
    {
        auto el = make_element<SingleElement>();
        foptics = el->get_front_optical_properties();
        foptics->reflectivity = 1.0;

        pos.set_values(5 * sin(i * PI / 2.0),
                       5 * cos(i * PI / 2.0),
                       0.0);
        el->set_origin(pos);
        vector_add(1.0, absorber->get_origin_global(),
                   -1.0, pos,
                   rvec);
        make_unit_vector(rvec);
        svec = sun->get_position();
        make_unit_vector(svec);
        vector_add(0.5, rvec, 0.5, svec, avec);

        vector_add(1.0, pos, 100.0, avec, aim);
        el->set_aim_vector(aim);

        // TODO: Set zrot as in python file?
        el->set_zrot(30.0 * i);
        // // Could also be set in radians
        // el->set_zrot_radians(PI / 3.0 * i);
        el->compute_coordinate_rotations();

        el->set_surface(make_surface<Flat>());
        el->set_aperture(make_aperture<Rectangle>(1.0, 1.95));

        std::stringstream name;
        name << "Reflector " << i;
        el->set_name(name.str());

        sd.add_element(el);
    }

    // Set parameters
    SimulationParameters &params = sd.get_simulation_parameters();
    params.number_of_rays = 100000;
    params.max_number_of_rays = params.number_of_rays * 100;
    params.include_optical_errors = false;
    params.include_sun_shape_errors = false;
    params.seed = 12345;

    NativeRunner runner;
    RunnerStatus sts = runner.initialize();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    // Setup runs but is not complete
    sts = runner.setup_simulation(&sd);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    // Run simulation runs but returns RunnerStatus::ERROR
    sts = runner.run_simulation();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);

    // TODO: Do some post processing tests here

    // const TSystem *sys = runner.get_system();
    // // auto ray_data = sys->AllRayData;
    // sys->AllRayData.Print();
}

TEST(TowerDemo, NativeRunnerWithErrors)
{
    SimulationData sd;

    // Sun
    auto sun = make_ray_source<Sun>();
    sun->set_position(0.0, 0.0, 100.0);
    sd.add_ray_source(sun);

    // Absorber -- Flat
    auto absorber = make_element<SingleElement>();
    absorber->set_origin(0.0, 0.0, 10.0);
    absorber->set_aim_vector(0.0, 5.0, 0.0);
    absorber->set_zrot(0.0);
    absorber->compute_coordinate_rotations();
    absorber->set_surface(make_surface<Flat>()); // surface(nullptr)
    absorber->set_aperture(make_aperture<Rectangle>(2.0, 2.0));
    OpticalProperties *foptics = absorber->get_front_optical_properties();
    foptics->my_type = InteractionType::REFLECTION;
    foptics->reflectivity = 0.0;
    absorber->set_name("Absorber");

    // Make stage 1 -- second stage -- these can be added to SimulationData
    // in any order but should be numbered in the desired order
    auto st1 = make_stage(1);
    // Origin is initialized to zero but set it explicitly
    st1->set_origin(0.0, 0.0, 0.0);
    // Set aim vector so stage and global coordinates are identical
    st1->set_aim_vector(0.0, 0.0, 1.0);
    st1->set_zrot(0.0);
    st1->compute_coordinate_rotations();
    st1->add_element(absorber);
    // Optional -- to help the user identify things
    st1->set_name("Stage 1--Absorber");

    // Make stage 0 -- this will be the first stage if the runner uses stages
    auto st0 = make_stage(0);
    st0->set_origin(0.0, 0.0, 0.0);
    st0->set_aim_vector(0.0, 0.0, 1.0);
    st0->set_zrot(0.0);
    st0->compute_coordinate_rotations();
    st0->set_name("Stage 0--Reflectors");

    Vector3d rvec, svec, avec;
    Vector3d aim, pos;

    for (int i = -1; i < 2; ++i)
    {
        auto el = make_element<SingleElement>();
        foptics = el->get_front_optical_properties();
        foptics->reflectivity = 1.0;

        pos.set_values(5 * sin(i * PI / 2.0),
                       5 * cos(i * PI / 2.0),
                       0.0);
        el->set_origin(pos);
        vector_add(1.0, absorber->get_origin_global(),
                   -1.0, pos,
                   rvec);
        make_unit_vector(rvec);
        svec = sun->get_position();
        make_unit_vector(svec);
        vector_add(0.5, rvec, 0.5, svec, avec);

        vector_add(1.0, pos, 100.0, avec, aim);
        el->set_aim_vector(aim);

        // TODO: Set zrot as in python file?
        el->set_zrot(30.0 * i);
        // // Could also be set in radians
        // el->set_zrot_radians(PI / 3.0 * i);
        el->compute_coordinate_rotations();

        el->set_surface(make_surface<Flat>());
        el->set_aperture(make_aperture<Rectangle>(1.0, 1.95));

        std::stringstream name;
        name << "Reflector " << i;
        el->set_name(name.str());

        st0->add_element(el);
    }

    // Stages must have all elements present before adding to SimulationData!!
    sd.add_stage(st0);
    sd.add_stage(st1);

    // Set parameters
    SimulationParameters &params = sd.get_simulation_parameters();
    params.number_of_rays = 100000;
    // params.number_of_rays = 10; // Above takes too long
    params.max_number_of_rays = params.number_of_rays * 100;
    params.include_optical_errors = true;
    params.include_sun_shape_errors = true;
    params.seed = 12345;

    NativeRunner runner;
    RunnerStatus sts = runner.initialize();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    // Setup runs but is not complete
    sts = runner.setup_simulation(&sd);
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);
    // Run simulation runs but returns RunnerStatus::ERROR
    sts = runner.run_simulation();
    EXPECT_EQ(sts, RunnerStatus::SUCCESS);

    // TODO: Do some post processing tests here

    // const TSystem *sys = runner.get_system();
    // // auto ray_data = sys->AllRayData;
    // sys->AllRayData.Print();
}

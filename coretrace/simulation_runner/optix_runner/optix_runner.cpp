
#include "simulation_runner/optix_runner/optix_runner.hpp"
#include "simulation_data/simulation_parameters.hpp"
#include "simulation_data/simulation_data.hpp"
#include "simulation_data/simulation_data_export.hpp"

using SolTrace::Runner::RunnerStatus;
using SolTrace::Runner::SimulationRunner;

using SolTrace::Result::SimulationResult;

OptixRunner::OptixRunner() : SimulationRunner(),
                             m_simdata(nullptr),
                             m_sys(10000) {}

OptixRunner::~OptixRunner()
{
}

RunnerStatus OptixRunner::initialize()
{
    // add elements to sys using data structure from SimulationData

    // set number of rays

    // set sun vector, and other sun properties

    return RunnerStatus::SUCCESS;
}

RunnerStatus OptixRunner::setup_simulation(const SimulationData *data)
{

    RunnerStatus sts;

    // this->simdata = data;

    this->setup_parameters(data);
    this->setup_sun(data);
    sts = this->setup_elements(data);

    m_sys.initialize();

    // std::cout << "Number of stages: " << this->tsys.StageList.size()
    //           << std::endl;

    return sts;
}

RunnerStatus OptixRunner::setup_parameters(const SimulationData *data)
{
    // Get Parameter data
    // TODO: Check that these parameters are used as expected
    const SimulationParameters &sim_params = data->get_simulation_parameters();
    m_sys.set_sun_points(sim_params.max_number_of_rays);

    // this->tsys.sim_errors_sunshape = sim_params.include_sun_shape_errors;
    // this->tsys.sim_errors_optical = sim_params.include_optical_errors;
    // this->tsys.sim_raycount = sim_params.number_of_rays;
    // this->tsys.sim_raymax = sim_params.max_number_of_rays;
    // this->tsys.seed = sim_params.seed;
    return RunnerStatus::SUCCESS;
}

RunnerStatus OptixRunner::setup_sun(const SimulationData *data)
{
    // Get RaySource data (this runner assumes there is only the Sun)
    assert(data->get_number_of_ray_sources() == 1);

    m_sys.set_sun_vector(ToVec3d(data->get_ray_source()->get_position()));

    //  TODO: sun angle and sun models

    return RunnerStatus::SUCCESS;
}

RunnerStatus OptixRunner::setup_elements(const SimulationData *data)
{
    for (auto iter = data->get_const_iterator();
         !data->is_at_end(iter);
         ++iter)
    {
        element_ptr el = iter->second;
        if (el->is_enabled())
        {
            auto optix_el = std::make_shared<OptixCSP::CspElement>();
            Vector3d origin = el->get_origin_global();
            OptixCSP::Vec3d origin_vec(origin[0], origin[1], origin[2]);
            optix_el->set_origin(ToVec3d(origin));
            optix_el->set_aim_point(ToVec3d(el->get_aim_vector_global()));

            // TODO: check zrot, radiance or degree here?

            std::cout << "adding elements " << el->get_name() << std::endl;
            std::cout << "Origin: " << origin[0] << ", " << origin[1] << ", " << origin[2] << std::endl;

            if (el->get_surface() != nullptr)
            {
                std::cout << "surface type: " << el->get_surface()->get_type() << std::endl;

                switch (el->get_surface()->get_type())
                {
                case SurfaceType::FLAT:
                {
                    auto surface = std::make_shared<OptixCSP::SurfaceFlat>();
                    optix_el->set_surface(surface);

                    break;
                }
                case SurfaceType::PARABOLA:
                {
                    auto surface = std::make_shared<OptixCSP::SurfaceParabolic>();
                    surface->set_curvature(0.02, 0.04);
                    optix_el->set_surface(surface);

                    break;
                }
                case SurfaceType::CYLINDER:
                {
                    auto surface = std::make_shared<OptixCSP::SurfaceCylinder>();
                    surface->set_half_height(2.);
                    surface->set_radius(1.);
                    optix_el->set_surface(surface);

                    break;
                }
                default:
                    std::cerr << "Unsupported surface type in OptixCSP" << std::endl;
                    break;
                }

                auto soltrace_aperture_type = el->get_aperture()->get_type();

                switch (soltrace_aperture_type)
                {

                case ApertureType::RECTANGLE:
                {

                    // TODO: still a placeholder now

                    auto aperture = std::make_shared<OptixCSP::ApertureRectangle>(1.5, 1.5);
                    optix_el->set_aperture(aperture);
                    break;
                }

                default:
                    std::cerr << "Unsupported aperture type in OptixCSP" << std::endl;
                    break;
                }

                m_sys.add_element(optix_el);
            }
            std::cout << "=====================================================" << std::endl;
        }
    }
    return RunnerStatus::SUCCESS;
}

RunnerStatus OptixRunner::update_simulation(const SimulationData *data)
{
    this->setup_simulation(data);
    // TODO: Implement this
    return RunnerStatus::SUCCESS;
}

RunnerStatus OptixRunner::run_simulation()
{
    m_sys.run();
    m_sys.write_hp_output("output.txt");
    return RunnerStatus::SUCCESS;
}

RunnerStatus OptixRunner::status_simulation(double *progress)
{
    // TODO: Implement this
    return RunnerStatus::SUCCESS;
}

RunnerStatus OptixRunner::cancel_simulation()
{
    // TODO: Implement this
    return RunnerStatus::SUCCESS;
}

RunnerStatus OptixRunner::report_simulation(SimulationResult *result,
                                            int level)
{
    // for different levels of reporting, populate result accordingly
    //
    return RunnerStatus::SUCCESS;
}

OptixCSP::Vec3d OptixRunner::ToVec3d(Vector3d v)
{

    OptixCSP::Vec3d vec(v[0], v[1], v[2]);
    return vec;
}

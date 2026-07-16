#include "simulation_runner/optix_runner/optix_runner.hpp"
#include "simulation_data/simulation_data_export.hpp"

#include <iostream>
#include <optical_properties.hpp>
#include <sstream>
#include <stdexcept>

using SolTrace::Runner::RunnerStatus;
using SolTrace::Runner::SimulationRunner;

using SolTrace::Result::SimulationResult;

using SolTrace::Data::optics_id;

OptixRunner::OptixRunner() : SimulationRunner(),
                             m_simdata(nullptr),
                             m_sys() {}

void OptixRunner::set_verbose(bool verbose)
{
    m_sys.set_verbose(verbose);
}

void OptixRunner::print_timing() const
{
    m_sys.print_timing();
}

void OptixRunner::set_max_ray_depth(uint_fast64_t depth)
{
    m_sys.set_max_ray_depth(depth);
}

void OptixRunner::set_batch_size(uint_fast64_t batch_size)
{
    m_sys.set_batch_size(batch_size);
}

uint_fast64_t OptixRunner::get_batch_size() const
{
    return m_sys.get_batch_size();
}

void OptixRunner::set_trim_excess_rays(bool trim)
{
    m_sys.set_trim_excess_rays(trim);
}

bool OptixRunner::get_trim_excess_rays() const
{
    return m_sys.get_trim_excess_rays();
}

uint64_t OptixRunner::get_N_run_iterations() const
{
    return m_sys.get_N_run_iterations();
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

    // Reset
    this->m_sys.reset();

    RunnerStatus sts;

    sts = this->setup_parameters(data);
    if (sts != RunnerStatus::SUCCESS)
        return sts;
    sts = this->setup_sun(data);
    if (sts != RunnerStatus::SUCCESS)
        return sts;
    sts = this->setup_elements(data);
    if (sts != RunnerStatus::SUCCESS)
        return sts;

    m_sys.initialize();

    // std::cout << "Number of stages: " << this->tsys.StageList.size()
    //           << std::endl;

    return sts;
}

RunnerStatus OptixRunner::setup_parameters(const SimulationData *data)
{
    // Get Parameter data
    const SimulationParameters &sim_params = data->get_simulation_parameters();

    m_sys.set_number_of_rays(sim_params.number_of_rays, sim_params.max_number_of_rays);
    m_sys.set_seed(static_cast<uint64_t>(sim_params.seed));

    m_sys.set_optical_errors(sim_params.include_optical_errors);
    m_sys.set_sun_shape_errors(sim_params.include_sun_shape_errors);

    return RunnerStatus::SUCCESS;
}

RunnerStatus OptixRunner::setup_sun(const SimulationData *data)
{
    // Get RaySource data (this runner assumes there is only the Sun)
    assert(data->get_number_of_ray_sources() == 1);

    // Verify that the only ray source is a Sun
    auto src = data->get_ray_source();
    auto sun = std::dynamic_pointer_cast<SolTrace::Data::Sun>(src);
    if (!sun)
    {
        // Needs to be a Sun
        return RunnerStatus::ERROR;
    }
    m_sys.set_sun(sun.get());

    // Check if sun shape is assigned
    if (data->get_simulation_parameters().include_sun_shape_errors)
    {
        const SolTrace::Data::SunShape shape = sun->get_shape();
        bool is_supported = false;
        for (auto supported_shape : OptixCSP::kSupportedSunshapes)
        {
            if (shape == supported_shape)
            {
                is_supported = true;
                break;
            }
        }
        if (!is_supported)
        {
            return RunnerStatus::ERROR;
        }
    }

    // Warn if Halton sampling is used with more rays than uint32_t can index,
    // since the Halton sequence index is truncated to 32 bits causing repeated positions.
    if (sun->get_gen_type() == SolTrace::Data::GenType::HALTON &&
        data->get_simulation_parameters().max_number_of_rays > static_cast<uint_fast64_t>(std::numeric_limits<uint32_t>::max()))
    {
        std::cerr << "Warning: max_number_of_rays exceeds 32-bit unsigned int maximum ("
                  << std::numeric_limits<uint32_t>::max()
                  << ") with Halton ray generation. Halton sequence positions will repeat after index "
                  << std::numeric_limits<uint32_t>::max() << "." << std::endl;
    }

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
            // Skip if element is not a single (i.e. stage, composite)
            if (el->is_single() == false)
                continue;

            auto optix_el = std::make_shared<OptixCSP::CspElement>();
            auto origin = el->get_origin_global();
            auto ap = el->get_aim_vector_global();
            // OptixCSP::Vec3d origin_vec(origin.x, origin.y, origin.z);
            optix_el->set_origin(ToVec3d(origin));
            optix_el->set_aim_point(ToVec3d(ap));
            optix_el->set_rotation_matrix(ToMatrix33d(el->get_local_to_global()));

            // Safely narrow element id to int32_t
            const auto id = el->get_id(); // int
            if (id < std::numeric_limits<int32_t>::min() || id > std::numeric_limits<int32_t>::max())
            {
                throw std::overflow_error("Element id out of int32_t range");
            }
            optix_el->set_id(static_cast<int32_t>(id));

            // Add optical properties
            auto opt_set = el->get_optical_property_set();

            if (opt_set == nullptr)
                throw std::runtime_error("Element has invalid optical property set.");

            DistributionType front_dist;
            double front_slope, front_spec;
            opt_set->get_errors(OpticalSide::Front, front_dist, front_slope, front_spec);

            OptixCSP::OpticalDistribution front_dist_optix = this->to_optical_distribution(front_dist);
            optix_el->set_optics_front(opt_set->get_interaction_type() == InteractionType::REFRACTION, 
                opt_set->get_reflectivity(OpticalSide::Front), opt_set->get_transmissivity(OpticalSide::Front), 
                front_slope, front_spec, front_dist_optix);

            DistributionType back_dist;
            double back_slope, back_spec;
            opt_set->get_errors(OpticalSide::Back, back_dist, back_slope, back_spec);

            OptixCSP::OpticalDistribution back_dist_optix = this->to_optical_distribution(back_dist);
            optix_el->set_optics_back(opt_set->get_interaction_type() == InteractionType::REFRACTION,
                opt_set->get_reflectivity(OpticalSide::Back), opt_set->get_transmissivity(OpticalSide::Back),
                back_slope, back_spec, back_dist_optix);

            if (m_sys.is_verbose())
            {
                std::cout << "adding elements " << el->get_name() << std::endl;
                auto origin = el->get_origin_global();
                std::cout << "Origin: (" << origin.x << ", " << origin.y << ", " << origin.z << ")" << std::endl;
            }

            if (el->get_surface() == nullptr)
            {
                throw std::runtime_error("Element must be assigned a surface.");
            }

            if (el->get_aperture() == nullptr)
            {
                throw std::runtime_error("Element must be assigned an aperture.");
            }

            if (m_sys.is_verbose())
            {
                std::cout << "surface type: " << el->get_surface()->get_type() << std::endl;
            }

            switch (el->get_surface()->get_type())
            {
            case SurfaceType::FLAT:
            {
                auto surface = std::make_shared<OptixCSP::SurfaceFlat>();
                assert(surface != nullptr);
                optix_el->set_surface(surface);

                break;
            }
            case SurfaceType::PARABOLA:
            {
                auto el_surface = std::dynamic_pointer_cast<Parabola>(el->get_surface());
                assert(el_surface != nullptr);
                double fx = el_surface->focal_length_x;
                double fy = el_surface->focal_length_y;

                double cx = 1. / (2. * fx);
                double cy = 1. / (2. * fy);

                auto optix_surface = std::make_shared<OptixCSP::SurfaceParabolic>();
                optix_surface->set_curvature(cx, cy);
                optix_el->set_surface(optix_surface);

                break;
            }
            case SurfaceType::CYLINDER:
            {
                auto el_surface = std::dynamic_pointer_cast<Cylinder>(el->get_surface());
                assert(el_surface != nullptr);
                auto el_aperture = std::dynamic_pointer_cast<Rectangle>(el->get_aperture());
                // assert(el_aperture != nullptr);
                if (el_aperture == nullptr)
                {
                    throw std::runtime_error("Cylinder surface type must have rectangular aperture.");
                }

                if (fabs(0.5 * el_aperture->x_length() - el_surface->radius) > 1e-6)
                {
                    throw std::runtime_error("Rectangle aperture has incorrect dimension for cylinder surface.");
                }

                auto surface = std::make_shared<OptixCSP::SurfaceCylinder>();
                surface->set_half_height(0.5 * el_aperture->y_length());
                surface->set_radius(el_surface->radius);
                optix_el->set_surface(surface);

                break;
            }
            default:
                // std::cerr << "Unsupported surface type in OptixCSP" << std::endl;
                throw std::runtime_error("Unsupported surface type in OptixCSP");
                break;
            }

            auto soltrace_aperture_type = el->get_aperture()->get_type();

            switch (soltrace_aperture_type)
            {
            case ApertureType::RECTANGLE:
            {
                auto el_aperture = std::dynamic_pointer_cast<Rectangle>(el->get_aperture());
                assert(el_aperture != nullptr);
                // TODO: account for x and y coord?
                // auto aperture = std::make_shared<OptixCSP::ApertureRectangle>(el_aperture->x_length(),
                // 							      el_aperture->y_length());
                auto aperture = std::make_shared<OptixCSP::ApertureRectangle>(el_aperture->x_length(),
                                                                              el_aperture->y_length(),
                                                                              el_aperture->x_coord(),
                                                                              el_aperture->y_coord());
                optix_el->set_aperture(aperture);
                break;
            }
            case ApertureType::ANNULUS:
            {
                auto el_aperture = std::dynamic_pointer_cast<Annulus>(el->get_aperture());
                assert(el_aperture != nullptr);
                auto aperture = std::make_shared<OptixCSP::ApertureAnnulus>(el_aperture->inner_radius, el_aperture->outer_radius, el_aperture->arc_angle * D2R);
                optix_el->set_aperture(aperture);
                break;
            }
            case ApertureType::CIRCLE:
            {
                auto el_aperture = std::dynamic_pointer_cast<Circle>(el->get_aperture());
                assert(el_aperture != nullptr);
                auto aperture = std::make_shared<OptixCSP::ApertureCircle>(0.5 * el_aperture->diameter);
                optix_el->set_aperture(aperture);
                break;
            }
            case ApertureType::EQUILATERAL_TRIANGLE:
            {
                auto el_aperture = std::dynamic_pointer_cast<EquilateralTriangle>(el->get_aperture());
                assert(el_aperture != nullptr);
                double r = 0.5 * el_aperture->circumscribe_diameter;

                OptixCSP::Vec3d p0(-sqrt(0.75) * r, -0.5 * r, 0.0);
                OptixCSP::Vec3d p1(sqrt(0.75) * r, -0.5 * r, 0.0);
                OptixCSP::Vec3d p2(0.0, r, 0.0);

                auto aperture = std::make_shared<OptixCSP::ApertureTriangle>(p0, p1, p2);
                optix_el->set_aperture(aperture);

                break;
            }
            case ApertureType::IRREGULAR_TRIANGLE:
            {
                auto el_aperture = std::dynamic_pointer_cast<IrregularTriangle>(el->get_aperture());
                assert(el_aperture != nullptr);


                OptixCSP::Vec3d p0(el_aperture->x1, el_aperture->y1, 0.0);
                OptixCSP::Vec3d p1(el_aperture->x2, el_aperture->y2, 0.0);
                OptixCSP::Vec3d p2(el_aperture->x3, el_aperture->y3, 0.0);

                // Ensure CCW winding (right-hand rule) required by ApertureTriangle
                // Aperture always lies in x-y-plane with positive z-axis corresponding to
                // the front of the geometric element so CCW winding always gives the front
                // of the element.
                const double signed_area = (p1[0] - p0[0]) * (p2[1] - p0[1]) - (p1[1] - p0[1]) * (p2[0] - p0[0]);
                if (signed_area < 0.0)
                    std::swap(p1, p2);

                auto aperture = std::make_shared<OptixCSP::ApertureTriangle>(p0, p1, p2);
                optix_el->set_aperture(aperture);

                break;
            }
            case ApertureType::IRREGULAR_QUADRILATERAL:
            {
                auto el_aperture = std::dynamic_pointer_cast<IrregularQuadrilateral>(el->get_aperture());
                assert(el_aperture != nullptr);

                OptixCSP::Vec3d p0(el_aperture->x1, el_aperture->y1, 0.0);
                OptixCSP::Vec3d p1(el_aperture->x2, el_aperture->y2, 0.0);
                OptixCSP::Vec3d p2(el_aperture->x3, el_aperture->y3, 0.0);
                OptixCSP::Vec3d p3(el_aperture->x4, el_aperture->y4, 0.0);

                // Ensure CCW winding using the shoelace signed-area formula.
                // For a simple (non-self-intersecting) quad the sign of 2*A tells
                // the winding without decomposing into triangles.
                const double signed_area2 = (p0[0] * p1[1] - p1[0] * p0[1])
                                          + (p1[0] * p2[1] - p2[0] * p1[1])
                                          + (p2[0] * p3[1] - p3[0] * p2[1])
                                          + (p3[0] * p0[1] - p0[0] * p3[1]);
                if (signed_area2 < 0.0)
                    std::swap(p1, p3); // reverse winding, keeping p0 and p2 fixed

                auto aperture = std::make_shared<OptixCSP::ApertureQuadrilateral>(p0, p1, p2, p3);
                optix_el->set_aperture(aperture);

                break;
            }
            case ApertureType::HEXAGON:
            {
                auto el_aperture = std::dynamic_pointer_cast<Hexagon>(el->get_aperture());
                assert(el_aperture != nullptr);
                auto aperture = std::make_shared<OptixCSP::ApertureHexagon>(el_aperture->radius_circumscribed_circle());
                optix_el->set_aperture(aperture);
                break;
            }
            default:
                // std::cerr << "Unsupported aperture type in OptixCSP" << std::endl;
                throw std::runtime_error("Unsupported aperture type in OptixRunner");
                break;
            }

            double xmin, xmax, ymin, ymax, zmin, zmax;
            el->get_aperture()->bounding_box(xmin, xmax, ymin, ymax);
            el->get_surface()->bounding_box(xmin, xmax, ymin, ymax, zmin, zmax);
            OptixCSP::Vec3d upper(xmax, ymax, zmax);
            OptixCSP::Vec3d lower(xmin, ymin, zmin);
            optix_el->set_bounding_box_local(lower, upper);

            if (m_sys.is_verbose())
            {
                std::cout << "Bounding Box Upper: " << optix_el->get_upper_bounding_box() << std::endl;
                std::cout << "Bounding Box Lower: " << optix_el->get_lower_bounding_box() << std::endl;
            }

            m_sys.add_element(optix_el);

            if (m_sys.is_verbose())
            {
                std::cout << "=====================================================" << std::endl;
            }
        }
    }
    return RunnerStatus::SUCCESS;
}

RunnerStatus OptixRunner::update_simulation(const SimulationData *data)
{
    return this->setup_simulation(data);
    // TODO: Implement this in a less lazy manner...
}

RunnerStatus OptixRunner::run_simulation()
{
    return run_simulation_core();
}

RunnerStatus OptixRunner::run_simulation_core()
{

    m_sys.run();

    return RunnerStatus::SUCCESS;
}

RunnerStatus OptixRunner::status_simulation(double *progress)
{
    // TODO: Implement this
    return RunnerStatus::SUCCESS;
}

// Temporary function to get hit points
RunnerStatus OptixRunner::get_hp_output(std::vector<float4> &hp_vec, std::vector<uint_fast64_t> &raynumber_vec,
                                        std::vector<int32_t> &element_id_vec)
{
    // for different levels of reporting, populate result accordingly
    //
    std::vector<uint8_t> hit_type_vec;
    m_sys.get_hp_output(hp_vec, raynumber_vec, element_id_vec, hit_type_vec);
    return RunnerStatus::SUCCESS;
}

SolTrace::Result::RayEvent hit_type_to_ray_event(OptixCSP::HitType hit_type)
{
    if (hit_type == OptixCSP::HitType::HIT_UNASSIGNED || hit_type == OptixCSP::HitType::HIT_UNKNOWN)
        return SolTrace::Result::RayEvent::UNKNOWN;

    return static_cast<SolTrace::Result::RayEvent>(hit_type);
}

RunnerStatus OptixRunner::report_simulation(SimulationResult *result,
                                            int level)
{
    // Declare results
    RunnerStatus retval = RunnerStatus::SUCCESS;
    std::map<unsigned int, SolTrace::Result::ray_record_ptr> ray_records;
    std::map<unsigned int, SolTrace::Result::ray_record_ptr>::iterator iter;

    // TODO: This should be redone without using these vectors and just using the
    // internal hit record vector
    // Get results from optixcsp
    std::vector<float4> hp_vec;
    std::vector<uint_fast64_t> raynumber_vec;
    std::vector<int32_t> element_id_vec;
    std::vector<uint8_t> hit_type_vec;
    m_sys.get_hp_output(hp_vec, raynumber_vec, element_id_vec, hit_type_vec);

    // Check sizes
    if (!(hp_vec.size() == raynumber_vec.size() && raynumber_vec.size() == element_id_vec.size() && element_id_vec.size() == hit_type_vec.size()))
    {
        return RunnerStatus::ERROR;
    }

    // Loop through data, populating ray records
    // Assumes ray data is grouped serially
    size_t ndata = hp_vec.size();
    // uint_fast64_t raynum_prev = -1;
    uint_fast64_t raynum = 0;
    SolTrace::Result::ray_record_ptr rec = nullptr;
    SolTrace::Result::interaction_ptr intr = nullptr;
    for (size_t ii = 0; ii < ndata; ++ii)
    {
        // Collect results for record
        raynum = raynumber_vec[ii];
        glm::dvec3 pos(hp_vec[ii].y, hp_vec[ii].z, hp_vec[ii].w); // x is depth
        glm::dvec3 cos(0.0);                                      // TODO: calculate directions
        int32_t element_id = element_id_vec[ii];
        uint8_t hit_type = hit_type_vec[ii];
        SolTrace::Result::RayEvent rev = hit_type_to_ray_event(static_cast<OptixCSP::HitType>(hit_type));

        // Make new ray record if necessary
        iter = ray_records.find(raynum);
        if (iter == ray_records.end())
        {
            rec = SolTrace::Result::make_ray_record(raynum);
            result->add_ray_record(rec);
            ray_records[raynum] = rec;
            assert(rev == SolTrace::Result::RayEvent::CREATE);
        }
        else
        {
            rec = iter->second;
        }

        // Make interaction record
        intr = SolTrace::Result::make_interaction_record(element_id, rev, pos, cos);
        rec->add_interaction_record(intr);
    }

    // Attach other results
    result->set_sun_ray_count(this->get_N_sun_rays());
    result->set_sun_dimensions(std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::quiet_NaN());
    result->set_sun_A_box(this->get_sun_plane_area());

    return RunnerStatus::SUCCESS;
}

RunnerStatus OptixRunner::cancel_simulation()
{
    RunnerStatus sts = RunnerStatus::ERROR;

    // TODO: Implement actual cancel

    return sts;
}

OptixCSP::Vec3d OptixRunner::ToVec3d(glm::dvec3 v)
{
    OptixCSP::Vec3d vec(v.x, v.y, v.z);
    return vec;
}

OptixCSP::Matrix33d OptixRunner::ToMatrix33d(const glm::dmat3 &mat)
{
    return OptixCSP::Matrix33d(
        mat[0][0], mat[1][0], mat[2][0],
        mat[0][1], mat[1][1], mat[2][1],
        mat[0][2], mat[1][2], mat[2][2]);
}

OptixCSP::OpticalDistribution OptixRunner::to_optical_distribution(SolTrace::Data::DistributionType dt)
{
    OptixCSP::OpticalDistribution od;
    if (dt == SolTrace::Data::DistributionType::NONE)
        od = OptixCSP::OpticalDistribution::OPT_NONE;
    else if (dt == SolTrace::Data::DistributionType::GAUSSIAN)
        od = OptixCSP::OpticalDistribution::OPT_GAUSSIAN;
    else if (dt == SolTrace::Data::DistributionType::PILLBOX)
        od = OptixCSP::OpticalDistribution::OPT_PILLBOX;
    else
    {
        std::stringstream ss;
        ss << "Unimplemented error distribution: "
           << SolTrace::Data::distribution_string(dt)
           << std::endl;

        throw std::invalid_argument(ss.str());
    }
    return od;
}

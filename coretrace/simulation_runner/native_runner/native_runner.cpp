
#include "native_runner.hpp"

#include <chrono>
#include <exception>
#include <map>
#include <mutex>
#include <thread>

// SimulationData headers
#include "composite_element.hpp"
#include "element.hpp"
#include "simulation_parameters.hpp"
#include "simulation_data.hpp"
#include "simulation_data_export.hpp"

// NativeRunner headers
#include "native_runner_types.hpp"
#include "trace.hpp"

namespace SolTrace::NativeRunner
{

    NativeRunner::NativeRunner() : SimulationRunner(),
                                   as_power_tower(false),
                                   number_of_threads(1)
    {
    }

    NativeRunner::~NativeRunner()
    {
    }

    RunnerStatus NativeRunner::initialize()
    {
        return RunnerStatus::SUCCESS;
    }

    RunnerStatus NativeRunner::setup_simulation(const SimulationData *data)
    {

        RunnerStatus sts;

        sts = this->setup_parameters(data);

        if (sts == RunnerStatus::SUCCESS)
            sts = this->setup_sun(data);

        if (sts == RunnerStatus::SUCCESS)
            sts = this->setup_elements(data);

        return sts;
    }

    RunnerStatus NativeRunner::setup_parameters(const SimulationData *data)
    {
        // Get Parameter data
        // TODO: Check that these parameters are used as expected
        const SimulationParameters &sim_params = data->get_simulation_parameters();
        this->tsys.sim_errors_sunshape = sim_params.include_sun_shape_errors;
        this->tsys.sim_errors_optical = sim_params.include_optical_errors;
        this->tsys.sim_raycount = sim_params.number_of_rays;
        this->tsys.sim_raymax = sim_params.max_number_of_rays;
        this->tsys.seed = sim_params.seed;
        return RunnerStatus::SUCCESS;
    }

    RunnerStatus NativeRunner::setup_sun(const SimulationData *data)
    {
        // TODO: This should throw an error...
        // Get RaySource data (this runner assumes there is only the Sun)
        assert(data->get_number_of_ray_sources() == 1);

        ray_source_ptr sun = data->get_ray_source();
        vector_copy(this->tsys.Sun.Origin, sun->get_position());
        this->tsys.Sun.ShapeIndex = sun->get_shape();

        // Set sunshape data
        switch (sun->get_shape())
        {
        case SunShape::GAUSSIAN:
            this->tsys.Sun.Sigma = sun->get_sigma();
            break;
        case SunShape::PILLBOX:
            this->tsys.Sun.Sigma = sun->get_half_width();
            break;
        case SunShape::LIMBDARKENED:
            this->tsys.Sun.MaxAngle = 4.65; // [mrad]
            this->tsys.Sun.MaxIntensity = 1.0;
            break;
        case SunShape::BUIE_CSR:
        {
            this->tsys.Sun.MaxAngle = 43.6; // [mrad]
            this->tsys.Sun.MaxIntensity = 1.0;
            double kappa, gamma;
            sun->calculate_buie_parameters(kappa, gamma);
            this->tsys.Sun.buie_kappa = kappa;
            this->tsys.Sun.buie_gamma = gamma;
            break;
        }
        case SunShape::USER_DEFINED:
        {
            std::vector<double> angle, intensity;
            sun->get_user_data(angle, intensity);
            int npoints = angle.size();

            // Set user data
            this->tsys.Sun.MaxAngle = 0;
            this->tsys.Sun.MaxIntensity = 0;

            this->tsys.Sun.SunShapeAngle.resize(2 * npoints - 1);
            this->tsys.Sun.SunShapeIntensity.resize(2 * npoints - 1);

            for (int i = 0; i < npoints; i++)
            {
                this->tsys.Sun.SunShapeAngle[npoints + i - 1] = angle[i];
                this->tsys.Sun.SunShapeIntensity[npoints + i - 1] = intensity[i];

                if (angle[i] > this->tsys.Sun.MaxAngle)
                    this->tsys.Sun.MaxAngle = angle[i];
                if (intensity[i] > this->tsys.Sun.MaxIntensity)
                    this->tsys.Sun.MaxIntensity = intensity[i];
            }
            // fill negative angle side of array -> I don't think we need this.
            // for (int i = 0; i < npoints - 1; i++)
            //{
            //    this->tsys.Sun.SunShapeAngle[i] = -angle[npoints - i - 1];
            //    this->tsys.Sun.SunShapeIntensity[i] = intensity[npoints - i - 1];
            //}
            break;
        }
        default:
            // TODO: add error
            break;
        }

        return RunnerStatus::SUCCESS;
    }

    RunnerStatus NativeRunner::setup_elements(const SimulationData *data)
    {
        // TODO: Improve error messages from this function.

        RunnerStatus sts = RunnerStatus::SUCCESS;
        auto my_map = std::map<int_fast64_t, tstage_ptr>();
        // int_fast64_t current_stage_id = -1;
        tstage_ptr current_stage = nullptr;
        int_fast64_t element_number = 1;
        bool element_found_before_stage = false;

        for (auto iter = data->get_const_iterator();
             !data->is_at_end(iter);
             ++iter)
        {
            element_ptr el = iter->second;
            if (el->is_enabled() && el->is_stage())
            {
                tstage_ptr stage = make_tstage(el, this->eparams);
                auto retval = my_map.insert(
                    std::make_pair(el->get_stage(), stage));

                // current_stage_id = stage->stage_id;

                // std::cout << "Created stage " << el->get_stage()
                //           << " with " << stage->ElementList.size() << " elements"
                //           << std::endl;

                if (retval.second == false)
                {
                    // TODO: Duplicate stage numbers. Need to make an error
                    // message.
                    sts = RunnerStatus::ERROR;
                }

                current_stage = stage;
                element_number = 1;
            }
            else if (el->is_enabled() && el->is_single())
            {
                if (current_stage == nullptr)
                {
                    // throw std::runtime_error("No stage to add element to");
                    element_found_before_stage = true;
                    continue;
                }
                else if (el->get_stage() != current_stage->stage_id)
                {
                    throw std::runtime_error(
                        "Element does not match current stage");
                }

                telement_ptr elem = make_telement(iter->second,
                                                  element_number,
                                                  this->eparams);
                ++element_number;
                current_stage->ElementList.push_back(elem);
            }
        }

        if (my_map.size() != 0 && element_found_before_stage)
        {
            throw std::runtime_error("Element found without a stage");
        }

        if (my_map.size() == 0)
        {
            // No stage elements found in the passed in data. However,
            // the runner requires stages. So make a single stage
            // and put everything there. Note that the coordinates are
            // set to correspond to global coordinates. This is necessary
            // so that the element coordinate setup in make_element are
            // correct.
            int_fast64_t element_number = 1;
            auto stage = make_tstage(this->eparams);
            stage->ElementList.reserve(data->get_number_of_elements());
            for (auto iter = data->get_const_iterator();
                 !data->is_at_end(iter);
                 ++iter)
            {
                element_ptr el = iter->second;
                if (el->is_enabled() && el->is_single())
                {
                    telement_ptr tel = make_telement(el,
                                                     element_number,
                                                     this->eparams);
                    stage->ElementList.push_back(tel);
                    ++element_number;
                }
            }
            my_map.insert(std::make_pair(0, stage));
        }

        // std::map (according to the documentation) is automatically
        // ordered by the keys so inserting into a map will sort the stages
        // and we can just transfer the pointers, in order, to the StageList
        // simply by pulling them out of the map.
        int_fast64_t last_stage_id = -1;
        for (auto iter = my_map.cbegin();
             iter != my_map.cend();
             ++iter)
        {
            assert(last_stage_id < iter->first);
            last_stage_id = iter->first;
            this->tsys.StageList.push_back(iter->second);
        }

        if (sts == RunnerStatus::SUCCESS)
        {
            // std::cout << "Setting ZAperture..." << std::endl;
            // Compute and set ZAperture field in each element
            bool success = set_aperture_planes(&this->tsys);
            sts = success ? RunnerStatus::SUCCESS : RunnerStatus::ERROR;
        }

        return sts;
    }

    RunnerStatus NativeRunner::update_simulation(const SimulationData *data)
    {
        // TODO: Do a more efficient implementation of this?
        this->tsys.ClearAll();
        this->setup_simulation(data);
        return RunnerStatus::SUCCESS;
    }

    RunnerStatus NativeRunner::run_simulation()
    {
        RunnerStatus sts = trace_native(
            &this->tsys,
            this->tsys.seed,
            this->tsys.sim_raycount,
            this->tsys.sim_raymax,
            this->tsys.sim_errors_sunshape,
            this->tsys.sim_errors_optical,
            this->as_power_tower);

        {
            // Hack to match up current state with return type
            std::lock_guard<std::mutex> lk(this->tsys.state_mutex);
            this->tsys.current_state = sts;
        }

        return sts;
    }

    RunnerStatus NativeRunner::status_simulation(double *progress)
    {
        RunnerStatus sts = RunnerStatus::ERROR;
        // Create isolated scope for lock guard
        {
            std::lock_guard<std::mutex> lk(this->tsys.state_mutex);
            sts = this->tsys.current_state;
            if (progress != nullptr)
            {
                *progress = this->tsys.progress;
            }
        }
        return sts;
    }

    RunnerStatus NativeRunner::cancel_simulation()
    {
        RunnerStatus sts = RunnerStatus::ERROR;

        // Create isolated scope for the lock
        {
            std::lock_guard<std::mutex> my_lock(this->tsys.state_mutex);
            this->tsys.cancel = true;
        }

        int count = 0;
        while (true)
        {
            ++count;
            std::this_thread::sleep_for(std::chrono::milliseconds(100));

            // Create isolated scope for the lock
            {
                std::lock_guard<std::mutex> my_lock(this->tsys.state_mutex);
                if (this->tsys.current_state != RunnerStatus::RUNNING)
                {
                    sts = this->tsys.current_state;
                    break;
                }
            }

            if (count > 30)
            {
                sts = RunnerStatus::TIMEOUT;
                break;
            }
        }

        return sts;
    }

    RunnerStatus NativeRunner::report_simulation(SolTrace::Result::SimulationResult *result,
                                                 int level)
    {
        RunnerStatus retval = RunnerStatus::SUCCESS;

        const TSystem *sys = this->get_system();
        // const TRayData ray_data = sys->AllRayData;
        const TRayData ray_data = sys->RayData;
        std::map<unsigned int, SolTrace::Result::ray_record_ptr> ray_records;
        std::map<unsigned int, SolTrace::Result::ray_record_ptr>::iterator iter;
        size_t ndata = ray_data.Count();

        bool sts;
        Vector3d point, cosines;
        int element;
        int stage;
        unsigned int raynum;

        telement_ptr el = nullptr;
        element_id elid;
        SolTrace::Result::ray_record_ptr rec = nullptr;
        SolTrace::Result::interaction_ptr intr = nullptr;
        SolTrace::Result::RayEvent rev;

        for (size_t ii = 0; ii < ndata; ++ii)
        {
            sts = ray_data.Query(ii,
                                 point.data,
                                 cosines.data,
                                 &element,
                                 &stage,
                                 &raynum,
                                 &rev);

            if (!sts)
            {
                retval = RunnerStatus::ERROR;
                break;
            }

            // std::cout << "ii: " << ii
            //           << "\npoint: " << point
            //           << "\ndirection: " << cosines
            //           << "\nelement: " << element
            //           << "\nstage: " << stage
            //           << "\nraynum: " << raynum
            //           << "\nevent: " << ray_event_string(rev)
            //           << std::endl;

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

            if (element > 0)
            {
                el = sys->StageList[stage - 1]->ElementList[element - 1];
                elid = el->sim_data_id;
            }
            else
            {
                elid = element;
            }

            intr = make_interaction_record(elid, rev, point, cosines);
            rec->add_interaction_record(intr);
        }

        return RunnerStatus::SUCCESS;
    }

    bool NativeRunner::set_aperture_planes(TSystem *tsys)
    {
        bool retval;

        for (auto iter = tsys->StageList.cbegin();
             iter != tsys->StageList.cend();
             ++iter)
        {
            retval = this->set_aperture_planes(*iter);
            if (!retval)
                break;
        }

        return retval;
    }

    bool NativeRunner::set_aperture_planes(tstage_ptr stage)
    {
        bool retval;

        for (auto eiter = stage->ElementList.begin();
             eiter != stage->ElementList.end();
             ++eiter)
        {
            retval = aperture_plane(*eiter);
            if (!retval)
                break;
        }

        return retval;
    }

    bool NativeRunner::aperture_plane(telement_ptr Element)
    {
        /*{Calculates the aperture plane of the element in element coord system.
        Applicable to rotationally symmetric apertures surfaces with small
        curvature: g, s, p, o, c, v, m, e, r, i.
          input - Element = Element record containing geometry of element
          output -
                 - Element.ZAperture  where ZAperture is the distance from
                   the origin to the plane.
        }*/

        Element->ZAperture =
            Element->icalc->compute_z_aperture(Element->aperture);

        return true;
    }

} // namespace SolTrace::NativeRunner

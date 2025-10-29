#include "simulation_result.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

#include "element.hpp"
#include "vector3d.hpp"

namespace SolTrace::Result
{

    using element_id = SolTrace::Data::element_id;
    using Vector3d = SolTrace::Data::Vector3d;

    const std::string &ray_event_string(const RayEvent rev)
    {
        auto item = REV_TO_STR.find(rev);
        if (item != REV_TO_STR.cend())
        {
            return item->second;
        }
        else
        {
            return REV_TO_STR.find(RayEvent::UNKNOWN)->second;
        }
    }

    InteractionRecord::InteractionRecord(
        element_id el,
        RayEvent rev,
        const Vector3d &location,
        const Vector3d &direction)
        // : index(-1),
        : element(el),
          event(rev),
          location(location),
          direction(direction)
    {
    }

    InteractionRecord::InteractionRecord(
        element_id el,
        RayEvent rev,
        double px,
        double py,
        double pz,
        double dx,
        double dy,
        double dz)
        : element(el),
          event(rev),
          location(px, py, pz),
          direction(dx, dy, dz)
    {
    }

    InteractionRecord::~InteractionRecord()
    {
    }

    std::ostream &operator<<(std::ostream &os, const InteractionRecord &rec)
    {
        os << "Element: " << rec.element
           << " Location: " << rec.location
           << " Direction: " << rec.direction
           << " Event: " << ray_event_string(rec.event);
        return os;
    }

    RayRecord::RayRecord(ray_id id) : id(id)
    {
    }

    RayRecord::~RayRecord()
    {
        this->interactions.clear();
        return;
    }

    void RayRecord::add_interaction_record(interaction_ptr ip)
    {
        // ip->index = this->interactions.size();
        this->interactions.push_back(ip);
        return;
    }

    void RayRecord::get_position(const interaction_ptr ip, Vector3d &pos)
    {
        pos = ip->location;
        return;
    }

    void RayRecord::get_position(int_fast64_t idx, Vector3d &pos)
    {
        const interaction_ptr ip0 = this->interactions[idx];
        return this->get_position(ip0, pos);
    }

    void RayRecord::get_direction(const interaction_ptr ip,
                                  Vector3d &cos)
    {
        cos = ip->direction;
        make_unit_vector(cos);
        return;
    }

    void RayRecord::get_direction(int_fast64_t idx,
                                  Vector3d &cos)
    {
        const interaction_ptr ip0 = this->interactions[idx];
        this->get_direction(ip0, cos);
        return;
    }

    const interaction_ptr &RayRecord::operator[](int_fast64_t idx) const
    {
        if (idx < 0 || idx >= this->interactions.size())
        {
            std::stringstream ss;
            ss << "RayRecord: Index " << idx
               << " is out of bounds [0, " << this->interactions.size()
               << "].";
            throw std::invalid_argument(ss.str());
        }
        return this->interactions[idx];
    }

    std::ostream &operator<<(std::ostream &os, const RayRecord &rec)
    {
        os << "Ray: " << rec.id
           << "\nInteractions:\n";
        for (uint_fast64_t k = 0; k < rec.interactions.size(); ++k)
        {
            os << k << " " << *rec.interactions[k] << "\n";
        }
        // os << "Last Cosines: " << rec.cos_last << "\n";
        return os;
    }

    SimulationResult::SimulationResult()
    {
        return;
    }

    SimulationResult::~SimulationResult()
    {
        this->ray_history.clear();
        return;
    }

    void SimulationResult::add_ray_record(ray_record_ptr rp)
    {
        this->ray_history.push_back(rp);
        return;
    }

    void SimulationResult::write_csv_file(std::string csv_name,
                                          int precision)
    {
        return this->write_csv_file(csv_name.c_str(), precision);
    }

    void SimulationResult::write_csv_file(const char *csv_name,
                                          int precision)
    {
        std::ofstream csv(csv_name);
        csv.precision(precision);
        csv << "Ray Number,Pos X,Pos Y,Pos Z,"
            << "Cos X,Cos Y,Cos Z,Element,Event\n";
        for (auto srit : this->ray_history)
        {
            for (auto cit : srit->interactions)
            {
                csv << srit->id << ","
                    << cit->location[0] << ","
                    << cit->location[1] << ","
                    << cit->location[2] << ","
                    << cit->direction[0] << ","
                    << cit->direction[1] << ","
                    << cit->direction[2] << ","
                    << cit->element << ","
                    << ray_event_string(cit->event) << "\n";
            }
        }
        csv.close();
        return;
    }

    const ray_record_ptr &SimulationResult::operator[](int_fast64_t idx) const
    {
        if (idx < 0 || idx >= this->ray_history.size())
        {
            std::stringstream ss;
            ss << "SimulationResult: Index " << idx
               << " is out of bounds [0, " << this->ray_history.size() - 1
               << "].";
            throw std::invalid_argument(ss.str());
        }
        return this->ray_history[idx];
    }

    // ray_record_ptr &SimulationResult::operator[](int_fast64_t idx)
    // {
    //     // TODO: Probably should do bounds checking...
    //     return this->ray_history[idx];
    // }

    std::ostream &operator<<(std::ostream &os, const SimulationResult &simres)
    {
        os << "Simulation Results -- " << simres.ray_history.size() << " Rays\n";
        for (uint_fast64_t k = 0; k < simres.ray_history.size(); ++k)
        {
            // os << "Ray: " << k << "\n"
            //    << *simres.ray_history[k];
            os << *simres.ray_history[k];
        }
        return os;
    }

} // namespace SolTrace::Result

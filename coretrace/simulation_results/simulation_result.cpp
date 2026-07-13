#include "simulation_result.hpp"

#define GLM_ENABLE_EXPERIMENTAL 1
#include <glm/gtx/io.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>

#include "element.hpp"
#include "records.hpp"

namespace SolTrace::Result
{
    using element_id = SolTrace::Data::element_id;

    RayHistoryResult::RayHistoryResult()
    {
        return;
    }

    RayHistoryResult::~RayHistoryResult()
    {
        this->ray_history.clear();
        this->element_view.clear();
        return;
    }

    void RayHistoryResult::add_ray_record(ray_record_ptr rp)
    {
        this->ray_history.push_back(rp);
        this->add_element_view(rp);
        return;
    }

    const element_record_ptr RayHistoryResult::get_element_record(element_id elid) const
    {
        auto iter = this->element_view.find(elid);
        element_record_ptr retval = nullptr;
        if (iter == this->element_view.end())
        {
            // Nothing in the record hit this element. Create an empty
            // record to return.
            retval = make_element_record(elid);
        }
        else
        {
            retval = iter->second;
        }
        return retval;
    }

    void RayHistoryResult::write_csv_file(std::string csv_name,
                                          int precision) const
    {
        return this->write_csv_file(csv_name.c_str(), precision);
    }

    void RayHistoryResult::write_csv_file(const char *csv_name,
                                          int precision) const
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

    const ray_record_ptr &RayHistoryResult::operator[](int_fast64_t idx) const
    {
        if (idx < 0 || idx >= this->ray_history.size())
        {
            std::stringstream ss;
            ss << "RayHistoryResult: Index " << idx
               << " is out of bounds [0, " << this->ray_history.size() - 1
               << "].";
            throw std::invalid_argument(ss.str());
        }
        return this->ray_history[idx];
    }

    std::ostream &operator<<(std::ostream &os, const RayHistoryResult &simres)
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

    void RayHistoryResult::add_element_view(const ray_record_ptr rp)
    {
        element_record_ptr erec;
        interaction_ptr ip;

        for (auto iter = rp->interactions.cbegin();
             iter != rp->interactions.cend();
             ++iter)
        {
            ip = *iter;
            auto eiter = this->element_view.find(ip->element);

            if (eiter == this->element_view.cend())
            {
                erec = make_element_record(ip->element);
                this->element_view.insert(std::make_pair(ip->element, erec));
            }
            else
            {
                erec = eiter->second;
            }

            erec->interactions.push_back(ip);
        }
        return;
    }

} // namespace SolTrace::Result

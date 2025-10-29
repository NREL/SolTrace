#ifndef SOLTRACE_SIMULATION_RESULT_H
#define SOLTRACE_SIMULATION_RESULT_H

#include <cstdint>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "element.hpp"
#include "vector3d.hpp"

namespace SolTrace::Result
{
    using ray_id = int_fast64_t;

    enum class RayEvent
    {
        CREATE = 1,
        ABSORB = 2,
        REFLECT = 3,
        TRANSMIT = 4,
        VIRTUAL = 5,
        EXIT = 6,
        UNKNOWN = 1000
    };

    const std::map<RayEvent, std::string> REV_TO_STR{
        {RayEvent::CREATE, "CREATE"},
        {RayEvent::ABSORB, "ABSORB"},
        {RayEvent::REFLECT, "REFLECT"},
        {RayEvent::TRANSMIT, "TRANSMIT"},
        {RayEvent::EXIT, "EXIT"},
        {RayEvent::UNKNOWN, "UNKNONWN"}};

    const std::string &ray_event_string(RayEvent rev);

    struct InteractionRecord
    {
        // int_fast64_t index;
        SolTrace::Data::element_id element;
        RayEvent event;
        SolTrace::Data::Vector3d location;
        SolTrace::Data::Vector3d direction;

        InteractionRecord(SolTrace::Data::element_id el, RayEvent rev,
                          const SolTrace::Data::Vector3d &location,
                          const SolTrace::Data::Vector3d &direction);
        InteractionRecord(SolTrace::Data::element_id el, RayEvent rev,
                          double px, double py, double pz,
                          double dx, double dy, double dz);
        ~InteractionRecord();

        friend std::ostream &operator<<(std::ostream &os,
                                        const InteractionRecord &rec);
    };

    using interaction_ptr = std::shared_ptr<InteractionRecord>;
    template <typename... Args>
    inline auto make_interaction_record(Args &&...args)
    {
        return std::make_shared<InteractionRecord>(std::forward<Args>(args)...);
    }

    struct RayRecord
    {
        ray_id id;
        std::vector<interaction_ptr> interactions;

        RayRecord(ray_id id);
        ~RayRecord();

        // Assumes that interactions are added in order
        void add_interaction_record(interaction_ptr ip);

        SolTrace::Data::element_id get_element(const interaction_ptr ip)
        {
            return ip->element;
        }
        SolTrace::Data::element_id get_element(int_fast64_t idx)
        {
            return get_element(this->interactions[idx]);
        }
        RayEvent get_event(const interaction_ptr ip)
        {
            return ip->event;
        }
        RayEvent get_event(int_fast64_t idx)
        {
            return get_event(this->interactions[idx]);
        }

        void get_position(const interaction_ptr ip, SolTrace::Data::Vector3d &pos);
        void get_position(int_fast64_t idx, SolTrace::Data::Vector3d &pos);
        void get_direction(const interaction_ptr ip, SolTrace::Data::Vector3d &cos);
        void get_direction(int_fast64_t idx, SolTrace::Data::Vector3d &cos);

        uint_fast64_t get_number_of_interactions()
        {
            return this->interactions.size();
        }

        const interaction_ptr &operator[](int_fast64_t idx) const;
        friend std::ostream &operator<<(std::ostream &os,
                                        const RayRecord &rec);
    };

    using ray_record_ptr = std::shared_ptr<RayRecord>;
    template <typename... Args>
    inline auto make_ray_record(Args &&...args)
    {
        return std::make_shared<RayRecord>(std::forward<Args>(args)...);
    }

    using RayRecordContainer = typename std::vector<ray_record_ptr>;

    class SimulationResult
    {
    public:
        SimulationResult();
        virtual ~SimulationResult();

        void add_ray_record(ray_record_ptr);
        uint_fast64_t get_number_of_records() const
        {
            return this->ray_history.size();
        }
        RayRecordContainer::const_iterator get_ray_record_iteratior()
        {
            return ray_history.cbegin();
        }
        bool is_at_end(RayRecordContainer::const_iterator citer)
        {
            return citer == this->ray_history.cend();
        }

        void write_csv_file(std::string csv_name, int precision=12);
        void write_csv_file(const char *csv_name, int precision=12);

        // Legacy stuff -- TODO:
        // void results_to_legacy_csv(std::string csv_name,
        //                            SimulationData *data);

        const ray_record_ptr &operator[](int_fast64_t idx) const;
        // ray_record_ptr &operator[](int_fast64_t idx);
        friend std::ostream &operator<<(std::ostream &os,
                                        const SimulationResult &simres);

    private:
        RayRecordContainer ray_history;
    };

} // namespace SolTrace::Result

#endif

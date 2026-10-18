#ifndef SOLTRACE_RECORDS_H
#define SOLTRACE_RECORDS_H

#include <cstdint>
#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include <element.hpp>
#include <glm/vec3.hpp>

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
        {RayEvent::VIRTUAL, "VIRTUAL"},
        {RayEvent::EXIT, "EXIT"},
        {RayEvent::UNKNOWN, "UNKNOWN"}};

    const std::string &ray_event_string(RayEvent rev);

    struct InteractionRecord
    {
        // int_fast64_t index;
        SolTrace::Data::element_id element;
        RayEvent event;
        glm::dvec3 location;
        glm::dvec3 direction;

        InteractionRecord(SolTrace::Data::element_id el, RayEvent rev,
                          const glm::dvec3 &location,
                          const glm::dvec3 &direction);
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
    using InteractionRecordContainer = typename std::vector<interaction_ptr>;

    struct RayRecord
    {
        ray_id id;
        InteractionRecordContainer interactions;

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

        void get_position(const interaction_ptr ip, glm::dvec3 &pos);
        void get_position(int_fast64_t idx, glm::dvec3 &pos);
        void get_direction(const interaction_ptr ip, glm::dvec3 &cos);
        void get_direction(int_fast64_t idx, glm::dvec3 &cos);

        uint_fast64_t get_number_of_interactions()
        {
            return this->interactions.size();
        }
        InteractionRecordContainer::const_iterator get_interaction_record_iterator() const
        { return interactions.cbegin(); }
        bool is_at_end(InteractionRecordContainer::const_iterator citer) const
        { return citer == this->interactions.cend(); }

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

    struct ElementRecord
    {
        SolTrace::Data::element_id elid;
        InteractionRecordContainer interactions;

        ElementRecord(SolTrace::Data::element_id eid);
        ~ElementRecord();

        void add_interaction_record(interaction_ptr ip);

        uint_fast64_t get_number_of_interactions()
        {
            return this->interactions.size();
        }

        const interaction_ptr &operator[](int_fast64_t idx) const;
        friend std::ostream &operator<<(std::ostream &os,
                                        const ElementRecord &rec);
    };

    using element_record_ptr = std::shared_ptr<ElementRecord>;
    template <typename... Args>
    inline auto make_element_record(Args &&...args)
    {
        return std::make_shared<ElementRecord>(std::forward<Args>(args)...);
    }

} // namespace SolTrace::Result

#endif

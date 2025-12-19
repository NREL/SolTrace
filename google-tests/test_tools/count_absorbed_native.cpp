#include "count_absorbed_native.h"

#include <map>

#include <native_runner_types.hpp>

uint_fast64_t count_absorbed_native(const SolTrace::NativeRunner::TRayData *ray_data)
{
    return count_event_native(ray_data,
                              SolTrace::Result::RayEvent::ABSORB);
}

uint_fast64_t count_event_native(const SolTrace::NativeRunner::TRayData *ray_data,
                                 SolTrace::Result::RayEvent event)
{
    uint_fast64_t count = 0;
    size_t n = ray_data->Count();
    for (size_t i = 0; i < n; i++)
    {
        double pos[3], cos[3];
        int elm, stage;
        uint_fast64_t ray;
        SolTrace::Result::RayEvent rev;
        if (ray_data->Query(i, pos, cos, &elm, &stage, &ray, &rev))
        {
            if (rev == event)
                ++count;
        }
    }
    return count;
}

void scan_events_native(const SolTrace::NativeRunner::TRayData *ray_data)
{
    size_t n = ray_data->Count();
    std::map<unsigned, SolTrace::Result::RayEvent> ray_end;
    for (size_t i = 0; i < n; i++)
    {
        double pos[3], cos[3];
        int elm, stage;
        uint_fast64_t ray;
        SolTrace::Result::RayEvent rev;
        if (ray_data->Query(i, pos, cos, &elm, &stage, &ray, &rev))
        {
            if (rev == SolTrace::Result::RayEvent::ABSORB ||
                rev == SolTrace::Result::RayEvent::EXIT)
            {
                if (ray_end.find(ray) == ray_end.cend())
                {
                    ray_end[ray] = rev;
                }
                else
                {
                    std::cout << "Ray number " << ray
                              << " terminated with status "
                              << SolTrace::Result::ray_event_string(ray_end[ray])
                              << " and "
                              << SolTrace::Result::ray_event_string(rev)
                              << std::endl;
                }
            }
        }
    }
    std::cout << "Scanned " << ray_end.size() << " rays." << std::endl;
    return;
}

#ifndef SOLTRACE_GROUP_RESULT_H
#define SOLTRACE_GROUP_RESULT_H

#include <cstdint>
#include <map>
#include <string>
#include <vector>
#include <nlohmann/json.hpp>

#include "records.hpp"

namespace SolTrace::Result
{
    struct GroupResult
    {
        int32_t group_id;
        // counts for each ray event type
        uint_fast64_t absorb_count, reflect_count, transmit_count, virtual_count;
        // counts for ray events where the previous element was the sun
        uint_fast64_t absorb_sun_previous, reflect_sun_previous, transmit_sun_previous, virtual_sun_previous;
        // counts for ray events based on the previous element's group
        // vector index is group number of previous element
        std::vector<uint_fast64_t> absorb_previous_group, reflect_previous_group, transmit_previous_group, virtual_previous_group;
        
        GroupResult() : 
            group_id(-1),
            absorb_count(0), reflect_count(0), transmit_count(0), virtual_count(0),
            absorb_sun_previous(0), reflect_sun_previous(0), transmit_sun_previous(0), virtual_sun_previous(0)
        {};
        GroupResult(int32_t group_id, size_t n_groups);

        void write_json(nlohmann::ordered_json& jnode) const;

        void increment(RayEvent rev, int32_t prev_group);
        
        // TODO: 
        // [] make function to calculate ungrouped counts by subtracting grouped and sun counted from totals.
        // [] probably add flux maps to groups
            /*
            HPM2D fluxGrid;
            std::vector<double> xValues, yValues;
            double binszx, binszy;
            double PeakFlux, PeakFluxUncertainty;
            double AveFlux, AveFluxUncertainty;
            double MinFlux, SigmaFlux, Uniformity;
            glm::dvec3 Centroid;
            double zScale;
            size_t NumberOfRays;
            int NotBinned;
            double max_neg_x_flux_err = 0;
            double max_pos_x_flux_err = 0;
            */
        // [] make fucntion for comparisions between GroupResult structs, get efficiency and stuff like that.
    };
}// namespace SolTrace::Result

#endif // SOLTRACE_GROUP_RESULT_H
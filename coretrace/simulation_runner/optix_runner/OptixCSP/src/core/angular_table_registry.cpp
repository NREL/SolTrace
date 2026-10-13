#include "angular_table_registry.h"

using namespace OptixCSP;
using SolTrace::Data::OpticalPropertySet;
using SolTrace::Data::OpticalSide;
using SolTrace::Data::optics_id;

namespace
{
void append_as_float(std::vector<float>& dest, const std::vector<double>& src)
{
    dest.reserve(dest.size() + src.size());
    for (double value : src)
        dest.push_back(static_cast<float>(value));
}
} // namespace

AngularTableRange AngularTableRegistry::intern(optics_id                 id,
                                               const OpticalPropertySet& optics,
                                               OpticalSide               side,
                                               bool use_refraction)
{
    const auto key = std::make_pair(id, side);
    const auto it  = m_cache.find(key);
    if (it != m_cache.end()) return it->second;

    const bool has_table = use_refraction
                               ? optics.uses_transmissivity_table(side)
                               : optics.uses_reflectivity_table(side);

    AngularTableRange range;
    if (has_table)
    {
        const std::vector<double>& cosines =
            use_refraction ? optics.get_transmissivity_cosine_cache(side)
                           : optics.get_reflectivity_cosine_cache(side);
        const std::vector<double> values =
            use_refraction ? optics.get_transmissivity_values(side)
                           : optics.get_reflectivity_values(side);

        range.offset = static_cast<uint32_t>(m_cos_pool.size());
        range.count  = static_cast<uint32_t>(cosines.size());

        append_as_float(m_cos_pool, cosines);
        append_as_float(m_value_pool, values);
    }

    m_cache.emplace(key, range);

    return range;
}

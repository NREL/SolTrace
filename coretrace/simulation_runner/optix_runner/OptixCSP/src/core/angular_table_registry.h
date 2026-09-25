#pragma once

#include <cstdint>
#include <map>
#include <utility>
#include <vector>

#include "optical_properties.hpp"

namespace OptixCSP
{
// Offset/count of a table within AngularTableRegistry's shared device pool.
struct AngularTableRange
{
    uint32_t offset = 0;
    uint32_t count  = 0;
};

// Collects angle-dependent reflectivity/transmissivity tables into a
// single device pool, keyed by optics_id so that every element sharing
// the same OpticalPropertySet reuses one entry. The table itself is only
// read out of `optics` on a cache miss, so sharing optics across many
// elements never costs more than one extraction.
class AngularTableRegistry
{
public:
    // Returns the pooled range for `id`'s `side` table, appending to the
    // shared pool only the first time that (id, side) pair is seen.
    // `use_refraction` selects whether the transmissivity or reflectivity
    // table is the one actually consulted at trace time. Returns an
    // empty range (offset=0, count=0) if no table is enabled.
    AngularTableRange intern(SolTrace::Data::optics_id                 id,
                             const SolTrace::Data::OpticalPropertySet& optics,
                             SolTrace::Data::OpticalSide               side,
                             bool use_refraction);

    const std::vector<float>& cos_pool() const { return m_cos_pool; }
    const std::vector<float>& value_pool() const { return m_value_pool; }

private:
    std::map<std::pair<SolTrace::Data::optics_id, SolTrace::Data::OpticalSide>,
             AngularTableRange>
        m_cache;

    std::vector<float> m_cos_pool;
    std::vector<float> m_value_pool;
};
} // namespace OptixCSP

#include "determine_interaction_type.hpp"

#include <cmath>
#include <sstream>

#include <glm/geometric.hpp>

#include <simulation_data_export.hpp>

#include "mtrand.hpp"
#include "native_runner_types.hpp"
#include "trace_logger.hpp"

namespace SolTrace::NativeRunner
{
    using SolTrace::Result::RayEvent;

    bool determine_interaction_type(trace_logger_ptr logger,
                                    int_fast64_t stage,
                                    unsigned thread_id,
                                    MTRand &myrng,
                                    const OpticalPropertySet *optics,
                                    const glm::dvec3 &LastDFXYZ,
                                    const glm::dvec3 &LastCosRaySurfElement,
                                    bool LastHitBackSide,
                                    RayEvent &rev)
    {
        bool good = true;
        rev = RayEvent::VIRTUAL;

        double TestValue;

        const OpticalSide side = LastHitBackSide == false ? OpticalSide::Front : OpticalSide::Back;

        // Only pay for the incident angle when a table lookup actually needs it.
        double IncidentCos = 1.0; // [mrad]
        if (optics->uses_angular_table(side))
        {
            const glm::dvec3 UnitLastDFXYZ = -glm::normalize(LastDFXYZ);
            IncidentCos = glm::dot(LastCosRaySurfElement, UnitLastDFXYZ);
        }

        switch (optics->get_interaction_type())
        {
        case InteractionType::REFRACTION:
            TestValue = optics->get_transmissivity(side, IncidentCos);
            rev = RayEvent::TRANSMIT;
            break;
        case InteractionType::REFLECTION:
            TestValue = optics->get_reflectivity(side, IncidentCos);
            rev = RayEvent::REFLECT;
            break;
        default:
            good = false;
            std::stringstream ss;
            ss << "Bad optical interaction."
               << " Type: " << static_cast<int>(optics->get_interaction_type())
               << " Stage: " << stage
               << " Thread: " << thread_id
               << "\n";
            logger->error_log(ss.str());
            return false;
        }

        // Apply MonteCarlo probability of absorption. Limited
        // for now, but can make more complex later on if desired
        if (TestValue <= myrng())
        {
            // ray was fully absorbed
            rev = RayEvent::ABSORB;
        }

        return good;
    }

} // namespace SolTrace::NativeRunner

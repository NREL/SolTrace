
#include "calculator_factory.hpp"

#include "cylinder_calculator.hpp"
#include "flat_calculator.hpp"
#include "native_runner_types.hpp"
#include "newton_calculator.hpp"
#include "parabola_calculator.hpp"
#include "sphere_calculator.hpp"
#include "simulation_data_export.hpp"
#include "surface.hpp"

#include <stdexcept>
#include <sstream>

namespace SolTrace::NativeRunner {

CalculatorFactory *CalculatorFactory::instance = nullptr;

calculator_ptr CalculatorFactory::make_calculator(
    aperture_ptr ap,
    surface_ptr surf,
    const ElementParameters &eparams)
{
    // Input validation
    if (surf == nullptr)
    {
        throw std::invalid_argument(
            "CalculatorFactory::make_calculator: Surface pointer cannot be null");
    }

    if (ap == nullptr)
    {
        throw std::invalid_argument(
            "CalculatorFactory::make_calculator: Aperture pointer cannot be null");
    }

    SurfaceType st = surf->get_type();
    calculator_ptr calc = nullptr;

    if (st == SurfaceType::PARABOLA)
    {
        calc = std::make_shared<ParabolaCalculator>(surf, ap);
    }
    else if (st == SurfaceType::FLAT)
    {
        calc = std::make_shared<FlatCalculator>(surf, ap);
    }
    else if (st == SurfaceType::CYLINDER)
    {
        calc = std::make_shared<CylinderCalculator>(surf, ap);
    }
    else if (st == SurfaceType::SPHERE)
    {
        calc = std::make_shared<SphereCalculator>(surf, ap);
    }
    else
    {
        std::stringstream ss;
        ss << "CalculatorFactory::make_calculator: Unsupported surface type: "
           << static_cast<int>(st);
        throw std::invalid_argument(ss.str());
    }

    // This should never happen but just in case...
    if (calc == nullptr)
    {
        std::stringstream ss;
        ss << "CalculatorFactory::make_calculator: Failed to create calculator for surface type "
           << static_cast<int>(st);
        throw std::runtime_error(ss.str());
    }

    return calc;
}

} // namespace SolTrace::NativeRunner

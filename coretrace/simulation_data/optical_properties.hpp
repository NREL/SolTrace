/**
 * @file optical_properties.hpp
 * @brief Optical properties definitions for materials
 *
 * Defines optical properties (reflectivity, transmissivity, refractive index)
 * and interaction types for optical surfaces and materials. Includes error
 * distribution parameters for modeling surface imperfections and optical
 * errors.
 */

#ifndef SOLTRACE_OPTICAL_PROPERTIES_H
#define SOLTRACE_OPTICAL_PROPERTIES_H

#include <container.hpp>
#include <cstdint>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <nlohmann/json.hpp>

#include "error_distributions.hpp"

namespace SolTrace::Data
{
using optics_id = std::int_fast64_t;

enum OPTICS_ID_TYPES
{
    OPTICS_ID_VIRTUAL    = -3,
    OPTICS_ID_UNASSIGNED = -2
};

struct OpticalPropertySet;
struct OpticalPropertySetReference;
using OpticalPropertySetContainer = Container<optics_id, OpticalPropertySet>;
using OpticalPropertySetResolver =
    std::function<OpticalPropertySetReference(const optics_id)>;

enum class InteractionType
{
    REFLECTION,
    REFRACTION,
    UNKNOWN
};

inline const std::map<InteractionType, std::string> InteractionTypeMap = {
    { InteractionType::REFLECTION, "REFLECTION" },
    { InteractionType::REFRACTION, "REFRACTION" },
    { InteractionType::UNKNOWN, "UNKNOWN" }
};

enum class OpticalSide
{
    Front,
    Back,
    Both
};

// A single (angle, value) sample of an angle-dependent reflectivity or
// transmissivity table. Angle is the incidence angle in [mrad].
struct AngularTablePoint
{
    double angle;
    double value;

    bool operator==(const AngularTablePoint& other) const
    { return this->angle == other.angle && this->value == other.value; }

    bool operator!=(const AngularTablePoint& other) const
    { return !(*this == other); }
};

struct OpticalPropertySetReference
{
    optics_id                               id = OPTICS_ID_UNASSIGNED;
    std::weak_ptr<const OpticalPropertySet> optical_property_set;
};

class OpticalPropertySet
{

    class OpticalPropertiesFace
    {
    public:
        DistributionType error_distribution_type;
        double           transmissivity;
        double           reflectivity;
        double           slope_error;       // [mrad]
        double           specularity_error; // [mrad]

        // Angle-dependent tables, keyed on incidence angle [mrad] and
        // interpolated linearly in cos(angle). Used instead of the
        // scalar reflectivity/transmissivity above when enabled.
        bool                           use_reflectivity_table;
        bool                           use_transmissivity_table;
        std::vector<AngularTablePoint> reflectivity_table;
        std::vector<AngularTablePoint> transmissivity_table;
        std::vector<double>            reflectivity_cache;
        std::vector<double>            transmissivity_cache;

        OpticalPropertiesFace()
            : error_distribution_type(DistributionType::UNKNOWN),
              transmissivity(0.0),
              reflectivity(0.0),
              slope_error(0.0),
              specularity_error(0.0),
              use_reflectivity_table(false),
              use_transmissivity_table(false)
        {
        }

        OpticalPropertiesFace(DistributionType dtype,
                              double           trans,
                              double           refl,
                              double           slope_err,
                              double           spec_err)
            : error_distribution_type(dtype),
              transmissivity(trans),
              reflectivity(refl),
              slope_error(slope_err),
              specularity_error(spec_err),
              use_reflectivity_table(false),
              use_transmissivity_table(false)
        { validate(); }

        OpticalPropertiesFace(const nlohmann::ordered_json& jnode);

        void validate() const;

        // Returns the scalar reflectivity, or the table value linearly
        // interpolated on incident_cosine (cosine of the incidence angle)
        // if a table is in use.
        double get_reflectivity(double incident_cosine) const;

        // Returns the scalar transmissivity, or the table value linearly
        // interpolated on incident_cosine (cosine of the incidence angle)
        // if a table is in use.
        double get_transmissivity(double incident_cosine) const;

        void update_table_cache();
        void cache_cosines(std::vector<double>&                  cache,
                           const std::vector<AngularTablePoint>& table);

        // TODO: What should the error settings be with the below?
        void write_json(nlohmann::ordered_json& jnode) const;

        bool operator==(const OpticalPropertiesFace& other) const;
        bool operator!=(const OpticalPropertiesFace& other) const;

        // cosine_cache holds std::cos(table[i].angle * kMradToRad), precomputed
        // whenever the table is set, so the hot lookup never calls std::acos.
        static double
        interpolate_table(const std::vector<double>&            cosine_cache,
                          const std::vector<AngularTablePoint>& table,
                          double incident_cosine);
    };

    OpticalPropertiesFace front;
    OpticalPropertiesFace back;

    InteractionType my_type;

    double refraction_index_front;
    double refraction_index_back;

    std::string my_name;

    void set_ideal_material(OpticalSide side)
    {
        if (side == OpticalSide::Front || side == OpticalSide::Both)
        {
            this->front.error_distribution_type = DistributionType::NONE;
            this->front.specularity_error       = 0.0;
            this->front.slope_error             = 0.0;
        }

        if (side == OpticalSide::Back || side == OpticalSide::Both)
        {
            this->back.error_distribution_type = DistributionType::NONE;
            this->back.specularity_error       = 0.0;
            this->back.slope_error             = 0.0;
        }

        return;
    }

public:
    OpticalPropertySet(InteractionType interaction_type,
                       double          refrac_front,
                       double          refrac_back,
                       std::string     name = "")
        : front(),
          back(),
          my_type(interaction_type),
          refraction_index_front(refrac_front),
          refraction_index_back(refrac_back),
          my_name(name) { };

    OpticalPropertySet(InteractionType interaction_type, std::string name = "")
        : OpticalPropertySet(interaction_type, 0, 0, name) { };

    OpticalPropertySet()
        : my_type(InteractionType::UNKNOWN),
          refraction_index_front(0.0),
          refraction_index_back(0.0),
          my_name("") { };

    OpticalPropertySet(const nlohmann::ordered_json& jnode);

    void set_properties(const OpticalSide side,
                        DistributionType  dtype,
                        double            trans,
                        double            refl,
                        double            slope_err,
                        double            spec_err)
    {
        auto set_face_props = [=](OpticalPropertiesFace& face)
        {
            face.error_distribution_type = dtype;
            face.transmissivity          = trans;
            face.reflectivity            = refl;
            face.slope_error             = slope_err;
            face.specularity_error       = spec_err;
        };

        if (side == OpticalSide::Front || side == OpticalSide::Both)
        {
            set_face_props(this->front);
            this->front.validate();
        }
        if (side == OpticalSide::Back || side == OpticalSide::Both)
        {
            set_face_props(this->back);
            this->back.validate();
        }

        return;
    }

    void set_interaction_type(const InteractionType type)
    { this->my_type = type; }

    void set_refraction_indices(double rfront, double rback)
    {
        this->refraction_index_front = rfront;
        this->refraction_index_back  = rback;
    }

    void set_reflectivity(const OpticalSide side, double refl)
    {
        if (side == OpticalSide::Front || side == OpticalSide::Both)
        {
            this->front.reflectivity = refl;
        }

        if (side == OpticalSide::Back || side == OpticalSide::Both)
        {
            this->back.reflectivity = refl;
        }
    }

    void set_transmissivity(const OpticalSide side, double trans)
    {
        if (side == OpticalSide::Front || side == OpticalSide::Both)
        {
            this->front.transmissivity = trans;
        }

        if (side == OpticalSide::Back || side == OpticalSide::Both)
        {
            this->back.transmissivity = trans;
        }
    }

    // Enables and sets an angle-dependent reflectivity table (angle in
    // [mrad], ascending) for the given side(s).
    void set_reflectivity_table(const OpticalSide                     side,
                                const std::vector<AngularTablePoint>& table)
    {
        if (side == OpticalSide::Front || side == OpticalSide::Both)
        {
            this->front.use_reflectivity_table = true;
            this->front.reflectivity_table     = table;
            this->front.validate();
            this->front.update_table_cache();
        }

        if (side == OpticalSide::Back || side == OpticalSide::Both)
        {
            this->back.use_reflectivity_table = true;
            this->back.reflectivity_table     = table;
            this->back.validate();
            this->back.update_table_cache();
        }
    }

    // Enables and sets an angle-dependent transmissivity table (angle in
    // [mrad], ascending) for the given side(s).
    void set_transmissivity_table(const OpticalSide                     side,
                                  const std::vector<AngularTablePoint>& table)
    {
        if (side == OpticalSide::Front || side == OpticalSide::Both)
        {
            this->front.use_transmissivity_table = true;
            this->front.transmissivity_table     = table;
            this->front.validate();
            this->front.update_table_cache();
        }

        if (side == OpticalSide::Back || side == OpticalSide::Both)
        {
            this->back.use_transmissivity_table = true;
            this->back.transmissivity_table     = table;
            this->back.validate();
            this->back.update_table_cache();
        }
    }

    // Toggles whether ray tracing uses the reflectivity table for the
    // given side(s), without discarding the table itself.
    void enable_reflectivity_table(const OpticalSide side)
    {
        if (side == OpticalSide::Front || side == OpticalSide::Both)
        {
            this->front.use_reflectivity_table = true;
            this->front.validate();
            this->front.update_table_cache();
        }

        if (side == OpticalSide::Back || side == OpticalSide::Both)
        {
            this->back.use_reflectivity_table = true;
            this->back.validate();
            this->back.update_table_cache();
        }
    }

    void disable_reflectivity_table(const OpticalSide side)
    {
        if (side == OpticalSide::Front || side == OpticalSide::Both)
        {
            this->front.use_reflectivity_table = false;
        }

        if (side == OpticalSide::Back || side == OpticalSide::Both)
        {
            this->back.use_reflectivity_table = false;
        }
    }

    // Toggles whether ray tracing uses the transmissivity table for the
    // given side(s), without discarding the table itself.
    void enable_transmissivity_table(const OpticalSide side)
    {
        if (side == OpticalSide::Front || side == OpticalSide::Both)
        {
            this->front.use_transmissivity_table = true;
            this->front.validate();
            this->front.update_table_cache();
        }

        if (side == OpticalSide::Back || side == OpticalSide::Both)
        {
            this->back.use_transmissivity_table = true;
            this->back.validate();
            this->back.update_table_cache();
        }
    }

    void disable_transmissivity_table(const OpticalSide side)
    {
        if (side == OpticalSide::Front || side == OpticalSide::Both)
        {
            this->front.use_transmissivity_table = false;
        }

        if (side == OpticalSide::Back || side == OpticalSide::Both)
        {
            this->back.use_transmissivity_table = false;
        }
    }

    void set_errors(const OpticalSide side,
                    DistributionType  dtype,
                    double            slope,
                    double            spec)
    {
        if (side == OpticalSide::Front || side == OpticalSide::Both)
        {
            this->front.error_distribution_type = dtype;
            this->front.slope_error             = slope;
            this->front.specularity_error       = spec;
            this->front.validate();
        }

        if (side == OpticalSide::Back || side == OpticalSide::Both)
        {
            this->back.error_distribution_type = dtype;
            this->back.slope_error             = slope;
            this->back.specularity_error       = spec;
            this->back.validate();
        }
    }

    void set_ideal_transmission()
    {
        this->my_type = InteractionType::REFRACTION;

        this->set_ideal_material(OpticalSide::Both);

        this->front.transmissivity = 1.0;
        this->front.reflectivity   = 0.0;

        this->back.transmissivity = 1.0;
        this->back.reflectivity   = 0.0;

        return;
    }

    void set_ideal_transmission(double refraction_index_front,
                                double refraction_index_back)
    {
        this->set_ideal_transmission();
        this->refraction_index_front = refraction_index_front;
        this->refraction_index_back  = refraction_index_back;
        return;
    }

    void
    set_ideal_one_sided_reflector(const OpticalSide side = OpticalSide::Front)
    {
        if (side == OpticalSide::Both)
        {
            throw std::invalid_argument("set_ideal_one_sided_reflector "
                                        "requires Front or Back, not Both.");
        }

        this->my_type = InteractionType::REFLECTION;

        this->set_ideal_material(OpticalSide::Both);

        this->front.reflectivity   = 0.0;
        this->front.transmissivity = 0.0;

        this->back.reflectivity   = 0.0;
        this->back.transmissivity = 0.0;

        if (side == OpticalSide::Front) { this->front.reflectivity = 1.0; }
        else
        {
            this->back.reflectivity = 1.0;
        }
    }

    void set_ideal_absorption(const OpticalSide side)
    {
        this->my_type = InteractionType::REFLECTION;
        this->set_ideal_material(side);

        if (side == OpticalSide::Front || side == OpticalSide::Both)
        {
            this->front.transmissivity = 0.0;
            this->front.reflectivity   = 0.0;
        }

        if (side == OpticalSide::Back || side == OpticalSide::Both)
        {
            this->back.transmissivity = 0.0;
            this->back.reflectivity   = 0.0;
        }

        return;
    }

    void set_ideal_reflection(const OpticalSide side)
    {
        this->my_type = InteractionType::REFLECTION;

        this->set_ideal_material(side);

        if (side == OpticalSide::Front || side == OpticalSide::Both)
        {
            this->front.transmissivity = 0.0;
            this->front.reflectivity   = 1.0;
        }

        if (side == OpticalSide::Back || side == OpticalSide::Both)
        {
            this->back.transmissivity = 0.0;
            this->back.reflectivity   = 1.0;
        }
        return;
    }

    const InteractionType get_interaction_type() const { return this->my_type; }

    const std::string get_name() const { return this->my_name; }

    void get_refraction_indices(double& rfront, double& rback) const
    {
        rfront = this->refraction_index_front;
        rback  = this->refraction_index_back;
    }

    double get_reflectivity(const OpticalSide side) const
    {
        switch (side)
        {
        case (OpticalSide::Front): return this->front.reflectivity;
        case (OpticalSide::Back): return this->back.reflectivity;
        default: return std::numeric_limits<double>::quiet_NaN();
        }
    }

    // Returns the reflectivity for the given side, using the
    // angle-dependent table (interpolated on incident_cosine, the cosine
    // of the incidence angle) if one is enabled for that side, otherwise
    // the scalar value.
    double get_reflectivity(const OpticalSide side,
                            double            incident_cosine) const
    {
        switch (side)
        {
        case (OpticalSide::Front):
            return this->front.get_reflectivity(incident_cosine);
        case (OpticalSide::Back):
            return this->back.get_reflectivity(incident_cosine);
        default: return std::numeric_limits<double>::quiet_NaN();
        }
    }

    double get_transmissivity(const OpticalSide side) const
    {
        switch (side)
        {
        case (OpticalSide::Front): return this->front.transmissivity;
        case (OpticalSide::Back): return this->back.transmissivity;
        default: return std::numeric_limits<double>::quiet_NaN();
        }
    }

    // Returns the transmissivity for the given side, using the
    // angle-dependent table (interpolated on incident_cosine, the cosine
    // of the incidence angle) if one is enabled for that side, otherwise
    // the scalar value.
    double get_transmissivity(const OpticalSide side,
                              double            incident_cosine) const
    {
        switch (side)
        {
        case (OpticalSide::Front):
            return this->front.get_transmissivity(incident_cosine);
        case (OpticalSide::Back):
            return this->back.get_transmissivity(incident_cosine);
        default: return std::numeric_limits<double>::quiet_NaN();
        }
    }

    // True if the given side has any angle-dependent reflectivity or
    // transmissivity table enabled.
    bool uses_angular_table(const OpticalSide side) const
    {
        switch (side)
        {
        case (OpticalSide::Front):
            return this->front.use_reflectivity_table ||
                   this->front.use_transmissivity_table;
        case (OpticalSide::Back):
            return this->back.use_reflectivity_table ||
                   this->back.use_transmissivity_table;
        default: return false;
        }
    }

    DistributionType get_error_distribution(const OpticalSide side) const
    {
        switch (side)
        {
        case (OpticalSide::Front): return this->front.error_distribution_type;
        case (OpticalSide::Back): return this->back.error_distribution_type;
        default: return DistributionType::UNKNOWN;
        }
    }

    double get_slope_error(const OpticalSide side) const
    {
        switch (side)
        {
        case (OpticalSide::Front): return this->front.slope_error;
        case (OpticalSide::Back): return this->back.slope_error;
        default: return std::numeric_limits<double>::quiet_NaN();
        }
    }

    double get_specularity_error(const OpticalSide side) const
    {
        switch (side)
        {
        case (OpticalSide::Front): return this->front.specularity_error;
        case (OpticalSide::Back): return this->back.specularity_error;
        default: return std::numeric_limits<double>::quiet_NaN();
        }
    }

    void get_errors(const OpticalSide side,
                    DistributionType& dtype,
                    double&           slope,
                    double&           spec) const
    {
        if (side == OpticalSide::Both)
        {
            dtype = DistributionType::UNKNOWN;
            slope = std::numeric_limits<double>::quiet_NaN();
            spec  = std::numeric_limits<double>::quiet_NaN();
            return;
        }

        auto& face = side == OpticalSide::Front ? this->front : this->back;

        dtype = face.error_distribution_type;
        slope = face.slope_error;
        spec  = face.specularity_error;
        return;
    }

    void write_json(nlohmann::ordered_json& jnode) const;

    bool operator==(const OpticalPropertySet& other) const;
    bool operator!=(const OpticalPropertySet& other) const;

    friend std::ostream& operator<<(std::ostream&             os,
                                    const OpticalPropertySet& op);

    friend std::ostream& operator<<(std::ostream&                os,
                                    const OpticalPropertiesFace& op);
};

} // namespace SolTrace::Data

#endif

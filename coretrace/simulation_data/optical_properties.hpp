/**
 * @file optical_properties.hpp
 * @brief Optical properties definitions for materials
 *
 * Defines optical properties (reflectivity, transmissivity, refractive index)
 * and interaction types for optical surfaces and materials. Includes error
 * distribution parameters for modeling surface imperfections and optical errors.
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

#include <nlohmann/json.hpp>

#include "error_distributions.hpp"

namespace SolTrace::Data
{
    using optics_id = std::int_fast64_t;

    enum OPTICS_ID_TYPES
    {
        OPTICS_ID_VIRTUAL = -3,
        OPTICS_ID_UNASSIGNED = -2        
    };

    struct OpticalPropertySet;
    struct OpticalPropertySetReference;
    using OpticalPropertySetContainer = Container<optics_id, OpticalPropertySet>;
    using OpticalPropertySetResolver = std::function<OpticalPropertySetReference(const optics_id)>;

    enum class InteractionType
    {
        REFLECTION,
        REFRACTION,
        UNKNOWN
    };

    inline const std::map<InteractionType, std::string> InteractionTypeMap =
    {
        {InteractionType::REFLECTION, "REFLECTION"},
        {InteractionType::REFRACTION, "REFRACTION"},
        {InteractionType::UNKNOWN, "UNKNOWN"}
    };

    enum class OpticalSide
    {
        Front,
        Back,
        Both
    };

    struct OpticalPropertySetReference
    {
        optics_id id = OPTICS_ID_UNASSIGNED;
        std::weak_ptr<const OpticalPropertySet> optical_property_set;
    };

    class OpticalPropertySet
    {

        class OpticalPropertiesFace
        {
        public:
            DistributionType error_distribution_type;
            double transmissivity;
            double reflectivity;
            double slope_error;                 // [mrad]
            double specularity_error;           // [mrad]

            OpticalPropertiesFace() : error_distribution_type(DistributionType::UNKNOWN),
                transmissivity(0.0),
                reflectivity(0.0),
                slope_error(0.0),
                specularity_error(0.0)
            {
            }

            OpticalPropertiesFace(DistributionType dtype,
                double trans, double refl,
                double slope_err, double spec_err)
                : error_distribution_type(dtype),
                transmissivity(trans),
                reflectivity(refl),
                slope_error(slope_err),
                specularity_error(spec_err)
            {
            }

            OpticalPropertiesFace(const nlohmann::ordered_json& jnode);

            // TODO: What should the error settings be with the below?
            void write_json(nlohmann::ordered_json& jnode) const;

            bool operator==(const OpticalPropertiesFace& other) const;
            bool operator!=(const OpticalPropertiesFace& other) const;

            
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
                this->front.specularity_error = 0.0;
                this->front.slope_error = 0.0;
            }

            if (side == OpticalSide::Back || side == OpticalSide::Both)
            {
                this->back.error_distribution_type = DistributionType::NONE;
                this->back.specularity_error = 0.0;
                this->back.slope_error = 0.0;
            }

            return;
        }

    public:

        OpticalPropertySet(InteractionType interaction_type, 
            double refrac_front, double refrac_back, 
            std::string name = "")
            : front(), back(),
            my_type(interaction_type), refraction_index_front(refrac_front),
            refraction_index_back(refrac_back),
            my_name(name)
        {
        };

        OpticalPropertySet(InteractionType interaction_type,
            std::string name = "")
            : OpticalPropertySet(interaction_type, 0, 0, name) 
        {
        };

        OpticalPropertySet()
            : my_type(InteractionType::UNKNOWN),
            refraction_index_front(0.0),
            refraction_index_back(0.0),
            my_name("")
        {
        };

        OpticalPropertySet(const nlohmann::ordered_json& jnode);

        void set_properties(const OpticalSide side, 
            DistributionType dtype,
            double trans, double refl,
            double slope_err, double spec_err)
        {
            auto set_face_props = [=](OpticalPropertiesFace& face)
            {
                face.error_distribution_type = dtype;
                face.transmissivity = trans;
                face.reflectivity = refl;
                face.slope_error = slope_err;
                face.specularity_error = spec_err;
            };

            if (side == OpticalSide::Front || side == OpticalSide::Both)
                set_face_props(this->front);
            if (side == OpticalSide::Back || side == OpticalSide::Both)
                set_face_props(this->back);

            return;
        }

        void set_interaction_type(const InteractionType type)
        {
            this->my_type = type;
        }

        void set_refraction_indices(double rfront, double rback)
        {
            this->refraction_index_front = rfront;
            this->refraction_index_back = rback;
        }

        void set_reflectivity(const OpticalSide side,
            double refl)
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

        void set_transmissivity(const OpticalSide side,
            double trans)
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

        void set_errors(const OpticalSide side,
            DistributionType dtype, double slope,
            double spec)
        {
            if (side == OpticalSide::Front || side == OpticalSide::Both)
            {
                this->front.error_distribution_type = dtype;
                this->front.slope_error = slope;
                this->front.specularity_error = spec;
            }

            if (side == OpticalSide::Back || side == OpticalSide::Both)
            {
                this->back.error_distribution_type = dtype;
                this->back.slope_error = slope;
                this->back.specularity_error = spec;
            }
        }

        void set_ideal_transmission()
        {
            this->my_type = InteractionType::REFRACTION;

            this->set_ideal_material(OpticalSide::Both);

            this->front.transmissivity = 1.0;
            this->front.reflectivity = 0.0;

            this->back.transmissivity = 1.0;
            this->back.reflectivity = 0.0;

            return;
        }

        void set_ideal_transmission(double refraction_index_front,
            double refraction_index_back)
        {
            this->set_ideal_transmission();
            this->refraction_index_front = refraction_index_front;
            this->refraction_index_back = refraction_index_back;
            return;
        }

        void set_ideal_one_sided_reflector(const OpticalSide side = OpticalSide::Front)
        {
            if (side == OpticalSide::Both)
            {
                throw std::invalid_argument("set_ideal_one_sided_reflector requires Front or Back, not Both.");
            }

            this->my_type = InteractionType::REFLECTION;

            this->set_ideal_material(OpticalSide::Both);

            this->front.reflectivity = 0.0;
            this->front.transmissivity = 0.0;

            this->back.reflectivity = 0.0;
            this->back.transmissivity = 0.0;

            if (side == OpticalSide::Front)
            {
                this->front.reflectivity = 1.0;
            }
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
                this->front.reflectivity = 0.0;
            }

            if (side == OpticalSide::Back || side == OpticalSide::Both)
            {
                this->back.transmissivity = 0.0;
                this->back.reflectivity = 0.0;
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
                this->front.reflectivity = 1.0;
            }

            if (side == OpticalSide::Back || side == OpticalSide::Both)
            {
                this->back.transmissivity = 0.0;
                this->back.reflectivity = 1.0;
            }
            return;
        }

        const InteractionType get_interaction_type() const
        {
            return this->my_type;
        }

        const std::string get_name() const
        {
            return this->my_name;
        }

        void get_refraction_indices(double& rfront, double& rback) const
        {
            rfront = this->refraction_index_front;
            rback = this->refraction_index_back;
        }

        double get_reflectivity(const OpticalSide side) const
        {
            switch (side)
            {
                case(OpticalSide::Front):
                    return this->front.reflectivity;
                case(OpticalSide::Back):
                    return this->back.reflectivity;
                default:
                    return std::numeric_limits<double>::quiet_NaN();
            }
        }

        double get_transmissivity(const OpticalSide side) const
        {
            switch (side)
            {
                case(OpticalSide::Front):
                    return this->front.transmissivity;
                case(OpticalSide::Back):
                    return this->back.transmissivity;
                default:
                    return std::numeric_limits<double>::quiet_NaN();
            }
        }

        DistributionType get_error_distribution(const OpticalSide side) const
        {
            switch (side)
            {
                case(OpticalSide::Front):
                    return this->front.error_distribution_type;
                case(OpticalSide::Back):
                    return this->back.error_distribution_type;
                default:
                    return DistributionType::UNKNOWN;
            }
        }

        double get_slope_error(const OpticalSide side) const
        {
            switch (side)
            {
                case(OpticalSide::Front):
                    return this->front.slope_error;
                case(OpticalSide::Back):
                    return this->back.slope_error;
                default:
                    return std::numeric_limits<double>::quiet_NaN();
            }
        }

        double get_specularity_error(const OpticalSide side) const
        {
            switch (side)
            {
                case(OpticalSide::Front):
                    return this->front.specularity_error;
                case(OpticalSide::Back):
                    return this->back.specularity_error;
                default:
                    return std::numeric_limits<double>::quiet_NaN();
            }
        }

        void get_errors(const OpticalSide side,
            DistributionType& dtype, double& slope,
            double& spec) const
        {
            if (side == OpticalSide::Both)
            {
                dtype = DistributionType::UNKNOWN;
                slope = std::numeric_limits<double>::quiet_NaN();
                spec = std::numeric_limits<double>::quiet_NaN();
                return;
            }

            auto& face = side == OpticalSide::Front ? this->front : this->back;

            dtype = face.error_distribution_type;
            slope = face.slope_error;
            spec = face.specularity_error;
            return;
        }

        void write_json(nlohmann::ordered_json& jnode) const;

        bool operator==(const OpticalPropertySet& other) const;
        bool operator!=(const OpticalPropertySet& other) const;

        friend std::ostream& operator<<(std::ostream& os,
            const OpticalPropertySet& op);

        friend std::ostream& operator<<(std::ostream& os,
            const OpticalPropertiesFace& op);
    };

} // namespace SolTrace::Data

#endif

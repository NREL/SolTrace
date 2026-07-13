#ifndef SOLTRACE_RESULT_SPEC_H
#define SOLTRACE_RESULT_SPEC_H

#include <memory>
#include <optional>
#include <set>

#include "simulation_result.hpp"

namespace SolTrace::Result
{
    /**
     * Abstract configuration base passed to setup_simulation().
     *
     * Carries the parameters the runner needs to configure its internal data
     * structures (e.g. GPU buffers, per-element filters) before tracing begins.
     * Also acts as a factory: create_result() produces the matching empty
     * SimulationResult that report_simulation() will populate.
     *
     * Subclasses add the fields relevant to their result type and override the
     * virtual query methods so runners can read configuration without casting.
     */
    struct ResultSpec
    {
        virtual ~ResultSpec() = default;

        virtual ResultType get_type() const = 0;

        /// Creates an empty result object of the appropriate concrete type.
        /// Call this after run_simulation() to obtain the container passed to
        /// report_simulation().
        virtual std::unique_ptr<SimulationResult> create_result() const = 0;

        /// Returns the element IDs whose hits should be recorded.
        /// std::nullopt means record all elements (default).
        virtual std::optional<std::set<SolTrace::Data::element_id>>
        get_element_filter() const
        {
            return std::nullopt;
        }

        /// Returns the ray event types that should be recorded.
        /// std::nullopt means record all event types (default).
        virtual std::optional<std::set<RayEvent>> get_event_filter() const
        {
            return std::nullopt;
        }

    protected:
        ResultSpec() = default;
    };

    /// Specification for a full per-ray interaction history result.
    struct RayHistorySpec : ResultSpec
    {
        /// If set, only hits on elements in this set will be recorded.
        std::optional<std::set<SolTrace::Data::element_id>> element_filter;

        /// If set, only interactions of these event types will be recorded.
        std::optional<std::set<RayEvent>> event_filter;

        ResultType get_type() const override { return ResultType::RAY_HISTORY; }

        std::unique_ptr<SimulationResult> create_result() const override
        {
            return std::make_unique<RayHistoryResult>();
        }

        std::optional<std::set<SolTrace::Data::element_id>>
        get_element_filter() const override
        {
            return element_filter;
        }

        std::optional<std::set<RayEvent>> get_event_filter() const override
        {
            return event_filter;
        }
    };

} // namespace SolTrace::Result

#endif

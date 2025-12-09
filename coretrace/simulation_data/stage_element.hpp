/**
 * @file stage_element.hpp
 * @brief Stage-based element organization
 *
 * Defines stage-based organization of optical elements,
 * allowing for multi-stage optical system modeling.
 * Stages provide logical grouping and sequential processing
 * of optical elements in ray tracing simulations.
 */

#ifndef SOLTRACE_STAGE_ELEMENT_H
#define SOLTRACE_STAGE_ELEMENT_H

#include <memory>
#include "composite_element.hpp"
#include "element.hpp"

namespace SolTrace::Data {

class StageElement: public CompositeElement
{
public:
    StageElement(int_fast64_t stage);
    StageElement(const nlohmann::ordered_json& jnode);
    ~StageElement();
    virtual bool is_stage() const override { return true; }
    virtual element_id add_element(element_ptr el);
    virtual void write_json(nlohmann::ordered_json& jnode) const override;
private:
};

using stage_ptr = std::shared_ptr<StageElement>;
template <typename... Args>
inline auto make_stage(Args &&...args)
{
    return make_element<StageElement>(std::forward<Args>(args)...);
}

} // namespace SolTrace::Data

#endif

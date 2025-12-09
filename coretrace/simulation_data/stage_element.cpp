#include "stage_element.hpp"

#include <sstream>

#include "element.hpp"
#include "json_helpers.hpp"


namespace SolTrace::Data {

StageElement::StageElement(int_fast64_t stage) : CompositeElement()
{
    this->set_stage(stage);
    return;
}

StageElement::StageElement(const nlohmann::ordered_json& jnode) : CompositeElement(jnode)
{
    // Check that it is a stage
    if (jnode.contains("is_stage") == false || jnode.at("is_stage").get<bool>() == false)
    {
        std::stringstream ss;
        ss << "JSON node is not a valid stage: ";
        if (jnode.contains("is_stage") == false) {
            ss << "missing 'is_stage' field.";
        } else {
            ss << "'is_stage' field is present but has value: " << jnode.at("is_stage") << ".";
        }
        throw std::invalid_argument(ss.str());
    }
    // Get and set stage number
    int stage = jnode.at("stage");
    this->set_stage(stage);
    
    // CompositeElement is responsible for populating elements
    // Its constructor has already been called
}

StageElement::~StageElement()
{
    return;
}


element_id StageElement::add_element(element_ptr el)
{
    element_id id = this->CompositeElement::add_element(el);
    if (Element::is_success(id))
    {
        el->set_stage(this->stage);
    }
    return id;
}

void StageElement::write_json(nlohmann::ordered_json& jnode) const
{
    // Write stage-specific json keys
    jnode["is_stage"] = true;
    jnode["is_single"] = false;

    // Call composite write_json
    this->CompositeElement::write_json(jnode);
}


} // namespace SolTrace::Data

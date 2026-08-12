#include "materialeditor.h"

#include "database.h"

namespace db {


MaterialEditor::MaterialEditor(QObject* parent)
    : QObject { parent },
      m_back_editor(new OpticalPropertiesObject(true, this)),
      m_front_editor(new OpticalPropertiesObject(false, this)),
      m_interaction_type_model(new QStringListModel(this)),
      m_distribution_type_model(new QStringListModel(this)) {

    build_options<SD::InteractionType>(*m_interaction_type_model);
    build_options<SD::DistributionType>(*m_distribution_type_model);
}

MaterialEditor::~MaterialEditor() = default;

void MaterialEditor::parameters_changed(entt::entity e) {
    if (this->m_current_group != e) return;
}

void MaterialEditor::set_new_database_connections(Database* ptr) {
    add_connection(connect(ptr->material_parameters.self(),
                           &ComponentAPIBase::changed,
                           this,
                           &MaterialEditor::parameters_changed));
}

void MaterialEditor::set(Database* database, entt::entity group) {
    observe(database);
    m_current_group = group;
    m_front_editor->set(database, group);
    m_back_editor->set(database, group);
    emit updated();
}

} // namespace db

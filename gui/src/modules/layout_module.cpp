#include "layout_module.h"


namespace SolTrace::GUI::App {

void LayoutModule::viewed_entity_changed() {
    child_model()->set_node(m_viewed_element);
    breadcrumb_model()->set_node(m_viewed_element);
}

void LayoutModule::edited_entity_changed() {
    instance_edit()->set(m_edited_element);

    if (!m_current_database) {
        set_edited_element_name(QString());
        return;
    }

    set_edited_element_name(m_current_database->name_of(m_edited_element));
}

void LayoutModule::reset(db::Database* database) {
    if (m_observed_database) {
        disconnect(
            m_observed_database->identity.self(), nullptr, this, nullptr);
    }

    m_observed_database = database;
    set_viewed_element({});
    set_edited_element({});

    if (!database) {
        set_edited_element_name(QString());
        return;
    }

    connect(database->identity.self(),
            &db::ComponentAPIBase::changed,
            this,
            &LayoutModule::identity_changed);

    set_edited_element_name(database->name_of(m_edited_element));
}

void LayoutModule::identity_changed(entt::entity entity) {
    if (!m_current_database) return;

    if (db::Entity(entity) != m_edited_element) return;

    set_edited_element_name(m_current_database->name_of(m_edited_element));
}

void LayoutModule::delete_edited_element() {
    if (!m_current_database) return;
    if (!m_edited_element.is_valid()) return;

    auto to_delete = m_edited_element;

    m_current_database->delete_element(to_delete);
    set_edited_element({});

    if (m_viewed_element == to_delete) { set_viewed_element({}); }
}

LayoutModule::LayoutModule(QObject* parent)
    : QObject(parent),
      m_status(new StatusComponent(this)),
      m_all_elements_model(new db::AllElementsModel(this)),
      m_root_elements_model(new db::RootElementsModel(this)),
      m_filtered_root_elements_model(new db::InstanceSortFilter(this)),
      m_child_model(new db::ChildModel(this)),
      m_filtered_child_model(new db::InstanceSortFilter(this)),
      m_breadcrumb_model(new db::BreadcrumbModel(this)),
      m_instance_edit(new db::AnInstanceEditor(this)),
      m_world_geometry_model(new db::WorldGeometryModel(this)) {

    m_filtered_root_elements_model->setSourceModel(m_root_elements_model);
    m_filtered_child_model->setSourceModel(m_child_model);

    connect(this,
            &LayoutModule::current_database_value_changed,
            m_all_elements_model,
            &db::AllElementsModel::reset);

    connect(this,
            &LayoutModule::current_database_value_changed,
            m_root_elements_model,
            &db::RootElementsModel::reset);

    connect(this,
            &LayoutModule::current_database_value_changed,
            this,
            &LayoutModule::reset);

    connect(this,
            &LayoutModule::current_database_value_changed,
            m_filtered_root_elements_model,
            &db::InstanceSortFilter::reset);

    connect(this,
            &LayoutModule::current_database_value_changed,
            m_filtered_child_model,
            &db::InstanceSortFilter::reset);

    connect(this,
            &LayoutModule::current_database_value_changed,
            m_child_model,
            &db::ChildModel::reset);

    connect(this,
            &LayoutModule::current_database_value_changed,
            m_breadcrumb_model,
            &db::BreadcrumbModel::reset);

    connect(this,
            &LayoutModule::current_database_value_changed,
            m_instance_edit,
            &db::AnInstanceEditor::reset);

    connect(m_instance_edit,
            &db::AnInstanceEditor::notify,
            this,
            &LayoutModule::notify);

    connect(this,
            &LayoutModule::current_database_value_changed,
            m_world_geometry_model,
            &db::WorldGeometryModel::reset);

    // Element changes

    connect(this,
            &LayoutModule::viewed_element_changed,
            this,
            &LayoutModule::viewed_entity_changed);

    connect(this,
            &LayoutModule::edited_element_changed,
            this,
            &LayoutModule::edited_entity_changed);
}

} // namespace SolTrace::GUI::App

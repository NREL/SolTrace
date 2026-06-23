#include "materials_module.h"

namespace SolTrace::GUI::App {

MaterialsModule::MaterialsModule(QObject* parent)
    : QObject(parent),
      m_status(new StatusComponent(this)),
      m_materials_list(new db::MaterialGroupsModel(this)),
      m_geometry_list(new db::GeometryGroupsModel(this)),
      m_material_edit(new db::MaterialEditor(this)),
      m_geometry_edit(new db::GeometryEditor(this)) {

    connect(this,
            &MaterialsModule::current_database_value_changed,
            this,
            &MaterialsModule::reset);

    connect(this,
            &MaterialsModule::current_database_value_changed,
            m_materials_list,
            &db::MaterialGroupsModel::reset);

    connect(this,
            &MaterialsModule::current_database_value_changed,
            m_geometry_list,
            &db::GeometryGroupsModel::reset);

    connect(this,
            &MaterialsModule::current_material_changed,
            this,
            &MaterialsModule::new_material_selected);

    connect(this,
            &MaterialsModule::current_geometry_changed,
            this,
            &MaterialsModule::new_geometry_selected);
}

void MaterialsModule::new_material_selected() {
    if (!m_current_database) {
        set_current_material_name(QString());
        return;
    }

    material_edit()->set(m_current_database, m_current_material);

    set_current_material_name(m_current_database->name_of(m_current_material));
}

void MaterialsModule::new_geometry_selected() {
    if (!m_current_database) {
        set_current_geometry_name(QString());
        return;
    }

    geometry_edit()->set(m_current_database, m_current_geometry);

    set_current_geometry_name(m_current_database->name_of(m_current_geometry));
}

void MaterialsModule::reset(db::Database* db) {
    if (m_database) {
        disconnect(m_database->identity.self(), nullptr, this, nullptr);
        disconnect(m_database, nullptr, this, nullptr);
    }

    m_database = db;

    set_current_geometry({});
    set_current_geometry_name({});
    set_current_material({});
    set_current_material_name({});

    m_material_edit->set(db, entt::null);
    m_geometry_edit->set(db, entt::null);

    if (!db) return;

    connect(db->identity.self(),
            &db::ComponentAPIBase::changed,
            this,
            [this](entt::entity e) {
                if (db::Entity(e) == m_current_material) {
                    set_current_material_name(
                        m_current_database->name_of(m_current_material));
                }
                if (db::Entity(e) == m_current_geometry) {
                    set_current_geometry_name(
                        m_current_database->name_of(m_current_geometry));
                }
            });
}

} // namespace SolTrace::GUI::App

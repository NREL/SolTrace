#pragma once

#include "data/database.h"
#include "editing/geometryeditor.h"
#include "editing/materialeditor.h"
#include "scene_models/group_models.h"
#include "support/qt_helpers.h"

#include <QObject>

namespace SolTrace::GUI::App {

/**
 * @class MaterialsModule
 * @brief Materials configuration module.
 *
 * Provides QML access to material and geometry group lists plus the editors for
 * the selected groups in the active GUI database.
 *
 * QML access pattern: App.materials.material_edit
 */
class MaterialsModule : public QObject {
    Q_OBJECT

    // TODO: A name module that always watches the name of an entity

    QPointer<SolTrace::GUI::Data::Database> m_database;

private slots:
    void new_material_selected();
    void new_geometry_selected();

    void reset(SolTrace::GUI::Data::Database*);

public:
    explicit MaterialsModule(QObject* parent = nullptr);

    QOBJECT_WRITABLE_PROPERTY(SolTrace::GUI::Data::Database, current_database)
    QOBJECT_WRITABLE_PROPERTY(SolTrace::GUI::Data::MaterialGroupsModel, materials_list)
    QOBJECT_WRITABLE_PROPERTY(SolTrace::GUI::Data::GeometryGroupsModel, geometry_list)

    QOBJECT_WRITABLE_PROPERTY(SolTrace::GUI::Data::MaterialEditor, material_edit);
    QOBJECT_WRITABLE_PROPERTY(SolTrace::GUI::Data::GeometryEditor, geometry_edit);

    Q_WRITABLE_PROPERTY(SolTrace::GUI::Data::Entity, current_material, { })
    Q_READONLY_PROPERTY(QString, current_material_name)

    Q_WRITABLE_PROPERTY(SolTrace::GUI::Data::Entity, current_geometry, { })
    Q_READONLY_PROPERTY(QString, current_geometry_name)
};

} // namespace SolTrace::GUI::App

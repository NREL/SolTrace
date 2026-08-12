#pragma once

#include "database/database.h"
#include "database/geometryeditor.h"
#include "database/materialeditor.h"
#include "database/models/group_models.h"
#include "module_common.h"
#include "utilities/qt_helpers.h"

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

    QPointer<db::Database> m_database;

private slots:
    void new_material_selected();
    void new_geometry_selected();

    void reset(db::Database*);

public:
    explicit MaterialsModule(QObject* parent = nullptr);

    QOBJECT_READONLY_PROPERTY(StatusComponent, status);
    QOBJECT_WRITABLE_PROPERTY(db::Database, current_database)
    QOBJECT_WRITABLE_PROPERTY(db::MaterialGroupsModel, materials_list)
    QOBJECT_WRITABLE_PROPERTY(db::GeometryGroupsModel, geometry_list)

    QOBJECT_WRITABLE_PROPERTY(db::MaterialEditor, material_edit);
    QOBJECT_WRITABLE_PROPERTY(db::GeometryEditor, geometry_edit);

    Q_WRITABLE_PROPERTY(db::Entity, current_material, { })
    Q_READONLY_PROPERTY(QString, current_material_name)

    Q_WRITABLE_PROPERTY(db::Entity, current_geometry, { })
    Q_READONLY_PROPERTY(QString, current_geometry_name)
};

} // namespace SolTrace::GUI::App

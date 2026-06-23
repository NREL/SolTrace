#pragma once

#include "backend.h"
#include "database/rootelementsmodel.h"
#include "module_common.h"
#include "utilities/notification.h"
#include "utilities/qt_helpers.h"
#include <QObject>

namespace SolTrace::GUI::App {


/**
 * @class GeometryModule
 * @brief Geometry configuration module.
 *
 * Provides QML access to the geometry database models owned by GeometryBackend.
 * Holds a non-owning QPointer reference to its backend slice.
 *
 * QML access pattern: App.geometry.backend.world_geometry_model
 */
class LayoutModule : public QObject {
    Q_OBJECT

    QPointer<db::Database> m_observed_database;

private slots:
    void viewed_entity_changed();
    void edited_entity_changed();
    void reset(db::Database*);
    void identity_changed(entt::entity);

public:
    explicit LayoutModule(QObject* parent = nullptr);

    QOBJECT_WRITABLE_PROPERTY(db::Database, current_database)

    QOBJECT_READONLY_PROPERTY(db::AllElementsModel, all_elements_model);

    QOBJECT_READONLY_PROPERTY(db::RootElementsModel, root_elements_model);
    QOBJECT_READONLY_PROPERTY(db::InstanceSortFilter,
                              filtered_root_elements_model);

    QOBJECT_WRITABLE_PROPERTY(db::ChildModel, child_model);
    QOBJECT_READONLY_PROPERTY(db::InstanceSortFilter, filtered_child_model);

    QOBJECT_WRITABLE_PROPERTY(db::BreadcrumbModel, breadcrumb_model);
    QOBJECT_WRITABLE_PROPERTY(db::AnInstanceEditor, instance_edit);
    QOBJECT_READONLY_PROPERTY(db::WorldGeometryModel, world_geometry_model);

    Q_WRITABLE_PROPERTY(db::Entity, viewed_element, { })
    Q_WRITABLE_PROPERTY(db::Entity, edited_element, { })
    Q_READONLY_PROPERTY(QString, edited_element_name);

    /// we need a selected element lists. need global pos and rot

    QOBJECT_READONLY_PROPERTY(StatusComponent, status);

public slots:
    void clear_viewed_element() { set_viewed_element({}); }
    void clear_edited_element() { set_edited_element({}); }
    void delete_edited_element();

signals:
    void notify(ANotification);
};


} // namespace SolTrace::GUI::App

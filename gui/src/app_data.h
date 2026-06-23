#pragma once
#include <QAbstractListModel>
#include <QHash>
#include <QJsonObject>
#include <QObject>
#include <QQmlEngine>
#include <QSharedPointer>

#include "logging.h"
#include "utilities/notification.h"
#include "utilities/qt_helpers.h"

#include <backend.h>

#include <modules/database_module.h>
#include <modules/documentation_module.h>
#include <modules/export_module.h>
#include <modules/flux_module.h>
#include <modules/intersections_module.h>
#include <modules/layout_module.h>
#include <modules/materials_module.h>
#include <modules/module_common.h>
#include <modules/simulation_module.h>
#include <modules/sun_module.h>
#include <modules/view_module.h>
#include <script/script.h>

/**
 * @namespace SolTrace::GUI::AppData
 * @brief AppDatalication layer between QML and the simulation backend.
 *
 * This namespace defines the GUI AppDatalication layer for SolTrace. It
 * mediates between the QML presentation layer and the simulation backend,
 * providing:
 *
 * - Domain-driven module decomposition (Sun, Tracing, Materials, Geometry,
 * etc.)
 * - Status lifecycle tracking per module
 * - Preset management for user-configurable parameters
 * - Inline documentation loaded from markdown files
 * - Non-owning references to backend services via QPointer
 *
 * Architecture:
 * @code
 *   QML  →  AppData (singleton)  →  Domain Modules  →  Backend (singleton)
 * @endcode
 *
 * Each domain module holds a QPointer to its corresponding backend slice,
 * constraining access and making dependencies explicit.
 */

namespace SolTrace::GUI::App {
/**
 * @class AppData
 * @brief QML singleton — top-level entrypoint for the AppDatalication layer.
 *
 * Provides a single, stable access point for all AppDatalication modules.
 * Registered as a QML singleton so all components share one instance.
 *
 * ## Initialization
 * After construction, wire backend references by calling the AppDataropriate
 * install methods on each module before the QML engine loads.
 *
 * @code
 *   // main.cpp
 *   auto* AppData = AppData::instance();
 *   AppData->docs().load();
 *   AppData->sun().backend() = QPointer(backend->sun_backend());
 *   // etc.
 * @endcode
 *
 */
class AppData : public QObject {
    Q_OBJECT
    QML_ELEMENT
    QML_SINGLETON

    void load_session();
    void save_session();
    void clear_session();

public:
    Q_READONLY_PROPERTY(QString, current_version_info);
    Q_READONLY_PROPERTY(QString, current_build_info);
    Q_READONLY_PROPERTY(bool, is_prerelease);

    static AppData* create(QQmlEngine* qmlEngine, QJSEngine*);

    explicit AppData(QObject*       parent,
                     QQmlEngine*    engine,
                     QString const& documentation_directory = "");

    ~AppData();

    Q_INVOKABLE void copy_build_info_to_clipboard();

    QOBJECT_WRITABLE_PROPERTY(db::Database, current_database)

    QOBJECT_WRITABLE_PROPERTY(LogList, log_list)

    QOBJECT_READONLY_PROPERTY(DatabaseModule, file_source)

    QOBJECT_READONLY_PROPERTY(DocumentationModule, docs)

    QOBJECT_READONLY_PROPERTY(SunModule, sun)

    QOBJECT_READONLY_PROPERTY(MaterialsModule, materials)

    QOBJECT_READONLY_PROPERTY(LayoutModule, layout)

    QOBJECT_READONLY_PROPERTY(ViewModule, view)

    QOBJECT_READONLY_PROPERTY(SimulationModule, simulation)

    QOBJECT_READONLY_PROPERTY(IntersectionsModule, intersections)

    QOBJECT_READONLY_PROPERTY(FluxModule, flux)

    QOBJECT_READONLY_PROPERTY(ExportModule, exporter)

    QOBJECT_READONLY_PROPERTY(Script::Script, script)

signals:
    void notification(ANotification);
    void new_results(db::SimulationResultPtr);
    void new_database(db::Database*);
};

} // namespace SolTrace::GUI::App

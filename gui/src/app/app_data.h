#pragma once
#include <QAbstractListModel>
#include <QHash>
#include <QJsonObject>
#include <QObject>
#include <QQmlEngine>
#include <QSharedPointer>
#include <QTranslator>

#include "app/logging.h"
#include "support/notification.h"
#include "support/qt_helpers.h"

#include <modules/database_module.h>
#include <modules/documentation_module.h>
#include <modules/export_module.h>
#include <modules/flux_module.h>
#include <modules/intersections_module.h>
#include <modules/layout_module.h>
#include <modules/materials_module.h>
#include <modules/simulation_module.h>
#include <modules/sun_module.h>
#include <modules/view_module.h>
#include <scripting/script.h>

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
 * - Inline documentation loaded from markdown files
 * - Session-level application settings
 *
 * Architecture:
 * @code
 *   QML  →  AppData (singleton)  →  Domain Modules  →  Backend (singleton)
 * @endcode
 *
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

    void apply_ui_locale(DocumentationModule::Locale locale);

public:
    Q_READONLY_PROPERTY(QString, current_version_info);
    Q_READONLY_PROPERTY(QString, current_build_info);
    Q_READONLY_PROPERTY(bool, is_prerelease);
    Q_WRITABLE_PROPERTY(bool, show_get_started_on_start, true);

    static AppData* create(QQmlEngine* qmlEngine, QJSEngine*);

    explicit AppData(QObject*       parent,
                     QQmlEngine*    engine,
                     QString const& documentation_directory = "");

    ~AppData();

    Q_INVOKABLE void copy_build_info_to_clipboard();
    Q_INVOKABLE void reset_all_settings();

    QOBJECT_WRITABLE_PROPERTY(SolTrace::GUI::Data::Database, current_database)

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

    QOBJECT_READONLY_PROPERTY(Scripting::Script, script)

private:
    QPointer<QQmlEngine> m_engine;
    QTranslator          m_ui_translator;
    bool                 m_ui_translator_installed = false;

signals:
    void notification(Support::ANotification);
    void new_results(SolTrace::GUI::Data::SimulationResultPtr);
    void new_database(SolTrace::GUI::Data::Database*);
};

} // namespace SolTrace::GUI::App

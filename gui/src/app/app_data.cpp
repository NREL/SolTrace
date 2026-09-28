#include "app/app_data.h"

#include "app_build_info.h"
#include "app/app_settings.h"

#include <QApplication>
#include <QClipboard>
#include <QGuiApplication>
#include <QLocale>

namespace SolTrace::GUI::App {

static QString build_info_string() {
    auto tag = BuildInfo::git_tag[0] == '\0'
                   ? QStringLiteral("(none)")
                   : QString::fromUtf8(BuildInfo::git_tag);

    return QStringLiteral("SolTrace build\n"
                          "Version: %1\n"
                          "Commit: %2\n"
                          "Describe: %3\n"
                          "Branch: %4\n"
                          "Tag: %5\n"
                          "Dirty: %6\n"
                          "Prerelease: %7")
        .arg(QString::fromUtf8(BuildInfo::version),
             QString::fromUtf8(BuildInfo::git_commit),
             QString::fromUtf8(BuildInfo::git_describe),
             QString::fromUtf8(BuildInfo::git_branch),
             tag,
             QString::fromUtf8(BuildInfo::git_dirty),
             BuildInfo::is_prerelease ? QStringLiteral("true")
                                      : QStringLiteral("false"));
}

void AppData::reset_all_settings() {
    AppSettings::clear_session();
    AppSettings::load_session(*this);
    emit notification(
        Support::ANotification::info("Settings reset to defaults."));
}

void AppData::apply_ui_locale(DocumentationModule::Locale locale) {
    if (m_ui_translator_installed) {
        qApp->removeTranslator(&m_ui_translator);
        m_ui_translator_installed = false;
    }

    switch (locale) {
    case DocumentationModule::Locale::EN:
        QLocale::setDefault(QLocale(QLocale::English));
        break;
    case DocumentationModule::Locale::ES:
        QLocale::setDefault(QLocale(QLocale::Spanish));
        if (m_ui_translator.load(QStringLiteral(":/i18n/soltrace_es.qm"))) {
            m_ui_translator_installed =
                qApp->installTranslator(&m_ui_translator);
        } else {
            qWarning() << "Unable to load UI translation"
                       << QStringLiteral(":/i18n/soltrace_es.qm");
        }
        break;
    }

    if (m_engine) { m_engine->retranslate(); }
}

AppData* AppData::create(QQmlEngine* qmlEngine, QJSEngine*) {
    return new AppData(nullptr, qmlEngine, "");
}

AppData::AppData(QObject*       parent,
                 QQmlEngine*    engine,
                 const QString& documentation_directory)
    : m_file_source(new DatabaseModule(this)),
      m_log_list(current_log_list()),
      m_view(new ViewModule(this)),
      m_docs(new DocumentationModule(this)),
      m_sun(new SunModule(this)),
      m_materials(new MaterialsModule(this)),
      m_layout(new LayoutModule(this)),
      m_simulation(new SimulationModule(this)),
      m_intersections(new IntersectionsModule(this)),
      m_flux(new FluxModule(engine, this)),
      m_exporter(new ExportModule(this)),
      m_script(new Scripting::Script(this)),
      m_engine(engine) {

    set_current_version_info(
        QString("%1 %2").arg(BuildInfo::version).arg(BuildInfo::git_commit));
    set_current_build_info(build_info_string());
    set_is_prerelease(BuildInfo::is_prerelease);


    connect(m_file_source,
            &DatabaseModule::current_database_value_changed,
            this,
            &AppData::set_current_database);

    connect(
        m_file_source, &DatabaseModule::notify, this, &AppData::notification);

    connect(
        m_simulation, &SimulationModule::notify, this, &AppData::notification);

    connect(m_layout, &LayoutModule::notify, this, &AppData::notification);

    connect(m_sun, &SunModule::notify, this, &AppData::notification);

    connect(m_flux, &FluxModule::notify, this, &AppData::notification);

    connect(m_exporter, &ExportModule::notify, this, &AppData::notification);

    connect(m_script, &Scripting::Script::notify, this, &AppData::notification);

    connect(m_docs, &DocumentationModule::locale_changed, this, [this] {
        apply_ui_locale(m_docs->locale());
    });

    connect(this,
            &AppData::current_database_value_changed,
            this,
            &AppData::new_database);

    connect(this,
            &AppData::new_database,
            m_simulation,
            &SimulationModule::set_current_database);

    connect(
        this, &AppData::new_database, m_sun, &SunModule::set_current_database);

    connect(this,
            &AppData::new_database,
            m_materials,
            &MaterialsModule::set_current_database);

    connect(this,
            &AppData::new_database,
            m_layout,
            &LayoutModule::set_current_database);

    connect(m_simulation,
            &SimulationModule::new_results,
            this,
            &AppData::new_results);

    connect(
        this, &AppData::new_database, m_script, &Scripting::Script::set_database);

    m_script->set_services(m_file_source, m_simulation);

    connect(qApp, &QCoreApplication::aboutToQuit, this, [this] {
        AppSettings::save_session(*this);
    });

    connect(this,
            &AppData::new_results,
            m_intersections,
            &IntersectionsModule::set_results);

    connect(this, &AppData::new_results, m_flux, &FluxModule::set_results);

    connect(
        this, &AppData::new_results, m_exporter, &ExportModule::set_results);

    connect(m_flux->pending_flux_maps(),
            &SolTrace::GUI::Data::PendingFluxMapModel::ready,
            m_exporter,
            &ExportModule::cache_flux_map);

    connect(m_simulation,
            &SimulationModule::edit_result_copy_requested,
            this,
            [this](SolTrace::GUI::Data::SimulationResultPtr result) {
                if (!m_file_source->append_clone(result)) return;

                m_view->set_workflow_phase(ViewModule::WorkflowPhase::Simulate);
                m_view->set_simulation_content_view(false);
            });

    AppSettings::load_session(*this);
    apply_ui_locale(m_docs->locale());

    m_file_source->load_new();
}

void AppData::copy_build_info_to_clipboard() {
    QGuiApplication::clipboard()->setText(current_build_info());
    emit notification(
        Support::ANotification::info("Copied build info to clipboard."));
}

AppData::~AppData() {
    AppSettings::save_session(*this);
}

} // namespace SolTrace::GUI::App

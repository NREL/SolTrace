#include "app_data.h"

#include "app_build_info.h"

#include <QApplication>
#include <QClipboard>
#include <QGuiApplication>
#include <QLocale>
#include <QSettings>

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

static DocumentationModule::Locale locale_from_setting(QVariant const& value) {
    auto const locale = value.toInt();
    if (locale == static_cast<int>(DocumentationModule::Locale::ES)) {
        return DocumentationModule::Locale::ES;
    }

    return DocumentationModule::Locale::EN;
}

void AppData::load_session() {
    QSettings s;

    s.beginGroup("Language");
    m_docs->set_locale(locale_from_setting(
        s.value("locale", static_cast<int>(DocumentationModule::Locale::EN))));
    s.endGroup();

    s.beginGroup("View");
    m_view->left_panel()->set_visible(
        s.value("show_left_panel", true).toBool());
    m_view->left_panel()->set_saved_visible(
        s.value("show_left_panel_saved", false).toBool());
    m_view->right_panel()->set_visible(
        s.value("show_right_panel", true).toBool());
    m_view->right_panel()->set_saved_visible(
        s.value("show_right_panel_saved", false).toBool());
    m_view->full_panel()->set_visible(
        s.value("show_settings_panel", false).toBool());
    m_view->left_panel()->update_size();
    m_view->right_panel()->update_size();

    m_view->left_panel()->set_width(s.value("left_panel_width", 550).toUInt());
    m_view->right_panel()->set_width(
        s.value("right_panel_width", 100).toUInt());

    m_view->set_inline_docs(
        s.value("inline_docs",
                s.value("left_panel_inline_docs", false).toBool() ||
                    s.value("right_panel_inline_docs", false).toBool())
            .toBool());

    m_view->set_workflow_phase(static_cast<ViewModule::WorkflowPhase>(
        s.value("workflow_phase", 0).toUInt()));

    m_view->set_configure_section(s.value("configure_section", 0).toUInt());
    m_view->set_simulate_section(s.value("simulate_section", 0).toUInt());
    m_view->set_analyze_section(s.value("analyze_section", 0).toUInt());

    m_view->set_sun_section(s.value("sun_section", 0).toUInt());
    m_view->set_right_panel_section(s.value("right_panel_section", 0).toUInt());

    m_view->full_panel()->set_mode(static_cast<FullPanelData::FullPanelMode>(
        s.value("full_panel_mode", 0).toUInt()));
    m_view->full_panel()->set_settings_section(
        s.value("settings_section", 0).toUInt());
    m_view->full_panel()->set_docs_section(s.value("docs_section", 0).toUInt());
    m_view->full_panel()->set_build_section(
        s.value("build_section", 0).toUInt());

    // Viewport
    auto* sim = m_view->sim();
    sim->set_camera(static_cast<SimulationViewState::Camera>(
        s.value("sim_camera", 1).toInt()));
    sim->set_perspective(static_cast<SimulationViewState::Perspective>(
        s.value("sim_perspective", 0).toInt()));
    sim->set_sun_viz(s.value("sim_sun_viz", true).toBool());
    sim->set_sun_viz_scale(s.value("sim_sun_viz_scale", 50.0).toDouble());
    sim->set_sun_color(
        s.value("sim_sun_color", QColor("yellow")).value<QColor>());
    sim->set_geometry_color(
        s.value("sim_geometry_color", QColor("white")).value<QColor>());
    m_view->set_ray_color(
        s.value("ray_color", QColor("white")).value<QColor>());
    sim->set_show_grid(s.value("sim_show_grid", true).toBool());
    s.endGroup();

    s.beginGroup("Sun");

    // Sun Position
    m_sun->ds_position()->set_from_calculator(
        s.value("ds_from_calculator", true).toBool());
    m_sun->ds_position()->set_x(s.value("ds_position_x", 1000.0).toDouble());
    m_sun->ds_position()->set_y(s.value("ds_position_y", 1000.0).toDouble());
    m_sun->ds_position()->set_z(s.value("ds_position_z", 1000.0).toDouble());

    m_sun->ps_position()->set_from_calculator(
        s.value("ps_from_calculator", true).toBool());
    m_sun->ps_position()->set_x(s.value("ps_position_x", 1000.0).toDouble());
    m_sun->ps_position()->set_y(s.value("ps_position_y", 1000.0).toDouble());
    m_sun->ps_position()->set_z(s.value("ps_position_z", 1000.0).toDouble());

    m_sun->set_type(static_cast<SunModule::Type>(s.value("type", 0).toInt()));

    // Calculator
    auto* calc = m_sun->calc_data();
    calc->set_calculator(static_cast<SolarCalculatorData::Calculator>(
        s.value("calculator", 0).toInt()));
    calc->set_latitude(s.value("latitude", 35.04).toDouble());
    calc->set_longitude(s.value("longitude", -105.10).toDouble());
    calc->set_year(s.value("year", 2026).toInt());
    calc->set_month(s.value("month", 3).toInt());
    calc->set_day(s.value("day", 20).toInt());
    calc->set_hour(s.value("hour", 12).toInt());
    calc->set_minute(s.value("minute", 0).toInt());
    calc->set_second(s.value("second", 0).toInt());
    calc->set_timezone_offset(s.value("timezone_offset", -7).toInt());
    calc->set_altitude(s.value("altitude", 1000).toDouble());
    calc->set_pressure(s.value("pressure", 1013.25).toDouble());
    calc->set_temperature(s.value("temperature", 20.0).toDouble());

    // Sun Shape
    m_sun->shape()->set_shape(
        static_cast<SunShape::Shape>(s.value("shape", 0).toDouble()));
    m_sun->shape()->set_sigma(s.value("sigma", 4.65).toDouble());
    m_sun->shape()->set_half_width(s.value("half_width", 4.65).toDouble());
    m_sun->shape()->set_csr(s.value("buie_csr", 0.1).toDouble());

    auto* cdist = m_sun->shape()->custom_distribution();
    if (s.contains("custom_shape"))
        cdist->set_variant_data(s.value("custom_shape").toList());
    else
        m_sun->shape()->reset_current_distribution();
    s.endGroup();
}

void AppData::save_session() {
    QSettings s;
    s.beginGroup("Language");
    s.setValue("locale", static_cast<int>(m_docs->locale()));
    s.endGroup();

    s.beginGroup("View");
    s.setValue("show_left_panel", m_view->left_panel()->visible());
    s.setValue("show_right_panel", m_view->right_panel()->visible());
    s.setValue("show_left_panel_saved", m_view->left_panel()->saved_visible());
    s.setValue("show_right_panel_saved",
               m_view->right_panel()->saved_visible());
    s.setValue("show_settings_panel", m_view->full_panel()->visible());

    s.setValue("left_panel_width", m_view->left_panel()->width());
    s.setValue("right_panel_width", m_view->right_panel()->width());

    s.setValue("inline_docs", m_view->inline_docs());

    s.setValue("workflow_phase", m_view->workflow_phase());

    s.setValue("configure_section", m_view->configure_section());
    s.setValue("simulate_section", m_view->simulate_section());
    s.setValue("analyze_section", m_view->analyze_section());

    s.setValue("sun_section", m_view->sun_section());
    s.setValue("right_panel_section", m_view->right_panel_section());

    s.setValue("full_panel_mode",
               static_cast<int>(m_view->full_panel()->mode()));
    s.setValue("settings_section", m_view->full_panel()->settings_section());
    s.setValue("docs_section", m_view->full_panel()->docs_section());
    s.setValue("build_section", m_view->full_panel()->build_section());

    // Viewport
    auto* sim = m_view->sim();
    s.setValue("sim_camera", static_cast<int>(sim->camera()));
    s.setValue("sim_perspective", static_cast<int>(sim->perspective()));
    s.setValue("sim_sun_viz", sim->sun_viz());
    s.setValue("sim_sun_viz_scale", sim->sun_viz_scale());
    s.setValue("sim_sun_color", sim->sun_color());
    s.setValue("sim_geometry_color", sim->geometry_color());
    s.setValue("ray_color", m_view->ray_color());
    s.setValue("sim_show_grid", sim->show_grid());

    s.endGroup();

    s.beginGroup("Sun");

    // Sun Position

    s.setValue("ds_from_calculator", m_sun->ds_position()->from_calculator());
    s.setValue("ds_position_x", m_sun->ds_position()->x());
    s.setValue("ds_position_y", m_sun->ds_position()->y());
    s.setValue("ds_position_z", m_sun->ds_position()->z());

    s.setValue("ps_from_calculator", m_sun->ps_position()->from_calculator());
    s.setValue("ps_position_x", m_sun->ps_position()->x());
    s.setValue("ps_position_y", m_sun->ps_position()->y());
    s.setValue("ps_position_z", m_sun->ps_position()->z());

    s.setValue("type", static_cast<int>(m_sun->type()));

    // Calculator
    auto* calc = m_sun->calc_data();
    s.setValue("calculator", static_cast<int>(calc->calculator()));
    s.setValue("latitude", calc->latitude());
    s.setValue("longitude", calc->longitude());
    s.setValue("year", calc->year());
    s.setValue("month", calc->month());
    s.setValue("day", calc->day());
    s.setValue("hour", calc->hour());
    s.setValue("minute", calc->minute());
    s.setValue("second", calc->second());
    s.setValue("timezone_offset", calc->timezone_offset());
    s.setValue("altitude", calc->altitude());
    s.setValue("pressure", calc->pressure());
    s.setValue("temperature", calc->temperature());

    // Sun Shape
    s.setValue("shape", static_cast<int>(m_sun->shape()->shape()));
    s.setValue("sigma", m_sun->shape()->sigma());
    s.setValue("half_width", m_sun->shape()->half_width());
    s.setValue("buie_csr", m_sun->shape()->csr());

    auto* cdist = m_sun->shape()->custom_distribution();
    s.setValue("custom_shape", cdist->variant_data());
    s.endGroup();
}

void AppData::clear_session() {
    QSettings s;
    s.clear();
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
      m_script(new Script::Script(this)),
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

    connect(m_script, &Script::Script::notify, this, &AppData::notification);

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
        this, &AppData::new_database, m_script, &Script::Script::set_database);

    m_script->set_services(m_file_source, m_simulation);

    connect(qApp, &QCoreApplication::aboutToQuit, this, &AppData::save_session);

    connect(this,
            &AppData::new_results,
            m_intersections,
            &IntersectionsModule::set_results);

    connect(this, &AppData::new_results, m_flux, &FluxModule::set_results);

    connect(
        this, &AppData::new_results, m_exporter, &ExportModule::set_results);

    connect(m_flux->pending_flux_maps(),
            &db::PendingFluxMapModel::ready,
            m_exporter,
            &ExportModule::cache_flux_map);

    connect(m_simulation,
            &SimulationModule::edit_result_copy_requested,
            this,
            [this](db::SimulationResultPtr result) {
                if (!m_file_source->append_clone(result)) return;

                m_view->set_workflow_phase(ViewModule::WorkflowPhase::Simulate);
                m_view->set_simulation_content_view(false);
            });

    load_session();
    apply_ui_locale(m_docs->locale());

    m_file_source->load_new();
}

void AppData::copy_build_info_to_clipboard() {
    QGuiApplication::clipboard()->setText(current_build_info());
    emit notification(ANotification::info("Copied build info to clipboard."));
}

AppData::~AppData() {
    save_session();
}

} // namespace SolTrace::GUI::App

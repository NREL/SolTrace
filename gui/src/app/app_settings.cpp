#include "app/app_settings.h"

#include "app/app_data.h"

#include <QSettings>

#include <algorithm>

namespace SolTrace::GUI::App {

static DocumentationModule::Locale locale_from_setting(QVariant const& value) {
    auto const locale = value.toInt();
    if (locale == static_cast<int>(DocumentationModule::Locale::ES)) {
        return DocumentationModule::Locale::ES;
    }

    return DocumentationModule::Locale::EN;
}

void AppSettings::load_session(AppData& app) {
    QSettings s;

    auto* docs = app.docs();
    auto* view = app.view();
    auto* sun  = app.sun();

    s.beginGroup("Language");
    docs->set_locale(locale_from_setting(
        s.value("locale", static_cast<int>(DocumentationModule::Locale::EN))));
    s.endGroup();

    s.beginGroup("View");
    view->left_panel()->set_visible(
        s.value("show_left_panel", true).toBool());
    view->right_panel()->set_visible(
        s.value("show_right_panel", true).toBool());
    view->left_panel()->update_size();
    view->right_panel()->update_size();

    view->left_panel()->set_width(s.value("left_panel_width", 550).toUInt());
    view->right_panel()->set_width(
        s.value("right_panel_width", 100).toUInt());

    view->set_inline_docs(
        s.value("inline_docs",
                s.value("left_panel_inline_docs", false).toBool() ||
                    s.value("right_panel_inline_docs", false).toBool())
            .toBool());

    app.set_show_get_started_on_start(
        s.value("show_get_started_on_start", true).toBool());

    auto workflow_phase =
        s.value("workflow_phase",
                static_cast<int>(ViewModule::WorkflowPhase::Load))
            .toInt();
    if (workflow_phase <= static_cast<int>(ViewModule::WorkflowPhase::Start) ||
        workflow_phase > static_cast<int>(ViewModule::WorkflowPhase::Analyze)) {
        workflow_phase = static_cast<int>(ViewModule::WorkflowPhase::Load);
    }
    view->set_workflow_phase(
        static_cast<ViewModule::WorkflowPhase>(workflow_phase));

    view->set_configure_section(s.value("configure_section", 0).toUInt());
    view->set_simulate_section(s.value("simulate_section", 0).toUInt());
    view->set_analyze_section(s.value("analyze_section", 0).toUInt());

    view->set_sun_section(s.value("sun_section", 0).toUInt());
    view->set_right_panel_section(s.value("right_panel_section", 0).toUInt());

    view->set_settings_section(
        std::clamp(s.value("settings_section", 0).toInt(), 0, 3));
    view->set_info_section(
        s.value("info_section", s.value("docs_section", 0)).toUInt());

    auto* sim = view->sim();
    sim->set_perspective(static_cast<SimulationViewState::Perspective>(
        s.value("sim_perspective", 0).toInt()));
    sim->set_sun_viz(s.value("sim_sun_viz", true).toBool());
    sim->set_sun_viz_scale(s.value("sim_sun_viz_scale", 50.0).toDouble());
    sim->set_sun_color(
        s.value("sim_sun_color", QColor("yellow")).value<QColor>());
    sim->set_geometry_color(
        s.value("sim_geometry_color", QColor("white")).value<QColor>());
    view->set_ray_color(s.value("ray_color", QColor("white")).value<QColor>());
    sim->set_show_grid(s.value("sim_show_grid", true).toBool());
    s.endGroup();

    s.beginGroup("Sun");

    sun->ds_position()->set_from_calculator(
        s.value("ds_from_calculator", true).toBool());
    sun->ds_position()->set_x(s.value("ds_position_x", 1000.0).toDouble());
    sun->ds_position()->set_y(s.value("ds_position_y", 1000.0).toDouble());
    sun->ds_position()->set_z(s.value("ds_position_z", 1000.0).toDouble());

    sun->ps_position()->set_from_calculator(
        s.value("ps_from_calculator", true).toBool());
    sun->ps_position()->set_x(s.value("ps_position_x", 1000.0).toDouble());
    sun->ps_position()->set_y(s.value("ps_position_y", 1000.0).toDouble());
    sun->ps_position()->set_z(s.value("ps_position_z", 1000.0).toDouble());

    sun->set_type(static_cast<SunModule::Type>(s.value("type", 0).toInt()));

    auto* calc = sun->calc_data();
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

    sun->shape()->set_shape(
        static_cast<SunShape::Shape>(s.value("shape", 0).toDouble()));
    sun->shape()->set_sigma(s.value("sigma", 4.65).toDouble());
    sun->shape()->set_half_width(s.value("half_width", 4.65).toDouble());
    sun->shape()->set_csr(s.value("buie_csr", 0.1).toDouble());

    auto* cdist = sun->shape()->custom_distribution();
    if (s.contains("custom_shape"))
        cdist->set_variant_data(s.value("custom_shape").toList());
    else
        sun->shape()->reset_current_distribution();
    s.endGroup();
}

void AppSettings::save_session(AppData const& app) {
    QSettings s;

    auto const* docs = app.docs();
    auto const* view = app.view();
    auto const* sun  = app.sun();

    s.beginGroup("Language");
    s.setValue("locale", static_cast<int>(docs->locale()));
    s.endGroup();

    s.beginGroup("View");
    s.setValue("show_left_panel", view->left_panel()->visible());
    s.setValue("show_right_panel", view->right_panel()->visible());
    s.setValue("left_panel_width", view->left_panel()->width());
    s.setValue("right_panel_width", view->right_panel()->width());

    s.setValue("inline_docs", view->inline_docs());
    s.setValue("show_get_started_on_start",
               app.show_get_started_on_start());

    auto workflow_phase = view->workflow_phase();
    if (workflow_phase == ViewModule::WorkflowPhase::Start) {
        workflow_phase = ViewModule::WorkflowPhase::Load;
    }
    s.setValue("workflow_phase", workflow_phase);

    s.setValue("configure_section", view->configure_section());
    s.setValue("simulate_section", view->simulate_section());
    s.setValue("analyze_section", view->analyze_section());

    s.setValue("sun_section", view->sun_section());
    s.setValue("right_panel_section", view->right_panel_section());

    s.setValue("settings_section", view->settings_section());
    s.setValue("info_section", view->info_section());

    auto* sim = view->sim();
    s.setValue("sim_perspective", static_cast<int>(sim->perspective()));
    s.setValue("sim_sun_viz", sim->sun_viz());
    s.setValue("sim_sun_viz_scale", sim->sun_viz_scale());
    s.setValue("sim_sun_color", sim->sun_color());
    s.setValue("sim_geometry_color", sim->geometry_color());
    s.setValue("ray_color", view->ray_color());
    s.setValue("sim_show_grid", sim->show_grid());
    s.endGroup();

    s.beginGroup("Sun");

    s.setValue("ds_from_calculator", sun->ds_position()->from_calculator());
    s.setValue("ds_position_x", sun->ds_position()->x());
    s.setValue("ds_position_y", sun->ds_position()->y());
    s.setValue("ds_position_z", sun->ds_position()->z());

    s.setValue("ps_from_calculator", sun->ps_position()->from_calculator());
    s.setValue("ps_position_x", sun->ps_position()->x());
    s.setValue("ps_position_y", sun->ps_position()->y());
    s.setValue("ps_position_z", sun->ps_position()->z());

    s.setValue("type", static_cast<int>(sun->type()));

    auto* calc = sun->calc_data();
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

    s.setValue("shape", static_cast<int>(sun->shape()->shape()));
    s.setValue("sigma", sun->shape()->sigma());
    s.setValue("half_width", sun->shape()->half_width());
    s.setValue("buie_csr", sun->shape()->csr());

    auto* cdist = sun->shape()->custom_distribution();
    s.setValue("custom_shape", cdist->variant_data());
    s.endGroup();
}

void AppSettings::clear_session() {
    QSettings s;
    s.clear();
}

} // namespace SolTrace::GUI::App

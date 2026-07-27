#include "sun_module.h"
#include <QClipboard>
#include <QDebug>
#include <QGuiApplication>
#include <QRegularExpression>
#include <QScopedValueRollback>
#include <algorithm>
#include <cmath>
#include <exception>
#include <numbers>
#include <optional>
#include <vector>

namespace SolTrace::GUI::App {

namespace {

QVector<SunShapePoint> normalized_radial_points(QVector<SunShapePoint> points) {
    for (auto& point : points) {
        point.angle = std::abs(point.angle);
    }

    std::sort(points.begin(), points.end(), [](auto const& a, auto const& b) {
        return a.angle < b.angle;
    });

    QVector<SunShapePoint> merged;
    for (auto const& point : points) {
        if (!merged.empty() &&
            std::abs(merged.back().angle - point.angle) < 1.0e-9) {
            merged.back().intensity =
                std::max(merged.back().intensity, point.intensity);
        } else {
            merged.push_back(point);
        }
    }

    return merged;
}

std::optional<SunShape::Shape> gui_shape_for_data_shape(Data::SunShape shape) {
    switch (shape) {
    case Data::SunShape::GAUSSIAN: return SunShape::Shape::Gaussian;
    case Data::SunShape::PILLBOX: return SunShape::Shape::Pillbox;
    case Data::SunShape::BUIE_CSR: return SunShape::Shape::Buie_CSR;
    case Data::SunShape::USER_DEFINED: return SunShape::Shape::Custom;
    case Data::SunShape::LIMBDARKENED: return SunShape::Shape::LimbDarkened;
    default: return std::nullopt;
    }
}

QVector<SunShapePoint>
points_from_user_data(std::vector<double> const& angles,
                      std::vector<double> const& intensities) {
    QVector<SunShapePoint> points;
    const auto             count = std::min(angles.size(), intensities.size());
    points.reserve(static_cast<qsizetype>(count));

    for (std::size_t i = 0; i < count; ++i) {
        points.push_back({
            .angle     = angles[i],
            .intensity = intensities[i],
        });
    }

    return normalized_radial_points(points);
}

} // namespace

SunModule::SunModule(QObject* parent)
    : QObject(parent),
      m_status(new StatusComponent(this)),
      m_shape(new SunShape(this)),
      m_ps_position(new SolarPositionData(this)),
      m_ds_position(new SolarPositionData(this)),
      m_calc_data(new SolarCalculatorData(this))

{
    connect(m_shape, &SunShape::changed, this, &SunModule::update_shape);

    /* connect(m_calc_data,
            &SolarCalculatorData::changed,
            this,
            &SunModule::update_position); */

    connect(this, &SunModule::type_changed, this, &SunModule::update_type);
    connect(this,
            &SunModule::current_database_changed,
            this,
            &SunModule::update_database_connections);
    connect(m_ps_position,
            &SolarPositionData::changed,
            this,
            &SunModule::write_position_to_database);
    connect(m_ds_position,
            &SolarPositionData::changed,
            this,
            &SunModule::write_position_to_database);

    update_type();
    update_position();
}


void SunModule::update_shape() {
    write_shape_to_database();
}

void SunModule::write_shape_to_database() {
    if (m_loading_from_database || m_writing_to_database) return;
    if (!m_current_database) return;

    qDebug() << Q_FUNC_INFO;

    QScopedValueRollback<bool> guard(m_writing_to_database, true);

    try {
        m_current_database->ray_source_resource.patch([this](auto& resource) {
            if (!resource.source) {
                resource.source = SD::make_ray_source<SD::Sun>();
            }

            resource.source->set_shape(
                m_shape->get_sunshape_data(),
                m_shape->sigma(),
                m_shape->half_width(),
                m_shape->csr(),
                m_shape->custom_distribution()->get_angle_data(),
                m_shape->custom_distribution()->get_intensity_data());
        });
    } catch (std::exception const& e) {
        qWarning() << "Unable to update sun shape:" << e.what();
        emit notify(
            ANotification::error(QString("Could not update the sun shape: %1")
                                     .arg(QString::fromUtf8(e.what()))));
    }
}

void SunModule::update_type() {
    if (m_type == Type::Directional) set_position(m_ds_position);
    else
        set_position(m_ps_position);

    write_position_to_database();
}

QString SunModule::update_position() {
    if (m_loading_from_database)
        return QStringLiteral("Currently loading from database");
    if (!m_position->from_calculator())
        return QStringLiteral("Not using calculator");

    qDebug() << Q_FUNC_INFO;

    try {
        m_calculator.set_method(selected_calculation_method());

        m_calculator.set_location(m_calc_data->latitude(),
                                  m_calc_data->longitude(),
                                  m_calc_data->timezone_offset(),
                                  m_calc_data->altitude());

        m_calculator.set_environment(m_calc_data->pressure(),
                                     m_calc_data->temperature());

        m_calculator.set_date(
            m_calc_data->year(), m_calc_data->month(), m_calc_data->day());
        m_calculator.set_time(
            m_calc_data->hour(), m_calc_data->minute(), m_calc_data->second());

        double x, y, z, azimuth, elevation;
        m_calculator.get_sun_vector(&x, &y, &z);
        m_calculator.get_azimuth_elevation(&azimuth, &elevation);
        {
            QScopedValueRollback<bool> guard(m_updating_calculated_position,
                                             true);
            m_position->set_x(x);
            m_position->set_y(y);
            m_position->set_z(z);
            m_position->set_azimuth(azimuth);
            m_position->set_elevation(elevation);
        }
        return write_position_to_database();
    } catch (std::exception const& e) {
        qWarning() << "Unable to calculate sun position:" << e.what();
        emit notify(ANotification::error(
            QString("Could not calculate the sun position: %1")
                .arg(QString::fromUtf8(e.what()))));
        return e.what();
    }

    return { };
}

QString SunModule::write_position_to_database() {
    if (m_loading_from_database || m_writing_to_database ||
        m_updating_calculated_position) {
        return { };
    }
    if (!m_current_database || !m_position) return { };

    qDebug() << Q_FUNC_INFO;

    QScopedValueRollback<bool> guard(m_writing_to_database, true);

    try {
        m_current_database->ray_source_resource.patch(
            [this](db::RaySourceResource& resource) {
                const bool created = !resource.source;
                if (created) {
                    resource.source = SD::make_ray_source<SD::Sun>();
                }

                if (created) {
                    resource.source->set_shape(
                        m_shape->get_sunshape_data(),
                        m_shape->sigma(),
                        m_shape->half_width(),
                        m_shape->csr(),
                        m_shape->custom_distribution()->get_angle_data(),
                        m_shape->custom_distribution()->get_intensity_data());
                }

                resource.source->set_position(
                    m_position->x(), m_position->y(), m_position->z());

                resource.source->set_gen_type(Data::GenType::RANDOM);
                resource.type = m_type == Type::PointSource
                                    ? db::RaySourceType::PointSource
                                    : db::RaySourceType::Directional;
            });
    } catch (std::exception const& e) {
        qWarning() << "Unable to update sun position:" << e.what();
        emit notify(ANotification::error(
            QString("Could not update the sun position: %1")
                .arg(QString::fromUtf8(e.what()))));

        return e.what();
    }

    return { };
}

void SunModule::update_database_connections() {
    for (auto const& connection : std::as_const(m_database_connections)) {
        QObject::disconnect(connection);
    }
    m_database_connections.clear();

    if (!m_current_database) return;

    m_database_connections.push_back(
        connect(m_current_database->ray_source_resource.self(),
                &db::ComponentAPIBase::changed,
                this,
                [this](entt::entity) {
                    if (m_writing_to_database) return;
                    load_from_database();
                }));

    m_database_connections.push_back(
        connect(m_current_database->ray_source_resource.self(),
                &db::ComponentAPIBase::removed,
                this,
                [this](entt::entity) {
                    if (m_writing_to_database) return;
                    load_from_database();
                }));

    load_from_database();
}

namespace {

constexpr double DEFAULT_SIGMA     = 4.65;
constexpr double DEFAULT_HALFWIDTH = 4.65;
constexpr double DEFAULT_CSR       = 0.1;

SD::Sun make_default_ray_source() {
    SD::Sun sun;
    sun.set_position(0.0, 0.0, 1.0);
    sun.set_shape(
        SD::SunShape::GAUSSIAN, DEFAULT_SIGMA, DEFAULT_HALFWIDTH, DEFAULT_CSR);
    sun.set_gen_type(SD::GenType::RANDOM);
    return sun;
}

SD::RaySource& default_ray_source() {
    // This is shared to avoid constructing a fallback on every scene change.
    // Treat the returned object as read-only. The mutable reference is required
    // because RaySource's parameter and user-data getters are not const.
    static auto ret = make_default_ray_source();
    return ret;
}

} // namespace

void SunModule::load_from_database() {
    qDebug() << Q_FUNC_INFO << m_current_database;
    if (!m_current_database) return;

    auto const* resource = m_current_database->ray_source_resource.get();
    auto const  source_type =
        resource ? resource->type : db::RaySourceType::Directional;

    if (!resource || !resource->source) {
        qDebug() << Q_FUNC_INFO << "No ray source in this database";
        load_from_ray_source(default_ray_source(), source_type);
        return;
    }

    load_from_ray_source(*resource->source, source_type);
}

void SunModule::load_from_ray_source(SD::RaySource&    ray_source,
                                     db::RaySourceType source_type) {
    qDebug() << Q_FUNC_INFO;
    // We always need a shape to work with here
    auto gui_shape = gui_shape_for_data_shape(ray_source.get_shape())
                         .value_or(SunShape::Shape::Gaussian);

    QScopedValueRollback<bool> guard(m_loading_from_database, true);

    set_type(source_type == db::RaySourceType::PointSource ? Type::PointSource
                                                           : Type::Directional);

    // A core Sun stores only the parameter used by its active shape and leaves
    // the other getters as NaN, so apply the complete GUI baseline explicitly.
    m_shape->set_sigma(DEFAULT_SIGMA);
    m_shape->set_half_width(DEFAULT_HALFWIDTH);
    m_shape->set_csr(DEFAULT_CSR);

    if (!std::isnan(ray_source.get_sigma())) {
        m_shape->set_sigma(ray_source.get_sigma());
    }

    if (!std::isnan(ray_source.get_half_width())) {
        m_shape->set_half_width(ray_source.get_half_width());
    }

    if (!std::isnan(ray_source.get_circumsolar_ratio())) {
        m_shape->set_csr(ray_source.get_circumsolar_ratio());
    }

    std::vector<double> user_angles;
    std::vector<double> user_intensities;

    ray_source.get_user_data(user_angles, user_intensities);

    m_shape->custom_distribution()->reset(
        points_from_user_data(user_angles, user_intensities));

    m_shape->set_shape(gui_shape);

    const auto& position = ray_source.get_position();
    const auto  length =
        std::sqrt(position.x * position.x + position.y * position.y +
                  position.z * position.z);

    double azimuth   = 0.0;
    double elevation = 0.0;
    if (length > 1.0e-12) {
        constexpr double radians_to_degrees = 180.0 / std::numbers::pi;
        azimuth = std::atan2(position.x, position.y) * radians_to_degrees;
        if (azimuth < 0.0) azimuth += 360.0;
        elevation = std::asin(std::clamp(position.z / length, -1.0, 1.0)) *
                    radians_to_degrees;
    }

    m_ds_position->set_from_calculator(false);
    m_ds_position->set_x(position.x);
    m_ds_position->set_y(position.y);
    m_ds_position->set_z(position.z);
    m_ds_position->set_azimuth(azimuth);
    m_ds_position->set_elevation(elevation);

    m_ps_position->set_from_calculator(false);
    m_ps_position->set_x(position.x);
    m_ps_position->set_y(position.y);
    m_ps_position->set_z(position.z);
    m_ps_position->set_azimuth(azimuth);
    m_ps_position->set_elevation(elevation);
}

Data::SolarPositionCalculationMethod
SunModule::selected_calculation_method() const {
    switch (m_calc_data->calculator()) {
    case SolarCalculatorData::Calculator::Legacy:
        return Data::SolarPositionCalculationMethod::LEGACY;
    case SolarCalculatorData::Calculator::Duffie:
        return Data::SolarPositionCalculationMethod::DUFFIE;
    case SolarCalculatorData::Calculator::SOLPOS:
        return Data::SolarPositionCalculationMethod::SOLPOS;
    case SolarCalculatorData::Calculator::SPA:
        return Data::SolarPositionCalculationMethod::SPA;
    }

    return Data::SolarPositionCalculationMethod::LEGACY;
}

QString SunModule::apply_calculator(int    calculator,
                                    double latitude,
                                    double longitude,
                                    int    year,
                                    int    month,
                                    int    day,
                                    int    hour,
                                    int    minute,
                                    int    second,
                                    int    timezone_offset,
                                    double altitude,
                                    double pressure,
                                    double temperature) {
    {
        QScopedValueRollback<bool> guard(m_loading_from_database, true);

        const auto clamped_calculator = std::clamp(calculator, 0, 3);
        m_calc_data->set_calculator(
            static_cast<SolarCalculatorData::Calculator>(clamped_calculator));
        m_calc_data->set_latitude(latitude);
        m_calc_data->set_longitude(longitude);
        m_calc_data->set_year(year);
        m_calc_data->set_month(month);
        m_calc_data->set_day(day);
        m_calc_data->set_hour(hour);
        m_calc_data->set_minute(minute);
        m_calc_data->set_second(second);
        m_calc_data->set_timezone_offset(timezone_offset);
        m_calc_data->set_altitude(altitude);
        m_calc_data->set_pressure(pressure);
        m_calc_data->set_temperature(temperature);
        m_ds_position->set_from_calculator(true);
    }

    if (m_type == Type::Directional) { set_position(m_ds_position); }

    return update_position();
}

SunShape::SunShape(QObject* parent)
    : QObject(parent),
      m_generated_distribution(new SunShapeModel(this)),
      m_custom_distribution(new SunShapeModel(this)) {

    // SunShape::shape_changed() -> SunShape::update_current_distribution()
    connect(this,
            &SunShape::shape_changed,
            this,
            &SunShape::update_current_distribution);

    // SunShape::*_changed() -> SunShape::regenerate()
    connect(this, &SunShape::shape_changed, this, &SunShape::regenerate);
    connect(this, &SunShape::sigma_changed, this, &SunShape::regenerate);
    connect(this, &SunShape::half_width_changed, this, &SunShape::regenerate);
    connect(this, &SunShape::csr_changed, this, &SunShape::regenerate);
    connect(m_custom_distribution,
            &SunShapeModel::changed,
            this,
            &SunShape::update_x_axis);

    // SunShape::*_changed() -> SunShape::changed()
    connect(this, &SunShape::shape_changed, this, &SunShape::changed);
    connect(this, &SunShape::sigma_changed, this, &SunShape::changed);
    connect(this, &SunShape::half_width_changed, this, &SunShape::changed);
    connect(this, &SunShape::csr_changed, this, &SunShape::changed);
    connect(m_custom_distribution,
            &SunShapeModel::changed,
            this,
            &SunShape::changed);

    // Initialization
    regenerate();
    update_current_distribution();
}

void SunShape::reset_current_distribution() {
    custom_distribution()->reset({
        { .angle = 0, .intensity = 1 },
        { .angle = 1, .intensity = 0.9 },
        { .angle = 2, .intensity = 0 },
    });
}

void SunShape::regenerate() {
    switch (m_shape) {
    case Shape::Gaussian: sample_gaussian(); break;
    case Shape::Pillbox: sample_pillbox(); break;
    case Shape::Buie_CSR: sample_buie(); break;
    case Shape::Custom: m_generated_distribution->clear(); break;
    case Shape::LimbDarkened: sample_limb_darkened(); break;
    }
    update_x_axis();
}

void SunShape::sample_gaussian() {
    // Point generation code referenced from app/src/sunshape
    // (SunShapeForm::UpdatePlot())

    int    num_points = 100;
    double theta_x    = 0;
    double theta_inc  = 3 * m_sigma / num_points;

    QVector<SunShapePoint> points;
    points.reserve(num_points);

    for (int i = 0; i < num_points; i++) {
        points.push_back({
            .angle     = theta_x,
            .intensity = 1.0 / exp(theta_x * theta_x / (2 * m_sigma * m_sigma)),
        });
        theta_x += theta_inc;
    }

    m_generated_distribution->reset(points);
}

void SunShape::sample_pillbox() {
    // Point generation code referenced from app/src/sunshape
    // (SunShapeForm::UpdatePlot())

    m_generated_distribution->reset({
        { .angle = 0, .intensity = 1 },
        { .angle = m_half_width, .intensity = 1 },
        { .angle = m_half_width, .intensity = 0 },
    });
}

void SunShape::sample_buie() {
    constexpr double min_csr_exclusive      = 0.0;
    constexpr double max_supported_csr      = 0.8;
    constexpr double solar_disk_radius_mrad = 4.65;
    constexpr double max_sample_angle_mrad  = 43.6;
    constexpr double sample_step_mrad       = 0.01;
    constexpr int    sample_count_estimate  = 4361;

    if (m_csr <= min_csr_exclusive || m_csr > max_supported_csr) {
        m_generated_distribution->reset();
        return;
    }

    // Buie's circumsolar model is driven by CSR, but the published intensity
    // equations use chi. These piecewise fits map the supported CSR range into
    // chi before deriving the aureole power-law terms below.
    const auto chi_from_csr = [](double csr) {
        if (csr > 0.145) {
            return -0.04419909985804843 +
                   csr * (1.401323894233574 +
                          csr * (-0.3639746714505299 +
                                 csr * (-0.9579768560161194 +
                                        1.1550475450828657 * csr)));
        }

        if (csr > 0.035) {
            return 0.022652077593662934 +
                   csr *
                       (0.5252380349996234 +
                        (2.5484334534423887 - 0.8763755326550412 * csr) * csr);
        }

        return 0.004733749294807862 +
               csr * (4.716738065192151 +
                      csr * (-463.506669149804 +
                             csr * (24745.88727411664 +
                                    csr * (-606122.7511711778 +
                                           5521693.445014727 * csr))));
    };

    const auto disk_intensity = [](double theta_mrad) {
        return std::cos(0.326 * theta_mrad) / std::cos(0.308 * theta_mrad);
    };

    const double chi   = chi_from_csr(m_csr);
    const double kappa = 0.9 * std::log(13.5 * chi) * std::pow(chi, -0.3);
    const double gamma = 2.2 * std::log(0.52 * chi) * std::pow(chi, 0.43) - 0.1;
    const double disk_edge_intensity = disk_intensity(solar_disk_radius_mrad);

    QVector<SunShapePoint> points;
    points.reserve(sample_count_estimate);

    // Store the non-negative radial profile. The graph mirrors it on demand
    // for display, while backend exports keep the radial profile unchanged.
    for (double theta = 0.0; theta <= max_sample_angle_mrad;
         theta += sample_step_mrad) {
        double intensity;

        if (theta <= solar_disk_radius_mrad) {
            intensity = disk_intensity(theta);
        } else {
            // Beyond the disk edge, Buie models the circumsolar aureole as a
            // power law. Cap it at the disk-edge value to avoid a
            // discontinuity.
            intensity = std::exp(kappa) * std::pow(theta, gamma);
            intensity = std::min(intensity, disk_edge_intensity);
        }

        points.push_back({
            .angle     = theta,
            .intensity = intensity,
        });
    }

    m_generated_distribution->reset(points);
}

void SunShape::sample_limb_darkened() {
    constexpr double disk_edge  = 4.65;
    constexpr int    num_points = 100;
    double           theta      = 0.0;
    double           theta_inc  = disk_edge / num_points;

    QVector<SunShapePoint> points;
    points.reserve(num_points + 1);

    for (int i = 0; i <= num_points; ++i) {
        points.push_back({
            .angle     = theta,
            .intensity = std::cos(0.326 * theta) / std::cos(0.308 * theta),
        });
        theta += theta_inc;
    }

    m_generated_distribution->reset(points);
}

void SunShape::update_x_axis() {
    QPointer<SunShapeModel> gdist = m_generated_distribution;
    QPointer<SunShapeModel> cdist = m_custom_distribution;

    switch (m_shape) {
    case Shape::Gaussian:
        gdist->set_x_axis_from(-3.3 * m_sigma);
        gdist->set_x_axis_to(3.3 * m_sigma);
        break;
    case Shape::Pillbox:
        gdist->set_x_axis_from(-3.3 * m_half_width);
        gdist->set_x_axis_to(3.3 * m_half_width);
        break;
    case Shape::Buie_CSR:
        gdist->set_x_axis_from(-20.0);
        gdist->set_x_axis_to(20.0);
        break;
    case Shape::LimbDarkened:
        gdist->set_x_axis_from(-1.3 * 4.65);
        gdist->set_x_axis_to(1.3 * 4.65);
        break;
    case Shape::Custom:
        // Code referenced from app/src/sunshape (SunShapeForm::UpdatePlot())
        if (cdist->count() >= 2) {
            double max_x = 0;
            for (int i = 0; i < cdist->count(); i++) {
                double angle = std::abs(cdist->get_at(i)->angle);
                if (angle > max_x) max_x = angle;
            }

            if (max_x == 0) {
                cdist->set_x_axis_from(-1.3);
                cdist->set_x_axis_to(1.3);
            } else {
                cdist->set_x_axis_from(-1.3 * max_x);
                cdist->set_x_axis_to(1.3 * max_x);
            }
        }
        break;
    }
}

void SunShape::update_current_distribution() {
    if (m_shape == Shape::Custom)
        set_current_distribution(m_custom_distribution);
    else
        set_current_distribution(m_generated_distribution);
}

SolarCalculatorData::SolarCalculatorData(QObject* parent) : QObject(parent) {
    connect(this,
            &SolarCalculatorData::calculator_changed,
            this,
            &SolarCalculatorData::changed);
    connect(this,
            &SolarCalculatorData::latitude_changed,
            this,
            &SolarCalculatorData::changed);
    connect(this,
            &SolarCalculatorData::longitude_changed,
            this,
            &SolarCalculatorData::changed);
    connect(this,
            &SolarCalculatorData::year_changed,
            this,
            &SolarCalculatorData::changed);
    connect(this,
            &SolarCalculatorData::month_changed,
            this,
            &SolarCalculatorData::changed);
    connect(this,
            &SolarCalculatorData::day_changed,
            this,
            &SolarCalculatorData::changed);
    connect(this,
            &SolarCalculatorData::hour_changed,
            this,
            &SolarCalculatorData::changed);
    connect(this,
            &SolarCalculatorData::minute_changed,
            this,
            &SolarCalculatorData::changed);
    connect(this,
            &SolarCalculatorData::second_changed,
            this,
            &SolarCalculatorData::changed);
    connect(this,
            &SolarCalculatorData::timezone_offset_changed,
            this,
            &SolarCalculatorData::changed);
    connect(this,
            &SolarCalculatorData::optional_solpos_fields_changed,
            this,
            &SolarCalculatorData::changed);
    connect(this,
            &SolarCalculatorData::interval_changed,
            this,
            &SolarCalculatorData::changed);
    connect(this,
            &SolarCalculatorData::optional_spa_fields_changed,
            this,
            &SolarCalculatorData::changed);
    connect(this,
            &SolarCalculatorData::dut1_changed,
            this,
            &SolarCalculatorData::changed);
    connect(this,
            &SolarCalculatorData::altitude_changed,
            this,
            &SolarCalculatorData::changed);
    connect(this,
            &SolarCalculatorData::pressure_changed,
            this,
            &SolarCalculatorData::changed);
    connect(this,
            &SolarCalculatorData::temperature_changed,
            this,
            &SolarCalculatorData::changed);
}

DateTime SolarCalculatorData::get_datetime_data() const {
    return DateTime { }; // TODO: stub
}

void SolarCalculatorData::set_spring() {
    set_month(3);
    set_day(20);
}

void SolarCalculatorData::set_summer() {
    set_month(6);
    set_day(21);
}

void SolarCalculatorData::set_fall() {
    set_month(9);
    set_day(22);
}

void SolarCalculatorData::set_winter() {
    set_month(12);
    set_day(21);
}

void SolarCalculatorData::set_dawn() {
    set_hour(6);
    set_minute(0);
    set_second(0);
}

void SolarCalculatorData::set_mid_morning() {
    set_hour(9);
    set_minute(0);
    set_second(0);
}

void SolarCalculatorData::set_noon() {
    set_hour(12);
    set_minute(0);
    set_second(0);
}

void SolarCalculatorData::set_mid_afternoon() {
    set_hour(15);
    set_minute(0);
    set_second(0);
}

void SolarCalculatorData::set_golden_hour() {
    set_hour(17);
    set_minute(0);
    set_second(0);
}

void SolarCalculatorData::set_dusk() {
    set_hour(19);
    set_minute(0);
    set_second(0);
}


SunShapeModel::SunShapeModel(QObject* parent) : StructTableModel(parent) {
    connect(this, &SunShapeModel::dataChanged, this, &SunShapeModel::changed);
    connect(this, &SunShapeModel::rowsInserted, this, &SunShapeModel::changed);
    connect(this, &SunShapeModel::rowsRemoved, this, &SunShapeModel::changed);
    connect(this, &SunShapeModel::modelReset, this, &SunShapeModel::changed);
}

std::vector<double> SunShapeModel::get_angle_data() {
    // Backend sun-shape data is radial-only. Views that need a symmetric curve
    // should mirror these non-negative points at presentation time.
    auto points = normalized_radial_points(m_records);

    std::vector<double> result;
    for (auto const& point : std::as_const(points)) {
        result.push_back(point.angle);
    }
    return result;
}

std::vector<double> SunShapeModel::get_intensity_data() {
    // Keep this paired with get_angle_data(): both export the same normalized,
    // non-negative radial profile to the backend.
    auto points = normalized_radial_points(m_records);

    std::vector<double> result;
    for (auto const& point : std::as_const(points)) {
        result.push_back(point.intensity);
    }
    return result;
}

QVariantList SunShapeModel::variant_data() {
    auto points = normalized_radial_points(m_records);

    QVariantList custom_shape;
    for (auto const& source : std::as_const(points)) {
        QVariantMap point;
        point["angle"]     = source.angle;
        point["intensity"] = source.intensity;
        custom_shape.append(point);
    }
    return custom_shape;
}

void SunShapeModel::set_variant_data(QVariantList data) {
    QVector<SunShapePoint> points;
    for (const auto& item : data) {
        QVariantMap point = item.toMap();
        points.push_back({
            .angle     = std::abs(point["angle"].toDouble()),
            .intensity = point["intensity"].toDouble(),
        });
    }
    reset(normalized_radial_points(points));
}

int SunShapeModel::count() const {
    return rowCount();
}

void SunShapeModel::append(double angle, double intensity) {
    StructTableModel::append({ std::abs(angle), intensity });
    emit countChanged();
}

void SunShapeModel::reset(QVector<SunShapePoint> points) {
    StructTableModel<SunShapePoint>::reset(points);
    emit countChanged();
}

void SunShapeModel::remove(int index) {
    if (index < 0 || index >= m_records.count()) return;
    remove_at(index);
    emit countChanged();
}

void SunShapeModel::clear() {
    reset();
}

void SunShapeModel::copy_to_clipboard() {
    QString text   = "Angle (mrad)\tIntensity\n";
    auto    points = normalized_radial_points(m_records);

    for (auto const& point : std::as_const(points)) {
        text += QString::number(point.angle) + "\t" +
                QString::number(point.intensity) + "\n";
    }
    QGuiApplication::clipboard()->setText(text);
}

void SunShapeModel::paste_from_clipboard() {
    QVariantList rows;
    QString      text = QGuiApplication::clipboard()->text();
    for (auto const& line : text.split('\n')) {
        if (line.trimmed() == "") continue;
        QStringList v = line.split(QRegularExpression("[\\t,]"));
        if (v.length() >= 2) {
            QVariantMap row;
            row["angle"]     = v[0].trimmed().toDouble();
            row["intensity"] = v[1].trimmed().toDouble();
            rows.append(row);
        }
    }
    set_variant_data(rows);
}

Data::SunShape SunShape::get_sunshape_data() const {
    switch (m_shape) {
    case Shape::Gaussian: return Data::SunShape::GAUSSIAN;
    case Shape::Pillbox: return Data::SunShape::PILLBOX;
    case Shape::Buie_CSR: return Data::SunShape::BUIE_CSR;
    case Shape::Custom: return Data::SunShape::USER_DEFINED;
    case Shape::LimbDarkened: return Data::SunShape::LIMBDARKENED;
    default: return Data::SunShape::UNKNOWN;
    }
}

SolarPositionData::SolarPositionData(QObject* parent) : QObject(parent) {
    connect(
        this, &SolarPositionData::x_changed, this, &SolarPositionData::changed);
    connect(
        this, &SolarPositionData::y_changed, this, &SolarPositionData::changed);
    connect(
        this, &SolarPositionData::z_changed, this, &SolarPositionData::changed);
    connect(this,
            &SolarPositionData::from_calculator_changed,
            this,
            &SolarPositionData::changed);
}


} // namespace SolTrace::GUI::App

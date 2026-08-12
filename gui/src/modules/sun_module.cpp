#include "sun_module.h"
#include <QDebug>
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

} // namespace SolTrace::GUI::App

#pragma once


#include "database/components.h"
#include "database/database.h"
#include "module_common.h"
#include "modules/sun/solar_calculator_data.h"
#include "modules/sun/solar_position_data.h"
#include "modules/sun/sun_shape.h"
#include "ray_source.hpp"
#include "solar_position_calculator.hpp"
#include "utilities/notification.h"
#include "utilities/qt_helpers.h"

#include <QDateTime>
#include <QMetaObject>
#include <QObject>
#include <QTimeZone>
#include <QVector>
#include <qqmlintegration.h>

namespace SolTrace::GUI::App {

class SunModule : public QObject {
    Q_OBJECT
    QML_ELEMENT

private:
    void update_database_connections();
    void load_from_database();

    // This should be const, but the library has non-const getters
    void load_from_ray_source(SD::RaySource&    ray_source,
                              db::RaySourceType source_type);

    void                                 write_shape_to_database();
    QString                              write_position_to_database();
    Data::SolarPositionCalculationMethod selected_calculation_method() const;

    Data::SolarPositionCalculator    m_calculator;
    QVector<QMetaObject::Connection> m_database_connections;
    bool                             m_loading_from_database        = false;
    bool                             m_writing_to_database          = false;
    bool                             m_updating_calculated_position = false;

public:
    explicit SunModule(QObject* parent = nullptr);

    QOBJECT_READONLY_PROPERTY(StatusComponent, status);
    QOBJECT_WRITABLE_PROPERTY(db::Database, current_database)

    QOBJECT_READONLY_PROPERTY(SunShape, shape)

    enum class Type { Directional, PointSource };
    enum class DirectionalPositionType { Calculator, Angle };

    Q_ENUM(Type)
    Q_ENUM(DirectionalPositionType)

    Q_WRITABLE_PROPERTY(Type, type, Type::Directional)
    Q_WRITABLE_PROPERTY(DirectionalPositionType,
                        ds_position_type,
                        DirectionalPositionType::Calculator)

    // Current sun position
    QOBJECT_WRITABLE_PROPERTY(SolarPositionData, position)

    // Point source position data
    QOBJECT_READONLY_PROPERTY(SolarPositionData, ps_position)

    // Directional sun position data (normalized direction vector)
    QOBJECT_READONLY_PROPERTY(SolarPositionData, ds_position)

    // Solar calculator fields (year, month, day, e.g.)
    QOBJECT_READONLY_PROPERTY(SolarCalculatorData, calc_data)

    Q_WRITABLE_PROPERTY(bool, sun_error_enabled, false)

public slots:
    QString apply_calculator(int    calculator,
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
                             double temperature);

    void    update_shape();
    void    update_type();
    QString update_position();

signals:
    void notify(ANotification);
};

} // namespace SolTrace::GUI::App

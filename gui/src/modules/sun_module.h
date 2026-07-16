#pragma once
#include "backend.h"
#include "module_common.h"
#include "utilities/notification.h"
#include "utilities/qt_helpers.h"
#include "utilities/structmodel.h"
#include <QDateTime>
#include <QMetaObject>
#include <QObject>
#include <QTimeZone>
#include <QVector>

namespace SolTrace::GUI::App {

struct SunShapePoint {
    double angle     = 0.0;
    double intensity = 0.0;

    RECORD_META(SunShapePoint, SM_EXPOSE_RW(angle), SM_EXPOSE_RW(intensity))
};

class SunShapeModel : public StructTableModel<SunShapePoint> {
    Q_OBJECT
public:
    explicit SunShapeModel(QObject* parent = nullptr);

    std::vector<double> get_angle_data();
    std::vector<double> get_intensity_data();

    Q_WRITABLE_PROPERTY(double, x_axis_from, -5)
    Q_WRITABLE_PROPERTY(double, x_axis_to, 5)
    Q_WRITABLE_PROPERTY(double, y_axis_from, 0)
    Q_WRITABLE_PROPERTY(double, y_axis_to, 1.2)

    QVariantList variant_data();
    void         set_variant_data(QVariantList data);

public slots:
    void reset(QVector<SunShapePoint> points = {});
    void append(double angle = 0.0, double intensity = 0.0);
    void remove(int index);
    void clear();

    void copy_to_clipboard();
    void paste_from_clipboard();

    int count() const;

signals:
    void countChanged();
    void changed();
};

class SunShape : public QObject {
    Q_OBJECT
    QML_ELEMENT

    QList<SunShapePoint> m_custom_shape;

    void regenerate();
    void sample_gaussian();
    void sample_pillbox();
    void sample_buie();
    void sample_limb_darkened();
    void update_x_axis();

    void update_current_distribution();

public:
    explicit SunShape(QObject* parent = nullptr);

    SolTrace::Data::SunShape get_sunshape_data() const;

    // Note that this is separate from SolTrace's SunShape enum to maintain
    // independence from backend modifications
    enum class Shape { Gaussian, Pillbox, Buie_CSR, Custom, LimbDarkened };
    Q_ENUM(Shape)
    Q_WRITABLE_PROPERTY(Shape, shape, Shape::Gaussian)

    // Gaussian
    Q_WRITABLE_PROPERTY(double, sigma, 4.65)

    // Pillbox
    Q_WRITABLE_PROPERTY(double, half_width, 4.65)

    // Buie
    Q_WRITABLE_PROPERTY(double, csr, 0.1)

    // Generated distribution (using sigma, half_width, csr)
    QOBJECT_READONLY_PROPERTY(SunShapeModel, generated_distribution)

    // Custom distribution (using user-defined points)
    QOBJECT_READONLY_PROPERTY(SunShapeModel, custom_distribution)

    // Current distribution
    QOBJECT_WRITABLE_PROPERTY(SunShapeModel, current_distribution)

    // Reset to default value
    void reset_current_distribution();

signals:
    void changed();
};

class SolarCalculatorData : public QObject {
    Q_OBJECT
    QML_ELEMENT
public:
    explicit SolarCalculatorData(QObject* parent = nullptr);

    DateTime get_datetime_data() const;

    // Calculator
    enum class Calculator { Legacy, Duffie, SOLPOS, SPA };
    Q_ENUM(Calculator)
    Q_WRITABLE_PROPERTY(Calculator, calculator, Calculator::Legacy)

    // Position
    Q_WRITABLE_PROPERTY(double, latitude, 35.04)
    Q_WRITABLE_PROPERTY(double, longitude, -105.10)

    // Date
    Q_WRITABLE_PROPERTY(int, year, 2026)
    Q_WRITABLE_PROPERTY(int, month, 3)
    Q_WRITABLE_PROPERTY(int, day, 20)

    // Time
    Q_WRITABLE_PROPERTY(int, hour, 12)
    Q_WRITABLE_PROPERTY(int, minute, 0)
    Q_WRITABLE_PROPERTY(int, second, 0)

    // Timezone offset in hours
    Q_WRITABLE_PROPERTY(int, timezone_offset, -7)

    // SOLPOS
    Q_WRITABLE_PROPERTY(bool, optional_solpos_fields, false)
    Q_WRITABLE_PROPERTY(int, interval, 1) ///< Averaging interval in seconds

    // SPA Optional fields
    Q_WRITABLE_PROPERTY(bool, optional_spa_fields, false)
    Q_WRITABLE_PROPERTY(double, dut1, 0.0)
    Q_WRITABLE_PROPERTY(double, altitude, 1000)
    Q_WRITABLE_PROPERTY(double, pressure, 1013.25)
    Q_WRITABLE_PROPERTY(double, temperature, 20.0)

signals:
    void changed();

public slots:
    void set_spring();
    void set_summer();
    void set_fall();
    void set_winter();

    void set_morning();
    void set_noon();
    void set_afternoon();
};

class SolarPositionData : public QObject {
    Q_OBJECT
public:
    explicit SolarPositionData(QObject* parent = nullptr);

    Q_WRITABLE_PROPERTY(double, x, 1000.0)
    Q_WRITABLE_PROPERTY(double, y, 1000.0)
    Q_WRITABLE_PROPERTY(double, z, 1000.0)

    Q_WRITABLE_PROPERTY(double, azimuth, 90)
    Q_WRITABLE_PROPERTY(double, elevation, 90)

    Q_WRITABLE_PROPERTY(bool, from_calculator, true)

signals:
    void changed();
};

class SunModule : public QObject {
    Q_OBJECT
    QML_ELEMENT

private:
    void update_database_connections();
    void load_from_database();

    // This should be const, but the library has non-const getters
    void load_from_ray_source(SD::RaySource& ray_source,
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

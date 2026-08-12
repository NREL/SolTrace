#pragma once

#include "datetime.hpp"
#include "utilities/qt_helpers.h"

#include <QObject>
#include <qqmlintegration.h>

namespace SolTrace::GUI::App {

/// Input fields for solar-position calculators.
class SolarCalculatorData : public QObject {
    Q_OBJECT
    QML_ELEMENT
public:
    explicit SolarCalculatorData(QObject* parent = nullptr);

    /// Convert UI fields to the library DateTime representation.
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
    /// Set date fields to spring equinox defaults.
    void set_spring();

    /// Set date fields to summer solstice defaults.
    void set_summer();

    /// Set date fields to fall equinox defaults.
    void set_fall();

    /// Set date fields to winter solstice defaults.
    void set_winter();

    /// Set time fields to a dawn preset.
    void set_dawn();

    /// Set time fields to a mid-morning preset.
    void set_mid_morning();

    /// Set time fields to a noon preset.
    void set_noon();

    /// Set time fields to a mid-afternoon preset.
    void set_mid_afternoon();

    /// Set time fields to a low-sun preset.
    void set_golden_hour();

    /// Set time fields to a dusk preset.
    void set_dusk();
};

} // namespace SolTrace::GUI::App

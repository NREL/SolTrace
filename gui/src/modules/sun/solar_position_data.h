#pragma once

#include "utilities/qt_helpers.h"

#include <QObject>

namespace SolTrace::GUI::App {

/// Shared position/direction fields for point and directional sun editing.
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

} // namespace SolTrace::GUI::App

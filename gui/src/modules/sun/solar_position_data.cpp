#include "modules/sun/solar_position_data.h"

namespace SolTrace::GUI::App {

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

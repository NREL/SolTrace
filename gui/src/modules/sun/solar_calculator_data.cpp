#include "modules/sun/solar_calculator_data.h"

namespace SolTrace::GUI::App {

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


} // namespace SolTrace::GUI::App

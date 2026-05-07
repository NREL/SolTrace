

#include <solar_position_calculator.hpp>
#include <basic_sun_position.hpp>
#include <solpos00.h>
#include <lib_irradproc.h>
#include <constants.hpp>

#include <cmath>
#include <stdexcept>



namespace SolTrace::Data {

SolarPositionCalculator::SolarPositionCalculator()
    : method(SolarPositionCalculationMethod::LEGACY),
        year(0), month(0), day(0),
        hour(0), minute(0), second(0),
        latitude(0.0), longitude(0.0), timeZone(0.0),
        dut1(0.0), altitude(0.0), pressure(1013.25), temperature(20.0),
        location_set(false), date_set(false), time_set(false), calculated(false),
        Azimuth(0.0), Zenith(0.0), Elevation(0.0),
        X(0.0), Y(0.0), Z(0.0)
{
}

void SolarPositionCalculator::set_method(SolarPositionCalculationMethod method) {
    this->calculated = false;
    this->method = method;
}

void SolarPositionCalculator::set_location(double latitude, double longitude, double timeZone, double altitude) {
    this->calculated = false;

    if (latitude > 90 || latitude < -90) {
        throw std::invalid_argument("Latitude must be between -90 and 90 degrees.");
    }
    this->latitude = latitude;

    if (longitude > 180 || longitude < -180) {
        throw std::invalid_argument("Longitude must be between -180 and 180 degrees.");
    }
    this->longitude = longitude;

    if (timeZone < -12 || timeZone > 14) {
        throw std::invalid_argument("Time zone must be between -12 and +14 hours.");
    }
    this->timeZone = timeZone;

    if (altitude < -440 || altitude > 8850) {
        throw std::invalid_argument("Altitude must be between -440 and 8850 meters.");
    }
    this->altitude = altitude;

    this->location_set = true;
}

void SolarPositionCalculator::set_datetime(int year, int month, int day, int hour, int minute, int second) {
    this->calculated = false;
    set_date(year, month, day);
    set_time(hour, minute, second);
}

void SolarPositionCalculator::set_date(int year, int month, int day) {
    this->calculated = false;

    // TODO: are there any limits on year?
    this->year = year;

    if (month < 1 || month > 12) {
        throw std::invalid_argument("Month must be between 1 and 12.");
    }
    this->month = month;
    if (day < 1 || day > 31) {
        throw std::invalid_argument("Day must be between 1 and 31.");
    }
    this->day = day;
        
    this->date_set = true;
}

void SolarPositionCalculator::set_time(int hour, int minute, int second) {
    this->calculated = false;

    if (hour < 0 || hour > 23) {
        throw std::invalid_argument("Hour must be between 0 and 23.");
    }
    this->hour = hour;
    if (minute < 0 || minute > 59) {
        throw std::invalid_argument("Minute must be between 0 and 59.");
    }
    this->minute = minute;
    if (second < 0 || second > 59) {
        throw std::invalid_argument("Second must be between 0 and 59.");
    }
    this->second = second;

    this->time_set = true;
}

void SolarPositionCalculator::set_environment(double pressure, double temperature) {
    this->calculated = false;
    if (pressure < 950 || pressure > 1050) {
        throw std::invalid_argument("Pressure must be between 950 and 1050 millibars.");
    }
    this->pressure = pressure;

    if (temperature < -50 || temperature > 60) {
        throw std::invalid_argument("Temperature must be between -50 and 60 degrees Celsius.");
    }
    this->temperature = temperature;
}

void SolarPositionCalculator::calculate_sun_position() {

    if (!this->location_set) {
        throw std::runtime_error("Location not set before calculating sun position.");
    }
    if (!this->date_set) {
        throw std::runtime_error("Date not set before calculating sun position.");
    }
    if (!this->time_set) {
        throw std::runtime_error("Time not set before calculating sun position.");
    }

    double azimuth = -1.0, zenith = -1.0;
    switch (this->method) {
        case SolarPositionCalculationMethod::LEGACY:
        {
            int doy = day_of_year(this->month, this->day);
            double hour = static_cast<double>(this->hour) + static_cast<double>(this->minute) / 60.0 + static_cast<double>(this->second) / 3600.0;
            legacy_sun_position(this->latitude, static_cast<double>(doy), hour, &azimuth, &zenith);
            break;
        }
        case SolarPositionCalculationMethod::DUFFIE:
        {
            int doy = day_of_year(this->month, this->day);
            double hour = static_cast<double>(this->hour) + static_cast<double>(this->minute) / 60.0 + static_cast<double>(this->second) / 3600.0;
            duffie_sun_position(this->latitude, this->longitude, this->timeZone, static_cast<double>(doy), hour, &azimuth, &zenith);
            break;
        }
        case SolarPositionCalculationMethod::SOLPOS:
        {
            // Call SOLPOS sun position calculation

            //Instantiate the solpos object
            struct posdata SP, * pdat;
            pdat = &SP;			//point to structure for convenience
            S_init(pdat);		//Initialize the values
            pdat->latitude = float(this->latitude);		//[deg] {float} North is positive
            pdat->longitude = float(this->longitude);	//[deg] {float} Degrees east. West is negative
            pdat->timezone = float(this->timeZone);		//[hr] {float} Time zone, east pos. west negative. Mountain -7, Central -6, etc..
            pdat->year = this->year;			    //[year] {int} 4-digit year
            pdat->month = this->month;				//[mo] {int} (1-12)
            pdat->day = this->day;			        //[day] {int} Day of the month
            //pdat->daynum = 0;				        //[day] {int} Day of the year
            pdat->function &= ~S_DOY;               // Use month and day input
            pdat->hour = this->hour;			    //[hr] {int} 0-23
            pdat->minute = this->minute;		    //[min] {int} 0-59
            pdat->second = this->second;			//[sec]	{int} 0-59
            pdat->interval = 0;						//[sec] {int} Measurement interval. See solpos documentation.

            long retcode = 0;			//Initialize with no errors
            retcode = S_solpos(pdat);	//Call the solar position algorithm
            S_decode(retcode, pdat);	//Check the return code

            azimuth = SP.azim;
            zenith = SP.zenetr;
            break;
        }
        case SolarPositionCalculationMethod::SPA_ORIGINAL:
        {
            double sun_results[9]; //vector to hold the results of solarpos function

            solarpos(this->year,
                this->month,
                this->day,
                this->hour,
                this->minute,
                this->latitude,
                this->longitude,
                this->timeZone,
                sun_results);

            azimuth = sun_results[0] * R2D; // Convert radians to degrees
            zenith = sun_results[1] * R2D;  // Convert radians to degrees
            break;
        }
        case SolarPositionCalculationMethod::SPA:
        {
            double sun_results[9]; //vector to hold the results of solarpos function
            solarpos_spa(this->year,
                this->month,
                this->day,
                this->hour,
                this->minute,
                this->second,
                this->latitude,
                this->longitude,
                this->timeZone,
                this->dut1,
                this->altitude,
                this->pressure,
                this->temperature,
                15, 180, sun_results);

            azimuth = sun_results[0] * R2D; // Convert radians to degrees
            zenith = sun_results[1] * R2D;  // Convert radians to degrees
            break;
        }
        default:
            // Handle unknown method
            break;
    }

    this->Azimuth = azimuth;
    this->Zenith = zenith;
    this->Elevation = 90.0 - zenith;

    this->X = std::sin(zenith * D2R) * std::sin(azimuth * D2R);
    this->Y = std::sin(zenith * D2R) * std::cos(azimuth * D2R);
    this->Z = std::cos(zenith * D2R);
    this->calculated = true;
}

void SolarPositionCalculator::get_azimuth_zenith(double* azimuth, double* zenith) {
    if (!this->calculated) {
        this->calculate_sun_position();
    }
    *azimuth = this->Azimuth;
    *zenith = this->Zenith;
}

void SolarPositionCalculator::get_azimuth_elevation(double* azimuth, double* elevation) {
    if (!this->calculated) {
        this->calculate_sun_position();
    }
    *azimuth = this->Azimuth;
    *elevation = this->Elevation;
}

void SolarPositionCalculator::get_sun_vector(double* sun_x, double* sun_y, double* sun_z) {
    if (!this->calculated) {
        this->calculate_sun_position();
    }
    *sun_x = this->X;
    *sun_y = this->Y;
    *sun_z = this->Z;
}

}



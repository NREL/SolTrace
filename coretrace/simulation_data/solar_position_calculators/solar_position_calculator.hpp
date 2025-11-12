
#ifndef SOLTRACE_SUNPOSITION_H
#define SOLTRACE_SUNPOSITION_H

namespace SolTrace::Data {

enum SolarPositionCalculationMethod {
    LEGACY,
    DUFFIE,
    SOLPOS,
    SPA_ORIGINAL,
    SPA
};

class SolarPositionCalculator {

    SolarPositionCalculationMethod method;

    // Time inputs
    int year;
    int month;
    int day;
    int hour;
    int minute;
    int second;

    // Location inputs
    double latitude;    // Degrees, north positive
    double longitude;   // Degrees, east positive
    double timeZone;    // Hours from UTC, west longitudes negative

    // Optional inputs
    double dut1;        // Seconds, fractional second difference between UTC and UT which is used to adjust UTC for earth's irregular rotation rate (http://maia.usno.navy.mil/ser7/ser7.dat) (-1 to 1 second), internally managed (NOT user-configuarble)
    double altitude;    // Meters, above sea level
    double pressure;    // Millibars
    double temperature; // Degrees Celsius, dry-bulb

    // Boolean checks
    bool location_set;  // Flag to indicate if location has been set
    bool date_set;      // Flag to indicate if date has been set
    bool time_set;      // Flag to indicate if time has been set
    bool calculated;    // Flag to indicate if calculation has been performed with current inputs



    // Outputs
    double Azimuth;     // Degrees
    double Zenith;      // Degrees
    double Elevation;   // Degrees

    double X;           // Sun vector X component
    double Y;           // Sun vector Y component
    double Z;           // Sun vector Z component

    void calculate_sun_position();


public:
    SolarPositionCalculator();
    //~SolarPositionCalculator();

    // Setters
    void set_method(SolarPositionCalculationMethod method);
    void set_location(double latitude, double longitude, double timeZone, double altitude = 0.0);
    void set_datetime(int year, int month, int day, int hour, int minute, int second = 0);
    void set_date(int year, int month, int day);
    void set_time(int hour, int minute, int second = 0);

    void set_environment(double pressure, double temperature);

    // Getters
    void get_azimuth_zenith(double* azimuth, double* zenith);
    void get_azimuth_elevation(double* azimuth, double* elevation);
    void get_sun_vector(double* sun_x, double* sun_y, double* sun_z);
    // TODO: overload sun_vector function to return an array type
};

} // namespace SolTrace::Data

#endif
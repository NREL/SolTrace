
#include "basic_sun_position.hpp"

namespace SolTrace::Data {

int st_sun_position(double lat, double day, double hour,
                    double *x, double *y, double *z)
{
    /*
    computes the sun vector xyz given arguments
    lat : [deg] latitude
    day : [] day of the year
    hour : [hour] solar time. 12.00 corresponds to sun at maximum elevation and does not necessarily match local time

    xyz coordinate system:
        x: +west
        y: +zenith
        z: +north
    */
    
    double Declination, HourAngle, Elevation, Azimuth;

    Declination = R2D * asin(0.39795 * cos(0.98563 * D2R * (day - 173)));
    HourAngle = 15 * (hour - 12);
    Elevation = R2D * asin(sin(Declination * D2R) * sin(lat * D2R) + cos(Declination * D2R) * cos(HourAngle * D2R) * cos(lat * D2R));
    Azimuth = R2D * acos((sin(D2R * Declination) * cos(D2R * lat) - cos(D2R * Declination) * sin(D2R * lat) * cos(D2R * HourAngle)) / cos(D2R * Elevation) + 0.0000000001);
    if (sin(HourAngle * D2R) > 0.0)
        Azimuth = 360 - Azimuth;

    // TODO: Update coordinate system
    *x = -sin(Azimuth * D2R) * cos(Elevation * D2R);
    *y = sin(Elevation * D2R);
    *z = cos(Azimuth * D2R) * cos(Elevation * D2R);

    return 1;
}

} // namespace SolTrace::Data
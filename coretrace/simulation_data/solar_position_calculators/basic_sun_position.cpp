
#include "basic_sun_position.hpp"
#include "constants.hpp"

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
    double Elevation, Azimuth, Zenith;

    legacy_sun_position(lat, day, hour, &Azimuth, &Zenith);
    Elevation = 90.0 - Zenith;

    // TODO: Update coordinate system
    *x = -sin(Azimuth * D2R) * cos(Elevation * D2R);
    *y = sin(Elevation * D2R);
    *z = cos(Azimuth * D2R) * cos(Elevation * D2R);

    return 1;
}

void legacy_sun_position(double lat, double day, double hour, double* azimuth, double* zenith)
{
    /*
    Computes the sun Azimuth and Zenith angles using the legacy SolTrace method. This method is similar to that presented in
    Duffie, J.A. and Beckman, W.A., "Solar Engineering of Thermal Processes," 4th Edition, Wiley, 2013. but not identical.
        lat : [deg] latitude
        day : [] day of the year
        hour : [hour] solar time. 12.00 corresponds to sun at maximum elevation and does not necessarily match local time

    Returns:
        azimuth : [deg] azimuth angle, measured east from north
        zenith : [deg] zenith angle

    */

    double Declination, HourAngle, Elevation, Azimuth;

    Declination = R2D * asin(0.39795 * cos(0.98563 * D2R * (day - 173)));
    HourAngle = 15 * (hour - 12);
    Elevation = R2D * asin(sin(Declination * D2R) * sin(lat * D2R) + cos(Declination * D2R) * cos(HourAngle * D2R) * cos(lat * D2R));
    Azimuth = R2D * acos((sin(D2R * Declination) * cos(D2R * lat) - cos(D2R * Declination) * sin(D2R * lat) * cos(D2R * HourAngle)) / cos(D2R * Elevation) + 0.0000000001);
    if (sin(HourAngle * D2R) > 0.0)
        Azimuth = 360 - Azimuth;

    *azimuth = Azimuth;
    *zenith = 90.0 - Elevation;
}

void duffie_sun_position(double lat, double lng, double tz, double day, double hour, double* azimuth, double* zenith)
{
    /*
    Computes the sun Azimuth and Zenith angles using the method presented in
    Duffie, J.A. and Beckman, W.A., "Solar Engineering of Thermal Processes," 4th Edition, Wiley, 2013.

    lat : [deg] latitude
    day : [] day of the year
    hour : [hour] solar time. 12.00 corresponds to sun at maximum elevation and does not necessarily match local time

    Returns:
    azimuth : [deg] azimuth angle, measured east from north
    zenith : [deg] zenith angle
    */
    // TODO: we could improve this by 
    //  - calculate a fractional day of year
    //  - Use More accurate equation for declination


    double B = (day - 1) * 360. / 365.; // [degrees] Equation 1.4.2
    double E = 229.2 * (0.000075 + 0.001868 * cos(D2R * B) - 0.032077 * sin(D2R * B) - 0.014615 * cos(D2R * 2 * B) - 0.04089 * sin(D2R * 2 * B)); // [minutes] Equation 1.5.3 
    double solar_time = hour + (lng / 15.0 - tz) + E / 60.0;  // [hour] solar time, Equation 1.5.2

    double HourAngle = 15 * (solar_time - 12);
    //double Declination = 23.45 * sin(D2R * (360.0 * (284 + day) / 365.0));
    double Declination = R2D * (0.006918 - 0.399912 * cos(D2R * B) + 0.070257 * sin(D2R * B) - 0.006758 * cos(D2R * 2 * B) + 0.000907 * sin(D2R * 2 * B) - 0.002697 * cos(D2R * 3 * B) + 0.00148 * sin(D2R * 3 * B));
    double Zenith = R2D * acos(cos(D2R * lat) * cos(D2R * Declination) * cos(D2R * HourAngle) + sin(D2R * lat) * sin(D2R * Declination));
    double ratio = (cos(D2R * Zenith) * sin(D2R * lat) - sin(D2R * Declination)) / (sin(D2R * Zenith) * cos(D2R * lat));

    double Azimuth;
    // adjusting numerical error
    if (ratio > 1.0 && (ratio - 1.0) < 1.e-6)
        Azimuth = 0.0;
    else if (ratio < -1.0 && (-1.0 - ratio) < 1.e-6)
        Azimuth = 180.0;
    else 
        Azimuth = R2D * fabs(acos(ratio));

    if (HourAngle < 0)
        Azimuth *= -1.0;
    Azimuth += 180.0;    // Shift to typical North reference (east positive angle)

    *azimuth = Azimuth;
    *zenith = Zenith;
}



} // namespace SolTrace::Data
#ifndef BASIC_SUN_H
#define BASIC_SUN_H

namespace SolTrace::Data {

int st_sun_position(double lat, double day, double hour,double *x, double *y, double *z);

void legacy_sun_position(double lat, double day, double hour, double *azimuth, double *zenith);

void duffie_sun_position(double lat, double lng, double tz, double day, double hour, double* azimuth, double* zenith);

} // namespace SolTrace::Data

#endif

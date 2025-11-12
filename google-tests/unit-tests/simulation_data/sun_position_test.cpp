#include <gtest/gtest.h>

#include <basic_sun_position.hpp>
#include <solpos00.h>
#include <lib_irradproc.h>
#include <solar_position_calculator.hpp>

#include <vector>
#include <limits>
#include <cmath>
//#include <iostream>
#include "common.hpp"

TEST(SolarPositionCalculator, LegacyDirectTest)
{
    double latitude = 35.0;
    double day_of_year = 172.0; // June 21
    double hour = 12.0; // Solar noon

    double x, y, z;

    int res = SolTrace::Data::st_sun_position(latitude, day_of_year, hour, &x, &y, &z);

    EXPECT_EQ(res, 1);

    double magnitude = std::sqrt(x * x + y * y + z * z);
    EXPECT_NEAR(magnitude, 1.0, 1e-6);

    // Expected values for June 21 at solar noon at 35 degrees latitude
    double expected_x = 0.0;
    double expected_y = 0.9797377;
    double expected_z = -0.20028;
    EXPECT_NEAR(x, expected_x, 1e-4);
    EXPECT_NEAR(y, expected_y, 1e-4);
    EXPECT_NEAR(z, expected_z, 1e-4);
}

TEST(SolarPositionCalculator, SolPosDirectTest)
{
	double hour = 12.0;
    double latitude = 35.0;
    double longitude = -106.0;
    double timezone = -7.0;
    int month = 6;
	int day_of_month = 21;
    int day_of_year = 172;

	//Instantiate the solpos object
	struct posdata SP, * pdat;
	pdat = &SP;			//point to structure for convenience
	S_init(pdat);		//Initialize the values

	//Calculate minutes/seconds
	double
		mins = 60. * (hour - floor(hour)),
		secs = 60. * (mins - floor(mins));

	pdat->latitude = float(latitude);		//[deg] {float} North is positive
	pdat->longitude = float(longitude);		//[deg] {float} Degrees east. West is negative
	pdat->timezone = float(timezone);		//[hr] {float} Time zone, east pos. west negative. Mountain -7, Central -6, etc..
	pdat->year = 2011;						//[year] {int} 4-digit year
	pdat->month = int(month);				//[mo] {int} (1-12)
	pdat->day = int(day_of_month);			//[day] {int} Day of the month
	pdat->daynum = day_of_year;				//[day] {int} Day of the year
	pdat->hour = int(hour + .0001);			//[hr] {int} 0-23
	pdat->minute = int(mins);				//[min] {int} 0-59
	pdat->second = int(secs);				//[sec]	{int} 0-59
	pdat->interval = 0;						//[sec] {int} Measurement interval. See solpos documentation.

	long retcode = 0;			//Initialize with no errors
	retcode = S_solpos(pdat);	//Call the solar position algorithm
	S_decode(retcode, pdat);	//Check the return code

	double azimuth = SP.azim;
	double zenith = SP.zenetr;

	double expected_azimuth = 173.417953;
    double expected_zenith = 11.630642;
    EXPECT_NEAR(azimuth, expected_azimuth, 1e-4);
    EXPECT_NEAR(zenith, expected_zenith, 1e-4);
}

TEST(SolarPositionCalculator, sunriseAndSunsetAtDifferentLocationsTest_lib_irradproc)
{
    /*locations to test:
    western hemisphere: Golden CO
    eastern hemisphere: Berlin Germany
    southern hemisphere: Lima Peru
    location near Greenwich meridian with negative longitude and positive time zone: Madrid Spain
    location near the international dateline with positive longitude and negative time zone: Lomaji, Fiji
    arctic circle: Kotzebue, Alaska
    arctic circle #2: Point Hope, Alaska
    arctic circle #3: Kotzebue, Alaska on the first day of continuous days
    */
    double e = 0.001;
    std::vector<double> latitudes = { 39.77, 52.5, -12.03, 40.43, -17.75, 66.9, 68.35, 66.9 };
    std::vector<double> longitudes = { -105.22, 13.3, -77.06, -3.72, -179.3, -162.6, -166.8, -162.6 };
    std::vector<double> time_zones = { -7, 1, -5, 1, 12, -9, -9, -9 };
    std::vector<double> sunrise_times = { 4.636, 3.849, 6.521, 5.833, 6.513, -100.0, 2.552, -100.0 };
    std::vector<double> sunset_times = { 19.455, 20.436, 17.814, 20.723, 17.449, 100.0, 25.885, 100.0 };
    std::vector<int> month = { 6, 6, 6, 6, 6, 6, 7, 6 };
    std::vector<int> day = { 21, 21, 21, 21, 21, 21, 14, 11 };


    double sun_results[9]; //vector to hold the results of solarpos function
    for (size_t i = 0; i < latitudes.size(); i++)
    {
        //run the solarpos function and check sunrise and sunset for each location
        solarpos(2010, month[i], day[i], 14, 30, latitudes[i], longitudes[i], time_zones[i], sun_results);
        EXPECT_NEAR((double)sun_results[4], sunrise_times[i], e) << "sunrise time for lat " << latitudes[i] << " long " << longitudes[i] << " failed\n";
        EXPECT_NEAR((double)sun_results[5], sunset_times[i], e) << "sunset time for lat " << latitudes[i] << " long " << longitudes[i] << " failed\n";
    }
}

TEST(SolarPositionCalculator, solarposTest_lib_irradproc) {

    double lat = 31.6340;
    double lon = 74.8723;
    double tz = 5.5;
    double year = 2017;
    double month = 7;
    double day = 19;

    double sun[9];
    std::vector<double> sunrise_times;
    std::vector<double> sunset_times;
    double e = 0.0001;
    /* Just before sunrise test case */
    solarpos(year, month, day, 4, 30, lat, lon, tz, sun);
    std::vector<double> solution = { 0.95662, 1.79457, -0.223771, 0.363938, 5.70882, 19.5183, 0.968276, 3.88646, 0 };
    sunrise_times.push_back(solution[4]);
    sunset_times.push_back(solution[5]);
    for (int i = 0; i < 9; i++) {
        EXPECT_NEAR((double)sun[i], solution[i], e) << "hourly before-sunrise case, parameter " << i << " fail\n";
    }
    solarpos(year, month, day, 5, 15, lat, lon, tz, sun);
    solution = { 1.0744, 1.65255, -0.0817513, 0.363839, 5.7091, 19.518, 0.96828, 4.63642, 0 };
    sunrise_times.push_back(solution[4]);
    sunset_times.push_back(solution[5]);
    for (int i = 0; i < 9; i++) {
        EXPECT_NEAR((double)sun[i], solution[i], e) << "15m before-sunrise case, parameter " << i << " fail\n";
    }

    /* Just after sunset test case */
    solarpos(year, month, day, 20, 30, lat, lon, tz, sun);
    solution = { 5.28748, 1.75391, -0.183117, 0.361807, 5.71544, 19.5131, 0.968361, 19.8857, 0 };
    sunrise_times.push_back(solution[4]);
    sunset_times.push_back(solution[5]);
    for (int i = 0; i < 9; i++) {
        EXPECT_NEAR((double)sun[i], solution[i], e) << "hourly after-sunset case, parameter " << i << " fail\n";
    }
    solarpos(year, month, day, 19, 45, lat, lon, tz, sun);
    solution = { 5.17431, 1.60864, -0.0378397, 0.361908, 5.71513, 19.5133, 0.968357, 19.1358, 0 };
    sunrise_times.push_back(solution[4]);
    sunset_times.push_back(solution[5]);
    for (int i = 0; i < 9; i++) {
        EXPECT_NEAR((double)sun[i], solution[i], e) << "15m after-sunrise case, parameter " << i << " fail\n";
    }
}

TEST(SolarPositionCalculator, solarposTest_lib_irradproc_other) {
    double lat = 31.6340;
    double lon = 74.8723;
    double tz = 5.5;
    double year = 2017;
    double month = 7;
    double day = 19;

    double sun[9];
    std::vector<double> sunrise_times;
    std::vector<double> sunset_times;
    double e = 0.0001;

    solarpos(year, month, day, 5, 30, lat, lon, tz, sun);
    std::vector<double> solution = { 1.11047, 1.6031, -0.0323028, 0.363806, 5.70924, 19.5179, 0.968281, 4.88641, 0 };
    sunrise_times.push_back(solution[4]);
    sunset_times.push_back(solution[5]);
    for (int i = 0; i < 9; i++) {
        EXPECT_NEAR((double)sun[i], solution[i], e) << "sunrise case, parameter " << i << " fail\n";
    }
}

TEST(SolarPositionCalculator, sunriseAndSunsetAtDifferentLocationsTest_spa_lib_irradproc) {
    /*locations to test:
    western hemisphere: Golden CO
    eastern hemisphere: Berlin Germany
    southern hemisphere: Lima Peru
    location near Greenwich meridian with negative longitude and positive time zone: Madrid Spain
    location near the international dateline with positive longitude and negative time zone: Lomaji, Fiji
    arctic circle: Kotzebue, Alaska
    arctic circle #2: Point Hope, Alaska
    arctic circle #3: Kotzebue, Alaska on the first day of continuous days
    */
    double e = 0.001;
    std::vector<double> latitudes = { 39.77, 52.5, -12.03, 40.43, -17.75, 66.9, 68.35, 66.9 };
    std::vector<double> longitudes = { -105.22, 13.3, -77.06, -3.72, -179.3, -162.6, -166.8, -162.6 };
    std::vector<double> time_zones = { -7, 1, -5, 1, 12, -9, -9, -9 };
    std::vector<double> sunrise_times = { 4.549, 3.726, 6.458, 5.745, 6.451, -100.0, 2.7847330, -100.0 };
    std::vector<double> sunset_times = { 19.541, 20.559, 17.877, 20.810, 17.514, 100.0, 25.4795, 100.0 };
    std::vector<int> month = { 6, 6, 6, 6, 6, 6, 7, 6 };
    std::vector<int> day = { 21, 21, 21, 21, 21, 21, 20, 11 };
    std::vector<int> alt = { 1730, 34, 154, 667, 0, 6, 2, 6 };

    double sun_results[9]; //vector to hold the results of solarpos function
    for (size_t i = 0; i < latitudes.size(); i++)
    {
        //run the solarpos function and check sunrise and sunset for each location
        solarpos_spa(2010, month[i], day[i], 14, 30, 0, latitudes[i], longitudes[i], time_zones[i], 0, alt[i], 0, 1016, 15, 180, sun_results);
        EXPECT_NEAR((double)sun_results[4], sunrise_times[i], e) << "sunrise time for lat " << latitudes[i] << " long " << longitudes[i] << " failed\n";
        EXPECT_NEAR((double)sun_results[5], sunset_times[i], e) << "sunset time for lat " << latitudes[i] << " long " << longitudes[i] << " failed\n";
    }
}

TEST(SolarPositionCalculator, atmos_refractionTest_spa_lib_irradproc) {
    //Test to check for atmospheric refraction correction occurring only if sun is above horizon
    double latitude = 31.6430;
    double longitude = 74.8723;
    double time_zone = 5.5;
    double elevation_angle = -.00175; //topocentric elevation angle corrected for atmospheric refraction (radians)
    double e = 0.001;
    //double sunset_time = 17.514;
    int month = 7;
    int day = 19;
    double sun_results[9];
    double alt = 0;
    solarpos_spa(2017, month, day, 5, 39, 0, latitude, longitude, time_zone, 0, 234, 1013.25, 15, latitude, 180, sun_results);
    EXPECT_NEAR((double)sun_results[2], elevation_angle, e) << "elevation angle for lat " << latitude << " long " << longitude << " failed\n";
}

TEST(SolarPositionCalculator, InputsValidation) {

    double azimuth, zenith;
    SolarPositionCalculator solar_position;

    // Test case: Attempt to get azimuth and zenith without setting location or datetime
    EXPECT_THROW(solar_position.get_azimuth_zenith(&azimuth, &zenith), std::runtime_error);

    // Test case: invalid latitude
    double latitude = 100.0; // Invalid latitude
    EXPECT_THROW(solar_position.set_location(latitude, -106.0, -7.0), std::invalid_argument);

    // Test case: invalid longitude
    latitude = 40.0;
    double longitude = 200.0; // Invalid longitude
    EXPECT_THROW(solar_position.set_location(latitude, longitude, -7.0), std::invalid_argument);

    // Test case: invalid timezone
    longitude = -105.0;
    double timezone = 15.0; // Invalid timezone
    EXPECT_THROW(solar_position.set_location(latitude, longitude, timezone), std::invalid_argument);

    // TEST case: Attempt to get azimuth and zenith without setting datetime
    timezone = -7.0;
    EXPECT_NO_THROW(solar_position.set_location(latitude, longitude, timezone));
    EXPECT_THROW(solar_position.get_azimuth_zenith(&azimuth, &zenith), std::runtime_error);

    // Test case: invalid month
    int month = 13; // Invalid month
    EXPECT_THROW(solar_position.set_datetime(2024, month, 15, 12, 0, 0), std::invalid_argument);

    // Test case: invalid day
    month = 6;
    int day = 32; // Invalid day
    EXPECT_THROW(solar_position.set_datetime(2024, month, day, 12, 0, 0), std::invalid_argument);

    // TEST case: Attempt to get azimuth and zenith without setting time
    day = 21;
    EXPECT_NO_THROW(solar_position.set_date(2024, month, day));
    EXPECT_THROW(solar_position.get_azimuth_zenith(&azimuth, &zenith), std::runtime_error);

    // Test case: invalid hour
    int hour = 25; // Invalid hour
    EXPECT_THROW(solar_position.set_datetime(2024, month, day, hour, 0, 0), std::invalid_argument);

    // Test case: invalid minute
    hour = 12;
    int minute = 60; // Invalid minute
    EXPECT_THROW(solar_position.set_datetime(2024, month, day, hour, minute, 0), std::invalid_argument);

    // Test case: invalid second
    minute = 30;
    int second = 60; // Invalid second
    EXPECT_THROW(solar_position.set_datetime(2024, month, day, hour, minute, second), std::invalid_argument);

    // Now set valid datetime
    second = 0;
    EXPECT_NO_THROW(solar_position.set_datetime(2024, month, day, hour, minute, second));

    // Now getting azimuth and zenith should not throw
    EXPECT_NO_THROW(solar_position.get_azimuth_zenith(&azimuth, &zenith));

    // Additional test case: invalid altitude
    double altitude = 9000.0; // Invalid altitude
    EXPECT_THROW(solar_position.set_location(latitude, longitude, timezone, altitude), std::invalid_argument);

    // Additional test case: invalid pressure
    double pressure = 900.0; // Invalid pressure
    EXPECT_THROW(solar_position.set_environment(pressure, 20.0), std::invalid_argument);

    // Additional test case: invalid temperature
    pressure = 1013.25;
    double temperature = -60.0; // Invalid temperature
    EXPECT_THROW(solar_position.set_environment(pressure, temperature), std::invalid_argument);
}

TEST(SolarPositionCalculator, CrossValidationTest) {
    // Example test case for SolarPositionCalculator
    double latitude = 40.0;
    double longitude = -105.0;
    double timezone = -7.0;
    int year = 2025;
    int month = 6;
    int day = 20;
    int hour = 12;

    SolarPositionCalculator solar_position;
    solar_position.set_location(latitude, longitude, timezone);
    solar_position.set_datetime(year, month, day, hour, 0, 0);
    solar_position.set_method(SolTrace::Data::SolarPositionCalculationMethod::SPA);
    
    double azimuth, elevation;
    solar_position.get_azimuth_elevation(&azimuth, &elevation);

    // Expected values - confirmed with https://gml.noaa.gov/grad/solcalc/azel.html
    double expected_azimuth = 178.61128380;
    double expected_elevation = 73.439035265;
    EXPECT_NEAR(azimuth, expected_azimuth, 1e-3);
    EXPECT_NEAR(elevation, expected_elevation, 1e-3);

    double zenith;
    solar_position.get_azimuth_zenith(&azimuth, &zenith);
    
    double expected_zenith = 90.0 - expected_elevation;
    EXPECT_NEAR(azimuth, expected_azimuth, 1e-3);
    EXPECT_NEAR(zenith, expected_zenith, 1e-3);

    double sun_x, sun_y, sun_z;
    solar_position.get_sun_vector(&sun_x, &sun_y, &sun_z);

    // Assumes x (+) -> east, y (+) -> north, and z (+) -> elevation
    double expected_sun_x = 0.006908;
    double expected_sun_y = -0.2849516;
    double expected_sun_z = 0.9585169;
    EXPECT_NEAR(sun_x, expected_sun_x, 1e-6);
    EXPECT_NEAR(sun_y, expected_sun_y, 1e-6);
    EXPECT_NEAR(sun_z, expected_sun_z, 1e-6);

    // Testing SPA original method
    solar_position.set_method(SolTrace::Data::SolarPositionCalculationMethod::SPA_ORIGINAL);
    solar_position.get_azimuth_elevation(&azimuth, &elevation);

    // Slight differences in solar angles
    expected_azimuth = 178.614529;
    expected_elevation = 73.443097;
    EXPECT_NEAR(azimuth, expected_azimuth, 1e-3);
    EXPECT_NEAR(elevation, expected_elevation, 1e-3);

    solar_position.get_azimuth_zenith(&azimuth, &zenith);

    expected_zenith = 90.0 - expected_elevation;
    EXPECT_NEAR(azimuth, expected_azimuth, 1e-3);
    EXPECT_NEAR(zenith, expected_zenith, 1e-3);

    expected_sun_x = 0.006890135;
    expected_sun_y = -0.284884148;
    expected_sun_z = 0.95853719;

    solar_position.get_sun_vector(&sun_x, &sun_y, &sun_z);
    EXPECT_NEAR(sun_x, expected_sun_x, 1e-6);
    EXPECT_NEAR(sun_y, expected_sun_y, 1e-6);
    EXPECT_NEAR(sun_z, expected_sun_z, 1e-6);

    // Testing SOLPOS method
    solar_position.set_method(SolTrace::Data::SolarPositionCalculationMethod::SOLPOS);
    solar_position.get_azimuth_elevation(&azimuth, &elevation);

    expected_azimuth = 178.6119537;
    expected_elevation = 73.430931;
    EXPECT_NEAR(azimuth, expected_azimuth, 1e-3);
    EXPECT_NEAR(elevation, expected_elevation, 1e-3);

    solar_position.get_azimuth_zenith(&azimuth, &zenith);

    expected_zenith = 90.0 - expected_elevation;
    EXPECT_NEAR(azimuth, expected_azimuth, 1e-3);
    EXPECT_NEAR(zenith, expected_zenith, 1e-3);

    expected_sun_x = 0.006909995;
    expected_sun_y = -0.28508729;
    expected_sun_z = 0.95847666;

    solar_position.get_sun_vector(&sun_x, &sun_y, &sun_z);
    EXPECT_NEAR(sun_x, expected_sun_x, 1e-5);
    EXPECT_NEAR(sun_y, expected_sun_y, 1e-6);
    EXPECT_NEAR(sun_z, expected_sun_z, 1e-6);

    // Testing Legacy method - less accurate
    solar_position.set_method(SolTrace::Data::SolarPositionCalculationMethod::LEGACY);
    solar_position.get_azimuth_elevation(&azimuth, &elevation);

    expected_azimuth = 179.9991897;
    expected_elevation = 73.435378;
    EXPECT_NEAR(azimuth, expected_azimuth, 1e-3);
    EXPECT_NEAR(elevation, expected_elevation, 1e-3);

    solar_position.get_azimuth_zenith(&azimuth, &zenith);

    expected_zenith = 90.0 - expected_elevation;
    EXPECT_NEAR(azimuth, expected_azimuth, 1e-3);
    EXPECT_NEAR(zenith, expected_zenith, 1e-3);

    expected_sun_x = 4.031894e-06;
    expected_sun_y = -0.28509658;
    expected_sun_z = 0.958498794;

    solar_position.get_sun_vector(&sun_x, &sun_y, &sun_z);
    EXPECT_NEAR(sun_x, expected_sun_x, 1e-6);
    EXPECT_NEAR(sun_y, expected_sun_y, 1e-6);
    EXPECT_NEAR(sun_z, expected_sun_z, 1e-6);

    // Testing Duffie method - less accurate -> similar to LEGACY
    solar_position.set_method(SolTrace::Data::SolarPositionCalculationMethod::DUFFIE);
    solar_position.get_azimuth_elevation(&azimuth, &elevation);

    expected_azimuth = 179.110386;
    expected_elevation = 73.43996784;
    EXPECT_NEAR(azimuth, expected_azimuth, 1e-3);
    EXPECT_NEAR(elevation, expected_elevation, 1e-3);

    solar_position.get_azimuth_zenith(&azimuth, &zenith);

    expected_zenith = 90.0 - expected_elevation;
    EXPECT_NEAR(azimuth, expected_azimuth, 1e-3);
    EXPECT_NEAR(zenith, expected_zenith, 1e-3);

    expected_sun_x = 0.00442455;
    expected_sun_y = -0.284985454;
    expected_sun_z = 0.958521629;

    solar_position.get_sun_vector(&sun_x, &sun_y, &sun_z);
    EXPECT_NEAR(sun_x, expected_sun_x, 1e-6);
    EXPECT_NEAR(sun_y, expected_sun_y, 1e-6);
    EXPECT_NEAR(sun_z, expected_sun_z, 1e-6);
}


TEST(SolarPositionCalculator, MultipleLocationsTimesTest) {

    /*locations to test:
    western hemisphere: Golden CO
    eastern hemisphere: Berlin Germany
    southern hemisphere: Chile
    southern hemisphere: Australia
    location near Greenwich meridian with negative longitude and positive time zone: Madrid Spain
    TODO: location near the international dateline with positive longitude and negative time zone: Lomaji, Fiji
    */

    std::vector<double> latitudes = { 39.77, 52.5, -22.44, -26.895, 40.43 };   // , -17.75};      // Would need to fix Duffie method.
    std::vector<double> longitudes = { -105.22, 13.3, -69.42, 130.38, -3.72 }; // , -179.3};
    std::vector<double> time_zones = { -7, 1, -5, 9, 0 };                      // , 12};

    std::vector<int> months = { 3, 6, 8, 11};
    std::vector<int> days = { 5, 15, 25};
    std::vector<int> hours = { 8, 12, 14};
    std::vector<int> minutes = { 0, 15, 30};

    SolarPositionCalculator solar_position;
    double latitude, longitude, timezone;
    double azi_spa, zen_spa;
    double azi_spa_original, zen_spa_original;
    double azi_solpos, zen_solpos;
    double azi_duffie, zen_duffie;
    double azi_legacy, zen_legacy;
    for (size_t i = 0; i < latitudes.size(); i++)
    {
        latitude = latitudes[i];
        longitude = longitudes[i];
        timezone = time_zones[i];
        solar_position.set_location(latitude, longitude, timezone);

        //std::cout << "Location: " << i << std::endl;
        //std::cout << "Latitude: " << latitude << std::endl;
        //std::cout << "Longitude: " << longitude << std::endl;
        //std::cout << "Timezone: " << timezone << std::endl;
        //std::cout << " " << std::endl;
        
        for (size_t j = 0; j < months.size(); j++)
        {
            int month = months[j];
            for (size_t k = 0; k < days.size(); k++)
            {
                int day = days[k];
                for (size_t l = 0; l < hours.size(); l++)
                {
                    int hour = hours[l];
                    for (size_t m = 0; m < minutes.size(); m++)
                    {
                        int minute = minutes[m];

                        //std::cout << "Date: " << month << "/" << day << "  Time: " << hour << ":" << minute << std::endl;

                        solar_position.set_datetime(2024, month, day, hour, minute, 0);
                        solar_position.set_method(SolTrace::Data::SolarPositionCalculationMethod::SPA);
                        EXPECT_NO_THROW(solar_position.get_azimuth_zenith(&azi_spa, &zen_spa));

                        // Compare other methods to SPA
                        solar_position.set_method(SolTrace::Data::SolarPositionCalculationMethod::SPA_ORIGINAL);
                        EXPECT_NO_THROW(solar_position.get_azimuth_zenith(&azi_spa_original, &zen_spa_original));
                        EXPECT_NEAR(azi_spa, azi_spa_original, 0.1);
                        EXPECT_NEAR(zen_spa, zen_spa_original, 0.6);

                        solar_position.set_method(SolTrace::Data::SolarPositionCalculationMethod::SOLPOS);
                        EXPECT_NO_THROW(solar_position.get_azimuth_zenith(&azi_solpos, &zen_solpos));
                        EXPECT_NEAR(azi_spa, azi_solpos, 0.1);
                        EXPECT_NEAR(zen_spa, zen_solpos, 0.6);

                        solar_position.set_method(SolTrace::Data::SolarPositionCalculationMethod::DUFFIE);
                        EXPECT_NO_THROW(solar_position.get_azimuth_zenith(&azi_duffie, &zen_duffie));
                        // Fixes issue of close to north.
                        if (azi_spa > 358.0 && azi_duffie < 2.0)
                            azi_duffie += 360.0; 
                        else if (azi_spa < 2.0 && azi_duffie > 358.0)
                            azi_duffie -= 360.0;
                        EXPECT_NEAR(azi_spa, azi_duffie, 2.5);
                        EXPECT_NEAR(zen_spa, zen_duffie, 1.0);  // 9.1

                        //solar_position.set_method(SolTrace::Data::SolarPositionCalculationMethod::LEGACY);
                        //EXPECT_NO_THROW(solar_position.get_azimuth_zenith(&azi_legacy, &zen_legacy));
                        //EXPECT_NEAR(azi_spa, azi_legacy, 10.0);
                        //EXPECT_NEAR(zen_spa, zen_legacy, 9.0);
                    }
                }
            }
        }
    }
    return;
}





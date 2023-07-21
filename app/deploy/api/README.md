## Usage: PySolTrace for Tracking Error Measurements (main-runsoltrace-iterate.py)

Built on top of the existing SolTrace Python API developed by Mike Wagner, this open-source tool feeds in time-series measurements of tracking error and sun positions and uses Monte-Carlo-based ray-tracing to generate flux distribution on the absorber tube and optical performance.

### Inputs
Site information: latitude, longitude, and altitude of the field measurement site


`optics_type` settings: 

`tracker_angle_input` settings

This tutorial will focus on the `field` setting. For this setting, the required input is a dataframe of timestamps and measured trough angle from each sensor. Example input shown below where `R1` denotes Row 1 and `DO`, `Mid`, and `SO` denote sensors in Row 1 at three lateral locations:

  | Index | R1_DO_Tilt_adjusted | R1_Mid_Tilt_adjusted | R1_SO_Tilt_adjusted |
  | ------------- | ------------- | ------------- | ------------- |
  | 2023-03-05 15:00:00  | 79.20174048386676 | 78.99083190802565 | 79.1755550947139
  | 2023-03-05 19:00:00  | 16.83815731068957 | 17.083237232132948 | 17.09770307440426
  | 2023-03-05 23:00:00  | -56.59770091821414 | -56.67091899881672 | -56.956761552518905

Trough angle sign convention: 0 degrees is flat, -90 is pointing west, +90 is pointing east

### Execution
The code then uses NREL's Solar Position Algorithm (https://pvlib-python.readthedocs.io/en/v0.4.2/generated/pvlib.solarposition.spa_python.html) to calculate the sun position (elevation and azimuth angles) at each timestamp in the input dataframe.


To calculate the optical performance and flux distribution, an instance of SolTrace via the Python API (https://github.com/brookeslawski/SolTrace/blob/develop/app/deploy/api/pysoltrace.py) is created for every time stamp (every row in the input dataframe) and every sensor location (every column in the input dataframe).


SolTrace models the ray interactions with each element in the CSP system using Monte-Carlo based ray-tracing (Wendelin, T. (2003). "SolTRACE: A New Optical Modeling Tool for Concentrating Solar Optics." Proceedings of the ISEC 2003: International Solar Energy Conference, 15-18 March 2003, Kohala Coast, Hawaii. New York: American Society of Mechanical Engineers, pp. 253-260; NREL Report No. CP-550-32866.). We use the back-end ray-tracing algorithm called `coretrace` through a Python API in a packaged called `pysoltrace`, which does not rely on a GUI.

### Outputs
For each SolTrace instance, a results dataframe is generated that includes the intercept factor, the flux distribution at the centerline (circumference of the absorber tube at the lateral mid-point), and the coefficient of variation. For example:

  R1_DO:
  | Index | intercept_factor | flux_centerline | coeff_var |
  | ------------- | ------------- | ------------- | ------------- |
  | 2023-03-05 15:00:00  | 0.862069 | [xx, ..., xx(nx)] | 3.248931
  | 2023-03-05 19:00:00  | 0.764706 | [xx, ..., xx(nx)] | xx
  | 2023-03-05 23:00:00  | 0.879518 | [xx, ..., xx(nx)] | xx


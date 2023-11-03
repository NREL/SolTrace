## Usage: PySolTrace for Timeseries Data of Tilt Angle Measurements in Parabolic Trough Systems

Built on top of the existing SolTrace Python API, this tool feeds in time-series measurements of trough tilt angle and sun positions and uses Monte-Carlo-based ray-tracing to generate flux distribution on the absorber tube and optical performance.


`example_run.py` is an illustrative example of this capability using field measurement data from the Nevada Solar One CSP plant [1].

### Inputs
Parabolic trough collector geometry: module length, aperture width, focal length, absorber tube diameter
Site information: latitude, longitude, and altitude of the field measurement site
Pandas DataFrame containing timeseries data of trough tilt angle for each location along the trough. 


The example uses `demo_field_data.p`, in which 'R4' denotes Row 4 of NSO's array and 'SO' and 'DO' are two locations at either end of the solar collector assembly.
  | Index | R4_SO_Tilt | R4_DO_Tilt | 
  | ------------- | ------------- | ------------- | 
  | 2023-01-15 16:00:00  | 75.975028 | 76.043902 
  | 2023-01-15 16:30:00  | 68.790124 | 68.915574
  | 2023-01-15 17:00:00  | 61.077125 | 61.154348
  | ... | ... | ...

 
Trough angle sign convention: 0 degrees is flat, -90 is pointing west, +90 is pointing east.

### Execution
The code then uses NREL's Solar Position Algorithm (https://pvlib-python.readthedocs.io/en/v0.4.2/generated/pvlib.solarposition.spa_python.html) to calculate the sun position (elevation and azimuth angles) at each timestamp in the input dataframe.


To calculate the optical performance and flux distribution, an instance of SolTrace via the Python API (https://github.com/brookeslawski/SolTrace/blob/develop/app/deploy/api/pysoltrace.py) is created for every time stamp (every row in the input dataframe) and every sensor location (every column in the input dataframe).


SolTrace models the ray interactions with each element in the CSP system using Monte-Carlo based ray-tracing (Wendelin, T. (2003). "SolTRACE: A New Optical Modeling Tool for Concentrating Solar Optics." Proceedings of the ISEC 2003: International Solar Energy Conference, 15-18 March 2003, Kohala Coast, Hawaii. New York: American Society of Mechanical Engineers, pp. 253-260; NREL Report No. CP-550-32866.). We use the back-end ray-tracing algorithm called `coretrace` through a Python API in a packaged called `pysoltrace`, which does not rely on a GUI.

### Outputs
For each sensor location (e.g. R4_SO), a results dataframe is generated that includes the intercept factor, the flux distribution at the centerline (circumference of the absorber tube at the lateral mid-point), and the coefficient of variation. For example:

  R4_SO:
  | Index | intercept_factor | flux_centerline | coeff_var |
  | ------------- | ------------- | ------------- | ------------- |
  | 2023-01-15 16:00:00  | 0.932097 | [xx, ..., xx(nx)] | 1.070096
  | 2023-01-15 16:30:00  | 0.930220 | [yy, ..., yy(nx)] | 1.068185
  | 2023-01-15 17:00:00  | 0.931526 | [zz, ..., zz(nx)] | 1.059448


### References
[1] National Renewable Energy Laboratory (NREL). (2021). Wind and Structural Loads on Parabolic Trough Solar Collectors at Nevada Solar One [data set].  Retrieved from https://dx.doi.org/10.25984/2001061.

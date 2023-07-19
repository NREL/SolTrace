# README
## Instructions for using PySolTrace for Tracking Error Measurements (main-runsoltrace-iterate.py)

Built on top of the existing SolTrace Python API developed by Mike Wagner, this open-source tool feeds in time-series measurements of tracking error and sun positions and uses Monte-Carlo-based ray-tracing to generate flux distribution on the absorber tube and optical performance.

### Inputs
Site information
Optics settings
Dataframe of timestamps and measured tracking error from each sensor. Example input shown below:

  | Index | R1_DO_Tilt_adjusted | R1_Mid_Tilt_adjusted | R1_SO_Tilt_adjusted |
  | ------------- | ------------- | ------------- | ------------- |
  | 2023-03-05 15:00:00  | 79.20174048386676 | 78.99083190802565 | 79.1755550947139
  | 2023-03-05 19:00:00  | 16.83815731068957 | 17.083237232132948 | 17.09770307440426
  | 2023-03-05 23:00:00  | -56.59770091821414 | -56.67091899881672 | -56.956761552518905

### Outputs


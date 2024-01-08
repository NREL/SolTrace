#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 11:11:23 2023

@author: bstanisl

This is an illustrative example using a snippet of field measurement data (`demo_field_data.py`) from the Nevada Solar One CSP plant [1].

[1] National Renewable Energy Laboratory (NREL). (2021). Wind and Structural Loads on Parabolic Trough Solar Collectors at Nevada Solar One [data set].  Retrieved from https://dx.doi.org/10.25984/2001061.

"""
import pandas as pd

from st_processing_functions import *
from run_pysoltrace_iterate import *


#% INPUTS ===========================================================================================

# define constant inputs                                                                                                                                                                                                                                                                                                                       
sunshape_flag = False
sfcerr_flag = False
optics_type = 'realistic' # 'ideal'
plot_rays = False # plot rays in ray-tracing simulation using plotly - this can be slow with more than 1e4 rays
save_pickle = False # save results as pickle file
n_hits = 1e4 #1e5 # 5e6 # 1e5 
nx = 30
ny = 30

# parabolic trough geometry definition ================================
# Nevada Solar One Trough Geometry (LS-2)
module_length = 12.0 # module length
aperture_width = 5.0 #5.77 # aperture width
focal_length = 1.49 #1.71 # focal length # this must be correct for results to make sense
d_abstube = 0.07 # diameter of absorber tube
ptc_pos = [0, 0, 0] # x, y, z
ptc_aim = [0, 0, 1] # x, y, z

# Nevada Solar One field site information
lat, lon = 35.8, -114.983 #coordinates of NSO
altitude = 543 #m

# running with field data timeseries =============================================
tracker_angle_input_mode = 'field' # 'validation' 'nominal'
sensorlocs = ['R4_DO'] #,'R4_SO'] # sensor location names
field_data = pd.read_pickle("./example_field_data.p") # previously prepared for this example

if __name__ == "__main__":
    results, df = run_soltrace_iterate(field_data, lat, lon, altitude, tracker_angle_input_mode, sensorlocs, 
                                       module_length, aperture_width, focal_length, d_abstube, 
                                       ptc_pos, ptc_aim, sunshape_flag, sfcerr_flag, 
                                       optics_type, plot_rays, save_pickle, n_hits, nx, ny)

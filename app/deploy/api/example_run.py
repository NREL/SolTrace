#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 11:11:23 2023

@author: bstanisl
"""
import os
os.chdir('/Users/bstanisl/OneDrive - NREL/Documents/seto-csp-project/SolTrace/s_SolTrace_gitclone_10_31_23/SolTrace/app/deploy/api/')
from os.path import exists
from pysoltrace import PySolTrace, Point
import random
import pandas as pd
import copy
import matplotlib.pyplot as plt
from pvlib import solarposition, tracking
import numpy as np
import math
import time
import pickle
# import plotly.express as px 
import plotly.graph_objects as go
import plotly.io as io
io.renderers.default='browser'

from st_processing_functions import *
from run_pysoltrace_iterate import *


#% INPUTS ===========================================================================================

# define constant inputs                                                                                                                                                                                                                                                                                                                       
sunshape_flag = False
sfcerr_flag = False
optics_type = 'realistic' # 'yang' 'realistic' # 'ideal'
plot_rays = False
save_pickle = False
n_hits = 1e3 # 5e6 # 1e5 #1e5 
nx = 30
ny = 30

# parabolic trough geometry definition ================================
# NSO Trough Geometry: using measurements from CAD file from Dave (aka LS-2)
module_length = 12.0 # module length
aperture_width = 5.0 #5.77 # aperture width
focal_length = 1.49 #1.71 # focal length # this must be correct for results to make sense
d_abstube = 0.07 # diameter of absorber tube
ptc_pos = [0, 0, 0] # x, y, z
ptc_aim = [0, 0, 1] # x, y, z

# NSO field site information
lat, lon = 35.8, -114.983 #coordinates of NSO
altitude = 543 #m

# running with field data timeseries =============================================
tracker_angle_input_mode = 'field' # 'validation' 'nominal'
sensorlocs = ['R4_DO','R4_SO'] #,'R1_Mid'] #,'R1_SO'] #,'R1_SO'] #['R1_SO','R1_Mid','R1_DO','R2_SO','R2_Mid','R2_DO','R4_SO','R4_Mid','R4_DO']
field_data = pd.read_pickle("/Users/bstanisl/OneDrive - NREL/Documents/seto-csp-project/SolTrace/s_SolTrace_gitclone_10_31_23/SolTrace/app/deploy/api/demo_field_data.p")

if __name__ == "__main__":
    results, df = run_soltrace_iterate(field_data, lat, lon, altitude, tracker_angle_input_mode, sensorlocs, 
                                       module_length, aperture_width, focal_length, d_abstube, 
                                       ptc_pos, ptc_aim, sunshape_flag, sfcerr_flag, 
                                       optics_type, plot_rays, save_pickle, n_hits, nx, ny)

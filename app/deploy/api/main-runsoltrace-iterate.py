#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file runs pysoltrace.py iteratively through time-series data 
and with field tilt angle data.

Notes:
    - may take too much memory to run in Spyder, so ideally will be launched 
    from a terminal, maybe on HPC resources

Created on Tue Feb 21 13:46:41 2023
@author: bstanisl
"""

import os
os.chdir('/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/')
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

from postprocessing_functions import *
global focal_len

# define constant inputs                                                                                                                                                                                                                                                                                                                       
sunshape_flag = False
sfcerr_flag = False

# NSO Trough Geometry: using measurements from CAD file from Dave
# l_c = 12.0 # module length
# a_w = 5.0 #5.77 # aperture width
# focal_len = 1.49 #1.71 # focal length # this must be correct for results to make sense
# d_abstube = 0.07 # diameter of absorber tube
# abs_height = focal_len - d_abstube/2. # pt on upper?? sfc of abs tube
# ptc_pos = [0, 0, 0] # x, y, z
# ptc_aim = [0, 0, 1] # x, y, z
# abs_aimz = focal_len*2. # 0. ??
# n_hits = 1e5 # 5e6 # 1e5 #1e5    

# Yang et al 2022 geometry
l_c = 7.8 # module length
a_w = 5.0 #5.77 # aperture width
focal_len = 1.84 #1.71 # focal length # this must be correct for results to make sense
d_abstube = 0.07 # diameter of absorber tube
abs_height = focal_len - d_abstube/2. # pt on upper?? sfc of abs tube
ptc_pos = [0, 0, 0] # x, y, z
ptc_aim = [0, 0, 1] # x, y, z
abs_aimz = focal_len*2. # 0. ??
n_hits = 1e5 # 5e6 # 1e5 #1e5    

# data output settings
# mesh definition for flux map
nx = 30
ny = 30
plotrays = False
save_pickle = False
# sampling_rate = 1 #hrs interval between sampling output

tracker_angle_input = 'validation' # 'nominal' # 'field'
sensorlocs = ['validation'] # ['nominal']
num_iters = 3
#sensorlocs = ['R2_Mid'] #,'R2_Mid','R4_Mid']
# sensorlocs = ['R1_SO','R1_Mid','R1_DO']
# sensorlocs = ['R2_SO','R2_Mid','R2_DO']

optics_type = 'yang' # 'realistic' # 'ideal'

#%% optics properties definition
if optics_type == 'realistic':
    refl_rho = 0.9 # 1. # trough reflectivity
    absr_rho = 0. # receiver reflectivity
    absr_alpha = 0.96 # 1. # receiver absorptivity
    # tau = 1. # transmittance of glass envelope
elif optics_type == 'ideal':
    refl_rho = 1. # trough reflectivity
    absr_rho = 0. # receiver reflectivity
    absr_alpha = 1. # receiver absorptivity
    # tau = 1. # transmittance of glass envelope
elif optics_type == 'yang':
    refl_rho = 0.93 # trough reflectivity
    absr_alpha = 0.96 # receiver absorptivity
    absr_rho = 1 - absr_alpha # receiver reflectivity (assumption in soltrace)
    refl_spec = 1.8 # [mrad] 0.2 == default
    # tau = 0.95 # transmittance of glass envelope
#%% load field data
if tracker_angle_input == 'field':
    year   = '*'
    month  = '12' # '01' # '*'
    day    = '16'
    fileres = '1min' # '1min' or '20Hz'
    outres = '0.5H'
    
    path = '/Users/bstanisl/Documents/seto-csp-project/NSO-field-data/' 
    field_data = load_field_data(path, year, month, day, fileres, outres)
# else:
#     sensorlocs = 'nominal'
#%% get sun positions from SPA directly through pvlib
tracker_angle_input = 'validation'
if tracker_angle_input == 'nominal':
    lat, lon = 35.8, -114.983 #coordinates of NSO
    times = pd.date_range('2022-12-16 15:36:00', '2022-12-17 00:00:00',
                          freq='3H') #, tz=tz)
    # times = pd.date_range('2022-12-16 19:31:00', '2022-12-16 19:40:00',
    #                       freq='0.5T') #, tz=tz)
    # times = pd.date_range('2022-12-16 15:36:00', '2022-12-16 20:00:00',
    #                       freq='4H') #, tz=tz)
    
    solpos = solarposition.get_solarposition(times, lat, lon, altitude=543) #, method='nrel_numba')
    # remove nighttime
    solpos = solpos.loc[solpos['apparent_elevation'] > 0, :]
    
    plt.figure(dpi=250)
    plt.plot(solpos.apparent_elevation,'ko', label='py-spa')
    # plt.plot(spa_sun_positions,'rx', label='SPA txt file')
    plt.ylabel('elevation angle [deg]')
    plt.xticks(rotation=45)
    plt.legend()

    # conclusion: python wrapper generates the same angles as the SPA website
    
    #% calc sun position based on sun vector
    # if tracker_angle_input == 'nominal':
    # solpos['sun_pos_x'] = 1000 * np.sin(np.radians(solpos.azimuth))
    # solpos['sun_pos_y'] = 1000 * np.sin(np.radians(solpos.azimuth))/np.tan(np.radians(solpos.azimuth)) #y
    # solpos['sun_pos_z'] = 1000 * np.tan(np.radians(solpos.apparent_elevation)) #z
    # solpos['sun_pos_x'] = 1000 * np.cos(np.radians(solpos.apparent_elevation)) * np.sin(np.radians(solpos.azimuth))
    # solpos['sun_pos_y'] = 1000 * np.cos(np.radians(solpos.apparent_elevation)) * np.cos(np.radians(solpos.azimuth))
    # solpos['sun_pos_z'] = 1000 * np.sin(np.radians(solpos.apparent_elevation)) #z
    # print(solpos['sun_pos_x'],solpos['sun_pos_y'],solpos['sun_pos_z'])

    [a, b, c] = get_aimpt_from_sunangles(solpos.apparent_elevation, solpos.azimuth)
    solpos['sun_pos_x'] = 1000 * a
    solpos['sun_pos_y'] = 1000 * b
    solpos['sun_pos_z'] = 1000 * c
    
elif tracker_angle_input == 'validation':
    # solpos = pd.DataFrame()
    solpos = pd.DataFrame([[0., 0., 100.]], columns=['sun_pos_x', 'sun_pos_y', 'sun_pos_z'])
    # solpos['sun_pos_x'] = 0. # 1000 * a
    # solpos['sun_pos_y'] = 0. # 1000 * b
    # solpos['sun_pos_z'] = 100. # 1000 * c
# plot_sun_position(solpos)

fig = plt.figure(dpi=250)
plt.plot(solpos['sun_pos_x'],solpos['sun_pos_z'],'ko')
# plt.plot(spa_sun_positions,'rx', label='SPA txt file')
plt.xlabel('sun position [x]')
plt.ylabel('sun position [z]')
# plt.xticks(rotation=45)

#%% calc nominal trough angles
if (tracker_angle_input == 'nominal') or (tracker_angle_input == 'field'):
    trough_angles = pd.DataFrame()
    trough_angles = sun_elev_to_trough_angles(solpos.apparent_elevation,solpos.azimuth)
    trough_angles = trough_angles.to_frame(name='nom_trough_angle')
    anglesdf = solpos.merge(trough_angles, left_index = True, right_index = True, how='inner')
else: # validation
    nom_trough_angle = 0.
    # trough_angles = pd.DataFrame()
    # trough_angles['nom_trough_angle'] = 90.
    # anglesdf = solpos.merge(trough_angles, left_index = True, right_index = True, how='inner')

#%% calculate trough angle deviation
if tracker_angle_input == 'field':
    # merge field data and nominal
    fulldata = anglesdf.merge(field_data, left_index = True, right_index = True, how='inner')
    
    #% calc angle deviation
    for sensorloc in sensorlocs:
        for column in fulldata.filter(regex='Tilt').columns:
            # absolute value
            fulldata['trough_angle_dev_{}'.format(column[0:6])] = abs(fulldata[column] - fulldata['nom_trough_angle'])
        plot_sun_trough_deviation_angles(fulldata, sensorloc)
elif tracker_angle_input == 'validation':
    # fulldata = pd.DataFrame()
    error = np.linspace(0,1,num_iters) # 0.05 #0.025 # [deg]
    data = nom_trough_angle + error
    trough_angles = pd.DataFrame(data, columns=['trough_angle'])
    fulldata = solpos.merge(trough_angles, how='cross') # repeat same sun position for all rows
    fulldata['nom_trough_angle'] = nom_trough_angle
    fulldata['trough_angle_dev'] = error
else: # 'nominal'
    fulldata = anglesdf

#%% main loop
circumf = math.pi*d_abstube
x = np.linspace(-circumf/2.,circumf/2., nx)
y = np.linspace(-l_c/2., l_c/2., ny)

results = {}

tstart = time.time()
# iterate through pandas dataframe at all sun positions
if __name__ == "__main__":
    # iterate over each row and sensor location
    for sensorloc in sensorlocs:
        print(sensorloc)
        
        # initialize
        intercept_factor = [] # np.ones(len(angles))*np.nan
        eta = [] # np.ones(len(angles))*np.nan
        flux_centerline_time = [] # np.ones((nx,len(angles)))*np.nan
        coeff_var = []
        
        for index, row in fulldata.iterrows():
            # set up simulation
            PT = PySolTrace()
            
            # Create two optics types - one for reflector, and one for absorber.
            opt_ref = PT.add_optic("Reflector")
            opt_ref.front.reflectivity = refl_rho # reflects all power
            opt_ref.front.spec_error = refl_spec
            opt_ref.back.reflectivity = refl_rho # reflects all power
            opt_ref.back.spec_error = refl_spec
        
            opt_abs = PT.add_optic("Absorber")
            opt_abs.front.reflectivity = absr_rho
            opt_abs.back.reflectivity = absr_rho
        
            # define sun
            sun = PT.add_sun()
            # Give sun position
            sun.position.x = row['sun_pos_x']
            sun.position.y = row['sun_pos_y']
            sun.position.z = row['sun_pos_z']
            
            # create stage for parabolic trough and absorber
            # (single stage)
            st = PT.add_stage()
            
            # set coordinate system for the PTC and absorber
            # stage origin is global origin = 0, 0, 0
            st.position = Point(0, 0, 0)
            
            # NOMINAL: stage aim towards sun: only take the x and z components of the sun position (trough doesn't adjust laterally)
            # st.aim = Point(sun_position[0], 0, sun_position[2]) #Point(1, 0, 1)
            
            # stage aim towards pt based on elev angle
            # stage_aim = get_aimpt_from_sunangles(row.apparent_elevation, row.azimuth) #, 10)
            # # stage_aim = get_aimpt_from_sunangles_pvlib(solpos.zenith[col], solpos.azimuth[col], focal_len)
            # st.aim = Point(stage_aim[0], 0, stage_aim[1])
            
            if tracker_angle_input == 'nominal':
                # stage aim as tracker angle
                stage_aim = get_aimpt_from_trough_angle(row.nom_trough_angle)
            elif tracker_angle_input == 'validation':
                stage_aim = get_aimpt_from_trough_angle(row.trough_angle)
            else: # 'field'
                # stage aim using actual tracker angle from field
                devkey = [col for col in fulldata.filter(regex='Tilt').columns if sensorloc in col]
                stage_aim = get_aimpt_from_trough_angle(row[devkey[0]])
            st.aim = Point(stage_aim[0], 0, stage_aim[1])
            
            # create parabolic trough element
            el = st.add_element()
            el.optic = opt_ref
            
            # # sets origin of element
            el.position = Point(ptc_pos[0], ptc_pos[1], ptc_pos[2])
            
            # the aim point - sets point to which the element z axis points to
            el.aim = Point(ptc_aim[0], ptc_aim[1], ptc_aim[2])
            
            # define PTC surface
            el.surface_parabolic(focal_len, float('infinity')) # focal_len_y should be 'inf', but that threw an error
            el.aperture_rectangle(a_w,l_c)
            
            # create absorber tube
            #sta = PT.add_stage()
            # single stage
            ela = st.add_element()
            
            #absorber element position
            abs_pos = Point(ptc_pos[0], ptc_pos[1], abs_height)
            ela.position = abs_pos
            ela.aim.z = abs_aimz
            
            # absorber optics
            ela.optic = opt_abs
            ela.surface_cylindrical(d_abstube/2.)
            ela.aperture_singleax_curve(0.,0.,l_c)
            
            # set simulation parameters
            PT.num_ray_hits = n_hits #10 # 1e5
            PT.max_rays_traced = PT.num_ray_hits*100 # 1000 #100 #PT.num_ray_hits*100
            PT.is_sunshape = sunshape_flag # True
            PT.is_surface_errors = sfcerr_flag # True
            
            # run pysoltrace at this instant in time / at this sun position
            # if __name__ == "__main__":
                
            # PT.write_soltrace_input_file('trough-singlestage-stagerotate.stinput')
            PT.run(10, False, 4)
            print("Num rays traced: {:d}".format(PT.raydata.index.size))
            
            # save data frame
            df = PT.raydata
            
            # calculate intercept factor
            intercept_factor.append(calc_intercept_factor(df))
            
            # calculate optical efficiency
            
            # calculate flux map
            ppr = PT.powerperray
            df_rec = generate_receiver_dataframe(df,d_abstube,focal_len)
            flux_st, c_v = compute_fluxmap(ppr,df_rec,d_abstube,l_c,nx,ny,plotflag=True)
            coeff_var.append(c_v)
            
            # save flux centerline for flux map in time
            flux_centerline_time.append(flux_st[:,int(ny/2)])
            
            # plot rays
            if plotrays==True:
                plot_rays_globalcoords(df, PT, st)
            
            #print(df.describe())
        
        #% plot time-varying variables
        # plot_time_series(nominaldf, fulldata, intercept_factor, flux_centerline_time, coeff_var, x)
    
        # combine results into a dataframe
        resultsdf = pd.DataFrame(list(zip(intercept_factor, flux_centerline_time, coeff_var)),
                    index = fulldata.index, 
                    columns =['intercept_factor', 'flux_centerline', 'coeff_var'])
        
        results[sensorloc] = resultsdf
        
        # read pickle file of nominal results if you haven't already
        print('reading nominal results dataframe for comparison')
        nominaldf = pickle.load(open('/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/nominal_12_16_22_1e5.p','rb'))
        
        #% compare nominal to actual
        plot_time_series_compare(nominaldf, fulldata, resultsdf, x, sensorloc)
    
    #%0 save variables to pickle file
    if save_pickle == True:
        pickle.dump([fulldata, results], open('/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/{}_{}.p'.format(tracker_angle_input,n_hits), 'wb'))
    
    if tracker_angle_input == 'field':
        plot_time_series_compare_sensors(nominaldf, fulldata, results, x, sensorlocs)
        
    # delete dataframe to save memory unless it's last iteration
    if (index != fulldata.index[-1]):
        del df
    else:
        tend = time.time()
        elapsed_time = tend - tstart
        print('Execution time: {:2f} seconds'.format(elapsed_time))
    #%% transforming from stage to global
    # plot_rays_globalcoords(df, PT, st)
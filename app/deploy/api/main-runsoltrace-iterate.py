#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Pseudo code for how to run pysoltrace iteratively through time-series data 
and with field data.

This file runs pysoltrace.py with specified inputs.

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
# import plotly.express as px 
import plotly.graph_objects as go
import plotly.io as io
io.renderers.default='browser'

from postprocessing_functions import *
global focal_len

# define constant inputs
refl_rho = 1. # trough reflectivity
absr_rho = 0. # receiver reflectivity
n_hits = 1e5 #1e5
sunshape_flag = False
sfcerr_flag = False

# NSO Trough Geometry: using SGX-2 from http://edge.rit.edu/edge/P15484/public/Detailed%20Design%20Documents/Solar%20Trough%20Preliminary%20analysis%20references/Parabolic%20Trough%20Technology.pdf
l_c = 12.0 # module length
a_w = 5.77 # aperture width
focal_len = 1.71 # focal length # this must be correct for results to make sense
d_abstube = 0.07 # diameter of absorber tube
abs_height = focal_len - d_abstube/2. # pt on upper?? sfc of abs tube
ptc_pos = [0, 0, 0] # x, y, z
ptc_aim = [0, 0, 1] # x, y, z
abs_aimz = focal_len*2. # 0. ??

# data output settings
# mesh definition for flux map
nx = 30
ny = 30
plotrays = False
# sampling_rate = 1 #hrs interval between sampling output

tracker_angle_input = 'field' # 'field'
sensorloc = 'R4_Mid'

#%% load field data
if tracker_angle_input == 'field':
    year   = '*'
    month  = '12' # '01' # '*'
    day    = '16'
    fileres = '1min' # '1min' or '20Hz'
    outres = '0.5H'
    
    path = '/Users/bstanisl/Documents/seto-csp-project/NSO-field-data/' 
    field_data = load_field_data(path, year, month, day, fileres, outres)


#%% get sun positions from SPA directly through pvlib
#tz = 'UTC'
lat, lon = 35.8, -114.983 #coordinates of NSO
times = pd.date_range('2022-12-16 15:36:00', '2022-12-17 00:00:00',
                      freq='1H') #, tz=tz)
# times = pd.date_range('2022-12-16 19:31:00', '2022-12-16 19:40:00',
#                       freq='0.5T') #, tz=tz)

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

#%% calc sun position based on sun vector
solpos['sun_pos_x'] = 1000 * np.sin(np.radians(solpos.azimuth))
solpos['sun_pos_y'] = 1000 * np.sin(np.radians(solpos.azimuth))/np.tan(np.radians(solpos.azimuth)) #y
solpos['sun_pos_z'] = 1000 * np.tan(np.radians(solpos.apparent_elevation)) #z

#%% calc nominal trough angles
trough_angles = pd.DataFrame()
trough_angles = sun_elev_to_trough_angles(solpos.apparent_elevation,solpos.azimuth)
trough_angles = trough_angles.to_frame(name='nom_trough_angle')
anglesdf = solpos.merge(trough_angles, left_index = True, right_index = True, how='inner')

#%% merge field data and nominal
if tracker_angle_input == 'field':
    fulldata = anglesdf.merge(field_data, left_index = True, right_index = True, how='inner')
    
    #% calc angle deviation
    for column in fulldata.filter(regex='Tilt').columns:
        # absolute value
        fulldata['trough_angle_dev_{}'.format(column[0:6])] = abs(fulldata[column] - fulldata['nom_trough_angle'])
    plot_sun_trough_deviation_angles(fulldata, sensorloc)

else: # 'nominal'
    fulldata = anglesdf

#%% main loop
circumf = math.pi*d_abstube
x = np.linspace(-circumf/2.,circumf/2., nx)
y = np.linspace(-l_c/2., l_c/2., ny)

# initialize intercept factor
intercept_factor = [] # np.ones(len(angles))*np.nan
eta = [] # np.ones(len(angles))*np.nan
flux_centerline_time = [] # np.ones((nx,len(angles)))*np.nan
coeff_var = []

tstart = time.time()
# iterate through pandas dataframe at all sun positions
#for col in range(sun_positions.shape[1]): # col iterates through through time
for index, row in fulldata.iterrows():
    # set up simulation
    PT = PySolTrace()
    
    # Create two optics types - one for reflector, and one for absorber.
    opt_ref = PT.add_optic("Reflector")
    opt_ref.front.reflectivity = refl_rho # reflects all power

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
    else: # 'field'
        # stage aim using actual tracker angle from field
        devkey = [col for col in fulldata.filter(regex='Tilt').columns if sensorloc in col]
        stage_aim = get_aimpt_from_trough_angle(row[devkey[0]])
    st.aim = Point(stage_aim[0], 0, stage_aim[1])
    
    # create parabolic trough element
    el = st.add_element()
    el.optic = opt_ref
    
    ptc_pos = [0, 0, 0] # x, y, z
    ptc_aim = [0, 0, 1] # x, y, z
    
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
    if __name__ == "__main__":
        
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
            # Data for a three-dimensional line
            # loc_x = df.loc_x.values
            # loc_y = df.loc_y.values
            # loc_z = df.loc_z.values

            # # Plotting with plotly
            # fig = go.Figure(data=go.Scatter3d(x=loc_x, y=loc_y, z=loc_z, mode='markers', marker=dict( size=1, color=df.stage, colorscale='bluered', opacity=0.8, ) ))

            # #for i in range(PT.raydata.index.size): # all rays
            # for i in range(50,100):
            # #for i in range(0,100000,500):
            #     dfr = df[df.number == i]
            #     ray_x = dfr.loc_x 
            #     ray_y = dfr.loc_y
            #     ray_z = dfr.loc_z
            #     raynum = dfr.number
            #     fig.add_trace(go.Scatter3d(x=ray_x, y=ray_y, z=ray_z, mode='lines', line=dict(color='black', width=0.5)))

            # fig.update_layout(showlegend=False)
            # fig.show()
            
            plot_rays_globalcoords(df, PT, st)
        
        #print(df.describe())
        
        # delete dataframe to save memory unless it's last iteration
        if (index != anglesdf.index[-1]):
            del df
        else:
            tend = time.time()
            elapsed_time = tend - tstart
            print('Execution time: {:2f} seconds'.format(elapsed_time))

#%% plot time-varying variables
# plot_time_series(nominaldf, fulldata, intercept_factor, flux_centerline_time, coeff_var, x)

#%% save variables to pickle file
# resultsdf = pd.DataFrame(list(zip(intercept_factor, flux_centerline_time, coeff_var)),
#             index = solpos.index, 
#             columns =['intercept_factor', 'flux_centerline', 'coeff_var'])
# pickle.dump(resultsdf, open('/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/nominal_12_16_22_1e5.p', 'wb'))

#%% read pickle file of nominal results
# nominaldf = pickle.load(open('/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/nominal_12_16_22_1e5.p','rb'))

#%% compare nominal to actual
# plot_time_series_compare(nominaldf, fulldata, intercept_factor, flux_centerline_time, coeff_var, x, sensorloc)

#%% plt sun position and PTC aim points
# xs = np.column_stack((origin, sun_positions))[0]/100
# ys = np.column_stack((origin, sun_positions))[1]/100
# zs = np.column_stack((origin, sun_positions))[2]/100

# #date_to_val = solpos.index.map(pd.Series(data=np.arange(len(solpos)-1), index=solpos.index.values).to_dict())

# fig = go.Figure(data=go.Scatter3d(x=xs, y=ys, z=zs, mode='markers')) #,marker_color=date_to_val))
# #fig.add_trace(go.Scatter3d(x=xs, y=ys, z=zs, mode='lines'))

# xaims, zaims = get_aimpt_from_sunangles(solpos.apparent_elevation, solpos.azimuth, focal_len)
# yaims = np.zeros((np.shape(xaims)))

# fig.add_trace(go.Scatter3d(x=xaims, y=yaims, z=zaims, mode='markers'))

# fig.update_layout() #showlegend=False)
# fig.show()

#%% transforming from stage to global
# plot_rays_globalcoords(df, PT, st)
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
from pvlib import solarposition
import numpy as np
import math
# import plotly.express as px 
import plotly.graph_objects as go
import plotly.io as io
io.renderers.default='browser'

from postprocessing_functions import *
global focal_len

# define constant inputs
refl_rho = 1. # trough reflectivity
absr_rho = 0. # receiver reflectivity
n_hits = 1e4 #1e5
sunshape_flag = False
sfcerr_flag = False

# NSO Trough Geometry: using SGX-2 from http://edge.rit.edu/edge/P15484/public/Detailed%20Design%20Documents/Solar%20Trough%20Preliminary%20analysis%20references/Parabolic%20Trough%20Technology.pdf
l_c = 12.0 # module length
a_w = 5.77 # aperture width
focal_len = 1.71 # focal length # this must be correct for results to make sense
d_abstube = 0.07 # diameter of absorber tube
abs_height = focal_len - d_abstube/2. # pt on upper?? sfc of abs tube
ptc_pos = [0, 0] # x, y 
ptc_aim = [0, 1] # x, z
abs_aimz = focal_len*2. # 0. ??

#%% add nominal trough angles from text file
# spa_sun_positions = get_sun_angles()
# tstart = '2022-12-16 15:00:00.000000' # '2022-12-16 00:00:00.000000-08:00' 
# tend = '2022-12-17 00:00:00.000000'
# spa_sun_positions = spa_sun_positions[tstart:tend]
# #%
# spa_sun_positions = spa_sun_positions.resample('4H').asfreq().interpolate() #1 minute

# plt.figure()
# plt.plot(spa_sun_positions,'.-')

#%% get sun positions from SPA directly through pvlib

#tz = 'UTC'
lat, lon = 35.8, -114.983 #coordinates of NSO
times = pd.date_range('2022-12-16 15:06:00', '2022-12-17 00:00:00',
                      freq='4H') #, tz=tz)
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

# now that we have the angles, find vector pointing toward the sun
sun_vectors = np.zeros((3,len(solpos.apparent_elevation)))
sun_vectors[0,:] = np.sin(np.radians(solpos.azimuth)) #x
sun_vectors[1,:] = np.sin(np.radians(solpos.azimuth))/np.tan(np.radians(solpos.azimuth)) #y
sun_vectors[2,:] = np.tan(np.radians(solpos.apparent_elevation)) #z

origin = np.zeros(3)

#%% 3d plot of sun vectors
# import plotly.express as px 
# import plotly.graph_objects as go
# import plotly.io as io
# io.renderers.default='browser'

# xs = np.column_stack((origin, sun_vectors))[0]
# ys = np.column_stack((origin, sun_vectors))[1]
# zs = np.column_stack((origin, sun_vectors))[2]

# #date_to_val = solpos.index.map(pd.Series(data=np.arange(len(solpos)-1), index=solpos.index.values).to_dict())

# fig = go.Figure(data=go.Scatter3d(x=xs, y=ys, z=zs, mode='markers')) #,marker_color=date_to_val))
# for i in range(len(solpos.apparent_elevation)):
#     xlines = [origin[0], sun_vectors[0,i]]
#     ylines = [origin[1], sun_vectors[1,i]]
#     zlines = [origin[2], sun_vectors[2,i]]
#     fig.add_trace(go.Scatter3d(x=xlines, y=ylines, z=zlines, mode='lines'))

# fig.update_layout(showlegend=False)
# fig.show()

#%% compute sun positions for soltrace input from sun vectors
sun_positions = 1000 * sun_vectors # multiplying by arbitrary distance for soltrace

# testing
#sun_ positions = np.array([0, 0, 100])[:,None]
# sun_positions = np.array([100, 0, 100])[:,None]

#%%
# data output settings
global nx
global ny
# mesh definition for flux map
nx = 30
ny = 30
circumf = math.pi*d_abstube
x = np.linspace(-circumf/2.,circumf/2., nx)
y = np.linspace(-l_c/2., l_c/2., ny)

plotrays = True
# sampling_rate = 1 #hrs interval between sampling output

# initialize intercept factor
intercept_factor = [] # np.ones(len(angles))*np.nan
eta = [] # np.ones(len(angles))*np.nan
flux_centerline_time = [] # np.ones((nx,len(angles)))*np.nan
coeff_var = []

# iterate through pandas dataframe at all sun positions
for col in range(sun_positions.shape[1]):
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
    sun_position = sun_positions[:,col]
    sun.position.x = sun_position[0] #2. #0.
    sun.position.y = sun_position[1]
    sun.position.z = sun_position[2]
    
    # create stage for parabolic trough and absorber
    # (single stage)
    st = PT.add_stage()
    
    # set coordinate system for the PTC and absorber
    # print('st position = {}'.format(st.position))
    # print('st aim = {}'.format(st.aim))
    
    # stage origin is global origin = 0, 0, 0
    st.position = Point(0, 0, 0)
    
    # only take the x and z components of the sun position (trough doesn't adjust laterally)
    st.aim = Point(sun_position[0], 0, sun_position[2]) #Point(1, 0, 1)
    # print('after st position = {}'.format(st.position))
    # print('after st aim = {}'.format(st.aim))
    # st.position.x = 0
    # st.position.y = 0
    
    # st.aim.x = 1
    # st.aim.y = 0
    # st.aim.z = 1
    
    # print(st.aim)
    
    # create parabolic trough element
    el = st.add_element()
    el.optic = opt_ref
    
    # # sets origin of element
    el.position.x = ptc_pos[0]
    el.position.y = ptc_pos[1]
    
    # the aim point - sets point to which the element z axis points to
    # print(el.aim)
    el.aim.x = ptc_aim[0] #0.01
    el.aim.z = ptc_aim[1]
    
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
        
        PT.write_soltrace_input_file('trough-singlestage-stagerotate.stinput')
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
            loc_x = df.loc_x.values
            loc_y = df.loc_y.values
            loc_z = df.loc_z.values

            # Plotting with plotly
            fig = go.Figure(data=go.Scatter3d(x=loc_x, y=loc_y, z=loc_z, mode='markers', marker=dict( size=1, color=df.stage, colorscale='bluered', opacity=0.8, ) ))

            #for i in range(PT.raydata.index.size): # all rays
            for i in range(50,100):
            #for i in range(0,100000,500):
                dfr = df[df.number == i]
                ray_x = dfr.loc_x 
                ray_y = dfr.loc_y
                ray_z = dfr.loc_z
                raynum = dfr.number
                fig.add_trace(go.Scatter3d(x=ray_x, y=ray_y, z=ray_z, mode='lines', line=dict(color='black', width=0.5)))

            fig.update_layout(showlegend=False)
            fig.show()
        
        #print(df.describe())
        
        # delete dataframe to save memory unless it's last iteration
        if (col != sun_positions.shape[1]-1):
            del df

#%% plot time-varying variables
# plot_time_series(solpos, intercept_factor, flux_centerline_time, coeff_var, x)

#%% coordinate rotation and transformation
# xs = np.column_stack((origin, sun_positions))[0]
# ys = np.column_stack((origin, sun_positions))[1]
# zs = np.column_stack((origin, sun_positions))[2]

# #date_to_val = solpos.index.map(pd.Series(data=np.arange(len(solpos)-1), index=solpos.index.values).to_dict())

# fig = go.Figure(data=go.Scatter3d(x=xs, y=ys, z=zs, mode='markers')) #,marker_color=date_to_val))
# fig.add_trace(go.Scatter3d(x=xs, y=ys, z=zs, mode='lines'))

# fig.update_layout(showlegend=False)
# fig.show()

#%% transforming from stage to global
# plot_rays_globalcoords(df, PT, st)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file runs pysoltrace.py iteratively through time-series data 
and with field tilt angle data.

Notes:
    - this code can be run in three modes:
        - tracker_angle_input = 'field' - sets tracker angle from field data over time
        - tracker_angle_input = 'nominal' - sets tracker angle to point the trough towards the sun over time
        - tracker_angle_input = 'validation' - sets tracker angle to vary error at a single sun position
    - may take too much memory to run in Spyder, so ideally will be launched 
    from a terminal, maybe on HPC resources

Created on Tue Feb 21 13:46:41 2023
@author: bstanisl
"""

import os
os.chdir('/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/')
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
global focal_len

#%% INPUTS ===========================================================================================

# define constant inputs                                                                                                                                                                                                                                                                                                                       
sunshape_flag = False
sfcerr_flag = False
optics_type = 'ideal' # 'yang' 'realistic' # 'ideal'
plot_rays = False
save_pickle = True
number_hits = 1e3 # 5e6 # 1e5 #1e5 

# parabolic trough geometry definition ================================
# NSO Trough Geometry: using measurements from CAD file from Dave (aka LS-2)
module_length = 12.0 # module length
aperture_width = 5.0 #5.77 # aperture width
focal_len = 1.49 #1.71 # focal length # this must be correct for results to make sense
d_abstube = 0.07 # diameter of absorber tube
abs_height = focal_len - d_abstube/2. # pt on upper?? sfc of abs tube
ptc_pos = [0, 0, 0] # x, y, z
ptc_aim = [0, 0, 1] # x, y, z
abs_aimz = focal_len*2. # 0. ??
critical_angle_error = 0.79 #[deg] from firstoptic validation dataset
lat, lon = 35.8, -114.983 #coordinates of NSO
altitude = 543 #m
save_path = '/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/'

# running with field data timeseries =============================================
# tracker_angle_input = 'field' # 'validation' 'nominal' # 'field'
# sensorlocs = ['R1_DO','R1_Mid'] #,'R1_SO'] #,'R1_SO'] #['R1_SO','R1_Mid','R1_DO','R2_SO','R2_Mid','R2_DO','R4_SO','R4_Mid','R4_DO']
# times = pd.date_range('2023-03-05 15:00:00', '2023-03-05 23:50:00',freq='4H') # in UTC
# field_data_path = '/Users/bstanisl/Documents/seto-csp-project/NSO-field-data/' 

# running with field data stats ======================================
# tracker_angle_input = 'stats'
# rows = [1, 2, 4]
# sensorlocs = ['DO','Mid','SO']
# adj_flag = False

# running characteristic median diurnal cycle from NSO ======================================
# tracker_angle_input = 'char' # 'validation' 'nominal' # 'field'
# rows = [1]
# sensorlocs = ['all']
# adj_flag = True

# running nominal =============================================
# tracker_angle_input = 'nominal'
# sensorlocs = ['nominal']
# times = pd.date_range('2023-01-01 15:00:00', '2023-01-01 23:50:00',freq='4H') # in UTC

# running for validation ===================================
tracker_angle_input = 'validation'
sensorlocs = ['validation']
num_iters = 1 #1 # number of trough dev angles to evaluate

# data output settings
# mesh discretization on absorber tube for flux map
nx = 30
ny = 30

#%% optics properties definition
refl_rho, absr_alpha, absr_rho, refl_spec = set_optics_props(optics_type)

#%% load field data
if tracker_angle_input == 'field':
    year   = '2023'  # !! fix this to not be hard coded in (just for testing)
    month  = '03' # '12' # '01' # '*'
    day    = '05' # '16'
    fileres = '1min' # '1min' or '20Hz'
    outres = '0.5H'
    
    field_data = load_field_data(field_data_path, year, month, day, fileres, outres)
    
    # field_data = pickle.load(open('/Users/bstanisl/Documents/seto-csp-project/NSO-field-data/field_data_pproc.p','rb'))
    
    # calculate adjusted tilt angle
    
    # select only necessary columns
    # selcols = list(field_data.filter(regex='Tilt_adj|trough_angle_').columns)
    # selcols = list(field_data.filter(regex='Tilt_adj').columns)
    selcols = list(field_data.filter(regex='Tilt$').columns)
    selcols.extend(['wspd_3m','wdir_3m',])
    field_data = field_data[selcols]
# else:
#     sensorlocs = 'nominal'


#% calculate sun positions from SPA directly through pvlib
if (tracker_angle_input == 'nominal') or (tracker_angle_input == 'field'):
    # sample field data at specified times
    field_data = field_data.loc[times]
    sunangles = get_trough_angles_py(times, lat, lon, altitude)
    field_data = field_data.merge(sunangles, left_index = True, right_index = True, how='inner')
    
    # [a, b, c] = get_aimpt_from_sunangles(field_data.apparent_elevation, field_data.azimuth)
    [a, b, c] = get_aimpt_from_sunangles(field_data.apparent_elevation, field_data.azimuth)
    field_data['sun_pos_x'] = 1000 * a
    field_data['sun_pos_y'] = 1000 * b
    field_data['sun_pos_z'] = 1000 * c
    
    fig = plt.figure(dpi=250)
    plt.plot(field_data['sun_pos_x'],field_data['sun_pos_z'],'ko')
    plt.xlabel('sun position [x]')
    plt.ylabel('sun position [z]')
    
    
elif tracker_angle_input == 'stats':
    if adj_flag:
        stats_data = pickle.load(open('/Users/bstanisl/Documents/seto-csp-project/NSO-field-data/tracker_error_stats_adj.p','rb'))
    else:
        stats_data = pickle.load(open('/Users/bstanisl/Documents/seto-csp-project/NSO-field-data/tracker_error_stats_no_adj.p','rb'))
    solpos = pd.DataFrame([[0., 0., 100.]], columns=['sun_pos_x', 'sun_pos_y', 'sun_pos_z'])
    
elif tracker_angle_input == 'char':
    mediandf = pickle.load(open('/Users/bstanisl/Documents/seto-csp-project/NSO-field-data/median_diurnal_cycle_R1_adj.p','rb'))
    mediandf['trough_angle_dev'] = abs(mediandf['trough_angle_dev'])
    # char_data = mediandf[(mediandf.trough_angle_dev > critical_angle_error) & (mediandf.trough_angle_dev < 1.5)]
    char_data = mediandf[mediandf.trough_angle_dev < 1.5]
    solpos = pd.DataFrame([[0., 0., 100.]], columns=['sun_pos_x', 'sun_pos_y', 'sun_pos_z'])
    
elif tracker_angle_input == 'validation':
    # if validating, sun position is directly overhead at arbitrary height of 100 m
    solpos = pd.DataFrame([[0., 0., 100.]], columns=['sun_pos_x', 'sun_pos_y', 'sun_pos_z'])
# plot_sun_position(solpos)

#%% calc nominal trough angles
# if (tracker_angle_input == 'nominal') or (tracker_angle_input == 'field'):
#     trough_angles = pd.DataFrame()
#     trough_angles = sun_elev_to_trough_angles(solpos.apparent_elevation,solpos.azimuth)
#     trough_angles = trough_angles.to_frame(name='nom_trough_angle')
#     anglesdf = solpos.merge(trough_angles, left_index = True, right_index = True, how='inner')
if tracker_angle_input == 'validation':
    nom_trough_angle = 0. # 0 degrees = flat, facing the sun directlly overhead

#%% calculate trough angle deviation
if tracker_angle_input == 'field':
    # merge field data and nominal
    #inputdata = anglesdf.merge(field_data, left_index = True, right_index = True, how='inner')
    # inputdata = solpos.merge(field_data, left_index = True, right_index = True, how='inner')
    inputdata = field_data
    
    #% calc angle deviation
    for sensorloc in sensorlocs:
        # for column in inputdata.filter(regex='Tilt').columns:
        #     # absolute value
        #     inputdata['trough_angle_dev_{}'.format(column[0:6])] = abs(inputdata[column] - inputdata['nom_trough_angle'])
        plot_sun_trough_deviation_angles(inputdata, sensorloc, adj_flag = False)

elif tracker_angle_input == 'stats':
    # select only the rows that match the senslocs requested in the inputs
    dfdata = stats_data.loc[(rows, sensorlocs),:]
    
    # plot that data
    plot_stats_deviation(dfdata)
    
    # concatenate the mean, mean+sigma, mean-sigma, and max into rows of a dataframe for parsing
    inputdata = pd.DataFrame()
    # dfdata = stats_data.loc[(rows, sensorlocs),:]
    colnames = ['absmean','absmean+std','absmean+2std','absmean-std','absmean-2std','absmax']
    for col in colnames:
        if '+std' in col:
            tmpdf = dfdata['absmean']+dfdata['absstd'] # captures 68% of the data
            tmpdf = tmpdf.to_frame()
            tmpdf.columns=[col]
        elif '+2std' in col:
            tmpdf = dfdata['absmean']+2*dfdata['absstd'] # captures 95% of the data
            tmpdf = tmpdf.to_frame()
            tmpdf.columns=[col]
        elif '-std' in col:
            tmpdf = dfdata['absmean']-dfdata['absstd']
            tmpdf = tmpdf.to_frame()
            tmpdf.columns=[col]
        elif '-2std' in col:
            tmpdf = dfdata['absmean']-2*dfdata['absstd'] # captures 95% of the data
            tmpdf = tmpdf.to_frame()
            tmpdf.columns=[col]
        else: #max?
            tmpdf = dfdata[col].to_frame()
        tmpdf.rename(columns={col: 'trough_angle'}, inplace=True)
        tmpdf['stat'] = col
        inputdata = pd.concat([inputdata, tmpdf])
        del tmpdf
    inputdata['nom_trough_angle'] = 0
    
    # combine with sun position
    for col in solpos.columns:
        inputdata[col] = solpos[col][0]
    
    inputdata = inputdata.set_index(['stat'],append=True)

elif tracker_angle_input == 'char':
    inputdata = char_data.rename(columns={'trough_angle_dev': 'trough_angle'})
    for col in solpos.columns:
        inputdata[col] = solpos[col][0]
    inputdata['nom_trough_angle'] = 0.

elif tracker_angle_input == 'validation':
    
    # array of tracking error values
    error = np.linspace(0,2,num_iters) # 0.05 #0.025 # [deg]
    #error = np.array([0.85])
    
    # define trough angle based on tracking error
    data = nom_trough_angle + error
    
    # create dataframe of trough angles
    trough_angles = pd.DataFrame(data, columns=['trough_angle'])
    inputdata = solpos.merge(trough_angles, how='cross') # repeat same sun position for all rows
    inputdata['nom_trough_angle'] = nom_trough_angle
    inputdata['trough_angle_dev'] = error
else: # 'nominal'
    inputdata = anglesdf

# print(inputdata)

#%% main loop of pysoltrace

# debugging
# inputdata = inputdata.iloc[0:3]
# print(inputdata)

circumf = math.pi*d_abstube
x = np.linspace(-circumf/2.,circumf/2., nx)
y = np.linspace(-module_length/2., module_length/2., ny)

results = {}

tstart = time.time()
# iterate through pandas dataframe at all sun positions
if __name__ == "__main__":
    # iterate over each row and sensor location
    for sensorloc in sensorlocs:
        print(sensorloc)
        
        dropcols = [col for col in inputdata.filter(regex='Tilt$').columns if sensorloc not in col]
        sensorinputdata = inputdata.drop(columns=dropcols)
        
        # initialize
        intercept_factor = [] # np.ones(len(angles))*np.nan
        flux_centerline_time = [] # np.ones((nx,len(angles)))*np.nan
        coeff_var = []
        
        for index, row in sensorinputdata.iterrows(): # for testing .iloc[0:1]
            # set up simulation
            PT = PySolTrace()
            
            # Create two optics types - one for reflector, and one for absorber.
            opt_ref = PT.add_optic("Reflector")
            opt_ref.front.reflectivity = refl_rho
            opt_ref.back.reflectivity = refl_rho
            
            if tracker_angle_input == 'validation':
                opt_ref.front.spec_error = refl_spec
                opt_ref.back.spec_error = refl_spec
                # ? these values seem to have no impact on intercept factor
        
            opt_abs = PT.add_optic("Absorber")
            opt_abs.front.reflectivity = absr_rho
            opt_abs.back.reflectivity = absr_rho
        
            # define sun
            sun = PT.add_sun()
            # Give sun position
            sun.position.x = row['sun_pos_x']
            sun.position.y = row['sun_pos_y']
            sun.position.z = row['sun_pos_z']
            print('sun position x = {}'.format(sun.position.x))
            print('sun position z = {}'.format(sun.position.z))
            
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
            elif (tracker_angle_input == 'validation') or (tracker_angle_input == 'char') or (tracker_angle_input == 'stats'):
                stage_aim = get_aimpt_from_trough_angle(row.trough_angle)
            else: # 'field'
                # stage aim using actual tracker angle from field
                # devkey = [col for col in inputdata.filter(regex='Tilt_adjusted').columns if sensorloc in col]
                devkey = [col for col in sensorinputdata.filter(regex='Tilt$').columns if sensorloc in col]
                stage_aim = get_aimpt_from_trough_angle(row[devkey[0]])
                print('stage aim = {}'.format(stage_aim))
            st.aim = Point(stage_aim[0], 0, stage_aim[1])
            
            # create parabolic trough element
            el = st.add_element()
            el.optic = opt_ref
            
            # # sets origin of element
            el.position = Point(ptc_pos[0], ptc_pos[1], ptc_pos[2])
            
            # the aim point - sets point to which the element z axis points to
            el.aim = Point(ptc_aim[0], ptc_aim[1], ptc_aim[2])
            
            # define parabolic trough surface
            el.surface_parabolic(focal_len, float('infinity')) # focal_len_y should be 'inf', but that threw an error
            el.aperture_rectangle(aperture_width,module_length)
            
            # create absorber tube
            # double stage
            # sta = PT.add_stage()
            # ela = sta.add_element()
            
            # single stage
            ela = st.add_element()
            
            #absorber element position
            abs_pos = Point(ptc_pos[0], ptc_pos[1], abs_height)
            ela.position = abs_pos
            ela.aim.z = abs_aimz
            
            # absorber optics
            ela.optic = opt_abs
            ela.surface_cylindrical(d_abstube/2.)
            ela.aperture_singleax_curve(0.,0.,module_length)
            
            # set simulation parameters
            PT.num_ray_hits = number_hits #10 # 1e5
            PT.max_rays_traced = PT.num_ray_hits*100 # 1000 #100 #PT.num_ray_hits*100
            PT.is_sunshape = sunshape_flag # True
            PT.is_surface_errors = sfcerr_flag # True
                
            # PT.write_soltrace_input_file('trough-singlestage-{}.stinput'.format(index))
            PT.run(10, False, 4)
            print("Num rays traced: {:d}".format(PT.raydata.index.size))
            
            # save data frame
            df = PT.raydata
            
            # calculate intercept factor
            intercept_factor.append(calc_intercept_factor(df))
                        
            # calculate flux map
            ppr = PT.powerperray
            df_rec = generate_receiver_dataframe(df,d_abstube,focal_len)
            flux_st, c_v = compute_fluxmap(ppr,df_rec,d_abstube,module_length,nx,ny,plotflag=True)
            coeff_var.append(c_v)
            
            # save flux centerline for flux map in time
            flux_centerline_time.append(flux_st[:,int(ny/2)])
            
            # plot rays
            if plot_rays==True:
                plot_rays_globalcoords(df, PT, st)
            
            #print(df.describe())
            # delete dataframe to save memory unless it's last iteration
            if (index != sensorinputdata.index[-1]):
                del df
            else:
                tend = time.time()
                elapsed_time = tend - tstart
                print('Execution time: {:2f} seconds'.format(elapsed_time))
        
        #% plot time-varying variables
        # plot_time_series(nominaldf, inputdata, intercept_factor, flux_centerline_time, coeff_var, x)
    
        # combine results into a dataframe
        resultsdf = pd.DataFrame(list(zip(intercept_factor, flux_centerline_time, coeff_var)),
                    index = sensorinputdata.index, 
                    columns =['intercept_factor', 'flux_centerline', 'coeff_var'])
        
        combineddf = sensorinputdata.merge(resultsdf, left_index = True, right_index = True, how='inner')
        results[sensorloc] = combineddf
        
        # read pickle file of nominal results if you haven't already
        if tracker_angle_input == 'field':
            print('reading nominal results dataframe for comparison')
            # nominaldf = pickle.load(open('/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/nominal_12_16_22_1e5.p','rb'))
            nomfn = '/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/nominal_{}_{}_{:.0E}hits_{}_optics.p'.format(sensorinputdata.index[0].month,inputdata.index[0].day,int(number_hits),optics_type)
            if exists(nomfn):
                tmp = pickle.load(open(nomfn,'rb'))
                nominaldf = tmp[1]['nominal']
            else:
                # if nominal has not been run already, then fill with NaNs
                nominaldf = resultsdf.copy()
                for col in nominaldf.columns:
                    nominaldf = nominaldf.assign(col=np.nan)
        
        #if (tracker_angle_input != 'stats') and (tracker_angle_input != 'char'):
            #% compare nominal to actual
            plot_time_series_compare(nominaldf, sensorinputdata, resultsdf, x, sensorloc)
    
    #% save variables to pickle file
    if save_pickle == True:
        if tracker_angle_input == 'field':
            pfn = '{}{}_{}_{}_{:.0E}hits_{}_optics.p'.format(save_path,tracker_angle_input,inputdata.index[0].month,inputdata.index[0].day,int(number_hits),optics_type)
        else:
            pfn = '{}{}_{:.0E}hits_{}_optics_2stage.p'.format(save_path,tracker_angle_input,int(number_hits),optics_type)
        pickle.dump([inputdata, results], open(pfn, 'wb'))
        
    if tracker_angle_input == 'field':
        plot_time_series_compare_sensors(nominaldf, inputdata, results, x, sensorlocs)
    elif tracker_angle_input == 'stats':
        plot_stats_intercept_factor(inputdata, resultsdf)

    #%% transforming from stage to global
    # plot_rays_globalcoords(df, PT, st)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 09:24:15 2023

@author: bstanisl
"""
import os
os.chdir('/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/')
# from os.path import exists
import glob
from pysoltrace import PySolTrace, Point
import pandas as pd
import matplotlib.pyplot as plt
from pvlib import solarposition, tracking
import numpy as np
import math
import time
import pickle
import plotly.graph_objects as go
import plotly.io as io
io.renderers.default='browser'

from st_processing_functions import *

def find_year_month_day(times):
    if times.year[0]==times.year[-1]:
        year = str(times.year[0])
    else:
        year = '*'
    
    if times.month[0]==times.month[-1]:
        month = times.strftime('%m')[0]
    else:
        month = '*'
        
    if times.day[0]==times.day[-1]:
        day = times.strftime('%d')[0]
    else:
        day = '*'
    
    return year,month,day

def run_soltrace_iterate(times, latitude, longitude, altitude, field_data_path, tracker_angle_input_mode, sensorlocs,
                         module_length, aperture_width, focal_length, absorber_diameter, ptc_position, ptc_aim, 
                         sunshape_flag=False, sfcerr_flag=False, optics_type='realistic', plot_rays=False,
                         save_pickle=False, number_hits=1e5, nx=30, ny=30, error_angles=np.array([0., 1.25, 2.5])):
    """
    Calculate the optical performance of a series of tracking error inputs of CSP Parabolic Trough systems using the SolTrace Python API [1].
    
    Applicability:
        - for a single trough
        
    
    Parameters
    ----------
    Input Dataframe with the following columns:
        *Data Series Inputs*
        times : pandas.DateTimeIndex in UTC
            e.g. times = pd.date_range('2023-03-05 15:00:00', '2023-03-05 23:50:00',freq='4H')
        latitude : float
        longitude : float
        altitude : float
        field_data_path : str
        tracker_angle_input_mode : str
        sensorlocs : list of str
        
        *Parabolic Trough Geometry*
        module_length : float
            length of parabolic trough / absorber tube along the tube direction
        aperture_width : float
        focal_length : float
        absorber_diameter : float
        ptc_position : list of 3 floats [x, y, z]
        ptc_aim : list of 3 floats [x, y, z]
        
        *SolTrace Inputs*
        sunshape_flag : boolean
        sfcerr_flag : boolean
        optics_type : str
        plot_rays : boolean
        number_hits : float
        nx : int
        ny : int
    
    Returns
    -------
    Dictionary of DataFrames where each DataFrame contains results from one sensor location with the following columns:
        input tilt angle [deg] (e.g. R1_DO_Tilt) - measured PTC tilt angle in transversal plane
        wspd_3m [m/s] - measured inflow wind speed at 3 m height
        wdir_3m [deg] - measured inflow wind direction at 3 m height
        apparent_zenith [deg]
        apparent_elevation [deg]
        azimuth [deg]
        nom_trough_angle [deg] - calculated from sun position
        sun_pos_x [m]
        sun_pos_y [m]
        sun_pos_z [m]
        track_err [deg]
        intercept_factor [-] - num of rays that hit absorber / num of rays that hit collector
        flux_centerline [W/m2] - heat flux values at lateral midpoint of absorber tube
        coeff_var [-] - a measure of heat distribution on absorber (np.std(flux_st)/np.mean(flux_st))
    
    References
    ----------
    [1] M. Wagner (2023) https://github.com/NREL/SolTrace/tree/develop/app/deploy/api
    """
    
    abs_height = focal_length - absorber_diameter/2. # pt on upper?? sfc of abs tube
    abs_aimz = focal_length*2. # 0. ??
    circumf = math.pi*absorber_diameter
    x = np.linspace(-circumf/2.,circumf/2., nx)
    y = np.linspace(-module_length/2., module_length/2., ny)
    
    refl_rho, absr_alpha, absr_rho, refl_spec = set_optics_props(optics_type)
    
    if tracker_angle_input_mode == 'field':
        
        year, month, day = find_year_month_day(times)
        # year   = '2023'  # !! fix this to not be hard coded in (just for testing)
        # month  = '03' # '12' # '01' # '*'
        # day    = '05' # '16'
        fileres = '1min' # '1min' or '20Hz'
        outres = '0.5H'
        
        field_data = load_field_data(field_data_path, year, month, day, fileres, outres)
    
        # calculate adjusted tilt angle
        
        # select only necessary columns
        # selcols = list(field_data.filter(regex='Tilt_adj|trough_angle_').columns)
        # selcols = list(field_data.filter(regex='Tilt_adj').columns)
        selcols = list(field_data.filter(regex='Tilt$').columns)
        selcols.extend(['wspd_3m','wdir_3m',])
        field_data = field_data[selcols]
        
        # sample field data at specified times
        field_data = field_data.loc[times]
        # else:
            
    #% calculate sun positions from SPA directly through pvlib
    if (tracker_angle_input_mode == 'nominal') or (tracker_angle_input_mode == 'field'):

        sunangles = get_trough_angles_py(times, latitude, longitude, altitude)
        
        if tracker_angle_input_mode == 'nominal':
            field_data = sunangles
        else:
            field_data = field_data.merge(sunangles, left_index = True, right_index = True, how='inner')
        
        [a, b, c] = get_aimpt_from_sunangles(field_data.apparent_elevation, field_data.azimuth)
        field_data['sun_pos_x'] = 1000 * a
        field_data['sun_pos_y'] = 1000 * b
        field_data['sun_pos_z'] = 1000 * c
        
        fig = plt.figure(dpi=250)
        plt.plot(field_data['sun_pos_x'],field_data['sun_pos_z'],'ko')
        plt.xlabel('sun position [x]')
        plt.ylabel('sun position [z]')
               
        inputdata = field_data
        for sensorloc in sensorlocs:
            plot_sun_trough_deviation_angles(field_data, sensorloc, adj_flag = False)
        
    if tracker_angle_input_mode == 'validation':
        # if validating, sun position is directly overhead at arbitrary height of 100 m
        print('in validation loop')
        #sun straight above
        solpos = pd.DataFrame([[0., 0., 100.]], columns=['sun_pos_x', 'sun_pos_y', 'sun_pos_z'])
        nom_trough_angle = 0. # 0 degrees = flat, facing the sun directlly overhead
    
        # array of tracking error values
        error = error_angles
        # if num_iters == 3:
        #     error = np.array([0., np.degrees(1.733333e-02), 2.5]) # for intercepts of 1, x, and 0
        # else:
        #     error = np.linspace(0,2.5,num_iters) # 0.05 #0.025 # [deg]
            #error = np.array([0.85])
        
        # define trough angle based on tracking error
        data = nom_trough_angle + error
        
        # create dataframe of trough angles
        trough_angles = pd.DataFrame(data, columns=['trough_angle'])
        inputdata = solpos.merge(trough_angles, how='cross') # repeat same sun position for all rows
        inputdata['nom_trough_angle'] = nom_trough_angle
        inputdata['trough_angle_dev'] = error
        
        # sun to the east
        solpos = pd.DataFrame([[100., 0., 0.]], columns=['sun_pos_x', 'sun_pos_y', 'sun_pos_z'])
        nom_trough_angle = 90. # facing directly east

        # define trough angle based on tracking error
        data = nom_trough_angle + error
        
        # create dataframe of trough angles
        trough_angles = pd.DataFrame(data, columns=['trough_angle'])
        inputdata2 = solpos.merge(trough_angles, how='cross') # repeat same sun position for all rows
        inputdata2['nom_trough_angle'] = nom_trough_angle
        inputdata2['trough_angle_dev'] = error
        
        # combined dataframe of sun overhead and also to the east
        inputdata = pd.concat([inputdata, inputdata2])
    
    print('input data ',inputdata)

    results = {}

    tstart = time.time()
    
    # iterate through pandas dataframe at all sun positions
    # iterate over each row and sensor location
    for sensorloc in sensorlocs:
        print(sensorloc)
        
        # remove cols for all other sensors
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
            
            if tracker_angle_input_mode == 'validation':
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
            sun.position.y = 0. # row['sun_pos_y']
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
            # # stage_aim = get_aimpt_from_sunangles_pvlib(solpos.zenith[col], solpos.azimuth[col], focal_length)
            # st.aim = Point(stage_aim[0], 0, stage_aim[1])
            
            if tracker_angle_input_mode == 'nominal':
                # stage aim as tracker angle
                stage_aim = get_aimpt_from_trough_angle(row.nom_trough_angle)
            elif (tracker_angle_input_mode == 'validation') or (tracker_angle_input_mode == 'char') or (tracker_angle_input_mode == 'stats'):
                stage_aim = get_aimpt_from_trough_angle(row.trough_angle)
            else: # 'field'
                # stage aim using actual tracker angle from field
                # devkey = [col for col in sensorinputdata.filter(regex='Tilt_adjusted').columns if sensorloc in col]
                devkey = [col for col in sensorinputdata.filter(regex='Tilt$').columns if sensorloc in col]
                stage_aim = get_aimpt_from_trough_angle(row[devkey[0]])
                print('stage aim = {}'.format(stage_aim))
            st.aim = Point(stage_aim[0], 0, stage_aim[1])
            
            # create parabolic trough element
            el = st.add_element()
            el.optic = opt_ref
            
            # # sets origin of element
            el.position = Point(ptc_position[0], ptc_position[1], ptc_position[2])
            
            # the aim point - sets point to which the element z axis points to
            el.aim = Point(ptc_aim[0], ptc_aim[1], ptc_aim[2])
            
            # define parabolic trough surface
            el.surface_parabolic(focal_length, float('infinity')) # focal_len_y should be 'inf', but that threw an error
            el.aperture_rectangle(aperture_width,module_length)
            
            # create absorber tube
            #sta = PT.add_stage()
            # single stage
            ela = st.add_element()
            
            #absorber element position
            abs_pos = Point(ptc_position[0], ptc_position[1], abs_height)
            ela.position = abs_pos
            ela.aim.z = abs_aimz
            
            # absorber optics
            ela.optic = opt_abs
            ela.surface_cylindrical(absorber_diameter/2.)
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
            df_rec = generate_receiver_dataframe(df,absorber_diameter,focal_length)
            flux_st, c_v = compute_fluxmap(ppr,df_rec,absorber_diameter,module_length,nx,ny,plotflag=True)
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
        if tracker_angle_input_mode == 'field':
            print('reading nominal results dataframe for comparison')
            # nominaldf = pickle.load(open('/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/nominal_12_16_22_1e5.p','rb'))
            nomfn = '/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/nominal_{}_{}_*hits_{}_optics.p'.format(sensorinputdata.index[0].month,sensorinputdata.index[0].day,optics_type)
            if len(glob.glob(nomfn)) > 0:
            # if exists(nomfn):
                tmp = pickle.load(open(glob.glob(nomfn)[-1],'rb'))
                nominaldf = tmp['nominal']
            else:
                # if nominal has not been run already, then fill with NaNs
                nominaldf = resultsdf.copy()
                for col in nominaldf.columns:
                    nominaldf = nominaldf.assign(col=np.nan)
        
        # if (tracker_angle_input_mode != 'stats') and (tracker_angle_input_mode != 'char'):
            #% compare nominal to actual
            plot_time_series_compare(nominaldf, sensorinputdata, resultsdf, x, sensorloc)
    
    #% save variables to pickle file
    if save_pickle == True:
        pickle.dump(results, open('/{}{}_{}_{}_{:.0E}hits_{}_optics.p'.format(save_path,tracker_angle_input_mode,inputdata.index[0].month,inputdata.index[0].day,int(number_hits),optics_type), 'wb'))
        # if tracker_angle_input_mode == 'field':
        #     pickle.dump(results, open('/{}{}_{}_{}_{:.0E}hits_{}_optics_test.p'.format(save_path,tracker_angle_input_mode,inputdata.index[0].month,inputdata.index[0].day,int(number_hits),optics_type), 'wb'))
        # if tracker_angle_input_mode == 'char':
        #     pickle.dump([inputdata, resultsdf], open('{}{}_{:.0E}hits_{}_optics.p'.format(save_path,tracker_angle_input_mode,int(number_hits),optics_type), 'wb'))
        # else:
        #     pickle.dump([inputdata, results], open('{}{}_{:.0E}hits_{}_optics.p'.format(save_path,tracker_angle_input_mode,int(number_hits),optics_type), 'wb'))
    
    if tracker_angle_input_mode == 'field':
        # plot_time_series_compare_sensors(nominaldf, inputdata, results, x, sensorlocs)
        plot_time_series_optical_results(results, nominaldf)
        plot_time_series_fluxmap_results(results, x, nominaldf)
    elif tracker_angle_input_mode == 'stats':
        plot_stats_intercept_factor(resultsdf)
    elif tracker_angle_input_mode == 'nominal':
        plot_time_series_compare_nominal(results, x)
    
    return results, df


#%% INPUTS ===========================================================================================

# define constant inputs                                                                                                                                                                                                                                                                                                                       
sunshape_flag = False
sfcerr_flag = False
optics_type = 'realistic' # 'yang' 'realistic' # 'ideal'
plot_rays = False
save_pickle = True
number_hits = 1e5 # 5e6 # 1e5 #1e5 

# parabolic trough geometry definition ================================
# NSO Trough Geometry: using measurements from CAD file from Dave (aka LS-2)
module_length = 12.0 # module length (doesn't matter much)
aperture_width = 5.0 #5.77 # aperture width
focal_len = 1.49 #1.71 # focal length # this must be correct for results to make sense
d_abstube = 0.07 # diameter of absorber tube
ptc_pos = [0, 0, 0] # x, y, z
ptc_aim = [0, 0, 1] # x, y, z
abs_aimz = focal_len*2. # 0. ??
critical_angle_error_min = 0.79 #[deg] from firstoptic validation dataset
critical_angle_error_max = np.degrees(2.363636e-02) #[deg] from firstoptic validation dataset

# field site definition ================================
lat, lon = 35.8, -114.983 #coordinates of Nevada Solar One
altitude = 543 #m
save_path = '/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/'

# running with field data timeseries =============================================
#path = 'smb://nrel.gov/shared/Wind-data/Restricted/Projects/NSO/Processed_data/'
tracker_angle_input = 'field' # 'validation' 'nominal' # 'field'
# sensorlocs = ['R1_Mid', 'R1_SO', 'R2_SO'] #,'R1_Mid','R1_SO'] #,'R1_SO'] 
sensorlocs = ['R1_DO','R1_Mid','R1_SO','R2_DO','R2_Mid','R2_SO','R4_DO','R4_Mid','R4_SO']
# sensorlocs = ['R1_SO','R1_Mid','R1_DO','R2_SO','R2_Mid','R2_DO','R4_SO','R4_Mid','R4_DO']
# times = pd.date_range('2023-03-05 15:00:00', '2023-03-05 23:50:00',freq='1H') # in UTC
tstart = '2023-01-15 16:00:00' # fulldata.index[0] # '2023-02-11 17:00:00'
tend = '2023-01-15 21:00:00' 
times = pd.date_range(tstart, tend, freq='0.5H') # in UTC
# field_data_path = '/Users/bstanisl/Documents/seto-csp-project/NSO-field-data/' CHANGE THIS TO SERVER
field_data_path = '/Volumes/Processed_data/'
error_angles = []

# running for validation ===================================
# times = [''] # in UTC
# field_data_path = '' 
# tracker_angle_input = 'validation'
# sensorlocs = ['validation']
# error_angles = np.concatenate((np.linspace(0.,critical_angle_error_min, 3), 
#                     np.linspace(critical_angle_error_min+.06, critical_angle_error_max-.05, 7),
#                     np.linspace(critical_angle_error_max, 3., 3)))
# error_angles = np.linspace(0,2.5,3) # 0.05 #0.025 # [deg]

# running  nominal (no tracking error) ===================================
# tracker_angle_input = 'nominal'
# sensorlocs = ['nominal']
# tstart = '2023-01-15 16:00:00' # fulldata.index[0] # '2023-02-11 17:00:00'
# tend = '2023-01-15 21:00:00' 
# times = pd.date_range(tstart, tend, freq='0.5H') # in UTC
# field_data_path = '/Volumes/Processed_data/'
# error_angles = []

# data output settings
# mesh discretization on absorber tube for flux map
nx = 30
ny = 30 #30

if __name__ == "__main__":
    results, df = run_soltrace_iterate(times, lat, lon, altitude, field_data_path, tracker_angle_input, sensorlocs, module_length, aperture_width, focal_len, d_abstube, ptc_pos, ptc_aim, sunshape_flag, sfcerr_flag, optics_type, plot_rays, save_pickle, number_hits, nx, ny, error_angles=error_angles)
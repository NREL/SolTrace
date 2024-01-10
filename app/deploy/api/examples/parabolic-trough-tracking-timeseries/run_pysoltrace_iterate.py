#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 25 09:24:15 2023

@author: bstanisl
"""
import sys
sys.path.insert(0, '../../')

import glob
from pysoltrace import PySolTrace, Point
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import math
import time
# import pickle
# import plotly.graph_objects as go
# import plotly.io as io
# io.renderers.default='browser'

from st_processing_functions import *

def run_soltrace_iterate(tilt_angle_data, latitude, longitude, altitude, tracker_angle_input_mode, sensorlocs,
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
        tilt_angle_data : pandas.DataFrame 
            with pandas.DateTimeIndex as index and 2+ columns of trough tilt angle
        latitude : float
            latitude of the field site
        longitude : float
            longitude of the field site
        altitude : float
            altitude of the field site
        tracker_angle_input_mode : str
            'field' or 'nominal' or 'validation' (should be 'field' for basic application)
        sensorlocs : list of str
            list of sensor locations (should appear in the tilt_angle_data tilt angle column names)
        
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
            'realistic' or 'ideal'
        plot_rays : boolean
        save_pickle : boolean
        number_hits : float
            dictates number of rays used in the SolTrace simulation
        nx : int
            number of grid cells in x to discretize absorber tube for flux map
        ny : int
            number of grid cells in y to discretize absorber tube for flux map
    Returns
    -------
    Dictionary of DataFrames where each DataFrame contains results from one sensor location with the following columns:
        input tilt angle [deg] (e.g. R1_DO_Tilt) - measured PTC tilt angle in transversal plane
        apparent_zenith [deg]
        apparent_elevation [deg]
        azimuth [deg]
        nom_trough_angle [deg] - calculated from sun position
        projected_sun_angle [deg] - projected onto east-west transversal plane (for single-axis tracking)
        sun_pos_x [m]
        sun_pos_y [m]
        sun_pos_z [m]
        intercept_factor [-] - num of rays that hit absorber / num of rays that hit collector
        flux_centerline [W/m2] - heat flux values at lateral midpoint of absorber tube
        coeff_var [-] - a measure of heat distribution on absorber (np.std(flux_st)/np.mean(flux_st))
    
    References
    ----------
    [1] M. Wagner (2023) https://github.com/NREL/SolTrace/tree/develop/app/deploy/api
    [2] Python Wrapper for NREL's Solar Position Algorithm https://pvlib-python.readthedocs.io/en/v0.4.2/generated/pvlib.solarposition.spa_python.html
    """
    
    #=======================================================================
    #% calculate geometric parameters of system  --------------------------------------
    #=======================================================================
    abs_height = focal_length - absorber_diameter/2. # pt on upper?? sfc of abs tube
    abs_aimz = focal_length*2. # 0. ??
    circumf = math.pi*absorber_diameter
    x = np.linspace(-circumf/2.,circumf/2., nx)
    y = np.linspace(-module_length/2., module_length/2., ny)
    
    #=======================================================================
    #% assign optical properties --------------------------------------
    #=======================================================================
    refl_rho, absr_alpha, absr_rho, refl_spec = set_optics_props(optics_type)
    
    #=======================================================================
    #% calculate sun positions  --------------------------------------
    #=======================================================================
    if (tracker_angle_input_mode == 'nominal') or (tracker_angle_input_mode == 'field'):
        times = tilt_angle_data.index
        
        # calculate sun positions from SPA directly through pvlib
        sunangles = get_trough_angles_py(times, latitude, longitude, altitude)
        
        tilt_angle_data = tilt_angle_data.merge(sunangles, left_index = True, right_index = True, how='inner')
        
        [a, b, c] = get_aimpt_from_sunangles(tilt_angle_data.apparent_elevation, tilt_angle_data.azimuth)
        tilt_angle_data['sun_pos_x'] = 1000 * a
        tilt_angle_data['sun_pos_y'] = 1000 * b
        tilt_angle_data['sun_pos_z'] = 1000 * c
        
        fig = plt.figure(dpi=250)
        plt.plot(tilt_angle_data['sun_pos_x'],tilt_angle_data['sun_pos_z'],'ko')
        plt.xlabel('sun position [x]')
        plt.ylabel('sun position [z]')
               
        inputdata = tilt_angle_data
        for sensorloc in sensorlocs:
            plot_sun_trough_deviation_angles(inputdata, sensorloc, adj_flag = False)

    elif tracker_angle_input_mode == 'validation':
        # if validating, sun position is directly overhead at arbitrary height of 100 m
        solpos = pd.DataFrame([[0., 0., 100.]], columns=['sun_pos_x', 'sun_pos_y', 'sun_pos_z'])
        nom_trough_angle = 0. # 0 degrees = flat, facing the sun directlly overhead
        
        # define trough angle based on tracking error
        data = nom_trough_angle + tilt_angle_data
        
        # create dataframe of trough angles
        trough_angles = pd.DataFrame(data, columns=['trough_angle'])
        inputdata = solpos.merge(trough_angles, how='cross') # repeat same sun position for all rows
        inputdata['nom_trough_angle'] = nom_trough_angle
        inputdata['trough_angle_dev'] = tilt_angle_data
        
        # sun to the east
        solpos = pd.DataFrame([[100., 0., 0.]], columns=['sun_pos_x', 'sun_pos_y', 'sun_pos_z'])
        nom_trough_angle = 90. # facing directly east

        # define trough angle based on tracking error
        data = nom_trough_angle + tilt_angle_data
        
        # create dataframe of trough angles
        trough_angles = pd.DataFrame(data, columns=['trough_angle'])
        inputdata2 = solpos.merge(trough_angles, how='cross') # repeat same sun position for all rows
        inputdata2['nom_trough_angle'] = nom_trough_angle
        inputdata2['trough_angle_dev'] = tilt_angle_data
        
        # combined dataframe of sun overhead and also to the east
        inputdata = pd.concat([inputdata, inputdata2])
    
    print('input data ',inputdata)
    
    #=======================================================================
    #% iterate over each row and sensor location  --------------------------------------
    #=======================================================================
    results = {}

    tstart = time.time()
    
    for sensorloc in sensorlocs:
        print(sensorloc)
        
        # remove cols for all other sensors
        dropcols = [col for col in inputdata.filter(regex='Tilt$').columns if sensorloc not in col]
        sensorinputdata = inputdata.drop(columns=dropcols)
        
        #=======================================================================
        #% iterate over each timestep in pandas dataframe --------------------------------------
        #=======================================================================
        
        # initialize
        intercept_factor = [] # np.ones(len(angles))*np.nan
        flux_centerline_time = [] # np.ones((nx,len(angles)))*np.nan
        coeff_var = []
        
        for index, row in sensorinputdata.iterrows(): # for testing .iloc[0:1]
            # set up simulation
            PT = PySolTrace()
            
            # Create two optics types - one for reflector, and one for absorber.
            opt_ref = PT.add_optic("Reflector")
            opt_ref.front.reflectivity = refl_rho # assign optical properties
            opt_ref.back.reflectivity = refl_rho # assign optical properties
            
            if tracker_angle_input_mode == 'validation':
                opt_ref.front.spec_error = refl_spec # assign optical properties
                opt_ref.back.spec_error = refl_spec # assign optical properties
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
            
            # create stage for parabolic trough and absorber
            # (single stage)
            st = PT.add_stage()
            
            # set coordinate system for the PTC and absorber
            # stage origin is global origin = 0, 0, 0
            st.position = Point(0, 0, 0)
            
            if tracker_angle_input_mode == 'nominal':
                # stage aim as tracker angle
                stage_aim = get_aimpt_from_trough_angle(row.nom_trough_angle)
            elif (tracker_angle_input_mode == 'validation') or (tracker_angle_input_mode == 'stats'):
                stage_aim = get_aimpt_from_trough_angle(row.trough_angle)
            else: # tracker_angle_input_mode == 'field'
                # stage aim using measured tracker angle from field
                devkey = [col for col in sensorinputdata.filter(regex='Tilt$').columns if sensorloc in col]
                stage_aim = get_aimpt_from_trough_angle(row[devkey[0]])
                # print('stage aim = {}'.format(stage_aim))
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
            # two stage
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
                
            # PT.write_soltrace_input_file('trough-singlestage-{}.stinput'.format(index)) # for comparison with SolTrace GUI results
            PT.run(10, False, 4)
            print("Num rays traced: {:d}".format(PT.raydata.index.size))
            
            # save data frame of ray info
            df = PT.raydata
            
            # calculate intercept factor
            intercept_factor.append(calc_intercept_factor(df))
                        
            # calculate flux map
            ppr = PT.powerperray
            df_rec = generate_receiver_dataframe(df,absorber_diameter,focal_length)
            flux_st, c_v = compute_fluxmap(ppr,df_rec,absorber_diameter,module_length,nx,ny,plotflag=True)
            
            # calculate coefficient of variation
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
    
        # combine results into a dataframe
        resultsdf = pd.DataFrame(list(zip(intercept_factor, flux_centerline_time, coeff_var)),
                    index = sensorinputdata.index, 
                    columns =['intercept_factor', 'flux_centerline', 'coeff_var'])
        
        combineddf = sensorinputdata.merge(resultsdf, left_index = True, right_index = True, how='inner')
        results[sensorloc] = combineddf
        
        # read pickle file of nominal results if you haven't already
        if tracker_angle_input_mode == 'field':
            print('reading nominal results dataframe for comparison')
            nomfn = './nominal_{}_{}_*hits_{}_optics.csv'.format(sensorinputdata.index[0].month,sensorinputdata.index[0].day,optics_type)
            if len(glob.glob(nomfn)) > 0: # if file exists
                nominaldf = read_results_csv(glob.glob(nomfn)[-1])
            else:
                # if nominal has not been run already, then fill with NaNs
                nominaldf = resultsdf.copy()
                for col in nominaldf.columns:
                    nominaldf = nominaldf.assign(col=np.nan)
        
            plot_time_series_compare(nominaldf, sensorinputdata, resultsdf, x, sensorloc)
    
    #% save variables to pickle file
    # if save_pickle == True:
    #     if tracker_angle_input_mode == 'field':
    #         pickle.dump(results, open('./{}_{}_{}_{:.0E}hits_{}_optics.p'.format(tracker_angle_input_mode,inputdata.index[0].month,inputdata.index[0].day,int(number_hits),optics_type), 'wb'))
    #     else:
    #         pickle.dump(results, open('./{}_{:.0E}hits_{}_optics.p'.format(tracker_angle_input_mode,int(number_hits),optics_type), 'wb'))

    if tracker_angle_input_mode == 'field':
        plot_time_series_optical_results(results, nominaldf)
        plot_time_series_fluxmap_results(results, x, nominaldf)
    elif tracker_angle_input_mode == 'validation':
        plot_validation_intercept_factor(results)
    elif tracker_angle_input_mode == 'nominal':
        plot_time_series_compare_nominal(results, x)
    
    return results, df
    
    
    
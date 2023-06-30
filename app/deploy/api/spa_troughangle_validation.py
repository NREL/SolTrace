#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 20 15:38:07 2023

@author: bstanisl

Compare SPA output to python wrapper version

"""
import os
os.chdir('/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pvlib import solarposition, tracking
from postprocessing_functions import *

lat, lon = 35.8, -114.983 #coordinates of NSO
freq = '10T'

def get_trough_angles_txt(tstart, tend):
    times = pd.date_range(tstart, tend, freq='10T')
    
    fn = '/Users/bstanisl/Documents/seto-csp-project/NSO-field-data/NREL_NSO_meas/trough_angles/sun_angles.txt'
    angles = pd.read_csv(fn, parse_dates={'UTC': [0, 1]}).set_index('UTC')  # ,nrows=200
    angles.iloc[:,-1] = angles.iloc[:,-1].where(angles.iloc[:,-1]>0)
    angles['nom_trough_angle'] = np.degrees( np.arctan2(np.sin(np.radians(angles.iloc[:,-1])), 
                                                    np.sin(np.radians(angles.iloc[:,-2])) ))
    angles.nom_trough_angle = angles.nom_trough_angle.where(angles.nom_trough_angle.isnull()==False, -30)
    angles = -angles + 90
    angles = angles[tstart:tend]
    return angles

def get_trough_angles_py(tstart, tend, lat, lon, freq):
    times = pd.date_range(tstart, tend, freq='10T')
    solpos = solarposition.get_solarposition(times, lat, lon, altitude=543) #, method='nrel_numba')
    # remove nighttime
    # solpos = solpos.loc[solpos['apparent_elevation'] > 0, :]

    angles = pd.DataFrame()
    angles = sun_elev_to_trough_angles(solpos.apparent_elevation,solpos.azimuth)
    angles = angles.to_frame(name='nom_trough_angle')
    anglesdf = solpos.merge(angles, left_index = True, right_index = True, how='inner')
    anglesdf.nom_trough_angle[anglesdf['apparent_elevation'] < 0] = 120
    return anglesdf

tstart = '2023-03-05 15:00:00' # fulldata.index[0] # '2022-12-17 07:59:06.500000' # '2022-12-16 00:00:00.000000-08:00' 
tend = '2023-03-06 00:00:00' # fulldata.index[-1]

angles = {}
angles['spa_txt'] = get_trough_angles_txt(tstart, tend)
angles['spa_py'] = get_trough_angles_py(tstart, tend)
#%%

plt.figure(dpi=250)
for key in angles.keys():
    plt.plot(angles[key].nom_trough_angle,'o-', label=key)
# plt.plot(spa_sun_positions,'rx', label='SPA txt file')
plt.ylabel('trough angle [deg]')
plt.xticks(rotation=45)
plt.legend()
#%%
plt.figure(dpi=250)
plt.plot(angles['spa_txt']['Top. elevation angle (corrected)'],'o-', label='spa_txt')
plt.plot(angles['spa_py']['apparent_elevation'],'o-', label='spa_py')
# plt.plot(spa_sun_positions,'rx', label='SPA txt file')
plt.ylabel('elevation angle [deg]')
plt.xticks(rotation=45)
plt.legend()
#%%
plt.figure(dpi=250)
plt.plot(angles['spa_txt']['Topocentric zenith angle'],'o-', label='spa_txt')
plt.plot(angles['spa_py']['apparent_zenith'],'o-', label='spa_py')
# plt.plot(spa_sun_positions,'rx', label='SPA txt file')
plt.ylabel('zenith angle [deg]')
plt.xticks(rotation=45)
plt.legend()

#%%
plt.figure(dpi=250)
plt.plot(angles['spa_txt']['Top. elevation angle (corrected)'],'o-', label='spa_txt elev')
plt.plot(angles['spa_py']['apparent_elevation'],'o-', label='spa_py elev')
plt.plot(angles['spa_txt']['Topocentric zenith angle'],'o-', label='spa_txt zenith')
plt.plot(angles['spa_py']['apparent_zenith'],'o-', label='spa_py zenith')
# plt.plot(spa_sun_positions,'rx', label='SPA txt file')
plt.ylabel('sun angles [deg]')
plt.xticks(rotation=45)
plt.legend()


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
post-processing the soltrace results for exploratory analysis

Created on Fri Mar 31 16:51:58 2023

@author: bstanisl
"""
import os
os.chdir('/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/')
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
from postprocessing_functions import *

figdir = '/Users/bstanisl/Documents/seto-csp-project/SolTrace/fig/'

tau = 1. # transmittance of glass envelope

l_c = 12.0 # module length
a_w = 5.0 #5.77 # aperture width
aperture_area = a_w * l_c
Ib = 1000 # W/m2

sensorlocs = ['R1_SO','R1_Mid','R1_DO']

optics_type = 'realistic' # 'ideal'

if optics_type == 'ideal':
    refl_rho = 0.9 # 1. # trough reflectivity
    absr_rho = 0. # receiver reflectivity
    absr_alpha = 0.96 # 1. # receiver absorptivity
    tau = 1. # transmittance of glass envelope
else:
    refl_rho = 1. # trough reflectivity
    absr_rho = 0. # receiver reflectivity
    absr_alpha = 1. # receiver absorptivity
    tau = 1. # transmittance of glass envelope

# read in results
tracker_angle_input = 'field' # 'field'
n_hits = 1e5 #1e5
fulldata,results = pickle.load(open('/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/{}_{}.p'.format(tracker_angle_input,n_hits), 'rb'))

#%%
plot_time_series_compare_sensors(nominaldf, fulldata, results, x, sensorlocs)

#%% calculate opt efficiency and heat flux
for sensorloc in sensorlocs:
    results[sensorloc]['eta'] = absr_alpha * tau * refl_rho * results[sensorloc].intercept_factor.values
    # results[sensorloc]['eta'] = results[sensorloc].intercept_factor.values
    results[sensorloc]['Q'] = results[sensorloc]['eta'].values * aperture_area * Ib

#%% plot outputs
fig, axs = plt.subplots(4,1,figsize=[9,7],dpi=250,sharex=True)

axs[0].plot(fulldata.apparent_elevation,'k.-')
axs[0].set_ylabel('sun elev. angle [deg]')

for sensorloc in sensorlocs:
    axs[1].plot(results[sensorloc].intercept_factor, '.-',label=sensorloc)
    axs[2].plot(results[sensorloc].eta, '.-',label=sensorloc)
    axs[3].plot(results[sensorloc].Q, '.-',label=sensorloc)

axs[1].set_ylabel('intercept factor')
axs[1].set_ylim([0, 1])
axs[2].set_ylabel('optical efficiency')
axs[2].set_ylim([0, 1])
axs[3].set_ylabel('Q [W]')

plt.legend()
 
for ax in axs:
    ax.tick_params(axis='x',labelrotation=30)

plt.tight_layout()

#%% validate against literature
valdata = {}

filedir = '/Users/bstanisl/Documents/seto-csp-project/SolTrace/StanekFig9.csv'
valdata['stanek-Q'] = pd.read_csv(filedir, header=None, names=['tracker_error', 'Q'])
valdata['stanek_aper_area'] = 1.*1.5

filedir = '/Users/bstanisl/Documents/seto-csp-project/SolTrace/YangFig12-Q-case4.csv'
valdata['yang-Q'] = pd.read_csv(filedir, header=None, names=['tracker_error', 'Q'])
#convert from mrad to degrees
valdata['yang-Q']['tracker_error'] = np.degrees(valdata['yang-Q']['tracker_error']/1000.)
valdata['yang_aper_area'] = 5.*7.8

filedir = '/Users/bstanisl/Documents/seto-csp-project/SolTrace/YangFig12-eta-case4.csv'
valdata['yang-eta'] = pd.read_csv(filedir, header=None, names=['tracker_error', 'eta'])
valdata['yang-eta']['tracker_error'] = np.degrees(valdata['yang-eta']['tracker_error']/1000.)

filedir = '/Users/bstanisl/Documents/seto-csp-project/SolTrace/YangFig10-intercept-case4.csv'
valdata['yang-intc'] = pd.read_csv(filedir, header=None, names=['tracker_error', 'intercept_factor'])
valdata['yang-intc']['tracker_error'] = np.degrees(valdata['yang-intc']['tracker_error']/1000.)

filedir = '/Users/bstanisl/Documents/seto-csp-project/SolTrace/firstOptic_Tracking.dat'
valdata['firstoptic'] = pd.read_csv(filedir, header=None, sep=' ', names=['tracker_error', 'intercept_factor'])
valdata['firstoptic']['tracker_error'] = np.degrees(valdata['firstoptic']['tracker_error'])


#%% plot validation of Q
fig = plt.figure(dpi=250)
plt.scatter(valdata['stanek-Q']['tracker_error'],valdata['stanek-Q']['Q']/valdata['stanek_aper_area'],color='k', label='Stanek et al. 2022')
plt.scatter(valdata['yang-Q']['tracker_error'],valdata['yang-Q']['Q']/valdata['yang_aper_area'],color='k',marker='^', label='Yang et al. 2022')
for sensorloc in sensorlocs:
    devkey = [col for col in fulldata.filter(regex='trough_angle_dev').columns if sensorloc in col]
    plt.scatter(fulldata[devkey],results[sensorloc].Q/aperture_area, label = sensorloc)
plt.xlabel('trough angle deviation [deg]')
plt.ylabel('$Q/A_a$ [W/m2]')
plt.legend()

#%% plot validation of opt efficiency

fig = plt.figure(dpi=250)
plt.scatter(valdata['yang-eta']['tracker_error'],valdata['yang-eta']['eta'],color='k',marker='^', label='Yang et al. 2022')
for sensorloc in sensorlocs:
    devkey = [col for col in fulldata.filter(regex='trough_angle_dev').columns if sensorloc in col]
    plt.scatter(fulldata[devkey],results[sensorloc].eta*100, label = sensorloc)
plt.xlabel('trough angle deviation [deg]')
plt.ylabel('$\eta$ [%]')
plt.legend()

#%% plot validation of intercept factor
plt.rcParams['font.family'] = 'serif'
plt.rcParams['mathtext.fontset'] = 'dejavuserif'
plt.rcParams['font.size'] = 14

resolution_value = 300

fig = plt.figure(dpi=resolution_value)
# plt.scatter(valdata['yang-intc']['tracker_error'],valdata['yang-intc']['intercept_factor'],color='k',marker='x', label='Yang et al. 2022')
plt.plot(valdata['firstoptic']['tracker_error'],valdata['firstoptic']['intercept_factor'],color='k', label='FirstOPTIC')
if tracker_angle_input == 'validation':
    plt.plot(fulldata['trough_angle_dev'],results[sensorloc].intercept_factor, 'kx', label = 'pysoltrace')
else:
    for sensorloc in sensorlocs:
        devkey = [col for col in fulldata.filter(regex='trough_angle_dev').columns if sensorloc in col]
        plt.scatter(fulldata[devkey],results[sensorloc].intercept_factor, color='r', label = sensorloc)
        # plt.scatter(fulldata[devkey],results[sensorloc].intercept_factor, label = sensorloc)
plt.xlabel('trough angle deviation [deg]')
plt.ylabel('intercept factor')
plt.legend()

plt.savefig('{}validation-firstoptic.pdf'.format(figdir), format='pdf', dpi=resolution_value)
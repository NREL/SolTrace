#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
post-processing the soltrace results for exploratory analysis

Created on Fri Mar 31 16:51:58 2023

@author: bstanisl
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import pickle
from postprocessing_functions import *

tau = 1. # transmittance of glass envelope

l_c = 12.0 # module length
a_w = 5.0 #5.77 # aperture width
aperture_area = a_w * l_c
Ib = 1000 # W/m2

# read in results
tracker_angle_input = 'field' # 'field'
n_hits = 1e5 #1e5
fulldata,results = pickle.load(open('/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/{}_{}.p'.format(tracker_angle_input,n_hits), 'rb'))

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
filedir = '/Users/bstanisl/Documents/seto-csp-project/SolTrace/StanekFig9.csv'
stanekdf = pd.read_csv(filedir, header=None, names=['tracker_error', 'Q'])
stanek_aper_area = 1.*1.5

fig = plt.figure(dpi=250)
plt.scatter(stanekdf['tracker_error'],stanekdf['Q']/stanek_aper_area,color='k', label='Stanek et al. 2022')
for sensorloc in sensorlocs:
    devkey = [col for col in fulldata.filter(regex='trough_angle_dev').columns if sensorloc in col]
    plt.scatter(fulldata[devkey],results[sensorloc].Q/aperture_area, label = sensorloc)
plt.xlabel('trough angle deviation [deg]')
plt.ylabel('$Q/A_a$ [W/m2]')
plt.legend()

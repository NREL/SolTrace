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

#%% plotting median diurnal cycle with intercept factor 
combineddf = resultsdf.merge(mediandf, left_index = True, right_index = True, how='outer')
combineddf['intercept_factor'] = combineddf['intercept_factor'].fillna(0)

tilt_col_list = [col for col in data.filter(regex='Tilt_adjusted').columns if srow in col]
markers = ['.','x','+']

plt.rcParams.update({'font.size': 14})
plt.rcParams.update({'font.family': 'serif'})
fig,axs = plt.subplots(5, 1, dpi=250, sharex=True, figsize=[10,12])
axs[0].set_title('median diurnal cycle')

for s, dfs in combineddf.groupby('season'):
    print(s)
    ind = dfs.index.get_level_values(0)
    color = assign_Color(s)
    print(color)
    
    axs[0].plot(ind, dfs.wspd_3m,'.-',linewidth=0.5,color=color,alpha=0.7)
    axs[1].plot(ind, dfs.wdir_3m,'.-',linewidth=0.5,color=color,alpha=0.7)
    #axs[2].plot(ind, dfs.Temp_3m,'.',color=color,alpha=0.7)
            
    axs[2].plot(ind, dfs.nom_trough_angle,"o", fillstyle='none', markeredgecolor='k')
    
    # print(dfs)

    for n, tc in enumerate(tilt_col_list):
        axs[2].plot(ind, dfs[tc], linestyle='-', linewidth=0.5, color=color, marker = markers[n]) #, label=column[:6])
    
        # calc_trough_dev = abs(dfs[tc] - dfs.nom_trough_angle)
        calc_trough_dev = dfs[tc] - dfs.nom_trough_angle
        axs[3].plot(ind, abs(calc_trough_dev), linestyle='-', linewidth=0.5, color=color, marker = markers[n]) #
        
    axs[4].plot(ind, dfs.intercept_factor, linestyle='-', linewidth=0.5, color=color, marker = 'o') #
    
axs[1].axhspan(225, 315, facecolor='0.8', alpha=0.9, label='west')
axs[3].axhline(critical_angle_error, color='0.6', label='critical angle \n deviation')
axs[3].axhline(0, color='k', label='')

axs[0].set_ylabel('wind speed at 3m \n [m/s]')      
axs[1].set_ylabel('wind dir at 3m \n [deg]')      
# axs[2].set_ylabel('air temperature at 3m \n [C]')      
# axs[2].set_ylabel('nom trough angle \n [deg]')
axs[2].set_ylabel('{} trough angle \n [deg]'.format(srow))
axs[3].set_ylabel('abs. val. trough \n angle deviation [deg]')
axs[4].set_ylabel('intercept factor \n $\lambda$')
axs[-1].set_xlabel('hour [UTC]')

axs[3].set_ylim([-2, 5])

axs[0].plot(np.nan, np.nan, label='spring', color='green')
axs[0].plot(np.nan, np.nan, label='summer', color='red')
# axs[0].plot(np.nan, np.nan, label='fall', color='orange')
axs[0].plot(np.nan, np.nan, label='winter', color='blue')

for n,column in enumerate(tilt_col_list):
    axs[2].plot(np.nan, np.nan, label=column[:6], color = 'k', marker=markers[n], linewidth=0.5)
    axs[3].plot(np.nan, np.nan, label=column[:6], color = 'k', marker=markers[n], linewidth=0.5)
axs[4].plot(np.nan, np.nan, label='R1 avg', color = 'k', marker='o', linewidth=0.5)

axs[0].legend(bbox_to_anchor=(1, 1.1), loc='upper left', fontsize=12)
axs[1].legend(bbox_to_anchor=(1, 1.1), loc='upper left', fontsize=12)
axs[2].legend(bbox_to_anchor=(1, 1.1), loc='upper left', fontsize=12)
axs[3].legend(bbox_to_anchor=(1, 1.1), loc='upper left', fontsize=12)
axs[4].legend(bbox_to_anchor=(1, 1.1), loc='upper left', fontsize=12)
    
#%%
noadjfn = '/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/stats_1E+05hits_Dec_June2023_realistic_optics_no_adj.p'
rawdf = pickle.load(open(noadjfn,'rb'))
inputdata = rawdf[0]
resultsdf = rawdf[1]['DO']
plot_stats_intercept_factor(resultsdf)
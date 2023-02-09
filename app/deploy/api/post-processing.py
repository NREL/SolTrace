#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 08:10:44 2023

@author: bstanisl
copied from lines 74+ in https://github.com/NREL/SolarPILOT/blob/develop/deploy/api/test_solarpilot_soltrace.py
and copied from Janna's 2/8/23 email
"""
import numpy as np
import math
import matplotlib.pyplot as plt

"""
soltrace numbering convention
stages:
    1 = collector
    2 = receiver
elements:
    1 = reflected
    -1 = absorbed (isfinal)
"""

#--- Get intersections for stage and element. Note: input stage/element is zero-indexed
def get_intersections(df, stage, elem = 'all'): #, isfinal = False):
# Note: stage/element numbering in SolTrace output starts from 1
# stage : integer, [0:PT.stages-1]
# elem  : string, ['all','absorbed','reflected'] 
    if elem == 'all': # all
        inds = np.where(df['stage'].values == stage)[0]
    
    elif elem == 'absorbed': # absorbed --> negative elem values only
        inds = np.where(np.logical_and(df['stage'].values == stage, df['element'].values < 0))[0]
    
    elif elem == 'reflected':# reflected --> positive elem values only
        inds = np.where(np.logical_and(df['stage'].values == stage, df['element'].values > 0))[0]
    return inds

#--- Get number of ray intersections with given element
def get_number_of_hits(df, stage, elem = 'all') :
    return len(get_intersections(df, stage, elem))

def get_power_per_ray(df):
    #PT = self.PT
    sunstats = PT.sunstats
    ppr = PT.dni * (sunstats['xmax']-sunstats['xmin']) * (sunstats['ymax']-sunstats['ymin']) / sunstats['nsunrays'] # power per ray [W]
    return ppr

#if __name__ == "__main__":
# get ray data
df = PT.raydata # equiv to get_ray_dataframe()
ppr = get_power_per_ray(df)

total_num_rays_traced = PT.raydata.index.size
print('checking that stage rays add up to {} rays traced: {} stage 1 rays + {} stage 2 rays = {} total'.format(
    total_num_rays_traced, get_number_of_hits(df,1), get_number_of_hits(df,2), get_number_of_hits(df,1) + get_number_of_hits(df,2))) 

n_coll_rays = get_number_of_hits(df,1,'all') # 'all' equivalent to PT.num_ray_hits
n_rcvr_rays = get_number_of_hits(df,2,'all') #just absorbed or all rays? 
intercept_factor = n_rcvr_rays/n_coll_rays # * PT.powerperray cancels out in numerator and denominator
print('intercept factor = {} = {}/{}'.format(intercept_factor,n_rcvr_rays_absorbed,n_coll_rays))

# compute optical efficiency (Eq 5, Zhu & Lewandowski 2012)
rho = opt_ref.reflectivity # reflector reflectance
tau = 0.963 # transmittance of receiver glass envelope -- from https://www.nrel.gov/docs/fy09osti/45633.pdf
alpha = 0.96 # average absorptance of receiver surface -- from https://www.nrel.gov/docs/fy09osti/45633.pdf
eta_optical = intercept_factor * rho * tau * alpha
print('optical efficiency = {}'.format(eta_optical))

#%% compute flux distribution
l_c = 12.0 # module length
a_w = 5.77 # aperture width
focal_len = 1.71 # focal length # this must be correct for results to make sense
d_abstube = 0.07 # diameter of absorber tube
abs_height = focal_len + d_abstube/2. # pt on upper?? sfc of abs tube

# copied from lines 74+ in https://github.com/NREL/SolarPILOT/blob/develop/deploy/api/test_solarpilot_soltrace.py
#dfr = df[df.stage==2]
#dfr = dfr[dfr.element==-1]  #absorbed rays
#df_ref = df[df.stage==1] # just reflector stage
#df_ref = dfr[dfr.element==-1]  #absorbed rays
df_rec = df[df.stage==2] # just receiver stage
df_rec = df_rec[df_rec.element==-1]  #absorbed rays

# creating mesh of receiver surface
d_rec = d_abstube
h_rec = abs_height
ny = 20 # how to determine ths?
nx = 20
#y_rec = np.arange(0, h_rec, h_rec/ny)
y_rec = np.arange(-l_c/2., l_c/2., l_c/ny)
x_rec = np.arange(0, np.pi*d_rec, np.pi*d_rec/nx)

Xr,Yr = np.meshgrid(x_rec, y_rec)

# finding rays that intersect a certain point?
tht = 0. #cp.data_get_number(spcxt, "solarfield.0.tht") # what is this?
df_rec['zpos'] = df_rec.loc_z - tht + h_rec/2. # what is this?
df_rec['cpos'] = (np.arctan2(df_rec.loc_x, df_rec.loc_y)+math.pi)*d_rec/2. # and this?

flux_st = np.zeros((ny,nx))
dx = d_rec*np.pi / nx # circumference of receiver / nx
dy = l_c / ny # y is vertical direction here?
anode = dx*dy # area of a node
ppr2 = PT.powerperray / anode *1e-3 #st.powerperray / anode *1e-3 # PT same as st?

# count flux at reflector - needed for intercept factor?

# count flux at receiver
for ind,ray in df_rec.iterrows():

    j = int(ray.cpos/dx)
    i = int(ray.zpos/dy)

    flux_st[i,j] += ppr
print(df_rec.describe())

plt.figure()
plt.title("Flux simulation from the SolTrace engine")
plt.contourf(Xr, Yr, flux_st, levels=25)
plt.colorbar()
plt.title(f"max flux {flux_st.max():.0f} kW/m2, mean flux {flux_st.mean():.1f}")
plt.show()

#title = '$\gamma = {}, \eta = {}$'.format(intercept_factor,eta_optical))

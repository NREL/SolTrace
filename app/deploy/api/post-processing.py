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
#df = PT.raydata # equiv to get_ray_dataframe()
ppr = get_power_per_ray(df)

total_num_rays_traced = PT.raydata.index.size
print('checking that stage rays add up to {} rays traced: {} stage 1 rays + {} stage 2 rays = {} total'.format(
    total_num_rays_traced, get_number_of_hits(df,1), get_number_of_hits(df,2), get_number_of_hits(df,1) + get_number_of_hits(df,2))) 

n_coll_rays = get_number_of_hits(df,1,'reflected') # 'all' equivalent to PT.num_ray_hits
n_rcvr_rays = get_number_of_hits(df,2,'all') #just absorbed or all rays? 
intercept_factor = n_rcvr_rays/n_coll_rays # * PT.powerperray cancels out in numerator and denominator
print('intercept factor = {} = {}/{}'.format(intercept_factor,n_rcvr_rays_absorbed,n_coll_rays))

# compute optical efficiency (Eq 5, Zhu & Lewandowski 2012)
rho = opt_ref.front.reflectivity # reflector reflectance
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

def create_xy_mesh_cyl(d,l,nx,ny):
    # assumes trough is located at 0,0,0
    #xmin = -d/2. # np.min(df.loc_x)
    #xmax = d/2. # np.max(df.loc_x)
    #ymin = -l/2. # np.min(df.loc_y)
    #ymax = l/2 # np.max(df.loc_y)
    #zmin = focal_len - d/2. # np.min(df.loc_z)
    #zmax = focal_len + d/2. #np.max(df.loc_z)
    #print(d, nx)
    x = np.linspace(-d/2., d/2., nx)
    y = np.linspace(-l/2., l/2., ny)
    #print('size of x = ',np.size(x))
    dx = x[1]-x[0]
    dy = y[1]-y[0]
    
    Xc,Yc = np.meshgrid(x,y,indexing='ij')
    print('size of Xc = ', np.size(Xc))
    return Xc,Yc,x,y,dx,dy

def create_polar_mesh_cyl(d,l,nc,ny):
    # creates mesh of unrolled cylinder surface
    circumf = math.pi*d/2.
    c = np.linspace(-circumf/4.,circumf/4., nc)
    y = np.linspace(-l/2., l/2., ny)
    dc = c[1]-c[0] # I think this is correct? dc*dy = area of the node b/c sfc is unrolled
    dy = y[1]-y[0]
    C,Y = np.meshgrid(c,y,indexing='ij')
    return C,Y,c,y,dc,dy
    
def convert_xy_polar_coords(d,x,zloc):
    # position on circumference becomes position in x
    # assumes axis of cylinder is y
    # uses s = r*theta and theta = atan(loc_x/r-loc_z)
    z = zloc - focal_len + d_rec/2. #reset height of tube to z=0
    r = d/2.
    cpos = r * np.arctan2(x,(r-z))
    #print('size of cpos = {}'.format(np.shape(cpos)))
    #print(cpos)
    #print(cpos.values)
    #print('c_pos = r * atan(x/(r-z))) = ')
    #print('{} = {} * atan({}/{})'.format(cpos.values[0],r,x.values[0],r-z.values[0]))
    return cpos

# copied from lines 74+ in https://github.com/NREL/SolarPILOT/blob/develop/deploy/api/test_solarpilot_soltrace.py
#dfr = df[df.stage==2]
#dfr = dfr[dfr.element==-1]  #absorbed rays
#df_ref = df[df.stage==1] # just reflector stage
#df_ref = dfr[dfr.element==-1]  #absorbed rays
df_rec = df[df.stage==2] # just receiver stage
df_rec = df_rec[df_rec.element==-1]  #absorbed rays should be -1? - shouldn't all rays be absorbed?
# is flux calculation based on incident or absorbed rays?

# creating mesh of receiver surface
d_rec = d_abstube
h_rec = abs_height
c_rec = np.pi*d_rec # circumference
ny = 20 # how to determine ths?
nx = 20

# xy coords
#Xc,Yc,x,y,dx,dy = create_xy_mesh_cyl(d_rec,l_c,nx,ny)
#anode = dx*dy # area of a node

# unrolled coords - position on circumference becomes position in x
Xc,Yc,x,y,dx,dy = create_polar_mesh_cyl(d_rec,l_c,nx,ny)
df_rec['ypos'] = df_rec.loc_y
df_rec['cpos'] = convert_xy_polar_coords(d_rec,df_rec.loc_x,df_rec.loc_z)
anode = dx*dy

#%%
flux_st = np.zeros((nx,ny))
#dx = d_rec*np.pi / nx # circumference of receiver / nx
#dy = l_c / ny # y is vertical direction here?

# is there a better way to do this?: seems as though anode shouldn't be in there
ppr = PT.powerperray / anode *1e-3 #st.powerperray / anode *1e-3 # PT same as st?

# count flux at receiver
for ind,ray in df_rec.iterrows():

    #i = int(ray.loc_x/dx)
    #j = int(ray.loc_y/dy)
    #print('ind,i,j = {},{},{}'.format(ind,i,j))
    
    #print(np.argmin(np.abs(x - ray.loc_x)))
    #print(np.where(np.abs(x - ray.loc_x) < tol)[0])
    # choose index of closest location to coordinates
    i = np.argmin(np.abs(x - ray.loc_x)) # int(np.where(np.abs(x - ray.loc_x) < tol)[0])
    j = np.argmin(np.abs(y - ray.loc_y)) # int(np.where(np.abs(y - ray.loc_y) < tol)[0])
    #print('ind,i,j = {},{},{}'.format(ind,i,j))
    
    flux_st[i,j] += ppr
print(df_rec.describe())

plt.figure(figsize=[3,4],dpi=250)
plt.title("Flux simulation from the SolTrace engine")
plt.contourf(Xc, Yc, flux_st, levels=25)
plt.colorbar()
plt.title(f"max flux {flux_st.max():.0f} kW/m2, mean flux {flux_st.mean():.1f}")
plt.xlabel('x [m]')
plt.ylabel('y [m]')
#plt.xlim([-0.5,0.5])
plt.show()

# %%

#y_rec = np.arange(0, h_rec, h_rec/ny)
y_rec = np.arange(-l_c/2., l_c/2., l_c/ny) # centered at y=0
x_rec = np.arange(0, c_rec, c_rec/nx)

Xr,Yr = np.meshgrid(x_rec, y_rec)


# finding rays that intersect a certain point?
tht = 0. #cp.data_get_number(spcxt, "solarfield.0.tht") # height of tower
df_rec['zpos'] = df_rec.loc_z - tht + h_rec/2. # vertical position of receiver
# unwrapped cylinder:
df_rec['cpos'] = (np.arctan2(df_rec.loc_x, df_rec.loc_y)+math.pi)*d_rec/2. # acompute circumferential position 

# count flux at reflector - needed for intercept factor?
for ind,ray in df_rec.iterrows():

    j = int(ray.cpos/dx)
    i = int(ray.zpos/dy)

    flux_st[i,j] += ppr
print(df_rec.describe())


#title = '$\gamma = {}, \eta = {}$'.format(intercept_factor,eta_optical))

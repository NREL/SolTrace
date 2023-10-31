#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 09:22:24 2023

@author: bstanisl
"""

import numpy as np
import math
import matplotlib.pyplot as plt
import pandas as pd
from st_processing_functions import *

"""
Validating flux map from post-processing.py
Against LS-3 parabolic trough standard case from SolTrace
First step: run NSO-trough.py with n_hits = 5e6
"""

#%% from pysoltrace
nx = 30
ny = 30
l_c = 10.0 # module length
a_w = 2.887*2 # aperture width
R = (1/0.292398) # from screenshot surface definition p-1/R
focal_len = R/2. #1.71 # focal length # this must be correct for results to make sense
r_abstube = 1/28.5714
d_abstube = r_abstube*2. #0.07 # diameter of absorber tube

PTppr = PT.powerperray
df_rec = generate_receiver_dataframe(df,d_abstube,focal_len)

#%
# flux_st, c_v = compute_fluxmap(ppr,df_rec,d_abstube,l_c,nx,ny,plotflag=False)
Xc,Yc,x,y,dx,dy,psi = create_polar_mesh_cyl(d_abstube,l_c,nx,ny)

flux_st_bs = np.zeros((nx,ny))
anode = dx*dy
ppr = PTppr / anode *1e-3

# count flux at receiver
for ind,ray in df_rec.iterrows():
    
    # i = int(np.where(np.abs(x - ray.loc_x) < tol)[0])
    # j = int(np.where(np.abs(y - ray.loc_y) < tol)[0])
    # choose index of closest location to coordinates
    i = np.argmin(np.abs(x - ray.cpos)) 
    j = np.argmin(np.abs(y - ray.ypos)) 
    
    flux_st_bs[i,j] += ppr    
print(ray)
# print(flux_st_bs)


#%% from GUI
flux_st = flux_st_bs
loc_val = '/Users/bstanisl/OneDrive - NREL/Documents/seto-csp-project/SolTrace/fluxmap.flx'
dfv = pd.read_csv(loc_val,skiprows=3,delimiter='   ')
# print(dfv)

# from screenshot
xmin = -0.11215496988813
xmax = 0.11215496988813
ymin = -5.099999
ymax = 5.099999
ppr = 0.00577356
nx = 30
ny = 30

gui_binszx = (xmax/1.01 - xmin/1.01)/nx
gui_binszy = (ymax/1.01 - ymin/1.01)/ny
gui_zscale = ppr/(gui_binszx*gui_binszy)/1000.

xv = np.linspace(xmin/1.01,xmax/1.01,nx) # reducing the extended xlimits from gui code
yv = np.linspace(ymin/1.01,ymax/1.01,ny)

vmax = np.max([dfv.max().max(), flux_st.max()])
vmin = np.min([dfv.min().min(), flux_st.min()])
# levels=15

# plt.figure(1,2,figsize=[6,4],dpi=250)
plt.rcParams['font.size'] = 14
fig, axs = plt.subplots(3,1,figsize=[6,10],dpi=300,sharex=True)
cf0 = axs[0].contourf(Xc, Yc, np.transpose(dfv), vmin=vmin, vmax=vmax, cmap='viridis') #, vmin=240, vmax=420)
# fig.colorbar(cf0, ax=axs[0])

cb_ax = fig.add_axes([1,.39,.04,.53])
fig.colorbar(cf0,orientation='vertical',cax=cb_ax,label='flux [kW/m2]')

axs[0].set_title(f"SolTrace GUI: \n max flux {dfv.max().max():.2f} kW/m2, mean flux {dfv.mean().mean():.2f}")

# cf1 = axs[1].pcolormesh(Xc, Yc, flux_st, cmap='viridis') #, shading='gouraud')
cf1 = axs[1].contourf(Xc, Yc, flux_st, vmin=vmin, vmax=vmax, cmap='viridis') #, vmin=240, vmax=420)
# fig.colorbar(cf1, ax=axs[1])
axs[1].set_title(f"pysoltrace: \n max flux {flux_st.max():.2f} kW/m2, mean flux {flux_st.mean():.2f}")

for ax in axs[:2]:
    ax.set_ylabel('y [m]')
# axs[1].set_xlabel('x [m]')
    #plt.xlim([-0.5,0.5])
# plt.show()

#% line plot
flux_centerline = flux_st[:,int(ny/2)]
flux_centerlinev = dfv.values[int(ny/2),:]

# plt.figure(figsize=[3,4],dpi=250)
#plt.title("Flux simulation from the SolTrace engine")
axs[2].plot(x,flux_centerlinev, 'k-', label='SolTrace GUI') #, vmin=240, vmax=420)
axs[2].plot(x,np.asarray(flux_centerline), 'kx', label='pysoltrace') #, vmin=240, vmax=420)
#plt.title(f"max flux {flux_st.max():.0f} kW/m2, mean flux {flux_st.mean():.1f}")
axs[2].set_title('Centerline Comparison')
axs[2].set_xlabel('x [m]')
axs[2].set_ylabel('flux at y=0 [kW/m2]')
axs[2].legend()
#plt.xlim([-0.5,0.5])
plt.tight_layout()
plt.show()

#%% debug plotting
dfr = df
dfr['ypos'] = dfr.loc_y
dfr['cpos'] = convert_xy_polar_coords(d_abstube,dfr.loc_x,dfr.loc_z,focal_len,gui_coords=True)
    
r = r_abstube
nbinsx = nx
nbinsy = ny

minx = -math.pi*r
maxx = math.pi*r

miny = -l_c/2.0
maxy = l_c/2.0

# total grid size
gridszx = 2.0 * math.pi * r
gridszy = maxy - miny

# bin size
binszx = gridszx/nbinsx
binszy = gridszy/nbinsy

RayCount = 0
fluxGrid = np.zeros((nbinsx,nbinsy))

# create grid of ray counts
for index,row in dfr.iterrows():
    x = row['cpos']
    y = row['ypos']
    # print(x)
    
    GridIncrementX = -1
    GridIncrementY = -1
    while ((minx + (GridIncrementX+1)*binszx) < x):
        GridIncrementX += 1
    while ((miny + (GridIncrementY+1)*binszy) < y):
        GridIncrementY += 1
    fluxGrid[GridIncrementX,GridIncrementY] += 1 #if ray falls inside a bin, increment count for that bin
    RayCount += 1 #increment ray intersection counter
print('RayCount = ',RayCount)

#%% calculate midpoints of bins
xValues = np.zeros(nbinsx)
yValues = np.zeros(nbinsy)
for i in range(nbinsx):
    xValues[i] = minx + binszx/2. + i*binszx
for i in range(nbinsy):
    yValues[i] = miny + binszy/2. + i*binszy

#% calculate peak, min, and avg fluxes
# zScale = PTppr/(binszx*binszy) # !! CHECK PTppr vs PowerPerRay
# GUIppr = 0.00577311
# PowerPerRay = PT.dni * (PT.sunstats['xmax'] - PT.sunstats['xmin'])
PowerPerRay = (PT.sunstats['xmax']-PT.sunstats['xmin'])*(PT.sunstats['ymax'] - PT.sunstats['ymin']) / PT.sunstats['nsunrays'] * PT.dni

zscale = PowerPerRay/(binszx*binszy)/PT.dni

PeakFlux = 0.
SumFlux = 0.
MinFlux = 1e99

for i in range(nbinsx):
    for j in range(nbinsy):
        z = fluxGrid[i,j]*zscale
        SumFlux += z

        if z > PeakFlux:
            PeakFlux = z
            NRaysInPeakFluxBin = fluxGrid[i,j]
        elif z < MinFlux:
            MinFlux = z
            NRaysInMinFluxBin = fluxGrid[i,j]
AveFlux = SumFlux/(nbinsx*nbinsy)

#%%
xx = np.linspace(xmin,xmax,nbinsx) # using extended x limits from GUI output
yy = np.linspace(ymin,ymax,nbinsy)
Xc,Yc = np.meshgrid(xx,yy,indexing='ij')
# Xc,Yc = np.meshgrid(xValues,yValues,indexing='ij')

vmax = np.max([dfv.max().max(), PeakFlux])
vmin = np.min([dfv.min().min(), MinFlux])
# levels=15

# plt.figure(1,2,figsize=[6,4],dpi=250)
fig, axs = plt.subplots(2,1,figsize=[6,8],dpi=300,sharex=True)
# cf0 = axs[0].contourf(xv, yv, dfv, levels=50, cmap='viridis') #, vmin=240, vmax=420)
cf0 = axs[0].pcolormesh(xv, yv, dfv/gui_zscale, cmap='viridis') #, vmin=240, vmax=420)
fig.colorbar(cf0, ax=axs[0])
axs[0].set_title(f"SolTrace GUI: \n max flux {dfv.max().max():.2f} kW/m2, mean flux {dfv.mean().mean():.2f}")

# cf1 = axs[1].contourf(Xc, Yc, fluxGrid*zscale, cmap='viridis') #, vmin=240, vmax=420)
# cf1 = axs[1].contourf(Xc, Yc, fluxGrid*zscale, levels=50, cmap='viridis') #, vmin=240, vmax=420)
cf1 = axs[1].pcolormesh(Xc, Yc, fluxGrid, cmap='viridis') #, vmin=240, vmax=420)
fig.colorbar(cf1, ax=axs[1])
axs[1].set_title(f"pysoltrace: \n peak flux {PeakFlux:.2f} kW/m2, avg flux {AveFlux:.2f}")

for ax in axs:
    ax.set_ylabel('y [m]')
axs[1].set_xlabel('x [m]')
    #plt.xlim([-0.5,0.5])
plt.show()

#%% line plot
flux_centerline = fluxGrid[:,int(ny/2)]*zscale
flux_centerlinev = dfv.values[int(ny/2),:]

plt.figure(figsize=[3,4],dpi=250)
#plt.title("Flux simulation from the SolTrace engine")
plt.plot(xv,flux_centerlinev, 'k.-', label='SolTrace GUI') #, vmin=240, vmax=420)
plt.plot(xValues,np.asarray(flux_centerline), 'r.-', label='pysoltrace') #, vmin=240, vmax=420)
#plt.title(f"max flux {flux_st.max():.0f} kW/m2, mean flux {flux_st.mean():.1f}")
plt.title('LS-3 Validation')
plt.xlabel('x [m]')
plt.ylabel('flux at y=0 [kW/m2]')
plt.legend(fontsize=8)
#plt.xlim([-0.5,0.5])
plt.show()

# %% flux_st = compute_fluxmap_solarpilot(ppr,df_rec,d_abstube,l_c,nx,ny)
# # df_rec['ypos'] = df_rec.loc_y
# # df_rec['cpos'] = convert_xy_polar_coords(d_rec,df_rec.loc_x,df_rec.loc_z,focal_len,gui_coords=True)

# # Xc,Yc,x_bs,y_bs,dx_bs,dy_bs,psi = create_polar_mesh_cyl(d_rec,focal_len,nx,ny)

# flux_st_sp = np.zeros((ny,nx))
# # dx = d_rec*np.pi / nx 
# # dy = h_rec / ny
# anode = dx*dy
# # ppr = PTppr / anode *1e-3

# for ind,ray in df_rec.iterrows():
#     i = int(ray.cpos/dx)
#     j = int(ray.ypos/dy)

#     flux_st_sp[i,j] += ppr
# print(ray)
# # print(flux_st_sp)

# #%%
# vmax = np.max([flux_st_bs.max(), flux_st_sp.max()])
# vmin = np.min([flux_st_bs.min(), flux_st_sp.min()])

# fig, axs = plt.subplots(2,1,figsize=[6,8],dpi=300,sharex=True)
# cf0 = axs[0].contourf(Xc, Yc, flux_st_bs, vmin=vmin, vmax=vmax, cmap='viridis') #, vmin=240, vmax=420)
# fig.colorbar(cf0, ax=axs[0])
# axs[0].set_title(f"bjs calculation: \n max flux {flux_st_bs.max():.2f} kW/m2, mean flux {flux_st_bs.mean():.2f}")

# # cf1 = axs[1].pcolormesh(Xc, Yc, flux_st, cmap='viridis') #, shading='gouraud')
# cf1 = axs[1].contourf(Xc, Yc, flux_st_sp, vmin=vmin, vmax=vmax, cmap='viridis') #, vmin=240, vmax=420)
# fig.colorbar(cf1, ax=axs[1])
# axs[1].set_title(f"solarpilot: \n max flux {flux_st_sp.max():.2f} kW/m2, mean flux {flux_st_sp.mean():.2f}")

# for ax in axs:
#     ax.set_ylabel('y [m]')
# axs[1].set_xlabel('x [m]')
#     #plt.xlim([-0.5,0.5])
# plt.show()
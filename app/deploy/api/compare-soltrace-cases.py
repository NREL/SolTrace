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

"""
Validating flux map from post-processing.py
Against LS-3 parabolic trough standard case from SolTrace
"""

loc_val = '/Users/bstanisl/Documents/seto-csp-project/SolTrace/fluxmap.flx'
dfv = pd.read_csv(loc_val,skiprows=3,delimiter='   ')
print(dfv)

# from screenshot
xmin = -0.11215496988813
xmax = 0.11215496988813
ymin = -5.099999
ymax = 5.099999
ppr = 0.00577356
nx = 30
ny = 30

xv = np.linspace(xmin,xmax,nx)
yv = np.linspace(ymin,ymax,ny)

plt.figure(figsize=[6,4],dpi=250)
plt.contourf(xv, yv, dfv, levels=15, cmap='viridis') #, vmin=240, vmax=420)
#plt.pcolormesh(Xc, Yc, flux_st, cmap='plasma', shading='gouraud')
plt.colorbar()
#plt.clim(240,420)
plt.title(f"SolTrace GUI: \n max flux {dfv.max().max():.2f} kW/m2, mean flux {dfv.mean().mean():.2f}")
plt.xlabel('x [m]')
plt.ylabel('y [m]')
#plt.xlim([-0.5,0.5])
plt.show()

#%% line plot
flux_centerline = flux_st[:,int(ny/2)]
flux_centerlinev = dfv.values[int(ny/2),:]

plt.figure(figsize=[3,4],dpi=250)
#plt.title("Flux simulation from the SolTrace engine")
plt.plot(xv,flux_centerlinev, 'k.-', label='SolTrace GUI') #, vmin=240, vmax=420)
plt.plot(x,np.asarray(flux_centerline), 'r.-', label='pysoltrace') #, vmin=240, vmax=420)
#plt.title(f"max flux {flux_st.max():.0f} kW/m2, mean flux {flux_st.mean():.1f}")
plt.title('LS-3 Validation')
plt.xlabel('x [m]')
plt.ylabel('flux at y=0 [kW/m2]')
plt.legend(fontsize=8)
#plt.xlim([-0.5,0.5])
plt.show()
"""
post-processing functions for output from pysoltrace
------------------------------------------------------------
soltrace numbering convention
stages:
    1 = collector
    2 = receiver
elements:
    1 = reflected
    -1 = absorbed (isfinal)
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.ion()
import math


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

def get_power_per_ray(PT,df):
    #PT = self.PT
    sunstats = PT.sunstats
    ppr = PT.dni * (sunstats['xmax']-sunstats['xmin']) * (sunstats['ymax']-sunstats['ymin']) / sunstats['nsunrays'] # power per ray [W]
    return ppr

def calc_intercept_factor(df):
    n_coll_rays = get_number_of_hits(df,1,'reflected') # 'all' equivalent to PT.num_ray_hits
    n_rcvr_rays = get_number_of_hits(df,2,'absorbed') #just absorbed or all rays? 
    intercept_factor = n_rcvr_rays/n_coll_rays # * PT.powerperray cancels out in numerator and denominator
    print('intercept factor = {} = {}/{}'.format(intercept_factor,n_rcvr_rays,n_coll_rays))

def create_xy_mesh_cyl(d,l,nx,ny):
    # assumes trough is located at 0,0,0
    #print(d, nx)
    x = np.linspace(-d/2., d/2., nx)
    y = np.linspace(-l/2., l/2., ny)
    #print('size of x = ',np.size(x))
    dx = x[1]-x[0]
    dy = y[1]-y[0]
    
    Xc,Yc = np.meshgrid(x,y,indexing='ij')
    print('size of Xc = ', np.size(Xc))
    return Xc,Yc,x,y,dx,dy

def create_polar_mesh_cyl(d,l,nx,ny):
    global dx
    global dy
    # creates mesh of unrolled cylinder surface
    # xmin = -0.11215496988813
    # xmax = 0.11215496988813
    # ymin = -5.099999
    # ymax = 5.099999
    # x = np.linspace(xmin, xmax, nx)
    # y = np.linspace(ymin, ymax, ny)
    # dx = x[1]-x[0]
    
    circumf = math.pi*d
    x = np.linspace(-circumf/2.,circumf/2., nx)
    y = np.linspace(-l/2., l/2., ny)
    dx = circumf/nx
    dy = y[1]-y[0]
    X,Y = np.meshgrid(x,y,indexing='ij')
    psi = np.linspace(-180.,180.,nx) # circumferential angle [deg]
    return X,Y,x,y,dx,dy,psi

def convert_xy_polar_coords(d,x,zloc,focal_len,gui_coords=False):
    # position on circumference becomes position in x
    # assumes axis of cylinder is y
    r = d/2.
    
    # assumes x=0 is at top of absorber tube
    # uses s = r*theta and theta = atan(loc_x/loc_z-focal_len)
    if gui_coords:
        #z = zloc - focal_len + d_rec/2. #reset height of tube to z=0
        cpos = r * np.arctan2(x,(zloc-focal_len))
        print('using gui coords')
    
    # assumes x=0 is at bottom of absorber tube
    # uses s = r*theta and theta = atan(loc_x/r-loc_z)
    else:
        z = zloc - focal_len + d/2. #reset height of tube to z=0
        cpos = r * np.arctan2(x,(r-z))
    #print('size of cpos = {}'.format(np.shape(cpos)))
    #print(cpos)
    #print(cpos.values)
    #print('c_pos = r * atan(x/(r-z))) = ')
    #print('{} = {} * atan({}/{})'.format(cpos.values[0],r,x.values[0],r-z.values[0]))
    return cpos

def generate_receiver_dataframe(df,d_rec,focal_len):
    # copied from lines 74+ in https://github.com/NREL/SolarPILOT/blob/develop/deploy/api/test_solarpilot_soltrace.py
    # assumes that receiver is last stage
    #df_rec = df[df.stage==2] # just receiver stage
    df_rec = df[df.stage==df.stage.unique()[-1]] # just receiver stage
    df_rec = df_rec[df_rec.element<0]  #absorbed rays should be negative? - shouldn't all rays be absorbed?
    df_rec['ypos'] = df_rec.loc_y
    df_rec['cpos'] = convert_xy_polar_coords(d_rec,df_rec.loc_x,df_rec.loc_z,focal_len,gui_coords=True)
    #print(df_rec.describe())
    return df_rec

def compute_fluxmap(PTppr,df_rec,d_rec,l_c,nx,ny,plotflag=False):
    Xc,Yc,x,y,dx,dy,psi = create_polar_mesh_cyl(d_rec,l_c,nx,ny)

    flux_st = np.zeros((nx,ny))
    anode = dx*dy
    ppr = PTppr / anode *1e-3

    # count flux at receiver
    for ind,ray in df_rec.iterrows():
        # choose index of closest location to coordinates
        i = np.argmin(np.abs(x - ray.cpos)) # int(np.where(np.abs(x - ray.loc_x) < tol)[0])
        j = np.argmin(np.abs(y - ray.ypos)) # int(np.where(np.abs(y - ray.loc_y) < tol)[0])
        
        flux_st[i,j] += ppr    

    # coeff of variation/dispersion from Zhang et al. 2022
    c_v = np.std(flux_st)/np.mean(flux_st)
    print('coeff of variation = {}'.format(c_v))

    if plotflag==True:
        #%% plot flux line 
        flux_centerline = np.array(flux_st[:,int(ny/2)])
        plt.figure(figsize=[3,4],dpi=250)
        plt.plot(x, flux_centerline, 'k.-') #, vmin=240, vmax=420)
        plt.title(f"max flux {flux_st.max():.0f} kW/m2, mean flux {flux_st.mean():.1f}")
        plt.xlabel('x [m]')
        plt.ylabel('flux at y=0 [kW/m2]')
        plt.savefig('flux-line.png')
        plt.show()

        #%% contour plot
        plt.figure(figsize=[6,4],dpi=250)
        plt.contourf(Xc, Yc, flux_st, levels=15, cmap='viridis') #, vmin=240, vmax=420)
        plt.colorbar()
        plt.title(f"pysoltrace: \n max flux {flux_st.max():.2f} kW/m2, mean flux {flux_st.mean():.2f}")
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
        plt.savefig('flux-map.png')
        plt.show()
    return flux_st
# %%

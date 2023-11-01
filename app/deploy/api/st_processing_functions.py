"""
pre- and post-processing functions for output from pysoltrace
------------------------------------------------------------
soltrace numbering convention
stages:
    1 = collector and receiver for single stage
elements:
    +N = reflected
    -N = absorbed (isfinal)
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.ion()
import math
from pvlib import solarposition, tracking
import glob
import plotly.graph_objects as go
import plotly.io as io
from scipy import interpolate
from os.path import exists
import pickle
# from scipy import stats as st
# from field_postprocessing_functions import assign_Color

io.renderers.default='browser'

def assign_Color(x):
    if x in 'Spring':
       return 'green'
    if x in 'Summer':
       return 'red'
    if x in 'Fall':
       return 'orange'
    else :
       return 'blue'

def set_optics_props(optics_type):
    if optics_type == 'realistic':
        refl_rho = 0.9 # 1. # trough reflectivity
        absr_alpha = 0.96 # 1. # receiver absorptivity
        absr_rho = 1 - absr_alpha #0. # receiver reflectivity
        refl_spec = 0.0 #None
        # tau = 1. # transmittance of glass envelope
    elif optics_type == 'ideal':
        refl_rho = 1. # trough reflectivity
        absr_rho = 0. # receiver reflectivity
        refl_spec = 0.0 # [mrad] 0.2 == default
        absr_alpha = 1. # receiver absorptivity
        # tau = 1. # transmittance of glass envelope
    elif optics_type == 'yang':
        refl_rho = 0.93 # trough reflectivity
        absr_alpha = 0.96 # receiver absorptivity
        absr_rho = 1 - absr_alpha # receiver reflectivity (assumption in soltrace)
        refl_spec = 1.8 # [mrad] 0.2 == default
        # tau = 0.95 # transmittance of glass envelope
    return refl_rho, absr_alpha, absr_rho, refl_spec

def get_trough_angles_py(times, lat, lon, alt):
    # tstart and tend in UTC
    #times = pd.date_range(tstart, tend, freq=inpfreq)

    solpos = solarposition.get_solarposition(times, lat, lon, altitude=alt) #, method='nrel_numba')

    angles = pd.DataFrame()
    angles = sun_elev_to_trough_angles(solpos.apparent_elevation,solpos.azimuth)
    angles = angles.to_frame(name='nom_trough_angle')
    anglesdf = solpos.merge(angles, left_index = True, right_index = True, how='inner')
    anglesdf['projected_sun_angle'] = anglesdf.nom_trough_angle
    anglesdf.nom_trough_angle[anglesdf['apparent_elevation'] < 0] = 120
    return anglesdf

def sun_elev_to_trough_angles(elev_angles, azimuth_angles):
    x, _, z = get_aimpt_from_sunangles(elev_angles, azimuth_angles)
    trough_angle = get_tracker_angle_from_aimpt(x,z)
    return trough_angle

def get_aimpt_from_sunangles(elev_angles, azimuth_angles):
    x = np.cos(np.radians(elev_angles))*np.sin(np.radians(azimuth_angles))
    y = np.cos(np.radians(elev_angles)) * np.cos(np.radians(azimuth_angles))
    z = np.sin(np.radians(elev_angles))
    return x,y,z

def get_tracker_angle_from_aimpt(x,z):
    tracker_angle = np.degrees(np.arctan2(x,z))
    return tracker_angle

def get_aimpt_from_trough_angle(trough_angle):
    x = np.sin(np.radians(trough_angle))
    z = np.cos(np.radians(trough_angle))
    return x,z

def get_aimpt_from_sunangles_pvlib(zenith, azimuth, factor):
    trough_angles = tracking.singleaxis(
        apparent_zenith=zenith,
        apparent_azimuth=azimuth,
        axis_tilt=0,
        axis_azimuth=180, # pointing east = negative
        max_angle=90,
        backtrack=False,  # for true-tracking
        gcr=0.5)  # irrelevant for true-tracking
    # unfinished
    x = 1
    z = 1
    return x,z

#--- Get intersections for stage and element. Note: input stage/element is zero-indexed
def get_intersections(df, stage, elem = 'all'): #, isfinal = False):
    # Note: stage/element numbering in SolTrace output starts from 1
    # stage : integer, [0:PT.stages-1]
    # elem  : string, ['all','absorbed','reflected']
    
    # if element value is negative, that ray was absorbed
    elems_negative = df.element.unique()[df.element.unique()<0]
    
    # if element value is positive, that ray was reflected
    elems_positive = df.element.unique()[df.element.unique()>0]
    
    if elem == 'all': # all
        inds = np.where(df['stage'].values == stage)[0]
    
    elif elem == 'absorbed': # absorbed --> negative elem values only
        # inds = np.where(np.logical_and(df['stage'].values == stage, df['element'].values < 0))[0]
        if elems_negative.size > 0: # rays have been absorbed
            # inds = np.where(np.logical_and(df['stage'].values == stage, df['element'].values == np.min(elems_negative)))[0]
            inds = df.loc[(df['stage'] == stage) & (df['element'] == np.min(elems_negative)) & (df['cos_z'] > 0.)].index.values # not counting rays from the sun that hit the absorber tube first ('cos_z' = -1)
        else: # all rays missed absorber
            inds = np.array([]) # return empty array with len zero (no intersections)
    
    elif elem == 'reflected':# reflected --> positive elem values only
        # inds = np.where(np.logical_and(df['stage'].values == stage, df['element'].values > 0))[0]
        # inds = np.where(np.logical_and(df['stage'].values == stage, df['element'].values == np.min(elems_positive)))[0]
        inds = df.loc[(df['stage'] == stage) & (df['element'] == np.min(elems_positive))].index.values
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
    stages = df.stage.unique()
    
    # needs fixing, not robust for more than 2 stages
    n_coll_rays = get_number_of_hits(df,np.min(stages),'reflected') # 'all' equivalent to PT.num_ray_hits
    n_rcvr_rays = get_number_of_hits(df,np.max(stages),'absorbed') #just absorbed or all rays? 
    
    # single stage
    # n_coll_rays = get_number_of_hits(df,stages[-1],'reflected') # 'all' equivalent to PT.num_ray_hits
    # n_rcvr_rays = get_number_of_hits(df,stages[-1],'absorbed') #just absorbed or all rays? 
    
    intercept_factor = n_rcvr_rays/n_coll_rays # * PT.powerperray cancels out in numerator and denominator
    print('intercept factor = {} = {}/{}'.format(intercept_factor,n_rcvr_rays,n_coll_rays))
    return intercept_factor

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
    # global dx
    # global dy
    # creates mesh of unrolled cylinder surface
    # xmin = -0.11215496988813
    # xmax = 0.11215496988813
    # ymin = -5.099999
    # ymax = 5.099999
    # x = np.linspace(xmin, xmax, nx)
    # y = np.linspace(ymin, ymax, ny)
    # dx = x[1]-x[0]
    
    circumf = math.pi*d
    # assuming x and y marks the center of the cell on the absorber tube
    # x = np.linspace(-circumf/2.,circumf/2., nx)
    # y = np.linspace(-l/2., l/2., ny)
    # dx = circumf/nx
    # dy = y[1]-y[0]
    
    # assuming x and y mark the edges of the cells on the absorber tube
    # and also ensures that all the cell sizes are the same (important for same anode in heat flux calc)
    dx = circumf/nx
    dy = l/ny
    
    if (nx % 2) == 0: #even
       xpositive = np.arange(0.+dx/2, circumf/2., dx)
       xnegative = -xpositive
    else: #odd
       xpositive = np.arange(0., circumf/2., dx)
       xnegative = -xpositive[1:] #without the initial zero to eliminate the -0., 0. repetition in the combined array
       
    if (ny % 2) == 0: #even
       ypositive = np.arange(0.+dy/2, l/2., dy)
       ynegative = -ypositive
    else: #odd
       ypositive = np.arange(0., l/2., dy)
       ynegative = -ypositive[1:] #without the initial zero to eliminate the -0., 0. repetition in the combined array
    
    y = np.sort(np.concatenate([ynegative, ypositive]))
    x = np.sort(np.concatenate([xnegative, xpositive]))
    
    # print('dx = {}, dy = {}'.format(dx,dy))
    # print('x = {}'.format(x))
    # print('y = {}'.format(y))
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
        cpos = r * np.arctan2(x,(zloc-focal_len)) # cpos=0 is top, cpos=circumference/2 is bottom
        # print('using gui coords')
    
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
    print('ppr = {}'.format(ppr))

    # count flux at receiver
    for ind,ray in df_rec.iterrows():
        # i = int(np.where(np.abs(x - ray.loc_x) < tol)[0])
        # j = int(np.where(np.abs(y - ray.loc_y) < tol)[0])
        # choose index of closest location to coordinates
        i = np.argmin(np.abs(x - ray.cpos)) 
        j = np.argmin(np.abs(y - ray.ypos)) 
        
        flux_st[i,j] += ppr    

    # coeff of variation/dispersion from Zhang et al. 2022
    c_v = np.std(flux_st)/np.mean(flux_st)
    print('coeff of variation = {}'.format(c_v))

    if plotflag==True:
        #% plot flux line 
        fig, axs = plt.subplots(1,2,figsize=[6,4],dpi=250)
        
        flux_centerline = np.array(flux_st[:,int(ny/2)])
        # plt.figure(figsize=[3,4],dpi=250)
        axs[0].plot(x, flux_centerline, 'k.-') #, vmin=240, vmax=420)
        # axs[0].set_title(f"max flux {flux_st.max():.2f} kW/m2, mean flux {flux_st.mean():.1f}")
        axs[0].set_xlabel('x [m]')
        axs[0].set_ylabel('flux at y=0 [kW/m2]')
        # plt.savefig('flux-line.png')
        # plt.show()

        #% contour plot
        # plt.figure(figsize=[6,4],dpi=250)
        # cf = axs[1].contourf(Xc, Yc, flux_st, levels=15, cmap='viridis') #, vmin=240, vmax=420)
        cf = axs[1].pcolormesh(Xc, Yc, flux_st, cmap='viridis') #, vmin=240, vmax=420)
        fig.colorbar(cf, ax=axs[1], label='flux [kW/m2]')
        axs[1].set_title(f"pysoltrace: \n max flux {flux_st.max():.2f} kW/m2, mean flux {flux_st.mean():.2f}")
        axs[1].set_xlabel('x [m]')
        axs[1].set_ylabel('y [m]')
        plt.tight_layout()
        plt.savefig('flux-map.png')
        plt.show()
        
    return flux_st, c_v

def compute_fluxmap_solarpilot(PTppr, dfr, d_rec, focal_len, nx, ny):
    # dfr = df[df.stage==df.stage.unique()[-1]] # just receiver stage
    # dfr = dfr[dfr.element<0]  #absorbed rays should be negative? - shouldn't all rays be absorbed?
    
    dfr['ypos'] = dfr.loc_y
    dfr['cpos'] = convert_xy_polar_coords(d_rec,dfr.loc_x,dfr.loc_z,focal_len,gui_coords=True)
    
    Xc,Yc,x,y,dx,dy,psi = create_polar_mesh_cyl(d_rec,focal_len,nx,ny)

    flux_st = np.zeros((ny,nx))
    # dx = d_rec*np.pi / nx 
    # dy = h_rec / ny
    anode = dx*dy
    ppr = PTppr / anode *1e-3

    for ind,ray in dfr.iterrows():

        j = int(ray.cpos/dx)
        i = int(ray.ypos/dy)

        flux_st[i,j] += ppr
    
    return flux_st

def compute_fluxmap_gui(PTppr, dfr, l_c, r, focal_len, nbinsx, nbinsy):
    # adapted to python from SolTrace https://github.com/NREL/SolTrace/blob/47b7abd41d281589700ab6199263df53183086f0/app/src/project.cpp#L1226
    
    # is partical cylinder
    # ZVal = r
    # x = Radius * asin(x / r);
    dfr['ypos'] = dfr.loc_y
    dfr['cpos'] = convert_xy_polar_coords(d_rec,dfr.loc_x,dfr.loc_z,focal_len,gui_coords=True)
        
    # r = d_rec/2.
    
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
    
    zScale = PTppr/(binszx*binszy) # !! CHECK PTppr vs PowerPerRay
    
    SumFlux = 0.
    RayCount = 0
    fluxGrid = np.zeros((nbinsx,nbinsy))
    
    GridIncrementX = -1
    GridIncrementY = -1
    for index,row in dfr.iterrows():
        x = row['cpos']
        y = row['ypos']
        
        while ((minx + (GridIncrementX+1)*binszx) < x):
            GridIncrementX += 1
        while ((miny + (GridIncrementY+1)*binszy) < y):
            GridIncrementY += 1
        fluxGrid[GridIncrementX,GridIncrementY] += 1 #if ray falls inside a bin, increment count for that bin
        RayCount += 1 #increment ray intersection counter
    print('RayCount = ',RayCount)
    
    return fluxGrid
    
    # for i in range(nbinsx):
    #     for j in range(nbinsy):
    #         z = fluxGrid[i,j]*zscale
    #         SumFlux = SumFlux + z

def plot_sun_trough_deviation_angles(fulldata, sensorloc, adj_flag = True):
    fig, axs = plt.subplots(3,1,figsize=[9,7],dpi=250,sharex=True)

    axs[0].plot(fulldata.apparent_elevation,'k.-')
    axs[0].set_ylabel('sun elev. angle [deg]')

    if adj_flag:
        devkey = [col for col in fulldata.filter(regex='Tilt_adjusted').columns if sensorloc in col][0]
    elif sensorloc == 'nominal':
        devkey = 'nom_trough_angle' # 'nominal'
    else:
        devkey = [col for col in fulldata.filter(regex='Tilt$').columns if sensorloc in col][0]
    
    axs[1].plot(fulldata.nom_trough_angle, '.-', label='nominal')
    axs[1].plot(fulldata[devkey], 'k.', label=devkey[0])
    axs[1].set_ylabel('trough_angle')
    axs[1].legend()
    
    # calculate tracking error
    track_error = fulldata[devkey] - fulldata['nom_trough_angle']
    axs[2].plot(track_error, '.-')
    axs[2].set_ylabel('deviation [deg]')
    
    for ax in axs:
        ax.tick_params(labelrotation=30)
    plt.tight_layout()
    plt.show()

def plot_stats_deviation(track_error_stats, critical_angle_error=0.79):
    sensor_locs = track_error_stats.index.unique(level=1).values
    fig,axs = plt.subplots(1,3,sharey=True,figsize=[9,3],dpi=250)
    for ax,sinputloc in zip(axs.ravel(),sensor_locs):
        rows = track_error_stats.index.unique(level=0).values
        ys = track_error_stats.loc[(rows,sinputloc),'absmean']
        stds = track_error_stats.loc[(rows,sinputloc),'absstd']
        ax.plot(rows, ys,'.-', label='$\overline{\epsilon}$')
        ax.fill_between(rows, ys-stds, ys+stds, color='C0', alpha=0.4, label='$\overline{\epsilon} + \sigma$')
        ax.fill_between(rows, ys-2*stds, ys+2*stds, color='C0', alpha=0.2, label='$\overline{\epsilon} + 2\sigma$')
        ax.plot(rows, track_error_stats.loc[(rows,sinputloc),'absmax'], 'k.', label='peak')
        # ax.plot(rows, track_error_stats.loc[(rows,sloc),'absmin'], 'k.', label='')
        ax.axhline(critical_angle_error, color='0.6', label='critical $\epsilon$')
        ax.axhline(0, color='0.5', linestyle=':')
        ax.set_xlabel('row')
        ax.set_title(sinputloc)
    
    axs[-1].legend(bbox_to_anchor=(1, 1.1), loc='upper left', fontsize=10)
    axs[0].set_ylabel('|trough angle deviation| ($\epsilon$) [deg]')
    plt.show()

def plot_stats_intercept_factor(results):
    # rows = inputdata.index.unique(level=0).values
    # sensorlocs = inputdata.index.unique(level=1).values
    # fig,axs = plt.subplots(1,3,sharey=True,figsize=[9,3],dpi=250)
    # for ax,sinputloc in zip(axs.ravel(),sensor_locs):
    #     rows = inputdata.index.unique(level=0).values
    #     ys = inputdata.loc[(rows,sinputloc),'absmean']
    #     stds = inputdata.loc[(rows,sinputloc),'absstd']
    #     ax.plot(rows, ys,'.-', label='$\overline{\epsilon}$')
    #     ax.fill_between(rows, ys-stds, ys+stds, color='C0', alpha=0.4, label='$\overline{\epsilon} + \sigma$')
    #     ax.fill_between(rows, ys-2*stds, ys+2*stds, color='C0', alpha=0.2, label='$\overline{\epsilon} + 2\sigma$')
    #     ax.plot(rows, track_error_stats.loc[(rows,sinputloc),'absmax'], 'k.', label='peak')
    #     # ax.plot(rows, track_error_stats.loc[(rows,sloc),'absmin'], 'k.', label='')
    #     ax.axhline(critical_angle_error, color='0.6', label='critical $\epsilon$')
    #     ax.axhline(0, color='0.5', linestyle=':')
    #     ax.set_xlabel('row')
    #     ax.set_title(sinputloc)
    
    # axs[-1].legend(bbox_to_anchor=(1, 1.1), loc='upper left', fontsize=10)
    # axs[0].set_ylabel('|trough angle deviation| ($\epsilon$) [deg]')
    
    resultsdf = results[list(results.keys())[0]]
    rows = resultsdf.index.unique(level=0).values
    sensorlocs = resultsdf.index.unique(level=1).values
    fig,axs = plt.subplots(1,3,figsize=[9,3],sharey=True,dpi=250)
    for ax,sloc in zip(axs.ravel(),sensorlocs):
        ys = resultsdf.loc[(rows,sloc,'absmean'),'intercept_factor']
        yp2std = resultsdf.loc[(rows,sloc,'absmean+2std'),'intercept_factor']
        ypstd = resultsdf.loc[(rows,sloc,'absmean+std'),'intercept_factor']
        ymstd = resultsdf.loc[(rows,sloc,'absmean-std'),'intercept_factor']
        maxs = resultsdf.loc[(rows,sloc,'absmax'),'intercept_factor']
        # stds = track_error_stats.loc[(rows,sloc),'absstd']
        ax.plot(rows, ys,'.-', label='$\overline{\gamma}$')
        ax.fill_between(rows, yp2std, ymstd, color='C0', alpha=0.2, label='$\overline{\gamma} + 2\sigma$')
        ax.fill_between(rows, ypstd, ymstd, color='C0', alpha=0.4, label='$\overline{\gamma} + \sigma$')
        ax.plot(rows, maxs, 'k.', label='$min(\gamma)$')
        ax.axhline(0,color='0.8',linestyle=':')
        # ax.plot(rows, track_error_stats.loc[(rows,sloc),'absmin'], 'k.', label='')
        # ax.axhline(0, color='0.5', linestyle=':')
        ax.set_xlabel('row')
        ax.set_title(sloc)

    axs[-1].legend(bbox_to_anchor=(1, 1.1), loc='upper left', fontsize=10)
    axs[0].set_ylabel('intercept factor ($\gamma$)')
    plt.show()

def plot_diurnal_cycle_optical_performance(mediandf, resultsdf, critical_angle_error):
    # from "char" mode of running SolTrace
    
    # merge dataframes
    combineddf = mediandf.merge(resultsdf, how='outer', left_index=True, right_index=True)
    # replace intercept factor values
    combineddf.loc[combineddf.trough_angle_dev > 1.5, 'intercept_factor'] = 0
    
    # get nominal value from wintertime nominal cycle
    nomfn = '/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/nominal_3_5_1E+05hits_realistic_optics.p'
    if exists(nomfn):
        tmp = pickle.load(open(nomfn,'rb'))
        nominaldf = tmp[1]['nominal']
    nominaldf['time'] = nominaldf.index.tz_localize('UTC').tz_convert('US/Pacific')
    nominaldf['time'] = nominaldf.time.dt.hour + nominaldf.time.dt.minute/60.

    f = interpolate.interp1d(nominaldf.time, nominaldf.intercept_factor, bounds_error=False, fill_value = nominaldf.intercept_factor.max())
    combineddf.loc[combineddf.trough_angle_dev < critical_angle_error, 'intercept_factor'] = f(combineddf.loc[combineddf.trough_angle_dev < critical_angle_error].index)

    # plotting
    plt.rcParams.update({'font.size': 14})
    markers = ['.','x','+']
    propcolors = ['orange','green','blue']
    colorsdict = {}
    for n,season in enumerate(combineddf.season.unique()):
        colorsdict[season] = propcolors[n]
        
    markersdict = {}
    for n,sensor in enumerate(combineddf.sensor.unique()):
        markersdict[sensor] = markers[n]
        
        
    tilt_col_list = combineddf.sensor.unique()
    
    fig,axs = plt.subplots(2, 1, dpi=250, sharex=True, figsize=[12,5])
    for season, d1 in combineddf.groupby('season'):
        # print(season)
        for sensor, d in d1.groupby('sensor'):
            axs[0].plot(d.trough_angle_dev, linestyle='', color=colorsdict[season], marker = markersdict[sensor]) #, label=column[:6])
            axs[1].plot(d.intercept_factor, linestyle='', color=colorsdict[season], marker = markersdict[sensor]) #, label=column[-6:])
    axs[0].axhline(critical_angle_error, color='k', label='critical angle \n deviation')

    # axs[0].set_ylabel('wind speed at 7m \n [m/s]')      
    # axs[1].set_ylabel('wind dir at 7m \n [deg]')      
    # axs[2].set_ylabel('nom trough angle \n [deg]')
    #axs[2].set_ylabel('{} trough angle \n [deg]'.format(srow))
    axs[0].set_ylabel('abs. val. trough \n angle deviation [deg]')
    axs[1].set_ylabel('intercept factor')

    axs[0].set_ylim([0, 1.4])

    axs[0].plot(np.nan, np.nan, label='spring', color='green')
    axs[0].plot(np.nan, np.nan, label='summer', color='red')
    axs[0].plot(np.nan, np.nan, label='fall', color='orange')
    axs[0].plot(np.nan, np.nan, label='winter', color='blue')

    for n,column in enumerate(tilt_col_list):
        axs[0].plot(np.nan, np.nan, label=column[:6], color = 'k', marker=markers[n], linewidth=0)
        axs[1].plot(np.nan, np.nan, label=column[:6], color = 'k', marker=markers[n], linewidth=0)

    axs[0].legend(bbox_to_anchor=(1, 1.1), loc='upper left', fontsize=9)
    axs[1].legend(bbox_to_anchor=(1, 1.1), loc='upper left', fontsize=9)
    plt.show()

def plot_time_series(solpos, intercept_factor, flux_centerline_time, c_v, x):
    fig, axs = plt.subplots(4,1,figsize=[9,7],dpi=250)

    axs[0].plot(solpos.apparent_elevation,'k.-')
    axs[0].set_ylabel('sun elev. angle [deg]')

    axs[1].plot(solpos.index, intercept_factor, '.-')
    axs[1].set_ylabel('intercept factor')
    
    axs[2].plot(solpos.index, c_v, '.-')
    axs[2].set_ylabel('coeff of variation')
     
    fluxcntr = np.array(flux_centerline_time).T
    cf = axs[3].contourf(solpos.index, x, fluxcntr, levels=100, cmap='turbo')
    axs[3].set_ylabel('x [m]')
    fig.colorbar(cf, ax=axs[3], label='flux at y=0')

    for ax in axs:
        ax.tick_params(labelrotation=30)

    plt.tight_layout()
    plt.show()

def plot_time_series_compare(nominaldf, inputsdf, outputsdf, x, sensorloc):
    #fig, axs = plt.subplots(5,1,figsize=[10,9],dpi=250)
    fig, axs = plt.subplot_mosaic("AE;BE;CF;DF",figsize=[12,7],dpi=250)

    # axs['A'].plot(inputsdf.apparent_elevation,'k.:')
    # axs['A'].set_ylabel('sun elev. angle [deg]')
    # axs['A'].set_title(sensorloc)
    
    axs['A'].plot(inputsdf.nom_trough_angle, 'k.-', label='nominal')
    if sensorloc == 'validation':
        axs['A'].plot(inputsdf.trough_angle, '.-', label='actual')
    else:
        devkey = [col for col in inputsdf.filter(regex='Tilt$').columns if sensorloc in col][0]
        axs['A'].plot(inputsdf[devkey],'.', label=sensorloc)
    axs['A'].set_ylabel('trough angle [deg]')
    
    if sensorloc == 'validation':
        axs['B'].plot(inputsdf.trough_angle_dev, '.-')
    else:
        # devkey = [col for col in inputsdf.filter(regex='trough_angle_dev').columns if sensorloc in col]
        track_error = inputsdf[devkey] - inputsdf['nom_trough_angle']
        axs['B'].plot(track_error,'.-')
    axs['B'].set_ylabel('trough angle \n deviation [deg]')

    if sensorloc == 'validation':
        axs['C'].plot(inputsdf.index, np.ones((len(inputsdf.index))), 'k.-', label='nominal')
        axs['C'].plot(outputsdf.index, outputsdf.intercept_factor, '.-', label=sensorloc)
    else:  
        axs['C'].plot(nominaldf.index, nominaldf.intercept_factor, 'k.-', label='nominal')
        axs['C'].plot(inputsdf.index, outputsdf.intercept_factor, '.-', label=sensorloc)
    axs['C'].set_ylabel('intercept factor')
    axs['C'].set_title('nominal avg = {:2f}, actual avg = {:2f}'.
                     format(nominaldf.intercept_factor.mean(),
                            np.mean(outputsdf.intercept_factor)))
    ymax = np.maximum(1., np.max(outputsdf.intercept_factor))
    axs['C'].set_ylim([0, ymax])
    axs['C'].legend()
    
    if sensorloc == 'validation':
        # axs['D'].plot(inputsdf.index, np.ones((len(inputsdf.index))), 'k.-', label='nominal')
        axs['D'].plot(outputsdf.index, outputsdf.coeff_var, '.-', label=sensorloc)
    else: 
        axs['D'].plot(nominaldf.index, nominaldf.coeff_var, 'k.-', label='nominal')
        axs['D'].plot(inputsdf.index, outputsdf.coeff_var, '.-', label=sensorloc)
    axs['D'].set_ylabel('coeff of variation')
    axs['D'].set_title('nominal avg = {:2f}, actual avg = {:2f}'.
                     format(nominaldf.coeff_var.mean(),
                            np.mean(outputsdf.coeff_var)))
    axs['D'].set_ylim([1, 6])
    axs['D'].legend()
    
    vmin = 0.0
    vmax = np.max(list(outputsdf.flux_centerline.values))
    levels = np.linspace(vmin,vmax,100)
    
    fluxcntr2 = np.stack(nominaldf.flux_centerline.values).T
    cf2 = axs['E'].contourf(nominaldf.index, x, fluxcntr2, 
                            levels=levels, cmap='turbo')
    axs['E'].set_ylabel('x [m]')
    axs['E'].set_title('nominal')
    fig.colorbar(cf2, ax=axs['E'], label='flux at y=0', extend='both')
     
    fluxcntr = np.stack(outputsdf.flux_centerline.values).T
    # print('inner function')
    # print('inputsdf.index = ',inputsdf.index)
    # print('x = ',x)
    # print('fluxcntr = ',fluxcntr)
    cf = axs['F'].contourf(inputsdf.index, x, fluxcntr, levels=levels, 
                           cmap='turbo')
    axs['F'].set_ylabel('x [m]')
    fig.colorbar(cf, ax=axs['F'], label='flux at y=0')
    axs['F'].set_title('actual')

    axs['D'].tick_params(labelrotation=30)
    axs['F'].tick_params(labelrotation=30)

    plt.tight_layout()
    plt.show()

def plot_time_series_compare_nominal(results, x):
    # nom trough angle
    # intercept factor
    # flux map
    resultsdf = results['nominal']
    
    fig, axs = plt.subplots(3,1,figsize=[8,8],dpi=250)
    # fig, axs = plt.subplot_mosaic("AE;BE;CF;DF",figsize=[12,7],dpi=250)

    axs[0].plot(resultsdf.nom_trough_angle,'.-k')
    axs[0].set_ylabel(r'$\beta_{nominal} \; [^\circ]$')
    
    axs[1].plot(resultsdf.intercept_factor, '.-k')
    axs[1].set_ylabel('$\lambda_{nominal} \; [-]$')
    
    vmin = 0.0
    vmax = np.max(list(resultsdf.flux_centerline.values))
    levels = np.linspace(vmin,vmax,100)
    
    fluxcntr2 = np.stack(resultsdf.flux_centerline.values).T
    cf2 = axs[2].contourf(resultsdf.index, x, fluxcntr2, 
                            levels=levels, cmap='turbo')
    axs[2].set_ylabel('x [m]')
    # axs['E'].set_title('nominal')
    fig.colorbar(cf2, ax=axs[2], label='flux at y=0', extend='both')
     
    
    # axs['A'].plot(inputsdf.apparent_elevation,'k.:')
    # axs['A'].set_ylabel('sun elev. angle [deg]')
    # axs['A'].set_title(sensorloc)
    
    # axs['A'].plot(inputsdf.nom_trough_angle, 'k.-', label='nominal')
    # if sensorloc == 'validation':
    #     axs['A'].plot(inputsdf.trough_angle, '.-', label='actual')
    # else:
    #     devkey = [col for col in inputsdf.filter(regex='Tilt$').columns if sensorloc in col][0]
    #     axs['A'].plot(inputsdf[devkey],'.', label=sensorloc)
    # axs['A'].set_ylabel('trough angle [deg]')
    
    # if sensorloc == 'validation':
    #     axs['B'].plot(inputsdf.trough_angle_dev, '.-')
    # else:
    #     # devkey = [col for col in inputsdf.filter(regex='trough_angle_dev').columns if sensorloc in col]
    #     track_error = inputsdf[devkey] - inputsdf['nom_trough_angle']
    #     axs['B'].plot(track_error,'.-')
    # axs['B'].set_ylabel('trough angle \n deviation [deg]')

    # if sensorloc == 'validation':
    #     axs['C'].plot(inputsdf.index, np.ones((len(inputsdf.index))), 'k.-', label='nominal')
    #     axs['C'].plot(outputsdf.index, outputsdf.intercept_factor, '.-', label=sensorloc)
    # else:  
    #     axs['C'].plot(nominaldf.index, nominaldf.intercept_factor, 'k.-', label='nominal')
    #     axs['C'].plot(inputsdf.index, outputsdf.intercept_factor, '.-', label=sensorloc)
    # axs['C'].set_ylabel('intercept factor')
    # axs['C'].set_title('nominal avg = {:2f}, actual avg = {:2f}'.
    #                  format(nominaldf.intercept_factor.mean(),
    #                         np.mean(outputsdf.intercept_factor)))
    # ymax = np.maximum(1., np.max(outputsdf.intercept_factor))
    # axs['C'].set_ylim([0, ymax])
    # axs['C'].legend()
    
    # if sensorloc == 'validation':
    #     # axs['D'].plot(inputsdf.index, np.ones((len(inputsdf.index))), 'k.-', label='nominal')
    #     axs['D'].plot(outputsdf.index, outputsdf.coeff_var, '.-', label=sensorloc)
    # else: 
    #     axs['D'].plot(nominaldf.index, nominaldf.coeff_var, 'k.-', label='nominal')
    #     axs['D'].plot(inputsdf.index, outputsdf.coeff_var, '.-', label=sensorloc)
    # axs['D'].set_ylabel('coeff of variation')
    # axs['D'].set_title('nominal avg = {:2f}, actual avg = {:2f}'.
    #                  format(nominaldf.coeff_var.mean(),
    #                         np.mean(outputsdf.coeff_var)))
    # axs['D'].set_ylim([1, 6])
    # axs['D'].legend()
    
    # vmin = 0.0
    # vmax = np.max(list(outputsdf.flux_centerline.values))
    # levels = np.linspace(vmin,vmax,100)
    
    # fluxcntr2 = np.stack(nominaldf.flux_centerline.values).T
    # cf2 = axs['E'].contourf(nominaldf.index, x, fluxcntr2, 
    #                         levels=levels, cmap='turbo')
    # axs['E'].set_ylabel('x [m]')
    # axs['E'].set_title('nominal')
    # fig.colorbar(cf2, ax=axs['E'], label='flux at y=0', extend='both')
     
    # fluxcntr = np.stack(outputsdf.flux_centerline.values).T
    # cf = axs['F'].contourf(inputsdf.index, x, fluxcntr, levels=levels, 
    #                        cmap='turbo')
    # axs['F'].set_ylabel('x [m]')
    # fig.colorbar(cf, ax=axs['F'], label='flux at y=0')
    # axs['F'].set_title('actual')

    # axs['D'].tick_params(labelrotation=30)
    # axs['F'].tick_params(labelrotation=30)

    plt.tight_layout()
    plt.show()


def transform_stage_to_global_coords(df, PT, st):
    locs_stage = np.array([df[k].values for k in ['loc_x','loc_y','loc_z']])
    #cos_stage = np.array([df[k].values for k in ['cos_x','cos_y','cos_z']])
    target = st #st
    euler = PT.util_calc_euler_angles([target.position.x, target.position.y, target.position.z], 
                                      [target.aim.x, target.aim.y, target.aim.z], target.zrot)
    T = PT.util_calc_transforms(euler)
    global_origin = np.array([0, 0, 0]).reshape((3,1))
    locs = np.matmul(T['rloctoref'][0:-1], locs_stage-global_origin)
    #locs_transform = PT.util_transform_to_ref(locs_stage, cos_stage, global_origin, T['rloctoref'])
    
    # create rotated dataframe
    col_list = ['loc_x', 'loc_y', 'loc_z', 'element', 'stage', 'number']
    dfr = df[col_list].copy()
    
    # replace coords with rotated coords
    dfr.loc[:,'loc_x'] = locs[0,:]
    dfr.loc[:,'loc_y'] = locs[1,:]
    dfr.loc[:,'loc_z'] = locs[2,:]
    return dfr

def plot_rays_globalcoords(df, PT, st):
    # Plotting with plotly
    dfr = transform_stage_to_global_coords(df, PT, st)
    
    # fig = go.Figure(data=go.Scatter3d(x=locs_stage[0], y=locs_stage[1], z=locs_stage[2], mode='markers', marker=dict( size=1, color='red', opacity=0.8, ) ))
    # fig.add_trace(go.Scatter3d(x=locs_global[0], y=locs_global[1], z=locs_global[2], mode='markers', marker=dict( size=1, color='black', opacity=0.8, ) ))

    # stage coords
    fig = go.Figure(data=go.Scatter3d(x=df.loc_x.values, y=df.loc_y.values, z=df.loc_z.values, mode='markers', marker=dict( size=1, color='red', opacity=0.1, ) ))
    # plot rays in stage coord sys
    for i in range(50,100):
        dfs = df[df.number == i]
        
        ray_x = dfs.loc_x 
        ray_y = dfs.loc_y
        ray_z = dfs.loc_z
        raynum = dfs.number
        fig.add_trace(go.Scatter3d(x=ray_x, y=ray_y, z=ray_z, mode='lines', line=dict(color='red', width=0.5)))
    
    # global coords
    fig.add_trace(go.Scatter3d(x=dfr.loc_x.values, y=dfr.loc_y.values, z=dfr.loc_z.values, mode='markers', marker=dict( size=1, color='black', opacity=0.8, ) ))
    # plot rays in global
    for i in range(50,100):
        dfs = dfr[dfr.number == i]
        
        ray_x = dfs.loc_x 
        ray_y = dfs.loc_y
        ray_z = dfs.loc_z
        raynum = dfs.number
        fig.add_trace(go.Scatter3d(x=ray_x, y=ray_y, z=ray_z, mode='lines', line=dict(color='black', width=0.5)))
    
    #fig.add_trace(go.Scatter3d(x=locs_transform[0], y=locs_transform[1], z=locs_transform[2], mode='markers', marker=dict( size=1, color='blue', opacity=0.8, ) ))
    fig.update_layout(showlegend=False)
    fig.show()
    
def plot_sun_position(solpos):
    #% 3d plot of sun vectors
    # origin = np.zeros((3,len(solpos['sun_pos_x'])))
    # xs = np.column_stack((origin[0,:], solpos['sun_pos_x']))
    # print(xs)
    # ys = np.column_stack((origin[1,:], solpos['sun_pos_y']))
    # zs = np.column_stack((origin[2,:], solpos['sun_pos_z']))
    # fig = go.Figure(go.Scatter3d(x=xs, y=ys, z=zs, mode='markers')) #,marker_color=date_to_val))
    fig = go.Figure(go.Scatter3d(x=solpos['sun_pos_x'], y=solpos['sun_pos_y'], 
                                 z=solpos['sun_pos_z'], mode='markers')) #,marker_color=date_to_val))
    # xaims, zaims = get_aimpt_from_sunangles(solpos.apparent_elevation, solpos.azimuth)
    # print(xaims)
    #yaims = np.zeros((np.shape(xaims)))
    #fig.add_trace(go.Scatter3d(x=xaims, y=yaims, z=zaims, mode='markers'))
    fig.update_layout(showlegend=False)
    fig.show()
    
def load_field_data(path, year, month, day, fileres, outres):
    inflow_files = sorted(glob.glob(path +'Met_masts_v0/Inflow_Mast_' + fileres + '_' + year + '-' + month + '-' + day + '_' + '*.pkl'))   #
    loads_files = sorted(glob.glob(path +'Loads_v1/Loads_' + fileres + '_' + year + '-' + month + '-' + day + '_' + '*.pkl'))

    inflow = pd.DataFrame()
    for datafile in inflow_files:
        #print(datafile)
        inflow = pd.concat( [inflow, pd.read_pickle(datafile)])
    #drop duplicates
    # inflow = inflow[inflow.index.drop_duplicates(keep='first')]
    # print(inflow)

    loads = pd.DataFrame()
    for datafile in loads_files:
        #print(datafile)
        tmpdf = pd.read_pickle(datafile)
        for col in tmpdf.columns:
            if 'R1_SO_tilt' in col:
                print('correcting {} for tilt vs Tilt in {}'.format(col,datafile))
                if tmpdf[col].isna().any()==False: # if the col contains no nans
                    # then just rename it
                    tmpdf = tmpdf.rename(columns={'R1_SO_tilt':'R1_SO_Tilt'})
                else:
                    print('need code to handle when there are nans: combine with other column')
        loads = pd.concat( [loads, tmpdf]) 
        # print(loads.keys())

    delta = loads.index[0]-inflow.index[0]
    if delta.microseconds > 0:
        print('rounding index to nearest full second')
        # round timestamp to nearest full second
        loads.index = loads.index.floor('T') # or should this be ceiling?
        #loads.index = loads.index.round('1S')

    #% merge into one dataframe
    # complete dataframe with inflow and masts and trough angles
    fulldata = inflow.merge(loads, left_index = True, right_index=True, how="inner") 
    # fulldata = fulldata.resample(outres).asfreq() #1 minute
    
    return fulldata

def plot_time_series_compare_sensors(nominaldf, inputsdf, results, x, sensorlocs):
    fig, axs = plt.subplot_mosaic("AE;BF;CG;DH",sharex=True,figsize=[12,7],dpi=250)

    axs['A'].plot(inputsdf.wspd_3m,'k.:')
    axs['A'].set_ylabel('wind speed \n [m/s]')
    
    axs['B'].plot(inputsdf.wdir_3m,'k.:')
    axs['B'].axhspan(225, 315, facecolor='0.8', alpha=0.9)
    axs['B'].set_ylabel('wind direction \n [deg]')
    # axs['B'].plot(inputsdf.nom_trough_angle, 'k-', label='nominal')
    # for sensorloc in sensorlocs:
    #     devkey = [col for col in inputsdf.filter(regex='Tilt').columns if sensorloc in col]
    #     axs['B'].plot(inputsdf[devkey],'.', label=sensorloc)
    # axs['B'].set_ylabel('trough angle \n [deg]')

    for sensorloc in sensorlocs:
        # devkey = [col for col in inputsdf.filter(regex='trough_angle_dev').columns if sensorloc in col]
        devkey = [col for col in inputsdf.filter(regex='Tilt$').columns if sensorloc in col]
        
        track_error = inputsdf[devkey[0]] - inputsdf['nom_trough_angle']
        axs['C'].plot(abs(track_error),'.-', label=sensorloc)
    axs['C'].set_ylabel('abs. val. trough \n angle deviation [deg]')
    # axs['C'].set_ylim([0, 1])

    axs['D'].plot(nominaldf.index, nominaldf.intercept_factor, 'k-', label='nominal')
    for sensorloc in sensorlocs:
        outputsdf = results[sensorloc]
        axs['D'].plot(inputsdf.index, outputsdf.intercept_factor, '.-', label=sensorloc)
    axs['D'].set_ylabel('intercept \n factor')
    # axs['D'].set_ylim([0.6, 1])
    axs['D'].legend(fontsize=7,loc='upper left')

    # axs['D'].plot(nominaldf.index, nominaldf.coeff_var, 'k.-', label='nominal')
    # for sensorloc in sensorlocs:
    #     outputsdf = results[sensorloc]
    #     axs['D'].plot(inputsdf.index, outputsdf.coeff_var, '.-', label=sensorloc)
    # axs['D'].set_ylabel('coeff of variation')
    # axs['D'].set_ylim([1, 6])
    # axs['D'].legend()

    vmin = 0.0
    vmax = 0.0 # just initializing
    for sensorloc in sensorlocs:
        outputsdf = results[sensorloc]
        tmpmax = np.max(list(outputsdf.flux_centerline.values))
        if tmpmax > vmax:
            vmax = tmpmax
    levels = np.linspace(vmin,vmax,100)

    fluxcntr2 = np.stack(nominaldf.flux_centerline.values).T
    cf2 = axs['E'].contourf(nominaldf.index, x, fluxcntr2, 
                            levels=levels, cmap='turbo')
    axs['E'].set_ylabel('x [m]')
    axs['E'].set_title('nominal')
    fig.colorbar(cf2, ax=axs['E'], label='flux at y=0', extend='both')

    cntraxs = [axs['F'],axs['G'],axs['H']]
    for n,sensorloc in enumerate(sensorlocs[0:3]):
        outputsdf = results[sensorloc]
        ax = cntraxs[n]
        fluxcntr = np.stack(outputsdf.flux_centerline.values).T
        # print('outer function')
        # print('inputsdf.index = ',inputsdf.index)
        # print('x = ',x)
        # print('fluxcntr = ',fluxcntr)
        cf = ax.contourf(inputsdf.index, x, fluxcntr, levels=levels, 
                                cmap='turbo')
        ax.set_xlim([inputsdf.index[0], inputsdf.index[-1]])
        ax.set_ylabel('x [m]')
        fig.colorbar(cf, ax=ax, label='flux at y=0')
        ax.set_title(sensorloc)

    axs['D'].tick_params(axis='x',labelrotation=30)
    axs['H'].tick_params(axis='x',labelrotation=30)

    plt.tight_layout()
    plt.show()

def plot_time_series_optical_results(results, nominaldf = None, critical_angle_error_min = 0.79, critical_angle_error_max = 1.3542636710518383):
    # extract row numbers
    rows = []
    for key in list(results.keys()):
        split_results = key.split('_',1)
        if split_results[0] not in rows:
            rows.append(split_results[0])
    
    fig, axs = plt.subplots(3*len(rows),1,figsize=[8,4*len(rows)],sharex=True,dpi=250)
    
    # n = 0, row = 1, ax = 0
    # n = 1, row = 2, ax = 3
    # n = 2, row = 4, ax = 6
    # nominal
    if nominaldf is not None:
        for n,row in enumerate(rows):
            axs[n*3].plot(nominaldf.intercept_factor*0., '-k', label='nominal') #zero tracking angle
            axs[n*3+1].plot(nominaldf.intercept_factor, '-k') #zero tracking angle
            axs[n*3+2].plot(nominaldf.coeff_var, '-k') #zero tracking angle
    
    
    for n,row in enumerate(rows):
        for column in [column for column in results.keys() if row in column]:
            tracking_error = abs(results[column]['{}_Tilt'.format(column)] - results[column].nom_trough_angle)
            
            # replace values with nans where non-operational
            tracking_error = tracking_error.where(tracking_error < critical_angle_error_max)
            # results[column] = results[column].where(tracking_error < critical_angle_error_max)
            
            axs[n*3].axhline(critical_angle_error_min, color='0.6')
            axs[n*3].plot(tracking_error, '.-', label=column)
            axs[n*3].set_ylim([-0.05,critical_angle_error_max])
            axs[n*3+1].plot(results[column].intercept_factor, '.-', label=column)
            axs[n*3+1].set_ylim([0,1])
            axs[n*3+2].plot(results[column].coeff_var, '.-', label=column)
            axs[n*3+2].set_ylim([1,2])
            # axs[0].plot(data[column],'.', markersize=markersize, label='') #, label=column[:6])
        axs[n*3].set_ylabel(r'$ |\epsilon| \; [^\circ]$')
        axs[n*3].legend(bbox_to_anchor=(1, 1.1), loc='upper left', fontsize=9)
    
        axs[n*3+1].set_ylabel(r'$ \lambda \; [-]$')   
        axs[n*3+2].set_ylabel(r'$ C_v \; [-]$')   

    plt.tight_layout()
    plt.show()
    
def plot_time_series_fluxmap_results(results, x, nominaldf = None):
    
    fig, axs = plt.subplots(1+len(results.keys()),1,figsize=[8,1.5+1.5*len(results.keys())],dpi=250)
    
    #
    # n = 0, row = 1, ax = 0
    # n = 1, row = 2, ax = 3
    # n = 2, row = 4, ax = 6
    
    xmin = np.min([results[list(results.keys())[0]].index[0],nominaldf.index[0]])
    xmax = np.max([results[list(results.keys())[0]].index[-1],nominaldf.index[-1]])
    
    # find vmax
    vmax = 0.
    for column in results.keys():
        tmpmax = np.max(list(results[column].flux_centerline.values))
        if tmpmax > vmax:
            vmax = tmpmax
    
    vmin = 0.0
    levels = np.linspace(vmin,vmax,100)
    
    # nominal
    if nominaldf is not None:
        fluxcntr = np.stack(nominaldf.flux_centerline.values).T
        cf = axs[0].contourf(nominaldf.index, x, fluxcntr, levels=levels, 
                               cmap='turbo')
        axs[0].set_ylabel('x [m]')
        axs[0].set_title('nominal')
    
    for ax,column in zip(axs[1:].ravel(),results.keys()):
            fluxcntr = np.stack(results[column].flux_centerline.values).T
            cf = ax.contourf(results[column].index, x, fluxcntr, levels=levels, 
                                   cmap='turbo')
            ax.set_ylabel('x [m]')
            ax.set_title(column)
            
    for ax in axs:
        fig.colorbar(cf, ax=ax) #, label='flux at y=0')
        ax.set_xlim([xmin, xmax])
        
    for ax in axs[:-1]:
        ax.xaxis.set_major_formatter(plt.NullFormatter())

    axs[-1].tick_params(labelrotation=30)
            
        #     axs[n*3].plot(tracking_error, '.-', label=column)
        #     axs[n*3+1].plot(results[column].intercept_factor, '.-', label=column)
        #     axs[n*3+2].plot(results[column].coeff_var, '.-', label=column)
        #     # axs[0].plot(data[column],'.', markersize=markersize, label='') #, label=column[:6])
        # axs[n*3].set_ylabel(r'$ |\epsilon| \; [^\circ]$')
        # axs[n*3].legend(bbox_to_anchor=(1, 1.1), loc='upper left', fontsize=9)
    
        # axs[n*3+1].set_ylabel(r'$ \lambda \; [-]$')   
        # axs[n*3+2].set_ylabel('coeff of var. [-]')   

    plt.tight_layout()
    plt.show()
    
def assign_marker(x):
    if 'DO' in x:
        marker = '.'
    elif 'Mid' in x:
        marker = 'x'
    else:
        marker = '+'
    return marker

def plot_time_series_median_optical(results, abs_val=False, critical_angle_error_min = 0.79, critical_angle_error_max = 1.3542636710518383):
    lw = 0.75
    ms = 8
    # combineddf = resultsdf.merge(mediandf, left_index = True, right_index = True, how='outer')
    # combineddf['intercept_factor'] = combineddf['intercept_factor'].fillna(0)
    # srows = ['R1', 'R2','R4']
    # markers = ['.','x','+']
    
    rows = []
    for key in list(results.keys()):
        split_results = key.split('_',1)
        if split_results[0] not in rows:
            rows.append(split_results[0])
            
    # find min and max vals for c_v
    ymin = 1.4
    ymax = 1.4 # just initializing
    for nk, key in enumerate(results.keys()):
        outputsdf = results[key]
        tmpmin = np.min(list(outputsdf.coeff_var.values))
        if tmpmin > ymin:
            ymin = tmpmin
        tmpmax = np.max(list(outputsdf.coeff_var.values))
        if tmpmax > ymax:
            ymax = tmpmax
    
    fig, axs = plt.subplots(3*len(rows),1,figsize=[8,4*len(rows)],sharex=True,dpi=250)

    for nk, key in enumerate(results.keys()):
        combineddf = results[key]
        
        marker = assign_marker(key)
        # print(key)
        # print(marker)
        # tilt_col_list = [key] # [col for col in combineddf.filter(regex='Tilt').columns if srow in col]


        tstart = combineddf.index[0] # '2022-12-17 07:59:06.500000' # '2022-12-16 00:00:00.000000-08:00' 
        tend = combineddf.index[-1]
        mediangroupbydf = combineddf.groupby(['timeofday','season']).median(numeric_only=True)

        
        # fig,axs = plt.subplots(4, 1, dpi=250, sharex=True, figsize=[10,10])
        # axs[0].set_title('median diurnal cycle from {} to {}'.format(tstart,tend))
        
        # for ax in axs:
        #     ax.axhline(0, color='k', label='')
        
        # for ax in axs[1:]:
        #     if abs_val:
        #         ax.set_ylabel(r'$|\epsilon| \; [^\circ]$')
        #         ax.axhline(critical_angle_error, color='0.6', label=r'$critical \; \epsilon$')
        #     else:
        #         ax.set_ylabel(r'$\epsilon \; [^\circ]$')
        #         ax.axhline(critical_angle_error, color='0.6', label=r'$critical \; \epsilon$')
        #         ax.axhline(-critical_angle_error, color='0.6', label='')

        for s, dfs in mediangroupbydf.groupby('season'):
            ind = dfs.index.get_level_values(0)
            color = assign_Color(s)
                    
            # axs[0].plot(ind, dfs.nom_trough_angle,"o", fillstyle='none', markeredgecolor='k')
            
            for i,srow in enumerate(rows):
                tilt_col_list = [col for col in combineddf.filter(regex='Tilt$').columns if srow in col]
                for n,tc in enumerate(tilt_col_list):
                    # axs[0].plot(ind, dfs[tc],  linestyle='-', linewidth=0.5, color=color, marker = markers[n]) #, label=column[:6])
                
                    if abs_val:
                        calc_trough_dev = abs(dfs[tc] - dfs.nom_trough_angle)
                    else:
                        calc_trough_dev = dfs[tc] - dfs.nom_trough_angle
                    axs[3*i].set_title(srow)
                    axs[3*i].axhline(critical_angle_error_min, color='0.6')
                    axs[3*i].plot(ind, calc_trough_dev,  linestyle='-', linewidth=lw, color=color, marker = marker, markersize = ms) #, label=column[-6:])
                    axs[3*i].set_ylabel(r'$ |\epsilon| \; [^\circ]$')
                    
                    
                    axs[1+3*i].plot(ind, dfs.intercept_factor,  linestyle='-', linewidth=lw, color=color, marker = marker, markersize = ms) #, label=column[-6:])
                    axs[1+3*i].set_ylabel(r'$ \lambda \; [-]$')   
                    axs[1+3*i].set_ylim([0, 1])
                    
                    axs[2+3*i].plot(ind, dfs.coeff_var,  linestyle='-', linewidth=lw, color=color, marker = marker, markersize = ms) #, label=column[-6:])
                    axs[2+3*i].set_ylabel(r'$ C_v \; [-]$')   
                    axs[2+3*i].set_ylim([ymin, ymax])
                
                if s == 'Fall' and nk == 0:
                    axs[3*i].plot(np.nan, np.nan, label='spring', color='green')
                    axs[3*i].plot(np.nan, np.nan, label='summer', color='red')
                    axs[3*i].plot(np.nan, np.nan, label='fall', color='orange')
                    axs[3*i].plot(np.nan, np.nan, label='winter', color='blue')
                    
                if s == 'Fall' and srow in key:
                    # print(tc)
                    # print(srow)
                    # print(key)
                    axs[1+3*i].plot(np.nan, np.nan, label=tc[:6], color = 'k', marker=marker, linewidth=0)
                    axs[2+3*i].plot(np.nan, np.nan, label=tc[:6], color = 'k', marker=marker, linewidth=0)
                


    
    for ax in axs:
        ax.legend(bbox_to_anchor=(1, 1.1), loc='upper left', fontsize=9)
    axs[-1].set_xlabel('hour of the day [UTC]')
    
    plt.tight_layout()
    
def plot_validation_intercept_factor(results):
    valdata = {}
    filedir = '/Users/bstanisl/Documents/seto-csp-project/SolTrace/firstOptic_Tracking.dat'
    valdata['firstoptic'] = pd.read_csv(filedir, header=None, sep=' ', names=['tracker_error', 'intercept_factor'])
    valdata['firstoptic']['tracker_error'] = np.degrees(valdata['firstoptic']['tracker_error'])
    
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['mathtext.fontset'] = 'dejavuserif'
    plt.rcParams['font.size'] = 14

    resolution_value = 300
    fulldata = results['validation']

    critical_angle_error_min = 0.79 #[deg] from firstoptic validation dataset
    critical_angle_error_max = np.degrees(2.363636e-02) #[deg] from firstoptic validation dataset

    fig = plt.figure(dpi=resolution_value)
    # plt.scatter(valdata['yang-intc']['tracker_error'],valdata['yang-intc']['intercept_factor'],color='k',marker='x', label='Yang et al. 2022')
    plt.plot(valdata['firstoptic']['tracker_error'],valdata['firstoptic']['intercept_factor'],color='k', label='FirstOPTIC')

    sensorloc='validation'
    plt.plot(fulldata['trough_angle_dev'],results[sensorloc].intercept_factor, 'kx', label = 'pysoltrace')
    # plt.scatter(fulldata[devkey],results[sensorloc].intercept_factor, color='r', label = sensorloc)
    # plt.scatter(fulldata[devkey],results[sensorloc].intercept_factor, label = sensorloc)
    plt.axvspan(critical_angle_error_min, critical_angle_error_max, color='0.7', alpha=0.6)
    plt.xlabel('tracking error ($\epsilon$) [$^\circ$]')
    plt.ylabel('intercept factor ($\gamma$) [-]')
    plt.legend()
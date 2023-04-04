"""
post-processing functions for output from pysoltrace
------------------------------------------------------------
soltrace numbering convention
stages:
    1 = collector
    2 = receiver
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
io.renderers.default='browser'

def get_trough_angles():
    fn = '/Users/bstanisl/Documents/seto-csp-project/NSO-field-data/NREL_NSO_meas/trough_angles/sun_angles.txt'
    angles = pd.read_csv(fn, parse_dates={'UTC': [0, 1]}).set_index('UTC')  # ,nrows=200
    angles.iloc[:,-1] = angles.iloc[:,-1].where(angles.iloc[:,-1]>0)
    angles['trough_angle'] = np.degrees( np.arctan2(np.sin(np.radians(angles.iloc[:,-1])), 
                                                    np.sin(np.radians(angles.iloc[:,-2])) ))
    angles.trough_angle = angles.trough_angle.where(angles.trough_angle.isnull()==False, -30)
    angles = -angles + 90
    return angles.trough_angle

def get_sun_angles():
    fn = '/Users/bstanisl/Documents/seto-csp-project/NSO-field-data/NREL_NSO_meas/trough_angles/sun_angles.txt'
    angles = pd.read_csv(fn, parse_dates={'UTC': [0, 1]}).set_index('UTC')  # ,nrows=200
    #angles.iloc[:,-1] = angles.iloc[:,-1].where(angles.iloc[:,-1]>0)
    return angles.iloc[:,-1]

def sun_elev_to_trough_angles(elev_angles, azimuth_angles):
    # trough_angles = np.degrees( np.arctan2(np.sin(np.radians(elev_angles)), np.sin(np.radians(azimuth_angles)) ))
    # print('trough angle = {:2f}'.format(trough_angles))
    # # print(trough_angles)
    # # trough_angles = trough_angles.where(trough_angles.isnull()==False, -30)
    # # print(trough_angles.where(trough_angles.isnull()==False, -30))
    # trough_angles = -trough_angles + 90
    # print('trough angle = {:2f}'.format(trough_angles))
    x, z = get_aimpt_from_sunangles(elev_angles, azimuth_angles)
    trough_angle = get_tracker_angle_from_aimpt(x,z)
    return trough_angle

def get_aimpt_from_sunangles(elev_angles, azimuth_angles):
    # trough_angles = sun_elev_to_trough_angles(elev_angles, azimuth_angles)
    # print('elev angle = {:2f}'.format(elev_angles))
    # print('azimuth angle = {:2f}'.format(azimuth_angles))
    # #print('trough angle = {:2f}'.format(trough_angles))
    # signed_elev_angles = 90 - trough_angles
    # x = factor * np.cos(np.radians(signed_elev_angles))
    # z = x * np.tan(np.radians(signed_elev_angles))
    x = np.cos(np.radians(elev_angles))*np.sin(np.radians(azimuth_angles))
    z = np.sin(np.radians(elev_angles))
    return x,z

def get_tracker_angle_from_aimpt(x,z):
    tracker_angle = np.degrees(np.arctan2(x,z))
    return tracker_angle

def get_aimpt_from_trough_angle(trough_angle):
    x = np.sin(np.radians(trough_angle))
    z = np.cos(np.radians(trough_angle))
    return(x,z)

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
    stages = df.stage.unique()
    
    # needs fixing, [0] and [-1] not robust for more than 2 stages
    n_coll_rays = get_number_of_hits(df,stages[0],'reflected') # 'all' equivalent to PT.num_ray_hits
    n_rcvr_rays = get_number_of_hits(df,stages[-1],'absorbed') #just absorbed or all rays? 
    
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
        #% plot flux line 
        flux_centerline = np.array(flux_st[:,int(ny/2)])
        plt.figure(figsize=[3,4],dpi=250)
        plt.plot(x, flux_centerline, 'k.-') #, vmin=240, vmax=420)
        plt.title(f"max flux {flux_st.max():.0f} kW/m2, mean flux {flux_st.mean():.1f}")
        plt.xlabel('x [m]')
        plt.ylabel('flux at y=0 [kW/m2]')
        plt.savefig('flux-line.png')
        plt.show()

        #% contour plot
        plt.figure(figsize=[6,4],dpi=250)
        plt.contourf(Xc, Yc, flux_st, levels=15, cmap='viridis') #, vmin=240, vmax=420)
        plt.colorbar()
        plt.title(f"pysoltrace: \n max flux {flux_st.max():.2f} kW/m2, mean flux {flux_st.mean():.2f}")
        plt.xlabel('x [m]')
        plt.ylabel('y [m]')
        plt.savefig('flux-map.png')
        plt.show()
    return flux_st, c_v

def plot_sun_trough_deviation_angles(fulldata, sensorloc):
    fig, axs = plt.subplots(3,1,figsize=[9,7],dpi=250,sharex=True)

    axs[0].plot(fulldata.apparent_elevation,'k.-')
    axs[0].set_ylabel('sun elev. angle [deg]')

    devkey = [col for col in fulldata.filter(regex='Tilt').columns if sensorloc in col]
    axs[1].plot(fulldata.nom_trough_angle, '.-', label='nominal')
    axs[1].plot(fulldata[devkey], 'k.', label=devkey[0])
    axs[1].set_ylabel('trough_angle')
    axs[1].legend()

    devkey = [col for col in fulldata.filter(regex='trough_angle_dev').columns if sensorloc in col]
    axs[2].plot(fulldata[devkey], '.-')
    axs[2].set_ylabel('deviation [deg]')
    
    for ax in axs:
        ax.tick_params(labelrotation=30)
    plt.tight_layout()

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

def plot_time_series_compare(nominaldf, inputsdf, outputsdf, x, sensorloc):
    #fig, axs = plt.subplots(5,1,figsize=[10,9],dpi=250)
    fig, axs = plt.subplot_mosaic("AE;BE;CF;DF",sharex=True,figsize=[12,7],dpi=250)

    axs['A'].plot(inputsdf.apparent_elevation,'k.:')
    axs['A'].set_ylabel('sun elev. angle [deg]')
    axs['A'].set_title(sensorloc)
    
    devkey = [col for col in inputsdf.filter(regex='trough_angle_dev').columns if sensorloc in col]
    axs['B'].plot(inputsdf[devkey],'r.-')
    axs['B'].set_ylabel('trough angle \n deviation [deg]')

    axs['C'].plot(nominaldf.index, nominaldf.intercept_factor, 'k.-', label='nominal')
    axs['C'].plot(inputsdf.index, outputsdf.intercept_factor, 'r.-', label=sensorloc)
    axs['C'].set_ylabel('intercept factor')
    axs['C'].set_title('nominal avg = {:2f}, actual avg = {:2f}'.
                     format(nominaldf.intercept_factor.mean(),
                            np.mean(outputsdf.intercept_factor)))
    axs['C'].set_ylim([0, 1])
    
    axs['D'].plot(nominaldf.index, nominaldf.coeff_var, 'k.-', label='nominal')
    axs['D'].plot(inputsdf.index, outputsdf.coeff_var, 'r.-', label=sensorloc)
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
    cf = axs['F'].contourf(inputsdf.index, x, fluxcntr, levels=levels, 
                           cmap='turbo')
    axs['F'].set_ylabel('x [m]')
    fig.colorbar(cf, ax=axs['F'], label='flux at y=0')
    axs['F'].set_title('actual')

    axs['D'].tick_params(labelrotation=30)
    axs['F'].tick_params(labelrotation=30)

    plt.tight_layout()

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
    
def load_field_data(path, year, month, day, fileres, outres):
    inflow_files = sorted(glob.glob(path +'Inflow_Mast_' + fileres + '_' + year + '-' + month + '-' + day + '_' + '*.pkl'))   #
    loads_files = sorted(glob.glob(path +'Loads_' + fileres + '_' + year + '-' + month + '-' + day + '_' + '*.pkl'))

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

    axs['A'].plot(inputsdf.wspd_7m,'k.:')
    axs['A'].set_ylabel('wind speed [m/s]')
    
    axs['B'].plot(inputsdf.nom_trough_angle, 'k-', label='nominal')
    for sensorloc in sensorlocs:
        devkey = [col for col in inputsdf.filter(regex='Tilt').columns if sensorloc in col]
        axs['B'].plot(inputsdf[devkey],'.', label=sensorloc)
    axs['B'].set_ylabel('trough angle [deg]')

    for sensorloc in sensorlocs:
        devkey = [col for col in inputsdf.filter(regex='trough_angle_dev').columns if sensorloc in col]
        axs['C'].plot(inputsdf[devkey],'.-', label=sensorloc)
    axs['C'].set_ylabel('trough angle \n deviation [deg]')

    axs['D'].plot(nominaldf.index, nominaldf.intercept_factor, 'k.-', label='nominal')
    for sensorloc in sensorlocs:
        outputsdf = results[sensorloc]
        axs['D'].plot(inputsdf.index, outputsdf.intercept_factor, '.-', label=sensorloc)
    axs['D'].set_ylabel('intercept factor')
    axs['D'].set_ylim([0, 1])
    axs['D'].legend()

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
        cf = ax.contourf(inputsdf.index, x, fluxcntr, levels=levels, 
                                cmap='turbo')
        ax.set_ylabel('x [m]')
        fig.colorbar(cf, ax=ax, label='flux at y=0')
        ax.set_title(sensorloc)

    axs['D'].tick_params(axis='x',labelrotation=30)
    axs['H'].tick_params(axis='x',labelrotation=30)

    plt.tight_layout()

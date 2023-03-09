import os
os.chdir('/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/')
from pysoltrace import PySolTrace, Point
import random
import pandas as pd
import copy
import matplotlib.pyplot as plt
plt.ion()
import numpy as np

from postprocessing_functions import *
global focal_len


# inputs

# geometry ===============================================================================
"""
Notes:
    the vector between the element (x,y,z) coordinate and the element (x,y,z) 
    aimpoint defines the z-axis of the element coordinate system. The apertures 
    and surface shapes are then all are defined relative to the element 
    coordinate system.
"""

# from Zhang et al. 2022
# l_c = 900*2 / 1000 # module length [m]
# a_w = 1300 / 1000 # aperture width
# focal_len = 270 / 1000 # 1.71 # focal length # this must be correct for results to make sense
# d_abstube = 40 / 1000 #??? not found in paper # diameter of absorber tube
# abs_height = focal_len + d_abstube/2. # pt on upper?? sfc of abs tube
# ptc_pos = [0, 0] # x, y 
# ptc_aim = [0, 1] # x, z
# abs_aimz = 0.

# # optical properties
# refl_rho = 1 # 0.935 # trough
# absr_rho = 0. # receiver tube

# IDEAL CASE: using SGX-2 from http://edge.rit.edu/edge/P15484/public/Detailed%20Design%20Documents/Solar%20Trough%20Preliminary%20analysis%20references/Parabolic%20Trough%20Technology.pdf
# l_c = 12.0 # module length
# a_w = 5.77 # aperture width
# focal_len = 1.71 # focal length # this must be correct for results to make sense
# d_abstube = 0.07 # diameter of absorber tube
# abs_height = focal_len - d_abstube/2. # pt on upper?? sfc of abs tube
# ptc_pos = [0, 0] # x, y 
# ptc_aim = [0, 1] # x, z
# abs_aimz = focal_len*2. # 0.

# optical properties
refl_rho = 1. # trough
absr_rho = 0. # receiver tube

# recreating LS-3 standard case from SolTrace GUI (from Janna)
# reflecter
l_c = 10.0 # module length
a_w = 2.887*2 # aperture width
R = (1/0.292398) # from screenshot surface definition p-1/R
focal_len = R/2. #1.71 # focal length # this must be correct for results to make sense
ptc_pos = [0, 0] # x, y 
ptc_aim = [0, 1] # x, z
# absorber
r_abstube = 1/28.5714
d_abstube = r_abstube*2. #0.07 # diameter of absorber tube
abs_height = focal_len + d_abstube/2. # pt on upper?? sfc of abs tube
abs_aimz = 0. # focal_len*2. # 0.

# sun
sun_position = [0., 0., 100.] # x, y, z

# sim inputs
n_hits = [int(1e2),int(1e3),int(1e4),int(1e5),int(1e6)] #,1e3,1e4,1e5]) # 1e7 #1e7 #10 # 1e5
sunshape_flag = False
sfcerr_flag = False

global nx
global ny
# mesh definition for flux map
nx = 30
ny = 30

plotrays = False

max_flux = np.zeros(np.shape(n_hits))

for n,s in enumerate(n_hits):
    print('n,s',n,s)

    # simulation setup ==========================================================================
    # Create API class instance
    PT = PySolTrace()
    
    # Create two optics types - one for reflector, and one for absorber.
    opt_ref = PT.add_optic("Reflector")
    opt_ref.front.reflectivity = refl_rho # reflects all power
    
    opt_abs = PT.add_optic("Absorber")
    opt_abs.front.reflectivity = absr_rho
    opt_abs.back.reflectivity = absr_rho
    
    # Sun
    sun = PT.add_sun()
    # Give sun an arbitrary position
    sun.position.x = sun_position[0] #2. #0.
    sun.position.y = sun_position[1]
    sun.position.z = sun_position[2]
    
    # Reflector stage
    st = PT.add_stage()
    
    # Create a parabolic trough at x,y position, reflecting to the receiver
    el = st.add_element()
    el.optic = opt_ref
    el.position.x = ptc_pos[0]
    el.position.y = ptc_pos[1]
    el.aim.x = ptc_aim[0] #0.01
    el.aim.z = ptc_aim[1]
    # rvec = (abs_pos - el.position).unitize()
    # svec = sun.position.unitize()
    # avec = (rvec + svec)/2.
    # # assign the aim vector. scale by a large number
    # el.aim = el.position + avec*100.
    # # compute surface z rotation to align with plane of the ground
    # el.zrot = PT.util_calc_zrot_azel(avec)
    # Set surface and aperture characteristics
    el.surface_parabolic(focal_len, float('infinity')) # focal_len_y should be 'inf', but that threw an error
    el.aperture_rectangle(a_w,l_c)
        
    # absorber stage
    # sta = PT.add_stage()
    
    #absorber element height
    abs_pos = Point(ptc_pos[0], ptc_pos[1], abs_height)
    #abs_pos = Point(0., 0., 1.745)
    
    ela = st.add_element()
    ela.position = abs_pos
    ela.aim.z = abs_aimz
    ela.optic = opt_abs
    ela.surface_cylindrical(d_abstube/2.)
    ela.aperture_singleax_curve(0.,0.,l_c)
    
    # set simulation parameters
    PT.num_ray_hits = s #n_hits #10 # 1e5
    PT.max_rays_traced = PT.num_ray_hits*100 # 1000 #100 #PT.num_ray_hits*100
    PT.is_sunshape = sunshape_flag # True
    PT.is_surface_errors = sfcerr_flag # True
    
    if __name__ == "__main__":
    
        #PT.write_soltrace_input_file('LS3trough.stinput')
        PT.run(10, False, 4)         #(seed, is point focus system?, number of threads)
        #PT.run(1, False, 1)         #(seed, is point focus system?, number of threads)
        
        print("Num rays traced: {:d}".format(PT.raydata.index.size))
    
        df = PT.raydata
        
        #print("Did any rays miss receiver? Element unique values = {}".format(df[df.stage==2].element.unique()))
    
        if plotrays==True:
            # Data for a three-dimensional line
            loc_x = df.loc_x.values
            loc_y = df.loc_y.values
            loc_z = df.loc_z.values
    
            # Plotting with plotly
            import plotly.express as px 
            import plotly.graph_objects as go
            import plotly.io as io
            io.renderers.default='browser'
    
            fig = go.Figure(data=go.Scatter3d(x=loc_x, y=loc_y, z=loc_z, mode='markers', marker=dict( size=1, color=df.stage, colorscale='bluered', opacity=0.8, ) ))
    
            #for i in range(PT.raydata.index.size): # all rays
            for i in range(50,100):
            #for i in range(0,100000,500):
                dfr = df[df.number == i]
                ray_x = dfr.loc_x 
                ray_y = dfr.loc_y
                ray_z = dfr.loc_z
                raynum = dfr.number
                fig.add_trace(go.Scatter3d(x=ray_x, y=ray_y, z=ray_z, mode='lines', line=dict(color='black', width=0.5)))
    
            fig.update_layout(showlegend=False)
            fig.show()
        
        #print(df.describe())
        ppr = PT.powerperray
        df_rec = generate_receiver_dataframe(df,d_abstube,focal_len)
        flux_st = compute_fluxmap(ppr,df_rec,d_abstube,l_c,nx,ny,plotflag=True)
        max_flux[n] = np.max(flux_st)
        del df
        
        
plt.figure(figsize=[3,4],dpi=250)
plt.semilogx(n_hits, max_flux, '.-')
plt.xlabel('number of hits')    
plt.ylabel('max flux [kW/m2]')
        
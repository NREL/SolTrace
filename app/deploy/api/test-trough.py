from pysoltrace import PySolTrace, Point
import random
import pandas as pd
import copy

# geometric inputs
# using SGX-2 from http://edge.rit.edu/edge/P15484/public/Detailed%20Design%20Documents/Solar%20Trough%20Preliminary%20analysis%20references/Parabolic%20Trough%20Technology.pdf
l_c = 12.0 # module length
a_w = 5.77 # aperture width
focal_len = 1.71 # focal length # this must be correct for results to make sense
d_abstube = 0.07 # diameter of absorber tube

# Create API class instance
PT = PySolTrace()

# Create two optics types - one for reflector, and one for absorber.
opt_ref = PT.add_optic("Reflector")
opt_ref.reflectivity = 1.

opt_abs = PT.add_optic("Absorber")
opt_abs.reflectivity = 0.

# Sun
sun = PT.add_sun()
# Give sun an arbitrary position
sun.position.x = 0.
sun.position.y = 0.
sun.position.z = 100.

# Reflector stage
st = PT.add_stage()

# Create a parabolic trough at x,y position, reflecting to the receiver
el = st.add_element()
el.optic = opt_ref
el.position.x = 0
el.position.y = 0
el.aim.x = 0 #0.05
el.aim.z = 1
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
sta = PT.add_stage()

#absorber element height
abs_height = focal_len + d_abstube/2. # pt on upper?? sfc of abs tube
abs_pos = Point(0., 0., abs_height)
#abs_pos = Point(0., 0., 1.745)

ela = sta.add_element()
ela.position = abs_pos
ela.aim.z = 0.
ela.optic = opt_abs
ela.surface_cylindrical(d_abstube/2.)
ela.aperture_singleax_curve(0.,0.,l_c)

# set simulation parameters
PT.num_ray_hits = 1e5
PT.max_rays_traced = PT.num_ray_hits*100
PT.is_sunshape = True 
PT.is_surface_errors = True

if __name__ == "__main__":

    PT.run(10, False, 4)         #(seed, is point focus system?, number of threads)
    
    print("Num rays traced: {:d}".format(PT.raydata.index.size))

    df = PT.raydata
    # Data for a three-dimensional line
    loc_x = df.loc_x.values
    loc_y = df.loc_y.values
    loc_z = df.loc_z.values


#     # Plotting with plotly
#     import plotly.express as px 
#     import plotly.graph_objects as go
#     import plotly.io as io
#     io.renderers.default='browser'

#     fig = go.Figure(data=go.Scatter3d(x=loc_x, y=loc_y, z=loc_z, mode='markers', marker=dict( size=1, color=df.stage, colorscale='bluered', opacity=0.8, ) ) )

# #    for i in range(50,100):
#     for i in range(0,100000,500):
#         dfr = df[df.number == i]
#         ray_x = dfr.loc_x 
#         ray_y = dfr.loc_y
#         ray_z = dfr.loc_z
#         raynum = dfr.number
#         fig.add_trace(go.Scatter3d(x=ray_x, y=ray_y, z=ray_z, mode='lines', line=dict(color='black', width=0.5)))

#     fig.update_layout(showlegend=False)
#     fig.show()
    
    # compute flux distribution
    import numpy as np
    import math
    import matplotlib.pyplot as plt
    
    # copied from lines 74+ in https://github.com/NREL/SolarPILOT/blob/develop/deploy/api/test_solarpilot_soltrace.py
    dfr = df[df.stage==2]
    dfr = dfr[dfr.element==-1]  #absorbed rays
    
    d_rec = d_abstube
    h_rec = abs_height
    ny = 100
    nx = 100
    y_rec = np.arange(0, h_rec, h_rec/ny)
    x_rec = np.arange(0, np.pi*d_rec, np.pi*d_rec/nx)

    Xr,Yr = np.meshgrid(x_rec, y_rec)
    
    tht = 0. #cp.data_get_number(spcxt, "solarfield.0.tht") # what is this?
    dfr['zpos'] = dfr.loc_z - tht + h_rec/2.
    dfr['cpos'] = (np.arctan2(dfr.loc_x, dfr.loc_y)+math.pi)*d_rec/2.
    
    flux_st = np.zeros((ny,nx))
    dx = d_rec*np.pi / nx # circumference of receiver / nx
    dy = h_rec / ny # y is vertical direction here?
    anode = dx*dy # area of a node
    ppr = PT.powerperray / anode *1e-3 #st.powerperray / anode *1e-3

    for ind,ray in dfr.iterrows():

        j = int(ray.cpos/dx)
        i = int(ray.zpos/dy)

        flux_st[i,j] += ppr
    print(dfr.describe())
    
    plt.figure()
    plt.title("Flux simulation from the SolTrace engine")
    plt.contourf(Xr, Yr, flux_st, levels=25)
    plt.colorbar()
    plt.title(f"max flux {flux_st.max():.0f} kW/m2, mean flux {flux_st.mean():.1f}")
    plt.show()
    
    
    
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
opt_ref.reflectivity = 1. # reflects all power

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
el.aim.x = 0 #0.01
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
PT.num_ray_hits = 10 # 1e5
PT.max_rays_traced = 100 #PT.num_ray_hits*100
PT.is_sunshape = False # True
PT.is_surface_errors = False # True

if __name__ == "__main__":

    #PT.run(10, False, 4)         #(seed, is point focus system?, number of threads)
    PT.run(1, False, 1)         #(seed, is point focus system?, number of threads)
    
    print("Num rays traced: {:d}".format(PT.raydata.index.size))

    df = PT.raydata
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

    for i in range(PT.raydata.index.size): # all rays
    #for i in range(50,100):
    #for i in range(0,100000,500):
        dfr = df[df.number == i]
        ray_x = dfr.loc_x 
        ray_y = dfr.loc_y
        ray_z = dfr.loc_z
        raynum = dfr.number
        fig.add_trace(go.Scatter3d(x=ray_x, y=ray_y, z=ray_z, mode='lines', line=dict(color='black', width=0.5)))

    fig.update_layout(showlegend=False)
    fig.show()
    
    
    
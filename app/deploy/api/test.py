from pysoltrace import PySolTrace, Point
import random
import pandas as pd
import copy

# def load_system(PT):
# def load_system(ii):
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
sun.position.x = 100.
sun.position.y = 0.
sun.position.z = 100.

# Reflector stage
st = PT.add_stage()

# Simple reflector at origin
# el = st.add_element()
# el.aim.z = 1
# el.optic = opt_ref
# el.surface_flat()
# el.aperture_rectangle(1,1)

#absorber element height
abs_pos = Point(0., 0., 10.)

# Create a heliostat at some random x,y position, reflecting to the receiver
for i in range(5):
    hpos = [random.uniform(-10,10), random.uniform(-10,10)]
    el = st.add_element()
    el.optic = opt_ref
    el.position.x = hpos[0]
    el.position.y = hpos[1]
    # calculate the vectors - receiver, sun, and aim
    rvec = (abs_pos - el.position).unitize()
    svec = sun.position.unitize()
    avec = (rvec + svec)/2.
    # assign the aim vector. scale by a large number
    el.aim = el.position + avec*100.
    # compute surface z rotation to align with plane of the ground
    el.zrot = PT.util_calc_zrot_azel(avec)
    # Set surface and aperture characteristics
    el.surface_flat()
    el.aperture_rectangle(0.5,1.0)
    

# absorber stage
sta = PT.add_stage()

ela = sta.add_element()
ela.position = abs_pos
ela.aim.z = 0.
ela.optic = opt_abs
ela.surface_flat()
ela.aperture_rectangle(5,5)  #target is 5x5 

# set simulation parameters
PT.num_ray_hits = 1e5
PT.max_rays_traced = PT.num_ray_hits*100
PT.is_sunshape = True 
PT.is_surface_errors = True



if __name__ == "__main__":

    PT.run(-1, True, 4)         #(seed, is point focus system?, number of threads)
    
    print("Num rays traced: {:d}".format(PT.raydata.index.size))

    df = PT.raydata
    # Data for a three-dimensional line
    loc_x = df.loc_x.values
    loc_y = df.loc_y.values
    loc_z = df.loc_z.values


    # Plotting with plotly
    import plotly.express as px 
    import plotly.graph_objects as go

    fig = go.Figure(data=go.Scatter3d(x=loc_x, y=loc_y, z=loc_z, mode='markers', marker=dict( size=1, color=df.stage, colorscale='bluered', opacity=0.8, ) ) )

    for i in range(50,100):
        dfr = df[df.number == i]
        ray_x = dfr.loc_x 
        ray_y = dfr.loc_y
        ray_z = dfr.loc_z
        raynum = dfr.number
        fig.add_trace(go.Scatter3d(x=ray_x, y=ray_y, z=ray_z, mode='lines', line=dict(color='black', width=0.5)))

    fig.update_layout(showlegend=False)
    fig.show()
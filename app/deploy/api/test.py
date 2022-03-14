from pysoltrace import PySolTrace, Point
import random

# Create API class instance
PT = PySolTrace()

# Create two optics types - one for reflector, and one for absorber.
opt_ref = PT.add_optic("Reflector")
opt_ref.reflectivity = 1.
opt_ref.Create()

opt_abs = PT.add_optic("Absorber")
opt_abs.reflectivity = 0.
opt_abs.Create()

# Sun
sun = PT.add_sun()
# Give sun an arbitrary position
sun.position.x = -20.
sun.position.y = 0
sun.position.z = 100.
sun.Create()

# Reflector stage
st = PT.add_stage()
st.Create()

# Simple reflector at origin
# st.add_elements()
# el = st.elements[-1]
# el.aim.z = 1
# el.optic = opt_ref
# el.Create()

#absorber element height
abs_pos = Point(0., 0., 10.)

# Create a heliostat at some random x,y position, reflecting to the receiver
for i in range(5):
    hpos = [random.uniform(-10,10), random.uniform(-10,10)]
    el = st.add_element()
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
    
    el.optic = opt_ref
    el.Create()

# absorber stage
sta = PT.add_stage()
sta.Create()

ela = sta.add_element()
ela.position = abs_pos
ela.aim.z = 0.
ela.optic = opt_abs
ela.aperture_rectangle(5,5)  #target is 5x5 
ela.Create()

# set simulation parameters
PT.num_ray_hits = 1e5
PT.max_rays_traced = PT.num_ray_hits*100
PT.is_sunshape = True 
PT.is_surface_errors = True
PT.run(-1, True)

df = PT.get_ray_dataframe()

PT.clear_context()


# ---------------------------------------------
# create plots of the ray data
# import matplotlib.pyplot as plt 
# fig,axs = plt.subplots(1,3, figsize=(12,4))
# dfe = df #[df.element < 0].reindex()      # rays that are absorbed have element<0
# dfes1 = dfe[df.stage == 1].reindex()    # rays in stage 1
# dfes2 = dfe[df.stage == 2].reindex()    # rays in stage 2

# # Scatter plot in each plane
# dfes1.plot('loc_x', 'loc_y', style='k.',ax=axs[0], markersize=0.05)
# dfes2.plot('loc_x', 'loc_y', style='r.',ax=axs[0], markersize=0.05)
# dfes1.plot('loc_x', 'loc_z', style='k.',ax=axs[1], markersize=0.05)
# dfes2.plot('loc_x', 'loc_z', style='r.',ax=axs[1], markersize=0.05)
# dfes1.plot('loc_y', 'loc_z', style='k.',ax=axs[2], markersize=0.05)
# dfes2.plot('loc_y', 'loc_z', style='r.',ax=axs[2], markersize=0.05)

# # Plot trace of selected rays, by number
# for i in range(50,70):
#     dfi = df[df.number == i]
#     dfi.plot('loc_x', 'loc_y', style='b-', ax=axs[0], linewidth=0.1)
#     dfi.plot('loc_x', 'loc_z', style='b-', ax=axs[1], linewidth=0.1)
#     dfi.plot('loc_y', 'loc_z', style='b-', ax=axs[2], linewidth=0.1)
# # no legend
# for ax in axs:
#     ax.get_legend().remove()

# plt.show()



# Data for a three-dimensional line
loc_x = df.loc_x.values
loc_y = df.loc_y.values
loc_z = df.loc_z.values


import plotly.express as px 
import plotly.graph_objects as go

# fig = px.scatter_3d(df, x='loc_x', y='loc_y', z='loc_z', color='stage', size='element') 

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
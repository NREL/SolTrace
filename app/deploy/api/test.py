from pysoltrace import PySolTrace

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

abs_z = 10.  #absorber element height

# Create a heliostat at some x,y,z position, reflecting to the receiver
for hpos in [[5,5], [5,-5], [-5,5], [-5,-5]]:   # picking some positions to test
    st.add_elements()
    el = st.elements[-1]
    el.position.x = hpos[0]
    el.position.y = hpos[1]
    # calculate the vectors - receiver, sun, and aim
    rvec = PT.util_calc_unitvect([-el.position.x,-el.position.x, abs_z])
    svec = PT.util_calc_unitvect([sun.position.x, sun.position.y, sun.position.z])
    avec = [(rvec[0]+svec[0])/2., (rvec[1]+svec[1])/2., (rvec[2]+svec[2])/2.]
    # assign the aim vector. scale by a large number
    el.aim.x = el.position.x + avec[0]*100.
    el.aim.y = el.position.y + avec[1]*100.
    el.aim.z = el.position.z + avec[2]*100.
    # compute surface z rotation to align with plane of the ground
    el.zrot = PT.util_calc_zrot_azel(avec)
    
    el.optic = opt_ref
    el.Create()

# absorber stage
sta = PT.add_stage()
sta.Create()

sta.add_elements()
ela = sta.elements[-1]
ela.position.z = abs_z
ela.aim.z = 0.
ela.optic = opt_abs
ela.aperture_rectangle(5,5)  #receiver is 5x5 
ela.Create()

# set simulation parameters
PT.num_ray_hits = 1e5
PT.max_rays_traced = PT.num_ray_hits*100
PT.is_sunshape = True 
PT.is_surface_errors = True
PT.run(-1, True)

df = PT.get_ray_dataframe()

PT.clear_context()

import matplotlib.pyplot as plt 
fig,axs = plt.subplots(1,3, figsize=(12,4))
df[df.element < 0].plot('loc_x', 'loc_y', style='k.',ax=axs[0])
df[df.element < 0].plot('loc_x', 'loc_z', style='k.',ax=axs[1])
df[df.element < 0].plot('loc_y', 'loc_z', style='k.',ax=axs[2])
plt.show()




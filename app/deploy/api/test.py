from pysoltrace import PySolTrace

import os 
print(os.getpid())

PT = PySolTrace()
opt_ref = PT.add_optic("Reflector")
opt_ref.reflectivity = 1.
opt_ref.Create()

opt_abs = PT.add_optic("Absorber")
opt_abs.reflectivity = 0.
opt_abs.Create()

# Sun
sun = PT.add_sun()
sun.position.z = 100.
sun.Create()

# Reflector stage
st = PT.add_stage()
st.Create()
st.add_elements()
el = st.elements[-1]
el.aim.z = 1
el.optic = opt_ref
el.Create()

# absorber stage
sta = PT.add_stage()
sta.Create()
sta.add_elements()
ela = sta.elements[-1]
ela.position.z = 1.
ela.aim.z = 0.
ela.optic = opt_abs
ela.Create()

ns = PT.get_num_stages()
nel = PT.stages[0].get_num_elements()

PT.run(10000, 1e6)

# nints = PT.get_num_intersections()
# ints = PT.get_intersect_locations()

df = PT.get_ray_dataframe()

xx=1

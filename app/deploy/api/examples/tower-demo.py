"""
This file demonstrates the modeling of a simple 3-mirror and tower system. 
The mirrors are placed at regular positions in a circle. The receiver can
be a flat surface or a cylindrical surface. The simulation runs in a few
seconds and generates a flux map and a trace plot (requires plotly).
"""

# Load the pysoltrace api from the parent directory ---
import sys
import os
sys.path.insert(1, os.path.join(sys.path[0], '..'))
# ----------

from pysoltrace import PySolTrace, Point
from math import sin,cos, pi


# Create API class instance
PT = PySolTrace()

# Create two optics types - one for reflector, and one for absorber.
opt_ref = PT.add_optic("Reflector")
opt_ref.front.reflectivity = 1.

opt_abs = PT.add_optic("Absorber")
opt_abs.front.reflectivity = 0.

# Sun
sun = PT.add_sun()
# Give sun an arbitrary position
sun.position.x = 0.
sun.position.y = 0.
sun.position.z = 100.

# Reflector stage
st = PT.add_stage()

#absorber element height
abs_pos = Point(0., 0., 10.)

# Create a heliostat at some random x,y position, reflecting to the receiver
for i in range(-1,2):
    hpos = [sin(i*pi/2)*5, cos(i*pi/2)*5]
    # hpos = [random.uniform(-10,10), random.uniform(-10,10)]
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
    el.aperture_rectangle(1.0,1.95)
    

sta = PT.add_stage()
    
# cylindrical absorber element. 
# r = 0.5
# ela = sta.add_element()
# abs_pos.y = -r
# ela.position = abs_pos
# ela.zrot = 0
# ela.aim = Point(0, 100, -r-50)
# ela.optic = opt_abs
# ela.surface_cylindrical(r)
# ela.aperture_singleax_curve(0,0,5)

# flat absorber element
ela = sta.add_element()
ela.position = abs_pos
ela.aim.z = 0.
ela.aim.x = 0.
ela.aim.y = 5.
ela.optic = opt_abs
ela.surface_flat()
ela.aperture_rectangle(2,2)  

# set simulation parameters
PT.num_ray_hits = 1e6
PT.max_rays_traced = PT.num_ray_hits*100
PT.is_sunshape = True 
PT.is_surface_errors = True

# Simulation needs to be inside of a __name__ guard to allow multi-threading
if __name__ == "__main__":

    # Run the configuration specified above
    PT.run(-1, True, 8)         #(seed, is point focus system?, number of threads)

    # Print a message after completion
    print("Num rays traced: {:d}".format(PT.raydata.index.size))
    
    # Generate a solatrace input file if desired
    PT.write_soltrace_input_file('simpletest.stinput')

    # Create a 3D trace plot of the simulation (requires plotly library)
    PT.plot_trace()

    # Plot the flux map on the last element of the last stage (the receiver in this case)
    PT.plot_flux(PT.stages[-1].elements[-1], nx=50, ny=50, levels=50)

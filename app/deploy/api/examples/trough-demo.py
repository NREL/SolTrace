"""
This file demonstrates the modeling of a simple parabolic trough collector. 
The model is used to predict loss as a function of tracking error by
varying the position of the sun while keeping the collector fixed.

The set of simulations take about 3 minutes to run.
"""

from pysoltrace import PySolTrace, Point

# Create API class instance
PT = PySolTrace()

# Create two optics types - one for reflector, and one for absorber.
opt_ref = PT.add_optic("Reflector")
opt_ref.front.reflectivity = 1.
opt_ref.front.slope_error = 2.

opt_abs = PT.add_optic("Absorber")
opt_abs.front.reflectivity = 0.
opt_abs.back.reflectivity = 0.

# Sun
sun = PT.add_sun()
# Give sun a position
sun.position.x = 0.
sun.position.y = 0.
sun.position.z = 100.

# Trough system parameters
W_ap = 2.887*2  #[m] aperture width
L_f = 1.71  #focal length
D_abs = 0.07 #[m] absorber tube diameter 

# Reflector stage
st = PT.add_stage()

el = st.add_element()
el.optic = opt_ref
el.position.x = 0
el.position.y = 0
el.position.z = -L_f

# assign the aim vector. scale by a large number
el.aim = Point(0, 0, 100)
el.zrot = 0.
# Set surface and aperture characteristics
el.surface_parabolic(L_f, float('inf'))
el.aperture_singleax_curve(-W_ap/2, W_ap/2, 10.)

# cylindrical absorber element. 
sta = PT.add_stage()
sta.is_tracethrough = True
ela = sta.add_element()
ela.position = Point(0., 0, -D_abs/2)
ela.zrot = 0
ela.aim = Point(0, 0, 100)
ela.optic = opt_abs
ela.surface_cylindrical(D_abs/2)
ela.aperture_singleax_curve(0,0,10)

# spillage catcher
stb = PT.add_stage()
els = stb.add_element()
els.position = Point(0., 0., 0.1)
els.aim = Point(0,0,0)
els.optic = opt_abs
els.surface_flat()
els.aperture_rectangle(W_ap, 10)

# set simulation parameters
PT.num_ray_hits = 5e5
PT.max_rays_traced = PT.num_ray_hits*100
PT.is_sunshape = True 
PT.is_surface_errors = True

# Simulation needs to be inside of a __name__ guard to allow multi-threading
if __name__ == "__main__":

    import numpy as np 
    import matplotlib.pyplot as plt 

    # Run a sweep over several surface errors and defocus angles
    ns = 4
    surferrs = np.linspace(1, 4, ns)
    nd = 15
    defocus = np.linspace(0, 4, nd)  #x offset / 100 equals angular displacement in rad
    f_inter = np.zeros((nd,ns))

    for j in range(len(surferrs)):  #iterate over surface errors

        PT.optics[0].front.slope_error = surferrs[j]  #assign surface error for collector

        for i in range(len(defocus)): #iterate over defocus angles
            
            PT.sun.position.x = defocus[i]  #assign the defocus angle

            # Run the configuration 
            PT.run(nthread=4)         
            
            #compute the intercept factor by examining the spillage
            df = PT.raydata 
            f_inter[i,j] = 1 - df[df.stage==3].size / df[df.stage==1].size

        plt.plot(defocus, f_inter[:,j], label=f'error={surferrs[j]:.1f} mrad')
    plt.legend()
    plt.xlabel("Defocusing angle (mrad)")
    plt.ylabel("Intercept factor")
    plt.savefig('trough-demo.png', dpi=200)
    plt.show()

    # # Generate a solatrace input file if desired
    # PT.write_soltrace_input_file('troughtest.stinput')

    # # Create a 3D trace plot of the simulation (requires plotly library)
    # PT.plot_trace()

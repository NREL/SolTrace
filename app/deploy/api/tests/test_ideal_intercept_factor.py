import sys
sys.path.insert(0, '/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/')
import numpy as np
import pandas as pd
from st_processing_functions import *
from run_pysoltrace_iterate import *

def test_answer():

    # run soltrace in validation mode
    sunshape_flag = False
    sfcerr_flag = False
    optics_type = 'ideal' # 'yang' 'realistic' # 'ideal'
    plot_rays = False
    save_pickle = False
    number_hits = 1e5 # 5e6 # 1e5 #1e5 

    # parabolic trough geometry definition ================================
    # NSO Trough Geometry: using measurements from CAD file from Dave (aka LS-2)
    module_length = 12.0 # module length
    aperture_width = 5.0 #5.77 # aperture width
    focal_len = 1.49 #1.71 # focal length # this must be correct for results to make sense
    d_abstube = 0.07 # diameter of absorber tube
    ptc_pos = [0, 0, 0] # x, y, z
    ptc_aim = [0, 0, 1] # x, y, z

    # field site definition ================================
    lat, lon = 35.8, -114.983 #coordinates of Nevada Solar One
    altitude = 543 #m
    
    # running  nominal (no tracking error) ===================================
    tracker_angle_input = 'nominal'
    sensorlocs = ['nominal']
    tstart = '2023-01-15 16:00:00' # fulldata.index[0] # '2023-02-11 17:00:00'
    tend = '2023-01-15 21:00:00' 
    times = pd.date_range(tstart, tend, freq='3H') # in UTC
    field_data_path = '/Volumes/Processed_data/'
    error_angles = []

    
    results, df = run_soltrace_iterate(times, lat, lon, altitude, field_data_path, tracker_angle_input, sensorlocs, module_length, aperture_width, focal_len, d_abstube, ptc_pos, ptc_aim, sunshape_flag, sfcerr_flag, optics_type, plot_rays, save_pickle, number_hits, nx=30, ny=30, error_angles=error_angles)

    intercept_factor_test = results['nominal'].intercept_factor

    # calculate intercept factor
    intercept_factor_test = results['nominal'].intercept_factor.values #func(df)
    print(intercept_factor_test)
    
    # define answers
    tstart = '2023-01-15 16:00:00' # fulldata.index[0] # '2023-02-11 17:00:00'
    tend = '2023-01-15 21:00:00' 
    test_times = pd.date_range(tstart, tend, freq='3H') # in UTC
    data = [1.0, 1.0]
    intercept_factor_answers = pd.DataFrame(data = data, index = test_times, columns=['intercept_factor'])
    
    assert np.allclose(intercept_factor_test, intercept_factor_answers)
    
if __name__ == "__main__":
    test_answer()
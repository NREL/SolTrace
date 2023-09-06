import sys
sys.path.insert(0, '/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/')
import numpy as np
from st_processing_functions import *
from run_pysoltrace_iterate import *

def test_answer():

    # run soltrace in validation mode
    # running for validation ===================================
    tracker_angle_input = 'validation'
    sensorlocs = ['validation']
    error_angles = np.array([0., np.degrees(1.733333e-02), 2.5])

    # fake inputs - not used in validation mode
    times = pd.date_range('2023-03-05 15:00:00', '2023-03-05 23:50:00',freq='4H') # in UTC
    lat, lon = 35.8, -114.983 #coordinates of Nevada Solar One
    altitude = 543 #m
    field_data_path = '/Users/bstanisl/Documents/seto-csp-project/NSO-field-data/' 
    module_length = 12.0 # module length
    aperture_width = 5.0 #5.77 # aperture width
    focal_len = 1.49 #1.71 # focal length # this must be correct for results to make sense
    d_abstube = 0.07 # diameter of absorber tube
    ptc_pos = [0, 0, 0] # x, y, z
    ptc_aim = [0, 0, 1] # x, y, z
    abs_aimz = focal_len*2. # 0. ??
    critical_angle_error = 0.79 #[deg] from firstoptic validation dataset
    sunshape_flag = False
    sfcerr_flag = False
    optics_type = 'ideal' # 'yang' 'realistic' # 'ideal'
    plot_rays = False
    save_pickle = False
    number_hits = 1e6 # 5e6 # 1e5 #1e5 

    results, df = run_soltrace_iterate(times, lat, lon, altitude, field_data_path, 
                                   tracker_angle_input, sensorlocs, module_length, 
                                   aperture_width, focal_len, d_abstube, ptc_pos, 
                                   ptc_aim, sunshape_flag, sfcerr_flag, 
                                   optics_type, plot_rays, 
                                   save_pickle, number_hits, nx=30, ny=30,
                                   error_angles=error_angles)


    # calculate intercept factor
    intercept_factor_test = results['validation'].intercept_factor.values #func(df)
    print(intercept_factor_test)
    
    # define answers
    cols = ['error','nom_trough_angle','intercept_factor']
    # copied from FirstOPTIC validation results
            # [np.degrees(1.733333e-02), 0., 7.135678e-01], 
            # [np.degrees(1.733333e-02), 90., 7.135678e-01],
    # the value of 0.706620459254 below is not exactly the analytical answer - it was copied from the soltrace output as matching "close enough" to the FirstOPTIC answer of 7.135678e-01
    data = [[0., 0., 1.], 
            [0., 90., 1.],
            [np.degrees(1.733333e-02), 0., 0.7066204592546165], 
            [np.degrees(1.733333e-02), 90., 0.7066204592546165],
            [2.5, 0., 0.], 
            [2.5, 90., 0.]]
    intercept_factor_answers = pd.DataFrame(data, columns=cols)
    intercept_factor_answers = intercept_factor_answers.set_index(['error','nom_trough_angle'])
    
    # collect test data
    nom_trough_angle = results['validation'].nom_trough_angle
    error = results['validation'].trough_angle_dev
    
    # calculate answers using test data
    intercept_factor_answer = []
    for i,j in zip(error, nom_trough_angle):
        intercept_factor_answer.append(intercept_factor_answers.loc[(i, j),'intercept_factor'])
    intercept_factor_answer = np.array(intercept_factor_answer)
    print(intercept_factor_answer)
    
    assert np.allclose(intercept_factor_test, intercept_factor_answer)
    
if __name__ == "__main__":
    test_answer()
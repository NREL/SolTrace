# content of test_sample.py
import pandas as pd
import numpy as np
import sys
sys.path.insert(0, '/Users/bstanisl/Documents/seto-csp-project/SolTrace/SolTrace/app/deploy/api/')
from st_processing_functions import get_trough_angles_py

def func(results):
    return results.projected_sun_angle.values

def test_answer():
    times = pd.date_range('2023-03-05 15:00:00', '2023-03-05 23:50:00',freq='4H') # in UTC
    lat, lon = 35.8, -114.983 #coordinates of Nevada Solar One
    altitude = 543 #m
    
    sunangles = get_trough_angles_py(times, lat, lon, altitude)
    answers = [79.385638, 16.929936, -56.151902]
    assert np.allclose(func(sunangles), answers)
    # for sloc in ['R1_DO']: # results.keys():
    #     assert func(sunangles) == answers
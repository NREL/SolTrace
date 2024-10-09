import sys
sys.path.append('../')
sys.path.append('../../../')
import numpy as np
import pandas as pd
from st_processing_functions import *

def f(nx,ny):
    
    r_abstube = 1/28.5714
    d_abstube = r_abstube*2. #0.07 # diameter of absorber tube
    R = (1/0.292398) # from screenshot surface definition p-1/R
    focal_len = R/2. #1.71 # focal length # this must be correct for results to make sense
    abs_height = focal_len + d_abstube/2. # pt on upper?? sfc of abs tube
    
    # define ray df
    # 3 ray intersections, one on leftmost pt of absorber, one on right most, one on lower
    
    # loc_x  loc_y  loc_z
    # -d_abs/2  0    0
    # 0         0    -d_abs/2
    # d_abs/2   0    0
    
    data = np.array([[-r_abstube, 0., focal_len], 
                    [0., 0., focal_len-r_abstube],
                    [r_abstube, 0., focal_len]])

    dfr = pd.DataFrame(data=data, columns=['loc_x','loc_y','loc_z'])
    
    ppr = 10.
    l_c = 1.0
    # nx = 5
    # ny = 4
    
    dfr['ypos'] = dfr.loc_y
    dfr['cpos'] = convert_xy_polar_coords(d_abstube,dfr.loc_x,dfr.loc_z,focal_len,gui_coords=True)
    print(dfr)
    
    flux_st, c_v = compute_fluxmap(ppr,dfr,d_abstube,l_c,nx,ny,plotflag=True)
    print(flux_st)

    # even ny

    return flux_st

def test_answer():

    nx = 5
    ny = 4
    if ny == 4:
        flux_answer = np.array([[0.,         0.,         0.,         0.        ],
                                [0.,         0.90945591, 0.,         0.        ],
                                [0.,         0.,         0.,         0.        ],
                                [0.,         0.90945591, 0.,         0.        ],
                                [0.,         0.90945591, 0.,         0.        ]])
    flux_st0 = f(nx,ny)       
    assert np.allclose(flux_st0, flux_answer)
    
    nx = 5
    ny = 3
    if ny == 3:
        flux_answer = np.array([[0.,         0.,         0.        ],
                                [0.,         0.68209193, 0.        ],
                                [0.,         0.,         0.        ],
                                [0.,         0.68209193, 0.        ],
                                [0.,         0.68209193, 0.        ]])
    flux_st1 = f(nx,ny)       
    assert np.allclose(flux_st1, flux_answer)
    
if __name__ == "__main__":
    test_answer()
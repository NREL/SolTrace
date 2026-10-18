import math, sys, os
sys.path.insert(1, os.path.join(sys.path[0], '../..'))

import numpy as np

from pysoltrace import PySolTrace as pst, Point

CONVERT_GLOBAL = np.array([[-1., 0., 0.],
                           [ 0., 0., 1.],
                           [ 0., 1., 0.]])

# NSTTF Coordinates (+X = West, +Y = Up, + Z = North)
# current values from g3p3_example_2mw_multifacet.stinput
# target_stage_center = np.array([40., 47.8177, 8.5])
target_stage_center = Point(40., 47.8177, 8.5)
# target_stage_aim = np.array([0, target_stage_center[1], 122])
target_stage_aim = Point(0, target_stage_center[1], 122)

target_euler = pst.util_calc_euler_angles(target_stage_center, target_stage_aim, 0)
target_transforms = pst.util_calc_transforms(target_euler)

# helio_stage_center = np.array([0., 0., 0.])
helio_stage_center = Point(0., 0., 0.)
# helio_stage_aim = np.array([0, 0, 1])
helio_stage_aim = Point(0, 0, 1)

helio_euler = pst.util_calc_euler_angles(helio_stage_center, helio_stage_aim, 0)
helio_transforms = pst.util_calc_transforms(helio_euler)

if __name__ == '__main__':
    print('\n\ndemo')
    print(target_euler)
    print(np.linalg.norm(target_euler))
    print(target_transforms['rloctoref'])
    print(helio_euler)
    print(np.linalg.norm(helio_euler))
    print(helio_transforms['rloctoref'])
    unstaged_pos = target_transforms['rreftoloc'] @ Point(0, 0, 8) + target_stage_center
    unstaged_aim = target_transforms['rreftoloc'] @ Point(0, 0, 9) + target_stage_center
    print(unstaged_pos)
    print(unstaged_aim)
    print(CONVERT_GLOBAL @ unstaged_pos)
    print(CONVERT_GLOBAL @ unstaged_aim)
import ctypes
from typing import Literal

from pysoltrace import dot_h
from pysoltrace.api.utils import st_function
from pysoltrace.api.dll import context
from pysoltrace.api.result.csv import csv

#############################################
# functions for SolTrace results management #
#############################################
class result(context):
    __slots__ = ('csv')
    
    def __init__(self, pdll, pcxt):
        super().__init__(pdll, pcxt)
        self.csv = csv(pdll, pcxt)

    @st_function
    def num(self):
        pcount = ctypes.c_uint64()
        code = self._pdll.st_num_intersections(self._pcxt, 
                                               ctypes.byref(pcount))
        return code, pcount.value

    def __len__(self):
        if not self._pdll or not self._pcxt:
            return 0
        return self.num()

    def __get_num(self):
        n = self.num()
        assert n > 0, 'No results to report.'
        return n

    @st_function    
    def locations(self):
        n = self.__get_num()
        loc_x, loc_y, loc_z = \
            ((ctypes.c_double * n)() for _ in range(3))
        code = self._pdll.st_locations(self._pcxt,
                                       loc_x, loc_y, loc_z)
        return code, loc_x[:n], loc_y[:n], loc_z[:n]

    @st_function
    def cosines(self):
        n = self.__get_num()
        cos_x, cos_y, cos_z = \
            ((ctypes.c_double * n)() for _ in range(3))
        code = self._pdll.st_cosines(self._pcxt,
                                     cos_x, cos_y, cos_z)
        return code, cos_x[:n], cos_y[:n], cos_z[:n]
    
    @st_function
    def elementmap(self):
        n = self.__get_num()
        element_map = (ctypes.c_uint64 * n)()
        code = self._pdll.st_elementmap(self._pcxt, element_map)
        return code, element_map[:n]

    @st_function
    def stagemap(self):
        n = self.__get_num()
        stage_map = (ctypes.c_uint64 * n)()
        code = self._pdll.st_stagemap(self._pcxt, stage_map)
        return code, stage_map[:n]

    @st_function
    def raynumbers(self):
        n = self.__get_num()
        ray_numbers = (ctypes.c_uint64 * n)()
        code = self._pdll.st_raynumbers(self._pcxt, ray_numbers)
        return code, ray_numbers[:n]

    @st_function
    def sun_stats(self):
        width, height, area = \
            (ctypes.c_double() for _ in range(3)) 
        nsunrays = ctypes.c_uint64()
        code = self._pdll.st_sun_stats(self._pcxt,
                                       ctypes.byref(width),
                                       ctypes.byref(height),
                                       ctypes.byref(area),
                                       ctypes.byref(nsunrays))
        return code, width.value, height.value, area.value, nsunrays.value

    @st_function
    def get(self):
        n = self.__get_num()
        loc_x, loc_y, loc_z, coz_x, coz_y, coz_z = \
            ((ctypes.c_double * n)() for _ in range(6))
        element_map, stage_map, ray_numbers= \
            ((ctypes.c_uint64 * n)() for _ in range(3))
        args = dot_h.args_results_data(loc_x, loc_y, loc_z,
                                       coz_x, coz_y, coz_z,
                                       element_map,
                                       stage_map,
                                       ray_numbers)
        code = self._pdll.st_get_results_data(self._pcxt, ctypes.byref(args))
        return code, { k: v[:n] for k, v in args.value.items() }

__all__ = ['result']
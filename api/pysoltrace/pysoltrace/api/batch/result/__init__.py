import ctypes

from pysoltrace import dot_h
from pysoltrace.api.batch.utils import batcher, generate_api_call
from pysoltrace.api.batch.result.csv import csv

# short names
_NUM = dot_h.st_api_call.CALL_ST_NUM_INTERSECTIONS
_LOC = dot_h.st_api_call.CALL_ST_LOCATIONS
_COS = dot_h.st_api_call.CALL_ST_COSINES
_EL  = dot_h.st_api_call.CALL_ST_ELEMENTMAP
_STA = dot_h.st_api_call.CALL_ST_STAGEMAP
_RAY = dot_h.st_api_call.CALL_ST_RAYNUMBERS
_SUN = dot_h.st_api_call.CALL_ST_SUN_STATS
_GET = dot_h.st_api_call.CALL_ST_GET_RESULTS_DATA

#############################################
# functions for SolTrace results management #
#############################################
class result(batcher):
    __slots__ = ('csv')

    def __init__(self, adder):
        super().__init__(adder)
        self.csv = csv(adder)

    def num(self):
        pcount = ctypes.c_uint64()
        call = generate_api_call(_NUM, ctypes.pointer(pcount))
        return self.adder(call, lambda: pcount.value)
    
    def locations(self, n: int):
        loc_x, loc_y, loc_z = \
            ((ctypes.c_double * n)() for _ in range(3))
        call = generate_api_call(_LOC, loc_x, loc_y, loc_z)
        return self.adder(call, lambda: (loc_x[:n], loc_y[:n], loc_z[:n]))

    def cosines(self, n: int):
        cos_x, cos_y, cos_z = \
            ((ctypes.c_double * n)() for _ in range(3))
        call = generate_api_call(_COS, cos_x, cos_y, cos_z)
        return self.adder(call, lambda: (cos_x[:n], cos_y[:n], cos_z[:n]))
    
    def elementmap(self, n: int):
        element_map = (ctypes.c_uint64 * n)()
        call = generate_api_call(_EL, element_map)
        return self.adder(call, lambda: element_map[:n])

    def stagemap(self, n: int):
        stage_map = (ctypes.c_uint64 * n)()
        call = generate_api_call(_STA, stage_map)
        return self.adder(call, lambda: stage_map[:n])

    def raynumbers(self, n: int):
        ray_numbers = (ctypes.c_uint64 * n)()
        call = generate_api_call(_RAY, ray_numbers)
        return self.adder(call, lambda: ray_numbers[:n])

    def sun_stats(self):
        width, height, area = \
            (ctypes.c_double() for _ in range(3)) 
        nsunrays = ctypes.c_uint64()
        call = generate_api_call(_SUN,
                                 ctypes.pointer(width),
                                 ctypes.pointer(height),
                                 ctypes.pointer(area),
                                 ctypes.pointer(nsunrays))
        return self.adder(call, lambda: (width.value, \
                                         height.value, \
                                         area.value, \
                                         nsunrays.value))

    def get(self, n: int):
        loc_x, loc_y, loc_z, coz_x, coz_y, coz_z = \
            ((ctypes.c_double * n)() for _ in range(6))
        element_map, stage_map, ray_numbers= \
            ((ctypes.c_uint64 * n)() for _ in range(3))
        args = dot_h.args_results_data(loc_x, loc_y, loc_z,
                                       coz_x, coz_y, coz_z,
                                       element_map,
                                       stage_map,
                                       ray_numbers)
        call = generate_api_call(_GET, ctypes.pointer(args))
        return self.adder(call, lambda: args.value)

__all__ = ['result']
import ctypes

from pysoltrace import dot_h, soltrace_structs as _STS
from pysoltrace.api.batch.utils import batcher, generate_api_call

# short names
_SET  = dot_h.st_api_call.CALL_ST_SET_SIMULATION_PARAMETERS
_RAYS = dot_h.st_api_call.CALL_ST_SIM_RAYS
_PT   = dot_h.st_api_call.CALL_ST_SIM_POWER_TOWER
_ERR  = dot_h.st_api_call.CALL_ST_SIM_ERRORS
_LOC  = dot_h.st_api_call.CALL_ST_SIM_LOCATION
_TOL  = dot_h.st_api_call.CALL_ST_SIM_TOLERANCE

##################################################
# functions for simulation parameters management #
##################################################
class parameters(batcher):
    def set(self, params: _STS.simulation_parameters) -> None:
        return self.adder(generate_api_call(_SET, ctypes.pointer(params.ctype)))

    def rays(self, raycount: int, maxcount: int) -> None:
        return self.adder(generate_api_call(_RAYS, raycount, maxcount))

    def power_tower(self, is_power_tower: bool) -> None:
        return self.adder(generate_api_call(_PT, is_power_tower))

    def errors(self, sun_shape: bool, optical: bool) -> None:
        return self.adder(generate_api_call(_ERR, sun_shape, optical))

    def location(self, latitude: float, longitude: float) -> None:
        return self.adder(generate_api_call(_LOC, latitude, longitude))

    def tolerance(self, tolerance: float) -> None:
        return self.adder(generate_api_call(_TOL, tolerance))

    # TODO:
    # def get(self) -> _STS.simulation_parameters:
    #     params = dot_h.args_simulation_parameters()
    #     code = self._pdll.st_get_simulation_parameters(self._pcxt,
    #                                                     ctypes.byref(params))
    #     return code, params.value

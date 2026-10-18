import ctypes
from typing import Literal

from pysoltrace import dot_h, soltrace_structs as _STS
from pysoltrace.api.utils import st_function
from pysoltrace.api.dll import context

##################################################
# functions for simulation parameters management #
##################################################
class parameters(context):
    @st_function
    def set(self, params: _STS.simulation_parameters) -> None:
        return self._pdll.st_set_simulation_parameters(self._pcxt,
                                                       ctypes.byref(params.ctype))

    @st_function
    def rays(self, raycount: int, maxcount: int) -> None:
        return self._pdll.st_sim_rays(self._pcxt, raycount, maxcount)

    @st_function
    def power_tower(self,
                    is_power_tower: bool) -> None:
        return self._pdll.st_sim_power_tower(self._pcxt, is_power_tower)

    @st_function
    def errors(self, sun_shape: bool, optical: bool) -> None:
        return self._pdll.st_sim_errors(self._pcxt, sun_shape, optical)

    @st_function
    def location(self, latitude: float, longitude: float) -> None:
        return self._pdll.st_sim_location(self._pcxt, latitude, longitude)

    @st_function
    def tolerance(self, tolerance: float) -> None:
        return self._pdll.st_sim_tolerance(self._pcxt, tolerance)

    @st_function
    def get(self) -> _STS.simulation_parameters:
        params = dot_h.args_simulation_parameters()
        code = self._pdll.st_get_simulation_parameters(self._pcxt,
                                                        ctypes.byref(params))
        return code, params.value

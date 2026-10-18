import ctypes
from typing import Literal

from pysoltrace import dot_h
from pysoltrace.api.utils import st_function
from pysoltrace.api.dll import context

class legacy(context):
    @st_function    
    def sim_params(self,
                   raycount: int,
                   maxcount: int,
                   include_dynamic_group: bool) -> None:
        return self._pdll.st_sim_params(self._pcxt, 
                                        raycount,
                                        maxcount,
                                        include_dynamic_group)
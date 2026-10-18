import ctypes
from pathlib import Path
from typing import Literal
import orjson

from pysoltrace import dot_h
from pysoltrace.api.utils import st_function
from pysoltrace.api.dll import context

##############################################################
# functions for simulation data management thru json strings #
##############################################################
class csv(context):
    @st_function
    def dump(self, filename: str, precision: int = 12) -> None:
        return self._pdll.st_write_results_csv(self._pcxt, filename.encode(), precision)

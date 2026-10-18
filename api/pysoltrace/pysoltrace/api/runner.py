import ctypes
from typing import Literal

from pysoltrace import dot_h
from pysoltrace.api.utils import st_function
from pysoltrace.api.dll import context, STATIC_DLL

############################################
# functions for SolTrace runner management #
############################################
class runner(context):
    @staticmethod
    @st_function
    def get_installed() -> dict[str, bool]:
        installed = ctypes.c_ubyte()
        code = STATIC_DLL.st_get_installed_runners(ctypes.byref(installed))
        return code, {
            dot_h.st_runner_type_t.NATIVE.name: bool(installed.value & (1 << dot_h.st_runner_type_t.NATIVE.value)),
            dot_h.st_runner_type_t.EMBREE.name: bool(installed.value & (1 << dot_h.st_runner_type_t.EMBREE.value)),
            dot_h.st_runner_type_t.OPTIX.name:  bool(installed.value & (1 << dot_h.st_runner_type_t.OPTIX.value)),
        }
    
    @staticmethod
    @st_function
    def is_installed(runner: int) -> bool:
        installed = ctypes.c_bool()
        code = STATIC_DLL.st_is_runner_installed(runner, ctypes.byref(installed))
        return code, installed.value

    @st_function
    def setup(self,
              runner_type: Literal[0, 1, 2],
              num_threads: int = 8,
              seeds:       list = None) -> None:
        num_seeds = 0
        # redefine seeds from list to C array
        _seeds = None
        if seeds and len(seeds):
            num_seeds = len(seeds)
            _seeds = (ctypes.c_uint * num_seeds)(*seeds)
        return self._pdll.st_sim_setup(self._pcxt,
                                      runner_type,
                                      num_threads,
                                      _seeds,
                                      num_seeds)

    @st_function
    def run(self) -> None:
        return self._pdll.st_sim_run_v2(self._pcxt)

    @st_function
    def report(self, level: int = 0) -> None:
        return self._pdll.st_sim_report(self._pcxt, level)
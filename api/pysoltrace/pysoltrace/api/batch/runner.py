import ctypes
from typing import Literal

from pysoltrace import dot_h
from pysoltrace.api.batch.utils import batcher, generate_api_call

# short names
_SET = dot_h.st_api_call.CALL_ST_SIM_SETUP
_RUN = dot_h.st_api_call.CALL_ST_SIM_RUN_V2
_REP = dot_h.st_api_call.CALL_ST_SIM_REPORT

############################################
# functions for SolTrace runner management #
############################################
class runner(batcher):
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
        return self.adder(generate_api_call(_SET,
                                            runner_type,
                                            num_threads,
                                            _seeds,
                                            num_seeds))

    def run(self) -> None:
        return self.adder(generate_api_call(_RUN))

    def report(self, level: int = 0) -> None:
        return self.adder(generate_api_call(_REP, level))
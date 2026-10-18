from pysoltrace import dot_h
from pysoltrace.api.batch.utils import batcher, generate_api_call

class legacy(batcher):
    def sim_params(self,
                   raycount: int,
                   maxcount: int,
                   include_dynamic_group: bool) -> None:
        return self.adder(generate_api_call(dot_h.st_api_call.CALL_ST_SIM_PARAMS,
                                            raycount,
                                            maxcount,
                                            include_dynamic_group))
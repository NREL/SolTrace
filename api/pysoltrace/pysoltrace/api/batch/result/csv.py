from pysoltrace import dot_h
from pysoltrace.api.batch.utils import batcher, generate_api_call

##############################################################
# functions for simulation data management thru json strings #
##############################################################
class csv(batcher):
    def dump(self, filename: str, precision: int = 12) -> None:
        return self.adder(generate_api_call(dot_h.st_api_call.CALL_ST_WRITE_RESULTS_CSV,
                                            filename.encode(), precision))

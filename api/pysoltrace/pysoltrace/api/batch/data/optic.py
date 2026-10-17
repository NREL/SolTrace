import ctypes

from pysoltrace import dot_h, soltrace_structs as _STS
from pysoltrace.api.batch.utils import batcher, generate_api_call

# short names
_NUM = dot_h.st_api_call.CALL_ST_NUM_OPTICS
_ADD = dot_h.st_api_call.CALL_ST_ADD_OPTICAL_PROPERIES_SET
_DEL = dot_h.st_api_call.CALL_ST_DELETE_OPTIC
_CL  = dot_h.st_api_call.CALL_ST_CLEAR_OPTICS

##################################################
# functions to add/remove/set optical properties #
##################################################
class optic(batcher):
    def num(self) -> int:
        pcount = ctypes.c_uint64()
        call = generate_api_call(_NUM, ctypes.pointer(pcount))
        return self.adder(call, lambda: pcount.value)

    def add(self,
            opt_set: _STS.optical_properties_set,
            front:   _STS.optical_properties_face,
            back:    _STS.optical_properties_face) -> int:
        opt_id = ctypes.c_uint64()
        call = generate_api_call(_ADD,
                                 ctypes.pointer(opt_set.ctype),
                                 ctypes.pointer(front.ctype),
                                 ctypes.pointer(back.ctype),
                                 ctypes.pointer(opt_id))
        return self.adder(call, lambda: opt_id.value)

    # TODO:
    # def get(self, optic_id: int) -> tuple[_STS.optical_properties_set, 
    #                                       _STS.optical_properties_face, 
    #                                       _STS.optical_properties_face]:
    #     opt_set = dot_h.args_optical_properties_set()
    #     front = dot_h.args_optical_properties_face()
    #     back = dot_h.args_optical_properties_face()
    #     code = self._pdll.st_get_optical_properties_set(self._pcxt,
    #                                                    optic_id,
    #                                                    ctypes.byref(opt_set),
    #                                                    ctypes.byref(front),
    #                                                    ctypes.byref(back))
    #     return code, opt_set.value, front.value, back.value

    def delete(self, idx: int) -> None:
        return self.adder(generate_api_call(_DEL, idx))

    def clear(self) -> None:
        return self.adder(generate_api_call(_CL))

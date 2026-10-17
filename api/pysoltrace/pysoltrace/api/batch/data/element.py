import ctypes

from pysoltrace import dot_h, soltrace_structs as _STS
from pysoltrace.api.batch.utils import batcher, generate_api_call
# from pysoltrace.api.utils import make_c_double_8

# short names
_NE = dot_h.st_api_call.CALL_ST_NUM_ELEMENTS
_AE = dot_h.st_api_call.CALL_ST_ADD_ELEMENT
_DE = dot_h.st_api_call.CALL_ST_DELETE_ELEMENT
_CE = dot_h.st_api_call.CALL_ST_CLEAR_ELEMENTS
_EL_EN = dot_h.st_api_call.CALL_ST_ELEMENT_ENABLED
_EL_VI = dot_h.st_api_call.CALL_ST_ELEMENT_VIRTUAL
_EL_XYZ = dot_h.st_api_call.CALL_ST_ELEMENT_XYZ
_EL_AIM = dot_h.st_api_call.CALL_ST_ELEMENT_AIM
_EL_Z = dot_h.st_api_call.CALL_ST_ELEMENT_ZROT
_EL_AP = dot_h.st_api_call.CALL_ST_ELEMENT_APERTURE
_EL_SU = dot_h.st_api_call.CALL_ST_ELEMENT_SURFACE
_EL_OP = dot_h.st_api_call.CALL_ST_ELEMENT_OPTIC
_EL_GR = dot_h.st_api_call.CALL_ST_ELEMENT_GROUP

###########################################
# functions to add/remove/modify elements #
###########################################
class element(batcher):
    def num(self) -> int:
        pcount = ctypes.c_uint64()
        call = generate_api_call(_NE,
                                 ctypes.pointer(pcount))
        return self.adder(call, lambda: pcount.value)

    def add(self,
            args: _STS.element,
            opt_id: int,
            geom: _STS.geometry) -> int:
        pid = ctypes.c_uint64()
        call = generate_api_call(_AE,
                                 ctypes.pointer(args.ctype),
                                 opt_id,
                                 ctypes.pointer(geom.ctype),
                                 ctypes.pointer(pid))
        return self.adder(call, lambda: pid.value)

    # TODO:
    # def get(self, id: int) -> tuple[_STS.element,
    #                                 int,
    #                                 list[float],
    #                                 list[float]]:
    #     args = dot_h.args_element()
    #     optic_id = ctypes.c_int64()
    #     a_params = make_c_double_8([0 for _ in range(8)])
    #     s_params = make_c_double_8([0 for _ in range(8)])
    #     code = self._pdll.st_get_element(self._pcxt,
    #                                     id,
    #                                     ctypes.byref(args),
    #                                     ctypes.byref(optic_id),
    #                                     ctypes.byref(a_params),
    #                                     ctypes.byref(s_params))
    #     return code, args.value, optic_id.value, a_params[:8], s_params[:8]

    def delete(self, idx: int) -> None:
        return self.adder(generate_api_call(_DE, idx))

    def clear(self) -> None:
        return self.adder(generate_api_call(_CE))
    
    ################################
    # functions to modify elements #
    ################################

    def enabled(self, idx: int, enabled_flag: bool) -> None:
        return self.adder(generate_api_call(_EL_EN, idx, enabled_flag))

    def virtual(self, idx: int, virtual_flag: bool) -> None:
        return self.adder(generate_api_call(_EL_VI, idx, virtual_flag))

    def xyz(self, idx: int, xyz: list[float]) -> None:
        assert len(xyz) == 3, f'Expected length 3, got {len(xyz)}'
        return self.adder(generate_api_call(_EL_XYZ, idx, *xyz))

    def aim(self, idx: int, aim: list[float]) -> None:
        assert len(aim) == 3, f'Expected length 3, got {len(aim)}'
        return self.adder(generate_api_call(_EL_AIM, idx, *aim))

    def zrot(self, idx: int, zrot: float) -> None:
        return self.adder(generate_api_call(_EL_Z, idx, zrot))

    def aperture(self,
                 idx: int,
                 ap: bytes,
                 params: list[float]) -> None:
        if isinstance(ap, str): ap = ap[0].encode()
        assert isinstance(ap, bytes) and len(ap) == 1, \
            "Aperture type must be a single character byte string."
        
        _params = _STS.make_c_double_8(params)
        return self.adder(generate_api_call(_EL_AP, idx, ap, _params))

    def surface(self,
                idx: int,
                surf: bytes,
                params: list[float]) -> None:
        if isinstance(surf, str): surf = surf[0].encode()
        assert isinstance(surf, bytes) and len(surf) == 1, \
            "Surface type must be a single character byte string."

        _params = _STS.make_c_double_8(params)
        return self.adder(generate_api_call(_EL_SU, idx, surf, _params))

    def optic(self, idx: int, opt_id: int) -> None:
        return self.adder(generate_api_call(_EL_OP, idx, opt_id))

    def group(self, idx: int, group: int) -> None:
        return self.adder(generate_api_call(_EL_GR, idx, group))
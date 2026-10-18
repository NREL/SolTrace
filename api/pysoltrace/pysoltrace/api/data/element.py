import ctypes

from pysoltrace import dot_h, soltrace_structs as _STS
from pysoltrace.api.utils import st_function
from pysoltrace.api.dll import context
# from pysoltrace.api.utils import make_c_double_8

###########################################
# functions to add/remove/modify elements #
###########################################
class element(context):
    @st_function
    def num(self) -> int:
        num_elements = ctypes.c_uint64()
        code = self._pdll.st_num_elements(self._pcxt, 
                                        ctypes.byref(num_elements))
        return code, num_elements.value

    def __len__(self):
        if not self._pdll or not self._pcxt:
            return 0
        return self.num()

    @st_function
    def add(self,
            args: _STS.element,
            opt_id: int,
            geom: _STS.geometry) -> int:
        pid = ctypes.c_uint64()
        code = self._pdll.st_add_element(self._pcxt,
                                         ctypes.byref(args.ctype),
                                         opt_id,
                                         ctypes.byref(geom.ctype),
                                         ctypes.byref(pid))
        return code, pid.value

    @st_function
    def get(self, id: int) -> tuple[_STS.element,
                                    int,
                                    list[float],
                                    list[float]]:
        args = dot_h.args_element()
        optic_id = ctypes.c_int64()
        geom = _STS.geometry().ctype
        code = self._pdll.st_get_element(self._pcxt,
                                         id,
                                         ctypes.byref(args),
                                         ctypes.byref(optic_id),
                                         ctypes.byref(geom))
        return code, args.value, optic_id.value, geom.value

    @st_function
    def delete(self, idx: int) -> None:
        return self._pdll.st_delete_element(self._pcxt, idx)

    @st_function
    def clear(self) -> None:
        return self._pdll.st_clear_elements(self._pcxt)
    
    ################################
    # functions to modify elements #
    ################################

    @st_function
    def enabled(self, idx: int, enabled_flag: bool) -> None:
        return self._pdll.st_element_enabled(self._pcxt, idx, enabled_flag)

    @st_function
    def virtual(self, idx: int, virtual_flag: bool) -> None:
        return self._pdll.st_element_virtual(self._pcxt, idx, virtual_flag)

    @st_function
    def xyz(self, idx: int, xyz: list[float]) -> None:
        assert len(xyz) == 3, f'Expected length 3, got {len(xyz)}'
        return self._pdll.st_element_xyz(self._pcxt, idx, *xyz)

    @st_function
    def aim(self, idx: int, aim: list[float]) -> None:
        assert len(aim) == 3, f'Expected length 3, got {len(aim)}'
        return self._pdll.st_element_aim(self._pcxt, idx, *aim)

    @st_function
    def zrot(self, idx: int, zrot: float) -> None:
        return self._pdll.st_element_zrot(self._pcxt, idx, zrot)

    @st_function
    def aperture(self,
                 idx: int,
                 ap: bytes,
                 params: list[float]) -> None:
        if isinstance(ap, str): ap = ap[0].encode()
        assert isinstance(ap, bytes) and len(ap) == 1, \
            "Aperture type must be a single character byte string."
        
        _params = _STS.make_c_double_8(params)
        return self._pdll.st_element_aperture(self._pcxt, idx, ap, _params)

    @st_function
    def surface(self,
                idx: int,
                surf: bytes,
                params: list[float]) -> None:
        if isinstance(surf, str): surf = surf[0].encode()
        assert isinstance(surf, bytes) and len(surf) == 1, \
            "Surface type must be a single character byte string."

        _params = _STS.make_c_double_8(params)
        return self._pdll.st_element_surface(self._pcxt, idx, surf, _params)

    @st_function
    def optic(self, idx: int, opt_id: int) -> None:
        return self._pdll.st_element_optic(self._pcxt, idx, opt_id)

    @st_function
    def group(self, idx: int, group: int) -> None:
        return self._pdll.st_element_group(self._pcxt, idx, group)
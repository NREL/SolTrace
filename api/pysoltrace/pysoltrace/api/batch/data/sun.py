import ctypes

from pysoltrace import dot_h, Point, soltrace_constants as _STC
from pysoltrace.api.batch.utils import batcher, generate_api_call

# short names
_ADD = dot_h.st_api_call.CALL_ST_ADD_SUN
_SHA = dot_h.st_api_call.CALL_ST_SUN_SHAPE
_XYZ = dot_h.st_api_call.CALL_ST_SUN_XYZ
_POS = dot_h.st_api_call.CALL_ST_SUN_POSITION
_UD  = dot_h.st_api_call.CALL_ST_SUN_USERDATA
_AZ  = dot_h.st_api_call.CALL_ST_GET_SUN_AZ_ZEN
_AE  = dot_h.st_api_call.CALL_ST_GET_SUN_AZ_EL
_VEC = dot_h.st_api_call.CALL_ST_GET_SUN_VECTOR

###################################
# functions to add/modify the sun #
###################################
class sun(batcher):
    def add(self,
            args: _STC.sun,
            angle: list[float] = [],
            intensity: list[float] = []) -> None:
        args.npoints = len(angle)
        assert len(angle) == args.npoints and \
               len(intensity) == args.npoints, \
            f'Expected {args.npoints} angle and intensity values,' \
            f'got {len(angle)} and {len(intensity)}'
        _angle     = (ctypes.c_double * len(angle))(*angle)
        _intensity = (ctypes.c_double * len(angle))(*intensity)
        return self.adder(generate_api_call(_ADD,
                                            ctypes.pointer(args.ctype),
                                            _angle,
                                            _intensity))

    # TODO:
    # def get(self) -> tuple[_STC.sun, list[float], list[float]]:
    #     args = dot_h.args_sun()
    #     angle = ctypes.c_double()
    #     intensity = ctypes.c_double()
    #     code = self._pdll.st_get_sun(self._pcxt,
    #                                 ctypes.byref(args),
    #                                 ctypes.byref(ctypes.pointer(angle)),
    #                                 ctypes.byref(ctypes.pointer(intensity)))
    #     return code, args.value, angle, intensity

    def shape(self,
              shape: bytes,
              sigma_halfwidth_csr: float) -> None:
        if isinstance(shape, str): shape = shape[0].encode()
        assert isinstance(shape, bytes) and len(shape) == 1, \
            "Sun shape type must be a single character byte string."
        return self.adder(generate_api_call(_SHA,
                                            shape,
                                            sigma_halfwidth_csr))

    def xyz(self, x: float, y: float, z: float) -> None:
        return self.adder(generate_api_call(_XYZ, x, y, z))

    def position(self, lat: float, day: float, hour: float) -> Point:
        px = ctypes.c_double()
        py = ctypes.c_double()
        pz = ctypes.c_double()
        call = generate_api_call(_POS,
                                 lat, day, hour,
                                 ctypes.pointer(px),
                                 ctypes.pointer(py),
                                 ctypes.pointer(pz))
        return self.adder(call, lambda: Point(px.value, py.value, pz.value))

    def userdata(self,
                 npoints:   int,
                 angle:     list[float],
                 intensity: list[float]) -> None:
        assert len(angle) == npoints and len(intensity) == npoints, \
            f'Expected {npoints} angle and intensity values,' \
            f'got {len(angle)} and {len(intensity)}'
        _angle     = (ctypes.c_double * npoints)(*angle)
        _intensity = (ctypes.c_double * npoints)(*intensity)
        return self.adder(generate_api_call(_UD,
                                         npoints,
                                         _angle,
                                         _intensity))

    def az_zen(self,
               calc: int,
               loc:  _STC.sun_location,
               dt:   _STC.sun_datetime) -> tuple[float, float]:
        az  = ctypes.c_double()
        zen = ctypes.c_double()
        call = generate_api_call(_AZ,
                                 calc,
                                 ctypes.pointer(loc.ctype),
                                 ctypes.pointer(dt.ctype),
                                 ctypes.pointer(az),
                                 ctypes.pointer(zen))
        return self.adder(call, lambda: (az.value, zen.value))

    def az_el(self,
              calc: int,
              loc:  _STC.sun_location,
              dt:   _STC.sun_datetime) -> tuple[float, float]:
        az = ctypes.c_double()
        el = ctypes.c_double()
        call = generate_api_call(_AE,
                                 calc,
                                 ctypes.pointer(loc.ctype),
                                 ctypes.pointer(dt.ctype),
                                 ctypes.pointer(az),
                                 ctypes.pointer(el))
        return self.adder(call, lambda: (az.value, el.value))

    def vector(self,
               calc: int,
               loc:  _STC.sun_location,
               dt:   _STC.sun_datetime) -> Point:
        sun_x = ctypes.c_double()
        sun_y = ctypes.c_double()
        sun_z = ctypes.c_double()
        call = generate_api_call(_VEC,
                                 calc,
                                 ctypes.pointer(loc.ctype),
                                 ctypes.pointer(dt.ctype),
                                 ctypes.pointer(sun_x),
                                 ctypes.pointer(sun_y),
                                 ctypes.pointer(sun_z))
        return self.adder(call, lambda: Point(sun_x.value,
                                              sun_y.value,
                                              sun_z.value))


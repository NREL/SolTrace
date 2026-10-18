import ctypes

from pysoltrace import dot_h, Point, soltrace_constants as _STC
from pysoltrace.api.utils import st_function
from pysoltrace.api.dll import context

###################################
# functions to add/modify the sun #
###################################
class sun(context):
    @st_function
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
        return self._pdll.st_add_sun(self._pcxt,
                                    ctypes.byref(args.ctype),
                                    _angle,
                                    _intensity)

    @st_function
    def get(self, npoints: int = 0) -> tuple[_STC.sun, list[float], list[float]]:
        args = dot_h.args_sun()
        expect_userdata = npoints > 0
        angle     = ctypes.cast((ctypes.c_double * npoints)(), ctypes.POINTER(ctypes.c_double)) \
                    if expect_userdata else ctypes.pointer(ctypes.c_double())
        intensity = ctypes.cast((ctypes.c_double * npoints)(), ctypes.POINTER(ctypes.c_double)) \
                    if expect_userdata else ctypes.pointer(ctypes.c_double())
        code = self._pdll.st_get_sun(self._pcxt,
                                    ctypes.byref(args),
                                    ctypes.byref(angle),
                                    ctypes.byref(intensity),
                                    npoints)
        rt = (args.value, angle[:npoints], intensity[:npoints]) if expect_userdata else (args.value, [], [])
        return code, *rt

    @st_function
    def shape(self,
              shape: bytes,
              sigma_halfwidth_csr: float) -> None:
        if isinstance(shape, str): shape = shape[0].encode()
        assert isinstance(shape, bytes) and len(shape) == 1, \
            "Sun shape type must be a single character byte string."
        return self._pdll.st_sun_shape(self._pcxt,
                                      shape,
                                      sigma_halfwidth_csr)

    @st_function
    def xyz(self, x: float, y: float, z: float) -> None:
        return self._pdll.st_sun_xyz(self._pcxt, x, y, z)

    @st_function
    def position(self, lat: float, day: float, hour: float) -> Point:
        px = ctypes.c_double()
        py = ctypes.c_double()
        pz = ctypes.c_double()
        code = self._pdll.st_sun_position(self._pcxt,
                                          lat, day, hour,
                                          ctypes.byref(px),
                                          ctypes.byref(py),
                                          ctypes.byref(pz))
        return code, Point(px.value, py.value, pz.value)

    @st_function
    def userdata(self,
                 npoints:   int,
                 angle:     list[float],
                 intensity: list[float]) -> None:
        assert len(angle) == npoints and len(intensity) == npoints, \
            f'Expected {npoints} angle and intensity values,' \
            f'got {len(angle)} and {len(intensity)}'
        _angle     = (ctypes.c_double * npoints)(*angle)
        _intensity = (ctypes.c_double * npoints)(*intensity)
        return self._pdll.st_sun_userdata(self._pcxt,
                                         npoints,
                                         _angle,
                                         _intensity)

    @st_function
    def az_zen(self,
               calc: int,
               loc:  _STC.sun_location,
               dt:   _STC.sun_datetime) -> tuple[float, float]:
        az  = ctypes.c_double()
        zen = ctypes.c_double()
        code = self._pdll.st_get_sun_az_zen(self._pcxt,
                                           calc,
                                           ctypes.byref(loc.ctype),
                                           ctypes.byref(dt.ctype),
                                           ctypes.byref(az),
                                           ctypes.byref(zen))
        return code, az.value, zen.value

    @st_function
    def az_el(self,
              calc: int,
              loc:  _STC.sun_location,
              dt:   _STC.sun_datetime) -> tuple[float, float]:
        az = ctypes.c_double()
        el = ctypes.c_double()
        code = self._pdll.st_get_sun_az_el(self._pcxt,
                                          calc,
                                          ctypes.byref(loc.ctype),
                                          ctypes.byref(dt.ctype),
                                          ctypes.byref(az),
                                          ctypes.byref(el))
        return code, az.value, el.value

    @st_function
    def vector(self,
               calc: int,
               loc:  _STC.sun_location,
               dt:   _STC.sun_datetime) -> Point:
        sun_x = ctypes.c_double()
        sun_y = ctypes.c_double()
        sun_z = ctypes.c_double()
        code = self._pdll.st_get_sun_vector(self._pcxt,
                                           calc,
                                           ctypes.byref(loc.ctype),
                                           ctypes.byref(dt.ctype),
                                           ctypes.byref(sun_x),
                                           ctypes.byref(sun_y),
                                           ctypes.byref(sun_z))
        return code, Point(sun_x.value, sun_y.value, sun_z.value)

import ctypes
from dataclasses import dataclass, field
from pysoltrace import dot_h

def make_c_double_8(params: list[float]) -> ctypes.Array:
    assert isinstance(params, list), f'Expected a list, got {type(params)}'
    assert len(params) <= 8,         f'Expected at most 8 elements, got {len(params)}'
    return (ctypes.c_double * 8)(*params)

# reexport structs in constants so it can be used in type hinting
class base_args:
    @property
    def ctype(self): return NotImplemented

@dataclass
class simulation_parameters(base_args):
    latitude:                 float
    longitude:                float
    number_of_rays:           int   = 1_000_000
    max_number_of_rays:       int   = 100_000_000
    tolerance:                float = .01
    include_sun_shape_errors: bool  = True
    include_optical_errors:   bool  = True
    as_power_tower:           bool  = False

    @property
    def ctype(self):
        return dot_h.args_simulation_parameters(self.number_of_rays,
                                                self.max_number_of_rays,
                                                self.tolerance,
                                                self.latitude,
                                                self.longitude,
                                                self.include_sun_shape_errors,
                                                self.include_optical_errors,
                                                self.as_power_tower)
@dataclass
class optical_properties_face(base_args):
    transmissivity:          float
    reflectivity:            float
    slope_error:             float
    specularity_error:       float
    error_distribution_type: str

    @property
    def ctype(self):
        return dot_h.args_optical_properties_face(self.transmissivity,
                                                 self.reflectivity,
                                                 self.slope_error,
                                                 self.specularity_error,
                                                 self.error_distribution_type)

@dataclass
class optical_properties_set(base_args):
    name:                   str
    refraction_index_front: float
    refraction_index_back:  float
    type:                   int

    @property
    def ctype(self):
        return dot_h.args_optical_properties_set(self.name,
                                                  self.refraction_index_front,
                                                  self.refraction_index_back,
                                                  self.type)
    
@dataclass
class element(base_args):
    xyz:          list[float]
    aim:          list[float]
    zrot:         float
    enabled_flag: bool
    virtual_flag: bool
    group:        int = -1

    @property
    def ctype(self):
        assert len(self.xyz) == 3, f'Expected length 3, got {len(self.xyz)}'
        assert len(self.aim) == 3, f'Expected length 3, got {len(self.aim)}'
        return dot_h.args_element(*self.xyz,
                                  *self.aim,
                                  self.zrot,
                                  self.enabled_flag,
                                  self.virtual_flag,
                                #   self.ap,
                                #   self.surf,
                                  self.group)

@dataclass
class geometry(base_args):
    a_params: list[float] = field(default_factory = lambda: [0 for _ in range(8)])
    s_params: list[float] = field(default_factory = lambda: [0 for _ in range(8)])
    ap:       bytes = b' '
    surf:     bytes = b' '

    @property
    def ctype(self):
        return dot_h.args_geometry(make_c_double_8(self.a_params),
                                   make_c_double_8(self.s_params),
                                   self.ap,
                                   self.surf)

@dataclass
class sun(base_args):
    npoints:             int    # TODO: make this default 0
    x:                   float
    y:                   float
    z:                   float
    sigma_halfwidth_csr: float
    shape:               bytes

    @property
    def ctype(self):
        return dot_h.args_sun(self.npoints,
                              self.x,
                              self.y,
                              self.z,
                              self.sigma_halfwidth_csr,
                              self.shape)

@dataclass
class sun_location(base_args):
    latitude:  float
    longitude: float
    timeZone:  float
    altitude:  float = 0

    @property
    def ctype(self):
        return dot_h.args_sun_location(self.latitude,
                                       self.longitude,
                                       self.timeZone,
                                       self.altitude)

# TODO: util to convert built-in datetime to st_datetime
@dataclass
class sun_datetime(base_args):
    year:   int
    month:  int
    day:    int
    hour:   int = 12
    minute: int = 0
    second: int = 0

    @property
    def ctype(self):
        return dot_h.args_sun_datetime(self.year,
                                       self.month,
                                       self.day,
                                       self.hour,
                                       self.minute,
                                       self.second)
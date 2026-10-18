"""
your home for creating and checking the compatability your soltrace json files

TODO: a checker for save_to_json
TODO: chedder, maybe:
        struct json_keys {
            static constexpr const char is_stage[] = "is_stage";
            static constexpr const char is_composite[] = "is_composite";
            static constexpr const char active[] = "active";
        };
        would need to check regex for parsing static/constexpr
        and a way to remap the ctypes as a type
"""

#################
# keys for json #
#################
K_IS_STAGE                 = 'is_stage'
K_IS_COMPOSITE             = 'is_composite'
K_ACTIVE                   = 'active'
K_VIRTUAL_FLAG             = 'virtual_flag'
K_IS_SINGLE                = 'is_single'
K_MY_NAME                  = 'my_name'
K_STAGE                    = 'stage'
K_ID                       = 'id'
K_POSITION                 = 'position'
K_AIM                      = 'aim'
K_ZROT                     = 'zrot'
K_IS_VIRTUAL               = 'is_virtual'
K_IS_MULTIHIT              = 'is_multihit'
K_IS_TRACETHROUGH          = 'is_tracethrough'
K_ORIGIN                   = 'origin'
K_INCLUDE_SUN_SHAPE_ERRORS = 'include_sun_shape_errors'
K_INCLUDE_OPTICAL_ERRORS   = 'include_optical_errors'
K_NUMBER_OF_RAYS           = 'number_of_rays'
K_MAX_NUMBER_OF_RAYS       = 'max_number_of_rays'
K_TOLERANCE                = 'tolerance'
K_LATITUDE                 = 'latitude'
K_LONGITUDE                = 'longitude'
K_SEED                     = 'seed'

# ray source keys
K_SOURCE_TYPE    = 'source_type'
K_MY_SHAPE       = 'my_shape'
K_SIGMA          = 'sigma'
K_HALF_WIDTH     = 'half_width'
K_CSR            = 'csr'
K_USER_ANGLE     = 'user_angle'
K_USER_INTENSITY = 'user_intensity'
K_GEN_TYPE       = 'gen_type'
K_POS            = 'pos'

# aperture keys
K_APERTURE_TYPE = 'aperture_type'
K_X_LENGTH = 'x_length'
K_Y_LENGTH = 'y_length'
K_X_COORD = 'x_coord'
K_Y_COORD = 'y_coord'
K_X1 = 'x1'
K_Y1 = 'y1'
K_X2 = 'x2'
K_Y2 = 'y2'
K_X3 = 'x3'
K_Y3 = 'y3'
K_X4 = 'x4'
K_Y4 = 'y4'

# surface keys
K_SURFACE_TYPE   = 'surface_type'
K_FOCAL_LENGTH_X = 'focal_length_x'
K_FOCAL_LENGTH_Y = 'focal_length_y'

# optical properties keys
K_ERROR_DISTRIBUTION_TYPE = 'error_distribution_type'
K_TRANSMISSIVITY          = 'transmissivity'
K_REFLECTIVITY            = 'reflectivity'
K_SLOPE_ERROR             = 'slope_error'
K_SPECULARITY_ERROR       = 'specularity_error'
K_MY_TYPE                 = 'my_type'
K_REFRACTION_INDEX_FRONT  = 'refraction_index_front'
K_REFRACTION_INDEX_BACK   = 'refraction_index_back'
K_FRONT                   = 'front'
K_BACK                    = 'back'

# element keys
K_MY_ID    = 'my_id'
K_GROUP    = 'group'
K_APERTURE = 'aperture'
K_SURFACE  = 'surface'
K_OPT_ID   = 'opt_id'

#############
# constants #
#############
enumify = lambda arr: { k: i for i, k in enumerate(arr) }
"""make a list of things an enum (ish)"""

"""sun shape"""
SUN_NONE         = 'NONE'
SUN_GAUSSIAN     = 'GAUSSIAN'
SUN_PILLBOX      = 'PILLBOX'
SUN_LIMBDARKENED = 'LIMBDARKENED'
SUN_BUIE_CSR     = 'BUIE_CSR'
SUN_USER_DEFINED = 'USER_DEFINED'
SUN_UNKNOWN      = 'UNKNOWN'

SUN_SHAPE_ENUM = enumify([
    SUN_NONE, SUN_GAUSSIAN, SUN_PILLBOX, SUN_LIMBDARKENED, SUN_BUIE_CSR, SUN_USER_DEFINED, SUN_UNKNOWN
])

RANDOM          = 'RANDOM'
HALTON          = 'HALTON'
RAY_GEN_UNKNOWN = 'UNKNOWN'

RAY_GENERATION_ENUM = enumify([RANDOM, HALTON, RAY_GEN_UNKNOWN])
RAY_GENERATION_ASSERT = [RAY_GENERATION_ENUM[RANDOM], RAY_GENERATION_ENUM[HALTON]]


"""aperture types"""
ANNULUS                       = 'ANNULUS'
CIRCLE                        = 'CIRCLE'
HEXAGON                       = 'HEXAGON'
RECTANGLE                     = 'RECTANGLE'
EQUILATERAL_TRIANGLE          = 'EQUILATERAL_TRIANGLE'
SINGLE_AXIS_CURVATURE_SECTION = 'SINGLE_AXIS_CURVATURE_SECTION'
IRREGULAR_TRIANGLE            = 'IRREGULAR_TRIANGLE'
IRREGULAR_QUADRILATERAL       = 'IRREGULAR_QUADRILATERAL'
APERTURE_UNKNOWN              = 'APERTURE_UNKNOWN'

APERTURE_ENUM = enumify([
    ANNULUS, CIRCLE, HEXAGON, RECTANGLE, 
    EQUILATERAL_TRIANGLE, SINGLE_AXIS_CURVATURE_SECTION, 
    IRREGULAR_TRIANGLE, IRREGULAR_QUADRILATERAL, APERTURE_UNKNOWN
])

APERTURE_JSON = {
    RECTANGLE: lambda x_l, y_l: {
        K_APERTURE_TYPE: RECTANGLE,
        K_X_LENGTH:      x_l,
        K_Y_LENGTH:      y_l,
        K_X_COORD:       -0.5 * x_l,
        K_Y_COORD:       -0.5 * y_l,
    },
    IRREGULAR_QUADRILATERAL: lambda x1, y1, x2, y2, x3, y3, x4, y4: {
        K_APERTURE_TYPE: IRREGULAR_QUADRILATERAL,
        K_X1: x1, K_Y1: y1, K_X2: x2, K_Y2: y2, K_X3: x3, K_Y3: y3, K_X4: x4, K_Y4: y4
    }
}

"""surface types"""
CONE                  = 'CONE'
CYLINDER              = 'CYLINDER'
FLAT                  = 'FLAT'
PARABOLA              = 'PARABOLA'
SPHERE                = 'SPHERE'
HYPER                 = 'HYPER'
GENERAL_SPENCER_MURTY = 'GENERAL_SPENCER_MURTY'
TORUS                 = 'TORUS'
SURFACE_UNKNOWN       = 'SURFACE_UNKNOWN'

SURFACE_ENUM = enumify([
    CONE, CYLINDER, FLAT, PARABOLA, SPHERE, 
    GENERAL_SPENCER_MURTY, HYPER, TORUS, SURFACE_UNKNOWN
])

SURFACE_JSON = {
    PARABOLA: lambda *args: {
            K_SURFACE_TYPE:   PARABOLA,
            K_FOCAL_LENGTH_X: args[0],
            K_FOCAL_LENGTH_Y: args[1],
        },
    FLAT: lambda *_: { K_SURFACE_TYPE: FLAT }
}

"""optical interaction types"""
REFLECTION          = 'REFLECTION'
REFRACTION          = 'REFRACTION'
INTERACTION_UNKNOWN = 'UNKNOWN'

OPTICAL_INTERACTIONS_ENUM = enumify([REFLECTION, REFRACTION, INTERACTION_UNKNOWN])

OPTICAL_INTERACTION_ASSERT = [
    OPTICAL_INTERACTIONS_ENUM[REFLECTION],
    OPTICAL_INTERACTIONS_ENUM[REFRACTION],
]

"""optical error types"""
ERROR_NONE    = 'NONE'
ERROR_GAUSSIAN      = 'GAUSSIAN'
ERROR_PILLBOX       = 'PILLBOX'
ERROR_DIFFUSE       = 'DIFFUSE'
ERROR_USER_DEFINED  = 'USER_DEFINED'
ERROR_UNKNOWN = 'UNKNOWN'

OPTICAL_ERRORS_ENUM = enumify([
    ERROR_NONE, ERROR_GAUSSIAN, ERROR_PILLBOX, ERROR_DIFFUSE, ERROR_USER_DEFINED, ERROR_UNKNOWN
])

OPTICAL_ERRORS_ASSERT = [
    OPTICAL_ERRORS_ENUM[ERROR_NONE],
    OPTICAL_ERRORS_ENUM[ERROR_GAUSSIAN],
    OPTICAL_ERRORS_ENUM[ERROR_PILLBOX],
    OPTICAL_ERRORS_ENUM[ERROR_DIFFUSE],
    OPTICAL_ERRORS_ENUM[ERROR_USER_DEFINED],
]
# end constants

##################
# default values #
##################
SUN_DEFAULT_USER_ANGLE = [
    0.0, 0.15, 0.3, 0.45, 0.6, 0.75, 0.9, 1.05, 1.2, 1.35, 1.5,
    1.65, 1.8, 1.95, 2.1, 2.25, 2.4, 2.55, 2.7, 2.85, 3.0, 3.15,
    3.3, 3.45, 3.6, 3.75, 3.9, 4.05, 4.2, 4.35, 4.5, 4.65, 4.8, 
    4.95, 5.1, 5.25, 5.4, 5.55, 5.7, 5.85, 6.0, 6.15, 6.3, 6.45, 
    6.6, 6.75, 6.9, 7.05, 7.2, 7.35, 7.5, 7.65, 7.8, 7.95,
]

SUN_DEFAULT_USER_INTENSITY = [
    1, 0.999872, 0.999485, 0.998837, 0.997923, 0.996734, 0.99526, 
    0.993487, 0.991399, 0.988976, 0.986193, 0.983019, 0.979417, 
    0.975345, 0.970747, 0.965558, 0.959697, 0.953063, 0.945528, 
    0.936933, 0.927069, 0.915665, 0.902358, 0.886653, 0.867855, 
    0.844965, 0.816477, 0.78003, 0.731687, 0.66436, 0.563875, 
    0.397159, 5.34414e-05, 5.07222e-05, 4.82164e-05, 4.59018e-05, 
    4.37589e-05, 4.17708e-05, 3.99224e-05, 3.82007e-05, 3.65941e-05, 
    3.50923e-05, 3.36861e-05, 3.23674e-05, 3.11289e-05, 2.9964e-05, 
    2.88669e-05, 2.78323e-05, 2.68554e-05, 2.59319e-05, 2.50579e-05, 
    2.42298e-05, 2.34443e-05, 2.26985e-05
]

DEFAULTS = {
    # simulation parameter defaults
    K_INCLUDE_SUN_SHAPE_ERRORS: True,
    K_INCLUDE_OPTICAL_ERRORS:   True,
    K_MAX_NUMBER_OF_RAYS:       100000,
    K_MAX_NUMBER_OF_RAYS:       1000000000,
    K_TOLERANCE:                1,
    K_SEED:                     608,

    # ray source defaults
    K_SOURCE_TYPE:    'Sun',
    K_SIGMA:          4.65,
    K_HALF_WIDTH:     4.65,
    K_CSR:            .05,
    K_USER_ANGLE:     [],
    K_USER_INTENSITY: [],
    K_GEN_TYPE:       HALTON,
    K_POS:            [.47239, -826.21038, 563.36151],

    # optical property defaults
    K_ERROR_DISTRIBUTION_TYPE: ERROR_GAUSSIAN,
    K_TRANSMISSIVITY:          0.0,
    K_REFLECTIVITY:            .95,
    K_SLOPE_ERROR:             1.0,
    K_SPECULARITY_ERROR:       0.1,
    K_MY_TYPE:                 REFLECTION,
    K_REFRACTION_INDEX_FRONT:  1.1,
    K_REFRACTION_INDEX_BACK:   1.1,

    # element defaults
    K_IS_SINGLE:    True,
    K_ACTIVE:       True,
    K_VIRTUAL_FLAG: False,
    K_MY_ID:        None,
    K_STAGE:        0,
}

def default_fallback(value, key):
    return value if value != None else DEFAULTS[key]

"""
removing staging from soltrace, so create constant global stage 
with +x, +y, and +z pointing east, north, and zenith, respectively.
"""
STAGE = {
    K_IS_STAGE:        True,
    K_IS_COMPOSITE:    True,
    K_ACTIVE:          True,
    K_VIRTUAL_FLAG:    False,
    K_IS_SINGLE:       False,
    K_MY_NAME:         'stage',
    K_STAGE:           0,
    K_ID:              0,
    K_POSITION:        [0, 0, 0],
    K_AIM:             [0, 0, 1],
    K_ZROT:            0,
    K_IS_VIRTUAL:      False,
    K_IS_MULTIHIT:     False,
    K_IS_TRACETHROUGH: False,
    K_ORIGIN:          [0, 0, 0],
}

def simulation_parameters_json(lat: float, long: float, **kwargs) -> dict:
    assert type(lat) == float, 'latitude must be a float'
    assert type(long) == float, 'longitude must be a float'

    for k, v in kwargs.items():
        if k in [K_INCLUDE_SUN_SHAPE_ERRORS, K_INCLUDE_OPTICAL_ERRORS]:
            assert type(v) == bool, f'{k} must be a boolean'
        elif k in [K_MAX_NUMBER_OF_RAYS, K_MAX_NUMBER_OF_RAYS, K_SEED]:
            assert type(v) == int, f'{k} must be an integer'
            assert v > 0, f'{k} must be positive'
        elif k == K_TOLERANCE:
            assert type(v) == float, 'tolerance must be a float'
            assert v > 0, 'tolerance must be positive'

    return {
        K_INCLUDE_SUN_SHAPE_ERRORS: DEFAULTS[K_INCLUDE_SUN_SHAPE_ERRORS],
        K_INCLUDE_OPTICAL_ERRORS:   DEFAULTS[K_INCLUDE_OPTICAL_ERRORS],
        K_MAX_NUMBER_OF_RAYS:       DEFAULTS[K_MAX_NUMBER_OF_RAYS],
        K_MAX_NUMBER_OF_RAYS:       DEFAULTS[K_MAX_NUMBER_OF_RAYS],
        K_TOLERANCE:                DEFAULTS[K_TOLERANCE],
        K_LATITUDE:                 lat,
        K_LONGITUDE:                long,
        K_SEED:                     DEFAULTS[K_SEED],
        **kwargs
    }

def _gaussian(_sigma: float = None, _gen_type: str = None) -> dict:
    _sigma    = default_fallback(_sigma, K_SIGMA)
    _gen_type = default_fallback(_gen_type, K_GEN_TYPE)

    assert type(_sigma) == float, 'Gaussian sigma must be float'
    assert _sigma >= 0,           f'Gaussian sigma ({_sigma}) must positive'
    assert RAY_GENERATION_ENUM[_gen_type] in RAY_GENERATION_ASSERT, f'invalid ray generation type ({_gen_type})'
    
    return {
        K_SOURCE_TYPE:    DEFAULTS[K_SOURCE_TYPE],
        K_MY_SHAPE:       SUN_GAUSSIAN,
        K_SIGMA:          _sigma,
        K_HALF_WIDTH:     DEFAULTS[K_HALF_WIDTH],
        K_CSR:            DEFAULTS[K_CSR], 
        K_USER_ANGLE:     DEFAULTS[K_USER_ANGLE],
        K_USER_INTENSITY: DEFAULTS[K_USER_INTENSITY],
        K_GEN_TYPE:       _gen_type,
        K_POS:            DEFAULTS[K_POS],
    }

def _pillbox(_half_width: float = None, _gen_type: str = None) -> dict:
    _half_width = default_fallback(_half_width, K_HALF_WIDTH)
    _gen_type   = default_fallback(_gen_type, K_GEN_TYPE)

    assert type(_half_width) == float, 'Pillbox half width must be float'
    assert _half_width >= 0,           f'Pillbox half width({_half_width}) must positive'
    assert RAY_GENERATION_ENUM[_gen_type] in RAY_GENERATION_ASSERT, f'invalid ray generation type ({_gen_type})'

    return {
        K_SOURCE_TYPE:    DEFAULTS[K_SOURCE_TYPE],
        K_MY_SHAPE:       SUN_PILLBOX,
        K_SIGMA:          4.65,
        K_HALF_WIDTH:     _half_width,
        K_CSR:            DEFAULTS[K_CSR], 
        K_USER_ANGLE:     DEFAULTS[K_USER_ANGLE],
        K_USER_INTENSITY: DEFAULTS[K_USER_INTENSITY],
        K_GEN_TYPE:       _gen_type,
        K_POS:            DEFAULTS[K_POS],
    }

def _buie_csr(_csr: float = None, _gen_type: str = None) -> dict:
    # TODO: come back and make sure default csr is a good number
    _csr      = default_fallback(_csr, K_CSR)
    _gen_type = default_fallback(_gen_type, K_GEN_TYPE)

    assert type(_csr) == float,      'Buie CSR must be float'
    assert _csr >= 0 and _csr <= .8, f'Buie CSR ({_csr}) must in range [0.0, 0.8]'
    assert RAY_GENERATION_ENUM[_gen_type] in RAY_GENERATION_ASSERT, f'invalid ray generation type ({_gen_type})'

    return {
        K_SOURCE_TYPE:    DEFAULTS[K_SOURCE_TYPE],
        K_MY_SHAPE:       SUN_BUIE_CSR,
        K_SIGMA:          4.65,
        K_HALF_WIDTH:     DEFAULTS[K_HALF_WIDTH],
        K_CSR:            _csr, 
        K_USER_ANGLE:     DEFAULTS[K_USER_ANGLE],
        K_USER_INTENSITY: DEFAULTS[K_USER_INTENSITY],
        K_GEN_TYPE:       _gen_type,
        K_POS:            DEFAULTS[K_POS],
    }

def _user_defined(_ua: list = None, _ui: list = None, _gen_type: str = None) -> dict:
    _ua       = _ua if _ua != None else SUN_DEFAULT_USER_ANGLE
    _ui       = _ui if _ui != None else SUN_DEFAULT_USER_INTENSITY
    _gen_type = default_fallback(_gen_type, K_GEN_TYPE)

    assert type(_ua) == list,    'user defined angles must be in a list'
    assert type(_ui) == list,    'user defined intensities must be in a list'
    assert len(_ua) == len(_ui), 'user defined angles and intensities must be same length'
    assert RAY_GENERATION_ENUM[_gen_type] in RAY_GENERATION_ASSERT, f'invalid ray generation type ({_gen_type})'

    return {
        K_SOURCE_TYPE:    DEFAULTS[K_SOURCE_TYPE],
        K_MY_SHAPE:       SUN_USER_DEFINED,
        K_SIGMA:          4.65,
        K_HALF_WIDTH:     DEFAULTS[K_HALF_WIDTH],
        K_CSR:            DEFAULTS[K_CSR], 
        K_USER_ANGLE:     _ua,
        K_USER_INTENSITY: _ui,
        K_GEN_TYPE:       _gen_type,
        K_POS:            DEFAULTS[K_POS],
    }

SUN_SHAPE_JSON = {
    SUN_GAUSSIAN: _gaussian,
    SUN_PILLBOX: _pillbox,
    SUN_BUIE_CSR: _buie_csr,
    SUN_USER_DEFINED: _user_defined,
}

def optical_side_json(_error_distribution_type: str   = None,
                      _transmissivity:          float = None,
                      _reflectivity:            float = None,
                      _slope_error:             float = None,
                      _specularity_error:       float = None) -> dict:
    # if argument is given set as the arg's value, otherwise default to heliostat optics
    _error_distribution_type = default_fallback(_error_distribution_type, K_ERROR_DISTRIBUTION_TYPE)
    _transmissivity          = default_fallback(_transmissivity, K_TRANSMISSIVITY)
    _reflectivity            = default_fallback(_reflectivity, K_REFLECTIVITY)
    _slope_error             = default_fallback(_slope_error, K_SLOPE_ERROR)
    _specularity_error       = default_fallback(_specularity_error, K_SPECULARITY_ERROR)

    # assertion checks
    assert OPTICAL_ERRORS_ENUM[_error_distribution_type] in OPTICAL_ERRORS_ASSERT, f'invalid optical error type {_error_distribution_type}'
    assert _transmissivity >= 0.0 and _transmissivity <= 1.0, f'transmissivity ({_transmissivity}) is out of bounds'
    assert _reflectivity >= 0.0 and _reflectivity <= 1.0,     f'reflectivity ({_reflectivity}) is out of bounds'
    assert type(_slope_error) == float,                       'slope error must be a float'
    assert type(_specularity_error) == float,                 'specularity error must be a float'

    return {
        K_ERROR_DISTRIBUTION_TYPE: _error_distribution_type,
        K_TRANSMISSIVITY:          _transmissivity,
        K_REFLECTIVITY:            _reflectivity,
        K_SLOPE_ERROR:             _slope_error,
        K_SPECULARITY_ERROR:       _specularity_error,
    }

class OpticalPropertyRegistry:
    def __init__(self):
        self.data = []
        self._current_item = 0

    def __contains__(self, key):
        """
        Allows usage of the 'in' keyword.
        - If int: checks if index is in bounds.
        - If not int: checks if an 'my_name' matches the key.
        """
        if isinstance(key, int):
            # Check if index is within the valid range of the list
            return 0 <= key < len(self.data)
        
        # Check if any dictionary has an 'my_name' matching the key
        return any(item.get("my_name") == key for item in self.data)

    def __getitem__(self, key):
        """
        Allows accessing the list by standard index (int) 
        or by the 'my_name' key of the dictionaries.
        """
        # If the key is an integer, access by list index
        if isinstance(key, int):
            return self.data[key]
        
        # Otherwise, search through the dictionaries for a matching 'my_name'
        for item in self.data:
            if item.get("my_name") == key:
                return item
                
        # If no match is found, raise a standard KeyError
        raise KeyError(f"No dictionary found with my_name: '{key}'")

    def __len__(self):
        return len(self.data)
    
    def __iter__(self):
        # reset iteration each time
        self._current_item = 0
        return self

    def __next__(self):
        if self._current_item < len(self):
            self._current_item += 1
            return self.data[self._current_item - 1]
        raise StopIteration

    def __repr__(self):
        return f"[\n{'    '}{'\n    '.join(uid for uid in [item.get("my_name") for item in self.data])}\n]"
    
    def add(self,
            my_name:                 str,
            front_side:              dict,
            back_side:               dict,
            _my_type:                str  = None,
            _refraction_index_front: float = None,
            _refraction_index_back:  float = None) -> str:
        # default values if args isn't given
        _my_type                = default_fallback(_my_type, K_MY_TYPE)
        _refraction_index_front = default_fallback(_refraction_index_front, K_REFRACTION_INDEX_FRONT)
        _refraction_index_back  = default_fallback(_refraction_index_back, K_REFRACTION_INDEX_BACK)

        # assertion checks for validity
        assert type(my_name) == str,                                             'my_name must be a string'
        assert type(front_side) == dict,                                         'front_side must be a dictionary'
        assert type(back_side) == dict,                                          'back_side must be a dictionary'
        assert OPTICAL_INTERACTIONS_ENUM[_my_type] in OPTICAL_INTERACTION_ASSERT, f'invalid optical interaction type {_my_type}'
        assert type(_refraction_index_front) == float,                           'refraction index front must be a float'
        assert type(_refraction_index_back) == float,                            'refraction index back must be a float'

        self.data.append({
            K_MY_NAME: my_name,
            K_MY_TYPE: _my_type,
            K_REFRACTION_INDEX_FRONT: _refraction_index_front,
            K_REFRACTION_INDEX_BACK: _refraction_index_back,
            K_FRONT: front_side,
            K_BACK: back_side,
        })
        return my_name
    
    def save_to_json(self):
        return { str(i): item for i, item in enumerate(self.data) }
    
    def index(self, key):
        """
        Returns the index of the dictionary with 'my_name' matching the key.
        Raises ValueError if no match is found.
        """
        for i, item in enumerate(self.data):
            if item.get("my_name") == key:
                return i
        raise ValueError(f"No dictionary found with my_name: '{key}'")
    
    def clear(self):
        self.data = []

"""genertic element dictionary constructor finish this after writing csv anaysis code"""
def element_json(my_name:          str,
                 origin:           list,
                 aim:              list,
                 zrot:             float,
                 aperture:         dict,
                 surface:          dict,
                 optical_ref:      int | str,
                 optical_registry: OpticalPropertyRegistry,
                 group:            int = -1) -> dict:
    assert type(my_name) == str,               "my_name must be a string"
    assert type(origin) == list,               "origin must be a list"
    assert type(aim) == list,                  "aim must be a list"
    assert type(zrot) == float,                "zrot must be a float"
    assert type(aperture) == dict,             "aperture must be a dictionary"
    assert type(surface) == dict,              "surface must be a dictionary"
    assert len(origin) == 3,                   "origin must be a 3D vector"
    assert len(aim) == 3,                      "aim must be a 3D vector"
    assert optical_ref in optical_registry,    f"optical_ref {optical_ref} not found in registry"
    assert type(group) == int and group >= -1, "group must be an integer greater than or equal to -1"

    return {
        K_IS_SINGLE:    DEFAULTS[K_IS_SINGLE],
        K_ACTIVE:       DEFAULTS[K_ACTIVE],
        K_VIRTUAL_FLAG: DEFAULTS[K_VIRTUAL_FLAG],
        K_MY_ID:        DEFAULTS[K_MY_ID],
        K_STAGE:        DEFAULTS[K_STAGE],
        K_MY_NAME:      my_name,
        K_GROUP:        group,
        K_ORIGIN:       origin,
        K_AIM:          aim,
        K_ZROT:         zrot,
        K_APERTURE:     aperture,
        K_SURFACE:      surface,
        K_OPT_ID:       optical_registry.index(optical_ref) if type(optical_ref) == str else optical_ref,
    }
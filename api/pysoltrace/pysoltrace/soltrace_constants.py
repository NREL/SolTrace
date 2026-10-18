import re

# chedder finds constants defined in stapi_v2.h and exposes
# them as dot_h. soltrace_constants holds constant information
# for the Python side of SolTrace, i.e. names of enum values
# that aren't from IntEnums, error/warning messages, or types
# defined in stapi_v2.h that are used in type hinting.
from enum import Enum, IntEnum
from pysoltrace import dot_h

# byte character / int enums for backwards compatibility
class optical_error_dist(Enum):
    GAUSSIAN     = b'g'
    PILLBOX      = b'p'
    DIFFUSE      = b'f'
    USER_DEFINED = b'd'

class sun_shape(Enum):
    GAUSSIAN     = b'g'
    PILLBOX      = b'p'
    BUIE_CSR     = b'b'
    USER_DEFINED = b'd'

class optical_interaction(IntEnum):
    REFRACTION = 1
    REFLECTION = 2

class aperture(Enum):
    CIRCLE                        = b'c'
    HEXAGON                       = b'h'
    EQUILATERAL_TRIANGLE          = b't'
    RECTANGLE                     = b'r'
    ANNULUS                       = b'a'
    SINGLE_AXIS_CURVATURE_SECTION = b'l'
    IRREGULAR_TRIANGLE            = b'i'
    IRREGULAR_QUADRILATERAL       = b'q'

class surface(Enum):
    SPHERE                = b's'
    PARABOLA              = b'p'
    HYPER                 = b'o'
    GENERAL_SPENCER_MURTY = b'g'
    FLAT                  = b'f'
    CONE                  = b'c'
    CYLINDER              = b't'
    TORUS                 = b'd'

######################################
# st_return_code name/message set up #
######################################

ST_RETURN_CODE_NAME = { code.value: code.name for code in dot_h.st_return_code }

# messages for return codes
ST_RETURN_CODE_ERROR_MSG = {
    dot_h.st_return_code.FAILURE:                                      '',
    dot_h.st_return_code.CANCEL:                                       'Simulation canceled.',
    dot_h.st_return_code.CONTEXT_NOT_FOUND:                            'Context pointer could not be cast as context.',
    dot_h.st_return_code.DATA_NOT_FOUND:                               'SimulationData pointer could not be found in context.',
    dot_h.st_return_code.RUNNER_NOT_FOUND:                             'SimulationRunner pointer could not be found in context.',
    dot_h.st_return_code.RESULT_NOT_FOUND:                             'SimulationResult pointer could not be found in context.',
    dot_h.st_return_code.INVALID_ARGUMENTS:                            'Invalid arguments for stapi_v2 function.',
    dot_h.st_return_code.DATA_INSERTION_FAILURE:                       'Unable to insert data into SimulationData.',
    dot_h.st_return_code.DATA_VALUE_NOT_FOUND:                         'Data required for function was not found in SimulationData',
    dot_h.st_return_code.RUNNER_INILIALIZE_FAILURE:                    'SimulationRunner could not be initialized.',
    dot_h.st_return_code.RUNNER_NUMBER_THREADS_SEEDS_MISMATCH_FAILURE: 'Number of threads requested and length of seeds list are not '
                                                                       'equal. Include the same number of seeds as threads requested.',
    dot_h.st_return_code.RUNNER_SETUP_ERROR:                           'A runtime error was thrown when setting up SimulationRunner.',
    dot_h.st_return_code.RUNNER_SETUP_FAILURE:                         'SimulationRunner could not be set up based on SimulationData provided.',
    dot_h.st_return_code.RUNNER_NOT_READY_TO_RUN:                      'SimulationRunner is not ready to run. Set up the SimulationRunner.',
    dot_h.st_return_code.RUNNER_NOT_READY_TO_REPORT:                   'SimulationRunner is not ready to report. Run the simulation.',
    dot_h.st_return_code.RESULT_NOT_REPORTED:                          'Reporting level does not include this information. Change to appropriate level.',
    dot_h.st_return_code.EXCEPTION:                                    'Exception raised. Check validity of the arguments passed to the '
                                                                       'function called or of the JSON against SolTrace schema version used.',
    dot_h.st_return_code.UKNOWN_BATCH_API_CALL_FAILURE:                'Unknown batch call received.'
}
ST_RETURN_CODE_WARNING_MSG = {
    dot_h.st_return_code.WARNING_FELLBACK_FROM_EMBREE:       'Requested EmbreeRunner, but is not installed. Fellback to NativeRunner.',
    dot_h.st_return_code.WARNING_FELLBACK_FROM_OPTIX:        'Requested OptixRunner, but is not installed. Fellback to NativeRunner.',
    dot_h.st_return_code.WARNING_ARGUMENT_IGNORED_BY_RUNNER: 'Requested a number of threads for OptixRunner. The argument does not apply '
                                                             'to this runner type and was ignored.',
    dot_h.st_return_code.WARNING_SUN_SHAPE_IGNORED:          'Requested invalid sun shape. Default sun created with a Gaussian '
                                                             'distribution with sigma = 4.65.',
    dot_h.st_return_code.WARNING_GROUP_IGNORED:              'Requested invalid group number. Setting element to ungrouped.',
    dot_h.st_return_code.WARNING_NOT_FOUND:                  'Requested an item that was not found.',
}

# used for extracting ctypes class name from str(ctypes obj)
_CTYPES_RE = re.compile(r"(?<=<class 'ctypes\.)[A-Za-z_0-9.<>]+(?='>)")


import ctypes, os, pathlib, sys
from pysoltrace import dot_h

# find cpu dll, overritten with gpu dll finder if found 
def find_dll() -> pathlib.Path:
    _path = pathlib.Path(__file__).parent.parent.resolve() / 'bin'
    assert _path.exists(), f"Could not find SolTrace DLL at {_path}"
                
    if sys.platform == "win32":
        _lib_name = "stapi_v2.dll"
    elif sys.platform == "darwin":
        _lib_name = "libstapi_v2.dylib"
    else:
        _lib_name = "libstapi_v2.so"

    return _path / _lib_name

_find_optix_dll = None
try:
    # search for nvsoltrace in production environment
    from nvsoltrace import find_dll as _nvfd
    _find_optix_dll = _nvfd
except (ImportError, ModuleNotFoundError):
    # search for nvsoltrace in development environment
    try:
        nvsoltrace_path = pathlib.Path(__file__).parent.parent.parent.parent / 'nvsoltrace' / 'nvsoltrace'
        sys.path.insert(1, str(nvsoltrace_path))
        import find_dll as _nvfd # pyright: ignore[reportMissingImports]
        _find_optix_dll = _nvfd.find_dll
    except (ImportError, ModuleNotFoundError): pass # silently fail searching

if _find_optix_dll != None:
    try:
        _find_optix_dll()
        find_dll = _find_optix_dll
    except AssertionError: pass

STATIC_FUNC_ARGS = ['args_st_get_installed_runners', 'args_st_is_runner_installed']

def __get_argtypes(args_struct: ctypes.Structure) -> list:
    return [dot_h.st_context_v2_t, *[args[1] for args in args_struct._fields_]]

def __get_dll(path: str = '') -> ctypes.CDLL:
    if not os.path.exists(path):
        raise FileNotFoundError(f'Could not find DLL at {path}')

    if sys.platform == "win32":
        os.add_dll_directory(str(path).rsplit(os.sep, maxsplit=1)[0])
        return ctypes.WinDLL(path, winmode=0)
    else: return ctypes.CDLL(path)

def setup_dll(path: str = ''):
    pdll = __get_dll(path)

    #############################################
    # functions for SolTrace context management #
    #############################################

    # warning that: "The function "POINTER" is deprecated ctypes.POINTER
    # with string" dot_h.st_context_v2_t is <class 'ctypes.c_void_p'>
    # which is not deprecated see ctypes entry for more information on
    # the following Python docs page here:
    # https://docs.python.org/3/whatsnew/3.14.html#new-deprecations 
    pdll.st_create_context.argtypes = [ctypes.POINTER(dot_h.st_context_v2_t),
                                       ctypes.CFUNCTYPE(ctypes.c_int,
                                                        ctypes.c_char_p,
                                                        ctypes.c_char_p)]
    pdll.st_create_context.restype  = dot_h.st_return_t

    pdll.st_reset_context.argtypes = [ctypes.c_void_p]
    pdll.st_reset_context.restype  = dot_h.st_return_t

    pdll.st_free_context.argtypes = [ctypes.c_void_p]
    pdll.st_free_context.restype  = dot_h.st_return_t
    
    pdll.st_locate_context.argtypes = [ctypes.c_void_p, ctypes.c_char_p]
    pdll.st_locate_context.restype  = dot_h.st_return_t

    ###############################################
    # creates ctypes function pointers for each   #
    # of the structs defined in stapi_v2.h        #
    # starting with "args_st" these functions do: #
    # - SolTrace data management                  #
    #   - thru json strings                       #
    #   - directly                                #
    #     - set simulation parameters             #
    #     - add/remove/set optical properties     #
    #     - add/remove/modify elements            #
    #     - add/modify the sun                    #
    # - input files for SolTrace writing          #
    # - SolTrace runner management                #
    # - SolTrace results management               #
    #   - thru writing files                      #
    #   - directly                                #
    ###############################################
    for k, v in vars(dot_h).items():
        if k.startswith("args_st_"):
            func_name = k.replace("args_", "")
            func = getattr(pdll, func_name)
            func.argtypes = __get_argtypes(v)
            func.restype  = dot_h.st_return_t
        
    ############################################
    # function for batching SolTrace API calls #
    ############################################
    
    pdll.st_batch.argtypes = [ctypes.c_void_p,
                              ctypes.POINTER(ctypes.c_void_p),
                              ctypes.c_uint,
                              ctypes.POINTER(ctypes.c_uint),
                              ctypes.c_bool]
    pdll.st_batch.restype  = dot_h.st_return_t
    return pdll

def setup_static_dll(path: str = ''):
    pdll = __get_dll(path)

    for f_args in STATIC_FUNC_ARGS:
        func_name = f_args.replace('args_', '')
        args_struct = getattr(dot_h, f_args)
        func = getattr(pdll, func_name)
        func.argtypes = [args[1] for args in args_struct._fields_]
        func.restype  = dot_h.st_return_t

    return pdll

STATIC_DLL = setup_static_dll(find_dll())

class context:
    __slots__ = ('_pdll', '_pcxt')
    
    def __init__(self, pdll, pcxt):
        self._pdll: ctypes.CDLL = pdll
        self._pcxt: int = pcxt
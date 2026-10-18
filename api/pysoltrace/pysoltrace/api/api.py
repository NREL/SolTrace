import atexit, os, pathlib, sys, warnings, weakref
from typing import Literal
from pathlib import Path
import orjson
import ctypes
from functools import partial
from colorama import just_fix_windows_console, Fore, Back, Style
just_fix_windows_console()

from pysoltrace.api.dll import context, setup_dll as _setup_dll, find_dll as _find_dll

# search for nvsoltrace hook into optix runner

from pysoltrace import soltrace_constants as _STC
from pysoltrace.api.utils import check_return_code, st_function, STAPIv2Exception, STAPIv2Warning
from pysoltrace.api.runner import runner
from pysoltrace.api.data import data
from pysoltrace.api.parameters import parameters
from pysoltrace.api.legacy import legacy
from pysoltrace.api.batch import batch
from pysoltrace.api.result import result

def free(dll, pcxt, testing: bool = False):
    code = dll.st_free_context(pcxt)
    if not testing:
        sys.stdout.write(f'Freed context ({pcxt:#x}) with code ({code}) from SolTrace DLL ({dll})\n')

@ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_char_p, ctypes.c_char_p)
def message_cb(loc, msg):
    sys.stdout.write(f"{Fore.MAGENTA}[stapi_v2] - Message callback triggered by ({loc.decode('utf-8')}){Style.RESET_ALL}: {msg.decode('utf-8')}\n")
    return 0

@ctypes.CFUNCTYPE(ctypes.c_int, ctypes.c_char_p, ctypes.c_char_p)
def testing_cb(loc, msg): return 0

#############################################################################
# STAPIv2 Class: wraps stapi_v2.{dll, so, dylib} with more Python-ish calls #
#############################################################################
class STAPIv2:
    __slots__ = ('parameters',
                 'data',
                 'runner',
                 'result',
                 'batch',
                 'legacy',
                 '__pcxt',
                 '__pdll',
                 '__benchmarking')
    
    def __init__(self, override_path: str = '', testing: bool = False, benchmarking: bool = False):
        self.__pdll, self.__pcxt = self.__create(override_path, testing)

        atexit.register(free, self.__pdll, self.__pcxt, testing)

        self.parameters = parameters(self.__pdll, self.__pcxt)
        self.data       = data(self.__pdll, self.__pcxt)
        self.runner     = runner(self.__pdll, self.__pcxt)
        self.result     = result(self.__pdll, self.__pcxt)
        self.batch      = batch(self.__pdll, self.__pcxt)
        self.legacy     = legacy(self.__pdll, self.__pcxt)

        self.__benchmarking = benchmarking

    def __repr__(self):
        rt = f'STAPIv2 Object at ({id(self.__pdll):#x})'
        rt += '\n' + '-' * len(rt) + '\n'

        runners = self.runner.get_installed()
        rt += f'Installed Runners: ' + ', '.join([k for k, v in runners.items() if v]) + '\n'
        max_key_len = max(len(str(k)) for k in self.__pdll.__dict__.keys())
        for k, v in self.__pdll.__dict__.items():
            rt += f'\n{str(k):<{max_key_len}}: {v}'
            if (isinstance(v, ctypes._CFuncPtr)):
                rt += f'\n{" " * max_key_len}: {[_STC._CTYPES_RE.search(str(arg)).group() for arg in v.argtypes]}'
        return rt

    def __create(self, override_path: str = '', testing: bool = False):
        dll_f = pathlib.Path(override_path).resolve() if len(override_path) else _find_dll()
        pdll = _setup_dll(dll_f)

        ppcxt = ctypes.c_void_p()
        code = pdll.st_create_context(ctypes.byref(ppcxt), message_cb if not testing else testing_cb)
        check_return_code(code)

        code = pdll.st_locate_context(ppcxt.value, bytes(dll_f.parent))
        check_return_code(code)

        return pdll, ppcxt.value

    def reset(self):
        code = self.__pdll.st_reset_context(self.__pcxt)
        check_return_code(code)
        if not self.__benchmarking:
            sys.stdout.write(f'Reset context ({self.__pcxt:#x})\n')

    @staticmethod
    def is_runner_installed(type: int) -> bool:
        return runner.is_installed(type)

    @staticmethod
    def get_installed_runners() -> dict[str, bool]:
        return runner.get_installed()
    
    def sneak(self): return self.__pdll, self.__pcxt, check_return_code
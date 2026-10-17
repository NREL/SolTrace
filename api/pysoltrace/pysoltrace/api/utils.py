import functools, warnings
from colorama import just_fix_windows_console, Fore, Back, Style
just_fix_windows_console()

from pysoltrace import soltrace_constants as _STC
from pysoltrace import dot_h

###################################
# exception and warning reporting #
###################################
class STAPIv2Exception(Exception):
    def __init__(self, code, name, msg) -> None:
        super().__init__(f'{Fore.RED}[stapi_v2] - Call returned with error code ({code}: {name}).{Style.RESET_ALL}\n  {msg}')
        self.code = code
        self.name = name
        self.message = msg

STAPI_V2_WARNING_PREFIX = '[stapi_v2] - Call returned with warning code'
STAPIv2Warning = lambda code, name, msg: warnings.warn(
    f'{Fore.YELLOW}{STAPI_V2_WARNING_PREFIX} ({code}: {name}).{Style.RESET_ALL}\n  {msg}',
    stacklevel=4
)

def _get_msg(code, msgs, batch_i = None):
    name = _STC.ST_RETURN_CODE_NAME[code]
    if batch_i == None: return name, msgs[code]
    else:
        batch_i_str = f'{Fore.RED}Batch call failed on call number {batch_i}.{Style.RESET_ALL}\n'
        return name, batch_i_str + msgs[code]

def check_return_code(code, batch_i = None):
    if code == dot_h.st_return_code.SUCCESS: return
    elif code > dot_h.st_return_code.SUCCESS and code < dot_h.st_return_code.WARNING_FELLBACK_FROM_EMBREE:
        raise STAPIv2Exception(code, *_get_msg(code, _STC.ST_RETURN_CODE_ERROR_MSG, batch_i))
    elif code >= dot_h.st_return_code.WARNING_FELLBACK_FROM_EMBREE and code < dot_h.st_return_code.RETURN_COUNT:
        STAPIv2Warning(code, *_get_msg(code, _STC.ST_RETURN_CODE_WARNING_MSG, batch_i))
    elif code >= dot_h.st_return_code.RETURN_COUNT:
        raise STAPIv2Exception(code, 'UNKNOWN', 'Unknown return code received.')

def st_function(func):
    @functools.wraps(func)
    def wrapper(*args, **kwargs):
        code_rt = func(*args, **kwargs)
        if isinstance(code_rt, tuple):
            check_return_code(code_rt[0])
            # if need to return a function to produce return args
            # assert len(rt) <= 1, f'Expected at most one return value, got {len(rt)}'
            # return rt[0]() if len(rt) == 1 else None

            return code_rt[1] if len(code_rt) == 2 else code_rt[1:] if len(code_rt) > 2 else None
        else: check_return_code(code_rt)
    return wrapper
